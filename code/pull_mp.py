"""
pull_mp.py — Week 1, stage 1: raw Materials Project pull for Paper 2
====================================================================
(Efficiency–coverage trade-off of scientific priors in materials discovery)

Pulls EVERY elasticity document plus matching summary documents and writes
an immutable, hashed raw snapshot. It applies NO scientific filters: the
metallicity, element-list, N>=2, Born-stability and reference-table
coverage filters all live in build_pool.py, so the candidate universe can be
audited and revised without re-pulling.

Verified against mp-api 0.46.5 / emmet-core 0.87.2 (field names, search
signatures, MPRester.db_version, MPRestError messages).

Procedure (pre-registered)
--------------------------
    1. install the pinned environment
    2. python pull_mp.py --inspect > inspect_output.txt   # archive this file
    3. python pull_mp.py                                  # full pull
    4. record mp_db_version + the four data-file SHA-256s from manifest.json
       in the pre-registration, and deposit the snapshot on Zenodo

Outputs (data/raw/mp_<db_version>/), all read-only after writing
-------
    elasticity.jsonl.gz   one record per RETURNED elasticity document,
                          structure removed. material_id uniqueness is
                          audited in manifest.json; duplicates are retained
                          and resolved downstream.
    structures.jsonl.gz   {material_id, structure} from the ELASTICITY doc
    summary.jsonl.gz      every returned summary document, as returned
    index.csv             flat convenience view (derived; never the source;
                          build_pool.py must read the .jsonl.gz files)
    manifest.json         acquisition metadata: db version, package versions,
                          counts, ID-integrity audit, SHA-256s

Determinism: records are sorted by material_id, keys are sorted, and gzip is
written with mtime=0, so the same acquisition against the same database
version and API state produces identical hashes for the FOUR DATA ARTIFACTS.
manifest.json holds acquisition metadata (timestamp, runtime) and is not
expected to hash identically across runs.

Immutability: an existing snapshot directory is never overwritten, and files
are made read-only. Read-only permissions only prevent accidents; the actual
guarantee is external (hashes in the pre-registration + Zenodo deposit).
Server-side changes made without a db-version bump cannot be detected by
any client and are a documented limitation.
"""

import argparse
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import platform
import re
import stat
import time
from collections import Counter
from importlib import metadata
from pathlib import Path

ELASTICITY_FIELDS = [
    "material_id", "formula_pretty", "composition", "composition_reduced",
    "chemsys", "elements", "nelements", "nsites", "volume", "density",
    "density_atomic", "symmetry", "structure", "order", "elastic_tensor",
    "bulk_modulus", "shear_modulus", "homogeneous_poisson",
    "universal_anisotropy", "state", "warnings", "deprecated",
    "deprecation_reasons", "fitting_method", "last_updated", "origins",
]
# Deliberately excluded: fitting_data (large raw stress/strain sets),
# compliance_tensor (derivable from elastic_tensor), sound_velocity,
# thermal_conductivity, debye_temperature (not used by the design).

SUMMARY_FIELDS = [
    "material_id", "formula_pretty", "is_metal", "band_gap", "density",
    "nelements", "elements", "symmetry", "energy_above_hull", "is_stable",
    "formation_energy_per_atom", "theoretical", "deprecated",
    "is_magnetic", "ordering",
]

SUMMARY_CHUNK = 500      # material_ids per summary request
RETRIES = 4


# ─────────────────────────────── helpers ────────────────────────────────

def _api_key():
    key = os.environ.get("MP_API_KEY", "").strip()
    if not key:
        raise SystemExit(
            "MP_API_KEY not set. In Colab: sidebar key icon -> add MP_API_KEY, then\n"
            "  from google.colab import userdata; import os\n"
            "  os.environ['MP_API_KEY'] = userdata.get('MP_API_KEY')")
    return key


def to_plain(obj):
    """Recursively convert API output (dicts or pydantic docs) to JSON types.
    Non-finite floats are preserved using Python's JSON representation
    (NaN/Infinity, which strict JSON parsers reject); build_pool.py must
    handle them explicitly when deciding eligibility."""
    if obj is None or isinstance(obj, (bool, int, str)):
        return obj
    if isinstance(obj, float):
        return obj
    if isinstance(obj, dict):
        return {str(k): to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_plain(v) for v in obj]
    if hasattr(obj, "model_dump"):                      # pydantic v2 doc
        return to_plain(obj.model_dump(mode="json"))
    if hasattr(obj, "as_dict"):                         # pymatgen object
        return to_plain(obj.as_dict())
    if hasattr(obj, "tolist"):                          # numpy
        return to_plain(obj.tolist())
    if hasattr(obj, "value") and not callable(obj.value):   # Enum
        return to_plain(obj.value)
    if isinstance(obj, (dt.datetime, dt.date)):
        return obj.isoformat()
    return str(obj)                                     # MPID etc.


# Messages mp-api 0.46.5 uses for errors that retrying cannot fix
# (invalid fields, HTTP 400 "server does not support", bad chunk arguments).
_PERMANENT = ("invalid fields requested", "does not support the request",
              "Chunk size must", "Number of chunks must")


def is_transient(e):
    """Retry only network/server failures. Programming errors (TypeError,
    KeyError, ...) and request-validation errors fail immediately."""
    from requests.exceptions import RequestException
    from mp_api.client.core import MPRestError
    if isinstance(e, RequestException):
        return True
    if not isinstance(e, MPRestError):
        return False
    msg = str(e)
    if any(p in msg for p in _PERMANENT):
        return False
    code = re.search(r"status code (\d{3})", msg)
    if code:                                            # HTTP error reply
        c = int(code.group(1))
        return c == 429 or c >= 500
    return True        # wrapped RequestException / timeout: free-form message


def with_retries(fn, what):
    for attempt in range(1, RETRIES + 1):
        try:
            return fn()
        except Exception as e:
            if attempt == RETRIES or not is_transient(e):
                raise
            wait = 10 * 2 ** (attempt - 1)
            print(f"  ! {what} failed ({type(e).__name__}: {e}); "
                  f"retry {attempt}/{RETRIES - 1} in {wait}s")
            time.sleep(wait)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_jsonl_gz(path, records):
    """Deterministic: sorted records, sorted keys, gzip mtime=0."""
    records = sorted(records, key=lambda r: r["material_id"])
    with open(path, "wb") as raw, \
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        for r in records:
            gz.write((json.dumps(r, sort_keys=True, separators=(",", ":"))
                      + "\n").encode("utf-8"))


def pkg_version(name):
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def _get(d, *keys):
    for k in keys:
        if not isinstance(d, dict):
            return None
        d = d.get(k)
    return d


def _num(x):
    return x if isinstance(x, (int, float)) and not isinstance(x, bool) else None


def get_db_version(mpr):
    """mp-api returns "" when the heartbeat is refused (HTTP 403)."""
    v = mpr.db_version
    if not v:
        raise SystemExit("Could not read the MP database version (empty "
                         "heartbeat response; check API key/access). Aborting.")
    return v


def checkpoint_fingerprint(db_version):
    """Everything that must be identical for a checkpoint to be reusable."""
    return {"mp_db_version": db_version,
            "script_sha256": sha256(__file__),
            "elasticity_fields": ELASTICITY_FIELDS,
            "summary_fields": SUMMARY_FIELDS,
            "elasticity_chunk_size": 1000,
            "summary_chunk": SUMMARY_CHUNK,
            "mp-api": pkg_version("mp-api"),
            "emmet-core": pkg_version("emmet-core")}


def open_checkpoint_dir(partial, db_version):
    """Create the checkpoint dir with its fingerprint, or verify an existing
    one. Refuses to resume if anything in the fingerprint differs."""
    meta = partial / "checkpoint_meta.json"
    want = checkpoint_fingerprint(db_version)
    if partial.exists():
        have = json.loads(meta.read_text()) if meta.exists() else None
        if have != want:
            diff = sorted(k for k in want if (have or {}).get(k) != want[k])
            raise SystemExit(f"Checkpoint {partial} was produced under a "
                             f"different state (differs in: {diff or 'no metadata'}). "
                             "Refusing to resume; delete it and re-run.")
        print(f"  resuming from verified checkpoint {partial}")
    else:
        partial.mkdir(parents=True)
        meta.write_text(json.dumps(want, indent=2))


def id_audit(ids_requested, docs):
    got = [d["material_id"] for d in docs]
    c = Counter(got)
    return {"docs": len(got), "unique_ids": len(c),
            "missing_ids": sorted(set(ids_requested) - set(c)),
            "unexpected_ids": sorted(set(c) - set(ids_requested)),
            "duplicate_ids": sorted(k for k, n in c.items() if n > 1)}


def structure_check(sd):
    """Round-trip a structure dict through pymatgen; report what it carries."""
    from pymatgen.core import Structure
    s = Structure.from_dict(sd)
    oxi = any(getattr(sp, "oxi_state", None) not in (None, 0)
              for site in s for sp in site.species)
    return {"formula": s.composition.reduced_formula, "nsites": len(s),
            "lattice_abc": [round(x, 4) for x in s.lattice.abc],
            "lattice_angles": [round(x, 3) for x in s.lattice.angles],
            "oxidation_states": oxi,
            "site_properties": sorted(s.site_properties)}


# ──────────────────────────────── fetch ─────────────────────────────────

def open_rester(key):
    from mp_api.client import MPRester
    # use_document_model=False -> plain dicts: no pydantic round-trip, and
    # server fields unknown to this emmet version are not dropped.
    return MPRester(key, use_document_model=False)


def fetch_elasticity(mpr, num_chunks=None, chunk_size=1000):
    return with_retries(
        lambda: mpr.materials.elasticity.search(
            fields=ELASTICITY_FIELDS, num_chunks=num_chunks,
            chunk_size=chunk_size),
        "elasticity search")


def fetch_summary(mpr, ids, partial_dir=None):
    out = []
    chunks = [ids[i:i + SUMMARY_CHUNK] for i in range(0, len(ids), SUMMARY_CHUNK)]
    for n, chunk in enumerate(chunks):
        ckpt = partial_dir / f"summary_{n:04d}.json" if partial_dir else None
        if ckpt and ckpt.exists():                    # resume after a dropout
            out.extend(json.loads(ckpt.read_text()))
            continue
        docs = with_retries(
            lambda c=chunk: mpr.materials.summary.search(
                material_ids=c, fields=SUMMARY_FIELDS),
            f"summary chunk {n + 1}/{len(chunks)}")
        recs = [to_plain(d) for d in docs]
        if ckpt:
            ckpt.write_text(json.dumps(recs))
        out.extend(recs)
        print(f"  summary chunk {n + 1}/{len(chunks)}: {len(recs)} docs")
    return out


# ─────────────────────────────── inspect ────────────────────────────────

def inspect_mode(key):
    """Fetch 3 docs per endpoint and show what the SERVER actually returns."""
    with open_rester(key) as mpr:
        print(f"MP database version: {get_db_version(mpr)}")
        el = [to_plain(d) for d in fetch_elasticity(mpr, num_chunks=1, chunk_size=3)]
        el = el[:3]
        ids = [d["material_id"] for d in el]
        sm = fetch_summary(mpr, ids)
    for name, docs, req in (("ELASTICITY", el, ELASTICITY_FIELDS),
                            ("SUMMARY", sm, SUMMARY_FIELDS)):
        got = set().union(*(d.keys() for d in docs)) if docs else set()
        print(f"\n{name}: {len(docs)} docs")
        print(f"  requested but NOT returned: {sorted(set(req) - got) or 'none'}")
        print(f"  returned but not requested: {sorted(got - set(req)) or 'none'}")
        for k in sorted(got):
            v = docs[0].get(k)
            s = json.dumps(v)[:110] if k != "structure" else "<structure dict>"
            print(f"    {k:<24}{type(v).__name__:<8}{s}")
    print("\nSTRUCTURE round-trip (pymatgen Structure.from_dict):")
    for d in el:
        try:
            print(f"  {d['material_id']}: {structure_check(d.get('structure'))}")
        except Exception as e:
            print(f"  {d['material_id']}: FAILED ({type(e).__name__}: {e})")
    print("\nCheck in particular: shear_modulus/bulk_modulus have voigt/reuss/vrh,"
          "\nelastic_tensor has ieee_format (6x6), state is a string, and"
          "\nsummary.is_metal is present. Nothing was written to disk.")


# ──────────────────────────────── full pull ─────────────────────────────

def flat_row(e, s, n_summary, n_elast):
    """s is the summary doc only when exactly one was returned for this id;
    summary columns are blank otherwise (see summary_n_docs)."""
    sym = e.get("symmetry") or {}
    return {
        "material_id": e["material_id"],
        "formula_pretty": e.get("formula_pretty"),
        "nelements": e.get("nelements"),
        "nsites": e.get("nsites"),
        "G_vrh": _num(_get(e, "shear_modulus", "vrh")),
        "G_voigt": _num(_get(e, "shear_modulus", "voigt")),
        "G_reuss": _num(_get(e, "shear_modulus", "reuss")),
        "K_vrh": _num(_get(e, "bulk_modulus", "vrh")),
        "density_elast": _num(e.get("density")),
        "spacegroup_number": sym.get("number"),
        "crystal_system": sym.get("crystal_system"),
        "state": e.get("state"),
        "n_warnings": len(e.get("warnings") or []),
        "deprecated_elast": e.get("deprecated"),
        "has_ieee_tensor": _get(e, "elastic_tensor", "ieee_format") is not None,
        "n_elasticity_docs_for_id": n_elast,
        "summary_n_docs": n_summary,
        "is_metal": (s or {}).get("is_metal"),
        "band_gap": _num((s or {}).get("band_gap")),
        "energy_above_hull": _num((s or {}).get("energy_above_hull")),
        "theoretical": (s or {}).get("theoretical"),
        "deprecated_summary": (s or {}).get("deprecated"),
    }


def descriptive_audit(rows):
    """Counts only — nothing here filters the snapshot."""
    n = len(rows)
    g = [r["G_vrh"] for r in rows]
    finite = [x for x in g if x is not None and math.isfinite(x)]
    print("\n" + "=" * 64 + "\n  DESCRIPTIVE AUDIT (no filtering applied)\n" + "=" * 64)
    print(f"  elasticity docs                 : {n}")
    print(f"  summary docs per id (0/1/2+)    : "
          f"{sum(r['summary_n_docs'] == 0 for r in rows)} / "
          f"{sum(r['summary_n_docs'] == 1 for r in rows)} / "
          f"{sum(r['summary_n_docs'] > 1 for r in rows)}")
    print(f"  state counts                    : {dict(Counter(r['state'] for r in rows))}")
    print(f"  with >=1 warning                : {sum(r['n_warnings'] > 0 for r in rows)}")
    print(f"  deprecated (elasticity doc)     : {sum(bool(r['deprecated_elast']) for r in rows)}")
    print(f"  missing ieee elastic tensor     : {sum(not r['has_ieee_tensor'] for r in rows)}")
    print(f"  G_vrh missing / non-finite      : {n - len(finite)}")
    print(f"  G_vrh <= 0                      : {sum(x <= 0 for x in finite)}")
    print(f"  G_vrh >= 600 (Paper-1 cutoff)   : {sum(x >= 600 for x in finite)}")
    print(f"  is_metal True / False / missing : "
          f"{sum(r['is_metal'] is True for r in rows)} / "
          f"{sum(r['is_metal'] is False for r in rows)} / "
          f"{sum(r['is_metal'] is None for r in rows)}")
    print(f"  nelements distribution          : "
          f"{dict(sorted(Counter(r['nelements'] for r in rows).items()))}")
    print(f"  is_metal & nelements>=2         : "
          f"{sum(r['is_metal'] is True and (r['nelements'] or 0) >= 2 for r in rows)}"
          "   (upper bound on the pool; element list and stability not yet applied)")


def full_pull(key, out_root):
    t0 = time.time()
    with open_rester(key) as mpr:
        db_start = get_db_version(mpr)
        print(f"MP database version: {db_start}")
        out = Path(out_root) / f"mp_{db_start}"
        if out.exists():
            raise SystemExit(f"{out} already exists. Refusing to overwrite an "
                             "existing snapshot. For a re-acquisition, use a "
                             "different --out-dir.")
        partial = Path(out_root) / f".partial_mp_{db_start}"
        Path(out_root).mkdir(parents=True, exist_ok=True)
        open_checkpoint_dir(partial, db_start)

        print("Pulling elasticity documents (all chunks)...")
        el_ckpt = partial / "elasticity.json"
        if el_ckpt.exists():
            elastic = json.loads(el_ckpt.read_text())
            print(f"  resumed {len(elastic)} docs from checkpoint")
        else:
            elastic = [to_plain(d) for d in fetch_elasticity(mpr)]
            el_ckpt.write_text(json.dumps(elastic))
        ids = sorted({d["material_id"] for d in elastic})
        if len(ids) != len(elastic):
            print(f"  ! {len(elastic) - len(ids)} duplicate material_ids in "
                  "elasticity docs (kept all; audited in manifest)")

        print(f"Pulling summary documents for {len(ids)} material_ids...")
        summary = fetch_summary(mpr, ids, partial_dir=partial)

    with open_rester(key) as mpr2:                     # release changed mid-pull?
        db_end = get_db_version(mpr2)
    if db_end != db_start:
        raise SystemExit(f"MP database changed during the pull ({db_start} -> "
                         f"{db_end}). Delete {partial} and re-run.")

    structures = [{"material_id": d["material_id"], "structure": d.get("structure")}
                  for d in elastic]
    elastic_ns = [{k: v for k, v in d.items() if k != "structure"} for d in elastic]
    el_audit = id_audit(ids, elastic)
    sm_audit = id_audit(ids, summary)
    s_count = Counter(d["material_id"] for d in summary)
    e_count = Counter(d["material_id"] for d in elastic)
    s_one = {d["material_id"]: d for d in summary if s_count[d["material_id"]] == 1}
    rows = [flat_row(e, s_one.get(e["material_id"]), s_count[e["material_id"]],
                     e_count[e["material_id"]]) for e in elastic_ns]

    print("Round-tripping structures through pymatgen...")
    struct_fail = []
    for d in structures:
        try:
            structure_check(d["structure"])
        except Exception as e:
            struct_fail.append({"material_id": d["material_id"],
                                "error": f"{type(e).__name__}: {e}"[:200]})

    out.mkdir(parents=True)
    write_jsonl_gz(out / "elasticity.jsonl.gz", elastic_ns)
    write_jsonl_gz(out / "structures.jsonl.gz", structures)
    write_jsonl_gz(out / "summary.jsonl.gz", summary)
    import csv
    rows_sorted = sorted(rows, key=lambda r: r["material_id"])
    with open(out / "index.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_sorted[0].keys()))
        w.writeheader()
        w.writerows(rows_sorted)

    returned = lambda docs: sorted(set().union(*(d.keys() for d in docs))) if docs else []
    manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "mp_db_version": db_start,
        "endpoint": "materials.elasticity + materials.summary",
        "versions": {"python": platform.python_version(),
                     "mp-api": pkg_version("mp-api"),
                     "emmet-core": pkg_version("emmet-core"),
                     "pymatgen": pkg_version("pymatgen")},
        "script_sha256": sha256(__file__),
        "fields_requested": {"elasticity": ELASTICITY_FIELDS, "summary": SUMMARY_FIELDS},
        "fields_returned": {"elasticity": returned(elastic), "summary": returned(summary)},
        "counts": {"elasticity_docs": len(elastic), "unique_material_ids": len(ids),
                   "summary_docs": len(summary)},
        "id_integrity": {"elasticity": el_audit, "summary": sm_audit},
        "structure_roundtrip_failures": struct_fail,
        "filters_applied": "none (raw snapshot)",
        "files": {},
        "elapsed_s": round(time.time() - t0, 1),
    }
    for name in ("elasticity.jsonl.gz", "structures.jsonl.gz",
                 "summary.jsonl.gz", "index.csv"):
        p = out / name
        manifest["files"][name] = {"sha256": sha256(p), "bytes": p.stat().st_size}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))

    for f in out.iterdir():                            # make snapshot read-only
        f.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    for f in partial.iterdir():
        f.unlink()
    partial.rmdir()

    descriptive_audit(rows)
    print(f"  summary missing / unexpected / duplicate ids: "
          f"{len(sm_audit['missing_ids'])} / {len(sm_audit['unexpected_ids'])} / "
          f"{len(sm_audit['duplicate_ids'])}")
    print(f"  duplicate elasticity ids        : {len(el_audit['duplicate_ids'])}")
    print(f"  structure round-trip failures   : {len(struct_fail)}")
    print(f"\nSnapshot written to {out}  ({manifest['elapsed_s']} s)")
    print("Record in the pre-registration: mp_db_version and the four file SHA-256s "
          "from manifest.json.")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--inspect", action="store_true",
                    help="fetch 3 docs per endpoint, print fields, write nothing")
    ap.add_argument("--out-dir", default="data/raw")
    a = ap.parse_args()
    key = _api_key()
    if a.inspect:
        inspect_mode(key)
    else:
        full_pull(key, a.out_dir)


if __name__ == "__main__":
    main()
