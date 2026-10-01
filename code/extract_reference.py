"""
extract_reference.py — Paper 2 elemental reference table (P_phys input)
=======================================================================
Extracts solid density, shear modulus G, Young's modulus E, Poisson number
and transverse sound velocity for every element from the "Mechanical
properties" tables (2.1-6B(b) ... 2.1-26B(b)) of

    W. Martienssen & H. Warlimont (eds.), Springer Handbook of Condensed
    Matter and Materials Data, Springer, 2005.

and applies the pre-registered internal-consistency screen (rule "C").

The handbook states (Sect. 2.1.3, printed p. 49) that most of its element
data are taken from Landolt-Boernstein, with additions from D'Ans-Lax and
the CRC Handbook. Values are therefore COMPILED data. The screen tests
internal consistency only; it is not evidence that a value is correct.

Pre-registered rule (frozen)
----------------------------
Poisson check, computable iff E, G and printed nu are all tabulated:
    nu_implied = E/(2G) - 1
    pass  iff  0 <= nu_implied <= 0.5  and  |nu_implied - nu_printed| <= 0.10
Sound check, computable iff density, v_t and G are all tabulated:
    pass  iff  |rho * v_t^2 - G| / G <= 0.15
Eligibility:
    - density or G not tabulated                  -> unavailable
    - G and density tabulated, no check computable -> unavailable
    - otherwise eligible iff AT LEAST ONE computable check passes
      (when only one check is computable, it alone decides).
No value is repaired, imputed or reconstructed (no rho*v_t^2, no E/2(1+nu),
no VRH from single-crystal constants).

Extraction notes
----------------
- Columns are keyed by ELEMENT NAME, never by the printed symbol:
  Table 2.1-26B(b) prints curium's symbol as "Cu". Asserted below.
- Footnote markers attached to values (e.g. Ca "7.85 a") are stripped from
  the number and kept; the footnote text is recorded.
- For an element appearing in a table and its "cont." table, the first
  tabulated value of each quantity (document order) is used.
- PDF page = printed page + 10 in this file.

Usage
-----
    python extract_reference.py HANDBOOK.pdf  [--out-dir DIR]

Outputs: elemental_reference_springer2005_v1.csv and
         elemental_reference_springer2005_v1.manifest.json
"""

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import subprocess
from pathlib import Path

VERSION = "v1"
PDF_FIRST, PDF_LAST, PAGE_OFFSET = 55, 170, 10
NU_TOL, SOUND_TOL = 0.10, 0.15

# name (as printed) -> symbol. Built in, to avoid a pymatgen dependency.
NAMES = {
 "Hydrogen": "H", "Helium": "He", "Lithium": "Li", "Beryllium": "Be", "Boron": "B",
 "Carbon": "C", "Nitrogen": "N", "Oxygen": "O", "Fluorine": "F", "Neon": "Ne",
 "Sodium": "Na", "Magnesium": "Mg", "Aluminium": "Al", "Aluminum": "Al",
 "Silicon": "Si", "Phosphorus": "P", "Phosphorous": "P", "Sulfur": "S",
 "Chlorine": "Cl", "Argon": "Ar", "Potassium": "K", "Calcium": "Ca",
 "Scandium": "Sc", "Titanium": "Ti", "Vanadium": "V", "Chromium": "Cr",
 "Manganese": "Mn", "Iron": "Fe", "Cobalt": "Co", "Nickel": "Ni", "Copper": "Cu",
 "Zinc": "Zn", "Gallium": "Ga", "Germanium": "Ge", "Arsenic": "As",
 "Selenium": "Se", "Bromine": "Br", "Krypton": "Kr", "Rubidium": "Rb",
 "Strontium": "Sr", "Yttrium": "Y", "Zirconium": "Zr", "Niobium": "Nb",
 "Molybdenum": "Mo", "Technetium": "Tc", "Ruthenium": "Ru", "Rhodium": "Rh",
 "Palladium": "Pd", "Silver": "Ag", "Cadmium": "Cd", "Indium": "In", "Tin": "Sn",
 "Antimony": "Sb", "Tellurium": "Te", "Iodine": "I", "Xenon": "Xe",
 "Cesium": "Cs", "Caesium": "Cs", "Barium": "Ba", "Lanthanum": "La",
 "Cerium": "Ce", "Praseodymium": "Pr", "Neodymium": "Nd", "Promethium": "Pm",
 "Samarium": "Sm", "Europium": "Eu", "Gadolinium": "Gd", "Terbium": "Tb",
 "Dysprosium": "Dy", "Holmium": "Ho", "Erbium": "Er", "Thulium": "Tm",
 "Ytterbium": "Yb", "Lutetium": "Lu", "Hafnium": "Hf", "Tantalum": "Ta",
 "Tungsten": "W", "Rhenium": "Re", "Osmium": "Os", "Iridium": "Ir",
 "Platinum": "Pt", "Gold": "Au", "Mercury": "Hg", "Thallium": "Tl", "Lead": "Pb",
 "Bismuth": "Bi", "Polonium": "Po", "Astatine": "At", "Radon": "Rn",
 "Francium": "Fr", "Radium": "Ra", "Actinium": "Ac", "Thorium": "Th",
 "Protactinium": "Pa", "Uranium": "U", "Neptunium": "Np", "Plutonium": "Pu",
 "Americium": "Am", "Curium": "Cm", "Berkelium": "Bk", "Californium": "Cf",
 "Einsteinium": "Es", "Fermium": "Fm", "Mendelevium": "Md", "Nobelium": "No",
 "Lawrencium": "Lr", "Darmstadtium": "Ds", "Roentgenium": "Rg", "Ununnilium": "Ds", "Rutherfordium": "Rf", "Dubnium": "Db",
 "Seaborgium": "Sg", "Bohrium": "Bh", "Hassium": "Hs", "Meitnerium": "Mt",
}

# Row label (prefix match on the joined, whitespace-normalised label) -> field
ROWS = [
    ("Density ρ, solid", "density_g_cm3"), ("Density , solid", "density_g_cm3"),
    ("Shear modulus G", "G_GPa"),
    ("Elastic modulus E", "E_GPa"),
    ("Poisson number", "nu_printed"),
    ("Sound velocity, solid, transverse", "v_t_m_s"),
    ("Modification", "modification"), ("Modiﬁcation", "modification"),
    ("Modifaction", "modification"),
]
UNIT_TOKENS = {"g/cm3", "GPa", "m/s", "1/K", "1/MPa", "1/TPa", "MPa", "N/m",
               "cm3 /mol", "mPa s", "pm"}

# Manual annotations from visual inspection of the rendered pages. These are
# descriptive metadata only; none of them changes eligibility.
PROVENANCE = {
    "Ca": "G footnote a: 'At 348 K' (not room temperature).",
    "Ce": "Handbook footnote: elastic data concern gamma-Ce.",
    "Tc": "Printed v_t = 50.6 m/s is a misprint (Tc is a stiff hcp metal); "
          "printed E (407 GPa) is identical to W's E.",
    "Ta": "Printed v_t inconsistent with G (rho*v_t^2 ~ 140 GPa); G agrees with E and nu.",
    "Tb": "Printed v_t inconsistent with G; G agrees with E and nu.",
    "Sm": "Printed E inconsistent with G and nu; G agrees with rho*v_t^2.",
    "Re": "Printed E inconsistent with nu; G agrees with rho*v_t^2.",
    "Pm": "No sound velocities or elastic constants tabulated; round values; "
          "no stable isotope. Provenance of G unclear.",
    "Hg": "Liquid at room temperature; no solid density or G tabulated.",
    "Cm": "Table 2.1-26B(b) prints the symbol as 'Cu' (misprint).",
}
VT_FLAG = {"Tc": "misprint", "Ta": "inconsistent_with_G", "Tb": "inconsistent_with_G"}
# Pre-registered sensitivity analysis (not the primary rule).
SENSITIVITY_EXCLUDE = {"Pm", "Tc"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def tokens(line):
    """(start_column, text) for runs separated by >=2 spaces."""
    return [(m.start(), m.group()) for m in re.finditer(r"\S+(?: \S+)*", line)]


def parse_number(s):
    """'7.85 a' -> (7.85, 'a'); '1.18 (280 K)' -> (1.18, '(280 K)'); '' -> (None, '')."""
    if not s:
        return None, ""
    s = s.replace(",", ".").replace("−", "-").strip()
    m = re.match(r"^(-?\d+(?:\.\d+)?)(.*)$", s.replace(" ", "", 1) if re.match(r"^\d+ \d{3}\b", s) else s)
    if not m:
        return None, s
    rest = m.group(2).strip()
    if re.match(r"^[×x]\s*10", rest):          # scientific notation: not a plain number
        return None, s
    return float(m.group(1)), rest


def extract(pdf):
    txt = subprocess.run(["pdftotext", "-f", str(PDF_FIRST), "-l", str(PDF_LAST),
                          "-layout", pdf, "-"], capture_output=True, text=True,
                         check=True).stdout
    L = txt.split("\n")
    page = []
    p = PDF_FIRST
    for line in L:
        if "\f" in line:
            p += line.count("\f")
        page.append(p)

    out = {}
    starts = [i for i, l in enumerate(L)
              if "Table 2.1-" in l and "Mechanical properties" in l]
    for s in starts:
        title = re.sub(r"\s+", " ", L[s].replace("\f", "")).strip()
        table_id = re.search(r"Table (2\.1-\S+?\(b\))", title).group(1)
        hdr = next(i for i in range(s, s + 15) if L[i].lstrip().startswith("Element name"))
        cols = []
        for pos, t in tokens(L[hdr]):
            if t == "Element name":
                continue
            name = re.sub(r"\s+[a-e]$", "", t)
            if name not in NAMES:
                raise SystemExit(f"Unknown element name {t!r} in {table_id}")
            cols.append((pos, NAMES[name]))
        first_col = cols[0][0]
        end = next((i for i in range(hdr + 1, len(L))
                    if "\f" in L[i] or "Table 2.1-" in L[i]), len(L))

        # logical rows: a line with a label and no values continues the previous label
        rows = []
        for i in range(hdr + 1, end):
            line = L[i]
            label = re.sub(r"\s+", " ", line[:max(first_col - 2, 0)]).strip()
            vals = [(pos, t) for pos, t in tokens(line) if pos >= first_col - 2]
            if not label and not vals:
                continue
            if label and not vals and rows and not rows[-1]["closed"]:
                rows[-1]["label"] += " " + label
                rows[-1]["closed"] = True
                continue
            rows.append({"label": label, "vals": vals, "line": i, "closed": False})

        footnotes = {}
        for i in range(hdr + 1, min(end + 4, len(L))):
            m = re.match(r"^\s*([a-e])\s+([A-Z].*\S)\s*$", L[i])
            if m and len(m.group(2)) < 120:
                footnotes.setdefault(m.group(1), m.group(2))

        for r in rows:
            field = next((f for pre, f in ROWS if r["label"].startswith(pre)), None)
            if field is None:
                continue
            for pos, t in r["vals"]:
                if t in UNIT_TOKENS or pos > cols[-1][0] + 14:
                    continue
                sym = min(cols, key=lambda c: abs(c[0] - pos))[1]
                d = out.setdefault(sym, {})
                if field in d:                      # first value wins
                    continue
                d[field] = t
                d[field + "_src"] = f"{table_id} p.{page[r['line']] - PAGE_OFFSET}"
                d.setdefault("pdf_page", page[r["line"]])
                d.setdefault("printed_page", page[r["line"]] - PAGE_OFFSET)
                d.setdefault("table", table_id)
                _, mark = parse_number(t)
                if field == "G_GPa" and mark in footnotes:
                    d["G_footnote"] = f"{mark}: {footnotes[mark]}"
        for pos, sym in cols:
            out.setdefault(sym, {}).setdefault("table", table_id)
            out[sym].setdefault("printed_page", page[hdr] - PAGE_OFFSET)
            out[sym].setdefault("pdf_page", page[hdr])
    return out


def screen(sym, d):
    rho, _ = parse_number(d.get("density_g_cm3"))
    G, _ = parse_number(d.get("G_GPa"))
    E, _ = parse_number(d.get("E_GPa"))
    nu, _ = parse_number(d.get("nu_printed"))
    vt, _ = parse_number(d.get("v_t_m_s"))

    nu_imp = E / (2 * G) - 1 if (E and G) else None
    pois = None
    if E and G and nu is not None:
        pois = bool(0 <= nu_imp <= 0.5 and abs(nu_imp - nu) <= NU_TOL)
    g_snd = rho * 1e3 * vt ** 2 / 1e9 if (rho and vt) else None
    snd = None
    if rho and vt and G:
        snd = bool(abs(g_snd - G) / G <= SOUND_TOL)

    if rho is None or G is None:
        status, eligible = "unavailable", False
    elif pois is None and snd is None:
        status, eligible = "unavailable_no_check_computable", False
    elif pois or snd:
        status, eligible = "passes", True
    elif pois is False and snd is False:
        status, eligible = "fails_both_checks", False
    else:
        status, eligible = "fails_only_computable_check", False

    temp_K = None
    m = re.search(r"At (\d+(?:\.\d+)?) K", d.get("G_footnote", ""))
    if m:
        temp_K = float(m.group(1))
    return {
        "element": sym,
        "table": d.get("table", ""),
        "printed_page": d.get("printed_page", ""),
        "pdf_page": d.get("pdf_page", ""),
        "modification": d.get("modification", ""),
        "density_g_cm3": rho if rho is not None else "",
        "G_GPa": G if G is not None else "",
        "E_GPa": E if E is not None else "",
        "nu_printed": nu if nu is not None else "",
        "v_t_m_s": vt if vt is not None else "",
        "density_raw": d.get("density_g_cm3", ""),
        "G_raw": d.get("G_GPa", ""),
        "G_footnote": d.get("G_footnote", ""),
        "nu_implied": round(nu_imp, 4) if nu_imp is not None else "",
        "G_rho_vt2_GPa": round(g_snd, 3) if g_snd is not None else "",
        "check_poisson": "" if pois is None else pois,
        "check_sound": "" if snd is None else snd,
        "status": status,
        "eligible": eligible,
        "eligible_sens_no_Pm_Tc": eligible and sym not in SENSITIVITY_EXCLUDE,
        "temperature_K": temp_K if temp_K is not None else "",
        "temperature_status": "stated_in_footnote" if temp_K else "not_stated",
        "vt_flag": VT_FLAG.get(sym, ""),
        "provenance_note": PROVENANCE.get(sym, ""),
    }


def assertions(rows):
    r = {x["element"]: x for x in rows}
    checks = {
        "Cu density is copper's (8.96), not curium's (13.51)": r["Cu"]["density_g_cm3"] == 8.96,
        "Cm density 13.51": r["Cm"]["density_g_cm3"] == 13.51,
        "Al G = 27.8": r["Al"]["G_GPa"] == 27.8,
        "Ca G = 7.85 with 348 K footnote": r["Ca"]["G_GPa"] == 7.85 and r["Ca"]["temperature_K"] == 348,
        "Nb G = 59.5 fails both checks": r["Nb"]["status"] == "fails_both_checks",
        "Hg unavailable": r["Hg"]["status"].startswith("unavailable"),
        "Pu modification alpha-Pu": "Pu" in r["Pu"]["modification"],
        "W G = 152": r["W"]["G_GPa"] == 152,
    }
    bad = [k for k, ok in checks.items() if not ok]
    if bad:
        raise SystemExit("ASSERTION FAILED: " + "; ".join(bad))
    return list(checks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf")
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()
    raw = extract(a.pdf)
    rows = sorted((screen(s, d) for s, d in raw.items()), key=lambda x: x["element"])
    passed = assertions(rows)

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / f"elemental_reference_springer2005_{VERSION}.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "source": "Martienssen & Warlimont (eds.), Springer Handbook of Condensed "
                  "Matter and Materials Data, Springer, 2005; tables 2.1-6B(b)..2.1-26B(b)",
        "source_pdf_sha256": sha256(a.pdf),
        "script_sha256": sha256(__file__),
        "pdftotext": subprocess.run(["pdftotext", "-v"], capture_output=True,
                                    text=True).stderr.splitlines()[0],
        "rule": {"nu_tol": NU_TOL, "sound_rel_tol": SOUND_TOL,
                 "nu_implied_range": [0, 0.5], "logic": "eligible iff >=1 computable check passes",
                 "no_check_computable": "unavailable", "repairs": "none"},
        "sensitivity_exclude": sorted(SENSITIVITY_EXCLUDE),
        "assertions_passed": passed,
        "n_elements": len(rows),
        "n_eligible": sum(r["eligible"] for r in rows),
        "csv_sha256": sha256(csv_path),
    }
    (out / f"elemental_reference_springer2005_{VERSION}.manifest.json").write_text(
        json.dumps(manifest, indent=2))
    print(f"{len(rows)} elements, {manifest['n_eligible']} eligible -> {csv_path}")
    print(f"csv sha256: {manifest['csv_sha256']}")


if __name__ == "__main__":
    main()
