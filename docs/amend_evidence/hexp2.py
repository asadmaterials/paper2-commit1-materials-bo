import json, sys
import numpy as np
from multiprocessing import Pool
import hexp
from hexp import pool, job, SPLITS
if __name__ == "__main__":
    jobs = []
    for sp in SPLITS:
        ids, Z, prow, P = pool(sp)
        G = np.array([float(prow[m]["G_vrh_GPa"]) for m in ids])
        D = np.sqrt(np.maximum(0, (Z ** 2).sum(1)[:, None] + (Z ** 2).sum(1)[None, :] - 2 * Z @ Z.T)); iu = np.triu_indices(len(ids), 1)
        cand = sorted([(abs(G[i] - G[j]), int(i), int(j)) for i, j in zip(*iu) if D[i, j] < 0.05], reverse=True)
        for _, i, j in cand[:2]:
            for lb in (1.0, 2.0):
                for beta in (0.0, 0.05):
                    jobs.append(("E2", sp, "Y1", 1, lb, beta, (i, j)))
            for lb in (0.01, 1.0, 2.0, 3.0):          # same block without the forced pair, for reference
                jobs.append(("E2ref", sp, "Y1", 1, lb, 0.0, None))
    for sp in SPLITS:
        for seed in range(1, 6):
            for obj in ("Y1", "Y2"):
                jobs.append(("E1", sp, obj, seed, 2.0, 0.0, None))
    print(len(jobs), flush=True)
    with Pool(2) as p, open("hexp2.jsonl", "w") as f:
        for r in p.imap_unordered(job, jobs, chunksize=1):
            f.write(json.dumps(r) + "\n"); f.flush()
    print("DONE", flush=True)
