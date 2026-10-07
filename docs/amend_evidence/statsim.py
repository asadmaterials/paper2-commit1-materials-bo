import numpy as np
from scipy import stats
rng = np.random.default_rng(7)
def sim(k, s, icc, mu, nsim=4000, B=1000):
    """2 types x k splits x s seeds; total variance 1; returns reject rates at alpha and alpha/12 for t and bootstrap."""
    ss, se = np.sqrt(icc), np.sqrt(1 - icc)
    y = mu + rng.normal(0, 0.3, (nsim, 2, 1, 1)) * 0 + ss * rng.normal(size=(nsim, 2, k, 1)) + se * rng.normal(size=(nsim, 2, k, s))
    m = y.mean(3)                                   # split means (nsim, 2, k)
    est = m.mean((1, 2))
    var = m.var(2, ddof=1).mean(1) / (2 * k)        # stratified: pooled within-type variance of split means
    t = est / np.sqrt(var); df = 2 * (k - 1)
    p_t = stats.t.sf(t, df)
    # hierarchical bootstrap as registered (splits within type, then seeds)
    p_b = np.empty(nsim)
    sp = rng.integers(0, k, (B, 2, k)); sd = rng.integers(0, s, (B, 2, k, s))
    kk = np.arange(2)[None, :, None, None]
    for i in range(nsim):
        bm = y[i][kk, sp[:, :, :, None], sd].reshape(B, -1).mean(1)
        p_b[i] = (1 + np.sum(bm <= 0)) / (B + 1)
    return [(p_t <= a).mean() for a in (0.05, 0.05 / 12)] + [(p_b <= a).mean() for a in (0.05, 0.05 / 12)]
print("design        ICC  mu   | t: a=.05  a=.05/12 | boot: a=.05  a=.05/12")
for k, s in ((3, 10), (10, 3), (10, 6)):
    for icc in (0.0, 0.2, 0.5):
        for mu in (0.0, 0.25, 0.5):
            r = sim(k, s, icc, mu, nsim=2000 if mu else 4000)
            print(f"{k:2d} splits x{s:2d}  {icc:.1f}  {mu:.2f} |   {r[0]:.3f}    {r[1]:.4f}  |     {r[2]:.3f}     {r[3]:.4f}", flush=True)
print("critical t one-sided at .05/12: df4 %.2f  df18 %.2f"%(stats.t.isf(0.05/12,4), stats.t.isf(0.05/12,18)))
