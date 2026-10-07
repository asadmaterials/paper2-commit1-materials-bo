import numpy as np
from scipy import stats
rng = np.random.default_rng(11)
def reject(sd_r, sd_c, n=200000, alpha=0.05/12):
    m = np.stack([rng.normal(0, sd_r, (n, 10)), rng.normal(0, sd_c, (n, 10))], 1)   # split means under H0
    v = m.var(2, ddof=1); est = m.mean((1, 2)); se = np.sqrt(v.sum(1) / 40)
    df = v.sum(1) ** 2 / ((v ** 2).sum(1) / 9)
    t = est / se
    return (stats.t.sf(t, 18) <= alpha).mean(), (stats.t.sf(t, df) <= alpha).mean()
print("nominal %.5f" % (0.05 / 12))
for a, b in ((1, 1), (0.2, 2), (0, 1)):
    print("between-split SD", a, b, "-> fixed 18 df: %.5f, Satterthwaite: %.5f" % reject(a if a else 1e-9, b))
# H2: probability that at least one of 6 betas is declared safe when the true margin is exactly 0 for all
n = 40000
out = {"unadj": 0, "adj": 0}
for rho in (0.0, 0.5, 0.9):
    shared = rng.normal(size=(n, 1, 2, 10)); own = rng.normal(size=(n, 6, 2, 10))
    m = np.sqrt(rho) * shared + np.sqrt(1 - rho) * own
    v = m.var(3, ddof=1); est = m.mean((2, 3)); se = np.sqrt(v.sum(2) / 40); df = v.sum(2) ** 2 / ((v ** 2).sum(2) / 9)
    un = (est + stats.t.ppf(0.95, df) * se < 0).any(1).mean(); ad = (est + stats.t.ppf(1 - 0.05 / 6, df) * se < 0).any(1).mean()
    print("correlation between betas %.1f: P(any safe) unadjusted %.3f, adjusted %.3f" % (rho, un, ad))
