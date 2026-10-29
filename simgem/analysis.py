"""
Gain extraction: Polya fits of the per-event gain distribution.

The Polya distribution of the avalanche size G,

    f(G; Gbar, theta) = (theta+1)^(theta+1) / Gamma(theta+1) * (G/Gbar)^theta
                        * exp(-(theta+1) G/Gbar) / Gbar,

is a Gamma distribution with shape k = theta+1 and scale Gbar/(theta+1), so
its CDF is the regularised incomplete gamma function. That gives exact bin
probabilities for integer-valued data (bins [n-1/2, n+1/2)), and the fit is a
binned maximum-likelihood fit conditional on the fit range: events below the
threshold (e.g. G = 0: primary lost, nothing reached the readout) and in the
far tail are excluded and the model is renormalised to the range. No scipy
needed; matplotlib only for the optional plot.

Reported numbers:
  mean_all       mean over ALL events (zeros included) = what a Fe-55 peak measures
  G0, theta      Polya parameters of the distribution above threshold
  efficiency     fraction of events above threshold
"""
import csv
import json
import math
import os

# ----------------------------------------------------------------- special functions
def gammainc_lower(a, x):
    """Regularised lower incomplete gamma P(a, x)."""
    if x <= 0:
        return 0.0
    if x < a + 1:                       # series
        s = term = 1.0 / a
        n = a
        for _ in range(1000):
            n += 1
            term *= x / n
            s += term
            if abs(term) < abs(s) * 1e-15:
                break
        return min(1.0, s * math.exp(-x + a * math.log(x) - math.lgamma(a)))
    return 1.0 - gammainc_upper(a, x)


def gammainc_upper(a, x):
    """Regularised upper incomplete gamma Q(a, x)."""
    if x <= 0:
        return 1.0
    if x < a + 1:
        return 1.0 - gammainc_lower(a, x)
    b, c, d = x + 1 - a, 1e300, 1.0 / (x + 1 - a)   # Lentz continued fraction
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-15:
            break
    return max(0.0, math.exp(-x + a * math.log(x) - math.lgamma(a)) * h)


def chi2_sf(x, k):
    """P(chi^2_k > x)."""
    return gammainc_upper(k / 2.0, x / 2.0) if x > 0 else 1.0


def polya_pdf(g, G0, theta):
    k = theta + 1.0
    if g <= 0:
        return 0.0 if k > 1 else (math.inf if k < 1 else k / G0)
    return math.exp(k * math.log(k) - math.lgamma(k) + (k - 1) * math.log(g / G0)
                    - k * g / G0) / G0


def polya_cdf(g, G0, theta):
    k = theta + 1.0
    return gammainc_lower(k, g * k / G0) if g > 0 else 0.0


# ----------------------------------------------------------------- minimiser
def nelder_mead(f, x0, step, tol=1e-10, max_iter=4000):
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += step[i]
        pts.append(p)
    vals = [f(p) for p in pts]
    for _ in range(max_iter):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol * (abs(vals[0]) + 1e-12):
            break
        cen = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        xr = [cen[j] + (cen[j] - pts[-1][j]) for j in range(n)]
        fr = f(xr)
        if fr < vals[0]:
            xe = [cen[j] + 2 * (cen[j] - pts[-1][j]) for j in range(n)]
            fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = [cen[j] + 0.5 * (pts[-1][j] - cen[j]) for j in range(n)]
            fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = [pts[0][j] + 0.5 * (pts[i][j] - pts[0][j]) for j in range(n)]
                    vals[i] = f(pts[i])
    return pts[0], vals[0]


def _hessian(f, x, h):
    n = len(x)
    H = [[0.0] * n for _ in range(n)]
    f0 = f(x)
    for i in range(n):
        for j in range(i, n):
            if i == j:
                xp, xm = list(x), list(x)
                xp[i] += h[i]
                xm[i] -= h[i]
                H[i][i] = (f(xp) - 2 * f0 + f(xm)) / h[i] ** 2
            else:
                pp, pm, mp, mm = list(x), list(x), list(x), list(x)
                pp[i] += h[i]; pp[j] += h[j]
                pm[i] += h[i]; pm[j] -= h[j]
                mp[i] -= h[i]; mp[j] += h[j]
                mm[i] -= h[i]; mm[j] -= h[j]
                H[i][j] = H[j][i] = (f(pp) - f(pm) - f(mp) + f(mm)) / (4 * h[i] * h[j])
    return H


def _inv2(H):
    a, b, c, d = H[0][0], H[0][1], H[1][0], H[1][1]
    det = a * d - b * c
    if det <= 0 or a <= 0:
        return None
    return [[d / det, -b / det], [-c / det, a / det]]


# ----------------------------------------------------------------- data
def read_column(path, column="effective_gain"):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if rows and column not in rows[0]:
        raise KeyError(f"{path} has no column '{column}' (has {list(rows[0])})")
    return [int(float(r[column])) for r in rows]


def find_events(target):
    """A gain_events.csv from: the file itself, a run folder, or a config
    (-> runs/<config name>/)."""
    if os.path.isfile(target) and target.endswith(".csv"):
        return target
    if os.path.isdir(target):
        p = os.path.join(target, "gain_events.csv")
        if os.path.exists(p):
            return p
    if target.endswith(".json"):
        from .config import run_name
        p = os.path.join(os.environ.get("SIMGEM_RUNS", "runs"), run_name(target), "gain_events.csv")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"no gain_events.csv found for '{target}'")


def _quantile(sorted_vals, q):
    if not sorted_vals:
        return 0
    i = min(len(sorted_vals) - 1, max(0, int(math.ceil(q * len(sorted_vals))) - 1))
    return sorted_vals[i]


# ----------------------------------------------------------------- the fit
def fit_polya(values, threshold=1, max_value=None, bin_width=None, min_events=20):
    """Binned ML Polya fit of integer gains in [threshold, max_value].
    Returns a dict (see module docstring); 'ok' False if the fit failed."""
    vals = sorted(int(v) for v in values)
    n_all = len(vals)
    res = {"ok": False, "events": n_all, "threshold": threshold, "warnings": []}
    if n_all == 0:
        res["warnings"].append("no events")
        return res
    mean_all = sum(vals) / n_all
    var_all = sum((v - mean_all) ** 2 for v in vals) / max(n_all - 1, 1)
    res.update({"mean_all": mean_all, "mean_all_err": math.sqrt(var_all / n_all),
                "rms_all": math.sqrt(var_all)})

    above = [v for v in vals if v >= threshold]
    res["efficiency"] = len(above) / n_all
    res["efficiency_err"] = math.sqrt(res["efficiency"] * (1 - res["efficiency"]) / n_all)
    if not above:
        res["warnings"].append("no events above threshold")
        return res
    hi = max_value if max_value is not None else max(_quantile(above, 0.999), threshold)
    fit_vals = [v for v in above if v <= hi]
    nfit = len(fit_vals)
    m = sum(fit_vals) / nfit
    v2 = sum((x - m) ** 2 for x in fit_vals) / max(nfit - 1, 1)
    res.update({"mean_above": m, "mean_above_err": math.sqrt(v2 / nfit), "fit_range": [threshold, hi],
                "events_in_fit": nfit})
    if nfit < min_events:
        res["warnings"].append(f"only {nfit} events in the fit range; no fit")
        return res
    if m < 5:
        res["warnings"].append(f"mean gain above threshold is only {m:.2g}: a continuous Polya "
                               "is a rough model for such small integer counts")

    # binning: unit bins for small gains, otherwise ~60 bins
    if bin_width is None:
        bin_width = max(1, int(math.ceil((hi - threshold + 1) / 60)))
    edges = [threshold - 0.5]
    while edges[-1] < hi + 0.5 - 1e-9:
        edges.append(edges[-1] + bin_width)
    nb = len(edges) - 1
    counts = [0] * nb
    for v in fit_vals:
        counts[min(int((v - edges[0]) // bin_width), nb - 1)] += 1

    lo_edge = max(edges[0], 0.0)

    def probs(G0, theta):
        cdf = [polya_cdf(max(e, 0.0), G0, theta) for e in edges]
        tot = cdf[-1] - polya_cdf(lo_edge, G0, theta)
        return [(cdf[i + 1] - cdf[i]) / tot for i in range(nb)], tot

    def nll_phys(G0, theta):
        if G0 <= 0 or theta <= -1:
            return 1e300
        p, tot = probs(G0, theta)
        if tot <= 0:
            return 1e300
        s = 0.0
        for n, pi in zip(counts, p):
            if n:
                if pi <= 0:
                    return 1e300
                s -= n * math.log(pi)
        return s

    # start from the moments of the (untruncated-ish) data
    th0 = max(min(m * m / v2 - 1 if v2 > 0 else 1.0, 50.0), -0.5)
    f = lambda u: nll_phys(math.exp(u[0]), math.exp(u[1]) - 1)
    u, _ = nelder_mead(f, [math.log(m), math.log(th0 + 1)], [0.2, 0.3])
    u, nll = nelder_mead(f, u, [0.05, 0.1])
    G0, theta = math.exp(u[0]), math.exp(u[1]) - 1

    cov = _inv2(_hessian(lambda x: nll_phys(x[0], x[1]), [G0, theta],
                         [max(1e-4 * G0, 1e-6), max(1e-3 * (theta + 1), 1e-5)]))
    p, _ = probs(G0, theta)
    mu = [nfit * pi for pi in p]
    # Baker-Cousins likelihood-ratio chi^2 (valid for small counts)
    chi2 = 2 * sum(mi - n + (n * math.log(n / mi) if n > 0 else 0.0) for n, mi in zip(counts, mu) if mi > 0)
    ndof = sum(1 for mi in mu if mi > 0) - 1 - 2
    res.update({
        "ok": cov is not None,
        "G0": G0, "theta": theta,
        "G0_err": math.sqrt(cov[0][0]) if cov else float("nan"),
        "theta_err": math.sqrt(cov[1][1]) if cov else float("nan"),
        "corr": (cov[0][1] / math.sqrt(cov[0][0] * cov[1][1])) if cov else float("nan"),
        "chi2": chi2, "ndof": ndof, "p_value": chi2_sf(chi2, ndof) if ndof > 0 else float("nan"),
        "relative_width": 1 / math.sqrt(theta + 1),
        "hist": {"edges": edges, "counts": counts, "expected": mu},
    })
    if cov is None:
        res["warnings"].append("fit did not converge to a proper minimum (no covariance)")
    return res


def fit_gain(values, max_value=None, bin_width=None, min_events=20):
    """Zero-inflated Polya fit of ALL events (no threshold) -> the gain.

    Model for the integer gain n of one primary electron:
        P(n = 0)  = (1 - eps) + eps * F(1/2)
        P(n >= 1) = eps * [F(n + 1/2) - F(n - 1/2)]
    with F the Polya CDF: a fraction 1 - eps of primaries gives nothing at
    the readout (attached, not collected, all electrons lost), the rest an
    avalanche of Polya(G0, theta). The likelihood factorises into
    Binomial(#zeros) x (Polya fit of the n >= 1 events), so (G0, theta) come
    from the n >= 1 fit and q = P(n >= 1) is the observed non-zero fraction.

    The gain is the mean of the model,
        gain = q * E[G | G >= 1/2] ~= eps * G0,
    i.e. the fitted counterpart of the mean over all events (what a current
    or Fe-55 peak measurement gives). G0 alone is NOT the gain.
    """
    r = fit_polya(values, threshold=1, max_value=max_value, bin_width=bin_width,
                  min_events=min_events)
    r["model"] = "zero-inflated Polya (all events, no threshold)"
    r["threshold"] = 0
    if not r.get("ok"):
        return r
    N, q = r["events"], r["efficiency"]
    G0, th, sG, sT, rho = r["G0"], r["theta"], r["G0_err"], r["theta_err"], r["corr"]

    def parts(g0, t):
        k = t + 1.0
        x = 0.5 * k / g0
        F0 = gammainc_lower(k, x)               # Polya mass below 1/2 (counted as 0)
        M0 = gammainc_lower(k + 1.0, x)         # its share of the first moment
        return F0, g0 * (1.0 - M0) / (1.0 - F0)  # P(G < 1/2), E[G | G >= 1/2]

    F0, cmean = parts(G0, th)
    dG, dT = 1e-4 * G0, 1e-4 * (th + 1)

    def grad(i):
        a = parts(G0 + dG, th)[i] - parts(G0 - dG, th)[i]
        b = parts(G0, th + dT)[i] - parts(G0, th - dT)[i]
        return a / (2 * dG), b / (2 * dT)

    def var(gg, gt):
        return gg * gg * sG * sG + gt * gt * sT * sT + 2 * gg * gt * rho * sG * sT

    var_q = q * (1 - q) / N
    gain = q * cmean
    gain_err = math.sqrt(cmean * cmean * var_q + q * q * var(*grad(1)))
    eps = q / (1.0 - F0)
    gF = grad(0)
    eps_err = math.sqrt(var_q / (1 - F0) ** 2 + (q / (1 - F0) ** 2) ** 2 * var(*gF))
    r.update({"gain": gain, "gain_err": gain_err, "survival": eps, "survival_err": eps_err,
              "polya_mass_below_half": F0,
              # resolution: relative variance of the per-primary gain (zeros included)
              "relative_variance_all": (r["rms_all"] / r["mean_all"]) ** 2 if r["mean_all"] > 0 else float("nan")})
    return r


def format_fit(r, title="Polya fit"):
    L = [title]
    zi = "gain" in r
    if zi:
        L.append(f"  GAIN (model mean)    : {r['gain']:.4g} +/- {r['gain_err']:.2g}"
                 f"   [zero-inflated Polya, all events]")
    if "mean_all" in r:
        L.append(f"  mean over all events : {r['mean_all']:.4g} +/- {r['mean_all_err']:.2g}"
                 f"   ({r['events']} events)")
    if zi:
        L.append(f"  survival eps         : {100 * r['survival']:.1f}% +/- {100 * r['survival_err']:.1f}%"
                 f"   (zero-gain events: {100 * (1 - r['efficiency']):.1f}%)")
    elif "efficiency" in r:
        L.append(f"  above threshold {r['threshold']:<4} : {100 * r['efficiency']:.1f}% "
                 f"(+/- {100 * r['efficiency_err']:.1f}%)")
    if r.get("ok"):
        L.append(f"  Polya G0             : {r['G0']:.4g} +/- {r['G0_err']:.2g}"
                 + ("   (avalanche mean; not the gain)" if zi else ""))
        L.append(f"  Polya theta          : {r['theta']:.3g} +/- {r['theta_err']:.2g}"
                 f"   (relative width {r['relative_width']:.2f})")
        L.append(f"  fit range            : {r['fit_range'][0]}..{r['fit_range'][1]}, "
                 f"chi2/ndof = {r['chi2']:.1f}/{r['ndof']}, p = {r['p_value']:.3g}")
    for w in r.get("warnings", []):
        L.append(f"  NOTE: {w}")
    return "\n".join(L)


def polya_curve(r, x, n_samples=2000, cl=0.68, seed=1):
    """Fitted Polya as expected counts per histogram bin at gains `x`, plus a
    `cl` band from sampling (G0, theta) with the fit covariance.
    Returns (y, y_low, y_high) as numpy arrays (numpy needed)."""
    import numpy as np
    x = np.asarray(x, dtype=float)
    e = r["hist"]["edges"]
    lo, hi, w = max(e[0], 0.0), e[-1], e[1] - e[0]
    n = r["events_in_fit"]

    def curve(G0, th):
        k = th + 1.0
        with np.errstate(divide="ignore", invalid="ignore"):
            pdf = np.where(x > 0, np.exp(k * math.log(k) - math.lgamma(k)
                                         + (k - 1) * np.log(np.maximum(x, 1e-300) / G0)
                                         - k * x / G0) / G0, 0.0)
        norm = polya_cdf(hi, G0, th) - polya_cdf(lo, G0, th)
        return n * w * pdf / norm

    y = curve(r["G0"], r["theta"])
    s0, s1, rho = r["G0_err"], r["theta_err"], r["corr"]
    cov = np.array([[s0 * s0, rho * s0 * s1], [rho * s0 * s1, s1 * s1]])
    rng = np.random.default_rng(seed)
    pars = rng.multivariate_normal([r["G0"], r["theta"]], cov, size=n_samples)
    pars = pars[(pars[:, 0] > 0) & (pars[:, 1] > -1)]
    ys = np.array([curve(g, t) for g, t in pars])
    q = 100 * (1 - cl) / 2
    ylo, yhi = np.percentile(ys, [q, 100 - q], axis=0)
    return y, ylo, yhi


def plot_fit(r, path, title=None, log=True):
    """Histogram + fitted Polya (+ pulls). Needs matplotlib; returns False if missing."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    if "hist" not in r:
        return False
    h = r["hist"]
    e, c, mu = h["edges"], h["counts"], h["expected"]
    x = [(e[i] + e[i + 1]) / 2 for i in range(len(c))]
    w = e[1] - e[0]
    fig, (ax, axp) = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    ax.errorbar(x, c, yerr=[math.sqrt(max(v, 1)) for v in c], fmt="k.", ms=4, label="simulation")
    ax.step(e, list(mu) + [mu[-1]], where="post", color="C3", lw=1.2,
            label=f"Polya: G0 = {r['G0']:.3g} ± {r['G0_err']:.2g}, θ = {r['theta']:.2g} ± {r['theta_err']:.2g}")
    if log:
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.5)
    ax.set_ylabel(f"events / {w:g}")
    ax.legend(fontsize=9)
    if title is None:
        title = (f"gain {r['gain']:.3g} ± {r['gain_err']:.2g} (mean of all events {r['mean_all']:.3g}), "
                 f"survival {100 * r['survival']:.0f}%" if "gain" in r else
                 f"mean (all events) {r['mean_all']:.3g}, {100 * r['efficiency']:.0f}% above {r['threshold']}")
    ax.set_title(title + f", χ²/ndf {r['chi2']:.0f}/{r['ndof']}", fontsize=10)
    pulls = [(n - m) / math.sqrt(m) if m > 0 else 0 for n, m in zip(c, mu)]
    axp.bar(x, pulls, width=w, color="0.6")
    axp.axhline(0, color="k", lw=0.5)
    axp.set_ylabel("pull")
    axp.set_xlabel({"effective_gain": "electrons at readout (effective gain)",
                    "absolute_gain": "electrons produced (absolute gain)"}.get(r.get("column"), r.get("column", "")))
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return True


def analyze_run(events_csv, outdir=None, column="effective_gain", threshold=None,
                max_value=None, bin_width=None, plot=True, quiet=False):
    """Fit one run's gain distribution; writes <outdir>/gain_fit.json (+ .png).
    threshold None: zero-inflated Polya of all events (the gain). A threshold
    (e.g. an electronics threshold in electrons) gives a plain Polya fit of
    the events above it, for comparison with single-electron spectra."""
    outdir = outdir or os.path.dirname(events_csv)
    vals = read_column(events_csv, column)
    if threshold is None:
        r = fit_gain(vals, max_value=max_value, bin_width=bin_width)
    else:
        r = fit_polya(vals, threshold=threshold, max_value=max_value, bin_width=bin_width)
    r["column"] = column
    r["source"] = os.path.abspath(events_csv)
    stem = "gain_fit" if column == "effective_gain" else f"gain_fit_{column}"
    with open(os.path.join(outdir, stem + ".json"), "w") as f:
        json.dump(r, f, indent=2)
    if plot and r.get("ok"):
        r["plot"] = os.path.join(outdir, stem + ".png") if plot_fit(r, os.path.join(outdir, stem + ".png")) else None
    if not quiet:
        print(format_fit(r, f"Polya fit of {column} ({events_csv})"))
        if r.get("plot"):
            print(f"  plot -> {r['plot']}")
    return r
