"""Polya fitting and scan-file expansion (pure Python; no solvers needed)."""
import json
import math
import os
import random
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from simgem import analysis as A     # noqa: E402
from simgem import sweep as S        # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "configs", "baseline.json")


class TestSpecialFunctions(unittest.TestCase):
    def test_gammainc(self):
        # P(1, x) = 1 - exp(-x); P(a, x) + Q(a, x) = 1
        for x in (0.1, 1.0, 3.0, 20.0):
            self.assertAlmostEqual(A.gammainc_lower(1.0, x), 1 - math.exp(-x), places=12)
        for a, x in ((0.3, 0.2), (2.5, 4.0), (40.0, 35.0)):
            self.assertAlmostEqual(A.gammainc_lower(a, x) + A.gammainc_upper(a, x), 1.0, places=12)
        self.assertAlmostEqual(A.chi2_sf(3.841459, 1), 0.05, places=5)

    def test_polya_normalised(self):
        G0, th = 40.0, 0.7
        s = sum(A.polya_pdf(g + 0.5, G0, th) for g in range(0, 2000))
        self.assertAlmostEqual(s, 1.0, places=3)
        self.assertAlmostEqual(A.polya_cdf(1e6, G0, th), 1.0, places=12)


class TestPolyaFit(unittest.TestCase):
    def _sample(self, G0, theta, n, seed):
        rng = random.Random(seed)
        k = theta + 1
        return [int(round(rng.gammavariate(k, G0 / k))) for _ in range(n)]

    def test_recovers_parameters(self):
        for G0, th, seed in ((500, 0.5, 1), (60, 1.5, 2), (2000, 0.2, 3)):
            r = A.fit_polya(self._sample(G0, th, 3000, seed))
            self.assertTrue(r["ok"], r["warnings"])
            self.assertLess(abs(r["G0"] - G0) / r["G0_err"], 4, (G0, r["G0"], r["G0_err"]))
            self.assertLess(abs(r["theta"] - th) / r["theta_err"], 4, (th, r["theta"], r["theta_err"]))
            self.assertGreater(r["p_value"], 1e-4)

    def test_zeros_and_threshold(self):
        vals = [0] * 700 + self._sample(100, 1.0, 300, 4)
        r = A.fit_polya(vals, threshold=1)
        self.assertAlmostEqual(r["efficiency"], 0.3, delta=0.01)
        self.assertAlmostEqual(r["mean_all"], 0.3 * 100, delta=6)
        self.assertLess(abs(r["G0"] - 100) / r["G0_err"], 4)

    def test_curve_matches_fit(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("numpy not available")
        r = A.fit_polya(self._sample(120, 0.8, 2000, 6))
        e = r["hist"]["edges"]
        x = np.linspace(e[0], e[-1], 20001)
        y, lo, hi = A.polya_curve(r, x, n_samples=200)
        w = e[1] - e[0]
        integral = (y[:-1] + y[1:]).sum() / 2 * (x[1] - x[0]) / w     # counts per bin -> events
        self.assertAlmostEqual(integral, r["events_in_fit"], delta=0.01 * r["events_in_fit"])
        self.assertTrue(np.all(lo <= y + 1e-9) and np.all(hi >= y - 1e-9))

    def test_zero_inflated_gain(self):
        rng = random.Random(8)
        for zf, G0, th in ((0.3, 1000, 0.5), (0.7, 5, 0.2), (0.0, 200, 1.0)):
            k = th + 1
            vals = [0 if rng.random() < zf else int(round(rng.gammavariate(k, G0 / k)))
                    for _ in range(3000)]
            r = A.fit_gain(vals)
            self.assertTrue(r["ok"])
            true_mean = sum(vals) / len(vals)
            # the fitted gain is the model's mean: it must agree with the data mean
            self.assertLess(abs(r["gain"] - true_mean) / r["gain_err"], 2, (zf, r["gain"], true_mean))
            if G0 > 50:
                self.assertAlmostEqual(r["survival"], 1 - zf, delta=0.03)
                self.assertLess(abs(r["G0"] - G0) / r["G0_err"], 4)

    def test_too_few_events(self):
        r = A.fit_polya([0, 0, 3, 5, 2])
        self.assertFalse(r["ok"])
        self.assertTrue(r["warnings"])

    def test_csv_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "gain_events.csv")
            with open(p, "w") as f:
                f.write("worker,event,absolute_gain,effective_gain\n")
                for i, v in enumerate(self._sample(80, 1.0, 400, 5)):
                    f.write(f"0,{i},{v + 3},{v}\n")
            self.assertEqual(A.find_events(d), p)
            r = A.analyze_run(p, plot=False, quiet=True)
            self.assertTrue(r["ok"])
            self.assertTrue(os.path.exists(os.path.join(d, "gain_fit.json")))


class TestScanSpec(unittest.TestCase):
    def spec(self, **kw):
        d = {"base": BASE, "_file": "x.json", "_name": "x"}
        d.update(kw)
        return d

    def test_grid_zip_and_ranges(self):
        sp = self.spec(axes=[{"name": "a", "path": "foils[*].dV_V", "range": {"start": 500, "stop": 560, "step": 20}},
                             {"name": "b", "path": "gas.pressure_bar", "values": [4, 8]}])
        _, pts = S.expand(sp, {})
        self.assertEqual(len(pts), 8)
        self.assertEqual(pts[1]["label"], "a=500_b=8")
        sp["mode"] = "zip"
        with self.assertRaises(S.ScanError):
            S.expand(sp, {})
        sp["axes"][1]["values"] = [1, 2, 3, 4]
        self.assertEqual(len(S.expand(sp, {})[1]), 4)
        self.assertEqual(S._expand_range({"start": 1, "stop": 100, "num": 3, "log": True}), [1, 10, 100.0])
        self.assertEqual(S._expand_range({"start": 0.1, "stop": 0.3, "step": 0.1}), [0.1, 0.2, 0.3])

    def test_linked_axis_derived_and_points(self):
        sp = self.spec(
            axes=[{"name": "p", "path": "gas.pressure_bar", "values": [4, 8]},
                  {"name": "s", "set": {"foils[0].dV_V": [500, 520], "foils[2].dV_V": [540, 520]}}],
            derived={"gaps[1].field_kV_cm": "base('gaps[1].field_kV_cm') * p / 8"},
            points=[{"label": "ref", "set": {"foils[1].dV_V": 555}}])
        with open(BASE) as f:
            base = json.load(f)
        _, pts = S.expand(sp, base)
        self.assertEqual(len(pts), 5)
        ov = dict(pts[0]["overrides"])
        self.assertEqual(ov["foils[0].dV_V"], 500)
        self.assertEqual(ov["foils[2].dV_V"], 540)
        self.assertAlmostEqual(ov["gaps[1].field_kV_cm"], 2.857 / 2)
        self.assertEqual(dict(pts[-1]["overrides"]), {"foils[1].dV_V": 555})   # literal
        self.assertEqual(pts[-1]["label"], "ref")

    def test_resolve_shares_fields_across_gases(self):
        sp = self.spec(axes=[{"name": "mix", "path": "gas.components",
                              "values": [{"ar": 98, "ch4": 2}, {"ar": 90, "co2": 10}]}],
                       set={"gas.penning.gas": None})
        _, pts, errors = S.resolve_points(sp)
        self.assertEqual(errors, [])
        self.assertEqual(pts[0]["field_hash"], pts[1]["field_hash"])
        self.assertEqual(pts[1]["cfg"]["gas"]["components"], {"ar": 90, "co2": 10})

    def test_invalid_points_reported(self):
        sp = self.spec(axes=[{"name": "d", "path": "foil.hole.inner_diameter_um", "values": [50, 90]}])
        _, _, errors = S.resolve_points(sp)
        self.assertEqual(len(errors), 1)
        self.assertIn("d=90", errors[0])

    def test_safe_eval(self):
        self.assertAlmostEqual(S.safe_eval("sqrt(p) * 2 + max(1, q)", {"p": 4, "q": 3}), 7)
        for bad in ("__import__('os')", "p.__class__", "open('x')", "[x for x in p]"):
            with self.assertRaises(S.ScanError, msg=bad):
                S.safe_eval(bad, {"p": 1})

    def test_example_scan_files(self):
        import glob
        files = sorted(glob.glob(os.path.join(ROOT, "scans", "**", "*.json"), recursive=True))
        self.assertGreater(len(files), 10)
        for path in files:
            fn = os.path.relpath(path, ROOT)
            if True:
                self.assertTrue(S.is_scan_file(path), fn)
                _, pts, errors = S.resolve_points(S.load_spec(path))
                self.assertEqual(errors, [], fn)
                self.assertGreater(len(pts), 0)
        self.assertFalse(S.is_scan_file(BASE))

    def test_parse_only(self):
        self.assertEqual(S.parse_only("0,3-5", 10), [0, 3, 4, 5])
        with self.assertRaises(S.ScanError):
            S.parse_only("12", 10)


if __name__ == "__main__":
    unittest.main()
