"""Fast checks of the pure-Python parts (no gmsh/Elmer/Garfield needed).

    python -m unittest discover tests
"""
import copy
import json
import math
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from simgem import config as C          # noqa: E402
from simgem import geometry as G        # noqa: E402
from simgem import misalignment, runner  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def raw(name="baseline"):
    with open(os.path.join(ROOT, "configs", f"{name}.json")) as f:
        return json.load(f)


class TestConfig(unittest.TestCase):
    def test_all_example_configs_resolve(self):
        for fn in os.listdir(os.path.join(ROOT, "configs")):
            if fn.endswith(".json"):
                C.load(os.path.join(ROOT, "configs", fn))

    def test_potentials_match_legacy_triple_gem(self):
        # values of the original SimGEM baseline (zero-thickness copper)
        st = C.Stack(C.resolve(raw("baseline_square")))
        legacy = {"drift": -3537.77, "gem1_top": -3509.97, "gem1_bot": -2959.97,
                  "gem2_top": -2360.00, "gem2_bot": -1810.00, "gem3_top": -1309.99,
                  "gem3_bot": -759.99, "readout": 0.0}
        for k, v in legacy.items():
            self.assertAlmostEqual(st.V[k], v, places=2)
        self.assertAlmostEqual(st.height, 7.45)

    def test_heights_with_copper(self):
        st = C.Stack(C.resolve(raw()))
        self.assertAlmostEqual(st.z["gem3_bot"], 2.1)
        self.assertAlmostEqual(st.z["gem3_top"], 2.16)
        self.assertAlmostEqual(st.kapton["gem3"][0], 2.105)
        self.assertAlmostEqual(st.kapton["gem3"][1], 2.155)
        self.assertAlmostEqual(st.z["drift"], 3 * 2.1 + 3 * 0.06 + 1.0)

    def test_any_number_of_foils(self):
        r = raw()
        r["foils"] = [{"name": "a", "dV_V": 400}, {"name": "b", "dV_V": 380}]
        r["gaps"] = r["gaps"][:3]
        st = C.Stack(C.resolve(r))
        self.assertEqual(st.electrodes(), ["drift", "a_top", "a_bot", "b_top", "b_bot", "readout"])
        with self.assertRaises(C.ConfigError):
            r["gaps"] = r["gaps"][:2]
            C.resolve(r)

    def test_validation_errors(self):
        bad = [("gas.penning.r", 1.5), ("gas.penning.gas", "xe"), ("foil.pattern", "tri"),
               ("foil.hole.inner_diameter_um", 80.0), ("simulation.start_depth_um", 5000.0)]
        for path, v in bad:
            r = raw()
            C.set_path(r, path, v)
            with self.assertRaises(C.ConfigError, msg=path):
                C.resolve(r)
        r = raw()
        r["gas"]["components"] = {"ar": 90, "co2": 8, "ch4": 2}
        r["gas"]["penning"]["r"] = "auto"
        with self.assertRaises(C.ConfigError):
            C.resolve(r)

    def test_overrides(self):
        r = raw()
        C.set_path(r, "foils[*].dV_V", 600)
        C.set_path(r, "gaps.1.thickness_mm", 1.0)
        C.set_path(r, "gas.penning.r", "auto")
        self.assertTrue(all(f["dV_V"] == 600 for f in r["foils"]))
        self.assertEqual(r["gaps"][1]["thickness_mm"], 1.0)
        self.assertEqual(C.parse_override("a.b=0.3"), ("a.b", 0.3))
        self.assertEqual(C.parse_override('a="x"'), ("a", "x"))
        self.assertEqual(C.parse_override("a=auto"), ("a", "auto"))

    def test_per_foil_override_merges_hole(self):
        r = raw()
        r["foils"][1]["hole"] = {"inner_diameter_um": 55.0}
        cfg = C.resolve(r)
        self.assertEqual(cfg["foils"][1]["hole"]["inner_diameter_um"], 55.0)
        self.assertEqual(cfg["foils"][1]["hole"]["outer_diameter_um"], 70.0)
        self.assertEqual(cfg["foils"][0]["hole"]["inner_diameter_um"], 50.0)

    def test_field_hash(self):
        cfg = C.resolve(raw())
        h = C.Stack(cfg).field_hash()
        c2 = copy.deepcopy(cfg)
        c2["gas"]["penning"]["r"] = 0.5
        c2["simulation"]["events"] = 7
        self.assertEqual(C.Stack(c2).field_hash(), h)        # gas/sim don't touch the field
        c3 = copy.deepcopy(cfg)
        c3["foils"][0]["dV_V"] = 551
        self.assertNotEqual(C.Stack(c3).field_hash(), h)
        self.assertNotEqual(C.Stack(cfg, {"gem2": (0.01, 0.0)}).field_hash(), h)

    def test_run_key(self):
        cfg = C.resolve(raw())
        k = runner.run_key(cfg, "abc", 1)
        c2 = copy.deepcopy(cfg)
        c2["gas"]["penning"]["r"] = 0.3
        self.assertNotEqual(runner.run_key(c2, "abc", 1), k)
        c2 = copy.deepcopy(cfg)
        c2["simulation"]["events"] = 5
        self.assertEqual(runner.run_key(c2, "abc", 1), k)


class TestPenningMatch(unittest.TestCase):
    def setUp(self):
        from simgem import penning
        self.P = penning
        self.r = raw("gorg")
        self.r["gas"]["penning"]["match"] = "table"

    def test_validation(self):
        C.resolve(self.r)
        for path, v in [("gas.penning.match", "micro"), ("gas.penning.r", None),
                        ("gas.penning.match_fields_kV_cm", [40, -1])]:
            r = copy.deepcopy(self.r)
            C.set_path(r, path, v)
            with self.assertRaises(C.ConfigError, msg=path):
                C.resolve(r)

    def test_fields_scale_with_pressure_and_key_ignores_voltage(self):
        cfg = C.resolve(self.r)
        f1 = self.P.match_fields(cfg)
        r = copy.deepcopy(self.r)
        r["gas"]["pressure_bar"] *= 2
        self.assertAlmostEqual(self.P.match_fields(C.resolve(r))[0], 2 * f1[0], places=0)
        r = copy.deepcopy(self.r)
        r["hv"]["V_gem"] += 100                  # a voltage scan calibrates once
        self.assertEqual(self.P.calibration_key(C.resolve(r)), self.P.calibration_key(cfg))
        r["gas"]["penning"]["match_fields_kV_cm"] = [50]
        self.assertNotEqual(self.P.calibration_key(C.resolve(r)), self.P.calibration_key(cfg))

    def test_effective_uses_cached_r_micro(self):
        import tempfile
        cfg = C.resolve(self.r)
        plain = C.resolve(raw("gorg"))
        self.assertIs(self.P.effective(plain), plain)                # no matching: unchanged
        with tempfile.TemporaryDirectory() as d:
            old, self.P.CACHE_DIR = self.P.CACHE_DIR, d
            try:
                with open(os.path.join(d, self.P.calibration_key(cfg) + ".json"), "w") as f:
                    json.dump({"r_table": 0.205, "r_micro": 0.235, "lambda": 0.0,
                               "penning_gas": "ar"}, f)
                eff = self.P.effective(cfg)
            finally:
                self.P.CACHE_DIR = old
        self.assertEqual(eff["gas"]["penning"]["r"], 0.235)
        self.assertEqual(eff["gas"]["penning"]["gas"], "ar")
        self.assertEqual(cfg["gas"]["penning"]["r"], "auto")       # input untouched
        args = runner.gas_args(eff)
        self.assertEqual(args[args.index("--penning-r") + 1], "0.235")
        self.assertEqual(args[args.index("--penning-gas") + 1], "ar")

    def test_shipped_calibration_found_from_any_directory(self):
        import tempfile
        cfg = C.resolve(self.r)
        with tempfile.TemporaryDirectory() as work, tempfile.TemporaryDirectory() as shipped:
            old = self.P.CACHE_DIR, self.P.SHIPPED_DIR
            self.P.CACHE_DIR, self.P.SHIPPED_DIR = work, shipped
            try:
                with open(os.path.join(shipped, self.P.calibration_key(cfg) + ".json"), "w") as f:
                    json.dump({"r_table": 0.205, "r_micro": 0.24, "lambda": 0.0,
                               "penning_gas": "ar"}, f)
                eff = self.P.effective(cfg)
            finally:
                self.P.CACHE_DIR, self.P.SHIPPED_DIR = old
            self.assertEqual(os.listdir(work), [])                  # nothing recomputed
        self.assertEqual(eff["gas"]["penning"]["r"], 0.24)


class TestCorrections(unittest.TestCase):
    def test_off_by_default(self):
        from simgem import corrections
        cfg = C.resolve(raw("gorg"))
        self.assertIsNone(cfg["corrections"])
        self.assertEqual(corrections.factor(cfg), 1.0)

    def test_literature_needs_matching_and_scales_per_foil(self):
        from simgem import corrections
        r = raw("gorg")
        r["corrections"] = "literature"
        with self.assertRaises(C.ConfigError):
            C.resolve(r)
        r["gas"]["penning"]["match"] = "table"
        cfg = C.resolve(r)
        L = corrections.LITERATURE
        k = (L["charging_up_per_foil"] * L["gain_scale_per_foil"]) ** 3
        self.assertAlmostEqual(corrections.factor(cfg), k)
        r["foils"] = r["foils"][:1]
        r["gaps"] = [r["gaps"][0], r["gaps"][-1]]
        r.pop("hv")
        r["foils"][0]["dV_V"] = 400
        r["gaps"][0]["field_kV_cm"] = 2.0
        r["gaps"][1]["field_kV_cm"] = 3.0
        self.assertAlmostEqual(corrections.factor(C.resolve(r)), k ** (1 / 3))

    def test_explicit_and_bad(self):
        from simgem import corrections
        r = raw()
        r["corrections"] = {"gain_scale_per_foil": 2.0}
        cfg = C.resolve(r)
        self.assertAlmostEqual(corrections.factor(cfg), 2.0 ** len(cfg["foils"]))
        for bad in ({"gain_scale_per_foil": -1}, {"scale": 2}, "lit"):
            r["corrections"] = bad
            with self.assertRaises(C.ConfigError, msg=str(bad)):
                C.resolve(r)

    def test_corrections_do_not_change_run_key(self):
        r = raw()
        k = runner.run_key(C.resolve(r), "abc", 1)
        r["corrections"] = {"gain_scale_per_foil": 2.0}
        self.assertEqual(runner.run_key(C.resolve(r), "abc", 1), k)


class TestHV(unittest.TestCase):
    # electrode potentials of the TECHTRA techtra_v1 divider at 2600 V (from
    # AGUILA/hardware/high_voltage/gem_divider.py, unloaded)
    TECHTRA_2600 = {"gem3_top": -2600.00, "gem3_bot": -2281.51, "gem2_top": -1702.45,
                    "gem2_bot": -1412.92, "gem1_top": -833.85, "gem1_bot": -579.06}

    def test_gorg_divider_matches_board_model(self):
        cfg = C.resolve(raw("gorg"))
        st = C.Stack(cfg)
        for k, v in self.TECHTRA_2600.items():
            self.assertAlmostEqual(st.V[k], v, places=2, msg=k)
        self.assertAlmostEqual(cfg["gaps"][0]["field_kV_cm"], 500 / 1.8 / 1e3, places=9)
        self.assertAlmostEqual(st.V["drift"], -2600 - 500 / 18.0, places=6)   # 1 mm of 500 V / 18 mm

    def test_shares_equal_resistances(self):
        r = raw("gorg")
        res = r["hv"].pop("divider_resistance_kohm")
        tot = sum(res.values())
        r["hv"]["divider"] = {k: v / tot for k, v in res.items()}
        st = C.Stack(C.resolve(r))
        for k, v in self.TECHTRA_2600.items():
            self.assertAlmostEqual(st.V[k], v, places=2, msg=k)

    def test_hv_scan_changes_field(self):
        r = raw("gorg")
        h1 = C.Stack(C.resolve(r)).field_hash()
        C.set_path(r, "hv.V_gem", 2650)
        self.assertNotEqual(C.Stack(C.resolve(r)).field_hash(), h1)

    def test_hv_validation(self):
        for path, v in (("hv.divider_resistance_kohm.gem1", -1), ("foils[0].dV_V", 300),
                        ("gaps[1].field_kV_cm", 2.0), ("hv.V_drift", 3100), ("hv.divider", {"gem1": 1.0})):
            r = raw("gorg")
            C.set_path(r, path, v)
            with self.assertRaises(C.ConfigError, msg=path):
                C.resolve(r)
        r = raw("gorg")
        del r["hv"]["V_drift_gap"]
        r["hv"]["V_drift"] = 3100.0                       # absolute drift voltage instead
        self.assertAlmostEqual(C.resolve(r)["gaps"][0]["field_kV_cm"], 500 / 1.8 / 1e3, places=9)


class TestLattice(unittest.TestCase):
    def test_cells(self):
        cfg = C.resolve(raw())
        st = C.Stack(cfg)
        self.assertAlmostEqual(st.cell[0], 0.14)
        self.assertAlmostEqual(st.cell[1], 0.14 * math.sqrt(3))

    def test_supercell(self):
        r = raw()
        r["foils"][1]["pitch_um"] = 280.0
        r["foils"][1]["hole"] = {"outer_diameter_um": 140.0, "inner_diameter_um": 100.0}
        st = C.Stack(C.resolve(r))
        self.assertAlmostEqual(st.cell[0], 0.28)
        r["foils"][1]["pitch_um"] = 150.0
        with self.assertRaises(C.ConfigError):
            C.Stack(C.resolve(r))

    def _open_area(self, foil, offset, cell, n=40000):
        """Monte-Carlo fraction of the cell inside a hole (inner radius)."""
        rng = random.Random(3)
        cs = G.hole_centres(foil, offset, cell)
        rI = G.hole_radii(foil)["kapton"][1]
        hit = 0
        for _ in range(n):
            x, y = rng.uniform(-cell[0] / 2, cell[0] / 2), rng.uniform(-cell[1] / 2, cell[1] / 2)
            hit += any((x - a) ** 2 + (y - b) ** 2 < rI * rI for a, b in cs)
        return hit / n

    def test_hole_density_independent_of_offset(self):
        cfg = C.resolve(raw())
        f = cfg["foils"][0]
        cell = C.Stack(cfg).cell
        expect = 2 * math.pi * 0.025 ** 2 / (cell[0] * cell[1])   # 2 holes per hex cell
        for off in [(0, 0), (0.07, 0.0404), (0.013, -0.1)]:
            self.assertAlmostEqual(self._open_area(f, off, cell), expect, delta=0.006)

    def test_wrap(self):
        f = C.resolve(raw())["foils"][0]
        cx, cy = 0.14, 0.14 * math.sqrt(3)
        for d in [(0.2, 0.3), (-0.5, 1.0), (0.07, -0.121)]:
            w = G.wrap_offset(f, d)
            self.assertTrue(-cx / 2 <= w[0] < cx / 2 and -cy / 2 <= w[1] < cy / 2)
            # same pattern: the wrapped offset differs by whole cell periods
            for a, b, c in ((d[0], w[0], cx), (d[1], w[1], cy)):
                k = (a - b) / c
                self.assertAlmostEqual(k, round(k), places=9)

    def test_safe_offsets_nudges_tangency(self):
        cfg = C.resolve(raw())
        # put a hole edge exactly on a wall: centre at x = Lx/2 - r_outer
        st = C.Stack(cfg, {"gem2": (0.07 - 0.035, 0.0)})
        off = G.safe_offsets(st, log=lambda *_: None)
        self.assertFalse(G._too_close(cfg["foils"][1], off["gem2"], st.cell))
        self.assertLess(abs(off["gem2"][0] - 0.035), 0.01)


class TestStatistics(unittest.TestCase):
    def test_chi2_sf(self):
        self.assertAlmostEqual(misalignment.chi2_sf(3.84, 1), 0.05, places=3)
        self.assertAlmostEqual(misalignment.chi2_sf(18.31, 10), 0.05, places=3)
        self.assertAlmostEqual(misalignment.chi2_sf(10, 10), 0.4405, places=3)

    def test_mc_statistics_recovers_geometry_spread(self):
        rng = random.Random(1)
        geoms = []
        for _ in range(400):
            true = 100 + rng.gauss(0, 3)            # geometry spread 3
            geoms.append({"effective_gain": true + rng.gauss(0, 4), "effective_gain_err": 4.0})
        st = misalignment.mc_statistics(geoms, {"events": 1, "effective_gain": 100, "effective_gain_err": 0.1})
        self.assertAlmostEqual(st["spread_geometry_rms"], 3, delta=0.6)
        self.assertLess(abs(st["shift_sigma"]), 3)

    def test_mc_statistics_null(self):
        rng = random.Random(2)
        geoms = [{"effective_gain": 50 + rng.gauss(0, 2), "effective_gain_err": 2.0} for _ in range(200)]
        st = misalignment.mc_statistics(geoms, None)
        self.assertGreater(st["p_value_no_geometry_effect"], 0.001)

    def test_chunk_sizes_cover_all_events(self):
        for total, jobs in ((10, 4), (3, 8), (1000, 6)):
            pending, chunks = total, []
            while pending:
                chunks.append(runner._chunk_size(pending, jobs))
                pending -= chunks[-1]
            self.assertEqual(sum(chunks), total)
            self.assertTrue(all(c >= 1 for c in chunks))
            self.assertGreaterEqual(chunks[0], chunks[-1])      # large chunks first

    def test_ncpu_respects_env_and_grid_slot(self):
        from unittest import mock
        from simgem.util import ncpu
        with mock.patch.dict(os.environ, {"SIMGEM_JOBS": "3", "OMP_NUM_THREADS": "8"}):
            self.assertEqual(ncpu(), 3)
        env = {k: v for k, v in os.environ.items()
               if k not in ("SIMGEM_JOBS", "OMP_NUM_THREADS", "_CONDOR_JOB_AD")}
        with mock.patch.dict(os.environ, dict(env, _CONDOR_SCRATCH_DIR="/x", OMP_NUM_THREADS="2"),
                             clear=True):
            self.assertEqual(ncpu(), 2)
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".ad", delete=False) as ad:
            ad.write('Owner = "someone"\nRequestCpus = 4\n')
        try:
            with mock.patch.dict(os.environ, dict(env, _CONDOR_JOB_AD=ad.name), clear=True):
                self.assertEqual(ncpu(), 4)
        finally:
            os.unlink(ad.name)

    def test_draw_offsets_deterministic(self):
        cfg = C.resolve(raw("misalign_mc"))
        a = misalignment.draw_offsets(cfg, cfg["misalignment_mc"], 3)
        b = misalignment.draw_offsets(cfg, cfg["misalignment_mc"], 3)
        c = misalignment.draw_offsets(cfg, cfg["misalignment_mc"], 4)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)


if __name__ == "__main__":
    unittest.main()
