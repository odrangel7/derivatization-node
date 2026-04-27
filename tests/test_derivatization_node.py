"""
Unit tests for derivatization_node.py
"""

import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

# Allow importing the module from the repo root without installing it.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from derivatization_node import (
    REACTIONS,
    AnnotationResult,
    Feature,
    annotate_features,
    list_reactions,
    load_features_csv,
    main,
    run_derivatization_node,
    save_results_csv,
    summarize_results,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_csv(rows, path):
    """Write a list-of-dicts to a CSV file."""
    if not rows:
        with open(path, "w") as fh:
            fh.write("feature_id,mz,rt,area\n")
        return
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _make_feature(fid="F1", mz=200.0, rt=1.0, area=1e6, charge=1):
    return Feature(feature_id=fid, mz=mz, rt=rt, area=area, charge=charge)


# ---------------------------------------------------------------------------
# Test: REACTIONS registry
# ---------------------------------------------------------------------------

class TestReactionsRegistry(unittest.TestCase):
    def test_all_reactions_have_required_keys(self):
        required = {"name", "functional_group", "mass_shift", "formula_delta", "adduct"}
        for rxn_id, rxn in REACTIONS.items():
            with self.subTest(reaction=rxn_id):
                self.assertTrue(
                    required.issubset(rxn.keys()),
                    f"Missing keys in {rxn_id}: {required - rxn.keys()}",
                )

    def test_mass_shifts_are_positive(self):
        for rxn_id, rxn in REACTIONS.items():
            with self.subTest(reaction=rxn_id):
                self.assertGreater(rxn["mass_shift"], 0)

    def test_list_reactions_returns_all(self):
        listed = list_reactions()
        self.assertEqual(len(listed), len(REACTIONS))
        ids = {r["id"] for r in listed}
        self.assertEqual(ids, set(REACTIONS.keys()))


# ---------------------------------------------------------------------------
# Test: Feature class
# ---------------------------------------------------------------------------

class TestFeature(unittest.TestCase):
    def test_neutral_mass_singly_charged(self):
        f = _make_feature(mz=201.007276, charge=1)
        # neutral mass = 201.007276 * 1 - 1 * 1.007276 = 200.0
        self.assertAlmostEqual(f.neutral_mass(), 200.0, places=4)

    def test_neutral_mass_doubly_charged(self):
        f = _make_feature(mz=101.007276, charge=2)
        # neutral mass = 101.007276 * 2 - 2 * 1.007276 = 200.0
        self.assertAlmostEqual(f.neutral_mass(), 200.0, places=4)


# ---------------------------------------------------------------------------
# Test: CSV I/O
# ---------------------------------------------------------------------------

class TestCSVIO(unittest.TestCase):
    def test_load_features_csv_basic(self):
        rows = [
            {"feature_id": "F1", "mz": "150.0500", "rt": "1.20", "area": "5000"},
            {"feature_id": "F2", "mz": "330.0929", "rt": "1.22", "area": "4500"},
        ]
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, newline=""
        ) as fh:
            path = fh.name
        try:
            _write_csv(rows, path)
            features = load_features_csv(path)
            self.assertEqual(len(features), 2)
            self.assertEqual(features[0].feature_id, "F1")
            self.assertAlmostEqual(features[0].mz, 150.0500)
            self.assertEqual(features[0].charge, 1)
        finally:
            os.unlink(path)

    def test_load_features_csv_with_charge(self):
        rows = [{"feature_id": "F1", "mz": "101.0", "rt": "2.0", "area": "100", "charge": "2"}]
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, newline=""
        ) as fh:
            path = fh.name
        try:
            _write_csv(rows, path)
            features = load_features_csv(path)
            self.assertEqual(features[0].charge, 2)
        finally:
            os.unlink(path)

    def test_load_features_csv_missing_column(self):
        rows = [{"feature_id": "F1", "mz": "150.0"}]  # missing rt, area
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, newline=""
        ) as fh:
            path = fh.name
        try:
            _write_csv(rows, path)
            with self.assertRaises(ValueError):
                load_features_csv(path)
        finally:
            os.unlink(path)

    def test_save_results_csv_empty(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as fh:
            path = fh.name
        try:
            save_results_csv([], path)
            content = Path(path).read_text()
            self.assertIn("No annotations", content)
        finally:
            os.unlink(path)

    def test_save_results_csv_roundtrip(self):
        feat_a = _make_feature("F1", mz=150.05, rt=1.2)
        feat_b = _make_feature("F2", mz=330.0929, rt=1.21)
        result = AnnotationResult(
            underivatized=feat_a,
            derivatized=feat_b,
            reaction_id="DNPH",
            mass_shift_expected=180.0429,
            mass_shift_observed=180.0429,
            delta_mz_da=0.0,
            delta_mz_ppm=0.0,
            delta_rt_min=0.01,
            confidence="High",
        )
        with tempfile.NamedTemporaryFile(
            suffix=".csv", delete=False
        ) as fh:
            path = fh.name
        try:
            save_results_csv([result], path)
            with open(path, newline="") as fh:
                reader = csv.DictReader(fh)
                rows = list(reader)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["underivatized_id"], "F1")
            self.assertEqual(rows[0]["derivatized_id"], "F2")
            self.assertEqual(rows[0]["reaction_id"], "DNPH")
        finally:
            os.unlink(path)


# ---------------------------------------------------------------------------
# Test: annotate_features – core algorithm
# ---------------------------------------------------------------------------

class TestAnnotateFeatures(unittest.TestCase):
    def _make_dnph_pair(self, base_mz=150.05, rt=1.2):
        """Return two features that form a valid DNPH pair."""
        deriv_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        return [
            _make_feature("F_underivatized", mz=base_mz, rt=rt),
            _make_feature("F_derivatized", mz=deriv_mz, rt=rt + 0.05),
        ]

    def test_dnph_pair_detected(self):
        features = self._make_dnph_pair()
        results = annotate_features(
            features,
            reaction_ids=["DNPH"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].reaction_id, "DNPH")
        self.assertEqual(results[0].underivatized.feature_id, "F_underivatized")
        self.assertEqual(results[0].derivatized.feature_id, "F_derivatized")

    def test_no_pair_outside_mz_tolerance(self):
        base_mz = 150.05
        # Place derivatized feature 1 Da away from expected (way outside tolerance)
        features = [
            _make_feature("F1", mz=base_mz, rt=1.0),
            _make_feature("F2", mz=base_mz + REACTIONS["DNPH"]["mass_shift"] + 1.0, rt=1.0),
        ]
        results = annotate_features(
            features,
            reaction_ids=["DNPH"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
        )
        self.assertEqual(len(results), 0)

    def test_no_pair_outside_rt_tolerance(self):
        base_mz = 150.05
        deriv_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        features = [
            _make_feature("F1", mz=base_mz, rt=1.0),
            _make_feature("F2", mz=deriv_mz, rt=3.0),  # 2 min apart
        ]
        results = annotate_features(
            features,
            reaction_ids=["DNPH"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
        )
        self.assertEqual(len(results), 0)

    def test_absolute_da_tolerance(self):
        base_mz = 150.05
        deriv_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        features = [
            _make_feature("F1", mz=base_mz, rt=1.0),
            _make_feature("F2", mz=deriv_mz, rt=1.0),
        ]
        results = annotate_features(
            features,
            reaction_ids=["DNPH"],
            mz_tolerance_da=0.01,
            rt_tolerance_min=0.3,
        )
        self.assertEqual(len(results), 1)

    def test_mutual_exclusive_tolerance_raises(self):
        with self.assertRaises(ValueError):
            annotate_features(
                [],
                mz_tolerance_da=0.01,
                mz_tolerance_ppm=5.0,
            )

    def test_unknown_reaction_raises(self):
        with self.assertRaises(ValueError):
            annotate_features([], reaction_ids=["UNKNOWN_RXN"])

    def test_multiple_reactions(self):
        base_mz = 200.0
        dnph_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        dansyl_mz = base_mz + REACTIONS["DANSYL"]["mass_shift"]
        features = [
            _make_feature("F_base", mz=base_mz, rt=1.0),
            _make_feature("F_dnph", mz=dnph_mz, rt=1.0),
            _make_feature("F_dansyl", mz=dansyl_mz, rt=1.0),
        ]
        results = annotate_features(
            features,
            reaction_ids=["DNPH", "DANSYL"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
        )
        rxn_ids = {r.reaction_id for r in results}
        self.assertIn("DNPH", rxn_ids)
        self.assertIn("DANSYL", rxn_ids)

    def test_multi_site_derivatization(self):
        """Double acetylation (2 × 42.0106 Da) should be detected."""
        base_mz = 200.0
        double_acetyl_mz = base_mz + REACTIONS["ACETYLATION"]["mass_shift"] * 2
        features = [
            _make_feature("F_base", mz=base_mz, rt=1.0),
            _make_feature("F_double", mz=double_acetyl_mz, rt=1.0),
        ]
        results = annotate_features(
            features,
            reaction_ids=["ACETYLATION"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
            max_derivatizations=3,
        )
        # At least one result corresponds to double derivatization
        double_results = [
            r for r in results
            if abs(r.mass_shift_expected - REACTIONS["ACETYLATION"]["mass_shift"] * 2) < 0.001
        ]
        self.assertTrue(double_results, "Double-site acetylation not detected")

    def test_confidence_high_for_exact_match(self):
        base_mz = 150.05
        deriv_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        features = [
            _make_feature("F1", mz=base_mz, rt=1.0),
            _make_feature("F2", mz=deriv_mz, rt=1.0),  # exact RT match
        ]
        results = annotate_features(
            features,
            reaction_ids=["DNPH"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].confidence, "High")

    def test_confidence_low_for_boundary_match(self):
        base_mz = 150.05
        deriv_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        features = [
            _make_feature("F1", mz=base_mz, rt=1.0),
            # RT just inside tolerance
            _make_feature("F2", mz=deriv_mz, rt=1.0 + 0.29),
        ]
        results = annotate_features(
            features,
            reaction_ids=["DNPH"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].confidence, "Low")


# ---------------------------------------------------------------------------
# Test: summarize_results
# ---------------------------------------------------------------------------

class TestSummarize(unittest.TestCase):
    def test_empty_summary(self):
        summary = summarize_results([])
        self.assertEqual(summary["total_annotations"], 0)
        self.assertEqual(summary["by_confidence"]["High"], 0)

    def test_summary_counts(self):
        feat_a = _make_feature("A", mz=150.05, rt=1.0)
        feat_b = _make_feature("B", mz=330.09, rt=1.0)
        r = AnnotationResult(
            underivatized=feat_a,
            derivatized=feat_b,
            reaction_id="DNPH",
            mass_shift_expected=180.0429,
            mass_shift_observed=180.04,
            delta_mz_da=0.0,
            delta_mz_ppm=0.0,
            delta_rt_min=0.0,
            confidence="High",
        )
        summary = summarize_results([r])
        self.assertEqual(summary["total_annotations"], 1)
        self.assertEqual(summary["by_reaction"]["DNPH"], 1)
        self.assertEqual(summary["by_confidence"]["High"], 1)


# ---------------------------------------------------------------------------
# Test: run_derivatization_node (CD integration entry point)
# ---------------------------------------------------------------------------

class TestRunDerivNode(unittest.TestCase):
    def test_basic_cd_integration(self):
        base_mz = 150.05
        deriv_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        feature_table = [
            {"feature_id": "F1", "mz": base_mz, "rt": 1.0, "area": 1e6},
            {"feature_id": "F2", "mz": deriv_mz, "rt": 1.0, "area": 9e5},
        ]
        output = run_derivatization_node(
            feature_table,
            reaction_ids=["DNPH"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
        )
        self.assertEqual(len(output), 1)
        self.assertEqual(output[0]["reaction_id"], "DNPH")

    def test_cd_integration_alternate_column_names(self):
        """CD may supply columns with different names."""
        base_mz = 150.05
        deriv_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        feature_table = [
            {"ID": "F1", "m/z": base_mz, "rt": 1.0, "area": 1e6},
            {"ID": "F2", "m/z": deriv_mz, "rt": 1.05, "area": 9e5},
        ]
        output = run_derivatization_node(
            feature_table,
            reaction_ids=["DNPH"],
            mz_tolerance_ppm=5.0,
            rt_tolerance_min=0.3,
        )
        self.assertEqual(len(output), 1)


# ---------------------------------------------------------------------------
# Test: CLI (main)
# ---------------------------------------------------------------------------

class TestCLI(unittest.TestCase):
    def _run(self, argv):
        return main(argv)

    def test_list_reactions(self):
        rc = self._run(["--list-reactions"])
        self.assertEqual(rc, 0)

    def test_version(self):
        with self.assertRaises(SystemExit) as cm:
            self._run(["--version"])
        self.assertEqual(cm.exception.code, 0)

    def test_missing_input_flag(self):
        with self.assertRaises(SystemExit) as cm:
            self._run(["--output", "out.csv"])
        self.assertNotEqual(cm.exception.code, 0)

    def test_full_pipeline(self):
        base_mz = 150.05
        deriv_mz = base_mz + REACTIONS["DNPH"]["mass_shift"]
        rows = [
            {"feature_id": "F1", "mz": base_mz, "rt": "1.0", "area": "1000000"},
            {"feature_id": "F2", "mz": deriv_mz, "rt": "1.05", "area": "900000"},
        ]
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, newline=""
        ) as fh:
            in_path = fh.name
        with tempfile.NamedTemporaryFile(
            suffix=".csv", delete=False
        ) as fh:
            out_path = fh.name
        try:
            _write_csv(rows, in_path)
            rc = self._run([
                "--input", in_path,
                "--output", out_path,
                "--reaction", "DNPH",
                "--mz-tolerance-ppm", "5",
                "--rt-tolerance", "0.3",
                "--summary",
            ])
            self.assertEqual(rc, 0)
            with open(out_path, newline="") as fh:
                reader = csv.DictReader(fh)
                result_rows = list(reader)
            self.assertEqual(len(result_rows), 1)
            self.assertEqual(result_rows[0]["reaction_id"], "DNPH")
        finally:
            os.unlink(in_path)
            os.unlink(out_path)

    def test_invalid_input_file(self):
        rc = self._run([
            "--input", "/nonexistent/path/features.csv",
            "--output", "/tmp/out.csv",
        ])
        self.assertEqual(rc, 1)

    def test_unknown_reaction_id_cli(self):
        rows = [{"feature_id": "F1", "mz": "150.0", "rt": "1.0", "area": "1000"}]
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, newline=""
        ) as fh:
            in_path = fh.name
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as fh:
            out_path = fh.name
        try:
            _write_csv(rows, in_path)
            rc = self._run([
                "--input", in_path,
                "--output", out_path,
                "--reaction", "TOTALLY_FAKE_REACTION",
            ])
            self.assertEqual(rc, 1)
        finally:
            os.unlink(in_path)
            os.unlink(out_path)


if __name__ == "__main__":
    unittest.main()
