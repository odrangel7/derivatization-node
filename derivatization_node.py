"""
Derivatization-Driven Annotation Processing Node
=================================================
Processing Node for Derivatization-Driven Annotation via Functional
Group-Specific Reactions Using Compound Discoverer Software.

This module implements the core algorithm for pairing underivatized and
derivatized features in a metabolomics feature table and annotating them
with their corresponding functional-group-specific reactions.

Usage (standalone CLI)
----------------------
    python derivatization_node.py --input features.csv --output annotated.csv \\
        --reaction DNPH --mz-tolerance 0.005 --rt-tolerance 0.3

Usage (Compound Discoverer scripted node)
-----------------------------------------
    Import and call ``run_derivatization_node()`` from the CD scripting
    entry point.  See README for the full integration guide.
"""

import argparse
import csv
import json
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import numpy as np
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:  # pragma: no cover
    _HAS_PANDAS = False

__version__ = "1.0.0"
__author__ = "Oscar Daniel Rangel Huerta"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Reaction registry
# ---------------------------------------------------------------------------

#: Built-in derivatization reactions.
#: Each entry contains:
#:   name            – human-readable reagent name
#:   functional_group – targeted functional group(s)
#:   mass_shift      – monoisotopic mass added by the reaction (Da)
#:   formula_delta   – elemental composition change (added – lost)
#:   adduct          – typical ion notation observed in positive mode
REACTIONS: Dict[str, Dict] = {
    "DNPH": {
        "name": "2,4-Dinitrophenylhydrazine (DNPH)",
        "functional_group": "Carbonyl (aldehydes, ketones)",
        "mass_shift": 180.0429,   # C6H3N3O4 added, H2O lost
        "formula_delta": "C6H3N3O4 – H2O",
        "adduct": "[M+DNPH−H2O+H]+",
    },
    "GIRARD_P": {
        "name": "Girard's Reagent P",
        "functional_group": "Carbonyl (aldehydes, ketones)",
        "mass_shift": 138.0793,   # C7H10N2O added, H2O lost
        "formula_delta": "C7H10N2O – H2O",
        "adduct": "[M+GirardP−H2O+H]+",
    },
    "DANSYL": {
        "name": "Dansyl Chloride (DNS-Cl)",
        "functional_group": "Primary and secondary amines, phenols",
        "mass_shift": 233.0567,   # C12H11NO2S added, HCl lost
        "formula_delta": "C12H11NO2S – HCl",
        "adduct": "[M+Dansyl+H]+",
    },
    "FMOC": {
        "name": "Fluorenylmethyloxycarbonyl Chloride (Fmoc-Cl)",
        "functional_group": "Primary amines",
        "mass_shift": 222.0681,   # C15H10O2 added, HCl lost
        "formula_delta": "C15H10O2 – HCl",
        "adduct": "[M+Fmoc+H]+",
    },
    "ACETYLATION": {
        "name": "Acetic Anhydride / Acetyl Chloride",
        "functional_group": "Hydroxyls, primary and secondary amines",
        "mass_shift": 42.0106,    # C2H2O added per site
        "formula_delta": "C2H2O",
        "adduct": "[M+Ac+H]+",
    },
    "NBD_CL": {
        "name": "4-Chloro-7-nitrobenzofurazan (NBD-Cl)",
        "functional_group": "Amines, thiols",
        "mass_shift": 166.0253,   # C6H4N3O3 added, HCl lost
        "formula_delta": "C6H4N3O3 – HCl",
        "adduct": "[M+NBD+H]+",
    },
    "OPA": {
        "name": "o-Phthalaldehyde (OPA)",
        "functional_group": "Primary amines (with thiol co-reagent)",
        "mass_shift": 116.0262,   # C8H4O added, H2O lost
        "formula_delta": "C8H4O – H2O",
        "adduct": "[M+OPA+H]+",
    },
    "TMS": {
        "name": "Trimethylsilylation (TMS)",
        "functional_group": "Hydroxyls, carboxylic acids, amines",
        "mass_shift": 72.0395,    # C3H8Si added, H lost per site
        "formula_delta": "C3H8Si – H",
        "adduct": "[M+TMS+H]+",
    },
}


def list_reactions() -> List[Dict]:
    """Return all built-in reactions as a list of dicts (for reporting)."""
    return [
        {"id": k, **v} for k, v in REACTIONS.items()
    ]


# ---------------------------------------------------------------------------
# Core data classes
# ---------------------------------------------------------------------------

@dataclass
class Feature:
    """A single chromatographic feature (peak) in the feature table."""

    feature_id: str
    mz: float
    rt: float
    area: float
    charge: int = 1

    # Annotation fields (populated by ``annotate_features``)
    reaction_id: Optional[str] = None
    functional_group: Optional[str] = None
    paired_feature_id: Optional[str] = None
    mass_shift_observed: Optional[float] = None
    mass_shift_expected: Optional[float] = None
    delta_mz_ppm: Optional[float] = None
    delta_rt_min: Optional[float] = None
    annotation_confidence: Optional[str] = None  # "High" | "Medium" | "Low"

    def neutral_mass(self) -> float:
        """Compute the neutral monoisotopic mass from observed m/z."""
        proton = 1.007276
        return self.mz * self.charge - self.charge * proton


@dataclass
class AnnotationResult:
    """Holds a paired annotation between an underivatized and a derivatized feature."""

    underivatized: Feature
    derivatized: Feature
    reaction_id: str
    mass_shift_expected: float
    mass_shift_observed: float
    delta_mz_da: float
    delta_mz_ppm: float
    delta_rt_min: float
    confidence: str  # "High" | "Medium" | "Low"

    def to_dict(self) -> Dict:
        return {
            "underivatized_id": self.underivatized.feature_id,
            "underivatized_mz": self.underivatized.mz,
            "underivatized_rt": self.underivatized.rt,
            "underivatized_area": self.underivatized.area,
            "derivatized_id": self.derivatized.feature_id,
            "derivatized_mz": self.derivatized.mz,
            "derivatized_rt": self.derivatized.rt,
            "derivatized_area": self.derivatized.area,
            "reaction_id": self.reaction_id,
            "reaction_name": REACTIONS[self.reaction_id]["name"],
            "functional_group": REACTIONS[self.reaction_id]["functional_group"],
            "mass_shift_expected_da": self.mass_shift_expected,
            "mass_shift_observed_da": self.mass_shift_observed,
            "delta_mz_da": self.delta_mz_da,
            "delta_mz_ppm": self.delta_mz_ppm,
            "delta_rt_min": self.delta_rt_min,
            "annotation_confidence": self.confidence,
            "adduct": REACTIONS[self.reaction_id]["adduct"],
        }


# ---------------------------------------------------------------------------
# Feature I/O helpers
# ---------------------------------------------------------------------------

_REQUIRED_COLS = {"feature_id", "mz", "rt", "area"}


def load_features_csv(path: str) -> List[Feature]:
    """Load features from a CSV file.

    Expected columns (case-insensitive): feature_id, mz, rt, area, [charge].

    Parameters
    ----------
    path:
        Path to the CSV file.

    Returns
    -------
    List of :class:`Feature` objects.
    """
    features: List[Feature] = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower().strip() for c in (reader.fieldnames or [])}
        missing = _REQUIRED_COLS - cols
        if missing:
            raise ValueError(
                f"Input CSV is missing required columns: {missing}. "
                f"Found columns: {set(reader.fieldnames or [])}"
            )
        for row in reader:
            row_lower = {k.lower().strip(): v for k, v in row.items()}
            charge_raw = row_lower.get("charge", "1").strip()
            try:
                charge = int(charge_raw) if charge_raw else 1
            except ValueError:
                logger.warning(
                    "Invalid charge value '%s' for feature '%s'; defaulting to 1.",
                    charge_raw,
                    row_lower.get("feature_id", "<unknown>"),
                )
                charge = 1
            features.append(
                Feature(
                    feature_id=row_lower["feature_id"].strip(),
                    mz=float(row_lower["mz"]),
                    rt=float(row_lower["rt"]),
                    area=float(row_lower["area"]),
                    charge=charge,
                )
            )
    logger.info("Loaded %d features from '%s'", len(features), path)
    return features


def save_results_csv(results: List[AnnotationResult], path: str) -> None:
    """Write annotation results to a CSV file.

    Parameters
    ----------
    results:
        List of :class:`AnnotationResult` objects.
    path:
        Output path for the CSV file.
    """
    if not results:
        logger.warning("No annotation results to write.")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("# No annotations found.\n")
        return

    rows = [r.to_dict() for r in results]
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d annotation results to '%s'", len(results), path)


# ---------------------------------------------------------------------------
# Core annotation algorithm
# ---------------------------------------------------------------------------

def _ppm_error(mz_observed: float, mz_expected: float) -> float:
    """Return the mass accuracy error in parts-per-million."""
    return abs(mz_observed - mz_expected) / mz_expected * 1e6


def _confidence_score(delta_mz_ppm: float, delta_rt_min: float,
                      mz_tol_ppm: float, rt_tol_min: float) -> str:
    """Assign a categorical confidence level to an annotation.

    Parameters
    ----------
    delta_mz_ppm:
        Observed mass accuracy error (ppm).
    delta_rt_min:
        Observed retention-time difference (min).
    mz_tol_ppm:
        User-specified m/z tolerance (ppm).
    rt_tol_min:
        User-specified RT tolerance (min).

    Returns
    -------
    ``"High"``, ``"Medium"``, or ``"Low"``.
    """
    mz_frac = delta_mz_ppm / mz_tol_ppm if mz_tol_ppm > 0 else 1.0
    rt_frac = delta_rt_min / rt_tol_min if rt_tol_min > 0 else 1.0
    combined = max(mz_frac, rt_frac)
    if combined <= 0.33:
        return "High"
    if combined <= 0.67:
        return "Medium"
    return "Low"


def annotate_features(
    features: List[Feature],
    reaction_ids: Optional[List[str]] = None,
    mz_tolerance_da: Optional[float] = None,
    mz_tolerance_ppm: Optional[float] = None,
    rt_tolerance_min: float = 0.3,
    max_derivatizations: int = 3,
) -> List[AnnotationResult]:
    """Pair underivatized features with their derivatized counterparts.

    For every feature in *features*, the function calculates the expected
    m/z of the derivatized adduct for each reaction and searches the feature
    table for a match within the specified tolerances.

    Parameters
    ----------
    features:
        Input feature list (all samples combined).
    reaction_ids:
        Subset of reaction keys from :data:`REACTIONS` to consider.
        Defaults to all built-in reactions.
    mz_tolerance_da:
        Absolute m/z matching window (Da).  Mutually exclusive with
        *mz_tolerance_ppm*.
    mz_tolerance_ppm:
        Relative m/z matching window (ppm).  If neither *mz_tolerance_da* nor
        *mz_tolerance_ppm* is provided, defaults to 5 ppm.
    rt_tolerance_min:
        Maximum allowed retention-time difference between the underivatized
        and derivatized feature (minutes).  Default 0.3 min.
    max_derivatizations:
        Maximum number of derivatization sites to consider (i.e. the reaction
        mass shift may be multiplied by 1 … *max_derivatizations*).

    Returns
    -------
    List of :class:`AnnotationResult` objects, one per matched pair.
    """
    if reaction_ids is None:
        reaction_ids = list(REACTIONS.keys())
    else:
        unknown = set(reaction_ids) - set(REACTIONS.keys())
        if unknown:
            raise ValueError(
                f"Unknown reaction ID(s): {unknown}. "
                f"Available: {list(REACTIONS.keys())}"
            )

    # Resolve tolerance settings
    if mz_tolerance_da is not None and mz_tolerance_ppm is not None:
        raise ValueError(
            "Specify either mz_tolerance_da or mz_tolerance_ppm, not both."
        )
    use_ppm = mz_tolerance_da is None
    if use_ppm:
        mz_tol_ppm = mz_tolerance_ppm if mz_tolerance_ppm is not None else 5.0
    else:
        mz_tol_ppm = 0.0  # will compute per feature

    results: List[AnnotationResult] = []
    # Track which features have already been used as derivatized partners
    used_derivatized: set = set()

    for rxn_id in reaction_ids:
        rxn = REACTIONS[rxn_id]
        base_shift = rxn["mass_shift"]

        for n_sites in range(1, max_derivatizations + 1):
            expected_shift = base_shift * n_sites

            for feat in features:
                # Predicted m/z of derivatized ion (positive mode, charge 1)
                expected_deriv_mz = feat.mz + expected_shift / feat.charge

                # Determine absolute tolerance in Da
                if use_ppm:
                    tol_da = mz_tol_ppm * expected_deriv_mz / 1e6
                else:
                    tol_da = mz_tolerance_da  # type: ignore[assignment]

                for candidate in features:
                    if candidate.feature_id == feat.feature_id:
                        continue
                    if candidate.feature_id in used_derivatized:
                        continue

                    delta_mz = abs(candidate.mz - expected_deriv_mz)
                    if delta_mz > tol_da:
                        continue

                    delta_rt = abs(candidate.rt - feat.rt)
                    if delta_rt > rt_tolerance_min:
                        continue

                    # Compute PPM error (always useful to report)
                    ppm_err = _ppm_error(candidate.mz, expected_deriv_mz)

                    # Compute effective ppm tolerance for confidence scoring
                    eff_ppm_tol = (
                        mz_tol_ppm
                        if use_ppm
                        else (tol_da / candidate.mz * 1e6)
                    )

                    confidence = _confidence_score(
                        ppm_err, delta_rt, eff_ppm_tol, rt_tolerance_min
                    )

                    mass_shift_obs = (
                        (candidate.mz - feat.mz) * feat.charge
                    )

                    result = AnnotationResult(
                        underivatized=feat,
                        derivatized=candidate,
                        reaction_id=rxn_id,
                        mass_shift_expected=expected_shift,
                        mass_shift_observed=round(mass_shift_obs, 6),
                        delta_mz_da=round(candidate.mz - expected_deriv_mz, 6),
                        delta_mz_ppm=round(ppm_err, 3),
                        delta_rt_min=round(delta_rt, 4),
                        confidence=confidence,
                    )
                    results.append(result)
                    used_derivatized.add(candidate.feature_id)
                    logger.debug(
                        "Paired %s ↔ %s  via %s (Δm/z=%.4f Da, %s confidence)",
                        feat.feature_id,
                        candidate.feature_id,
                        rxn_id,
                        result.delta_mz_da,
                        confidence,
                    )

    logger.info(
        "Annotation complete: %d pair(s) found across %d reaction(s).",
        len(results),
        len(reaction_ids),
    )
    return results


# ---------------------------------------------------------------------------
# Summary utilities
# ---------------------------------------------------------------------------

def summarize_results(results: List[AnnotationResult]) -> Dict:
    """Produce a summary dictionary from a list of annotation results.

    Parameters
    ----------
    results:
        Output of :func:`annotate_features`.

    Returns
    -------
    Dictionary with total counts, per-reaction counts, and confidence
    distribution.
    """
    if not results:
        return {
            "total_annotations": 0,
            "by_reaction": {},
            "by_confidence": {"High": 0, "Medium": 0, "Low": 0},
        }

    by_reaction: Dict[str, int] = {}
    by_confidence: Dict[str, int] = {"High": 0, "Medium": 0, "Low": 0}

    for r in results:
        by_reaction[r.reaction_id] = by_reaction.get(r.reaction_id, 0) + 1
        by_confidence[r.confidence] = by_confidence.get(r.confidence, 0) + 1

    return {
        "total_annotations": len(results),
        "by_reaction": by_reaction,
        "by_confidence": by_confidence,
    }


# ---------------------------------------------------------------------------
# Compound Discoverer integration entry point
# ---------------------------------------------------------------------------

def run_derivatization_node(
    feature_table: List[Dict],
    reaction_ids: Optional[List[str]] = None,
    mz_tolerance_ppm: float = 5.0,
    rt_tolerance_min: float = 0.3,
    max_derivatizations: int = 3,
) -> List[Dict]:
    """Entry point for Compound Discoverer scripted processing nodes.

    This function accepts the feature table as a list of dictionaries (as
    returned by the CD scripting API) and returns the annotation results in
    the same format.

    Parameters
    ----------
    feature_table:
        List of feature dicts.  Each dict must contain at least the keys
        ``feature_id`` (or ``"ID"``), ``mz`` (or ``"MolecularWeight"``/
        ``"m/z"``), ``rt``, and ``area`` (or ``"Area"``).
    reaction_ids:
        Reaction IDs to apply (default: all built-in reactions).
    mz_tolerance_ppm:
        m/z matching tolerance in ppm (default 5).
    rt_tolerance_min:
        RT matching tolerance in minutes (default 0.3).
    max_derivatizations:
        Maximum number of reactive sites per molecule (default 3).

    Returns
    -------
    List of annotation-result dicts suitable for injection back into CD.
    """
    # Normalize the dict keys coming from CD's API
    features: List[Feature] = []
    for row in feature_table:
        low = {k.lower().strip(): v for k, v in row.items()}
        fid = (
            str(low.get("feature_id") or low.get("id") or low.get("compound_id", ""))
        )
        mz_val = float(
            low.get("mz") or low.get("m/z") or low.get("molecularweight", 0)
        )
        rt_val = float(low.get("rt") or low.get("retentiontime", 0))
        area_val = float(low.get("area") or low.get("peakarea", 0))
        charge_val = int(low.get("charge", 1) or 1)
        features.append(
            Feature(
                feature_id=fid,
                mz=mz_val,
                rt=rt_val,
                area=area_val,
                charge=charge_val,
            )
        )

    annotation_results = annotate_features(
        features,
        reaction_ids=reaction_ids,
        mz_tolerance_ppm=mz_tolerance_ppm,
        rt_tolerance_min=rt_tolerance_min,
        max_derivatizations=max_derivatizations,
    )
    return [r.to_dict() for r in annotation_results]


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="derivatization_node",
        description=(
            "Derivatization-Driven Annotation Processing Node – annotate MS "
            "feature tables using functional group-specific reaction mass shifts."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(
            [
                "Available reaction IDs:",
                *[
                    f"  {k:15s} {v['functional_group']}"
                    for k, v in REACTIONS.items()
                ],
                "",
                "Examples:",
                "  python derivatization_node.py --input features.csv "
                "--output annotated.csv --reaction DNPH",
                "  python derivatization_node.py --input features.csv "
                "--output annotated.csv --mz-tolerance-ppm 5 --rt-tolerance 0.2",
                "  python derivatization_node.py --list-reactions",
            ]
        ),
    )
    parser.add_argument(
        "--input", "-i",
        metavar="CSV",
        help="Input feature table (CSV).  Required columns: feature_id, mz, rt, area.",
    )
    parser.add_argument(
        "--output", "-o",
        metavar="CSV",
        help="Output annotation table (CSV).",
    )
    parser.add_argument(
        "--reaction", "-r",
        metavar="ID",
        action="append",
        dest="reactions",
        help=(
            "Reaction ID(s) to apply.  May be specified multiple times.  "
            "Default: all built-in reactions."
        ),
    )
    tol_group = parser.add_mutually_exclusive_group()
    tol_group.add_argument(
        "--mz-tolerance-da",
        type=float,
        metavar="DA",
        help="Absolute m/z matching tolerance in Daltons.",
    )
    tol_group.add_argument(
        "--mz-tolerance-ppm",
        type=float,
        default=5.0,
        metavar="PPM",
        help="Relative m/z matching tolerance in ppm (default: 5).",
    )
    parser.add_argument(
        "--rt-tolerance",
        type=float,
        default=0.3,
        metavar="MIN",
        help="Retention-time matching tolerance in minutes (default: 0.3).",
    )
    parser.add_argument(
        "--max-sites",
        type=int,
        default=3,
        metavar="N",
        help="Maximum number of derivatization sites per compound (default: 3).",
    )
    parser.add_argument(
        "--list-reactions",
        action="store_true",
        help="Print all available reaction IDs and exit.",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print a summary of annotation results to stdout.",
    )
    parser.add_argument(
        "--version", "-v",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point.  Returns 0 on success, non-zero on error."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.list_reactions:
        header = f"{'ID':<15} {'Functional Group':<40} {'Mass Shift (Da)':<18} Reagent"
        print(header)
        print("-" * len(header))
        for rxn in list_reactions():
            print(
                f"{rxn['id']:<15} {rxn['functional_group']:<40} "
                f"{rxn['mass_shift']:<18.4f} {rxn['name']}"
            )
        return 0

    if not args.input:
        parser.error("--input is required (unless --list-reactions is used).")
    if not args.output:
        parser.error("--output is required (unless --list-reactions is used).")

    try:
        features = load_features_csv(args.input)
    except (FileNotFoundError, ValueError) as exc:
        logger.error("Failed to load input file: %s", exc)
        return 1

    try:
        results = annotate_features(
            features,
            reaction_ids=args.reactions,
            mz_tolerance_da=args.mz_tolerance_da,
            mz_tolerance_ppm=(
                args.mz_tolerance_ppm if args.mz_tolerance_da is None else None
            ),
            rt_tolerance_min=args.rt_tolerance,
            max_derivatizations=args.max_sites,
        )
    except ValueError as exc:
        logger.error("Annotation error: %s", exc)
        return 1

    save_results_csv(results, args.output)

    if args.summary:
        summary = summarize_results(results)
        print(json.dumps(summary, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())
