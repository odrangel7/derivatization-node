# derivatization-node

**Processing Node for Derivatization-Driven Annotation via Functional Group-Specific Reactions Using Compound Discoverer Software**

This repository accompanies the manuscript:

> *Derivatization-Driven Annotation via Functional Group-Specific Reactions Using Compound Discoverer Software*  
> Oscar Daniel Rangel Huerta *et al.*

---

## Overview

`derivatization_node.py` is a Python processing node that annotates metabolomics feature tables by detecting mass-shift patterns introduced by functional group-specific chemical derivatization reagents.

The algorithm:

1. Reads a feature table (m/z, retention time, area) from a single or combined sample run.
2. For each specified derivatization reaction, calculates the expected m/z of the derivatized adduct for every feature.
3. Searches the table for a matching feature within user-defined m/z and retention-time tolerances.
4. Assigns a confidence level (High / Medium / Low) based on how closely the match fits within the tolerances.
5. Outputs an annotated CSV table and an optional summary.

The module can be used:
- **Standalone** via the command-line interface (CLI).
- **Integrated into Compound Discoverer** as a scripted processing node via the `run_derivatization_node()` entry point.

---

## Supported Derivatization Reactions

| ID           | Reagent                                         | Targeted Functional Group              | Mass Shift (Da) |
|--------------|-------------------------------------------------|----------------------------------------|-----------------|
| `DNPH`       | 2,4-Dinitrophenylhydrazine                      | Carbonyl (aldehydes, ketones)          | 180.0429        |
| `GIRARD_P`   | Girard's Reagent P                              | Carbonyl (aldehydes, ketones)          | 138.0793        |
| `DANSYL`     | Dansyl Chloride (DNS-Cl)                        | Primary/secondary amines, phenols     | 233.0567        |
| `FMOC`       | Fluorenylmethyloxycarbonyl Chloride (Fmoc-Cl)  | Primary amines                         | 222.0681        |
| `ACETYLATION`| Acetic Anhydride / Acetyl Chloride              | Hydroxyls, primary/secondary amines   | 42.0106         |
| `NBD_CL`     | 4-Chloro-7-nitrobenzofurazan (NBD-Cl)           | Amines, thiols                         | 166.0253        |
| `OPA`        | o-Phthalaldehyde (OPA)                         | Primary amines (with thiol co-reagent) | 116.0262        |
| `TMS`        | Trimethylsilylation (TMS)                       | Hydroxyls, carboxylic acids, amines   | 72.0395         |

To list all reactions from the command line:

```bash
python derivatization_node.py --list-reactions
```

Custom reactions can be added by extending the `REACTIONS` dictionary in `derivatization_node.py`.

---

## Requirements

- Python ≥ 3.9
- numpy ≥ 1.24
- pandas ≥ 2.0 *(optional – used only if you import `pandas` directly)*

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Installation

No installation is required. Copy `derivatization_node.py` to your project directory and run it directly.

---

## Input Format

The input feature table must be a **CSV file** with at minimum the following columns (case-insensitive):

| Column       | Description                                       |
|--------------|---------------------------------------------------|
| `feature_id` | Unique identifier for the feature (string)        |
| `mz`         | Observed m/z value (float)                        |
| `rt`         | Retention time in minutes (float)                 |
| `area`       | Peak area / intensity (float)                     |
| `charge`     | *(Optional)* Ion charge state, default `1` (int)  |

Example `features.csv`:

```
feature_id,mz,rt,area
F001,150.0550,1.20,1500000
F002,330.0979,1.21,1200000
F003,201.1234,2.45,980000
F004,381.1663,2.47,870000
```

---

## Usage – Command-Line Interface

### Annotate with all built-in reactions

```bash
python derivatization_node.py \
    --input features.csv \
    --output annotated.csv \
    --mz-tolerance-ppm 5 \
    --rt-tolerance 0.3
```

### Annotate with a specific reaction

```bash
python derivatization_node.py \
    --input features.csv \
    --output annotated.csv \
    --reaction DNPH \
    --mz-tolerance-ppm 5 \
    --rt-tolerance 0.3
```

### Apply multiple reactions and print a summary

```bash
python derivatization_node.py \
    --input features.csv \
    --output annotated.csv \
    --reaction DNPH \
    --reaction DANSYL \
    --mz-tolerance-ppm 5 \
    --rt-tolerance 0.3 \
    --summary
```

### Use absolute m/z tolerance (Da) instead of ppm

```bash
python derivatization_node.py \
    --input features.csv \
    --output annotated.csv \
    --mz-tolerance-da 0.005 \
    --rt-tolerance 0.3
```

### All options

```
usage: derivatization_node [--input CSV] [--output CSV]
                           [--reaction ID] ...
                           [--mz-tolerance-ppm PPM | --mz-tolerance-da DA]
                           [--rt-tolerance MIN]
                           [--max-sites N]
                           [--list-reactions] [--summary] [--version]

optional arguments:
  --input,  -i CSV         Input feature table CSV (required)
  --output, -o CSV         Output annotation table CSV (required)
  --reaction, -r ID        Reaction ID(s) to apply (repeatable; default: all)
  --mz-tolerance-ppm PPM   m/z matching window in ppm (default: 5)
  --mz-tolerance-da DA     m/z matching window in Da (mutually exclusive with ppm)
  --rt-tolerance MIN       RT matching window in minutes (default: 0.3)
  --max-sites N            Maximum reactive sites per compound (default: 3)
  --list-reactions         Print all available reaction IDs and exit
  --summary                Print a JSON summary to stdout after annotation
  --version, -v            Show version and exit
```

---

## Output Format

The output CSV contains one row per annotated feature pair with the following columns:

| Column                  | Description                                             |
|-------------------------|---------------------------------------------------------|
| `underivatized_id`      | Feature ID of the underivatized compound                |
| `underivatized_mz`      | m/z of the underivatized feature                        |
| `underivatized_rt`      | RT of the underivatized feature (min)                   |
| `underivatized_area`    | Peak area of the underivatized feature                  |
| `derivatized_id`        | Feature ID of the derivatized compound                  |
| `derivatized_mz`        | m/z of the derivatized feature                          |
| `derivatized_rt`        | RT of the derivatized feature (min)                     |
| `derivatized_area`      | Peak area of the derivatized feature                    |
| `reaction_id`           | Reaction key (e.g., `DNPH`)                             |
| `reaction_name`         | Human-readable reagent name                             |
| `functional_group`      | Targeted functional group                               |
| `mass_shift_expected_da`| Expected mass shift (Da)                                |
| `mass_shift_observed_da`| Observed mass shift (Da)                                |
| `delta_mz_da`           | Residual m/z error (Da)                                 |
| `delta_mz_ppm`          | Residual m/z error (ppm)                                |
| `delta_rt_min`          | RT difference between paired features (min)             |
| `annotation_confidence` | `High`, `Medium`, or `Low`                              |
| `adduct`                | Ion notation (e.g., `[M+DNPH−H2O+H]+`)                 |

---

## Confidence Levels

| Level    | Criterion                                                                 |
|----------|---------------------------------------------------------------------------|
| `High`   | Both Δm/z and ΔRT are within the innermost third of the tolerance window  |
| `Medium` | Both are within the middle third                                          |
| `Low`    | Either is in the outer third (but still within tolerance)                 |

---

## Usage – Compound Discoverer Integration

To integrate this node into Compound Discoverer's scripted processing node framework:

1. Add `derivatization_node.py` to your CD scripting path.
2. In your CD scripting node, call:

```python
from derivatization_node import run_derivatization_node

# `feature_table` is the list of feature dicts supplied by CD
results = run_derivatization_node(
    feature_table=feature_table,        # list of dicts from CD
    reaction_ids=["DNPH", "DANSYL"],    # reactions to apply
    mz_tolerance_ppm=5.0,
    rt_tolerance_min=0.3,
    max_derivatizations=3,
)
# `results` is a list of annotation dicts ready to inject back into CD
```

`run_derivatization_node()` automatically normalises common Compound Discoverer column-name variants (`ID`, `m/z`, `MolecularWeight`, `RetentionTime`, `PeakArea`, etc.).

---

## Running Tests

```bash
pip install pytest numpy pandas
python -m pytest tests/ -v
```

---

## Repository Structure

```
derivatization-node/
├── derivatization_node.py      # Main processing node script
├── requirements.txt            # Python dependencies
├── tests/
│   └── test_derivatization_node.py  # Unit tests
├── LICENSE                     # MIT License
└── README.md                   # This file
```

---

## Citation

If you use this software in your research, please cite:

> Rangel Huerta OD *et al.* *Derivatization-Driven Annotation via Functional Group-Specific Reactions Using Compound Discoverer Software.* *(manuscript in preparation)*

---

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
