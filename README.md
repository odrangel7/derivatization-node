# derivatization-node

> **Processing Node for Derivatization-Driven Annotation via Functional Group-Specific Reactions Using Compound Discoverer Software**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Compound%20Discoverer-blue)]()
[![Status](https://img.shields.io/badge/status-in%20development-orange)]()

---

## Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Node Configuration & Parameters](#node-configuration--parameters)
6. [Workflow Integration](#workflow-integration)
7. [Input & Output](#input--output)
8. [Usage Example](#usage-example)
9. [Derivatization Reactions Supported](#derivatization-reactions-supported)
10. [Troubleshooting](#troubleshooting)
11. [Contributing](#contributing)
12. [Citation](#citation)
13. [License](#license)
14. [Contact](#contact)

---

## Overview

<!-- TODO: Provide a detailed description of the node's purpose and scientific rationale.
     Example:
     This scripted processing node extends Thermo Fisher Scientific's Compound Discoverer (CD)
     software to enable derivatization-driven annotation of small molecules. By applying
     functional group-specific chemical derivatization reactions in silico, the node predicts
     expected mass shifts and generates candidate annotations for detected features, thereby
     improving coverage and confidence in untargeted metabolomics workflows.
-->

_Fill in a paragraph describing the scientific background, the problem this node solves, and how it fits into a Compound Discoverer workflow._

---

## Features

<!-- TODO: List the key capabilities of the node. Example bullet points:
- Supports multiple derivatization reagents (e.g., DNPH, dansyl chloride, TMS, …)
- Automatically calculates expected mass shifts for each reagent/functional-group pair
- Integrates directly into Compound Discoverer as a Scripted Node
- Outputs annotated compound tables with derivatization evidence scores
- Configurable reaction specificity and mass tolerance
-->

- [ ] Feature 1 — _describe it_
- [ ] Feature 2 — _describe it_
- [ ] Feature 3 — _describe it_

---

## Requirements

| Requirement | Version / Notes |
|---|---|
| Thermo Fisher Compound Discoverer | ≥ X.X <!-- TODO: specify minimum CD version --> |
| Python / R / other runtime | X.X <!-- TODO: specify language & version --> |
| Operating System | Windows <!-- TODO: confirm OS support --> |
| Additional libraries / packages | _list here_ |

<!-- TODO: Add any hardware requirements (RAM, disk space) if relevant. -->

---

## Installation

> **TODO:** Replace the steps below with the actual installation procedure.

### 1. Clone or download this repository

```bash
git clone https://github.com/odrangel7/derivatization-node.git
```

### 2. Install dependencies

```bash
# TODO: add the real command, e.g.:
# pip install -r requirements.txt
```

### 3. Register the Scripted Node in Compound Discoverer

1. Open Compound Discoverer.
2. Navigate to **Administration → Scripted Node Manager** (or the equivalent menu path).
3. Click **Add** and browse to the script file (`<!-- TODO: script filename -->`).
4. Configure the node parameters (see [Node Configuration & Parameters](#node-configuration--parameters)).
5. Save and restart CD if prompted.

---

## Node Configuration & Parameters

<!-- TODO: Describe every parameter exposed in the CD node panel. Use the table below as a template. -->

| Parameter | Type | Default | Description |
|---|---|---|---|
| `reagent` | String | — | Derivatization reagent name (e.g., `DNPH`) |
| `mass_tolerance_ppm` | Float | `5.0` | Mass tolerance in ppm for matching derivatized adducts |
| `min_score` | Float | `0.0` | Minimum annotation score threshold |
| _add more_ | | | |

---

## Workflow Integration

<!-- TODO: Insert a diagram or step-by-step description of where this node sits in a typical
     Compound Discoverer workflow (e.g., after feature detection, before gap filling). -->

```
[Data Import] → [Feature Detection] → [Derivatization Annotation Node] → [Compound Consolidation] → [Results]
                                              ↑
                                         (this node)
```

_Describe any upstream nodes whose output this node depends on, and any downstream nodes that consume its output._

---

## Input & Output

### Input

<!-- TODO: Describe the expected input format (e.g., mzXML, mzML, CD result file). -->

| File / Table | Format | Description |
|---|---|---|
| Raw data or CD result | `.raw` / `.cdResult` | Instrument data or previous CD result |
| _add more_ | | |

### Output

<!-- TODO: Describe all output tables, columns, and files produced by the node. -->

| Output | Format | Description |
|---|---|---|
| Annotated compound table | CD compound table | Includes new columns: `Derivatization_Reagent`, `Derivatization_Score`, … |
| _add more_ | | |

---

## Usage Example

<!-- TODO: Walk through a minimal end-to-end example. Include screenshots or command-line
     invocations if the node can also be run standalone. -->

```
# TODO: provide a concrete example
```

---

## Derivatization Reactions Supported

<!-- TODO: List each reagent/reaction supported, the functional groups it targets, and
     the expected mass shift(s). Example table: -->

| Reagent | Target Functional Group | Mass Shift (Da) | Reference |
|---|---|---|---|
| DNPH | Aldehydes & Ketones | +180.0424 | _cite source_ |
| Dansyl chloride | Amines & Phenols | +233.0567 | _cite source_ |
| _add more_ | | | |

---

## Troubleshooting

<!-- TODO: Document common errors and their solutions. -->

| Symptom | Likely Cause | Solution |
|---|---|---|
| Node fails to load in CD | Incompatible CD version | Check [Requirements](#requirements) |
| No annotations produced | Mass tolerance too strict | Increase `mass_tolerance_ppm` |
| _add more_ | | |

---

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository.
2. Create a feature branch: `git checkout -b feature/my-feature`.
3. Commit your changes: `git commit -m "Add my feature"`.
4. Push to your fork and open a Pull Request.

Please open an issue first to discuss significant changes.

---

## Citation

<!-- TODO: Add the citation once the work is published. Example:
Rangel Huerta, O.D. et al. (2026). "Derivatization-Driven Annotation via Functional
Group-Specific Reactions Using Compound Discoverer." *Journal Name*, vol(issue), pages.
DOI: https://doi.org/...
-->

_Citation information will be added upon publication._

---

## License

This project is licensed under the [MIT License](LICENSE) © 2026 Oscar Daniel Rangel Huerta.

---

## Contact

**Oscar Daniel Rangel Huerta**
<!-- TODO: add email, ORCID, institution, or other contact details -->

_GitHub: [@odrangel7](https://github.com/odrangel7)_
