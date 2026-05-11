# Derivatization-node

**Processing Node for Derivatization-Driven Annotation via Functional Group-Specific Reactions Using Compound Discoverer Software**

## Overview
The **Expected Compounds 2.0** node is an R-based Scripting Node designed to integrate seamlessly with Thermo Fisher's Compound Discoverer software. It facilitates the discovery and annotation of derivatized compounds by detecting mass shifts associated with specific functional group reactions (e.g., Thiols, DiThiols, or Custom user-defined mass shifts).

By comparing potential products against precursors found in control samples, the node effectively filters out noise and highlights chemically and biologically relevant derivatization pairs.

## Key Features
* **Automated Mass Shift Detection:** Identifies potential products of derivatization reactions based on specific, user-defined mass shifts.
* **Control Sample Filtering:** Uses Control samples to locate original precursor compounds, ensuring that matched products are derived from actual sample components rather than background noise.
* **Precursor-Product Pair Linking:** Automatically generates paired links (e.g., `Thiol_1`) between the suspected precursor and the resulting derivatized product.
* **Flexible Configuration:** Easily customize mass shifts and mass matching tolerances directly from the Compound Discoverer user interface.

## Prerequisites
* **Compound Discoverer** (Thermo Fisher) installed.
* **R (v4.0 or higher)** installed on the machine running Compound Discoverer.
* Required R packages: `data.table`, `RJSONIO`, `dplyr`, `stringr`, `tidyr`. The script will attempt to install these automatically if they are missing.

## Installation & Setup
1. Clone or download this repository to your local machine.
2. Open the `node.json` file in a text editor.
3. Update the `ExecutablePath` to point to your local `Rscript.exe` (e.g., `"C:\\Program Files\\R\\R-4.5.2\\bin\\Rscript.exe"`).
4. Update the `ExecutableCommandLineArguments` to point to the absolute path of the `ExpectedCompounds_Node_Paper.R` script on your machine.
5. Register the `node.json` file in Compound Discoverer as a custom Scripting Node.

## Data for testing
Data presented in the manuscript can be downloaded at **insert the link when published**

## Usage
1. Add the **"Expected Compounds 2.0"** node to your Compound Discoverer workflow as a Post-Processing node.
2. Ensure it is connected to upstream nodes that output the **Compounds**, **Compounds per File**, and **Input Files** tables.
3. Select the node in your workflow to configure its parameters:
   * **Tolerance:** Mass tolerance in Daltons (Da) for matching precursors to products. (Default: `0.005`)
   * **Thiol Mass:** Mass shift for Thiol reactions. Set to `0` to disable.
   * **DiThiol Mass:** Mass shift for DiThiol reactions. Set to `0` to disable.
   * **Custom 1:** Define any custom mass shift for alternative derivatization reactions. Set to `0` to disable.

## Output
Upon successful execution, the node will append new data to your Compound Discoverer result tables:
* **Boolean Flags:** New columns indicating whether a compound is a product (e.g., `Thiol = True`) or a precursor (e.g., `Suspected_precursor_Thiol = True`).
* **Pairs Column:** A `Pairs` column detailing specific Precursor-Product linkages, allowing for easy sorting and filtering of successfully derivatized components.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
