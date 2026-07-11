# DASS-anoikis

**DASS-anoikis** provides code and workflows for analyzing anoikis-state heterogeneity and anoikis-associated chemoresistance across cancer types. This repository accompanies our study on pan-cancer anoikis dynamics, single-cell resolution of ovarian cancer chemoresistance, FN1-mediated microenvironmental remodeling, and the development of **DASS**: **Deep-learning Anoikis Scoring System**.

## Overview

Anoikis is a detachment-induced form of apoptosis that plays important roles in tumor progression, metastasis, and therapeutic resistance. However, how heterogeneous anoikis states contribute to chemoresistance across cancer types remains incompletely understood.

In this project, we profiled anoikis-related pathway activity in **9,310 tumor samples across 30 solid cancer types** and identified six anoikis-associated subgroups spanning a continuum from anoikis-active to anoikis-repressed states. These subgroups showed distinct prognostic outcomes, immune microenvironment features, and multi-omic alterations.

Single-cell RNA-seq analysis of ovarian cancer further resolved this anoikis continuum at cellular resolution and revealed a non-linear association between anoikis activity and chemotherapy response. Both anoikis-activated and anoikis-repressed states were enriched in chemoresistant populations, whereas intermediate anoikis states were associated with chemosensitivity. Chemoresistant states were characterized by extracellular matrix remodeling and stromal-immune crosstalk centered on **FN1 signaling**.

This repository includes code for reproducing key analyses and applying **DASS** to quantify anoikis dynamics from transcriptomic data.

## Key Features

* Pan-cancer anoikis pathway activity profiling
* Identification of anoikis-associated molecular subgroups
* Survival, immune microenvironment, and multi-omic characterization
* Single-cell mapping of anoikis-state heterogeneity in ovarian cancer
* Chemotherapy response-associated anoikis state analysis
* FN1-centered stromal-immune communication analysis
* Deep learning-based anoikis scoring using DASS
* Example workflows for applying DASS to bulk or single-cell transcriptomic data

## Repository Structure

```text
DASS-anoikis/
├── data/                     # Example input files or metadata templates
├── scripts/                  # Main analysis scripts
│   ├── 01_preprocessing/
│   ├── 02_anoikis_scoring/
│   ├── 03_pancancer_subgroups/
│   ├── 04_survival_analysis/
│   ├── 05_immune_analysis/
│   ├── 06_single_cell_analysis/
│   ├── 07_cell_communication/
│   └── 08_DASS_model/
├── models/                   # Trained DASS model files
├── results/                  # Output tables and figures
├── notebooks/                # Example notebooks
├── docs/                     # Additional documentation
├── requirements.txt          # Python dependencies
├── environment.yml           # Conda environment file
└── README.md
```

## Installation

Clone this repository:

```bash
git clone https://github.com/your-username/DASS-anoikis.git
cd DASS-anoikis
```

Create the required environment:

```bash
conda env create -f environment.yml
conda activate dass-anoikis
```

Alternatively, install Python dependencies manually:

```bash
pip install -r requirements.txt
```

## Input Data

DASS can be applied to transcriptomic expression matrices from bulk RNA-seq, microarray, or single-cell RNA-seq data.

The expected input format is a gene expression matrix:

```text
Gene      Sample1    Sample2    Sample3
FN1       12.31      9.84       15.02
BCL2      7.42       8.11       6.93
CASP3     5.33       4.98       6.21
...
```

For single-cell data, the input can be a gene-by-cell or cell-by-gene matrix, depending on the selected script. Please check the corresponding example script before running the analysis.

## Quick Start

### 1. Calculate anoikis pathway activity

```bash
python scripts/02_anoikis_scoring/calculate_anoikis_activity.py \
  --input data/example_expression_matrix.tsv \
  --output results/anoikis_activity_scores.tsv
```

### 2. Apply DASS to quantify anoikis dynamics

```bash
python scripts/08_DASS_model/run_DASS.py \
  --input data/example_expression_matrix.tsv \
  --model models/DASS_model.pt \
  --output results/DASS_scores.tsv
```

### 3. Perform subgroup classification

```bash
Rscript scripts/03_pancancer_subgroups/classify_anoikis_subgroups.R \
  --input results/DASS_scores.tsv \
  --output results/anoikis_subgroups.tsv
```

## Main Analyses

### Pan-cancer anoikis heterogeneity

The pan-cancer analysis identifies anoikis-associated subgroups across 30 solid cancer types and evaluates their clinical, immune, and molecular characteristics.

Main steps include:

1. Anoikis pathway activity estimation
2. Consensus clustering or subgroup classification
3. Survival analysis
4. Immune infiltration analysis
5. Multi-omic feature comparison

### Single-cell anoikis continuum

The single-cell workflow maps anoikis-state heterogeneity at cellular resolution, especially in ovarian cancer. It can be used to examine anoikis activity across malignant cells, stromal cells, and immune cell populations.

Main steps include:

1. Single-cell quality control and preprocessing
2. Cell type annotation
3. Anoikis score calculation
4. Chemotherapy response-associated state analysis
5. Visualization of anoikis-active, intermediate, and anoikis-repressed states

### FN1-centered microenvironmental remodeling

The FN1-related analysis investigates stromal-immune communication and extracellular matrix remodeling associated with chemoresistant anoikis states.

Main steps include:

1. FN1 expression and pathway activity analysis
2. Cell-cell communication analysis
3. ECM remodeling signature analysis
4. Association with chemotherapy response

## DASS: Deep-learning Anoikis Scoring System

DASS is designed to quantify anoikis dynamics from transcriptomic data. It provides a sample-level or cell-level anoikis score that can be used to investigate anoikis-state transitions, tumor microenvironment remodeling, and therapy response.

The DASS output includes:

```text
Sample      DASS_score      Anoikis_state
Sample1     0.82            Anoikis-active
Sample2     0.45            Intermediate
Sample3     0.13            Anoikis-repressed
```

## Citation

If you use this repository or DASS in your research, please cite:

```text
[Author names]. Anoikis-state heterogeneity shapes chemoresistance across cancer types through FN1-mediated microenvironmental remodeling. [Journal/Preprint], [Year].
```

## Contact

For questions, suggestions, or bug reports, please open an issue in this repository or contact the corresponding author.

## License

This repository is released for academic and research use. Please refer to the `LICENSE` file for details.
