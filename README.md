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

## Contact

For questions, suggestions, or bug reports, please open an issue in this repository or contact the corresponding author.
