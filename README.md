# mapTALL

Transcriptomic subtyping of T-cell acute lymphoblastic leukaemia (T-ALL) using scPred
and projection onto the St Jude reference dataset.

---

## Overview

**mapTALL** is an R-based workflow for assigning transcriptomic subtypes to bulk
RNA-seq T-ALL samples by projecting them onto the St Jude reference cohort
published by **Pölönen et al. (2024)** and applying **scPred**.

Although scPred was originally developed for single-cell RNA-seq, mapTALL applies
the same supervised classification framework to bulk RNA-seq by treating each
sample as a “cell” in a Seurat object.

The workflow is designed to be:
- transparent
- modular
- conservative in interpretation
- suitable for public reuse

mapTALL is a script-based workflow rather than a packaged R library.

---

## What the script does

The GitHub script performs the following steps:

### 1. Load inputs
- A pre-built St Jude Seurat reference object with embedded scPred models
- A user-provided Seurat object containing test cohort RNA-seq data

### 2. scPred classification
Three independent scPred models are applied:

1. **Reviewed.subtype**
   - Primary transcriptomic subtype assignment
   - Applied to all samples

2. **Reviewed.genetic.subtype**
   - Applied **only** to samples whose predicted subtype belongs to a
     genetic-eligible family
   - Genetic subtypes are not defined for all transcriptomic families

3. **meta_parent**
   - Higher-level grouping of subtypes
   - Applied to all samples

---

## Genetic subtype logic (important)

Genetic subtype prediction is conditional on subtype applicability.

### Genetic-eligible subtype families

By default, genetic subtype scPred is only applied to samples predicted as:

- TAL1_AB-like
- TAL1_DP-like
- TLX3
- ETP-like
- NKX2-1

This behaviour mirrors the biological scope of the genetic subtype labels in the
Pölönen dataset.

### Interpretation of outputs

| Value | Meaning |
|------|--------|
| NA | Genetic subtype is not applicable for this subtype family |
| unassigned | Genetic subtype is applicable, but confidence threshold was not met |
| Label | High-confidence genetic subtype assignment |

This distinction is intentional and preserved in the final output tables.

---

## Outputs

The script writes the following key outputs:

### 1. Summary table
`Results_mapTALL_summary.csv`

Contains, per sample:
- subtype prediction, no-rejection class, and max probability
- genetic subtype prediction (or NA if not applicable)
- meta_parent prediction
- risk category (if enabled)

### 2. Annotated Seurat object
`test.scpred.obj.rds`

Includes all predictions and probabilities stored in `@meta.data`.

### 3. Metadata table
`test.scpred.meta.csv`

A flat table of Seurat metadata for downstream analysis.

### 4. Session information
A dated `sessionInfo_YYYY-MM-DD.txt` file for reproducibility.

---

## Risk categories

Risk categories are derived from subtype and, where applicable, genetic subtype
predictions. Genetic subtypes are only used when biologically meaningful for the
predicted subtype family.

---

## Requirements

- R ≥ 4.2
- Required packages:
  - Seurat
  - scPred
  - tidyverse
  - data.table
- Optional:
  - GSVA
  - PROGENy

The St Jude reference dataset is not included and must be obtained from the
original publication under appropriate data access conditions.

---

## Citation

If you use mapTALL, please cite:

Pölönen et al. (2024)  
Integrative multiomic classification of childhood T-ALL

Alquicira-Hernandez et al. (2019)  
scPred: accurate supervised cell-type classification

Please also cite this repository if the workflow is used in published work.

---

## Notes for users

- Genetic subtype predictions should not be interpreted independently of
  transcriptomic subtype.
- NA values for genetic subtype are expected and meaningful.
- This workflow prioritises conservative interpretation over maximal output.

---

## Contact

Questions, issues, or suggestions should be raised via GitHub Issues.