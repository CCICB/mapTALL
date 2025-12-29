# mapTALL

Transcriptomic subtyping of T-cell acute lymphoblastic leukaemia (T-ALL) using scPred
and projection onto the St Jude reference dataset.

---

## Overview

**mapTALL** is an R-based workflow for assigning transcriptomic subtypes to bulk
RNA-seq T-ALL samples by projecting them into the feature space defined by the
St Jude reference cohort published by **Pölönen et al. (2024)** and applying
supervised classification using **scPred**.

Although scPred was originally developed for single-cell RNA-seq, mapTALL applies
the same framework to bulk RNA-seq by treating each sample as a “cell” in a
Seurat object.

The workflow is designed to be:
- transparent
- modular
- conservative in interpretation
- suitable for public reuse

mapTALL is provided as a script-based workflow rather than an R package.

---

## Reference data

### St Jude reference cohort

The reference dataset is derived from the publicly released transcriptomic data
accompanying **Pölönen et al. (2024)**. The full reference object is constructed
internally using the authors’ published subtype annotations and dimensional
reduction strategy.

For public distribution, a **lean reference object** is provided:

- contains RNA **counts** and **normalised data** for those genes
- retains the **PCA reduction** defining the scPred feature space
- includes essential annotation fields:
  - `Reviewed.subtype`
  - `Reviewed.genetic.subtype`
  - `meta_parent`

The lean reference object is created after adding a parent classification to the St Jude annotations, followed by removal of bulky components (e.g. scaled data,
neighbour graphs) for public distribution.

### Reference object

Due to file size constraints, the St Jude reference Seurat object used by
mapTALL is **not stored directly in this GitHub repository**.

The reference object is hosted on the Open Science Framework (OSF):

https://osf.io/adht8/

Users must download the reference file:

- `SJ_with_parent.lean.obj.rds`

and place it in a local `Data/` directory before running the mapTALL script.

---

## What the script does

The mapTALL run script performs the following steps:

### 1. Load inputs
- A pre-built **St Jude reference Seurat object**
- A user-provided test cohort (bulk RNA-seq counts)

### 2. Pre-processing
- Intersects test genes with the reference gene set
- Applies **ComBat-seq** batch adjustment to the shared gene set
- Constructs a Seurat object for the test cohort

### 3. scPred classification
Three scPred classifiers are applied sequentially:

1. **Reviewed.subtype**
   - Primary transcriptomic subtype classification
   - Applied to all samples

2. **Reviewed.genetic.subtype**
   - Applied only to samples whose predicted subtype belongs to a
     genetic-eligible family

3. **meta_parent**
   - Higher-level grouping of transcriptomic subtypes


---

## Genetic subtype logic

Genetic subtype prediction is conditional on subtype applicability.

### Genetic-eligible subtype families

By default, genetic subtype classification is only attempted for samples
predicted as:

- TAL1_AB-like  
- TAL1_DP-like  
- TLX3  
- ETP-like  
- NKX2-1  

This reflects the biological scope of the genetic subtype annotations in the
Pölönen dataset.

### Interpretation of outputs

| Value | Meaning |
|------|--------|
| `NA` | Genetic subtype is not applicable for this subtype family |
| `unassigned` | Genetic subtype is applicable, but confidence threshold was not met |
| Label | High-confidence genetic subtype assignment |

This distinction is intentional and preserved in the output tables.

---

## Outputs

The script produces the following key outputs:

### 1. Summary table  
`Results_mapTALL_summary.csv`

Per sample:
- transcriptomic subtype prediction, no-rejection class, and max probability
- genetic subtype prediction (or NA where not applicable)
- meta-parent classification
- risk category (if enabled)

### 2. Annotated Seurat object  
`test.scpred.obj.rds`

Contains all predictions and probabilities in `@meta.data`.

### 3. Metadata table  
`test.scpred.meta.csv`

Flat table of Seurat metadata for downstream analysis.

### 4. Session information  
A dated `sessionInfo_YYYY-MM-DD.txt` file for reproducibility.

---


## Requirements

- R ≥ 4.2
- Required packages:
  - Seurat (v5)
  - scPred
  - tidyverse
  - data.table
- Optional:
  - GSVA
  - PROGENy

---

## Citation

If you use mapTALL, please cite:

Pölönen et al. (2024)  
*The genomic basis of childhood T-lineage acute lymphoblastic leukaemia*

Alquicira-Hernandez et al. (2019)  
*scPred: accurate supervised method for cell-type classification from single-cell RNA-seq data*

### mapTALL workflow

If mapTALL or the parent-level classification (`meta_parent`) introduced in this
repository is used in analysis or publication, please cite this repository:

Corley SM. *mapTALL: Transcriptomic classification of T-ALL using scPred and the
St Jude reference cohort*. GitHub repository.  
https://github.com/CCICB/mapTALL

A manuscript describing mapTALL is in preparation.
---


## Contact

Questions, issues, or suggestions should be raised via GitHub Issues.