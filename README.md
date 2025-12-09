# mapTALL
### Transcriptomic Subtyping of T-ALL by Projection onto the St. Jude Reference Using scPred

**mapTALL** is an R-based workflow for assigning transcriptomic subtypes to T-cell acute lymphoblastic leukaemia (T-ALL) samples.  
It projects local cohort RNA-seq data onto the large St. Jude (SJ) reference dataset published by **Pölönen et al. (2024)** and performs subtype classification using **scPred**.

Although scPred was originally designed for single-cell data, mapTALL applies the same supervised-learning approach to **bulk RNA-seq**, treating each bulk sample as a “cell” within a Seurat object.

mapTALL is a **transparent, modular analysis workflow**, not a standalone R package.

---

## 🔍 Key Features

- Uses the curated **SJ T-ALL reference (N = 1335)** with Reviewed subtype and genetic subtype labels  
- Supervised classification via **scPred** using:
  - 300 variable genes defined by Pölönen  
  - 50 principal components  
  - Probability-based class assignment
- Multi-level prediction:
  - **Reviewed.subtype** (e.g., TAL1-DP-like, LMO2/LYL1, TLX1, TLX3, ETP)
  - **Reviewed.genetic.subtype** (e.g., STIL-TAL1, HOXA-MLL, TLX3-RAG)
  - **meta_parent** grouping (TAL1, TLX, NKX2, HOXA_MLL, ETP, etc.)
- Built-in confidence metrics:
  - `scpred_prediction`
  - `scpred_no_rejection`
  - `scpred_max` (max probability)
- QC module categorises samples as:
  - **High confidence**
  - **Medium confidence**
  - **Discordant / Ambiguous**
- Optional downstream analyses:
  - TAL1 **DP–AB axis** scoring  
  - **GSVA** and **PROGENy** pathway activity  
  - Drug-response signatures (e.g., dasatinib sensitivity)

---

## 📦 Requirements

- **R ≥ 4.2**
- Required packages:
  - `Seurat`
  - `scPred`
  - `tidyverse`
  - `data.table`
  - `ggplot2`, `patchwork`
- Optional:
  - `GSVA`
  - `progeny`

**Dataset availability:**  
The St. Jude reference dataset is *not* included and must be downloaded from the Pölönen et al. 2024 repository under the appropriate data access conditions.
