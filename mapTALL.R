# ============================================================
# mapTALL_run.R — test subtyping

#
# What users provide:
#   1) A raw counts matrix (genes x samples) with:
#        - rownames = gene symbols/IDs matching the SJ reference
#        - colnames = sample IDs using ONLY letters/numbers/_ (underscores)
#          (no spaces, no dots, no slashes). Example: "ZCC_001_TALL"
#   2) (Optional) a sample metadata CSV with a "sample_id" column that
#      matches the counts colnames exactly.
#
# What this script does:
#   1) Load SJ reference object (bundled with repo)
#   2) Extract reference counts directly from ref.obj 
#   3) Intersect genes between reference and test
#   4) Batch-correct SJ + test with ComBat_seq (counts-based)
#   5) Build test Seurat object and run scPred predictions:
#        - Reviewed.subtype
#        - Reviewed.genetic.subtype
#        - meta_parent
#   6) Compute ONE risk category per sample 
#        - subtype-only risk where appropriate (e.g. SPI1)
#        - genetic risk only used when subtype supports that family
#   7) Write outputs (tables + RDS)
# ============================================================

# ---------------------------
# 0) Libraries
# ---------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(sva)       # ComBat_seq
  library(scPred)
  library(Matrix)
  library(readr)
})

# ---------------------------
# 1) User settings 
# ---------------------------

# Reference object (ship in repo)
REF_RDS <- "Data/SJ_with_parent.obj.rds"

# Test counts (genes x samples), CSV or TSV
TEST_COUNTS_FILE <- "Input/Example.TALL.counts.csv"
TEST_COUNTS_FORMAT <- "csv"   # "csv" or "tsv"

# Optional metadata (must have sample_id column matching counts colnames)
TEST_META_FILE <- "NA" #"input/test_meta.csv"  # set to NA if not used
TEST_META_SAMPLE_COL <- "sample_id"

# Output directories
OUT_DIR <- "output"
OUT_OBJ <- file.path(OUT_DIR, "objects")
OUT_TAB <- file.path(OUT_DIR, "tables")

# Batch correction
DO_COMBAT_SEQ <- TRUE

# Seurat/UMAP settings (reasonable defaults; not critical for scPred)
NPCS <- 50
DIMS <- 1:50
UMAP_SEED <- 42

# scPred settings (from your current pipeline)
SCPRED_SIG <- 1
SCPRED_THRESHOLD <- 0.55
SCPRED_SEED <- 44

BAD_TRUTH_LABELS <- c("", "noAlloc", "unassigned", "NA", "NA.", "None/Other")



# Genetic subtype prediction behaviour
# - Default (recommended): Only run genetic-subtype scPred for subtype families where
#   genetic.subtype labels are intended to be interpreted.
#   This ensures:
#     NA        = not applicable for that predicted subtype family
#     unassigned = applicable, but confidence threshold not reached
#
# - Optional: also run the genetic-subtype model across ALL samples and export
#   the top class/probability as exploratory-only fields (clearly labelled).
RUN_GENETIC_ALL <- FALSE

# Subtype families for which a genetic.subtype prediction is considered applicable
# (edit this list if you expand your reference annotations)
GENETIC_ELIGIBLE_SUBTYPES <- c("TAL1_AB-like", "TAL1_DP-like", "TLX3", "ETP-like", "NKX2-1")
# ---------------------------
# 2) Risk mapping (simple + guarded)
# ---------------------------

risk_categories <- list(
  "Very High Risk" = c("KMT2A_ETP-like", "SPI1", "LMO2_GD-like"),
  "High Risk"      = c("MLLT10_ETP-like", "NUP98_ETP-like", "Rare_ETP-like", "TAL1_AB-like_Other", "NKX2-5"),
  "Low Risk"       = c("HOXA13_ETP-like", "MED12_ETP-like", "ETV6_ETP-like", "NUP214_ETP-like",
                       "TAL1_DP-like_JAK", "TAL1_DP-like_LEF1/LYL1", "TAL1_DP-like_Other",
                       "TAL1_AB-like_Notch_wt", "TAL1_AB-like_Loss_6q", "TLX3_Immature",
                       "BCL11B", "MLLT10", "STAG2_LMO2", "NUP98", "TME-enriched"),
  "Very Low Risk"  = c("NKX2-1", "TLX1", "TLX3_DP-like", "KMT2A", "NUP214", "HOXA9_TCR",
                       "TAL1_DP-like_RPL10", "ZFP36L2_ETP-like")
)

# Subtypes whose risk is defined by subtype alone (do not require genetic subtype)
subtype_driven <- c("SPI1", "BCL11B", "NKX2-1", "TLX1", "KMT2A", "MLLT10", "STAG2_LMO2", "TME-enriched", "LMO2_GD-like")

# Families where risk is defined by genetic.subtype, but ONLY if subtype supports it
genetic_driven_families <- list(
  "TAL1_AB-like" = c("TAL1_AB-like_Notch_wt", "TAL1_AB-like_Loss_6q", "TAL1_AB-like_Other"),
  "TAL1_DP-like" = c("TAL1_DP-like_JAK", "TAL1_DP-like_LEF1/LYL1", "TAL1_DP-like_Other", "TAL1_DP-like_RPL10"),
  "ETP-like"     = c("KMT2A_ETP-like", "MLLT10_ETP-like", "NUP98_ETP-like", "Rare_ETP-like",
                     "HOXA13_ETP-like", "MED12_ETP-like", "ETV6_ETP-like", "NUP214_ETP-like", "ZFP36L2_ETP-like"),
  "TLX3"         = c("TLX3_Immature", "TLX3_DP-like")
)

risk_of_label <- function(label) {
  if (is.na(label) || label == "unassigned") return(NA_character_)
  for (rk in names(risk_categories)) {
    if (label %in% risk_categories[[rk]]) return(rk)
  }
  NA_character_
}

categorize_risk <- function(subtype_pred, genetic_pred) {
  # 1) If subtype itself is a subtype-driven risk label, use subtype only.
  if (!is.na(subtype_pred) && subtype_pred %in% subtype_driven) {
    return(risk_of_label(subtype_pred))
  }

  # 2) If subtype is unassigned, do not assign risk based only on genetic subtype.
  if (is.na(subtype_pred) || subtype_pred == "unassigned") {
    return(NA_character_)
  }

  # 3) If subtype is a genetic-driven family, only accept genetic labels within that family.
  if (subtype_pred %in% names(genetic_driven_families)) {
    if (!is.na(genetic_pred) && genetic_pred %in% genetic_driven_families[[subtype_pred]]) {
      return(risk_of_label(genetic_pred))
    }
    # fall back to subtype risk if defined
    return(risk_of_label(subtype_pred))
  }

  # 4) Fallback: subtype risk if defined
  risk_of_label(subtype_pred)
}

# ---------------------------
# 3) Helpers
# ---------------------------

dir_create <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)

read_counts <- function(path, format = c("csv", "tsv")) {
  format <- match.arg(format)
  if (!file.exists(path)) stop("Counts file not found: ", path)
  if (format == "csv") {
    x <- read.csv(path, header = TRUE, row.names = 1, check.names = FALSE)
  } else {
    x <- read.delim(path, header = TRUE, row.names = 1, check.names = FALSE)
  }
  as.matrix(x)
}

as_integer_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[is.na(x)] <- 0
  x <- round(x)
  storage.mode(x) <- "integer"
  x
}

assert_has_cols <- function(df, cols, df_name = "data.frame") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) stop(df_name, " missing required column(s): ", paste(missing, collapse = ", "))
}

check_sample_ids <- function(ids) {
  # enforce simple sample IDs: letters, numbers, underscore only
  bad <- ids[!grepl("^[A-Za-z0-9_]+$", ids)]
  if (length(bad) > 0) {
    stop(
      "Sample IDs must contain only letters/numbers/underscore.\n",
      "Bad IDs (first 10 shown): ", paste(head(bad, 10), collapse = ", "), "\n",
      "Please fix your counts colnames (and metadata sample_id) to match this rule."
    )
  }
  invisible(TRUE)
}

norm_label <- function(x) {
  ifelse(is.na(x), NA_character_, trimws(x))
}


# Seurat v5 compatibility 

`%||%` <- function(x, y) if (!is.null(x)) x else y

GetAssayData <- function(object, assay = NULL, slot = NULL, layer = NULL, ...) {
  # If "data"/"counts"/"scale.data" is passed as the 2nd positional arg,
  # treat it as a layer/slot name (old Seurat v3/v4 behaviour).
  if (!is.null(assay) && assay %in% c("data", "counts", "scale.data") &&
      is.null(slot) && is.null(layer)) {
    layer <- assay
    assay <- Seurat::DefaultAssay(object)
  }
  
  # Seurat v5 prefers 'layer'; SeuratObject supports 'slot' for compatibility.
  if (!is.null(layer) && is.null(slot)) slot <- layer
  
  SeuratObject::GetAssayData(
    object = object,
    assay  = assay %||% Seurat::DefaultAssay(object),
    slot   = slot,
    ...
  )
}

# ---------------------------
# 4) Create output dirs
# ---------------------------
dir_create(OUT_DIR); dir_create(OUT_OBJ); dir_create(OUT_TAB)

# ---------------------------
# 5) Load reference and extract counts
# ---------------------------
stopifnot(file.exists(REF_RDS))
ref.obj <- readRDS(REF_RDS)
DefaultAssay(ref.obj) <- "RNA"
ref.counts <- GetAssayData(ref.obj, slot = "counts")
ref.counts <- as(ref.counts, "dgCMatrix")

message("Reference loaded: ", REF_RDS)
message("Reference samples: ", ncol(ref.obj), " | Reference genes: ", nrow(ref.obj))

# ---------------------------
# 6) Load test counts + metadata
# ---------------------------
test.counts <- read_counts(TEST_COUNTS_FILE, TEST_COUNTS_FORMAT)
check_sample_ids(colnames(test.counts))

test.meta <- NULL
if (!is.na(TEST_META_FILE) && file.exists(TEST_META_FILE)) {
  test.meta <- read.csv(TEST_META_FILE, header = TRUE, check.names = FALSE)
  assert_has_cols(test.meta, TEST_META_SAMPLE_COL, "test.meta")
  test.meta[[TEST_META_SAMPLE_COL]] <- as.character(test.meta[[TEST_META_SAMPLE_COL]])
  check_sample_ids(test.meta[[TEST_META_SAMPLE_COL]])
  rownames(test.meta) <- test.meta[[TEST_META_SAMPLE_COL]]
  message("Test metadata loaded: ", TEST_META_FILE)
}

# Align counts to metadata if provided
if (!is.null(test.meta)) {
  keep <- intersect(colnames(test.counts), rownames(test.meta))
  test.counts <- test.counts[, keep, drop = FALSE]
  test.meta <- test.meta[keep, , drop = FALSE]
}

# ---------------------------
# 7) Harmonise genes and ComBat-seq
# ---------------------------
common_genes <- intersect(rownames(ref.counts), rownames(test.counts))
if (length(common_genes) < 5000) {
  warning("Only ", length(common_genes), " common genes. Check gene identifiers.")
}
message("Common genes used: ", length(common_genes))

ref.ed  <- ref.counts[common_genes, , drop = FALSE]
test.ed <- test.counts[common_genes, , drop = FALSE]

combined <- cbind(ref.ed, test.ed)
combined <- as_integer_matrix(combined)

batch.labels <- c(rep("SJ", ncol(ref.ed)), rep("TEST", ncol(test.ed)))

if (DO_COMBAT_SEQ) {
  message("Running ComBat_seq on SJ + TEST.")
  adj <- ComBat_seq(combined, batch = batch.labels)
  adj <- as_integer_matrix(adj)
} else {
  message("Skipping batch correction (DO_COMBAT_SEQ = FALSE).")
  adj <- combined
}

adj_test <- adj[, (ncol(ref.ed) + 1):ncol(adj), drop = FALSE]

saveRDS(adj_test, file.path(OUT_OBJ, "adj_test_counts.rds"))

# ---------------------------
# 8) Build test Seurat object
# ---------------------------
test.obj <- CreateSeuratObject(counts = adj_test, project = "mapTALL", min.cells = 0, min.features = 0)

if (!is.null(test.meta)) {
  test.obj <- AddMetaData(test.obj, test.meta)
}

test.obj <- NormalizeData(test.obj) %>% ScaleData()
VariableFeatures(test.obj) <- VariableFeatures(ref.obj)

test.obj <- RunPCA(test.obj, npcs = NPCS, seed.use = UMAP_SEED) %>%
  RunUMAP(dims = DIMS, seed.use = UMAP_SEED)

saveRDS(test.obj, file.path(OUT_OBJ, "test.obj.rds"))
write_csv(test.obj@meta.data %>% tibble::rownames_to_column("sample_id"),
          file.path(OUT_TAB, "test.obj.meta.csv"))

# ---------------------------
# 9) scPred predictions
# ---------------------------
threshold <- SCPRED_THRESHOLD
seed <- SCPRED_SEED
sig <- SCPRED_SIG

DefaultAssay(ref.obj)  <- "RNA"
DefaultAssay(test.obj) <- "RNA"

ref.filt.obj <- subset(ref.obj, subset = !(Reviewed.genetic.subtype %in% BAD_TRUTH_LABELS))

run_one <- function(ref_train, test_obj, truth_col, prefix) {
  ref_fs <- getFeatureSpace(ref_train, truth_col, sig = sig, reduction = "pca")
  ref_fs <- trainModel(ref_fs, model = "svmRadial", seed = seed)
  out <- scPredict(test_obj, ref_fs, threshold = threshold, seed = seed)@meta.data

  out <- out %>%
    rownames_to_column("sample_id") %>%
    rename(
      !!paste0(prefix, "_scpred_max")          := scpred_max,
      !!paste0(prefix, "_scpred_prediction")   := scpred_prediction,
      !!paste0(prefix, "_scpred_no_rejection") := scpred_no_rejection
    )
  out
}

# --- scPred/Seurat v5 workaround: scPred expects an assay named "data"
if (!("data" %in% names(Seurat::Assays(ref.obj)))) {
  ref.obj[["data"]] <- ref.obj[["RNA"]]
}
if (!("data" %in% names(Seurat::Assays(test.obj)))) {
  test.obj[["data"]] <- test.obj[["RNA"]]
}

# 1) SUBTYPE first (needed to decide if genetic-subtype prediction is applicable)
sub_md <- run_one(ref.obj, test.obj, "Reviewed.subtype", "subtype")

# Attach subtype predictions to test.obj so we can subset for genetic prediction
sub_md_rn <- sub_md %>% column_to_rownames("sample_id")
test.obj  <- AddMetaData(test.obj, sub_md_rn[["subtype_scpred_prediction"]],
                         col.name = "subtype_scpred_prediction")

# 2) GENETIC subtype prediction (applicable-only, recommended default)
test_gen_eligible <- subset(
  test.obj,
  subset = subtype_scpred_prediction %in% GENETIC_ELIGIBLE_SUBTYPES
)

if (nrow(test_gen_eligible@meta.data) > 0) {
  gen_md <- run_one(ref.filt.obj, test_gen_eligible, "Reviewed.genetic.subtype", "genetic")
} else {
  gen_md <- tibble::tibble(
    sample_id = character(),
    genetic_scpred_max = numeric(),
    genetic_scpred_prediction = character(),
    genetic_scpred_no_rejection = character()
  )
}

# 2b) Optional: run the genetic-subtype model on ALL samples (exploratory only)
# These fields are clearly labelled and should NOT be interpreted as a genetic subtype
# assignment unless the sample is within GENETIC_ELIGIBLE_SUBTYPES.
if (isTRUE(RUN_GENETIC_ALL)) {
  gen_md_all <- run_one(ref.filt.obj, test.obj, "Reviewed.genetic.subtype", "genetic_model_all")
} else {
  gen_md_all <- tibble::tibble(sample_id = character())
}

# 3) PARENT grouping on all samples
par_md <- run_one(ref.obj, test.obj, "meta_parent", "parent")
preds <- sub_md %>%
  left_join(gen_md, by = "sample_id") %>%
  left_join(gen_md_all, by = "sample_id") %>%
  left_join(par_md, by = "sample_id")

write_csv(preds, file.path(OUT_TAB, "Results_all_predictions.csv"))

# Attach prediction columns to test.obj
to_add <- preds %>% column_to_rownames("sample_id")
for (cn in colnames(to_add)) {
  test.obj <- AddMetaData(test.obj, to_add[[cn]], col.name = cn)
}

# ---------------------------
# 10) Risk category (single column)
# ---------------------------
test.obj$subtype_scpred_prediction <- norm_label(test.obj$subtype_scpred_prediction)
test.obj$genetic_scpred_prediction <- norm_label(test.obj$genetic_scpred_prediction)

test.obj$risk_category <- mapply(
  categorize_risk,
  test.obj$subtype_scpred_prediction,
  test.obj$genetic_scpred_prediction
)

# ---------------------------
# 11) Write compact summary table
# ---------------------------

preds <- readr::read_csv(
  file.path(OUT_TAB, "Results_all_predictions.csv"),
  show_col_types = FALSE
)

risk_df <- test.obj@meta.data %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::select(sample_id, risk_category)

preds2 <- preds %>%
  dplyr::left_join(risk_df, by = "sample_id")

# ---------------------------
# Enforce NA for non-applicable genetic subtype predictions
# (NA = not applicable; "unassigned" = applicable but below confidence threshold)
# ---------------------------

GENETIC_ELIGIBLE_SUBTYPES <- c("TAL1_AB-like", "TAL1_DP-like", "TLX3", "ETP-like", "NKX2-1")

is_genetic_applicable <- (preds2$subtype_scpred_prediction %in% GENETIC_ELIGIBLE_SUBTYPES) |
  (preds2$subtype_scpred_no_rejection %in% GENETIC_ELIGIBLE_SUBTYPES)

preds2 <- preds2 %>%
  dplyr::mutate(
    genetic_scpred_prediction   = dplyr::if_else(is_genetic_applicable, genetic_scpred_prediction,   NA_character_),
    genetic_scpred_no_rejection = dplyr::if_else(is_genetic_applicable, genetic_scpred_no_rejection, NA_character_),
    genetic_scpred_max          = dplyr::if_else(is_genetic_applicable, genetic_scpred_max,          NA_real_)
  )


summary_tab <- preds2 %>%
  dplyr::select(
    sample_id,
    
    subtype_prediction = subtype_scpred_prediction,
    subtype_noRej      = subtype_scpred_no_rejection,
    subtype_max        = subtype_scpred_max,
    
    genetic_prediction = genetic_scpred_prediction,
    genetic_noRej      = genetic_scpred_no_rejection,
    genetic_max        = genetic_scpred_max,
    
    parent_prediction  = parent_scpred_prediction,
    parent_noRej       = parent_scpred_no_rejection,
    parent_max         = parent_scpred_max,
    
    risk_category
  )

readr::write_csv(summary_tab, file.path(OUT_TAB, "Results_mapTALL_summary.csv"))
# Final outputs
saveRDS(test.obj, file.path(OUT_OBJ, "test.scpred.obj.rds"))
write_csv(test.obj@meta.data %>% tibble::rownames_to_column("sample_id"),
          file.path(OUT_TAB, "test.scpred.meta.csv"))

message("Risk category counts:")
print(table(test.obj$risk_category, useNA = "ifany"))

session_file <- file.path(OUT_DIR, paste0("sessionInfo_", format(Sys.Date(), "%Y-%m-%d"), ".txt"))
sink(session_file)
cat("Session info written on:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")
print(sessionInfo())
sink()

message("✅ mapTALL_run_simple complete. See: ", OUT_DIR)