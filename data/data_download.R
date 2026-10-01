# ============================================================
#  TCGA UCEC Data Download & CSV Export
#  Uses: TCGAbiolinks (Bioconductor)
#  Output: CSV files ready for Python/pandas
# ============================================================

# ── 0. Install / load packages ───────────────────────────────
library(TCGAbiolinks)
library(SummarizedExperiment)

# All CSVs will land here — change to your preferred path
OUTPUT_DIR <- "./data/raw/TCGA_UCEC_data"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# TCGAbiolinks' own download cache (raw GDC files before GDCprepare
# processes them) — nested under data/raw/ so it's not a separate
# top-level folder, and .gitignore already excludes data/raw/ wholesale
GDC_CACHE_DIR <- "./data/raw/GDCdata"
dir.create(GDC_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Output directory:", OUTPUT_DIR, "\n")
cat("GDC cache directory:", GDC_CACHE_DIR, "\n")


# ============================================================
#  1. RNA-seq (STAR counts)
#     Gene × Sample count matrix — primary model input
# ============================================================
cat("\n[1/5] Querying RNA-seq...\n")

query_rna <- GDCquery(
  project           = "TCGA-UCEC",
  data.category     = "Transcriptome Profiling",
  data.type         = "Gene Expression Quantification",
  workflow.type     = "STAR - Counts",
  sample.type       = "Primary Tumor"   # exclude normals
)

GDCdownload(query_rna, method = "api", files.per.chunk = 20, directory = GDC_CACHE_DIR)

rna_se <- GDCprepare(query_rna, directory = GDC_CACHE_DIR)   # returns SummarizedExperiment

# --- Extract count matrix (genes × samples) ---
count_mat <- assay(rna_se, "unstranded")   # raw STAR counts

# Convert to data.frame with gene_id as first column
count_df <- as.data.frame(count_mat)
count_df <- cbind(gene_id = rownames(count_df), count_df)

write.csv(count_df,
          file      = file.path(OUTPUT_DIR, "ucec_rnaseq_counts.csv"),
          row.names = FALSE,
          quote     = FALSE)

cat("  -> Saved ucec_rnaseq_counts.csv\n",
    "     Shape:", nrow(count_df), "genes x", ncol(count_df) - 1, "samples\n")

# --- Also save TPM (normalized) if available ---
if ("tpm_unstrand" %in% assayNames(rna_se)) {
  tpm_mat <- assay(rna_se, "tpm_unstrand")
  tpm_df  <- as.data.frame(tpm_mat)
  tpm_df  <- cbind(gene_id = rownames(tpm_df), tpm_df)
  write.csv(tpm_df,
            file      = file.path(OUTPUT_DIR, "ucec_rnaseq_tpm.csv"),
            row.names = FALSE,
            quote     = FALSE)
  cat("  -> Saved ucec_rnaseq_tpm.csv\n")
}

# --- Gene metadata (Ensembl ID → gene name, biotype, etc.) ---
gene_meta <- as.data.frame(rowData(rna_se))
gene_meta <- cbind(gene_id = rownames(gene_meta), gene_meta)
write.csv(gene_meta,
          file      = file.path(OUTPUT_DIR, "ucec_gene_metadata.csv"),
          row.names = FALSE,
          quote     = FALSE)
cat("  -> Saved ucec_gene_metadata.csv\n")

# --- Sample metadata (barcodes, submitter IDs, etc.) ---
sample_meta_rna <- as.data.frame(colData(rna_se))
sample_meta_rna <- cbind(sample_id = rownames(sample_meta_rna), sample_meta_rna)

# colData often contains list-type columns (nested data from GDC API).
# Flatten them all to semicolon-separated strings so write.csv can handle them.
for (col in names(sample_meta_rna)) {
  if (is.list(sample_meta_rna[[col]])) {
    sample_meta_rna[[col]] <- sapply(sample_meta_rna[[col]], function(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      paste(x, collapse = ";")
    })
  }
}

write.csv(sample_meta_rna,
          file      = file.path(OUTPUT_DIR, "ucec_sample_metadata_rna.csv"),
          row.names = FALSE,
          quote     = FALSE)
cat("  -> Saved ucec_sample_metadata_rna.csv\n")

# ============================================================
#  2. Clinical data (survival, grade, stage, histology, age)
# ============================================================
cat("\n[2/5] Downloading clinical data...\n")

# TCGAbiolinks provides two sources — GDC and indexed clinical
# We get both and merge for maximum completeness

# (a) GDC clinical XML — richest source
clinical_gdc <- GDCquery_clinic(project = "TCGA-UCEC", type = "clinical")
write.csv(clinical_gdc,
          file      = file.path(OUTPUT_DIR, "ucec_clinical_gdc.csv"),
          row.names = FALSE,
          quote     = FALSE)
cat("  -> Saved ucec_clinical_gdc.csv (", nrow(clinical_gdc), "rows )\n")

# (b) Indexed clinical — survival endpoints
query_clin <- GDCquery(
  project       = "TCGA-UCEC",
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)
GDCdownload(query_clin, directory = GDC_CACHE_DIR)
clinical_tab <- GDCprepare(query_clin, directory = GDC_CACHE_DIR)

# GDCprepare returns a list of tables for BCR Biotab
# The main patient table is usually named "clinical_patient_ucec"
if (is.list(clinical_tab)) {
  for (tbl_name in names(clinical_tab)) {
    clean_name <- gsub("[^a-zA-Z0-9_]", "_", tbl_name)
    out_path   <- file.path(OUTPUT_DIR, paste0("ucec_clinical_", clean_name, ".csv"))
    write.csv(clinical_tab[[tbl_name]],
              file      = out_path,
              row.names = FALSE,
              quote     = FALSE)
    cat("  -> Saved ucec_clinical_", clean_name, ".csv\n")
  }
} else {
  write.csv(clinical_tab,
            file      = file.path(OUTPUT_DIR, "ucec_clinical_supplement.csv"),
            row.names = FALSE,
            quote     = FALSE)
  cat("  -> Saved ucec_clinical_supplement.csv\n")
}


# ============================================================
# 3. Molecular subtype labels
# ============================================================
cat("\n[3/5] Downloading molecular subtype labels...\n")

library(TCGAbiolinks)

subtype_df <- TCGAquery_subtype(tumor = "UCEC")

if (!is.null(subtype_df) && nrow(subtype_df) > 0) {
  
  write.csv(
    subtype_df,
    file = file.path(OUTPUT_DIR, "ucec_molecular_subtypes.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  cat(" -> Saved ucec_molecular_subtypes.csv (", nrow(subtype_df), "patients)\n")
  cat("    Columns:", paste(names(subtype_df), collapse = ", "), "\n")
  
} else {
  
  cat(" -> TCGAquery_subtype returned empty\n")
  
}


# ============================================================
#  4. Copy Number Variation (CNV) segments
#     Useful for validating CN-High subtype assignment
# ============================================================
cat("\n[4/5] Downloading CNV segment data...\n")

query_cnv <- GDCquery(
  project       = "TCGA-UCEC",
  data.category = "Copy Number Variation",
  data.type     = "Copy Number Segment",
  sample.type   = "Primary Tumor"
)

GDCdownload(query_cnv, method = "api", files.per.chunk = 20, directory = GDC_CACHE_DIR)
cnv_data <- GDCprepare(query_cnv, directory = GDC_CACHE_DIR)

write.csv(cnv_data,
          file      = file.path(OUTPUT_DIR, "ucec_cnv_segments.csv"),
          row.names = FALSE,
          quote     = FALSE)
cat("  -> Saved ucec_cnv_segments.csv (", nrow(cnv_data), "rows )\n")


# ============================================================
#  5. Somatic mutations (MAF format)
#     For TMB calculation and POLE/MSI validation
# ============================================================
cat("\n[5/5] Downloading somatic mutations (MAF)...\n")

query_maf <- GDCquery(
  project           = "TCGA-UCEC",
  data.category     = "Simple Nucleotide Variation",
  data.type         = "Masked Somatic Mutation",
  access            = "open"
)

GDCdownload(query_maf, method = "api", files.per.chunk = 10, directory = GDC_CACHE_DIR)
maf_data <- GDCprepare(query_maf, directory = GDC_CACHE_DIR)

write.csv(maf_data,
          file      = file.path(OUTPUT_DIR, "ucec_somatic_mutations.csv"),
          row.names = FALSE,
          quote     = FALSE)
cat("  -> Saved ucec_somatic_mutations.csv (", nrow(maf_data), "mutations )\n")


# ============================================================
#  6. Build master sample key
#     Maps TCGA barcodes across all data types — essential for
#     aligning RNA-seq columns with clinical rows in Python
# ============================================================
cat("\nBuilding master sample key...\n")

rna_barcodes <- data.frame(
  sample_barcode = colnames(count_mat),
  patient_barcode = substr(colnames(count_mat), 1, 12),
  stringsAsFactors = FALSE
)

write.csv(
  rna_barcodes,
  file = file.path(OUTPUT_DIR, "ucec_sample_key.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat(" -> Saved ucec_sample_key.csv\n")


# ============================================================
#  7. Session info — for reproducibility
# ============================================================
cat("\nSaving session info...\n")
sink(file.path(OUTPUT_DIR, "session_info.txt"))
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("DONE. All files saved to:", OUTPUT_DIR, "\n")
cat("============================================================\n")
cat("\nFiles created:\n")
for (f in list.files(OUTPUT_DIR, pattern = "\\.csv$")) {
  fpath <- file.path(OUTPUT_DIR, f)
  fsize <- round(file.size(fpath) / 1e6, 1)
  cat(sprintf("  %-55s %6.1f MB\n", f, fsize))
}
