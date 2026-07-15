library(data.table)

maf <- fread(
  file.path("data/raw/TCGA_UCEC_data/ucec_somatic_mutations.csv"),
  sep = ",",
  fill = TRUE,
  data.table = FALSE,
  showProgress = TRUE
)

maf <- maf[, c("Hugo_Symbol",
               "Tumor_Sample_Barcode",
               "Variant_Classification")]
maf