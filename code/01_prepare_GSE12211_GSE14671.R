source(file.path("code", "R", "00_setup.R"))


dir.create(file.path("results", "rds"), recursive = TRUE, showWarnings = FALSE)



gse12211_rds <- file.path("results", "rds", "GSE12211_eset.rds")

if (!file.exists(gse12211_rds)) {
  message("[01] Downloading GSE12211 from GEO...")
  gse_list_12211 <- GEOquery::getGEO("GSE12211", GSEMatrix = TRUE)
  eset12211 <- if (is.list(gse_list_12211)) gse_list_12211[[1]] else gse_list_12211
  saveRDS(eset12211, gse12211_rds)
} else {
  message("[01] Loading cached GSE12211_eset.rds ...")
  eset12211 <- readRDS(gse12211_rds)
}

ph12211 <- Biobase::pData(eset12211)


title_vec <- as.character(ph12211$title)
timepoint <- ifelse(
  grepl("before",  title_vec, ignore.case = TRUE), "pre",
  ifelse(grepl("7 days", title_vec, ignore.case = TRUE), "day7", NA)
)

patient_id <- stringr::str_extract(title_vec, "Patient [0-9]+")
patient_id <- gsub("Patient ", "", patient_id)

ph12211$timepoint  <- timepoint
ph12211$patient_id <- patient_id

message("[01] GSE12211 timepoints:")
print(table(ph12211$timepoint, useNA = "ifany"))


expr12211_gene <- exprs_by_gene(eset12211)


hippo_score_12211 <- compute_score(expr12211_gene, hippo_kinase_genes)
yap_score_12211   <- compute_score(expr12211_gene, yap_tead_output_genes)


df12211_scores <- tibble::tibble(
  sample      = colnames(expr12211_gene),
  hippo_score = as.numeric(hippo_score_12211),
  yap_score   = as.numeric(yap_score_12211)
) %>%
  dplyr::left_join(
    ph12211 %>% dplyr::select(geo_accession, timepoint, patient_id),
    by = c("sample" = "geo_accession")
  )

saveRDS(df12211_scores, file.path("results", "rds", "GSE12211_scores.rds"))
message("[01] Saved results/rds/GSE12211_scores.rds")


gse14671_rds <- file.path("results", "rds", "GSE14671_eset.rds")

if (!file.exists(gse14671_rds)) {
  message("[01] Downloading GSE14671 from GEO...")
  gse_list_14671 <- GEOquery::getGEO("GSE14671", GSEMatrix = TRUE)
  eset14671 <- if (is.list(gse_list_14671)) gse_list_14671[[1]] else gse_list_14671
  saveRDS(eset14671, gse14671_rds)
} else {
  message("[01] Loading cached GSE14671_eset.rds ...")
  eset14671 <- readRDS(gse14671_rds)
}

ph14671 <- Biobase::pData(eset14671)

title_vec <- as.character(ph14671$title)
response_group <- ifelse(
  grepl("NonResponder", title_vec, ignore.case = TRUE),
  "NR",
  "R"
)

ph14671$response_group <- response_group

message("[01] GSE14671 response groups:")
print(table(ph14671$response_group, useNA = "ifany"))

expr14671_gene <- exprs_by_gene(eset14671)

hippo_score_14671 <- compute_score(expr14671_gene, hippo_kinase_genes)
yap_score_14671   <- compute_score(expr14671_gene, yap_tead_output_genes)

df14671_scores <- tibble::tibble(
  sample      = colnames(expr14671_gene),
  hippo_score = as.numeric(hippo_score_14671),
  yap_score   = as.numeric(yap_score_14671)
) %>%
  dplyr::left_join(
    ph14671 %>% dplyr::select(geo_accession, response_group),
    by = c("sample" = "geo_accession")
  )

saveRDS(df14671_scores, file.path("results", "rds", "GSE14671_scores.rds"))
message("[01] Saved results/rds/GSE14671_scores.rds")

