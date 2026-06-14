#!/usr/bin/env Rscript
# =============================================================================
# 02_tcga_analysis.R — TCGA-LUAD KEAP1 / Plunder Score analysis
#
# Outputs (results/):
#   FigA_KEAP1-MUT_correlation_heatmap.pdf/png
#   FigB_KEAP1_MUT_vs_WT_boxplot.pdf/png
#   FigC_PlunderScore_KM_OS.pdf/png
#
# Usage:
#   Rscript scripts/02_tcga_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(showtext)
  library(sysfonts)
})

# --- paths -------------------------------------------------------------------
args_ok <- TRUE
script_dir <- tryCatch({
  dirname(normalizePath(sub("--file=([^ ]+)", "\\1", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])))
}, error = function(e) NULL)

if (is.null(script_dir) || script_dir == ".") {
  script_dir <- file.path(getwd(), "scripts")
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
logs_dir    <- file.path(project_dir, "logs")
cache_dir   <- file.path(data_dir, "processed")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

status_file <- file.path(logs_dir, "02_analysis.status")
log_file    <- file.path(logs_dir, "02_analysis.log")

TARGET_GENES <- c("SLC7A11", "GGT1", "SLC1A5", "ABCC1", "ABCC2", "ABCC3", "DPEP1")

log_msg <- function(...) {
  msg <- paste(c(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ...), collapse = " ")
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

set_status <- function(step, detail = "") {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", step,
                 if (nzchar(detail)) paste0(": ", detail) else "")
  writeLines(line, status_file)
  log_msg(step, detail)
}

# --- theme -------------------------------------------------------------------
font_add("Arial", regular = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
         bold = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")
showtext_auto()

theme_pub <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    legend.background = element_rect(fill = "white", color = NA)
  )

save_dual <- function(plot_obj, prefix, w = 7, h = 6) {
  ggsave(paste0(prefix, ".pdf"), plot_obj, width = w, height = h, bg = "white")
  ggsave(paste0(prefix, ".png"), plot_obj, width = w, height = h, bg = "white", dpi = 300)
  log_msg("Saved:", prefix, "(pdf/png)")
}

patient_id <- function(x) substr(x, 1, 12)

# --- Step 1: sample mapping (GDC file_id -> TCGA barcode) --------------------
set_status("STEP 1/6", "Building sample ID mapping")
map_cache <- file.path(cache_dir, "file_sample_map.tsv")

if (!file.exists(map_cache)) {
  stop("Missing ", map_cache, " — run sample mapping first.")
} else {
  file_map <- fread(map_cache)
}

# --- Step 2: load expression (7 genes, Primary Tumor only) -------------------
set_status("STEP 2/6", "Loading STAR-Counts expression for target genes")
manifest <- fread(file.path(data_dir, "manifests", "luad_star_counts.manifest.tsv"))
star_dir <- file.path(data_dir, "star_counts")

expr_list <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  fid <- manifest$id[i]
  fname <- manifest$filename[i]
  fpath <- file.path(star_dir, fid, fname)
  if (!file.exists(fpath)) next
  dt <- fread(fpath, skip = 1, select = c("gene_name", "tpm_unstranded"))
  dt <- dt[gene_name %in% TARGET_GENES]
  if (nrow(dt) == 0) next
  sid <- file_map[file_id == fid, sample_id]
  if (length(sid) != 1) next
  vec <- setNames(dt$tpm_unstranded, dt$gene_name)
  vec <- vec[TARGET_GENES]
  names(vec) <- TARGET_GENES
  expr_list[[i]] <- c(sample_id = sid, vec)
  if (i %% 100 == 0) log_msg(sprintf("  read %d/%d count files", i, nrow(manifest)))
}

expr_mat <- as.data.table(do.call(rbind, expr_list))
expr_mat <- expr_mat[substr(sample_id, 14, 15) == "01"]  # Primary Tumor only
expr_mat <- expr_mat[!duplicated(sample_id)]
rownames_expr <- expr_mat$sample_id
expr_mat[, sample_id := NULL]
expr_mat <- as.matrix(expr_mat)
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- rownames_expr
colnames(expr_mat) <- TARGET_GENES
log_msg("Primary tumor samples:", nrow(expr_mat))

log_expr <- log2(expr_mat + 1)

# --- Step 3: KEAP1 mutation status from MAF ----------------------------------
set_status("STEP 3/6", "Parsing KEAP1 mutations from MAF files")
nonsilent <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site",
  "Nonstop_Mutation"
)
maf_files <- list.files(
  file.path(data_dir, "somatic_maf"), pattern = "\\.maf\\.gz$",
  recursive = TRUE, full.names = TRUE
)
maf_parts <- lapply(maf_files, function(f) {
  fread(f, skip = "#", select = c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"))
})
maf_dt <- rbindlist(maf_parts, fill = TRUE)
keap1_mut_samples <- unique(maf_dt[
  Hugo_Symbol == "KEAP1" & Variant_Classification %in% nonsilent,
  Tumor_Sample_Barcode
])
keap1_mut_patients <- unique(patient_id(keap1_mut_samples))
log_msg("KEAP1-MUT patients:", length(keap1_mut_patients))

sample_status <- data.table(
  sample_id = rownames_expr,
  patient_id = patient_id(rownames_expr),
  keap1_status = ifelse(patient_id(rownames_expr) %in% keap1_mut_patients,
                        "KEAP1-MUT", "KEAP1-WT")
)
fwrite(sample_status, file.path(results_dir, "sample_keap1_status.tsv"), sep = "\t")

# --- Step 4: Fig A — Spearman correlation heatmap (KEAP1-MUT) ----------------
set_status("STEP 4/6", "Figure A: correlation heatmap")
mut_samples <- sample_status[keap1_status == "KEAP1-MUT", sample_id]
mut_samples <- intersect(mut_samples, rownames(log_expr))
mut_expr <- log_expr[mut_samples, , drop = FALSE]

cor_r <- cor(mut_expr, method = "spearman")
cor_p <- matrix(NA, nrow = length(TARGET_GENES), ncol = length(TARGET_GENES),
                dimnames = list(TARGET_GENES, TARGET_GENES))
for (i in seq_along(TARGET_GENES)) {
  for (j in seq_along(TARGET_GENES)) {
    if (i == j) { cor_p[i, j] <- NA; next }
    ct <- cor.test(mut_expr[, i], mut_expr[, j], method = "spearman", exact = FALSE)
    cor_p[i, j] <- ct$p.value
  }
}

sig_labels <- matrix("", nrow = length(TARGET_GENES), ncol = length(TARGET_GENES),
                     dimnames = list(TARGET_GENES, TARGET_GENES))
for (i in seq_along(TARGET_GENES)) {
  for (j in seq_along(TARGET_GENES)) {
    if (i == j) next
    p <- cor_p[i, j]
    sig_labels[i, j] <- if (is.na(p)) "" else format.pval(p, digits = 2, eps = 0.001)
  }
}

figA_prefix <- file.path(results_dir, "FigA_KEAP1-MUT_correlation_heatmap")
pdf(paste0(figA_prefix, ".pdf"), width = 7, height = 6, bg = "white")
pheatmap(
  cor_r,
  display_numbers = sig_labels,
  number_color = "black",
  fontsize_number = 8,
  color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
  border_color = "grey90",
  main = paste0("Spearman Correlation (KEAP1-MUT, n=", nrow(mut_expr), ")"),
  fontsize = 11,
  family = "Arial"
)
dev.off()

png(paste0(figA_prefix, ".png"), width = 7, height = 6, units = "in", res = 300, bg = "white")
pheatmap(
  cor_r, display_numbers = sig_labels, number_color = "black",
  fontsize_number = 8,
  color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
  border_color = "grey90",
  main = paste0("Spearman Correlation (KEAP1-MUT, n=", nrow(mut_expr), ")"),
  fontsize = 11, family = "Arial"
)
dev.off()
log_msg("Saved Figure A")

# --- Step 5: Fig B — Boxplot MUT vs WT ---------------------------------------
set_status("STEP 5/6", "Figure B: expression boxplot")
plot_df <- as.data.table(log_expr, keep.rownames = "sample_id") |>
  merge(sample_status[, .(sample_id, keap1_status)], by = "sample_id") |>
  pivot_longer(all_of(TARGET_GENES), names_to = "gene", values_to = "log2_TPM")

plot_df$gene <- factor(plot_df$gene, levels = TARGET_GENES)
plot_df$keap1_status <- factor(plot_df$keap1_status, levels = c("KEAP1-WT", "KEAP1-MUT"))

pval_df <- plot_df |>
  group_by(gene) |>
  summarise(
    p_value = wilcox.test(log2_TPM ~ keap1_status, exact = FALSE)$p.value,
    .groups = "drop"
  ) |>
  mutate(
    p_label = format.pval(p_value, digits = 2, eps = 0.001),
    x = 1.5, y = Inf
  )

figB <- ggplot(plot_df, aes(x = keap1_status, y = log2_TPM, fill = keap1_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.35, color = "grey30") +
  geom_text(data = pval_df, aes(x = x, y = y, label = p_label),
            inherit.aes = FALSE, vjust = 1.4, size = 3.5, family = "Arial") +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_npg() +
  labs(
    title = "Predatory Metabolism Genes: KEAP1-MUT vs KEAP1-WT",
    subtitle = "Primary Tumor | log2(TPM+1) | Wilcoxon exact P",
    x = NULL, y = expression(log[2](TPM + 1)), fill = "Group"
  ) +
  theme_pub +
  theme(strip.text = element_text(face = "italic", family = "Arial"))

save_dual(figB, file.path(results_dir, "FigB_KEAP1_MUT_vs_WT_boxplot"), w = 12, h = 7)

# --- Step 6: Fig C — Plunder Score KM survival --------------------------------
set_status("STEP 6/6", "Figure C: Kaplan-Meier OS by Plunder Score")
clinical <- fread(file.path(data_dir, "clinical", "tcga_luad_clinical.tsv"))
clinical_pt <- clinical[sample_type == "Primary Tumor"]
clinical_pt[, patient_id := substr(sample_submitter_id, 1, 12)]

surv_clin <- clinical_pt[, .(
  vital_status = if (any(vital_status == "Dead")) "Dead" else "Alive",
  days_to_death = suppressWarnings(max(as.numeric(days_to_death), na.rm = TRUE)),
  days_to_follow = suppressWarnings(max(as.numeric(days_to_last_follow_up), na.rm = TRUE))
), by = patient_id]
surv_clin[, days_to_death := as.numeric(ifelse(is.infinite(days_to_death), NA, days_to_death))]
surv_clin[, days_to_follow := as.numeric(ifelse(is.infinite(days_to_follow), NA, days_to_follow))]
surv_clin[, os_time := NA_real_]
surv_clin[vital_status == "Dead" & !is.na(days_to_death), os_time := days_to_death]
surv_clin[is.na(os_time) & !is.na(days_to_follow), os_time := days_to_follow]
surv_clin[, os_event := as.integer(vital_status == "Dead")]

plunder <- data.table(
  sample_id = rownames(log_expr),
  patient_id = patient_id(rownames(log_expr)),
  plunder_score = rowMeans(log_expr)
)
plunder_pt <- plunder[, .(plunder_score = mean(plunder_score)), by = patient_id]

surv_df <- merge(plunder_pt, surv_clin, by = "patient_id", all.x = TRUE)
surv_df <- surv_df[!is.na(os_time) & os_time > 0]
med <- median(surv_df$plunder_score, na.rm = TRUE)
surv_df[, plunder_group := factor(
  ifelse(plunder_score >= med, "High", "Low"),
  levels = c("Low", "High")
)]

fit <- survfit(Surv(os_time, os_event) ~ plunder_group, data = surv_df)
logrank_p <- survdiff(Surv(os_time, os_event) ~ plunder_group, data = surv_df)
p_os <- 1 - pchisq(logrank_p$chisq, df = 1)

figC <- ggsurvplot(
  fit, data = surv_df,
  pval = FALSE, conf.int = TRUE, risk.table = TRUE,
  palette = c("#00A087", "#E64B35"),
  legend.title = "Plunder Score",
  legend.labs = c("Low", "High"),
  xlab = "Time (days)", ylab = "Overall Survival Probability",
  title = "TCGA-LUAD OS by Plunder Score (Primary Tumor)",
  font.main = c(12, "black", "bold"),
  font.x = c(11, "black", "plain"),
  font.y = c(11, "black", "plain"),
  font.tickslab = c(10, "black", "plain"),
  font.legend = c(10, "black", "plain"),
  ggtheme = theme_pub,
  tables.theme = theme_pub
)
figC$plot <- figC$plot +
  annotate("text", x = max(surv_df$os_time, na.rm = TRUE) * 0.55, y = 0.95,
           label = paste0("Log-rank P = ", format.pval(p_os, digits = 3, eps = 0.001)),
           size = 4)

figC_prefix <- file.path(results_dir, "FigC_PlunderScore_KM_OS")
pdf(paste0(figC_prefix, ".pdf"), width = 8, height = 8, bg = "white")
print(figC)
dev.off()
png(paste0(figC_prefix, ".png"), width = 8, height = 8, units = "in", res = 300, bg = "white")
print(figC)
dev.off()
log_msg("Saved Figure C")

fwrite(plunder_pt, file.path(results_dir, "plunder_scores_by_patient.tsv"), sep = "\t")
set_status("DONE", "All figures saved to results/")
