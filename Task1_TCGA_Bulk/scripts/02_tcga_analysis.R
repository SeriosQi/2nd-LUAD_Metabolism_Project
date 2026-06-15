#!/usr/bin/env Rscript
# =============================================================================
# 02_tcga_analysis.R — TCGA-LUAD KEAP1 / 9-gene predatory metabolism analysis
#
# Target genes (9):
#   SLC7A11 | ABCC1, ABCC2, ABCC3 | GGT1, GGT5, DPEP1 | SLC1A4 (ASCT1), SLC1A5 (ASCT2)
#
# Outputs (results/):
#   FigA_KEAP1-MUT_correlation_heatmap.pdf/png
#   FigB_KEAP1_MUT_vs_WT_boxplot.pdf/png
#   FigC_KM_OS_by_gene_panels.pdf/png          — 9-gene KM overview (3×3)
#   FigC_KM_OS_{GENE}.pdf/png                   — per-gene KM + risk table
#   gene_survival_logrank.tsv                     — log-rank P for each gene
#   sample_keap1_status.tsv
#   patient_gene_expression.tsv
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
  library(patchwork)
})

script_dir <- tryCatch({
  dirname(normalizePath(sub("--file=([^ ]+)", "\\1",
    commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])))
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

# Hugo symbols in GENCODE / STAR counts (ASCT1=SLC1A4, ASCT2=SLC1A5; no separate ASCT rows)
TARGET_GENES <- c(
  "SLC7A11",
  "ABCC1", "ABCC2", "ABCC3",
  "GGT1", "GGT5", "DPEP1",
  "SLC1A4", "SLC1A5"
)
GENE_LABELS <- c(
  SLC7A11 = "SLC7A11",
  ABCC1   = "ABCC1",
  ABCC2   = "ABCC2",
  ABCC3   = "ABCC3",
  GGT1    = "GGT1",
  GGT5    = "GGT5",
  DPEP1   = "DPEP1",
  SLC1A4  = "SLC1A4 (ASCT1)",
  SLC1A5  = "SLC1A5 (ASCT2)"
)

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

font_add("Arial", regular = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
         bold = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf")
showtext_auto()

theme_pub <- theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black", linewidth = 0.4),
    legend.background = element_rect(fill = "white", color = NA)
  )

save_dual <- function(plot_obj, prefix, w = 7, h = 6) {
  ggsave(paste0(prefix, ".pdf"), plot_obj, width = w, height = h, bg = "white")
  ggsave(paste0(prefix, ".png"), plot_obj, width = w, height = h, bg = "white", dpi = 300)
  log_msg("Saved:", prefix, "(pdf/png)")
}

patient_id <- function(x) substr(x, 1, 12)

build_surv_clin <- function(clinical_path) {
  clinical <- fread(clinical_path)
  clinical_pt <- clinical[sample_type == "Primary Tumor"]
  clinical_pt[, patient_id := substr(sample_submitter_id, 1, 12)]
  surv <- clinical_pt[, .(
    vital_status = if (any(vital_status == "Dead")) "Dead" else "Alive",
    days_to_death = suppressWarnings(max(as.numeric(days_to_death), na.rm = TRUE)),
    days_to_follow = suppressWarnings(max(as.numeric(days_to_last_follow_up), na.rm = TRUE))
  ), by = patient_id]
  surv[, days_to_death := as.numeric(ifelse(is.infinite(days_to_death), NA, days_to_death))]
  surv[, days_to_follow := as.numeric(ifelse(is.infinite(days_to_follow), NA, days_to_follow))]
  surv[, os_time := NA_real_]
  surv[vital_status == "Dead" & !is.na(days_to_death), os_time := days_to_death]
  surv[is.na(os_time) & !is.na(days_to_follow), os_time := days_to_follow]
  surv[, os_event := as.integer(vital_status == "Dead")]
  surv
}

# --- Step 1: sample mapping --------------------------------------------------
set_status("STEP 1/6", "Loading sample ID mapping")
map_cache <- file.path(cache_dir, "file_sample_map.tsv")
if (!file.exists(map_cache)) stop("Missing ", map_cache)
file_map <- fread(map_cache)

# --- Step 2: expression matrix (9 genes, Primary Tumor) ----------------------
set_status("STEP 2/6", "Loading STAR-Counts for 9 target genes")
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
  vec <- setNames(rep(NA_real_, length(TARGET_GENES)), TARGET_GENES)
  vec[dt$gene_name] <- dt$tpm_unstranded
  expr_list[[length(expr_list) + 1]] <- c(sample_id = sid, vec)
  if (i %% 100 == 0) log_msg(sprintf("  read %d/%d count files", i, nrow(manifest)))
}

expr_mat <- as.data.table(do.call(rbind, expr_list))
expr_mat <- expr_mat[substr(sample_id, 14, 15) == "01"]
expr_mat <- expr_mat[!duplicated(sample_id)]
rownames_expr <- expr_mat$sample_id
expr_mat[, sample_id := NULL]
expr_mat <- as.matrix(expr_mat)
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- rownames_expr
colnames(expr_mat) <- TARGET_GENES
log_msg("Primary tumor samples:", nrow(expr_mat))

log_expr <- log2(expr_mat + 1)

# patient-level mean expression (for survival)
expr_pt <- as.data.table(log_expr, keep.rownames = "sample_id")
expr_pt[, patient_id := patient_id(sample_id)]
expr_pt <- expr_pt[, lapply(.SD, mean, na.rm = TRUE), by = patient_id,
                   .SDcols = TARGET_GENES]
fwrite(expr_pt, file.path(results_dir, "patient_gene_expression.tsv"), sep = "\t")

# --- Step 3: KEAP1 status ------------------------------------------------------
set_status("STEP 3/6", "Parsing KEAP1 mutations from MAF")
nonsilent <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site",
  "Nonstop_Mutation"
)
maf_files <- list.files(file.path(data_dir, "somatic_maf"), pattern = "\\.maf\\.gz$",
                        recursive = TRUE, full.names = TRUE)
maf_parts <- lapply(maf_files, function(f) {
  fread(f, skip = "#", select = c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"))
})
maf_dt <- rbindlist(maf_parts, fill = TRUE)
keap1_mut_patients <- unique(patient_id(maf_dt[
  Hugo_Symbol == "KEAP1" & Variant_Classification %in% nonsilent,
  Tumor_Sample_Barcode
]))
log_msg("KEAP1-MUT patients:", length(keap1_mut_patients))

sample_status <- data.table(
  sample_id = rownames_expr,
  patient_id = patient_id(rownames_expr),
  keap1_status = ifelse(patient_id(rownames_expr) %in% keap1_mut_patients,
                        "KEAP1-MUT", "KEAP1-WT")
)
fwrite(sample_status, file.path(results_dir, "sample_keap1_status.tsv"), sep = "\t")

# --- Step 4: Fig A — correlation heatmap (KEAP1-MUT) ---------------------------
set_status("STEP 4/6", "Figure A: 9-gene correlation heatmap (KEAP1-MUT)")
mut_samples <- intersect(
  sample_status[keap1_status == "KEAP1-MUT", sample_id],
  rownames(log_expr)
)
mut_expr <- log_expr[mut_samples, , drop = FALSE]
cor_r <- cor(mut_expr, method = "spearman")
cor_p <- matrix(NA, 9, 9, dimnames = list(TARGET_GENES, TARGET_GENES))
sig_labels <- matrix("", 9, 9, dimnames = list(TARGET_GENES, TARGET_GENES))
for (i in seq_along(TARGET_GENES)) {
  for (j in seq_along(TARGET_GENES)) {
    if (i == j) next
    ct <- cor.test(mut_expr[, i], mut_expr[, j], method = "spearman", exact = FALSE)
    cor_p[i, j] <- ct$p.value
    sig_labels[i, j] <- format.pval(ct$p.value, digits = 2, eps = 0.001)
  }
}

label_mat <- matrix(GENE_LABELS[TARGET_GENES], nrow = 1)
dimnames(cor_r) <- list(GENE_LABELS[TARGET_GENES], GENE_LABELS[TARGET_GENES])
dimnames(sig_labels) <- dimnames(cor_r)

figA_prefix <- file.path(results_dir, "FigA_KEAP1-MUT_correlation_heatmap")
draw_heatmap <- function() {
  pheatmap(
    cor_r,
    display_numbers = sig_labels,
    number_color = "black",
    fontsize_number = 7,
    color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
    border_color = "grey90",
    main = paste0("Spearman Correlation (KEAP1-MUT, n=", nrow(mut_expr), ")"),
    fontsize = 10,
    family = "Arial",
    angle_col = 45
  )
}
pdf(paste0(figA_prefix, ".pdf"), width = 9, height = 8, bg = "white")
draw_heatmap()
dev.off()
png(paste0(figA_prefix, ".png"), width = 9, height = 8, units = "in", res = 300, bg = "white")
draw_heatmap()
dev.off()
log_msg("Saved Figure A")

# --- Step 5: Fig B — MUT vs WT boxplot -----------------------------------------
set_status("STEP 5/6", "Figure B: 9-gene boxplot (KEAP1-MUT vs WT)")
plot_df <- as.data.table(log_expr, keep.rownames = "sample_id") |>
  merge(sample_status[, .(sample_id, keap1_status)], by = "sample_id") |>
  pivot_longer(all_of(TARGET_GENES), names_to = "gene", values_to = "log2_TPM")
plot_df$gene <- factor(plot_df$gene, levels = TARGET_GENES,
                       labels = GENE_LABELS[TARGET_GENES])
plot_df$keap1_status <- factor(plot_df$keap1_status, levels = c("KEAP1-WT", "KEAP1-MUT"))

pval_df <- plot_df |>
  group_by(gene) |>
  summarise(p_value = wilcox.test(log2_TPM ~ keap1_status, exact = FALSE)$p.value,
            .groups = "drop") |>
  mutate(p_label = format.pval(p_value, digits = 2, eps = 0.001), x = 1.5, y = Inf)

figB <- ggplot(plot_df, aes(x = keap1_status, y = log2_TPM, fill = keap1_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.35, color = "grey30") +
  geom_text(data = pval_df, aes(x = x, y = y, label = p_label),
            inherit.aes = FALSE, vjust = 1.4, size = 3.2, family = "Arial") +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_npg() +
  labs(
    title = "Predatory Metabolism Genes: KEAP1-MUT vs KEAP1-WT",
    subtitle = "Primary Tumor | log2(TPM+1) | Wilcoxon P",
    x = NULL, y = expression(log[2](TPM + 1)), fill = "Group"
  ) +
  theme_pub +
  theme(strip.text = element_text(size = 9, family = "Arial"))

save_dual(figB, file.path(results_dir, "FigB_KEAP1_MUT_vs_WT_boxplot"), w = 12, h = 10)

# --- Step 6: Fig C — per-gene KM OS (patient-level) ----------------------------
set_status("STEP 6/6", "Figure C: per-gene KM overall survival")
surv_clin <- build_surv_clin(file.path(data_dir, "clinical", "tcga_luad_clinical.tsv"))

surv_results <- vector("list", length(TARGET_GENES))
km_plots <- vector("list", length(TARGET_GENES))
names(km_plots) <- TARGET_GENES

for (g in TARGET_GENES) {
  df <- merge(expr_pt[, .(patient_id, expr = get(g))], surv_clin, by = "patient_id")
  df <- df[!is.na(os_time) & os_time > 0 & !is.na(expr)]
  med <- median(df$expr, na.rm = TRUE)
  df[, expr_group := factor(ifelse(expr >= med, "High", "Low"), levels = c("Low", "High"))]

  fit <- survfit(Surv(os_time, os_event) ~ expr_group, data = df)
  lr  <- survdiff(Surv(os_time, os_event) ~ expr_group, data = df)
  p_lr <- 1 - pchisq(lr$chisq, df = 1)

  surv_results[[g]] <- data.table(
    gene = g,
    gene_label = GENE_LABELS[g],
    n_patients = nrow(df),
    n_events = sum(df$os_event),
    median_expr_cutoff = med,
    logrank_p = p_lr,
    n_low = sum(df$expr_group == "Low"),
    n_high = sum(df$expr_group == "High")
  )

  fig_g <- ggsurvplot(
    fit, data = df,
    pval = FALSE, conf.int = TRUE, risk.table = TRUE,
    palette = c("#00A087", "#E64B35"),
    legend.title = "Expression",
    legend.labs = c("Low", "High"),
    xlab = "Time (days)", ylab = "OS Probability",
    title = paste0(GENE_LABELS[g], " — TCGA-LUAD OS"),
    font.main = c(11, "black", "bold"),
    font.x = c(10, "black", "plain"),
    font.y = c(10, "black", "plain"),
    font.tickslab = c(9, "black", "plain"),
    font.legend = c(9, "black", "plain"),
    ggtheme = theme_pub,
    tables.theme = theme_pub,
    risk.table.height = 0.25
  )
  fig_g$plot <- fig_g$plot +
    annotate("text", x = max(df$os_time, na.rm = TRUE) * 0.55, y = 0.95,
             label = paste0("Log-rank P = ", format.pval(p_lr, digits = 3, eps = 0.001)),
             size = 3.5)

  prefix_g <- file.path(results_dir, paste0("FigC_KM_OS_", g))
  pdf(paste0(prefix_g, ".pdf"), width = 8, height = 8, bg = "white")
  print(fig_g)
  dev.off()
  png(paste0(prefix_g, ".png"), width = 8, height = 8, units = "in", res = 300, bg = "white")
  print(fig_g)
  dev.off()
  log_msg("Saved KM:", g)

  km_plots[[g]] <- fig_g$plot +
    labs(title = GENE_LABELS[g]) +
    theme(plot.title = element_text(size = 9, face = "bold", family = "Arial"))
}

surv_tab <- rbindlist(surv_results)
fwrite(surv_tab, file.path(results_dir, "gene_survival_logrank.tsv"), sep = "\t")

panel <- wrap_plots(km_plots, ncol = 3) +
  plot_annotation(
    title = "TCGA-LUAD OS by Target Gene Expression (patient-level, median split)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", family = "Arial"))
  )
save_dual(panel, file.path(results_dir, "FigC_KM_OS_by_gene_panels"), w = 14, h = 14)

# remove legacy plunder outputs if present
legacy <- c(
  "FigC_PlunderScore_KM_OS.pdf", "FigC_PlunderScore_KM_OS.png",
  "plunder_scores_by_patient.tsv"
)
for (f in legacy) {
  p <- file.path(results_dir, f)
  if (file.exists(p)) file.remove(p)
}

set_status("DONE", "Phase 3 complete — 9 genes, per-gene KM survival")
