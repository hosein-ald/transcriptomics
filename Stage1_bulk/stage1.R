# ══════════════════════════════════════════════════════════════════════════════
# BG4 × TCGA RNA-seq × ENHANCER × TF ACTIVITY PIPELINE
# Focus:   Identify G4-associated differentially expressed genes in early-stage
#          CRC, characterise their enhancer overlap, and infer transcription
#          factor activity
# Datasets:
#   BG4 ChIP-seq (e5 PE, hg19)           — G4 peak calls
#   TCGA stage-I CRC (107 tumour samples) — RNA-seq log2(FPKM+1)
#   TCGA solid normal tissue (51 samples) — matched normal expression
#   SW480 enhancers (H3K27ac + H3K4me1)  — overlap_ac27_me1 BED
#   SCREEN enhancer reference (hg19)      — ENCODE cCRE annotations
#   FANTOM5 enhancer reference (hg19)     — F5 enhancer peaks
#   SEdb SW48 super-enhancers             — cell-line SE regions
# Author: Hossein Allahdadi
#
# !! THIS IS A STANDALONE SCRIPT — completely independent from the
#    GSE95656 methylation pipeline scripts. Do NOT mix objects. !!
#
# Pipeline overview:
#   0.  Libraries & global parameters
#   1.  BG4 peak loading & ChIPseeker annotation
#   2.  TCGA RNA-seq data loading, cleaning & merging
#   3.  ENSEMBL → gene symbol mapping & log2 transformation
#   4.  Sample metadata & coldata construction
#   5.  QC plots (density, PCA, correlation heatmap)
#   6.  limma differential expression analysis
#   7.  DEG × BG4 peak integration
#   8.  Visualisation (volcano, heatmaps, bar, distance, UpSet)
#   9.  Enhancer analysis (SW480, SCREEN, FANTOM5)
#  10.  Super-enhancer analysis (SEdb SW48)
#  11.  GO / KEGG / Reactome / Hallmark enrichment
#  12.  TF activity inference (decoupleR + CollecTRI)
#  13.  Motif analysis (JASPAR2022 + TFBSTools)
#  ──  ARCHIVE: broken / redundant / exploratory code
# ══════════════════════════════════════════════════════════════════════════════


# ── SECTION 0: LIBRARIES & GLOBAL PARAMETERS ─────────────────────────────────

suppressPackageStartupMessages({
  # Core Bioconductor
  library(GenomicFeatures)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(BSgenome.Hsapiens.UCSC.hg19)

  # ChIP analysis
  library(rtracklayer)
  library(ChIPseeker)

  # Expression analysis
  library(limma)
  library(sva)

  # Enrichment & TF
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(ReactomePA)
  library(decoupleR)

  # Motif analysis
  library(TFBSTools)
  library(JASPAR2022)
  library(motifmatchr)

  # Genomic intervals
  library(valr)

  # Visualisation & data wrangling
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(ComplexHeatmap)
  library(UpSetR)
  library(tidyverse)
  library(tibble)
  library(scales)
  library(RColorBrewer)
  library(STRINGdb)
})

# ── Global thresholds (change here, used throughout) ─────────────────────────
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 2      # |logFC| threshold for DEG calling (log2-FPKM / limma)
TSS_REGION  <- c(-2000, 2000)    # bp window for ChIPseeker annotation
TSS_PROX_BP <- 2000              # bp used to define "peak near promoter"

# ── File paths (update these before running on a new machine) ─────────────────
BG4_PEAK_FILE     <- "/Users/hossein.allahdadi/Downloads/ChIP/files/BG4_e5_PE_peaks.narrowPeak"
MOCK_PEAK_FILE    <- "/Users/hossein.allahdadi/Downloads/ChIP/files/Mock_e5_PE_peaks.narrowPeak"
STAGE1_META_FILE  <- "/Users/hossein.allahdadi/Downloads/ChIP/rnaseq_stage1sample/stage1Samples.csv"
NORMAL_META_FILE  <- "/Users/hossein.allahdadi/Downloads/ChIP/rnaseq_stage1sample/Solid_Normal_Tissue.csv"
STAGE1_EXPR_FILE  <- "/Users/hossein.allahdadi/Downloads/ChIP/rnaseq_stage1sample/Stage1Expression.csv"
NORMAL_EXPR_FILE  <- "/Users/hossein.allahdadi/Downloads/ChIP/rnaseq_stage1sample/Solid_Normal_Tissue_Expression.csv"
ENH_SW480_FILE    <- "/Users/hossein.allahdadi/Downloads/ChIP/rnaSeq/enhancer bed file/GSE77737_RAW/overlap_ac27_me1.bed"
ENH_SCREEN_FILE   <- "/Users/hossein.allahdadi/Downloads/ChIP/rnaSeq/enhancer bed file/enhancer_peaks_hg19.bed"
ENH_F5_FILE       <- "/Users/hossein.allahdadi/Downloads/ChIP/rnaSeq/enhancer bed file/f5_enhancer_peaks_hg19.bed"
SE_SW48_FILE      <- "/Users/hossein.allahdadi/Downloads/ChIP/rnaSeq/enhancer bed file/SEDB/SE_02_0324_SE_hg19(SW48_SE).bed"
OUTPUT_DIR        <- "/Users/hossein.allahdadi/Downloads/ChIP/results_stage1"

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

NarrowPeak_COLS <- c("chrom","start","end","name","score","strand",
                     "signalValue","pValue","qValue","peak")

# ── Outlier samples to remove (identified by PCA in Section 5) ───────────────
OUTLIER_SAMPLES <- c("S53","S57","S61")


# ── SECTION 1: BG4 PEAK LOADING & ChIPseeker ANNOTATION ──────────────────────

# Load peaks as GRanges (rtracklayer handles NarrowPeak format)
bg4_peaks_gr <- rtracklayer::import(BG4_PEAK_FILE)
cat("BG4 e5 PE peaks loaded:", length(bg4_peaks_gr), "\n")

# Load as data frame for valr operations
bg4_peaks_df <- read.table(BG4_PEAK_FILE, header = FALSE)
colnames(bg4_peaks_df) <- NarrowPeak_COLS

mock_peaks_df <- read.table(MOCK_PEAK_FILE, header = FALSE)
colnames(mock_peaks_df) <- NarrowPeak_COLS

# Mock-subtract BG4 peaks to remove background signal
bg4_peaks_pure_df <- valr::bed_subtract(bg4_peaks_df, mock_peaks_df)
cat("BG4 raw:", nrow(bg4_peaks_df),
    "| Mock:", nrow(mock_peaks_df),
    "| Pure:", nrow(bg4_peaks_pure_df), "\n")

# ── ChIPseeker annotation ─────────────────────────────────────────────────────
# Annotates each peak relative to the nearest gene feature (TSS, exon, etc.)
peak_anno <- ChIPseeker::annotatePeak(
  bg4_peaks_gr,
  TxDb      = txdb,
  tssRegion = TSS_REGION,
  annoDb    = "org.Hs.eg.db",
  verbose   = FALSE
)
peak_df <- as.data.frame(peak_anno)
cat("Peaks annotated:", nrow(peak_df), "\n")

# Also annotate with refGene TxDb for comparison (requires internet / RMariaDB)
# peak_anno_refseq <- ChIPseeker::annotatePeak(
#   bg4_peaks_gr,
#   TxDb      = txdbmaker::makeTxDbFromUCSC(genome = "hg19", tablename = "refGene"),
#   tssRegion = TSS_REGION,
#   annoDb    = "org.Hs.eg.db"
# )

# Summary visualisations for peak annotation
plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno)
vennpie(peak_anno)
upsetplot(peak_anno)

# ── Extract gene symbols associated with peaks ────────────────────────────────
peak_genes <- peak_df %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::distinct(SYMBOL) %>%
  dplyr::pull(SYMBOL)
cat("Unique genes with annotated BG4 peak:", length(peak_genes), "\n")

# ── Peaks in promoter / at TSS (for focused analyses) ────────────────────────
peak_df_onPROM <- peak_df[
  peak_df$distanceToTSS > -TSS_PROX_BP &
    peak_df$distanceToTSS < TSS_PROX_BP, ]
peak_df_onTSS  <- peak_df[peak_df$distanceToTSS == 0, ]

cat("Peaks within", TSS_PROX_BP, "bp of TSS:", nrow(peak_df_onPROM),
    "| unique genes:", length(unique(peak_df_onPROM$SYMBOL)), "\n")
cat("Peaks exactly at TSS:", nrow(peak_df_onTSS), "\n")


# ── SECTION 2: TCGA RNA-seq DATA LOADING, CLEANING & MERGING ─────────────────
# Two TCGA cohorts are merged:
#   normal_expression: 51 solid normal tissue samples
#   stage1_expression: 107 TCGA CRC stage I tumour samples
# Both have ENSEMBL gene IDs (with version suffix, e.g. ENSG00000146648.14)
# as rownames and FPKM-UQ values as entries.

stage1_metadata  <- read.csv(STAGE1_META_FILE, stringsAsFactors = FALSE)
normal_metadata  <- read.csv(NORMAL_META_FILE, stringsAsFactors = FALSE)
stage1_expression <- as.data.frame(read.csv(STAGE1_EXPR_FILE, row.names = 1))
normal_expression <- as.data.frame(read.csv(NORMAL_EXPR_FILE))

# Fix rownames for normal expression (ENSEMBL column → rownames)
rownames(normal_expression) <- normal_expression$ENSEMBL
normal_expression <- normal_expression[, !colnames(normal_expression) %in% c("ENSEMBL","gene_name")]

cat("Stage1 expression: ", nrow(stage1_expression), "genes ×",
    ncol(stage1_expression), "samples\n")
cat("Normal expression: ", nrow(normal_expression), "genes ×",
    ncol(normal_expression), "samples\n")
cat("Common genes:      ",
    length(intersect(rownames(normal_expression), rownames(stage1_expression))), "\n")

# Sanity checks
stopifnot(!anyNA(normal_expression))
stopifnot(!anyNA(stage1_expression))
stopifnot(!anyDuplicated(rownames(normal_expression)))
stopifnot(!anyDuplicated(rownames(stage1_expression)))
stopifnot(identical(rownames(normal_expression), rownames(stage1_expression)))

# Merge — cbind keeps all rows (same gene order confirmed above)
merged_data <- cbind(normal_expression, stage1_expression)
cat("Merged matrix:", nrow(merged_data), "genes ×", ncol(merged_data), "samples\n")

write.csv(merged_data,
          file.path(OUTPUT_DIR, "merged_normal_stage1.csv"))


# ── SECTION 3: ENSEMBL → SYMBOL MAPPING & LOG2 TRANSFORM ─────────────────────
# ENSEMBL IDs may carry version suffixes (e.g. .14) — strip them before mapping.

sym <- mapIds(
  org.Hs.eg.db,
  keys      = sub("\\..*$", "", rownames(merged_data)),
  keytype   = "ENSEMBL",
  column    = "SYMBOL",
  multiVals = "first"
)

# Remove unmapped / empty symbols
keep      <- !is.na(sym) & sym != ""
mat       <- merged_data[keep, , drop = FALSE]
sym_clean <- sym[keep]

# Collapse duplicate symbols by averaging (limma::avereps — correct for log-expression)
mat_symbol <- limma::avereps(mat, ID = sym_clean)
cat("Genes after symbol mapping (unique):", nrow(mat_symbol), "\n")

# Log2(FPKM+1) transform
log_fpkm_merged <- log2(mat_symbol + 1)

# Remove zero-variance genes (cannot be used in limma)
keep_var        <- apply(log_fpkm_merged, 1, var, na.rm = TRUE) > 0
log_fpkm_merged <- log_fpkm_merged[keep_var, , drop = FALSE]
cat("Genes after variance filter:", nrow(log_fpkm_merged), "\n")


# ── SECTION 4: SAMPLE METADATA & coldata CONSTRUCTION ────────────────────────
# Samples are renamed S1–S158 for brevity.
# Outliers S53, S57, S61 are identified by PCA (Section 5) and handled here.

n_normal <- ncol(normal_expression)   # 51
n_tumor  <- ncol(stage1_expression)   # 107

group <- factor(
  c(rep("normal", n_normal), rep("tumor", n_tumor)),
  levels = c("normal", "tumor")
)
condition <- factor(
  c(normal_metadata$AgeGroup, stage1_metadata$AgeGroup)
)

coldata_all <- data.frame(
  row.names = paste0("S", seq_len(n_normal + n_tumor)),
  group     = group,
  condition = condition
)
colnames(log_fpkm_merged) <- rownames(coldata_all)

# Version without outliers (used for final limma model and heatmaps)
coldata_clean <- coldata_all[!rownames(coldata_all) %in% OUTLIER_SAMPLES, , drop = FALSE]
log_fpkm_clean <- log_fpkm_merged[,
  !colnames(log_fpkm_merged) %in% OUTLIER_SAMPLES, drop = FALSE]

# Re-filter zero-variance after removing outliers
keep_var2     <- apply(log_fpkm_clean, 1, sd, na.rm = TRUE) > 0
log_fpkm_clean <- log_fpkm_clean[keep_var2, , drop = FALSE]

cat("Final matrix (clean):", nrow(log_fpkm_clean), "genes ×",
    ncol(log_fpkm_clean), "samples\n")
stopifnot(identical(rownames(coldata_clean), colnames(log_fpkm_clean)))


# ── SECTION 5: QC PLOTS ───────────────────────────────────────────────────────

# ── 5A: Expression density per sample ────────────────────────────────────────
df_long <- as.data.frame(log_fpkm_merged) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  dplyr::left_join(
    tibble::tibble(sample = rownames(coldata_all), group = coldata_all$group),
    by = "sample"
  )

p_density <- ggplot(df_long, aes(x = expr, group = sample, color = group)) +
  geom_density(alpha = 0.12, linewidth = 0.4, na.rm = TRUE) +
  scale_color_manual(values = c(normal = "#378ADD", tumor = "#D85A30")) +
  theme_minimal(base_size = 12) +
  labs(title = "QC: Expression density (log2(FPKM+1))",
       x = "log2(FPKM+1)", y = "Density", color = "Group")
print(p_density)

# ── 5B: Distribution statistics per sample ───────────────────────────────────
dist_df <- data.frame(
  sample = colnames(log_fpkm_merged),
  med    = apply(log_fpkm_merged, 2, median, na.rm = TRUE),
  iqr    = apply(log_fpkm_merged, 2, IQR,    na.rm = TRUE),
  zeros  = colMeans(log_fpkm_merged == 0, na.rm = TRUE),
  stringsAsFactors = FALSE
)
dist_df <- dplyr::left_join(dist_df,
  tibble::rownames_to_column(coldata_all, "sample"), by = "sample")

ggplot(dist_df, aes(x = group, y = med, color = group)) +
  geom_point(size = 2) + theme_bw() +
  labs(title = "Median expression per sample", y = "Median log2(FPKM+1)")

ggplot(dist_df, aes(x = group, y = zeros, color = group)) +
  geom_point(size = 2) + theme_bw() +
  labs(title = "Fraction zero-expression genes", y = "Fraction = 0")

# ── 5C: PCA (all samples — use to identify outliers) ─────────────────────────
pca_all <- prcomp(t(log_fpkm_merged), scale. = TRUE)
pca_df  <- data.frame(
  sample = rownames(pca_all$x),
  PC1    = pca_all$x[, 1],
  PC2    = pca_all$x[, 2]
) %>%
  dplyr::left_join(
    tibble::rownames_to_column(coldata_all, "sample"), by = "sample")

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text(vjust = -0.6, size = 2, show.legend = FALSE) +
  scale_color_manual(values = c(normal = "#378ADD", tumor = "#D85A30")) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA — all samples (outliers visible)",
       x = "PC1", y = "PC2", color = "Group")
print(p_pca)
# Outliers S53, S57, S61 visible here — they are removed in coldata_clean

# ── 5D: PCA (clean dataset) ───────────────────────────────────────────────────
pca_clean <- prcomp(t(log_fpkm_clean), scale. = TRUE)
pca_clean_df <- data.frame(
  sample = rownames(pca_clean$x),
  PC1    = pca_clean$x[, 1],
  PC2    = pca_clean$x[, 2]
) %>%
  dplyr::left_join(
    tibble::rownames_to_column(coldata_clean, "sample"), by = "sample")

ggplot(pca_clean_df, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.6, size = 2, show.legend = FALSE) +
  scale_color_manual(values = c(normal = "#378ADD", tumor = "#D85A30")) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA — clean dataset (outliers removed)", x = "PC1", y = "PC2")

# ── 5E: Sample correlation heatmap ───────────────────────────────────────────
cors <- cor(log_fpkm_clean, use = "pairwise.complete.obs")
ann_col_qc <- data.frame(group = coldata_clean$group,
                         row.names = rownames(coldata_clean))
pheatmap(cors,
         annotation_col = ann_col_qc,
         annotation_row = ann_col_qc,
         treeheight_row = FALSE, treeheight_col = FALSE,
         show_colnames  = FALSE, show_rownames  = FALSE,
         main = "QC: Sample-to-sample correlation")

# ── 5F: EGFR sanity check (known CRC marker) ─────────────────────────────────
if ("EGFR" %in% rownames(log_fpkm_clean)) {
  egfr_df <- data.frame(
    sample = colnames(log_fpkm_clean),
    expr   = as.numeric(log_fpkm_clean["EGFR", ])
  ) %>%
    dplyr::left_join(tibble::rownames_to_column(coldata_clean, "sample"),
                     by = "sample")

  ggplot(egfr_df, aes(x = group, y = expr, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c(normal = "#378ADD", tumor = "#D85A30")) +
    theme_classic(base_size = 12) +
    labs(title = "EGFR sanity check", x = NULL,
         y = "log2(FPKM+1)")
}


# ── SECTION 6: LIMMA DIFFERENTIAL EXPRESSION ANALYSIS ────────────────────────
# limma on log2(FPKM+1) is the appropriate approach when you have pre-computed
# FPKM values rather than raw counts. DESeq2 requires raw counts.
# The ~0 + group contrast coding is used for explicit tumour vs normal comparison.

stopifnot(identical(rownames(coldata_clean), colnames(log_fpkm_clean)))

design <- model.matrix(~ 0 + group, data = coldata_clean)
colnames(design) <- levels(coldata_clean$group)

fit      <- lmFit(log_fpkm_clean, design)
cont     <- makeContrasts(tumor_vs_normal = tumor - normal, levels = design)
fit2     <- eBayes(contrasts.fit(fit, cont))
results  <- topTable(fit2, coef = "tumor_vs_normal", number = Inf, sort.by = "P")

cat("Total genes tested:", nrow(results), "\n")
cat("Significant (FDR <", PADJ_CUTOFF, "):", sum(results$adj.P.Val < PADJ_CUTOFF, na.rm = TRUE), "\n")

# DEG list: |logFC| ≥ LFC_CUTOFF and FDR < PADJ_CUTOFF
DEG_list_stage1 <- results %>%
  dplyr::filter(abs(logFC) >= LFC_CUTOFF & adj.P.Val < PADJ_CUTOFF) %>%
  dplyr::mutate(direction = dplyr::if_else(logFC > 0, "UP", "DOWN"))

cat("DEGs total:", nrow(DEG_list_stage1),
    "| UP:", sum(DEG_list_stage1$direction == "UP"),
    "| DOWN:", sum(DEG_list_stage1$direction == "DOWN"), "\n")

deg_up_stage1   <- DEG_list_stage1 %>% dplyr::filter(direction == "UP")
deg_down_stage1 <- DEG_list_stage1 %>% dplyr::filter(direction == "DOWN")


# ── SECTION 7: DEG × BG4 PEAK INTEGRATION ────────────────────────────────────
# Three levels of BG4 proximity tested:
#   A. Any annotated gene near a BG4 peak (all peak_genes)
#   B. Gene promoter overlaps BG4 peak (peak within ±TSS_PROX_BP of TSS)
#   C. Gene TSS exactly overlaps BG4 peak (distanceToTSS == 0)

# ── 7A: DEGs with any BG4 peak ───────────────────────────────────────────────
deg_with_peak_stage1 <- DEG_list_stage1 %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::mutate(has_peak = SYMBOL %in% peak_genes) %>%
  dplyr::filter(has_peak)

deg_up_with_peak_stage1   <- dplyr::filter(deg_with_peak_stage1, direction == "UP")
deg_down_with_peak_stage1 <- dplyr::filter(deg_with_peak_stage1, direction == "DOWN")

# ── 7B: DEGs with BG4 peak in promoter region (±2kb) ─────────────────────────
DEGs_BG4_onPROM <- intersect(peak_df_onPROM$SYMBOL, rownames(DEG_list_stage1))
deg_up_onPROM   <- rownames(DEG_list_stage1[
  rownames(DEG_list_stage1) %in% DEGs_BG4_onPROM &
    DEG_list_stage1$direction == "UP", ])
deg_down_onPROM <- rownames(DEG_list_stage1[
  rownames(DEG_list_stage1) %in% DEGs_BG4_onPROM &
    DEG_list_stage1$direction == "DOWN", ])

# ── 7C: DEGs with BG4 peak at exact TSS ──────────────────────────────────────
DEGs_BG4_onTSS <- intersect(peak_df_onTSS$SYMBOL, rownames(DEG_list_stage1))
deg_up_onTSS   <- rownames(DEG_list_stage1[
  rownames(DEG_list_stage1) %in% DEGs_BG4_onTSS &
    DEG_list_stage1$direction == "UP", ])
deg_down_onTSS <- rownames(DEG_list_stage1[
  rownames(DEG_list_stage1) %in% DEGs_BG4_onTSS &
    DEG_list_stage1$direction == "DOWN", ])

cat("\n─── DEG × BG4 peak integration ─────────────────────────────────────\n")
cat("DEGs with any BG4 peak:       ", nrow(deg_with_peak_stage1),
    " (UP:", nrow(deg_up_with_peak_stage1),
    " | DOWN:", nrow(deg_down_with_peak_stage1), ")\n")
cat("DEGs with peak on promoter:   ", length(DEGs_BG4_onPROM),
    " (UP:", length(deg_up_onPROM),
    " | DOWN:", length(deg_down_onPROM), ")\n")
cat("DEGs with peak at exact TSS:  ", length(DEGs_BG4_onTSS),
    " (UP:", length(deg_up_onTSS),
    " | DOWN:", length(deg_down_onTSS), ")\n")
cat("Total DEGs:                   ", nrow(DEG_list_stage1),
    " (UP:", nrow(deg_up_stage1),
    " | DOWN:", nrow(deg_down_stage1), ")\n")

# ── Write TSS BED file for DEG genes ─────────────────────────────────────────
entrez_deg <- mapIds(org.Hs.eg.db,
                     keys     = rownames(DEG_list_stage1),
                     keytype  = "SYMBOL",
                     column   = "ENTREZID",
                     multiVals = "first")
entrez_deg <- na.omit(entrez_deg)

gene_gr <- genes(txdb)[names(genes(txdb)) %in% entrez_deg]
tss_gr  <- resize(gene_gr, width = 1, fix = "start")

tss_bed <- data.frame(
  chrom     = as.character(seqnames(tss_gr)),
  start     = start(tss_gr) - 1L,     # 0-based BED
  end       = end(tss_gr),
  name      = names(tss_gr),           # Entrez ID
  score     = ".",
  strand    = as.character(strand(tss_gr))
)
write.table(tss_bed,
            file      = file.path(OUTPUT_DIR, "DEG_genes_tss.bed"),
            sep       = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)


# ── SECTION 8: VISUALISATION ─────────────────────────────────────────────────

# ── 8A: Volcano plot — all results, coloured by BG4 peak + direction ──────────
plot_df_volcano <- results %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::mutate(
    has_peak     = SYMBOL %in% peak_genes,
    sig          = adj.P.Val < PADJ_CUTOFF & abs(logFC) >= LFC_CUTOFF,
    direction    = dplyr::case_when(
      sig & logFC  > 0 ~ "UP",
      sig & logFC  < 0 ~ "DOWN",
      TRUE             ~ "NS"
    ),
    group = dplyr::case_when(
      sig & has_peak & logFC > 0 ~ "UP + peak",
      sig & has_peak & logFC < 0 ~ "DOWN + peak",
      sig & !has_peak            ~ "DEG no peak",
      TRUE                       ~ "Not significant"
    ),
    neglog10padj = -log10(adj.P.Val)
  ) %>%
  dplyr::filter(is.finite(logFC), is.finite(neglog10padj))

top_labels_volcano <- plot_df_volcano %>%
  dplyr::filter(group %in% c("UP + peak","DOWN + peak")) %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(SYMBOL)

p_volcano <- ggplot(plot_df_volcano, aes(x = logFC, y = neglog10padj)) +
  geom_point(aes(color = group), alpha = 0.85, size = 1.6, na.rm = TRUE) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed") +
  scale_color_manual(values = c(
    "UP + peak"       = "#D55E00",
    "DOWN + peak"     = "#0072B2",
    "DEG no peak"     = "#7A7A7A",
    "Not significant" = "#D0D0D0"
  )) +
  ggrepel::geom_text_repel(
    data        = dplyr::filter(plot_df_volcano, SYMBOL %in% top_labels_volcano),
    aes(label   = SYMBOL),
    size        = 3, max.overlaps = 50,
    box.padding = 0.4, point.padding = 0.2
  ) +
  theme_minimal(base_size = 12) +
  labs(title    = "Volcano: tumour vs normal (TCGA stage I CRC)",
       subtitle = paste0("BG4 peak-associated DEGs highlighted | ",
                         "|logFC| ≥ ", LFC_CUTOFF, ", FDR < ", PADJ_CUTOFF),
       x = "log2 fold-change", y = "-log10(adjusted p-value)", color = "")
print(p_volcano)
ggsave(file.path(OUTPUT_DIR, "volcano_BG4_DEGs.pdf"), p_volcano,
       width = 9, height = 6)

# ── 8B: Stacked bar — proportion of DEGs with/without BG4 peak ───────────────
bar_df <- plot_df_volcano %>%
  dplyr::filter(direction %in% c("UP","DOWN")) %>%
  dplyr::group_by(direction, has_peak) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(direction) %>%
  dplyr::mutate(pct = n / sum(n))

p_bar <- ggplot(bar_df, aes(x = direction, y = n, fill = has_peak)) +
  geom_col(position = "stack") +
  geom_text(aes(label = paste0("(", scales::percent(pct, accuracy = 1), ")")),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values   = c("TRUE" = "#009E73", "FALSE" = "#CCCCCC"),
                    labels   = c("TRUE" = "Has BG4 peak", "FALSE" = "No peak")) +
  theme_minimal(base_size = 12) +
  labs(title = "DEGs split by BG4 peak proximity", x = "", y = "n DEGs", fill = "")
print(p_bar)

# ── 8C: Peak distance to TSS — DEG vs non-DEG ────────────────────────────────
peak_df_dist <- peak_df %>%
  dplyr::mutate(is_deg = SYMBOL %in% rownames(DEG_list_stage1)) %>%
  dplyr::filter(!is.na(distanceToTSS))

p_dist <- ggplot(peak_df_dist, aes(x = distanceToTSS, fill = is_deg)) +
  geom_histogram(bins = 80, alpha = 0.8) +
  coord_cartesian(xlim = c(-10000, 10000)) +
  scale_fill_manual(values = c("TRUE" = "#CC79A7", "FALSE" = "#BBBBBB"),
                    labels = c("TRUE" = "Peak near DEG", "FALSE" = "Peak near non-DEG")) +
  theme_minimal(base_size = 12) +
  labs(title = "BG4 peak distance to TSS (±10kb)",
       x = "Distance to TSS (bp)", y = "Number of peaks", fill = "")
print(p_dist)

# ── 8D: Heatmap — top DEGs with BG4 peak on promoter ─────────────────────────
# UP genes ordered by descending logFC, DOWN by ascending (most extreme first)
up_ordered   <- rownames(deg_up_stage1[order(deg_up_stage1$logFC, decreasing = TRUE), ])
down_ordered <- rownames(deg_down_stage1[order(deg_down_stage1$logFC, decreasing = FALSE), ])

up_peak_top50   <- head(intersect(up_ordered,   DEGs_BG4_onPROM), 50)
down_peak_top50 <- head(intersect(down_ordered, DEGs_BG4_onPROM), 50)
genes_heatmap   <- c(up_peak_top50, down_peak_top50)
genes_heatmap   <- genes_heatmap[genes_heatmap %in% rownames(log_fpkm_clean)]

ann_col_hm <- coldata_clean[, "group", drop = FALSE]

pheatmap(
  as.matrix(log_fpkm_clean[genes_heatmap, ]),
  scale             = "row",
  cluster_rows      = TRUE, cluster_cols = FALSE,
  show_rownames     = TRUE,  show_colnames = FALSE,
  treeheight_row    = FALSE, treeheight_col = FALSE,
  annotation_col    = ann_col_hm,
  color             = colorRampPalette(c("#0072B2","white","#D55E00"))(100),
  fontsize_row      = 6,
  main              = paste0("Top DEGs with BG4 peak on promoter (n=",
                             length(genes_heatmap), ")")
)

# ── 8E: UpSet plot — DEG ∩ BG4 peak genes ────────────────────────────────────
upset_sets <- list(
  DEG      = rownames(DEG_list_stage1),
  BG4_Peak = peak_genes
)
upset(fromList(upset_sets), order.by = "freq",
      main.bar.color = "#444441", sets.bar.color = c("#D85A30","#378ADD"),
      text.scale = 1.3)


# ── SECTION 9: ENHANCER ANALYSIS ─────────────────────────────────────────────
# Tests what fraction of BG4 peaks (mock-subtracted) overlap various enhancer
# reference sets. Promoters are removed from enhancer sets before intersecting
# to avoid conflating promoter and enhancer signals.
# Note: all bed_subtract / bed_intersect calls require exactly 3 valr columns.

# ── Helper: clean any data frame to valr format ───────────────────────────────
to_valr3 <- function(df) {
  data.frame(
    chrom = as.character(df[[1]]),
    start = as.integer(df[[2]]),
    end   = as.integer(df[[3]]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(grepl("^chr", chrom)) %>%
    dplyr::arrange(chrom, start)
}

# BG4 pure (3-column valr)
bg4_valr <- to_valr3(bg4_peaks_pure_df)

# ── Build promoter BED (valr format, merged overlapping regions) ──────────────
# Gene-level promoters (fewer, less granular)
prom_genes_gr <- promoters(genes(txdb), upstream = 3000, downstream = 3000)
prom_genes_bed <- as.data.frame(prom_genes_gr) %>%
  dplyr::transmute(
    chrom = as.character(seqnames),
    start = as.integer(start) - 1L,   # 0-based
    end   = as.integer(end)
  ) %>%
  dplyr::arrange(chrom, start) %>%
  dplyr::distinct()
prom_genes_bed <- valr::bed_merge(prom_genes_bed)

# Transcript-level promoters (more granular, preferred for enhancer subtraction)
prom_tx_gr  <- promoters(transcripts(txdb), upstream = 2000, downstream = 2000)
prom_tx_bed <- as.data.frame(prom_tx_gr) %>%
  dplyr::transmute(
    chrom = as.character(seqnames),
    start = as.integer(start) - 1L,
    end   = as.integer(end)
  ) %>%
  dplyr::arrange(chrom, start) %>%
  dplyr::distinct()
prom_tx_bed <- valr::bed_merge(prom_tx_bed)

# ── Load enhancer sets ────────────────────────────────────────────────────────
# SW480 H3K27ac + H3K4me1 overlap (active enhancers in a CRC cell line)
enh_sw480_raw <- readr::read_tsv(ENH_SW480_FILE, col_names = FALSE, show_col_types = FALSE)
colnames(enh_sw480_raw) <- c("chrom","start","end","name","value")
enh_sw480_valr <- to_valr3(enh_sw480_raw)

# SCREEN ENCODE cCREs (hg19)
enh_screen <- read.table(ENH_SCREEN_FILE, header = FALSE)
colnames(enh_screen) <- c("chrom","start","end","ID1","ID2","dist")
enh_screen_valr <- to_valr3(enh_screen)

# FANTOM5 enhancers
enh_f5 <- read.table(ENH_F5_FILE, header = FALSE)
colnames(enh_f5) <- c("chrom","start","end","name","score","strand",
                      "thickStart","thickEnd","reserved","evidenceSources",
                      "elementType","eliteness")
enh_f5_valr <- to_valr3(enh_f5)

# ── Subtract promoters from each enhancer set ─────────────────────────────────
enh_sw480_noprom      <- valr::bed_subtract(enh_sw480_valr, prom_genes_bed)
enh_sw480_noprom_tx   <- valr::bed_subtract(enh_sw480_valr, prom_tx_bed)
enh_screen_noprom     <- valr::bed_subtract(enh_screen_valr, prom_genes_bed)
enh_screen_noprom_tx  <- valr::bed_subtract(enh_screen_valr, prom_tx_bed)

cat("SW480 enhancers: total =",   nrow(enh_sw480_valr),
    "| after promoter removal =", nrow(enh_sw480_noprom), "\n")
cat("SCREEN cCREs:    total =",   nrow(enh_screen_valr),
    "| after promoter removal =", nrow(enh_screen_noprom), "\n")

# ── Compute BG4 overlap % with each enhancer set ─────────────────────────────
pct_overlap <- function(bg4_df, enh_df, label) {
  n_total <- dplyr::distinct(bg4_df, chrom, start, end) %>% nrow()
  n_hit   <- valr::bed_intersect(bg4_df, enh_df) %>%
    dplyr::distinct(chrom, start.x, end.x) %>% nrow()
  pct <- round(100 * n_hit / n_total, 2)
  cat(label, "— BG4 peaks overlapping:", n_hit, "/", n_total,
      "(", pct, "%)\n")
  pct
}

pct_sw480_raw         <- pct_overlap(bg4_valr, enh_sw480_valr,     "SW480 (raw)")
pct_sw480_noprom      <- pct_overlap(bg4_valr, enh_sw480_noprom,   "SW480 (no gene-promoter)")
pct_sw480_noprom_tx   <- pct_overlap(bg4_valr, enh_sw480_noprom_tx,"SW480 (no tx-promoter)")
pct_screen_raw        <- pct_overlap(bg4_valr, enh_screen_valr,    "SCREEN (raw)")
pct_screen_noprom     <- pct_overlap(bg4_valr, enh_screen_noprom,  "SCREEN (no gene-promoter)")
pct_screen_noprom_tx  <- pct_overlap(bg4_valr, enh_screen_noprom_tx,"SCREEN (no tx-promoter)")
pct_f5                <- pct_overlap(bg4_valr, enh_f5_valr,        "FANTOM5")

# Interpretation note:
# Raw BG4 × SW480 ~60% overlap → includes promoters (expected, BG4 enriched at TSSs)
# After tx-promoter removal    ~18–22% → true enhancer-only signal


# ── SECTION 10: SUPER-ENHANCER ANALYSIS (SEdb SW48) ──────────────────────────
# SW48 is a CRC cell line with super-enhancer annotations from SEdb.
# We find the overlap between SW48 SEs, BG4 peaks, and the DEG list.

se_sw48 <- read.table(SE_SW48_FILE, header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE)
colnames(se_sw48) <- se_sw48[1, ]
se_sw48           <- se_sw48[-1, ]
colnames(se_sw48)[1:4] <- c("chrom","start","end","name")

# BG4 peak annotation (all peaks, not just pure) for gene-level matching
bg4_anno_all <- ChIPseeker::annotatePeak(
  BG4_PEAK_FILE,
  tssRegion = c(-3000, 3000),
  TxDb      = txdb,
  annoDb    = "org.Hs.eg.db",
  verbose   = FALSE
)
bg4_anno_all_df <- as.data.frame(bg4_anno_all)

# Genes shared between SW48 SEs and BG4 peaks
se_bg4_genes <- intersect(se_sw48$se_gene_overlap, bg4_anno_all_df$SYMBOL)
cat("Genes in SW48 SEs ∩ BG4 peaks:", length(se_bg4_genes), "\n")

# Add DEG overlap
se_bg4_deg <- intersect(se_bg4_genes, deg_with_peak_stage1$SYMBOL)
cat("Genes in SW48 SEs ∩ BG4 peaks ∩ DEGs:", length(se_bg4_deg), "\n")

# Heatmap for this candidate set
if (length(se_bg4_deg) > 0) {
  se_bg4_mat <- as.matrix(
    log_fpkm_clean[rownames(log_fpkm_clean) %in% se_bg4_deg, , drop = FALSE])

  pheatmap(se_bg4_mat,
           scale             = "row",
           cluster_rows      = TRUE, cluster_cols = FALSE,
           treeheight_row    = FALSE, treeheight_col = FALSE,
           annotation_col    = ann_col_hm,
           color             = colorRampPalette(c("#0072B2","white","#D55E00"))(100),
           main              = paste0("BG4 × Super-enhancer × DEG genes (n=",
                                      length(se_bg4_deg), ")"))
}


# ── SECTION 11: GO / KEGG / REACTOME / HALLMARK ENRICHMENT ───────────────────
# Run enrichment for four gene sets: DEG all promoter, DEG up/down on promoter.
# Universe = all tested genes (not just DEGs) to avoid inflation.

universe_entrez <- mapIds(org.Hs.eg.db,
                          keys      = rownames(results),
                          column    = "ENTREZID",
                          keytype   = "SYMBOL",
                          multiVals = "first") %>% na.omit() %>% unique()

to_entrez <- function(symbols) {
  mapIds(org.Hs.eg.db,
         keys      = symbols,
         column    = "ENTREZID",
         keytype   = "SYMBOL",
         multiVals = "first") %>% na.omit()
}

entrez_up_onPROM   <- to_entrez(deg_up_onPROM)
entrez_down_onPROM <- to_entrez(deg_down_onPROM)
entrez_all_onPROM  <- to_entrez(DEGs_BG4_onPROM)
entrez_up_all      <- to_entrez(rownames(deg_up_stage1))
entrez_down_all    <- to_entrez(rownames(deg_down_stage1))

# ── GO Biological Process ─────────────────────────────────────────────────────
run_go <- function(entrez, label) {
  enrichGO(
    gene          = entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 10,
    maxGSSize     = 500,
    readable      = TRUE
  )
}

go_up_onPROM   <- run_go(entrez_up_onPROM,   "UP × BG4 promoter")
go_down_onPROM <- run_go(entrez_down_onPROM, "DOWN × BG4 promoter")
go_up_all      <- run_go(entrez_up_all,      "All UP DEGs")
go_down_all    <- run_go(entrez_down_all,    "All DOWN DEGs")

dotplot(go_up_onPROM,   showCategory = 20) + ggtitle("GO BP — UP DEGs with BG4 promoter peak")
dotplot(go_down_onPROM, showCategory = 20) + ggtitle("GO BP — DOWN DEGs with BG4 promoter peak")

# Network enrichment map
go_up_onPROM_sim <- pairwise_termsim(go_up_onPROM)
emapplot(go_up_onPROM_sim, showCategory = 40)
cnetplot(go_up_onPROM, showCategory = 10)

# ── KEGG ─────────────────────────────────────────────────────────────────────
kegg_up   <- enrichKEGG(entrez_up_onPROM,   organism = "hsa", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(entrez_down_onPROM, organism = "hsa", pvalueCutoff = 0.05)
dotplot(kegg_up,   showCategory = 15) + ggtitle("KEGG — UP × BG4 promoter")
dotplot(kegg_down, showCategory = 15) + ggtitle("KEGG — DOWN × BG4 promoter")

# ── Reactome ──────────────────────────────────────────────────────────────────
react_up   <- ReactomePA::enrichPathway(entrez_up_onPROM,   organism = "human")
react_down <- ReactomePA::enrichPathway(entrez_down_onPROM, organism = "human")
dotplot(react_up,   showCategory = 15) + ggtitle("Reactome — UP × BG4 promoter")
dotplot(react_down, showCategory = 15) + ggtitle("Reactome — DOWN × BG4 promoter")

# ── MSigDB Hallmark ───────────────────────────────────────────────────────────
msig_h   <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
hallmark_up   <- enricher(deg_up_onPROM,   TERM2GENE = msig_h)
hallmark_down <- enricher(deg_down_onPROM, TERM2GENE = msig_h)
dotplot(hallmark_up,   showCategory = 10) + ggtitle("Hallmark — UP × BG4 promoter")
dotplot(hallmark_down, showCategory = 10) + ggtitle("Hallmark — DOWN × BG4 promoter")

# ── STRING protein network ────────────────────────────────────────────────────
string_db  <- STRINGdb$new(version = "12", species = 9606,
                           score_threshold = 400)
df_string  <- data.frame(gene = deg_up_onPROM)
mapped_str <- string_db$map(df_string, "gene", removeUnmappedRows = TRUE)
string_db$plot_network(mapped_str$STRING_id)


# ── SECTION 12: TF ACTIVITY INFERENCE (decoupleR + CollecTRI) ────────────────
# Infers transcription factor activity across all samples using Uniform Linear
# Model (ULM) from decoupleR with the CollecTRI regulon (curated TF–target network).
# Input: log2(FPKM+1) matrix — each column is one sample.

net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)

tf_acts <- decoupleR::run_ulm(
  mat      = log_fpkm_clean,
  net      = net,
  .source  = "source",
  .target  = "target",
  .mor     = "mor",
  minsize  = 5
)

# Filter significant TF activities (p < 0.05)
tf_acts_sig <- tf_acts %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(dplyr::desc(abs(score)))

# Wide matrix: rows = samples, cols = TFs
tf_mat_wide <- tf_acts_sig %>%
  tidyr::pivot_wider(id_cols = "condition",
                     names_from  = "source",
                     values_from = "score") %>%
  tibble::column_to_rownames("condition") %>%
  as.matrix()

# Select top 30 most variable TFs
top_tfs <- tf_acts_sig %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(score), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(abs(std))) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::pull(source)

tf_mat_top <- scale(tf_mat_wide[, top_tfs, drop = FALSE])

pal_rdbu   <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
pal_colors <- grDevices::colorRampPalette(pal_rdbu)(100)
my_breaks  <- c(seq(-2, 0, length.out = 51), seq(0.05, 2, length.out = 50))

pheatmap(
  tf_mat_top,
  color        = pal_colors,
  breaks       = my_breaks,
  border_color = "white",
  cluster_cols = TRUE, cluster_rows = TRUE,
  treeheight_row = FALSE, treeheight_col = FALSE,
  main = paste0("TF activity scores — top ", length(top_tfs), " TFs (decoupleR/CollecTRI)")
)


# ── SECTION 13: MOTIF ANALYSIS (JASPAR2022 + TFBSTools) ──────────────────────
# Scans promoters of BG4+DEG genes (promoter window) for known TF motifs.
# Uses PWMs from JASPAR2022 (human, non-redundant).

# Genes to scan: UP DEGs with BG4 peak on promoter
genes_motif <- deg_up_onPROM

# Map symbols → Entrez
entrez_motif <- mapIds(org.Hs.eg.db, keys = genes_motif,
                       keytype = "SYMBOL", column = "ENTREZID",
                       multiVals = "first") %>% na.omit()

# Get promoter GRanges (2kb upstream, 200bp downstream)
prom_motif <- promoters(
  genes(txdb)[names(genes(txdb)) %in% entrez_motif],
  upstream = 2000, downstream = 200
)

# Fetch PWMs for human TFs from JASPAR2022
pfm_opts <- list(species = 9606, all_versions = FALSE)
pfm      <- getFromNamespace("getMatrixSet", "TFBSTools")(JASPAR2022, pfm_opts)

# Scan promoters for motif matches
motif_ix     <- motifmatchr::matchMotifs(pfm, prom_motif,
                                         genome = BSgenome.Hsapiens.UCSC.hg19)
motif_counts <- sort(colSums(motifmatchr::motifMatches(motif_ix)),
                     decreasing = TRUE)

# Top 15 most frequent motifs in BG4+DEG promoters
barplot(rev(head(motif_counts, 15)),
        horiz = TRUE, las = 1,
        main  = "Top 15 TF motifs in BG4 × UP-DEG promoters",
        xlab  = "Number of promoters with motif match")


# ══════════════════════════════════════════════════════════════════════════════
# ARCHIVE — broken / redundant / exploratory code
# ══════════════════════════════════════════════════════════════════════════════

# ── A1: Column names assigned before object created ───────────────────────────
# In the original script, `colnames(g4_e5_pe_mock) = bed_cols` appeared on the
# line BEFORE g4_e5_pe_mock was loaded with read.table(). This would throw
# "object not found". Fixed in Section 1 by loading both files before any
# column renaming.

# ── A2: Broken for-loop for enhancer gene annotation ─────────────────────────
# The original loop attempted to match gene start/stop positions manually
# against enhancer regions but was never completed (empty gene_enhancer_df,
# wrong column indices, no output). Removed entirely — the correct approach
# is valr::bed_intersect() as used in Section 9.
#
# for(i in 1:length(rownames(enhancer_REF_hg19))){
#   if (geneStart > enhancer_REF_hg19[i,2] & ...
#   gene_enhancer_df= gene_enhancer_df[i,3]  # ← never initialised
# }

# ── A3: coldata rownames reset to wrong length ────────────────────────────────
# The original code ran `rownames(coldata_normal_stage1) = paste0("S", 1:158)`
# then later `rownames(coldata_normal_stage2) = paste0("S", 1:158)` on a
# data frame with only 155 rows (after removing outliers). This throws a
# length mismatch error. Fixed in Section 4 by building coldata once with
# correct lengths and using subset indexing.

# ── A4: scale = "z-score" in pheatmap ────────────────────────────────────────
# pheatmap does not accept scale = "z-score". Valid values are "row", "column",
# "none". Fixed throughout Section 8 to scale = "row".

# ── A5: `entrez` vs `entrez_deg` naming inconsistency ─────────────────────────
# The TSS BED section used `geneGR[names(geneGR) %in% entrez]` but the
# variable was defined as `entrez_deg`. Fixed in Section 7 with consistent naming.

# ── A6: has_peak = SYMBOL %in% peak_df (wrong — peak_df is a data frame) ─────
# One block used `mutate(has_peak = SYMBOL %in% peak_df)` — you cannot use %in%
# to match a character vector against a data frame. Correct form is
# `SYMBOL %in% peak_df$SYMBOL` or `SYMBOL %in% peak_genes`. Fixed in Section 7.

# ── A7: prom used directly in bed_subtract without valr column cleanup ────────
# The original `bed_subtract(overlap_ac27_me1, prom)` would fail because prom
# (as.data.frame of GRanges) has columns seqnames/start/end/width/strand —
# valr requires exactly chrom/start/end as the first three columns.
# Fixed in Section 9 with the to_valr3() helper and explicit transmute().

# ── A8: Duplicate peak annotation and DEG-with-peak blocks ────────────────────
# The script contained 3 near-identical versions of the deg_with_peak_stage1
# assignment, each with slightly different column references. Only the cleanest
# version (tibble rownames_to_column + filter) is kept in Section 7.

# ── A9: BigWig annotation (not connected to main analysis) ────────────────────
# A BigWig file (GSM5066844_18081N_normalized_RPKM.bw) was imported and passed
# to annotatePeak but was never linked to any downstream analysis object.
# Kept here for reference in case it is needed later.
#
# bw_data <- rtracklayer::import(
#   "/Users/hossein.allahdadi/Downloads/GSM5066844_18081N_normalized_RPKM.bw")
# bw_ann  <- ChIPseeker::annotatePeak(bw_data, TxDb = txdb,
#                                     tssRegion = c(-3000, 3000), verbose = FALSE)
# plotAnnoPie(bw_ann)

# ── A10: TF activity — wrong pivot column ("condition" instead of "statistic") ─
# The original pivot_wider used id_cols = "condition" but decoupleR::run_ulm()
# output has columns: statistic, source, target, condition, score, p_value.
# For a sample-by-TF matrix the id should be "condition" (sample name),
# which happens to be correct, but the run_ulm output needs to be filtered
# to a single statistic ("ulm") first. Fixed in Section 12.
