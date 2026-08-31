# ================================================================================
# Figure 5 and Supplementary Figure S6
# CRISPRa-mediated LINE-1 reactivation in osimertinib-resistant cells
# ================================================================================

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
library(scales)
library(ggExtra)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

# -------------------------------------------------------------------------------
# Differential LINE-1 expression after CRISPRa-mediated reactivation
# -------------------------------------------------------------------------------

OR_VPR_count <- read.table("H1975_OR_VPR_merge_geneTE.count", row.names = 1)
colnames(OR_VPR_count) <- c("sgCtrl_1","sgCtrl_2","sgL1_1","sgL1_2")
condition <- factor(rep(c("sgCtrl","sgL1"), c(2, 2)))
colData <- data.frame(row.names = colnames(OR_VPR_count), condition)

OR_VPR_count_dds <- DESeqDataSetFromMatrix(countData = OR_VPR_count,
                                 colData = colData ,
                                 design = ~ condition)

keep <- rowSums(counts(OR_VPR_count_dds)) >= 10
OR_VPR_count_dds <- OR_VPR_count_dds[keep, ]
OR_VPR_count_dds <- DESeq(OR_VPR_count_dds)
scaling_factors <- sizeFactors(OR_VPR_count_dds)

normalized <- as.data.frame(counts(OR_VPR_count_dds, normalized = TRUE))
write.table(normalized, file ="H1975_OR_VPR_diffGeneTE_norm.txt", sep ="\t", quote = FALSE)
OR_VPR_count_res <- results(OR_VPR_count_dds)
write.table(OR_VPR_count_res, file ="H1975_OR_VPR_diffGeneTE.txt", sep ="\t", col.names = TRUE, quote = FALSE)

L1_count <- read.table("H1975_OR_VPR_L16k.count", row.names = 1)
colnames(L1_count) <- c("sgCtrl_1","sgCtrl_2","sgL1_1","sgL1_2")
condition <- factor(rep(c("sgCtrl","sgL1"), c(2, 2)))
colData <- data.frame(row.names = colnames(L1_count), condition)

L1_count_dds <- DESeqDataSetFromMatrix(countData = L1_count,
                                 colData = colData ,
                                 design = ~ condition)
sizeFactors(L1_count_dds) <- scaling_factors

L1_count_dds <- DESeq(L1_count_dds)
normalized <- as.data.frame(counts(L1_count_dds, normalized = TRUE))
write.table(normalized, file ="H1975_OR_VPR_diffL1_genescale_norm.txt", sep ="\t", quote = FALSE)

L1_count_res <- results(L1_count_dds)
write.table(L1_count_res, file ="H1975_OR_VPR_diffL1_genescale.txt", sep ="\t", col.names = TRUE, quote = FALSE)

# -------------------------------------------------------------------------------
# Volcano plot of CRISPRa-responsive LINE-1 loci
# -------------------------------------------------------------------------------

H1975_OR_VPR_L1 <- read.table("H1975_OR_VPR_diffL1_genescale.txt", header = TRUE)
H1975_OR_VPR_L1 <- na.omit(H1975_OR_VPR_L1)
H1975_OR_VPR_L1$class <- ifelse(
  abs(H1975_OR_VPR_L1$log2FoldChange) >= 0.5 & H1975_OR_VPR_L1$padj <= 0.05 ,
  ifelse(H1975_OR_VPR_L1$log2FoldChange >= 0.5 & H1975_OR_VPR_L1$padj <= 0.05,'Up','Down'),'None')
table(H1975_OR_VPR_L1$class)

ggplot(H1975_OR_VPR_L1, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = class), alpha = 0.8, size = 3) +
  scale_color_manual(values = c("royalblue4","grey","firebrick2"), labels = c("Down (2022)","None","Up (35)")) +
  geom_vline(xintercept = c(-0.5, 0.5), lty = 3, color ='black', lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 3, color ='black', lwd = 0.8) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-5, 5, 5)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 60, 20)) +
  theme(axis.title.y = element_text(vjust = 2, size = 28)) +
  theme(axis.text.y = element_text(vjust = 1, size = 25)) +
  theme(axis.title.x = element_text(vjust = 2, size = 28)) +
  theme(axis.text.x = element_text(vjust = 1, size = 25)) +
  theme(legend.text = element_text(size = 22)) +
  labs(x ="Log2FC(sgLINE-1/sgCtrl)", y ="-Log10(padj)") +
  theme(plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle("LINE-1s")

ggsave("H1975_OR_VPR_L1_volcano.pdf", width = 13, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Subfamily distribution of CRISPRa-responsive LINE-1 loci
# -------------------------------------------------------------------------------

OR_VPR_LINE1 <- read.table("H1975_OR_VPR_diffL1_genescale.txt", header = TRUE)
OR_VPR_LINE1$class <- ifelse(
  is.na(OR_VPR_LINE1$padj),"None",
  ifelse(
    abs(OR_VPR_LINE1$log2FoldChange) >= 0.5 & OR_VPR_LINE1$padj <= 0.05,
    ifelse(OR_VPR_LINE1$log2FoldChange >= 0.5 & OR_VPR_LINE1$padj <= 0.05,'Up','Down'),
'None'
  )
)

OR_VPR_LINE1_select <- OR_VPR_LINE1 %>%
  mutate(family = sub("_dup\\d+$","", LINE1)) %>%
  filter(family %in% c("L1HS", paste0("L1PA", 2:7)))

OR_VPR_LINE1_select2 <- OR_VPR_LINE1_select %>%
  group_by(family, class) %>%
  summarise(count = n()) %>%
  group_by(family) %>%
  mutate(percentage = count / sum(count) * 100)

OR_VPR_LINE1_select2$family <- factor(OR_VPR_LINE1_select2$family,
                                      levels = c("L1PA7","L1PA6","L1PA5","L1PA4","L1PA3","L1PA2","L1HS"))

ggplot(OR_VPR_LINE1_select2, aes(x = family, y = percentage, fill = class)) +
  geom_bar(stat ="identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", percentage)),
            position = position_stack(vjust = 0.5), size = 4) +
  labs(x =" ", y ="Percentage", fill ="Category") +
  scale_fill_manual(values = c("Up" ="firebrick2","Down" ="royalblue4","None" ="grey")) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  theme(axis.title.y = element_text(vjust = 0.8, size = 28)) +
  theme(axis.text.y = element_text(size = 28)) +
  theme(axis.title.x = element_text(vjust = 0.8, size = 28)) +
  theme(axis.text.x = element_text(size = 25)) +
  theme(legend.text = element_text(size = 25)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  coord_flip()

ggsave("H1975 OR_VPR L1 family percentage.pdf", width = 11, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Comparison of resistance-associated and CRISPRa-associated LINE-1 changes
# -------------------------------------------------------------------------------

H1975_OR_diffL1 <- read.table("H1975_OR_diffL1_genescale.txt", header = TRUE)
H1975_OR_VPR_diffL1 <- read.table("H1975_OR_VPR_diffL1_genescale.txt", header = TRUE)

H1975_OR_VPR_diffL1_merge <- merge(
  H1975_OR_diffL1,
  H1975_OR_VPR_diffL1,
  by ="LINE1",
  suffixes = c("_OR","_VPR")
)
H1975_OR_VPR_diffL1_merge_logFC <- na.omit(H1975_OR_VPR_diffL1_merge[, c(3, 10)])

cor_result <- cor.test(
  H1975_OR_VPR_diffL1_merge_logFC[, 1],
  H1975_OR_VPR_diffL1_merge_logFC[, 2],
  method ="pearson"
)

pdf("H1975_OR_VPR_diffL1_merge_logFC Scatterplot.pdf", 8, 7)
par(mar = c(5, 5, 3, 3))
palette <- colorRampPalette(c("blue","yellow","red"))
smoothScatter(H1975_OR_VPR_diffL1_merge_logFC,
              colramp = palette,
              main ="LINE-1 loci",
              xlab ="Log2FC (H1975OR/H1975)",
              ylab ="Log2FC (sgLINE-1/sgCtrl)",
              cex.main = 1.5,
              cex.lab = 1.5,
              cex.axis = 1.5,
              font.main = 2,
              font.lab = 2,
              xaxs ="i",
              yaxs ="i")

abline(h = 0, v = 0, col ="gray40", lty = 2, lwd = 1.5)
x_val <- H1975_OR_VPR_diffL1_merge_logFC[, 1]
y_val <- H1975_OR_VPR_diffL1_merge_logFC[, 2]
total_pts <- length(x_val)

if(total_pts > 0) {
  q1_per <- round(sum(x_val > 0 & y_val > 0) / total_pts * 100, 1)
  q2_per <- round(sum(x_val < 0 & y_val > 0) / total_pts * 100, 1)
  q3_per <- round(sum(x_val < 0 & y_val < 0) / total_pts * 100, 1)
  q4_per <- round(sum(x_val > 0 & y_val < 0) / total_pts * 100, 1)

  usr <- par("usr") #
  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]

  text(usr[2] - x_range * 0.08, usr[4] - y_range * 0.06, paste0(q1_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[1] + x_range * 0.08, usr[4] - y_range * 0.06, paste0(q2_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[1] + x_range * 0.08, usr[3] + y_range * 0.06, paste0(q3_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[2] - x_range * 0.08, usr[3] + y_range * 0.06, paste0(q4_per,"%"), cex = 1.4, font = 2, col ="black")
}

dev.off()

# -------------------------------------------------------------------------------
# Neighboring gene expression around reactivated LINE-1 loci
# -------------------------------------------------------------------------------
ORVPR_upL1_geneFC <- read.csv("H1975 OR VPR upL1 neibor diffGene Log2FC.csv", header = TRUE)
ORVPR_upL1_geneFC_rep <- melt(ORVPR_upL1_geneFC)
ORVPR_upL1_geneFC_rep <- na.omit(ORVPR_upL1_geneFC_rep)
group_counts <- ORVPR_upL1_geneFC_rep %>%
  group_by(variable) %>%
  summarize(count = n(), max_value = 6)

ggplot(ORVPR_upL1_geneFC_rep, aes(x = variable, y = value)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, fill ="ghostwhite") +
  geom_hline(yintercept = 0, linetype = 1, color ='red') +
  scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 1, 1)) +
  theme_bw() +theme(panel.grid = element_blank(), legend.position ="none") +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position ="none") +
  theme(axis.title.y = element_text(vjust = 2, size = 28)) +
  theme(axis.text.y = element_text(vjust = 1, size = 25)) +
  theme(axis.title.x = element_text(vjust = 2, size = 28)) +
  theme(axis.text.x = element_text(vjust = 0.6, size = 25, angle = 25)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) +
  labs(title =" ",
       x ="Distance to upregulate LINE-1 (kb)", y ="Gene expression change")

ggsave("gene expression of upregulate LINE1 OR_VPR.pdf", width = 12, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Gene set enrichment analysis after LINE-1 reactivation
# -------------------------------------------------------------------------------

hallmark_gmt <- read.gmt("/RStudio/gmt/human/halmark.all.v2024.1.Hs.symbols.gmt")

OR_VPR_all <- read.table("H1975_OR_VPR_diffGene_symbol.txt", header = TRUE)
OR_VPR_expr <- setNames(OR_VPR_all$stat, OR_VPR_all$symbol)
OR_VPR_expr <- sort(OR_VPR_expr, decreasing = TRUE)
head(OR_VPR_expr)

OR_VPR_hallmark <- GSEA(
  geneList = OR_VPR_expr,
  TERM2GENE = hallmark_gmt,
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod ="BH",
  verbose = FALSE
)
write.csv(OR_VPR_hallmark@result,"OR_VPR_hallmark_gsea.csv", row.names = FALSE)

pdf("OR_VPR_IFN_hallmark_gsea.pdf", 15, 11)
gseaplot2(OR_VPR_hallmark, 1 , pvalue_table = FALSE)
dev.off()

OR_VPR_gsea_data <- read.csv("OR_VPR_hallmark_gsea.csv")
plot_data <- OR_VPR_gsea_data %>%
  mutate(
    log_pvalue = -log10(pvalue),
    color_group = ifelse(NES > 0,"Positive","Negative"),
    bubble_size = abs(NES),
    y_position = 1
  )

ggplot(plot_data, aes(x = NES, y = reorder(ID, NES))) +
  geom_point(aes(
    size = bubble_size,
    color = color_group,
    alpha = log_pvalue)) +
  scale_color_manual(
    values = c("Positive" ="firebrick2","Negative" ="royalblue4"),
    name ="NES Direction") +
  scale_alpha_continuous(
    range = c(0.5, 1),
    name ="-log10(p-value)") +
  scale_size_continuous(
    range = c(4, 8),
    name ="|NES|") +
  labs(
    x ="Normalized Enrichment Score (NES)",
    y ="Pathway",
    title ="GSEA Enrichment Analysis") +
  theme_minimal() +
  theme(
    legend.position ="right",
    plot.title = element_text(hjust = 0.5, face ="bold")) +
  geom_vline(xintercept = 0, linetype ="dashed", color ="gray50")

ggsave("H1975 OR VPR hallmark gsea.pdf", width = 13, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Interferon-response genes near reactivated LINE-1 loci
# -------------------------------------------------------------------------------

# IFN should contain the interferon-response gene symbols used for highlighting.
diff_genes <- read.csv("OR LINE1 VPR diffGene.csv", header = TRUE)
diff_genes <- diff_genes %>% filter(baseMean > 1)
distance_data <- read.table("gene_nearest_L1_distance_full.txt", header = TRUE)

diff_genes_clean <- diff_genes %>%
  select(symbol, ENSEMBL, log2FoldChange, padj, baseMean) %>%
  rename(GeneName = symbol, ensembl_id = ENSEMBL)

combined_data <- distance_data %>%
  inner_join(diff_genes_clean, by ="ensembl_id") %>%
  filter(distance != -1) %>%
  mutate(
    group = ifelse(GeneName %in% IFN,
"IFN Genes","Other Genes")
  )

label_data <- combined_data %>%
  filter(group =="IFN Genes") %>%
  group_by(GeneName) %>%
  slice(1) %>%
  ungroup()

main_plot <- combined_data %>%
  arrange(desc(group)) %>%
  ggplot(aes(x = distance, y = log2FoldChange,
             color = group, group = group,
             alpha = group, size = group)) +
  geom_point(stroke = 0.1) +
  scale_alpha_manual(values = c("IFN Genes" = 1,"Other Genes" = 0.3)) +
  scale_size_manual(values = c("IFN Genes" = 3,"Other Genes" = 0.8)) +
  scale_color_manual(values = c("IFN Genes" ="#A83030","Other Genes" ="grey")) +

  geom_text_repel(
    data = label_data,
    aes(label = GeneName),
    size = 5,
    min.segment.length = 0,
    segment.size = 0.3,
    segment.color ="black",
    segment.alpha = 0.5,
    box.padding = 0.5,
    max.overlaps = 20,
    force = 1,
    direction ="both") +

 geom_vline(xintercept = -500000, color ="black", alpha = 0.5, linetype ="dashed") +
 geom_vline(xintercept = 500000, color ="black", alpha = 0.5, linetype ="dashed") +
 geom_hline(yintercept = 0, color ="black", alpha = 0.5, linetype ="dashed") +

labs(
    x ="Distance to nearest full-length LINE-1 (Mb)",
    y ="Gene Expression Changes
        log2FC(sgLINE-1/sgCtrl)",
    title =" ") +

scale_y_continuous(limits = c(-1.5, 3)) +
scale_x_continuous(
    limits = c(-4000000, 4000000),
    breaks = seq(-4000000, 4000000, 1000000),
    labels = function(x) x / 1000000
  ) +

theme_bw() +
theme(
    plot.title = element_text(hjust = 0.5, face ="bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position ="bottom",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color ="black"),
    axis.title = element_text(color ="black"),
    panel.border = element_rect(color ="black", fill = NA, linewidth = 0.5),
    plot.margin = margin(1, 1, 1, 1,"cm")
  )

final_plot <- ggMarginal(
  main_plot,
  type ="density",
  margins ="x",
  groupColour = TRUE,
  groupFill = TRUE,
  size = 4)

ggsave("L1_distance_vs_expression_scatter.pdf", final_plot, width = 10, height = 10, dpi = 300)

ifn_boxplot <- combined_data %>%
  arrange(desc(group)) %>%
  ggplot(aes(x = group, y = log2FoldChange, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("IFN Genes" ="#A83030","Other Genes" ="grey70")) +
  geom_jitter(data = filter(combined_data, group =="IFN Genes"),
              width = 0.2, size = 1, color ="#A83030", alpha = 0.5) +

  stat_compare_means(
    comparisons = list(c("IFN Genes","Other Genes")),
    method ="wilcox.test",
    label ="p.format",
    size = 3
  ) +

  labs(x ="", y ="log2 Fold Change") +
  theme_minimal() +
  theme(
    legend.position ="none",
    axis.text.x = element_text(angle = 45, hjust = 1, color ="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color ="black"),
    axis.ticks = element_line(color ="black")
  )

ggsave("IFN gene expression boxplot.pdf", ifn_boxplot, width = 4, height = 5)

IFN_genes_exp <- read.csv("IFN_genes_exp.csv", header = TRUE)
rownames(IFN_genes_exp) <- IFN_genes_exp$symbol
IFN_genes_exp <- IFN_genes_exp[, 15:18]
IFN_genes_exp_scale <- t(apply(IFN_genes_exp, 1, scale))
colnames(IFN_genes_exp_scale) <- colnames(IFN_genes_exp)

set.seed(1234)
pdf("H1975 OR VPR IFN gene heatmap.pdf", 6, 8)
exp_ht <- Heatmap(IFN_genes_exp_scale,
                  cluster_rows = TRUE ,
                  cluster_columns = FALSE,
                  column_names_gp = gpar(fontsize = 14),
                  column_names_rot = 60,
                  show_row_names = TRUE,
                  width = ncol(IFN_genes_exp_scale)*grid::unit(13.5,"mm"),
                  height = nrow(IFN_genes_exp_scale)*grid::unit(6,"mm"),
                  col = colorRamp2(c(-1, 0, 1), c("royalblue4","white","firebrick2")),
                  heatmap_legend_param = list(
                    title ="Expression Level
                    Row Z-score",
                    legend_width = grid::unit(2,"cm"),
                    at = c(-1, 1), labels = c("-1","1"),
                    title_position ="lefttop-rot"),
)

draw(exp_ht, heatmap_legend_side ="left")

dev.off()
