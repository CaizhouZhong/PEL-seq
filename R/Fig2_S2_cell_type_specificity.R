# ================================================================================
# Figure 2 and Supplementary Figure S2
# Cell-type-selective LINE-1 transcription and neighboring gene associations
# ================================================================================

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(reshape2)
library(ggpubr)

gene_count <- read.table("Celltye_merge_geneTE.count", row.names = 1)
colnames(gene_count) <- c("BEAS2B_r1","BEAS2B_r2","H1975_r1","H1975_r2","PC9_r1","PC9_r2",
"HCC827_r1","HCC827_r2","H1_r1","H1_r2","H9_r1","H9_r2")
condition <- factor(rep(c("BEAS2B","H1975","PC9","HCC827","H1","H9"), c(2, 2, 2, 2, 2, 2)))
colData <- data.frame(row.names = colnames(gene_count), condition)
gene_count_dds <- DESeqDataSetFromMatrix(countData = gene_count,
                                         colData = colData ,
                                         design = ~ condition)
keep <- rowSums(counts(gene_count_dds)) >= 10
gene_count_dds <- gene_count_dds[keep, ]
gene_count_dds <- DESeq(gene_count_dds)
cell_size_factors <- sizeFactors(gene_count_dds)

# Apply size factors derived from the combined gene/TE count matrix.
L1_count <- read.table("merge_L16k.count", row.names = 1)
colnames(L1_count) <- c("BEAS2B_r1","BEAS2B_r2","H1975_r1","H1975_r2","PC9_r1","PC9_r2",
"HCC827_r1","HCC827_r2","H1_r1","H1_r2","H9_r1","H9_r2")
condition <- factor(rep(c("BEAS2B","H1975","PC9","HCC827","H1","H9"), c(2, 2, 2, 2, 2, 2)))
colData <- data.frame(row.names = colnames(L1_count), condition)

scaling_factors <- cell_size_factors
L1_count_dds <- DESeqDataSetFromMatrix(countData = L1_count,
                                 colData = colData ,
                                 design = ~ condition)
sizeFactors(L1_count_dds) <- scaling_factors
L1_count_dds <- DESeq(L1_count_dds)
sizeFactors(L1_count_dds)
normalized <- as.data.frame(counts(L1_count_dds, normalized = TRUE))
write.table(normalized, file ="Nor_cancer_hESC_L16k_norm_genescale.txt", sep ="\t", quote = FALSE)

L16k_B_H1975 <- results(L1_count_dds, contrast = c("condition","H1975","BEAS2B"))
write.table(L16k_B_H1975, file ="BEAS2B_H1975_diffL1_genescale.txt", sep ="\t", col.names = TRUE, quote = FALSE)

L16k_B_PC9 <- results(L1_count_dds, contrast = c("condition","PC9","BEAS2B"))
write.table(L16k_B_PC9, file ="BEAS2B_PC9_diffL1_genescale.txt", sep ="\t", col.names = TRUE, quote = FALSE)

L16k_B_HCC827 <- results(L1_count_dds, contrast = c("condition","HCC827","BEAS2B"))
write.table(L16k_B_HCC827, file ="BEAS2B_HCC827_diffL1_genescale.txt", sep ="\t", col.names = TRUE, quote = FALSE)

L16k_B_H1 <- results(L1_count_dds, contrast = c("condition","H1","BEAS2B"))
write.table(L16k_B_H1, file ="BEAS2B_H1_diffL1_genescale.txt", sep ="\t", col.names = TRUE, quote = FALSE)

L16k_B_H9 <- results(L1_count_dds, contrast = c("condition","H9","BEAS2B"))
write.table(L16k_B_H9, file ="BEAS2B_H9_diffL1_genescale.txt", sep ="\t", col.names = TRUE, quote = FALSE)

# -------------------------------------------------------------------------------
# PCA of locus-specific LINE-1 expression
# -------------------------------------------------------------------------------
vsd <- vst(L1_count_dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup ="condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData,"percentVar"))

pdf("L1_count_dds_vsd_PCA.pdf", width = 7.5, height = 5)

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 7) +
  scale_color_manual(
    values = c(
"BEAS2B" ="steelblue",
"H1975" ="maroon3",
"PC9" ="coral",
"HCC827" ="red",
"H1" ="olivedrab",
"H9" ="seagreen2"
    )
  ) +
  xlab(paste0("PC1: ", percentVar[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar[2],"% variance")) +
  theme_bw() +
  theme(
    legend.title = element_blank()
  )

dev.off()

# -------------------------------------------------------------------------------
# Cell-type-specific LINE-1 expression heatmap
# -------------------------------------------------------------------------------

mergeLINE1_norm <- read.table("Nor_cancer_hESC_L16k_norm_genescale.txt", header = TRUE, row.names = 1)
row_totals <- rowSums(mergeLINE1_norm)
cutoff <- quantile(row_totals, probs = 0.30)
mergeLINE1_norm <- mergeLINE1_norm[row_totals >= cutoff, ]
mergeLINE1_norm_scale <- t(apply(mergeLINE1_norm, 1, scale))
colnames(mergeLINE1_norm_scale) <- colnames(mergeLINE1_norm)

set.seed(1234)
wss <- numeric(10)
for (k in 1:10) {
  wss[k] <- sum(kmeans(mergeLINE1_norm_scale, centers = k)$withinss)
}

plot(1:10, wss, type ="b",
     xlab ="Number of Clusters (k)",
     ylab ="Total Within-cluster SS",
     main ="Elbow Method for Choosing k")

for (k in 2:8) {
  wss[k] <- sum(kmeans(mergeLINE1_norm_scale, centers = k, nstart = 50)$withinss)
}
plot(2:8, wss[2:8], type ="b")

set.seed(1234)
pdf("BEAS2B cancer hESCs locus_LINE1 heatmap.pdf", 10, 12)
exp_ht <- Heatmap(mergeLINE1_norm_scale,
                  cluster_rows = FALSE ,
                  cluster_columns = FALSE,
                  column_names_gp = gpar(fontsize = 14),
                  column_names_rot = 60,
                  show_row_names = FALSE,
                  km = 5,
                  width = ncol(mergeLINE1_norm_scale)*grid::unit(12,"mm"),
                  height = nrow(mergeLINE1_norm_scale)*grid::unit(0.02,"mm"),
                  col = colorRamp2(c(-1, 0, 1), c("royalblue4","white","firebrick2")),
                  heatmap_legend_param = list(
                    title ="Expression Level
                    Row Z-score",
                    legend_width = grid::unit(2,"cm"),
                    at = c(-1, 1), labels = c("-1","1"),
                    title_position ="lefttop-rot"),
)

exp_ht_drawn <- draw(exp_ht, heatmap_legend_side ="left")
dev.off()

row_clusters <- row_order(exp_ht_drawn)
cluster_L1 <- lapply(row_clusters, function(idx) rownames(mergeLINE1_norm_scale)[idx])
names(cluster_L1) <- paste0("Cluster_", names(cluster_L1))
sapply(cluster_L1, length)

output_dir <-"L1_cluster"
dir.create(output_dir, showWarnings = FALSE)

for (name in names(cluster_L1)) {
  file_path <- file.path(output_dir, paste0(name,".txt"))
  writeLines(cluster_L1[[name]], file_path)
}

# -------------------------------------------------------------------------------
# Intersections of upregulated LINE-1 loci across cell types
# -------------------------------------------------------------------------------

Beas2b_H1975_L16k <- read.table("BEAS2B_H1975_diffL1_genescale.txt", header = TRUE)
Beas2b_H1975_L16k$class <- ifelse(
  abs(Beas2b_H1975_L16k$log2FoldChange) >= 1 & Beas2b_H1975_L16k$padj <= 0.05 ,
  ifelse(Beas2b_H1975_L16k$log2FoldChange >= 1 & Beas2b_H1975_L16k$padj <= 0.05,'Up','Down'),'None')
table(Beas2b_H1975_L16k$class)
Beas2b_H1975_L16k_up <- subset(Beas2b_H1975_L16k, class%in%'Up')

Beas2b_PC9_L16k <- read.table("BEAS2B_PC9_diffL1_genescale.txt", header = TRUE)
Beas2b_PC9_L16k$class <- ifelse(
  abs(Beas2b_PC9_L16k$log2FoldChange) >= 1 & Beas2b_PC9_L16k$padj <= 0.05 ,
  ifelse(Beas2b_PC9_L16k$log2FoldChange >= 1 & Beas2b_PC9_L16k$padj <= 0.05,'Up','Down'),'None')
table(Beas2b_PC9_L16k$class)
Beas2b_PC9_L16k_up <- subset(Beas2b_PC9_L16k, class%in%'Up')

Beas2b_HCC827_L16k <- read.table("BEAS2B_HCC827_diffL1_genescale.txt", header = TRUE)
Beas2b_HCC827_L16k$class <- ifelse(
  abs(Beas2b_HCC827_L16k$log2FoldChange) >= 1 & Beas2b_HCC827_L16k$padj <= 0.05 ,
  ifelse(Beas2b_HCC827_L16k$log2FoldChange >= 1 & Beas2b_HCC827_L16k$padj <= 0.05,'Up','Down'),'None')
table(Beas2b_HCC827_L16k$class)
Beas2b_HCC827_L16k_up <- subset(Beas2b_HCC827_L16k, class%in%'Up')

Beas2b_H1_L16k <- read.table("BEAS2B_H1_diffL1_genescale.txt", header = TRUE)
Beas2b_H1_L16k$class <- ifelse(
  abs(Beas2b_H1_L16k$log2FoldChange) >= 1 & Beas2b_H1_L16k$padj <= 0.05 ,
  ifelse(Beas2b_H1_L16k$log2FoldChange >= 1 & Beas2b_H1_L16k$padj <= 0.05,'Up','Down'),'None')
table(Beas2b_H1_L16k$class)
Beas2b_H1_L16k_up <- subset(Beas2b_H1_L16k, class%in%'Up')

Beas2b_H9_L16k <- read.table("BEAS2B_H9_diffL1_genescale.txt", header = TRUE)
Beas2b_H9_L16k$class <- ifelse(
  abs(Beas2b_H9_L16k$log2FoldChange) >= 1 & Beas2b_H9_L16k$padj <= 0.05 ,
  ifelse(Beas2b_H9_L16k$log2FoldChange >= 1 & Beas2b_H9_L16k$padj <= 0.05,'Up','Down'),'None')
table(Beas2b_H9_L16k$class)
Beas2b_H9_L16k_up <- subset(Beas2b_H9_L16k, class%in%'Up')

listInput <- list(
  H1975_L16k_up = Beas2b_H1975_L16k_up$LINE1,
  PC9_L16k_up = Beas2b_PC9_L16k_up$LINE1,
  HCC827_L16k_up = Beas2b_HCC827_L16k_up$LINE1,
  H1_L16k_up = Beas2b_H1_L16k_up$LINE1,
  H9_L16k_up = Beas2b_H9_L16k_up$LINE1)

pdf("BEAS2B cancer hESCs uprelate L1 upsetR.pdf", 15, 10)
upset(fromList(listInput), order.by ="freq",
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label ="Intersect LINE-1 loci",
      sets.x.label ="Number of LINE-1 loci",
      text.scale = c(3, 4, 3, 2, 3, 4),
      queries = list(
        list(query = intersects,
             params = list("H1975_L16k_up","PC9_L16k_up","HCC827_L16k_up"),
             color ="red3",
             active = TRUE),
        list(query = intersects,
             params = list("H1_L16k_up","H9_L16k_up"),
             color ="springgreen",
             active = TRUE)))
dev.off()

# -------------------------------------------------------------------------------
# Total full-length LINE-1 abundance across cell types
# -------------------------------------------------------------------------------
total_line1_raw <- read.csv(
"celtype_norm_count.csv",
  header = TRUE,
  check.names = FALSE
)

if (all(c("Sample","Value") %in% colnames(total_line1_raw))) {
  total_LINE1 <- total_line1_raw[, c("Sample","Value")]
} else if (nrow(total_line1_raw) == 1) {
  total_LINE1 <- data.frame(
    Sample = colnames(total_line1_raw),
    Value = as.numeric(total_line1_raw[1, ])
  )
} else {
  stop("celtype_norm_count.csv must contain Sample/Value columns or one row of sample values.")
}

total_LINE1$Group <- sub("-\\d+$","", total_LINE1$Sample)
total_LINE1$Group <- factor(total_LINE1$Group, levels = c("BEAS2B","H1975","PC9","HCC827","H1","H9"))
print(total_LINE1)

ggplot(total_LINE1, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  stat_boxplot(geom ="errorbar", width = 0.2) +
  scale_fill_manual(values = c("BEAS2B" ="steelblue","H1975" ="maroon3","PC9" ="coral",
"HCC827" ="red","H1" ="olivedrab","H9" ="seagreen2")) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_y_continuous(
    limits = c(100000, 500000),
    labels = scales::scientific) +
  theme(axis.title.y = element_text(vjust = 0.8, size = 25)) +
  theme(axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = 0.8, size = 25)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(legend.position ="none") +
  labs(x =" ", y ="DESeq2 Normalized Counts") +
  ggtitle(" ")
ggsave("BEAS2B cancer hESCs FL LINE1.pdf", width = 13, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Number of upregulated LINE-1 loci
# -------------------------------------------------------------------------------
upL1_count <- data.frame(c(677, 625, 1386, 5486, 5254))
upL1_count$group <- c("H1975","PC9","HCC827","H1","H9")
colnames(upL1_count) <- c("upL1","group")
upL1_count$group <- factor(upL1_count$group, levels = c("H1975","PC9","HCC827","H1","H9"))

ggplot(upL1_count, aes(x = group, y = upL1, fill = group)) +
  geom_bar(stat ="identity", width = 0.7, alpha = 0.8) +

  geom_text(
    aes(label = upL1),
    vjust = -0.5,
    size = 6,
    fontface ="bold"
  ) +

  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(
    limits = c(0, 6000),
    expand = expansion(mult = c(0, 0.1)),
    breaks = seq(0, 6000, 2000)
  ) +
  scale_fill_manual(values = c("H1975" ="indianred","PC9" ="coral2",
"HCC827" ="firebrick1","H1" ="olivedrab","H9" ="seagreen")) +
  labs(title =" ", x ="", y ="Upregulated LINE-1 Loci") +
  theme(axis.title.y = element_text(vjust = 0.8, size = 25),
        axis.title.x = element_text(vjust = 0.8, size = 25),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.position ="none")
ggsave("BEAS2B cancer hESCs uprelate L1 loci.pdf", width = 13, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Functional enrichment of LINE-1-associated clusters
# -------------------------------------------------------------------------------

cluster125_L1_enrichment <- read.csv("cluster125 L1 terms enrichment.csv")
cluster125_L1_enrichment$annotation_x <- -5
ggplot(cluster125_L1_enrichment, aes(x = reorder(cluster125_L1_enrichment$Description, abs(LogP)), y = -(LogP))) +
  geom_bar(stat ="identity", width = 0.7, fill ="tomato2") +
  coord_flip() +
  geom_text(aes(label = Description, y = annotation_x), size = 7) +
  labs(x = NULL, y ="LogP", title =" ") +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 5)) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 25),
        legend.position ="none",
        panel.grid = element_blank(),
        strip.text = element_blank(),
        panel.spacing = grid::unit(0.3,"cm")  )
ggsave("cluster125_L1_enrichment.pdf", width = 15, height = 10, dpi = 300)

cluster34_L1_enrichment <- read.csv("cluster34 L1 terms enrichment.csv")
cluster34_L1_enrichment$annotation_x <- -5
ggplot(cluster34_L1_enrichment, aes(x = reorder(cluster34_L1_enrichment$Description, abs(LogP)), y = -(LogP))) +
  geom_bar(stat ="identity", width = 0.7, fill ="tomato2") +
  coord_flip() +
  geom_text(aes(label = Description, y = annotation_x), size = 7) +
  labs(x = NULL, y ="LogP", title =" ") +
  scale_y_continuous(limits = c(-10, 8), breaks = seq(-5, 5, 5)) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 25),
        legend.position ="none",
        panel.grid = element_blank(),
        strip.text = element_blank(),
        panel.spacing = grid::unit(0.3,"cm")  )
ggsave("cluster34_L1_enrichment.pdf", width = 15, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Correlation between LINE-1 and nearby gene expression
# -------------------------------------------------------------------------------
L1_tss2kb_gene <- read.table("Gene_L16k_tss2kb.txt")
colnames(L1_tss2kb_gene) <- c("LINE1","symbol","ENSEMBL")
L1_tss2kb_gene_L1 <- B_H1975_L16k[match(L1_tss2kb_gene$LINE1, B_H1975_L16k$LINE1), ]
L1_tss2kb_gene_gene <- B_H1975_DEG[match(L1_tss2kb_gene$ENSEMBL, B_H1975_DEG$ENSEMBL), ]
L1_tss2kb_gene_L1gene_merge <- cbind(L1_tss2kb_gene_L1, L1_tss2kb_gene_gene)
L1_tss2kb_gene_L1gene_merge <- na.omit(L1_tss2kb_gene_L1gene_merge)
write.csv(L1_tss2kb_gene_L1gene_merge, file ="B_H1975_L1_tss2kb_gene_merge.csv", quote = FALSE, row.names = FALSE)
L1_tss2kb_gene_L1gene_logFC <- L1_tss2kb_gene_L1gene_merge[, c(3, 12)]
L1_tss2kb_gene_L1gene_logFC <- na.omit(L1_tss2kb_gene_L1gene_logFC)
cor(L1_tss2kb_gene_L1gene_logFC)

pdf("B_H1975_L1_tss2kb_gene cor.pdf", 8, 8)
par(mar = c(5, 5, 3, 3))
palette <- colorRampPalette(c("blue","yellow","red"))
smoothScatter(L1_tss2kb_gene_L1gene_logFC,
              colramp = palette,
              main ="H1975/BEAS2B",
              xlab ="Log2FC of LINE1 ",
              ylab ="Log2FC of Gene",
              cex.main = 1.5,
              cex.lab = 1.5,
              cex.axis = 1.5,
              font.main = 2,
              font.lab = 2,
              xaxs ="i",
              yaxs ="i")

abline(h = 0, v = 0, col ="gray40", lty = 2, lwd = 1.5)

if(ncol(L1_tss2kb_gene_L1gene_logFC) >= 2) {
  x <- L1_tss2kb_gene_L1gene_logFC[, 1]
  y <- L1_tss2kb_gene_L1gene_logFC[, 2]

  cor_test <- cor.test(x, y, method ="pearson", use ="complete.obs")
  rho <- round(cor_test$estimate, 3)
  p_value <- round(cor_test$p.value, 4)

  legend("topleft",
         legend = c(paste0("R = ", rho),
                    paste0("p = ", ifelse(p_value < 0.0001,"< 0.0001", p_value))),
         bty ="n",
         cex = 1.5,
         text.col ="darkred")
}

dev.off()

# -------------------------------------------------------------------------------
# Distance-dependent correlation between LINE-1 and neighboring genes
# -------------------------------------------------------------------------------

distances <- c("2kb","100kb","200kb","300kb","400kb","500kb")

analyze_distance_correlation <- function(l1_diff, gene_diff, prefix, summary_file) {
  results <- purrr::map_dfr(distances, function(dist) {
    message("Processing distance: ", dist)

    l1_gene_data <- read.table(
      paste0("Gene_L16k_tss", dist,".txt"),
      header = FALSE
    )
    colnames(l1_gene_data) <- c("LINE1","symbol","ENSEMBL")

    l1_matched <- l1_diff[match(l1_gene_data$LINE1, l1_diff$LINE1), ]
    gene_matched <- gene_diff[match(l1_gene_data$ENSEMBL, gene_diff$ENSEMBL), ]
    merged_data <- na.omit(cbind(l1_matched, gene_matched))

    write.csv(
      merged_data,
      file = paste0(prefix,"_L1_tss", dist,"_gene_merge.csv"),
      quote = FALSE,
      row.names = FALSE
    )

    log_fc_data <- na.omit(merged_data[, c(3, 12)])
    if (nrow(log_fc_data) < 3) {
      return(data.frame(Distance = dist, Pearson_R = NA_real_, P_value = NA_real_))
    }

    cor_test <- cor.test(
      log_fc_data[, 1],
      log_fc_data[, 2],
      method ="pearson"
    )

    data.frame(
      Distance = dist,
      Pearson_R = unname(cor_test$estimate),
      P_value = cor_test$p.value
    )
  })

  write.csv(results, summary_file, row.names = FALSE)
  results
}

B_H1975_L16k <- read.table("BEAS2B_H1975_diffL1_genescale.txt", header = TRUE)
B_H1975_DEG <- read.table("Lung_cell_BEAS2B_H1975_DEG_symbol.txt", header = TRUE)
B_PC9_L16k <- read.table("BEAS2B_PC9_diffL1_genescale.txt", header = TRUE)
B_PC9_DEG <- read.table("Lung_cell_BEAS2B_PC9_DEG_symbol.txt", header = TRUE)
B_HCC827_L16k <- read.table("BEAS2B_HCC827_diffL1_genescale.txt", header = TRUE)
B_HCC827_DEG <- read.table("Lung_cell_BEAS2B_HCC827_DEG_symbol.txt", header = TRUE)

df_h1975 <- analyze_distance_correlation(
  B_H1975_L16k,
  B_H1975_DEG,
"B_H1975",
"B_H1975 Summary_L1_Gene_Correlations.csv"
)

df_pc9 <- analyze_distance_correlation(
  B_PC9_L16k,
  B_PC9_DEG,
"B_PC9",
"B_PC9 Summary_L1_Gene_Correlations.csv"
)

df_hcc827 <- analyze_distance_correlation(
  B_HCC827_L16k,
  B_HCC827_DEG,
"B_HCC827",
"B_HCC827 Summary_L1_Gene_Correlations.csv"
)

plot_data <- bind_rows(
  mutate(df_h1975, CellLine ="H1975"),
  mutate(df_pc9, CellLine ="PC9"),
  mutate(df_hcc827, CellLine ="HCC827")
)

plot_data$Distance <- factor(plot_data$Distance, levels = distances)
plot_data$CellLine <- factor(
  plot_data$CellLine,
  levels = c("H1975","PC9","HCC827")
)

ggplot(
  plot_data,
  aes(x = Distance, y = Pearson_R, group = CellLine, color = CellLine)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4.5) +
  scale_y_continuous(
    limits = c(0, max(plot_data$Pearson_R, na.rm = TRUE) * 1.15)
  ) +
  scale_color_manual(
    values = c(
"H1975" ="maroon3",
"PC9" ="coral",
"HCC827" ="red"
    )
  ) +
  labs(
    x ="Genomic Distance to LINE-1 TSS",
    y ="Pearson Correlation (R)",
    color ="Cell Line"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 13, face ="bold"),
    axis.text = element_text(size = 12, color ="black"),
    legend.title = element_text(size = 12, face ="bold"),
    legend.text = element_text(size = 11),
    legend.position ="top"
  )

ggsave("LINE1 Gene Cor Distance Normal cancer.pdf", width = 7, height = 5)

# -------------------------------------------------------------------------------
# Validation in public datasets
# -------------------------------------------------------------------------------
Public_geneTE <- read.table("Public_tumor_merge_geneTE.count", row.names = 1)
condition <- factor(rep(c("control","tumor"), c(14, 14)))
colData <- data.frame(row.names = colnames(Public_geneTE), condition)
Public_geneTE_dds <- DESeqDataSetFromMatrix(countData = Public_geneTE,
                                         colData = colData ,
                                         design = ~ condition)
keep <- rowSums(counts(Public_geneTE_dds)) >= 20
Public_geneTE_dds <- Public_geneTE_dds[keep, ]
Public_geneTE_dds <- DESeq(Public_geneTE_dds)
sizeFactors(Public_geneTE_dds)
geneTE_factor <- sizeFactors(Public_geneTE_dds)

normalized <- as.data.frame(counts(Public_geneTE_dds, normalized = TRUE))
write.table(normalized, file ="Public_geneTE_all_norm.txt", sep ="\t", quote = FALSE)
Public_geneTE_res <- results(Public_geneTE_dds)
write.table(Public_geneTE_res, file ="Public_Tumor_geneTE_diff.txt", sep ="\t", col.names = TRUE, quote = FALSE)

PublicL1_count <- read.table("Public_L16k.count", row.names = 1)
condition <- factor(rep(c("control","tumor"), c(14, 14)))
colData <- data.frame(row.names = colnames(PublicL1_count), condition)

scaling_factors <- geneTE_factor
PublicL1_count_dds <- DESeqDataSetFromMatrix(countData = PublicL1_count,
                                 colData = colData ,
                                 design = ~ condition)
sizeFactors(PublicL1_count_dds) <- scaling_factors
PublicL1_count_dds <- DESeq(PublicL1_count_dds)
sizeFactors(PublicL1_count_dds)

normalized <- as.data.frame(counts(PublicL1_count_dds, normalized = TRUE))
write.table(normalized, file ="Public_L16k_norm_genescale.txt", sep ="\t", quote = FALSE)
PublicL1_res <- results(PublicL1_count_dds)
write.table(PublicL1_res, file ="Public_Tumor_L16k_diff.txt", sep ="\t", col.names = TRUE, quote = FALSE)

Norm_L16k <- read.csv("Public_Tumor_L16k_sum.csv", header = TRUE)
Norm_L16k <- data.frame(Norm_L16k)
Norm_L16k$group <- rep(c("Normal","Tumor"), c(14, 14))
ggplot(Norm_L16k,
       aes(x = group, y = Norm_L16k, fill = group)) +
  stat_boxplot(geom ="errorbar", width = 0.2) +
  geom_boxplot(width = 0.4) +
  scale_fill_manual(values = c("Normal" ="steelblue","Tumor" ="indianred")) +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.y = element_text(size = 38)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 28)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), legend.position ="none") +
  scale_y_continuous(limits = c(9000, 36000), breaks = seq(10000, 30000, 10000)) +
  labs(x =" ", y ="Normalized counts")
ggsave("GSE81089 L16k Normalized DESeq2 count .pdf", width = 12, height = 10, dpi = 300)

Norm_L16k_normal <- subset(Norm_L16k, group %in%"Normal")
Norm_L16k_tumor <- subset(Norm_L16k, group %in%"Tumor")
t.test(Norm_L16k_normal$Norm_L16k, Norm_L16k_tumor$Norm_L16k)

GSE81089_L1_norm <- read.table("Public_L16k_norm_genescale.txt", header = TRUE)
GSE81089_L15534 <- subset(GSE81089_L1_norm, LINE1 %in%"L1PA3_dup5534")
GSE81089_L15534_melt <- melt(GSE81089_L15534)
GSE81089_L15534_melt$group <- rep(c("Normal","Tumor"), c(14, 14))
ggplot(GSE81089_L15534_melt,
       aes(x = group, y = value, fill = group)) +
  stat_boxplot(geom ="errorbar", width = 0.2) +
  geom_boxplot(width = 0.4) +
  scale_fill_manual(values = c("Tumor" ="indianred","Normal" ="aquamarine3")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust = 0.8, size = 30)) +
  theme(axis.text.y = element_text(vjust = 1, size = 28)) +
  theme(axis.title.x = element_text(vjust = 0.8, size = 30)) +
  theme(axis.text.x = element_text(size = 28)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), legend.position ="none") +
  labs(x =" ", y ="Normalized Counts")
ggsave("GSE81089 LINE1_dup5534 Normalized DESeq2 Counts.pdf", width = 12, height = 10, dpi = 300)

GSE81089_ZFPM2_norm <- read.table("Public_Tumor_geneTE_diff.txt", header = TRUE)
GSE81089_ZFPM2 <- subset(GSE81089_ZFPM2_norm, gene %in% c("ZFPM2-AS1"))
GSE81089_ZFPM2_melt <- melt(GSE81089_ZFPM2)
GSE81089_ZFPM2_melt$group <- rep(c("Normal","Tumor"), c(14, 14))
ggplot(GSE81089_ZFPM2_melt,
       aes(x = group, y = value, fill = group)) +
  stat_boxplot(geom ="errorbar", width = 0.2) +
  geom_boxplot(width = 0.4) +
  scale_fill_manual(values = c("Tumor" ="indianred","Normal" ="aquamarine3")) +
  theme_bw() +
  theme(axis.title.y = element_text(vjust = 0.8, size = 30)) +
  theme(axis.text.y = element_text(vjust = 1, size = 28)) +
  theme(axis.title.x = element_text(vjust = 0.8, size = 30)) +
  theme(axis.text.x = element_text(size = 28)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28), legend.position ="none") +
  labs(x =" ", y ="Normalized Counts")
ggsave("GSE81089 ZFPM2-AS1 Normalized DESeq2 Counts.pdf", width = 12, height = 10, dpi = 300)

ZFPM2_L1_norm <- subset(GSE81089_L1_norm, LINE1 %in% c("L1PA3_dup5534"))
rownames(ZFPM2_L1_norm) <- ZFPM2_L1_norm$LINE1
ZFPM2_L1_norm_t <- t(ZFPM2_L1_norm)
ZFPM2_L1_norm_t <- ZFPM2_L1_norm_t[-1, ]
ZFPM2_L1_norm_t <- apply(ZFPM2_L1_norm_t, 2, as.numeric)
ZFPM2_L1_norm_merge <- cbind(ZFPM2_L1_norm_t, GSE81089_ZFPM2_melt)

ggplot(ZFPM2_L1_norm_merge, aes(x = ZFPM2_L1_norm_merge$L1PA3_dup5534, y = ZFPM2_L1_norm_merge$value)) +
  geom_point(aes(color = group), size = 3) +
  theme_bw() +
  scale_color_manual(values = c("Tumor" ="indianred","Normal" ="aquamarine3")) +
  theme(axis.title.y = element_text(vjust = 0.8, size = 28)) +
  theme(axis.text.y = element_text(size = 26)) +
  theme(axis.title.x = element_text(vjust = 0.8, size = 28)) +
  theme(axis.text.x = element_text(size = 26)) +
  theme(plot.title = element_text(hjust = 0.5, size = 30)) +
  labs( x ="L1PA3_dup5534 norm",
        y ="ZFPM2-AS1 norm") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 30)) +
  geom_smooth(method ='lm', color ='steelblue', fill ='lightgrey') +
  stat_cor(method ="spearman", size = 8,
           label.x = 0,
           label.y = 400 )
ggsave("Correlation between ZFPM2-AS1 and L1PA3_dup5534.pdf", width = 12, height = 10, dpi = 300)
