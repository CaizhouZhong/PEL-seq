# ================================================================================
# Figure 4 and Supplementary Figures S4-S5
# LINE-1 remodeling during acquisition of osimertinib resistance
# ================================================================================

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(reshape2)
library(ggpubr)
library(ggalluvial)
library(UpSetR)

# -------------------------------------------------------------------------------
# Differential Gene and LINE-1 expression in H1975 osimertinib-resistant cells
# -------------------------------------------------------------------------------

OR_count <- read.table("H1975_OR_merge_geneTE.count", row.names = 1)
colnames(OR_count) <- c("WT1","WT2","OR1","OR2")
condition <- factor(rep(c("control","OR"), c(2, 2)))
colData <- data.frame(row.names = colnames(OR_count), condition)

OR_count_dds <- DESeqDataSetFromMatrix(countData = OR_count,
                                 colData = colData ,
                                 design = ~ condition)

keep <- rowSums(counts(OR_count_dds)) >= 10
OR_count_dds <- OR_count_dds[keep, ]
OR_count_dds <- DESeq(OR_count_dds)
scaling_factors <- sizeFactors(OR_count_dds)

normalized <- as.data.frame(counts(OR_count_dds, normalized = TRUE))
write.table(normalized, file ="H1975_OR_diffGeneTE_norm.txt", sep ="\t", quote = FALSE)
OR_count_res <- results(OR_count_dds)
write.table(OR_count_res, file ="H1975_OR_diffGeneTE.txt", sep ="\t", col.names = TRUE, quote = FALSE)

L1_count <- read.table("H1975_OR_L16k.count", row.names = 1)
colnames(L1_count) <- c("WT1","WT2","OR1","OR2")
condition <- factor(rep(c("control","OR"), c(2, 2)))
colData <- data.frame(row.names = colnames(L1_count), condition)

L1_count_dds <- DESeqDataSetFromMatrix(countData = L1_count,
                                 colData = colData ,
                                 design = ~ condition)
sizeFactors(L1_count_dds) <- scaling_factors

L1_count_dds <- DESeq(L1_count_dds)
normalized <- as.data.frame(counts(L1_count_dds, normalized = TRUE))
write.table(normalized, file ="H1975_OR_diffL1_genescale_norm.txt", sep ="\t", quote = FALSE)

L1_count_res <- results(L1_count_dds)
write.table(L1_count_res, file ="H1975_OR_diffL1_genescale.txt", sep ="\t", col.names = TRUE, quote = FALSE)

# -------------------------------------------------------------------------------
# Figure 4: Volcano plot of differential LINE-1 expression
# -------------------------------------------------------------------------------

H1975_OR_L1 <- read.table("H1975_OR_diffL1_genescale.txt", header = TRUE)
H1975_OR_L1 <- na.omit(H1975_OR_L1)
H1975_OR_L1$class <- ifelse(
  abs(H1975_OR_L1$log2FoldChange) >= 0.5 & H1975_OR_L1$padj <= 0.05 ,
  ifelse(H1975_OR_L1$log2FoldChange >= 0.5 & H1975_OR_L1$padj <= 0.05,'Up','Down'),'None')
table(H1975_OR_L1$class)

ggplot(H1975_OR_L1, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = class), alpha = 0.8, size = 3) +
  scale_color_manual(values = c("royalblue4","grey","firebrick2"), labels = c("Down (210)","None","Up (90)")) +
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
  labs(x ="Log2FC(H1975 OR/H1975)", y ="-Log10(padj)") +
  theme(plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle("LINE-1s")

ggsave("H1975_OR_L1_volcano.pdf", width = 13, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Differential epigenomic signal at LINE-1 5' UTRs
# -------------------------------------------------------------------------------
read_count <- read.table("H1975_OR_L16k_tss1kb_H3K27ac.count", row.names = 1, header = FALSE)
colnames(read_count) <- c("WT1","WT2","OR1","OR2")
condition <- factor(rep(c("control","OR"), c(2, 2)))
colData <- data.frame(row.names = colnames(read_count), condition)
sample_factors <- read.table("H3K27ac_factor.txt", header = TRUE)
scaling_factors <- 1/sample_factors
read_count_dds <- DESeqDataSetFromMatrix(countData = read_count,
                                 colData = colData ,
                                 design = ~ condition)
sizeFactors(read_count_dds) <- scaling_factors
sizeFactors(read_count_dds)
read_count_dds <- DESeq(read_count_dds)

# Save normalized counts.
normalized <- as.data.frame(counts(read_count_dds, normalized = TRUE))
write.table(normalized, file ="H1975_OR_H3K27ac_L16k_tss1kb_norm.txt", sep ="\t", quote = FALSE)

# Save differential counts.
OR_diffpeak_L16k <- results(read_count_dds)
write.table(OR_diffpeak_L16k, file ="H1975_OR_H3K27ac_diffL1_1kb.txt", sep ="\t", col.names = TRUE, quote = FALSE)

# -------------------------------------------------------------------------------
# Association between LINE-1 expression and epigenomic changes in resistant cells
# -------------------------------------------------------------------------------

H1975_OR_L16k <- read.table("H1975_OR_all_diffL1_genescale.txt", header = TRUE)
H1975_OR_L16k <- na.omit(H1975_OR_L16k)
H1975_OR_L16k_H3K27ac <- read.table("H1975_OR_H3K27ac_diffL1_1kb.txt", header = TRUE)
H1975_OR_L16k_H3K27ac_merge <- H1975_OR_L16k %>%
  inner_join(H1975_OR_L16k_H3K27ac, by ="LINE1", suffix = c("_exp","_27ac"))
H1975_OR_L16k_H3K27ac_merge_logFC <- H1975_OR_L16k_H3K27ac_merge[, c(9, 3)]
cor(H1975_OR_L16k_H3K27ac_merge_logFC)
cor_result <- cor.test(
  H1975_OR_L16k_H3K27ac_merge_logFC[, 1],
  H1975_OR_L16k_H3K27ac_merge_logFC[, 2],
  method ="pearson"
)

pdf("H1975_OR_L16k_H3K27ac Scatterplot.pdf", 8, 8)
par(mar = c(5, 5, 3, 3))
palette <- colorRampPalette(c("blue","yellow","red"))
smoothScatter(H1975_OR_L16k_H3K27ac_merge_logFC,
              colramp = palette,
              main ="LINE-1 loci",
              xlab ="H3K27ac Changes (log2FC)",
              ylab ="Expression Changes (log2FC)",
              cex.main = 1.5,
              cex.lab = 1.5,
              cex.axis = 1.5,
              font.main = 2,
              font.lab = 2,
              xaxs ="i",
              yaxs ="i")

abline(h = 0, v = 0, col ="gray40", lty = 2, lwd = 1.5)

x_val <- H1975_OR_L16k_H3K27ac_merge_logFC[, 1]
y_val <- H1975_OR_L16k_H3K27ac_merge_logFC[, 2]
total_pts <- length(x_val)

if(total_pts > 0) {
  q1_per <- round(sum(x_val > 0 & y_val > 0) / total_pts * 100, 1)
  q2_per <- round(sum(x_val < 0 & y_val > 0) / total_pts * 100, 1)
  q3_per <- round(sum(x_val < 0 & y_val < 0) / total_pts * 100, 1)
  q4_per <- round(sum(x_val > 0 & y_val < 0) / total_pts * 100, 1)

  usr <- par("usr")
  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]

  text(usr[2] - x_range * 0.08, usr[4] - y_range * 0.06, paste0(q1_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[1] + x_range * 0.08, usr[4] - y_range * 0.06, paste0(q2_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[1] + x_range * 0.08, usr[3] + y_range * 0.06, paste0(q3_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[2] - x_range * 0.08, usr[3] + y_range * 0.06, paste0(q4_per,"%"), cex = 1.4, font = 2, col ="black")
}

dev.off()

# -------------------------------------------------------------------------------
# Association between LINE-1 expression and DNA methylation changes in resistant cells
# -------------------------------------------------------------------------------

ctrl1 <- read.table("H1975-1_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
ctrl2 <- read.table("H1975-2_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)

treat1 <- read.table("H1975_OR-1_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
treat2 <- read.table("H1975_OR-2_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)

mCG_ctrl <- bind_rows(ctrl1, ctrl2) %>%
  group_by(LINE1) %>%
  summarize(mean_ctrl = mean(as.numeric(mCG), na.rm = TRUE), .groups ="drop")

mCG_treat <- bind_rows(treat1, treat2) %>%
  group_by(LINE1) %>%
  summarize(mean_treat = mean(as.numeric(mCG), na.rm = TRUE), .groups ="drop")

methylation_delta <- inner_join(mCG_ctrl, mCG_treat, by ="LINE1") %>%
  mutate(diffMeth = mean_treat - mean_ctrl) %>%
  select(LINE1, diffMeth)
methylation_delta <- subset(methylation_delta, diffMeth > -0.4)

# Read differential LINE-1 expression results.
H1975_OR_L16k <- read.table("H1975_OR_diffL1_genescale.txt", header = TRUE)
H1975_OR_L16k <- na.omit(H1975_OR_L16k)
exp_df <- H1975_OR_L16k %>% select(LINE1, log2FoldChange)

plot_data <- inner_join(methylation_delta, exp_df, by ="LINE1") %>%
  filter(!is.na(diffMeth) & !is.na(log2FoldChange))

x_val <- plot_data$diffMeth
y_val <- plot_data$log2FoldChange
total_pts <- length(x_val)
cor_result <- cor.test(
  x_val,
  y_val,
  method ="pearson"
)

pdf("H1975_OR_L16k_mCG Scatterplot.pdf", 8, 8)
par(mar = c(5, 5, 3, 3))
palette <- colorRampPalette(c("blue","yellow","red"))
smoothScatter(x = x_val,
              y = y_val,
              colramp = palette,
              main ="Full-length LINE-1",
              xlab ="Methylation Changes (Δβ)",
              ylab ="Expression Changes (log2FC)",
              cex.main = 1.5,
              cex.lab = 1.5,
              cex.axis = 1.5,
              font.main = 2,
              font.lab = 2,
              xaxs ="i",
              yaxs ="i")

abline(h = 0, v = 0, col ="gray40", lty = 2, lwd = 1.5)

if(total_pts > 0) {
  q1_per <- round(sum(x_val > 0 & y_val > 0) / total_pts * 100, 1)
  q2_per <- round(sum(x_val < 0 & y_val > 0) / total_pts * 100, 1)
  q3_per <- round(sum(x_val < 0 & y_val < 0) / total_pts * 100, 1)
  q4_per <- round(sum(x_val > 0 & y_val < 0) / total_pts * 100, 1)

  usr <- par("usr")
  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]

  text(usr[2] - x_range * 0.08, usr[4] - y_range * 0.06, paste0(q1_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[1] + x_range * 0.08, usr[4] - y_range * 0.06, paste0(q2_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[1] + x_range * 0.08, usr[3] + y_range * 0.06, paste0(q3_per,"%"), cex = 1.4, font = 2, col ="black")
  text(usr[2] - x_range * 0.08, usr[3] + y_range * 0.06, paste0(q4_per,"%"), cex = 1.4, font = 2, col ="black")
}

dev.off()

# -------------------------------------------------------------------------------
# Integrated multi-omic state transitions at LINE-1 loci
# -------------------------------------------------------------------------------

H1975_OR_L16k <- read.table("H1975_OR_diffL1_genescale.txt", header = TRUE)
H1975_OR_L16k_H3K27ac <- read.table("H1975_OR_H3K27ac_diffL1_1kb.txt", header = TRUE)
H1975_OR_L16k_H3K9me3 <- read.table("H1975_OR_H3K9me3_diffL1_1kb.txt", header = TRUE)
H1975_OR_L16k_ssDNA <- read.table("H1975_OR_KAS-seq_diffL1_1kb.txt", header = TRUE)

df_rna <- H1975_OR_L16k %>% select(LINE1, log2FC_RNA = log2FoldChange)
df_k9 <- H1975_OR_L16k_H3K9me3 %>% select(LINE1, log2FC_K9me3 = log2FoldChange)
df_k27 <- H1975_OR_L16k_H3K27ac %>% select(LINE1, log2FC_K27ac = log2FoldChange)
df_ssdna <- H1975_OR_L16k_ssDNA %>% select(LINE1, log2FC_ssDNA = log2FoldChange)

L1_multi_omics <- df_rna %>%
  inner_join(df_k9, by ="LINE1") %>%
  inner_join(df_k27, by ="LINE1") %>%
  inner_join(df_ssdna, by ="LINE1") %>%
  na.omit()

alluvial_binary_df <- L1_multi_omics %>%
  mutate(
    K27ac_Status = ifelse(log2FC_K27ac > 0,"K27ac_Up","K27ac_Down"),
    K9me3_Status = ifelse(log2FC_K9me3 > 0,"K9me3_Up","K9me3_Down"),
    ssDNA_Status = ifelse(log2FC_ssDNA > 0,"ssDNA_Up","ssDNA_Down"),
    RNA_Status = ifelse(log2FC_RNA > 0,"RNA_Up","RNA_Down")
  ) %>%
  group_by(K27ac_Status, K9me3_Status, ssDNA_Status, RNA_Status) %>%
  tally(name ="Freq")

alluvial_binary_df <- alluvial_binary_df %>%
  mutate(
    Highlight_Group = case_when(
      K27ac_Status =="K27ac_Down" & K9me3_Status =="K9me3_Up" & RNA_Status =="RNA_Down" ~"Core_Mechanism",
      RNA_Status =="RNA_Down" ~"Other_RNA_Down",
      TRUE ~"Other_Stream"
    )
  )

ggplot(alluvial_binary_df, aes(y = Freq,
                               axis1 = K27ac_Status,
                               axis2 = K9me3_Status,
                               axis3 = ssDNA_Status,
                               axis4 = RNA_Status)) +
  geom_alluvium(aes(fill = Highlight_Group), width = 1/8, alpha = 0.7, knot.prop = 0.4) +
  geom_stratum(width = 1/8, fill ="#F8F9FA", color ="gray40", linewidth = 0.5) +
  geom_text(stat ="stratum", aes(label = after_stat(stratum)), size = 3.8, fontface ="bold") +
  scale_x_discrete(limits = c("H3K27ac","H3K9me3","ssDNA","RNA-seq"), expand = c(.06, .06)) +
  scale_fill_manual(values = c(
"Core_Mechanism" ="#08519C",
"Other_RNA_Down" ="#6BAED6",
"Other_Stream" ="#D9D9D9"
  )) +
  theme_minimal() +
  labs(
    title ="LINE-1 Epigenetic Modulation",
    y ="Number of Full-length LINE-1 Loci",
    fill ="Functional Stream"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face ="bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color ="gray40", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color ="black", size = 10),
    axis.text.x = element_text(face ="bold", size = 12, color ="black"),
    legend.position ="right"
  )

ggsave("H1975 OR L1 alluvial plot.pdf", width = 10, height = 7, dpi = 300)

# -------------------------------------------------------------------------------
# Neighboring gene expression around downregulated LINE-1 loci
# -------------------------------------------------------------------------------

downL1_geneFC <- read.csv("H1975 OR downL1 neighborgene logFC.csv", header = TRUE)
downL1_geneFC_rep <- melt(downL1_geneFC)
downL1_geneFC_rep <- na.omit(downL1_geneFC_rep)
downgroup_counts <- downL1_geneFC_rep %>%
  group_by(variable) %>%
  summarize(count = n(), max_value = 4)

ggplot(downL1_geneFC_rep, aes(x = variable, y = value)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, fill ="ghostwhite") +
  geom_hline(yintercept = 0, linetype = 1, color ='red') +
  scale_y_continuous(limits = c(-3, 2), breaks = seq(-2, 2, 2)) +
  theme_bw() +theme(panel.grid = element_blank(), legend.position ="none") +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5), legend.position ="none") +
  theme(axis.title.y = element_text(vjust = 2, size = 28)) +
  theme(axis.text.y = element_text(vjust = 1, size = 25)) +
  theme(axis.title.x = element_text(vjust = 2, size = 28)) +
  theme(axis.text.x = element_text(vjust = 0.6, size = 25, angle = 25)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25)) +
  labs(title ="Neighboring gene expression around downregulated LINE-1 loci",
       x ="Distance to downregulated LINE-1 (kb)", y ="Gene expression change")
ggsave("H1975 OR downregulate LINE1 gene expression changes.pdf", width = 12, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Subfamily distribution of downregulated LINE-1 loci
# -------------------------------------------------------------------------------

H1975_OR_L16k <- read.table("H1975_OR_diffL1_genescale.txt", header = TRUE)
H1975_OR_L16k <- na.omit(H1975_OR_L16k)
H1975_OR_L16k$class <- ifelse(
  abs(H1975_OR_L16k$log2FoldChange) >= 0.5 & H1975_OR_L16k$padj <= 0.05 ,
  ifelse(H1975_OR_L16k$log2FoldChange >= 0.5 & H1975_OR_L16k$padj <= 0.05,'Up','Down'),'None')
table(H1975_OR_L16k$class)
H1975_OR_L16k_down <- subset(H1975_OR_L16k, class %in%"Down")

PC9_OR_L16k <- read.table("PC9_OR_diffL1_genescale.txt", header = TRUE)
PC9_OR_L16k <- na.omit(PC9_OR_L16k)
PC9_OR_L16k$class <- ifelse(
  abs(PC9_OR_L16k$log2FoldChange) >= 0.5 & PC9_OR_L16k$padj <= 0.05 ,
  ifelse(PC9_OR_L16k$log2FoldChange >= 0.5 & PC9_OR_L16k$padj <= 0.05,'Up','Down'),'None')
table(PC9_OR_L16k$class)
PC9_OR_L16k_down <- subset(PC9_OR_L16k, class %in%"Down")

HCC827_OR_L16k <- read.table("HCC827_OR_diffL1_genescale.txt", header = TRUE)
HCC827_OR_L16k <- na.omit(HCC827_OR_L16k)
HCC827_OR_L16k$class <- ifelse(
  abs(HCC827_OR_L16k$log2FoldChange) >= 0.5 & HCC827_OR_L16k$padj <= 0.05 ,
  ifelse(HCC827_OR_L16k$log2FoldChange >= 0.5 & HCC827_OR_L16k$padj <= 0.05,'Up','Down'),'None')
table(HCC827_OR_L16k$class)
HCC827_OR_L16k_down <- subset(HCC827_OR_L16k, class %in%"Down")

target_subfamilies <- c("L1HS", paste0("L1PA", 2:7))
l1_raw <- read.table("L1.gt6k.bed", header = FALSE, sep ="\t", stringsAsFactors = FALSE)
l1_background <- data.frame(
  LINE1 = l1_raw$V4,
  Subfamily = l1_raw$V5,
  stringsAsFactors = FALSE
) %>%
  filter(Subfamily %in% target_subfamilies)

bg_counts <- l1_background %>%
  count(Subfamily, name ="Total_Background")

down_data_list <- list(
"H1975_OR" = H1975_OR_L16k_down,
"PC9_OR" = PC9_OR_L16k_down,
"HCC827_OR" = HCC827_OR_L16k_down
)

calculate_down_percentage <- function(cell_name, down_df) {
  down_df %>%
    inner_join(l1_background, by ="LINE1") %>%
    count(Subfamily, name ="Down_Count") %>%
    complete(Subfamily = target_subfamilies, fill = list(Down_Count = 0)) %>%
    left_join(bg_counts, by ="Subfamily") %>%
    mutate(
      CellLine = cell_name,
      Percentage = (Down_Count / Total_Background) * 100
    )
}

final_stats <- map2_dfr(names(down_data_list), down_data_list, calculate_down_percentage) %>%
  mutate(Subfamily = factor(Subfamily, levels = target_subfamilies))
write.csv(final_stats,"LINE1_down_regulated_percentage_summary.csv", row.names = FALSE)

ggplot(final_stats, aes(x = Subfamily, y = Percentage, fill = CellLine)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, color ="white", size = 0.3) +
  scale_fill_manual(values = c("H1975_OR" ="#e74c3c","PC9_OR" ="#3498db","HCC827_OR" ="#2ecc71")) +
  labs(
    x ="LINE-1 Subfamilies",
    y ="Percentage (%)",
    title ="Proportion of Downregulated LINE-1 Loci",
    fill ="Cell Line"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic() +
  theme(
    plot.title = element_text(face ="bold", size = 13, hjust = 0.5, color ="black"),
    axis.title = element_text(face ="bold", size = 11, color ="black"),
    axis.text = element_text(face ="bold", size = 10, color ="black"),
    legend.title = element_text(face ="bold", size = 10),
    legend.position ="top"
  )

ggsave("OR downregulate LINE1 percentage.pdf", width = 8, height = 6, dpi = 300)

listInput <- list(
  H1975_OR_downL1 = H1975_OR_L16k_down$LINE1,
  PC9_OR_downL1 = PC9_OR_L16k_down$LINE1,
  HCC827_OR_downL1 = HCC827_OR_L16k_down$LINE1)

pdf("OR downregulate LINE1 upsetR.pdf", 12, 8)
upset(fromList(listInput), order.by ="freq",
      point.size = 2.5,
      line.size = 1,
      mainbar.y.label ="Intersect LINE-1 loci",
      sets.x.label ="Number of Genes",
      text.scale = c(3, 4, 3, 2, 3, 4),
      queries = list(
        list(query = intersects,
             params = list("H1975_OR_downL1","PC9_OR_downL1","HCC827_OR_downL1"),
             color ="red3",
             active = TRUE)))
dev.off()

# -------------------------------------------------------------------------------
# Distance-dependent correlation between LINE-1 and neighboring genes
# -------------------------------------------------------------------------------

H1975_OR_L16k <- read.table("H1975_OR_diffL1_genescale.txt", header = TRUE)
H1975_OR_DEG <- read.csv("H1975 OR RNA-seq DEG.csv", header = TRUE)

distances <- c("2kb","100kb","200kb","300kb","400kb","500kb")
cor_results_summary <- data.frame(Distance = character(),
                                  Pearson_R = numeric(),
                                  P_value = numeric())

for (dist in distances) {

  input_file <- paste0("Gene_L16k_tss", dist,".txt")
  L1_gene_data <- read.table(input_file, header = FALSE)
  colnames(L1_gene_data) <- c("LINE1","symbol","ENSEMBL")
  L1_matched <- H1975_OR_L16k[match(L1_gene_data$LINE1, H1975_OR_L16k$LINE1), ]
  gene_matched <- H1975_OR_DEG[match(L1_gene_data$ENSEMBL, H1975_OR_DEG$ENSEMBL), ]
  merge_data <- cbind(L1_matched, gene_matched)
  merge_data <- na.omit(merge_data)
  csv_name <- paste0("H1975_OR_L1_tss", dist,"_gene_merge.csv")
  write.csv(merge_data, file = csv_name, quote = FALSE, row.names = FALSE)

  logFC_data <- na.omit(merge_data[, c(3, 12)])

  pdf_name <- paste0("H1975_OR_L1_tss", dist,"_gene_cor.pdf")
  pdf(pdf_name, 8, 8)

  par(mar = c(5, 5, 3, 3))
  palette_func <- colorRampPalette(c("blue","yellow","red"))

  smoothScatter(logFC_data,
                colramp = palette_func,
                main = paste0("LINE-1 TSS ±", dist),
                xlab ="LINE-1 log2FC",
                ylab ="Gene log2FC",
                cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,
                font.main = 2, font.lab = 2,
                xaxs ="i", yaxs ="i")

  abline(h = 0, v = 0, col ="gray40", lty = 2, lwd = 1.5)

  if(nrow(logFC_data) > 2) {
    x_val <- logFC_data[, 1]
    y_val <- logFC_data[, 2]

    cor_test <- cor.test(x_val, y_val, method ="pearson", use ="complete.obs")
    rho <- round(cor_test$estimate, 3)
    p_val <- cor_test$p.value

    cor_results_summary <- rbind(cor_results_summary,
                                 data.frame(Distance = dist, Pearson_R = rho, P_value = p_val))

    p_label <- ifelse(p_val < 0.0001,"< 0.0001", round(p_val, 4))
    legend("topleft",
           legend = c(paste0("R = ", rho), paste0("p ", p_label)),
           bty ="n", cex = 1.5, text.col ="darkred")

    total_pts <- length(x_val)
    q1_per <- round(sum(x_val > 0 & y_val > 0) / total_pts * 100, 1)
    q2_per <- round(sum(x_val < 0 & y_val > 0) / total_pts * 100, 1)
    q3_per <- round(sum(x_val < 0 & y_val < 0) / total_pts * 100, 1)
    q4_per <- round(sum(x_val > 0 & y_val < 0) / total_pts * 100, 1)

    usr <- par("usr")
    x_range <- usr[2] - usr[1]
    y_range <- usr[4] - usr[3]

    text(usr[2] - x_range * 0.1, usr[4] - y_range * 0.08, paste0(q1_per,"%"), cex = 1.4, font = 2, col ="black")
    text(usr[1] + x_range * 0.1, usr[4] - y_range * 0.25, paste0(q2_per,"%"), cex = 1.4, font = 2, col ="black")
    text(usr[1] + x_range * 0.1, usr[3] + y_range * 0.08, paste0(q3_per,"%"), cex = 1.4, font = 2, col ="black")
    text(usr[2] - x_range * 0.1, usr[3] + y_range * 0.08, paste0(q4_per,"%"), cex = 1.4, font = 2, col ="black")
  }

  dev.off()
}

print(cor_results_summary)
write.csv(cor_results_summary,"H1975_OR Summary_L1_Gene_Correlations.csv", row.names = FALSE)

df_h1975 <- data.frame(cor_results_summary)
df_h1975$CellLine <-"H1975"

if (!exists("df_pc9") || !exists("df_hcc827")) {
  stop("Generate df_pc9 and df_hcc827 using the same correlation workflow before the combined plot.")
}

plot_data <- rbind(df_h1975, df_hcc827, df_pc9)
plot_data$Distance <- factor(plot_data$Distance,
                             levels = c("2kb","100kb","200kb","300kb","400kb","500kb"))

get_stars <- function(p) {
  ifelse(p < 0.0001,"****",
         ifelse(p < 0.001,"***",
                ifelse(p < 0.01,"**",
                       ifelse(p < 0.05,"*","ns"))))
}
plot_data$Significance <- get_stars(plot_data$P_value)

ggplot(plot_data, aes(x = Distance, y = Pearson_R, group = CellLine, color = CellLine)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  geom_text(aes(label = Significance), vjust = -0.6, show.legend = FALSE, fontface ="bold", size = 4.5) +
  scale_y_continuous(limits = c(0, max(plot_data$Pearson_R) * 1.2)) +
  scale_color_manual(values = c(
"H1975" ="#e74c3c",
"HCC827" ="#2ecc71",
"PC9" ="#3498db"
  )) +
  theme_classic() +
  labs(
    title ="Correlation",
    x ="Genomic Distance to LINE-1 TSS",
    y ="Pearson Correlation Coefficient (R)",
    color ="Cell Line"
  ) +
  theme(
    plot.title = element_text(size = 16, face ="bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color ="gray30", hjust = 0.5),
    axis.title = element_text(size = 13, face ="bold"),
    axis.text = element_text(size = 12, color ="black"),
    legend.title = element_text(size = 12, face ="bold"),
    legend.text = element_text(size = 11),
    legend.position ="top"
  )

ggsave("LINE1_Gene_Cor_Distance_Decay.pdf", width = 7, height = 5.5)
