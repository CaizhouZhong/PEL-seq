# ================================================================================
# Figure 1 and Supplementary Figure S1
# PEL-seq performance and locus-level LINE-1 expression
# ================================================================================

library(tidyverse)
library(reshape2)

# -------------------------------------------------------------------------------
# LINE-1 read composition in public sequencing datasets
# -------------------------------------------------------------------------------

Public_count <- read.csv("Public LINE1 count.csv", header = TRUE)
Public_long <- melt(Public_count, id.vars ="sample")
Public_percent <- Public_long %>%
  group_by(sample) %>%
  mutate(percent = value / sum(value) * 100)
Public_percent$sample <- factor(Public_percent$sample, levels = c("RNA-seq_r1","RNA-seq_r2","RNA-seq_r3","ONT-seq_r1","ONT-seq_r2","ONT-seq_r3"))

ggplot(Public_percent, aes(x = sample, y = percent, fill = variable)) +
  geom_bar(stat ="identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", percent)),
            position = position_stack(vjust = 0.5), size = 7) +
  labs(x =" ", y ="Percentage", fill ="Category") +
  scale_fill_manual(values = c('aquamarine2','skyblue','snow2')) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  theme(axis.title.y = element_text(vjust = 2, size = 28)) +
  theme(axis.text.y = element_text(vjust = 1, size = 28)) +
  theme(axis.title.x = element_text(vjust = 2, size = 28)) +
  theme(axis.text.x = element_text(vjust = 0.6, size = 25)) +
  theme(legend.text = element_text(size = 25))
ggsave("Public LINE-1 percentage.pdf", width = 15, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# LINE-1 composition of PEL-seq enrichment
# -------------------------------------------------------------------------------

PEL_seq_count <- read.csv("PEL-seq L1 celltype.csv", header = TRUE)
PEL_seq_long <- melt(PEL_seq_count, id.vars ="sample")
PEL_seq_percent <- PEL_seq_long %>%
  group_by(sample) %>%
  mutate(percent = value / sum(value) * 100)
PEL_seq_percent$sample <- factor(PEL_seq_percent$sample, levels = c("H1975_r1","H1975_r2","PC9_r1","PC9_r2","HCC827_r1","HCC827_r2"))

ggplot(PEL_seq_percent, aes(x = sample, y = percent, fill = variable)) +
  geom_bar(stat ="identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", percent)),
            position = position_stack(vjust = 0.5), size = 7) +
  labs(x =" ", y ="Percentage", fill ="Category") +
  scale_fill_manual(values = c('aquamarine2','skyblue','snow2')) +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  theme(axis.title.y = element_text(vjust = 2, size = 28)) +
  theme(axis.text.y = element_text(vjust = 1, size = 28)) +
  theme(axis.title.x = element_text(vjust = 2, size = 28)) +
  theme(axis.text.x = element_text(vjust = 0.6, size = 25)) +
  theme(legend.text = element_text(size = 25))
ggsave("PEL-seq LINE-1 percentage.pdf", width = 15, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Detection of full-length LINE-1 loci by subfamily
# -------------------------------------------------------------------------------

target_subfamilies <- c("L1PA7","L1PA6","L1PA5","L1PA4","L1PA3","L1PA2","L1HS")
tpm_data <- read.table("H1975_L16k.TPM.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>%
  mutate(LINE1 = sub("_[^_]+$","", ID)) %>%
  filter(LINE1 %in% target_subfamilies)

subf_stats <- tpm_data %>%
  group_by(LINE1) %>%
  summarise(
    Total_Loci = n(),
    Det_Rep1 = sum(TPM1 >= 3),
    Det_Rep2 = sum(TPM2 >= 3),
    .groups ="drop"
  ) %>%
  mutate(
    Mean_Det = (Det_Rep1 + Det_Rep2) / 2,
    Min_Det = pmin(Det_Rep1, Det_Rep2),
    Max_Det = pmax(Det_Rep1, Det_Rep2),
    Pct_Det = round((Mean_Det / Total_Loci) * 100),
    Pct_Label = paste0(Pct_Det,"%"),
    Text_X = Mean_Det / 2,
    LINE1 = factor(LINE1, levels = target_subfamilies)
  )

ggplot(subf_stats, aes(x = LINE1)) +
  geom_col(
    aes(y = Total_Loci, fill ="Not detected"),
    width = 0.65 ) +
  geom_col(
    aes(y = Mean_Det, fill ="Detected"),
    width = 0.65 ) +
  geom_errorbar(
    aes(ymin = Min_Det, ymax = Max_Det),
    width = 0.22,
    color ="black",
    size = 0.8 ) +
  geom_text(
    aes(y = Text_X, label = Pct_Label),
    color ="black",
    size = 5.5,
    fontface ="plain" ) +
  scale_fill_manual(
    values = c("Detected" ="mediumturquoise","Not detected" ="paleturquoise"),
    breaks = c("Detected","Not detected") ) +
  labs(
    x ="",
    y ="Count",
    fill ="",
    title ="Full-length LINE-1 detected by PEL-seq") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color ="black", fill = NA, size = 1.2),
    plot.title = element_text(hjust = 0.5, size = 22, face ="plain", margin = margin(b = 15)),
    axis.text.y = element_text(size = 20, color ="black", face ="plain"),
    axis.title.x = element_text(size = 22, color ="black", face ="plain", vjust = -0.5),
    axis.text.x = element_text(size = 18, color ="black", face ="plain"),
    legend.position = c(0.82, 0.90),
    legend.background = element_blank(),
    legend.key.size = grid::unit(0.8,"cm"),
    legend.text = element_text(size = 18, color ="black")
  )

ggsave("H1975_full_length_LINE1_detected.pdf", width = 9, height = 8, dpi = 300)

# -------------------------------------------------------------------------------
# LINE-1 expression by subfamily
# -------------------------------------------------------------------------------

L1_tpm <- read.table("H1975_L16k.TPM.txt", header = TRUE)
L1_tpm$Expression <- (L1_tpm$TPM1+L1_tpm$TPM2)/2
L1_tpm$subfamily <- sub("_dup\\d+","", L1_tpm$LINE1)
L1_tpm_subfamily <- subset(L1_tpm, subfamily %in% c("L1HS","L1PA2","L1PA3","L1PA4","L1PA5","L1PA6","L1PA7"))
L1_tpm_subfamily$subfamily <- factor(L1_tpm_subfamily$subfamily, levels = c("L1PA7","L1PA6","L1PA5","L1PA4","L1PA3","L1PA2","L1HS"))
ggplot(L1_tpm_subfamily, aes(x = subfamily, y = log2(Expression+1))) +
  stat_boxplot(geom ="errorbar", width = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill ="mediumturquoise") +
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank()) +
  theme(axis.title.y = element_text(vjust = 0.8, size = 35)) +
  theme(axis.text.y = element_text(size = 32)) +
  theme(axis.title.x = element_text(vjust = 0.8, size = 35)) +
  theme(axis.text.x = element_text(size = 32)) +
  labs(x =" ", y ="Log2(TPM+1)") +
  theme(plot.title = element_text(hjust = 0.5, size = 35)) +
  ggtitle("H1975")
ggsave("H1975 L1 family log2TPM boxplot.pdf", width = 13, height = 10, dpi = 300)

# -------------------------------------------------------------------------------
# Read-length distribution
# -------------------------------------------------------------------------------
pdf("Reads_length_density.pdf", 13, 12)
Reads_length <- read.table("Sample_Reads_length.txt", header = TRUE)
plot(density(Reads_length$length), col ="gray", lwd = 2, main ="Reads length distribution",
     xlim = c(0, 4000), cex.lab = 2, cex.axis = 1.5, cex.main = 2,
     xlab ="Reads length (bp)", ylab ="Density")
dev.off()

# -------------------------------------------------------------------------------
# Mapping-quality distribution
# -------------------------------------------------------------------------------
pdf("Reads_Mapping_distr.pdf", 13, 12)
Reads_Mapping <- read.table("Sample_Reads_Mapping.txt", header = TRUE)
hist(
  Reads_Mapping$MAPQ,
  breaks = 60,
  col ="gray",
  xlim = c(0, 60),
  cex.lab = 2,
  cex.axis = 1.5,
  cex.main = 2,
  main ="Read mapping quality",
  xlab ="MAPQ",
  ylab ="Frequency"
)
dev.off()

# -------------------------------------------------------------------------------
# Contribution of highly expressed LINE-1 loci
# -------------------------------------------------------------------------------

L1_tpm <- read.csv("H1975_L16k.TPM.txt", header = TRUE)
L1_tpm$Expression <- (L1_tpm$L1_1+L1_tpm$L1_2)/2
L1_tpm <- L1_tpm[order(L1_tpm$Expression, decreasing = TRUE), ]
L1_tpm$subfamily <- sub("_dup\\d+","", L1_tpm$LINE1)
L1_tpm_L1HS <- L1_tpm[grepl("^L1HS_", L1_tpm$LINE1), ]
L1_tpm_L1HS$per <- L1_tpm_L1HS$Expression / sum(L1_tpm_L1HS$Expression) * 100
L1_tpm_L1PA2 <- L1_tpm[grepl("^L1PA2_", L1_tpm$LINE1), ]
L1_tpm_L1PA2$per <- L1_tpm_L1PA2$Expression / sum(L1_tpm_L1PA2$Expression) * 100
L1_tpm_L1PA3 <- L1_tpm[grepl("^L1PA3_", L1_tpm$LINE1), ]
L1_tpm_L1PA3$per <- L1_tpm_L1PA3$Expression / sum(L1_tpm_L1PA3$Expression) * 100
L1_tpm_L1PA4 <- L1_tpm[grepl("^L1PA4_", L1_tpm$LINE1), ]
L1_tpm_L1PA4$per <- L1_tpm_L1PA4$Expression / sum(L1_tpm_L1PA4$Expression) * 100
L1_tpm_L1PA5 <- L1_tpm[grepl("^L1PA5_", L1_tpm$LINE1), ]
L1_tpm_L1PA5$per <- L1_tpm_L1PA5$Expression / sum(L1_tpm_L1PA5$Expression) * 100
L1_tpm_L1PA6 <- L1_tpm[grepl("^L1PA6_", L1_tpm$LINE1), ]
L1_tpm_L1PA6$per <- L1_tpm_L1PA6$Expression / sum(L1_tpm_L1PA6$Expression) * 100
L1_tpm_L1PA7 <- L1_tpm[grepl("^L1PA7_", L1_tpm$LINE1), ]
L1_tpm_L1PA7$per <- L1_tpm_L1PA7$Expression / sum(L1_tpm_L1PA7$Expression) * 100

L1_top100 <- cbind(
  L1_tpm_L1HS[1:100, ],
  L1_tpm_L1PA2[1:100, ],
  L1_tpm_L1PA3[1:100, ],
  L1_tpm_L1PA4[1:100, ],
  L1_tpm_L1PA5[1:100, ],
  L1_tpm_L1PA6[1:100, ],
  L1_tpm_L1PA7[1:100, ]
)
L1_top100_sub <- L1_top100[, c(4, 8, 12, 16, 20, 24, 28)]
L1_top100_sub_fal <- rbind(L1_top100_sub, 100-colSums(L1_top100_sub))

colnames(L1_top100_sub_fal) <- c("L1HS","L1PA2","L1PA3","L1PA4","L1PA5","L1PA6","L1PA7")
L1_top100_sub_fal_long <- melt(L1_top100_sub_fal)
L1_top100_sub_fal_long$variable <- factor(L1_top100_sub_fal_long$variable, levels = c("L1PA7","L1PA6","L1PA5","L1PA4","L1PA3","L1PA2","L1HS"))
labels <- paste0("C", 1:101)
L1_top100_sub_fal_long$rank <- rep(labels, times = 7)
L1_top100_sub_fal_long$rank <- factor(L1_top100_sub_fal_long$rank, levels = labels)
base_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
"#8c564b","#e377c2","snow","#bcbd22","#17becf")
generate_gradient <- function(start_color, end_color, n) {
  colorRampPalette(c(start_color, end_color))(n)
}
gradient_colors <- unlist(lapply(seq_along(base_colors)[-length(base_colors)], function(i) {
  generate_gradient(base_colors[i], base_colors[i + 1], 11)
}))
colors <- c(gradient_colors,"grey","grey")
final_colors <- rep(colors, times = 7)

ggplot(L1_top100_sub_fal_long, aes(x = variable, y = value, fill = rank)) +
  geom_bar(stat ="identity", width = 0.65, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = colors) +
  labs(x =" ", y ="Percentage") +
  ggtitle("Cumulative Distribution") +
  theme_classic() +
  theme(
    axis.title.y = element_text(vjust = 0.8, size = 35),
    axis.text.y = element_text(size = 32),
    axis.title.x = element_text(vjust = 0.8, size = 35),
    axis.text.x = element_text(size = 32),
    plot.title = element_text(hjust = 0.5, size = 32),
    legend.position ='none') +
  coord_cartesian(expand = FALSE) +
  coord_flip()
ggsave("Cumulative Distribution of LINE1 Family Expression H1975.pdf", width = 16, height = 12, dpi = 300)
