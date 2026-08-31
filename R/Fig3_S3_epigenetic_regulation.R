# ================================================================================
# Figure 3 and Supplementary Figure S3
# Epigenetic features of full-length LINE-1 loci
# ================================================================================

library(tidyverse)
library(ggpubr)

# -------------------------------------------------------------------------------
# PEL-seq expression across LINE-1 subfamilies
# -------------------------------------------------------------------------------

target_subfamilies <- c("L1HS", paste0("L1PA", 2:7))

tpm_data <- read.table("H1975_L16k.TPM.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>%
  mutate(Signal = (TPM1 + TPM2) / 2) %>%
  mutate(Log2_Signal = log2(Signal + 1)) %>%
  select(LINE1 = LINE1, Log2_Signal)

top100_log <- read.table("mergeL1_top100.sort.bed", header = FALSE, sep ="\t", stringsAsFactors = FALSE) %>%
  select(LINE1 = V4, Subfamily = V5) %>%
  filter(Subfamily %in% target_subfamilies) %>%
  inner_join(tpm_data, by ="LINE1") %>%
  mutate(Subfamily = factor(Subfamily, levels = target_subfamilies))

all_l1_log <- read.table("L1.gt6k.bed", header = FALSE, sep ="\t", stringsAsFactors = FALSE) %>%
  select(LINE1 = V4, Subfamily = V5) %>%
  filter(Subfamily %in% target_subfamilies) %>%
  inner_join(tpm_data, by ="LINE1") %>%
  mutate(Subfamily = factor(Subfamily, levels = target_subfamilies))

ggplot(top100_log, aes(x = Subfamily, y = Log2_Signal, fill = Subfamily)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, color ="#2c3e50", size = 0.6) +
  scale_fill_viridis_d(option ="plasma", end = 0.85) +
  labs(
    x ="LINE-1 Subfamilies",
    y ="PEL-seq(log2(TPM + 1))",
    title ="H1975 (Top 100 loci)"
  ) +
  theme_classic() +
  theme(
    legend.position ="none",
    plot.title = element_text(face ="bold", size = 13, hjust = 0.5, color ="black"),
    axis.title = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(face ="bold", size = 11, color ="black"),
    axis.text.y = element_text(size = 11, color ="black")
  )

ggsave("H1975_L1_Top100_Log2_boxplot.pdf", width = 7.5, height = 5, dpi = 300)

ggplot(all_l1_log, aes(x = Subfamily, y = Log2_Signal, fill = Subfamily)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, color ="#2c3e50", size = 0.6) +
  scale_fill_viridis_d(option ="plasma", end = 0.85) +
  labs(
    x ="LINE-1 Subfamilies",
    y ="PEL-seq (log2(TPM + 1))",
    title ="H1975 (All loci)"
  ) +
  theme_classic() +
  theme(
    legend.position ="none",
    plot.title = element_text(face ="bold", size = 13, hjust = 0.5, color ="black"),
    axis.title = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(face ="bold", size = 11, color ="black"),
    axis.text.y = element_text(size = 11, color ="black")
  )

ggsave("H1975_L1_All_Background_Log2_boxplot.pdf", width = 7.5, height = 5, dpi = 300)

# -------------------------------------------------------------------------------
# Epigenomic signal across LINE-1 subfamilies
# -------------------------------------------------------------------------------

matrix_file <-"H1975_H3K27ac_L16k_5UTR.gz"
matrix_data <- read.table(
  gzfile(matrix_file),
  header = FALSE,
  sep ="\t",
  stringsAsFactors = FALSE,
  comment.char ="@"
)

matrix_data$Mean_Signal <- rowMeans(matrix_data[, 7:ncol(matrix_data)], na.rm = TRUE)

target_subfamilies <- c("L1HS", paste0("L1PA", 2:7))
processed_data <- matrix_data %>%
  select(LINE1 = V4, Signal = Mean_Signal) %>%
  mutate(
    Subfamily = str_extract(LINE1,"^L1HS|^L1PA[2-7]")
  ) %>%
  filter(!is.na(Subfamily)) %>%
  mutate(
    Subfamily = factor(Subfamily, levels = target_subfamilies)
  )

plot_data_clean <- processed_data

print(table(plot_data_clean$Subfamily))

ggplot(plot_data_clean, aes(x = Subfamily, y = Signal, fill = Subfamily)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, color ="#2c3e50", size = 0.6) +
  scale_fill_viridis_d(option ="plasma", end = 0.85) +
  labs(
    x ="LINE-1 Subfamilies",
    y ="H3K27ac Signal (5' UTR)",
    title ="H1975 (All loci)"
  ) +
  theme_classic() +
  theme(
    legend.position ="none",
    plot.title = element_text(face ="bold", size = 14, hjust = 0.5, color ="black"),
    axis.title = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(face ="bold", size = 11, color ="black"),
    axis.text.y = element_text(size = 11, color ="black"),
    axis.line = element_line(size = 0.6, color ="#2c3e50")
  )

ggsave("L1_subfamilies_H3K27ac_5UTR_boxplot.pdf", width = 7.5, height = 5, dpi = 300)

matrix_file <-"H1975_H3K27ac_topL16k_5UTR.gz"

matrix_data <- read.table(
  gzfile(matrix_file),
  header = FALSE,
  sep ="\t",
  stringsAsFactors = FALSE,
  comment.char ="@"
)

matrix_data$Mean_Signal <- rowMeans(matrix_data[, 7:ncol(matrix_data)], na.rm = TRUE)
target_subfamilies <- c("L1HS", paste0("L1PA", 2:7))

processed_data <- matrix_data %>%
  select(LINE1 = V4, Signal = Mean_Signal) %>%
  mutate(
    Subfamily = str_extract(LINE1,"^L1HS|^L1PA[2-7]")
  ) %>%
  filter(!is.na(Subfamily)) %>%
  mutate(
    Subfamily = factor(Subfamily, levels = target_subfamilies)
  )

plot_data_clean <- processed_data
print(table(plot_data_clean$Subfamily))

ggplot(plot_data_clean, aes(x = Subfamily, y = Signal, fill = Subfamily)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, color ="#2c3e50", size = 0.6) +
  scale_fill_viridis_d(option ="plasma", end = 0.85) +
  labs(
    x ="LINE-1 Subfamilies",
    y ="H3K27ac Signal (5' UTR)",
    title ="H1975 (Top 100 loci)"
  ) +
  theme_classic() +
  theme(
    legend.position ="none",
    plot.title = element_text(face ="bold", size = 14, hjust = 0.5, color ="black"),
    axis.title = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(face ="bold", size = 11, color ="black"),
    axis.text.y = element_text(size = 11, color ="black"),
    axis.line = element_line(size = 0.6, color ="#2c3e50")
  )

ggsave("topL1_subfamilies_H3K27ac_5UTR_boxplot.pdf", width = 7.5, height = 5, dpi = 300)

# -------------------------------------------------------------------------------
# DNA methylation across LINE-1 subfamilies
# -------------------------------------------------------------------------------

target_subfamilies <- c("L1HS", paste0("L1PA", 2:7))

mcg_rep1 <- read.table("H1975-1_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
mcg_rep2 <- read.table("H1975-2_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
mcg_averaged <- mcg_rep1 %>%
  inner_join(mcg_rep2, by = c("LINE1","Subfamily"), suffix = c("_rep1","_rep2")) %>%
  mutate(Mean_mCG = (mCG_rep1 + mCG_rep2) / 2) %>%
  select(LINE1, Subfamily, Mean_mCG)

top100_mcg <- read.table("mergeL1_top100.sort.bed", header = FALSE, sep ="\t", stringsAsFactors = FALSE) %>%
  select(LINE1 = V4) %>%
  inner_join(mcg_averaged, by ="LINE1") %>%
  filter(Subfamily %in% target_subfamilies) %>%
  mutate(Subfamily = factor(Subfamily, levels = target_subfamilies))

all_l1_mcg <- read.table("L1.gt6k.bed", header = FALSE, sep ="\t", stringsAsFactors = FALSE) %>%
  select(LINE1 = V4) %>%
  inner_join(mcg_averaged, by ="LINE1") %>%
  filter(Subfamily %in% target_subfamilies) %>%
  mutate(Subfamily = factor(Subfamily, levels = target_subfamilies))

ggplot(top100_mcg, aes(x = Subfamily, y = Mean_mCG, fill = Subfamily)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, color ="#2c3e50", size = 0.6) +
  scale_fill_viridis_d(option ="plasma", end = 0.85) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x ="LINE-1 Subfamilies",
    y ="mCG Level (5' UTR)",
    title ="H1975 (Top 100 loci)"
  ) +
  theme_classic() +
  theme(
    legend.position ="none",
    plot.title = element_text(face ="bold", size = 13, hjust = 0.5, color ="black"),
    axis.title = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(face ="bold", size = 11, color ="black"),
    axis.text.y = element_text(size = 11, color ="black")
  )

ggsave("H1975_L1_Top100_mCG_boxplot.pdf", width = 7.5, height = 5, dpi = 300)

ggplot(all_l1_mcg, aes(x = Subfamily, y = Mean_mCG, fill = Subfamily)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, color ="#2c3e50", size = 0.6) +
  scale_fill_viridis_d(option ="plasma", end = 0.85) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x ="LINE-1 Subfamilies",
    y ="mCG Level (5' UTR)",
    title ="H1975 (All loci)"
  ) +
  theme_classic() +
  theme(
    legend.position ="none",
    plot.title = element_text(face ="bold", size = 13, hjust = 0.5, color ="black"),
    axis.title = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(face ="bold", size = 11, color ="black"),
    axis.text.y = element_text(size = 11, color ="black")
  )

ggsave("H1975_L1_All_mCG_boxplot.pdf", width = 7.5, height = 5, dpi = 300)

# -------------------------------------------------------------------------------
# Association between LINE-1 expression and Epigenomic signal
# -------------------------------------------------------------------------------

matrix_data <- read.table(
  gzfile("H1975_H3K27ac_L16k_5UTR.gz"),
  header = FALSE,
  sep ="\t",
  stringsAsFactors = FALSE,
  comment.char ="@"
)

matrix_data$H3K27ac_Mean <- rowMeans(matrix_data[, 7:ncol(matrix_data)], na.rm = TRUE)

H3K27ac_df <- data.frame(
  LINE1 = matrix_data$V4,
  Log2_H3K27ac = log2(matrix_data$H3K27ac_Mean + 1),
  stringsAsFactors = FALSE
)

exp_data <- read.table("H1975_L16k.TPM.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>%
  mutate(Log2_TPM = log2(((TPM1 + TPM2) / 2) + 1)) %>%
  .[, c("LINE1","Log2_TPM")] %>%
  dplyr::rename(LINE1 = LINE1)

all_loci_data <- H3K27ac_df %>%
  dplyr::inner_join(exp_data, by ="LINE1")

ggplot(all_loci_data, aes(x = Log2_H3K27ac, y = Log2_TPM)) +
  geom_point(color ="#2c3e50", alpha = 0.5, size = 1.2) +
  geom_smooth(method ="lm", color ="#2c3e50", fill ="#bdc3c7", alpha = 0.2, size = 0.8) +
  stat_cor(method ="pearson", fontface ="bold", label.x.npc ="left", label.y.npc ="top") +
  labs(
    x ="H3K27ac Signal (log2(Mean + 1))",
    y ="PEL-seq (log2(TPM + 1))",
    title ="H1975 (All loci)"  ) +
    theme_classic() +
    theme(
    plot.title = element_text(face ="bold", size = 12, hjust = 0.5),
    axis.title = element_text(face ="bold", size = 11),
    axis.text = element_text(face ="bold", size = 10)
    )

ggsave("H1975_L1_all_loci_H3K27ac_Expression.pdf", width = 5.5, height = 5, dpi = 300)

# -------------------------------------------------------------------------------
# Association between LINE-1 expression and DNA methylation
# -------------------------------------------------------------------------------

mcg_rep1 <- read.table("H1975-1_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>%
  dplyr::select(LINE1, mCG_1 = mCG)
mcg_rep2 <- read.table("H1975-2_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>%
  dplyr::select(LINE1, mCG_2 = mCG)
mcg_combined <- mcg_rep1 %>%
  dplyr::inner_join(mcg_rep2, by ="LINE1") %>%
  mutate(Mean_mCG = (mCG_1 + mCG_2) / 2)

exp_data <- read.table("H1975_L16k.TPM.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>%
  mutate(Log2_TPM = log2(((TPM1 + TPM2) / 2) + 1)) %>%
  .[, c("LINE1","Log2_TPM")] %>%
  dplyr::rename(LINE1 = LINE1)

all_loci_data <- mcg_combined %>%
  dplyr::inner_join(exp_data, by ="LINE1")

ggplot(all_loci_data, aes(x = Mean_mCG, y = Log2_TPM)) +
  geom_point(color ="#2c3e50", alpha = 0.5, size = 1.2) +
  geom_smooth(method ="lm", color ="#2c3e50", fill ="#bdc3c7", alpha = 0.2, size = 0.8) +
  stat_cor(method ="pearson", fontface ="bold", label.x.npc ="left", label.y.npc ="top") +
  labs(
    x ="mCG Level (mean ratio)",
    y ="PEL-seq (log2(TPM + 1))",
    title ="LINE-1 loci") +
   theme_classic() +
   theme(
    plot.title = element_text(face ="bold", size = 12, hjust = 0.5),
    axis.title = element_text(face ="bold", size = 11),
    axis.text = element_text(face ="bold", size = 10)
   )
ggsave("H1975_L1_all_loci_mCG_Expression.pdf", width = 5.5, height = 5, dpi = 300)

# -------------------------------------------------------------------------------
# Highly and lowly expressed LINE-1 loci
# -------------------------------------------------------------------------------

subfamilies <- c("L1HS","L1PA2","L1PA3")
base_dir <-"./"

expr_groups <- map_df(subfamilies, function(subf) {
  high_file <- file.path(base_dir, paste0(subf,"_highExpre_L1.bed"))
  low_file <- file.path(base_dir, paste0(subf,"_lowExpre_L1.bed"))
  df_high <- data.frame()
  df_low <- data.frame()
  if (file.exists(high_file)) {
    df_high <- read.table(high_file, header = FALSE, sep ="\t", stringsAsFactors = FALSE) %>%
      select(LINE1 = V4, Subfamily = V5) %>%
      mutate(Group ="highExpre_L1")
  }
  if (file.exists(low_file)) {
    df_low <- read.table(low_file, header = FALSE, sep ="\t", stringsAsFactors = FALSE) %>%
      select(LINE1 = V4, Subfamily = V5) %>%
      mutate(Group ="lowExpre_L1")
  }

  bind_rows(df_high, df_low)
})

tpm_file <- file.path(base_dir,"H1975_L16k.TPM.txt")

h1975_tpm <- read.table(tpm_file, header = TRUE, sep ="\t", stringsAsFactors = FALSE) %>%
  mutate(Signal = (TPM1 + TPM2) / 2) %>%
  select(LINE1 = LINE1, Signal)

filtered_data <- expr_groups %>%
  inner_join(h1975_tpm, by ="LINE1") %>%
  mutate(
    Subfamily = factor(Subfamily, levels = subfamilies),
    Group = factor(Group, levels = c("highExpre_L1","lowExpre_L1"))
  )
all_subf_data <- filtered_data
print(table(filtered_data$Subfamily, filtered_data$Group))

ggplot(filtered_data, aes(x = Group, y = Signal, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("highExpre_L1" ="#e74c3c","lowExpre_L1" ="#2c3e50")) +
  facet_wrap(~Subfamily, scales ="free_y", nrow = 1) +
  labs(
    x = NULL,
    y ="PEL-seq (TPM)",
    title =" "
  ) +
  theme_classic() +
  theme(
    legend.position ="none",
    strip.background = element_blank(),
    strip.text = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(hjust = 0.5)
  )

ggsave("H1975_L1_subfamilies_TPM_boxplot.pdf", width = 7, height = 5, dpi = 300)

p_value_table <- filtered_data %>%
  group_by(Subfamily) %>%
  summarise(
    p_value = wilcox.test(Signal ~ Group)$p.value
  )
print(p_value_table)
write.csv(p_value_table, file ="H1975_L1_subfamilies_TPM_p_value.csv", quote = FALSE, row.names = FALSE)

# -------------------------------------------------------------------------------
# Epigenomic Signal of highly and lowly expressed LINE-1 loci
# -------------------------------------------------------------------------------

subfamilies <- c("L1HS","L1PA2","L1PA3")
all_subf_data <- data.frame()

for (subf in subfamilies) {
  file_name <- paste0(subf,"_H3K27ac_5UTR.gz")

  if (file.exists(file_name)) {
    con <- gzfile(file_name,"rt")
    header_line <- readLines(con, n = 1)
    close(con)
    match <- regmatches(header_line, regexec('"group_boundaries":\\[([0-9,]+)\\]', header_line))
    if (length(match[[1]]) < 2) {
      warning(paste("File not found:", file_name))
      next
    }

    boundaries <- as.numeric(strsplit(match[[1]][2],",")[[1]])

    con_data <- gzfile(file_name,"rt")
    mat_data <- read.table(con_data, header = FALSE, comment.char ="@", stringsAsFactors = FALSE)
    close(con_data)

    mean_signals <- rowMeans(mat_data[, 7:ncol(mat_data)], na.rm = TRUE)
    high_idx <- 1:boundaries[2]
    low_idx <- (boundaries[2] + 1):boundaries[3]

    high_df <- data.frame(Signal = mean_signals[high_idx], Group ="highExpre_L1", Subfamily = subf)
    low_df <- data.frame(Signal = mean_signals[low_idx], Group ="lowExpre_L1", Subfamily = subf)

    all_subf_data <- rbind(all_subf_data, high_df, low_df)

  } else {
    warning(paste("File not found:", file_name))
  }
}

if (nrow(all_subf_data) == 0) {
  stop("No valid epigenomic signal files were found.")
}

all_subf_data$Subfamily <- factor(all_subf_data$Subfamily, levels = subfamilies)

filtered_data <- all_subf_data

ggplot(filtered_data, aes(x = Group, y = Signal, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("highExpre_L1" ="#e74c3c","lowExpre_L1" ="#2c3e50")) +
  facet_wrap(~Subfamily, scales ="free_y", nrow = 1) +
  labs(x = NULL, y ="H3K27ac Signal (5' UTR)", title =" ") +
  theme_classic() +
  theme(
    legend.position ="none",
    strip.background = element_blank(),
    strip.text = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(hjust = 0.5)
  )
ggsave("H1975_L1_subfamilies_5UTR_H3K27ac_from.pdf", width = 7, height = 5, dpi = 300)

p_value_table <- all_subf_data %>%
  group_by(Subfamily) %>%
  summarise(
    p_value = wilcox.test(Signal ~ Group)$p.value
  )
print(p_value_table)
write.csv(p_value_table, file ="H1975_L1_subfamilies_5UTR_H3K27ac_p_value.csv", quote = FALSE, row.names = FALSE)

# -------------------------------------------------------------------------------
# DNA methylation of highly versus lowly expressed LINE-1 loci
# -------------------------------------------------------------------------------

ctrl1 <- read.table("H1975-1_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)
ctrl2 <- read.table("H1975-2_L1_mCG.txt", header = TRUE, sep ="\t", stringsAsFactors = FALSE)

mCG_ctrl <- bind_rows(ctrl1, ctrl2) %>%
  group_by(LINE1) %>%
  summarize(mCG = mean(as.numeric(mCG), na.rm = TRUE), .groups ="drop") %>%
  mutate(Group ="H1975")

read_and_label <- function(file_path, subfamily, group_name) {
  read_tsv(file_path, col_names = FALSE, show_col_types = FALSE) %>%
    select(LINE1 = X4) %>%
    mutate(Subfamily = subfamily, Group = group_name)
}

all_beds <- bind_rows(
  read_and_label("L1HS_highExpre_L1.bed","L1HS","highExpre_L1"),
  read_and_label("L1HS_lowExpre_L1.bed","L1HS","lowExpre_L1"),
  read_and_label("L1PA2_highExpre_L1.bed","L1PA2","highExpre_L1"),
  read_and_label("L1PA2_lowExpre_L1.bed","L1PA2","lowExpre_L1"),
  read_and_label("L1PA3_highExpre_L1.bed","L1PA3","highExpre_L1"),
  read_and_label("L1PA3_lowExpre_L1.bed","L1PA3","lowExpre_L1")
)

mCG_ctrl_prep <- mCG_ctrl %>%
  rename(CellLine = Group, Signal = mCG)

all_subf_data <- all_beds %>%
  inner_join(mCG_ctrl_prep, by ="LINE1")

ggplot(all_subf_data, aes(x = Group, y = Signal, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c("highExpre_L1" ="#e74c3c","lowExpre_L1" ="#2c3e50")) +
  facet_wrap(~Subfamily, scales ="free_y", nrow = 1) +
  labs(x = NULL, y ="mCG Level (5' UTR)", title =" ") +
  theme_classic() +
  theme(
    legend.position ="none",
    strip.background = element_blank(),
    strip.text = element_text(face ="bold", size = 12, color ="black"),
    axis.text.x = element_text(hjust = 0.5)
  )

ggsave("H1975_L1_subfamilies_5UTR_mCG_boxplot.pdf", width = 7, height = 5, dpi = 300)

p_value_table <- all_subf_data %>%
  group_by(Subfamily) %>%
  summarise(
    p_value = wilcox.test(Signal ~ Group)$p.value
  )

print(p_value_table)
write.csv(p_value_table, file ="H1975_L1_subfamilies_5UTR_mCG_p_value.csv", quote = FALSE, row.names = FALSE)
