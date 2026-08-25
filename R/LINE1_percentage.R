library(ggplot2)
library(dplyr)
library(reshape2)


Public_count <- read.csv("Public LINE1 count.csv",header = T)
Public_long <- melt(Public_count, id.vars = "sample")
Public_percent <- Public_long %>%
  group_by(sample) %>%
  mutate(percent = value / sum(value) * 100)
Public_percent$sample <- factor(Public_percent$sample,levels = c("RNA-seq_r1", "RNA-seq_r2", "RNA-seq_r3","ONT-seq_r1", "ONT-seq_r2", "ONT-seq_r3"))

ggplot(Public_percent, aes(x = sample, y = percent, fill = variable)) +
  geom_bar(stat = "identity",width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", percent)), 
            position = position_stack(vjust = 0.5), size = 7) +
  labs(x = " ", y = "Percentage", fill = "Category") +
  scale_fill_manual(values = c( 'aquamarine2','skyblue','snow2'))+
  theme_bw()+
  theme(legend.title = element_blank(),panel.grid = element_blank())+
  theme(axis.title.y=element_text(vjust=2, size=28))+
  theme(axis.text.y=element_text(vjust=1,size=28))+
  theme(axis.title.x=element_text(vjust=2, size=28))+
  theme(axis.text.x=element_text(vjust=0.6,size=25))+
  theme(legend.text = element_text(size = 25))
ggsave("Public LINE-1 percentage.pdf",width = 15,height = 10,dpi = 300) 



PEL_seq_count <- read.csv("PEL-seq L1 celltype.csv",header = T)
PEL_seq_long <- melt(PEL_seq_count, id.vars = "sample")
PEL_seq_percent <- PEL_seq_long %>%
  group_by(sample) %>%
  mutate(percent = value / sum(value) * 100)
PEL_seq_percent$sample <- factor(PEL_seq_percent$sample,levels = c("H1975-1","H1975-2","PC9-1","PC9-2","HCC827-1","HCC827-2"))

ggplot(PEL_seq_percent, aes(x = sample, y = percent, fill = variable)) +
  geom_bar(stat = "identity",width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", percent)), 
            position = position_stack(vjust = 0.5), size = 7) +
  labs(x = " ", y = "Percentage", fill = "Category") +
  scale_fill_manual(values = c( 'aquamarine2','skyblue','snow2'))+
  theme_bw()+
  theme(legend.title = element_blank(),panel.grid = element_blank())+
  theme(axis.title.y=element_text(vjust=2, size=28))+
  theme(axis.text.y=element_text(vjust=1,size=28))+
  theme(axis.title.x=element_text(vjust=2, size=28))+
  theme(axis.text.x=element_text(vjust=0.6,size=25))+
  theme(legend.text = element_text(size = 25))
ggsave("PEL-seq LINE-1 percentage.pdf",width = 15,height = 10,dpi = 300) 



