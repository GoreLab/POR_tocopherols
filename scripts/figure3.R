#Figure 3 plot for POR paper
library(tidyverse)
library(readxl)

#load in DESEQ output
#have to run DESeq2.R to obtain .csv - can speak the visualization parts
#go to line 296 and run from there
deseq <- read.csv("data/all_DEG_deseq2_res_ALL.csv")
#load in supplemental dataset
genes <- read_xlsx('data/Supplemental Dataset S1.xlsx')


split_parts <- sapply(strsplit(deseq$comparison, "_deseq2"), `[`, 1)

# Further split the first parts into genotype and treatment
genotype_treatment <- t(sapply(split_parts, function(x) strsplit(x, "_", fixed = TRUE)[[1]]))

# Create a data frame
df <- data.frame(
  Genotype = genotype_treatment[, 1],
  Treatment = apply(genotype_treatment,1, function(x) paste(x[-1], collapse = "_")),
  stringsAsFactors = FALSE
)

#naming for figure
deseq$geno <- df$Genotype
deseq$treatment_compare <- df$Treatment

#merge datasets
df  <- merge(genes,deseq,by.x="RefGen_v4.Gene.ID" ,by.y="gene")

#get means of log fold change per gene per treatment comparison
df_treatment_means <- df %>%
  group_by(treatment_compare, RefGen_v4.Gene.ID, A.Priori.Candidate.Gene.Pathway) %>%
  summarise(mean(log2FoldChange))

#Plotting format
df <- as.data.frame(df_treatment_means)
colnames(df)[3] <- "pathway" 
df$pathway <- str_to_title(df$pathway)

df$pathway  <- str_replace_all(df$pathway , "Ipp Synthesis", "IPP Synthesis")
df$pathway  <- str_replace_all(df$pathway , "3,8-Divinyl-Chlorophyllide Biosynthesis I \\(Aerobic, Light-Dependent\\)", "3,8-Divinyl-Chlorophyllide Biosynthesis I") 

df$treatment_compare  <- str_replace_all(df$treatment_compare , "dark_vs_control", "Control vs. Dark") 
df$treatment_compare  <- str_replace_all(df$treatment_compare , "dark_vs_light", "Light vs. Dark") 
df$treatment_compare  <- str_replace_all(df$treatment_compare , "light_vs_control", "Control vs. Light") 
df2 <- df[df$pathway %in% c("Chlorophyll Cycle", "Chlorophyll Degradation", "Tocochromanol Pathway","3,8-Divinyl-Chlorophyllide Biosynthesis I", "Prenyl Group Synthesis", "IPP Synthesis"),]
colnames(df2)[1] <- "Treatment Comparison" 
colnames(df2)[4] <- "Log Fold Change"
colnames(df2)[2] <- "gene"
df2$`Treatment Comparison` <- factor(df2$`Treatment Comparison`, levels = c("Control vs. Dark","Light vs. Dark","Control vs. Light"))
#highlight por genes
highlight_genes <- c("Zm00001d032576", "Zm00001d013937","Zm00001d046909")
df2 <- df2 %>%
  mutate(highlight = if_else(gene %in% highlight_genes, as.character(gene), 'Other'),
         point_size = if_else(gene %in% highlight_genes, 1.05, 1)) 
df2$highlight <- str_replace_all(df2$highlight ,"Zm00001d032576", "POR1")
df2$highlight <- str_replace_all(df2$highlight ,"Zm00001d013937", "POR2")
df2$highlight <- str_replace_all(df2$highlight ,"Zm00001d046909", "VTE2")

# plot and save
ggplot(df2, aes(x = `Treatment Comparison`, y = `Log Fold Change`) )+
  geom_boxplot(width = 0.5,outlier.shape = NA) +
  geom_point(aes(color = highlight,size = point_size),position = position_jitter(width = 0.2))+
  scale_color_manual(values = c('POR1' = "#882255", 'POR2' = "#E69F00","VTE2"="#56B4E9",   'Other' = 'grey'),breaks = c("POR1", "POR2", "VTE2")) +
  
  facet_wrap(~pathway) +
  
  labs( x = "Treatment Comparisons", y = "Gene Expression Log Fold Change") +
  theme_minimal() +
  labs(color = 'Gene name')+
  guides(size = "none")+
  theme(strip.text.x = element_text(angle = 0, size = 10, face = "bold"),axis.text.x = element_text(size = 12),legend.text = element_text(face = "italic"))
ggsave("Figure3v2_03_11_24.png",width = 10, height =5,units = "in" )
