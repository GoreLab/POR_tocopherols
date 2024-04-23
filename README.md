PREPRINT:
This repository and website documents all analyses, summary, tables and figures associated with the following PREPRINT: 
"Tocopherol concentration in maize grain is dependent on the chlorophyll biosynthesis within the embryo"

DIRECTORIES:

data: contains data relevant to analysis and visualizations
"figure2_data.xlsx" -  metabolomics data harvested from 24 Days after pollination embryos from the 2018 physiological experiment
"figure5a_data.xlsx" - metabolomics data harvested from mature kernels selected from CRISPR/Cas9 mutants in the 2019 field experiment
"figure5b_data.csv" - metabolomics data harvested from mature kernels selected from CRISPR/Cas9 mutants in the 2020 Greenhouse experiment 
"figure6_data.xlsx" - metabolomics data harvested from 24 Days after pollination embryos selected from selfed hemizygous por1 and homozygous por2 mutant in the 2020 Greenhouse experiment 
"figureS2.xlsx" -  metabolomics data harvested from mature kernels from the 2018 physiological experiment
"figureS3.xlsx" -  metabolomics data harvested from mature kernels from the 2019 physiological experiment
"RNAseq_counts.csv" - RNAseq counts data collected from 24 Days after pollination embryos from the 2018 physiological experiment - used for DESeq2 analysis

scripts: 
"1.light_dark_analysis.Rmd" - Analysis of the physiological experiment. Obtains BLUEs and performs outlier detection, ANOVA, and Tukey's HSD.
                            - input: "figure2_data.xlsx", "figureS2.xlsx", "figureS3.xlsx"
"2.Figure2_visualization.Rmd" - Script to output Figure 2.
                              - input: output of "1.light_dark_analysis.Rmd" 
"3.DESeq2_analysis.R" - Performs DESeq2 analysis on all genotypes from the 2018 physiological experiment
                      - input: "RNAseq_counts.csv"
"4.Figure3_visualization.R" - Script to output Figure 3
                            - input: output of "3.DESeq2_analysis.R" 
"5.Figure5_analysis_visualization.Rmd" - Analysis and visualization of mature kernels metabolomics from the 2019 field experiment & 2020 Greenhouse experiment 
                                       - input: "figure5a_data.xlsx", "figure5b_data.csv"
"6.Figure6_analysis_visualization.Rmd" - Analysis and visualization of 24 DAP embryo metabolomics from the 2020 Greenhouse experiment 
                                       - input: "figure6_data.xlsx"
                                  
                   
                            
