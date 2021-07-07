setwd("~/Desktop/L&D/expression")

library(gdata)
library(ggplot2)
library(plyr)
library(reshape2)
library(lme4)
library(asreml)
library(scales)
require(multcomp)
library(lmerTest)
library(car)
library(emmeans)

# combined candidate genes (to include vte7 and samt1)-------------------------------------------------------------
raw_genes=read.csv('../raw/light_dark_rna_fpkms_CPM_DEG_18june19_FPKM.csv')#FPKM
##get a priori genes
genes=read.xls('../tocochromanol_all_candidate_genes_combined_DW_20190624.xlsx')
raw_genes=merge(genes,raw_genes,by.x="RefGen_v4.Gene.ID" ,by.y="gene")
names(raw_genes)[1]='Gene_ID'

raw_genes_re=raw_genes[,-c(2:11)]
raw_genes_re.m=melt(raw_genes_re)

##reorder genes
raw_genes_re.m$Gene_ID=factor(raw_genes_re.m$Gene_ID,levels=unique(raw_genes$Gene_ID))

##reorder genotypes
all_geno=names(raw_genes[,-c(1:11)])
all_geno_sorted=all_geno[c(grep('B73',all_geno),grep('b97',all_geno),grep('M37W',all_geno),
                           grep('NC358',all_geno),grep('Ki11',all_geno),grep('MS71',all_geno),
                           grep('OH7B',all_geno))]

raw_genes_re.m$variable=factor(raw_genes_re.m$variable,levels=all_geno_sorted)

#format data for model fitting
data=data.frame(t(raw_genes_re))
names(data)=as.character(unlist(data[1,]))
data=data[-1,]

data$Rep=gsub('.*_.*_R','',row.names(data))
data$Treatment=gsub('.*_(.*)_.*','\\1',row.names(data))
data$genotype=gsub('_.*_.*','',row.names(data))

data[,c(1:130)] <- sapply(data[,c(1:130)], as.character)
data[,c(1:130)] <- sapply(data[,c(1:130)], as.numeric)
data[,c(1:130)] = data[,c(1:130)]+1e-9 #add a small constant to avoid 0 value before log
data[,c(1:130)] <- sapply(data[,c(1:130)], log2)#log2 transformation

genes=names(data)[1:130]

genes_pval=data.frame(matrix(nrow=0,ncol=3))
names(genes_pval)=c('genotype','treatment','genotype:treatment')
tukey_pval=data.frame(matrix(nrow=0,ncol=3))
names(tukey_pval)=c('dark_control','light_control','light_dark')
tukey_pval_sliced=data.frame(matrix(nrow=0,ncol=21))
names(tukey_pval_sliced)=c("B73_dark_control","B73_light_control","B73_light_dark","b97_dark_control","b97_light_control",
                           "b97_light_dark","Ki11_dark_control","Ki11_light_control","Ki11_light_dark","M37W_dark_control","M37W_light_control",
                           "M37W_light_dark","MS71_dark_control","MS71_light_control","MS71_light_dark","NC358_dark_control","NC358_light_control",
                           "NC358_light_dark","OH7B_dark_control","OH7B_light_control","OH7B_light_dark")

sink('genes_cand_genes_anova_tukey.txt')
for (i in genes){
  ##generate anova table
  cat(paste('#############################################################\n############ anova table for gene ',i,' ############\n#############################################################\n',sep=''))
  mod<-lmer(data[,i] ~ genotype + (1|Rep) + (1|Rep:genotype) + Treatment + genotype:Treatment , data) #singular fit
  print(Anova(mod,test='F'))
  cat('\n\n')
  #anova
  a=Anova(mod,test='F')
  genes_pval[i,]=a$`Pr(>F)`
  #tukey
  b=emmeans(mod, list(pairwise ~ Treatment), adjust = "tukey",lmer.df = "satterthwaite")
  b=as.data.frame(b$`pairwise differences of Treatment`)
  print(b)
  tukey_pval[i,]=b$p.value
  
  #nested tukey
  curr_pval=c()
  for(j in unique(data$genotype)){
    data_curr=data[which(data$genotype==j),]
    if(min(data_curr[,i])==max(data_curr[,i])){curr_pval=c(curr_pval,NA,NA,NA)}#no need to fit model if all values are the same
    else{
      mod<-lmer(data_curr[,i] ~ (1|Rep) + Treatment , data_curr) #singular fit
      c=emmeans(mod, list(pairwise ~ Treatment), adjust = "tukey",lmer.df = "satterthwaite")
      c=as.data.frame(c$`pairwise differences of Treatment`)
      print(c)
      curr_pval=c(curr_pval,c$p.value)
    }
  }
  tukey_pval_sliced[i,]=curr_pval
}

sink()
genes_pval$Gene_ID=row.names(genes_pval)
genes_pval=merge(raw_genes[,1:8],genes_pval,by='Gene_ID')

tukey_pval$Gene_ID=row.names(tukey_pval)
tukey_pval=merge(raw_genes[,1:8],tukey_pval,by='Gene_ID')

tukey_pval_sliced$Gene_ID=row.names(tukey_pval_sliced)
tukey_pval_sliced=merge(raw_genes[,1:8],tukey_pval_sliced,by='Gene_ID')

write.table(genes_pval,'1.genes_anova_pval.txt',quote=F,row.names = F,sep='\t')
write.table(tukey_pval,'1.genes_tukey_pval.txt',quote=F,row.names = F,sep='\t')
write.table(tukey_pval_sliced,'1.genes_tukey_sliced_pval.txt',quote=F,row.names = F,sep='\t')


