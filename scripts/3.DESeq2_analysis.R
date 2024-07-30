#Light and Dark DEG 
#Author: Di Wu, Sam Herr, Mike Gore


#Load DESeq2
library(DESeq2)

# Processing the input data for DESeq2: Formating the input data ----------------------------------


#Read in the data
filename <- "../data/RNAseq_counts.csv" #file containing the read count matrix


cts <- as.matrix(read.csv(filename,row.names="gene", header = TRUE))
cts=cts[,-c(1:4)]#retain only counts, remove all other gene info
rownames(cts)
dim(cts)
cts_num=t(apply(cts, 1,function(x) as.numeric(as.character(x))))
colnames(cts_num)=colnames(cts)
dim(cts_num)
# Processing the input data for DESeq2: Removing Uninformative Rows -------------------------------


#Remove rows where there all samples don't have any counts.
#We're using the 'apply' function, which applys the user-specified function to either rows or columns
ind=c()

for (i in 1:nrow(cts_num)){
  if(all(cts_num[i,] == 0)){ind=c(ind,i)}
}
countsNonZero <- cts_num[-ind,]

#How many rows are left?
nrow(countsNonZero)
colnames(countsNonZero)


# Generating the coldata table -------------------------------------------------------------------

coldata=data.frame(matrix(colnames(countsNonZero)))
coldata['rep']=gsub('.*_.*_R','',as.character(coldata[,1]))
coldata['treatment']=gsub('.*_(.*)_.*','\\1',as.character(coldata[,1]))
coldata['genotype']=gsub('_.*_.*','',as.character(coldata[,1]))
rownames(coldata)=coldata[,1]
coldata$geno_trt <- paste(coldata$genotype, coldata$treatment, sep = "_")
coldata$geno_rep <- paste(coldata$genotype, coldata$rep, sep = "_")
rownames(coldata)

write.csv(coldata, "coldata.csv")



#Testing if the row names in the coldata table are exactly the same as the counts matrix
all(rownames(coldata) == colnames(countsNonZero))

# DESeq2 ---------------------------------------------------------------------------------------



#Building the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countsNonZero, 
                              colData=coldata, 
                              design= ~  geno_trt) #Design for DESeq

#Look at data structure of the DESeqDataSet object
str(dds)


# DESeq2: Sample QÇ ---------------------------------------------------------------

#Normalizing the counts for variance
vsd <- rlog(dds, blind = F) #have to use rlog to normalize sequencing data vst wont normalize for sequencing depth
vsd$matrix.colnames.countsNonZero..
#Extract the norm counts
norm.counts <- assay(vsd)
head(norm.counts[,1:6])



#save the normalized counts to a matrix
write.csv(norm.counts, file = "count_matrix_RLOG_blindF.csv", row.names = T, col.names = T)
#Calculate the Pearson's Correlation Coefficients for the normalized counts
counts.cor <- cor(norm.counts,method = "pearson")

head(counts.cor[,1:6])

#Load pheatmap package
library(pheatmap)

#Generate heatmap of the pcc
pdf('Correlation_normalized_counts.pdf')
pheatmap(counts.cor)
dev.off()

#Generating the PCA plot. Notice that the function takes the DESeq2 object 
#obtained after running the 'rlog' function, NOT the normalized counts directly. 
#See ?plotPCA for more infomration.
#below increase readability with labels 
library("RColorBrewer") 
library("ggplot2")
library("ggrepel")
pcaData <- plotPCA(vsd,intgroup=c("genotype", "treatment"), ntop = 5)

pdf('PCA.pdf')
plotPCA(vsd,intgroup=c("genotype","treatment"), ntop =20000)
pcaData + geom_label_repel(size = 2, aes(label = name)) # this is to label points on PCA, we have seen PCA. The default above just has hard to differentiate colors
dev.off()

# DESeq2: Calculate the test statistics -----------------------------------

#Calculate the dispersion parameters and the test statistics
dds <- DESeq(dds)

#Let's plot the dispersion to see how DESeq2 shrunk the values
pdf('dispersion.pdf')
plotDispEsts(dds)
dev.off()

#Lets look at the coefficient names in our object. 
resultsNames(dds)

# We are extracting the results for each of the contrasts below: 
dds.res.B73_light <- results(dds, contrast = c("geno_trt", "B73_light", "B73_control"), 
                             alpha = 0.05, format = "DataFrame")
dds.res.B73_dark <- results(dds, contrast = c("geno_trt", "B73_dark", "B73_control"), 
                            alpha = 0.05, format = "DataFrame")
dds.res.B73_dl <- results(dds, contrast = c("geno_trt", "B73_dark", "B73_light"), 
                          alpha = 0.05, format = "DataFrame")
dds.res.OH7B_light <- results(dds, contrast = c("geno_trt", "OH7B_light", "OH7B_control"), 
                              alpha = 0.05, format = "DataFrame")
dds.res.OH7B_dark <- results(dds, contrast = c("geno_trt", "OH7B_dark", "OH7B_control"), 
                             alpha = 0.05, format = "DataFrame")
dds.res.OH7B_dl <- results(dds, contrast = c("geno_trt", "OH7B_dark", "OH7B_light"), 
                           alpha = 0.05, format = "DataFrame")
dds.res.M37W_light <- results(dds, contrast = c("geno_trt", "M37W_light", "M37W_control"), 
                              alpha = 0.05, format = "DataFrame")
dds.res.M37W_dark <- results(dds, contrast = c("geno_trt", "M37W_dark", "M37W_control"), 
                             alpha = 0.05, format = "DataFrame")
dds.res.M37W_dl <- results(dds, contrast = c("geno_trt", "M37W_dark", "M37W_light"), 
                           alpha = 0.05, format = "DataFrame")
dds.res.MS71_light <- results(dds, contrast = c("geno_trt", "MS71_light", "MS71_control"), 
                              alpha = 0.05, format = "DataFrame")
dds.res.MS71_dark <- results(dds, contrast = c("geno_trt", "MS71_dark", "MS71_control"), 
                             alpha = 0.05, format = "DataFrame")
dds.res.MS71_dl <- results(dds, contrast = c("geno_trt", "MS71_dark", "MS71_light"), 
                           alpha = 0.05, format = "DataFrame")
dds.res.b97_light <- results(dds, contrast = c("geno_trt", "b97_light", "b97_control"), 
                             alpha = 0.05, format = "DataFrame")
dds.res.b97_dark <- results(dds, contrast = c("geno_trt", "b97_dark", "b97_control"), 
                            alpha = 0.05, format = "DataFrame")
dds.res.b97_dl <- results(dds, contrast = c("geno_trt", "b97_dark", "b97_light"), 
                          alpha = 0.05, format = "DataFrame")
dds.res.Ki11_light <- results(dds, contrast = c("geno_trt", "Ki11_light", "Ki11_control"), 
                              alpha = 0.05, format = "DataFrame")
dds.res.Ki11_dark <- results(dds, contrast = c("geno_trt", "Ki11_dark", "Ki11_control"), 
                             alpha = 0.05, format = "DataFrame")
dds.res.Ki11_dl <- results(dds, contrast = c("geno_trt", "Ki11_dark", "Ki11_light"), 
                           alpha = 0.05, format = "DataFrame")
dds.res.NC358_light <- results(dds, contrast = c("geno_trt", "NC358_light", "NC358_control"), 
                               alpha = 0.05, format = "DataFrame")
dds.res.NC358_dark <- results(dds, contrast = c("geno_trt", "NC358_dark", "NC358_control"), 
                              alpha = 0.05, format = "DataFrame")
dds.res.NC358_dl <- results(dds, contrast = c("geno_trt", "NC358_dark", "NC358_light"), 
                            alpha = 0.05, format = "DataFrame")

head(dds.res.B73_light)

#calculates LFC that are adjusted to account for low expression levels that can lead to inflated LFCs
dds.res.B73_light <- lfcShrink(dds, contrast=c("geno_trt", "B73_light", "B73_control"), res=dds.res.B73_light,type='normal')
dds.res.B73_dark <- lfcShrink(dds, contrast=c("geno_trt", "B73_dark", "B73_control"), res=dds.res.B73_dark,type='normal')
dds.res.B73_dl <- lfcShrink(dds, contrast=c("geno_trt", "B73_dark", "B73_light"), res=dds.res.B73_dl,type='normal')
dds.res.OH7B_light <- lfcShrink(dds, contrast=c("geno_trt", "OH7B_light", "OH7B_control"), res=dds.res.OH7B_light,type='normal')
dds.res.OH7B_dark <- lfcShrink(dds, contrast=c("geno_trt", "OH7B_dark", "OH7B_control"), res=dds.res.OH7B_dark,type='normal')
dds.res.OH7B_dl <- lfcShrink(dds, contrast=c("geno_trt", "OH7B_dark", "OH7B_light"), res=dds.res.OH7B_dl,type='normal')
dds.res.M37W_light <- lfcShrink(dds, contrast=c("geno_trt", "M37W_light", "M37W_control"), res=dds.res.M37W_light,type='normal')
dds.res.M37W_dark <- lfcShrink(dds, contrast=c("geno_trt", "M37W_dark", "M37W_control"), res=dds.res.M37W_dark,type='normal')
dds.res.M37W_dl <- lfcShrink(dds, contrast=c("geno_trt", "M37W_dark", "M37W_light"), res=dds.res.M37W_dl,type='normal')
dds.res.MS71_light <- lfcShrink(dds, contrast=c("geno_trt", "MS71_light", "MS71_control"), res=dds.res.MS71_light,type='normal')
dds.res.MS71_dark <- lfcShrink(dds, contrast=c("geno_trt", "MS71_dark", "MS71_control"), res=dds.res.MS71_dark,type='normal')
dds.res.MS71_dl <- lfcShrink(dds, contrast=c("geno_trt", "MS71_dark", "MS71_light"), res=dds.res.MS71_dl,type='normal')
dds.res.b97_light <- lfcShrink(dds, contrast=c("geno_trt", "b97_light", "b97_control"), res=dds.res.b97_light,type='normal')
dds.res.b97_dark <- lfcShrink(dds, contrast=c("geno_trt", "b97_dark", "b97_control"), res=dds.res.b97_dark,type='normal')
dds.res.b97_dl <- lfcShrink(dds, contrast=c("geno_trt", "b97_dark", "b97_light"), res=dds.res.b97_dl,type='normal')
dds.res.Ki11_light <- lfcShrink(dds, contrast=c("geno_trt", "Ki11_light", "Ki11_control"), res=dds.res.Ki11_light,type='normal')
dds.res.Ki11_dark <- lfcShrink(dds, contrast=c("geno_trt", "Ki11_dark", "Ki11_control"), res=dds.res.Ki11_dark,type='normal')
dds.res.Ki11_dl <- lfcShrink(dds, contrast=c("geno_trt", "Ki11_dark", "Ki11_light"), res=dds.res.Ki11_dl,type='normal')
dds.res.NC358_light <- lfcShrink(dds, contrast=c("geno_trt", "NC358_light", "NC358_control"), res=dds.res.NC358_light,type='normal')
dds.res.NC358_dark <- lfcShrink(dds, contrast=c("geno_trt", "NC358_dark", "NC358_control"), res=dds.res.NC358_dark,type='normal')
dds.res.NC358_dl <- lfcShrink(dds, contrast=c("geno_trt", "NC358_dark", "NC358_light"), res=dds.res.NC358_dl,type='normal')

# We're going to combine them all into a list so they are easier to work
# with
dds.res <- list(`B73_light` = dds.res.B73_light, `B73_dark` = dds.res.B73_dark, `B73_dl` = dds.res.B73_dl, 
                `OH7B_light` = dds.res.OH7B_light, `OH7B_dark` = dds.res.OH7B_dark, `OH7B_dl` = dds.res.OH7B_dl, 
                `M37W_light` = dds.res.M37W_light, `M37W_dark` = dds.res.M37W_dark,  `M37W_dl` = dds.res.M37W_dl, 
                `MS71_light` = dds.res.MS71_light, `MS71_dark` = dds.res.MS71_dark, `MS71_dl` = dds.res.MS71_dl, 
                `b97_light` = dds.res.b97_light, `b97_dark` = dds.res.b97_dark,  `b97_dl` = dds.res.b97_dl, 
                `Ki11_light` = dds.res.Ki11_light, `Ki11_dark` = dds.res.Ki11_dark, `Ki11_dl` = dds.res.Ki11_dl, 
                `NC358_light` = dds.res.NC358_light, `NC358_dark` = dds.res.NC358_dark, `NC358_dl` = dds.res.NC358_dl )

# The next step is filtering based on lfc and adj. p-value. First, we will
# make a function to remove NAs and then filter out genes with an adj.
# p-value > 0.05 and -2 < lfc > 2.
dds.filter <- function(x) {
  x <- x[!is.na(x$padj), ]
  y <- x[((x$padj < 0.05)), ]
  return(y)
}

# Then apply function to each element of the list
dds.res.flt <- lapply(dds.res, dds.filter)

x <- rbindlist(dds.res.flt)

# Look at some summary statistics for the results after filtering
summary(dds.res.flt)



#write ALL UNFILTERED results to a file for later use and further analysis 
write.csv(dds.res$`B73_light`, file = "B73_light_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`B73_dark`, file = "B73_dark_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`B73_dl`, file = "B73_dark_vs_light_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`OH7B_light`, file = "OH7B_light_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`OH7B_dark`, file = "OH7B_dark_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`OH7B_dl`, file = "OH7B_dark_vs_light_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`M37W_light`, file = "M37W_light_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`M37W_dark`, file = "M37W_dark_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`M37W_dl`, file = "M37W_dark_vs_light_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`MS71_light`, file = "MS71_light_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`MS71_dark`, file = "MS71_dark_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`MS71_dl`, file = "MS71_dark_vs_light_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`b97_light`, file = "B97_light_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`b97_dark`, file = "B97_dark_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`b97_dl`, file = "B97_dark_vs_light_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`Ki11_light`, file = "Ki11_light_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`Ki11_dark`, file = "Ki11_dark_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`Ki11_dl`, file = "Ki11_dark_vs_light_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`NC358_light`, file = "NC358_light_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`NC358_dark`, file = "NC358_dark_vs_control_deseq2_res_ALL.csv", quote = F)
write.csv(dds.res$`NC358_dl`, file = "NC358_dark_vs_light_deseq2_res_ALL.csv", quote = F)



files=list.files(pattern='_deseq2_res_ALL.csv')
files <- files[-1]
all_res=data.frame(matrix(nrow=0,ncol=8))
for (i in files){
  f=read.csv(i)
  names(f)[1]='gene'
  f$comparison=gsub('_deseq2_res.csv','',i)
  names(all_res)=names(f)
  all_res=rbind.data.frame(all_res,f)
}
#for figure3
write.csv(all_res,"all_DEG_deseq2_res_ALL.csv", quote = F,row.names = F)


# DESeq2: Plot the Log Fold Changes Against the Mean Counts ---------------

#Creat a function to plot horizontal lines at lfc of -2 and 2
drawLines <- function() abline(h=c(-2,2),col="dodgerblue",lwd=2)

#Plot the results
pdf('DEG_fold_change_with_mean_counts.pdf',width=10,height=5)
#Set the graphical parameters such that the plotting window will plot graphs together
par(mfrow=c(1,3))
plotMA(dds.res$`B73_light`,alpha = 0.05,main="B73_light_vs_control"); drawLines()
plotMA(dds.res$`B73_dark`,alpha = 0.05,main="B73_dark_vs_control"); drawLines()
plotMA(dds.res$`B73_dl`,alpha = 0.05,main="B73_dark_vs_light"); drawLines()
plotMA(dds.res$`OH7B_light`,alpha = 0.05,main="OH7B_light_vs_control"); drawLines()
plotMA(dds.res$`OH7B_dark`,alpha = 0.05,main="OH7B_dark_vs_control"); drawLines()
plotMA(dds.res$`OH7B_dl`,alpha = 0.05,main="OH7B_dark_vs_light"); drawLines()
plotMA(dds.res$`M37W_light`,alpha = 0.05,main="M37W_light_vs_control"); drawLines()
plotMA(dds.res$`M37W_dark`,alpha = 0.05,main="M37W_dark_vs_control"); drawLines()
plotMA(dds.res$`M37W_dl`,alpha = 0.05,main="M37W_dark_vs_light"); drawLines()
plotMA(dds.res$`MS71_light`,alpha = 0.05,main="MS71_light_vs_control"); drawLines()
plotMA(dds.res$`MS71_dark`,alpha = 0.05,main="MS71_dark_vs_control"); drawLines()
plotMA(dds.res$`MS71_dl`,alpha = 0.05,main="MS71_dark_vs_light"); drawLines()
plotMA(dds.res$`b97_light`,alpha = 0.05,main="B97_light_vs_control"); drawLines()
plotMA(dds.res$`b97_dark`,alpha = 0.05,main="B97_dark_vs_control"); drawLines()
plotMA(dds.res$`b97_dl`,alpha = 0.05,main="B97_dark_vs_light"); drawLines()
plotMA(dds.res$`Ki11_light`,alpha = 0.05,main="Ki11_light_vs_control"); drawLines()
plotMA(dds.res$`Ki11_dark`,alpha = 0.05,main="Ki11_dark_vs_control"); drawLines()
plotMA(dds.res$`Ki11_dl`,alpha = 0.05,main="Ki11_dark_vs_light"); drawLines()
plotMA(dds.res$`NC358_light`,alpha = 0.05,main="NC358_light_vs_control"); drawLines()
plotMA(dds.res$`NC358_dark`,alpha = 0.05,main="NC358_dark_vs_control"); drawLines()
plotMA(dds.res$`NC358_dl`,alpha = 0.05,main="NC358_dark_vs_light"); drawLines()
dev.off()
#Return the graphical parameters to normal
par(mfrow=c(1,1))





