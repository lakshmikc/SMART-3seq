#Description: The sva package contains functions for removing batch effects and other unwanted variation in high-throughput
#experiment. Specifically, the sva package contains functions for the identifying and building surrogate variables for
#high-dimensional data sets. Surrogate variables are covariates constructed directly from high-dimensional data (like gene
#expression/RNA sequencing/methylation/brain imaging data) that can be used in subsequent analyses to adjust for unknown,
#unmodeled, or latent sources of noise. The sva package can be used to remove artifacts in three ways: (1) identifying and
##estimating surrogate variables for unknown sources of variation in high-throughput experiments (Leek and Storey 2007 PLoS
#Genetics,2008 PNAS), (2) directly removing known batch effects using ComBat (Johnson et al. 2007 Biostatistics) and (3) removing
#batch effects with known control probes (Leek 2014 biorXiv). Removing batch effects and using surrogate variables in
#differential expression analysis have been shown to reduce dependence, stabilize error rate estimates, and improve
#reproducibility, see (Leek and Storey 2007 PLoS Genetics, 2008 PNAS or Leek et al. 2011 Nat. Reviews Genetics).
#https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
#Batch correction: ComBat


library(DESeq2)
library(Rtsne)
library(ggplot2)
library(cluster)
library(gtools)
library(pheatmap)
library(pcaMethods)
library(reshape2)
library(rgl)

library(SCnorm)
library(limma)

library(sva)
library(edgeR)

setwd("/Users/Lakshmi/Dropbox (SBG)/SPARC-Data-Acquisition/Downstream_Analysis/PIG_Data/Batch-correction")
prefix <- "/Users/Lakshmi"
source ("/Users/Lakshmi/Dropbox (SBG)/SPARC-Data-Acquisition/Downstream_Analysis/PIG_Data/Batch-correction/help_combat_seq.R")
source ("/Users/Lakshmi/Dropbox (SBG)/SPARC-Data-Acquisition/Downstream_Analysis/PIG_Data/Batch-correction/ComBat_seq.R")
source(paste0(prefix,"/Dropbox (SBG)/SHR Work/NTS_analysis/qPCR_data_Functions11.R"))


full.all <- read.table("PR1643_raw_RNAseqData_gene_names_rnd3.txt", sep="\t", header=T, row.names = 1)
full <- full.all[which(rowSums(full.all)>0),]
batch47 <- full[,1:47]
batch95 <- full[,48:142]

new.set <- as.matrix(cbind(batch47,batch95))
col.data <- data.frame(batch=c(rep("Batch47",47), rep("Batch95",95)));
rownames(col.data) <- colnames(new.set)
annot_samp <- col.data
annot_cols <- NA


#############################
#Batch correction: ComBat_seq
#############################

counts <- as.matrix(new.set)
batch <- as.factor(c(rep("b47", 47), rep("b95", 95)))

group <- NULL
full_mod <- FALSE
covar_mod <- NULL
shrink <- FALSE
shrink.disp <- NULL
gene.subset.n <- NULL

adjusted<- ComBat_seq(counts, batch=batch)

write.table(adjusted,"FullSet_142_adjusted.txt", sep="\t", quote=F)

#Total percentage of zeros in samples
perc_zeroes <- colSums(adjusted==0)/nrow(adjusted)*100

#Ratio of total counts 
ratio_counts_sample <- colSums(adjusted)/nrow(adjusted)
ratio_counts_gene <- rowSums(adjusted)/ncol(adjusted)

adjusted.rlog.med <- t(apply(adjusted.rlog.count,1,function(x)(x-median(x)))); 
rownames(adjusted.rlog.med) <- rownames(adjusted.rlog.count)
write.table(adjusted.rlog.med,"PR1643-69-adjusted-rlog-med.txt", sep="\t", quote=F, col.names = NA, row.names = T)




#Reading the adjusted filtered rlog values

adjusted <- read.table("PR1643_adjusted_raw_90samps_15kgenes.txt", sep="\t", header=T, row.names = 1)
adjusted.rlog <- read.table("PR1643_adjusted_rlog_90samps_15kgenes.txt", sep="\t", header=T, row.names = 1)
rownames(adjusted.rlog) <- rownames(adjusted)

ragp.raw <- as.matrix(read.table("PR1643_adjusted_raw_90samps_15kgenes.txt"), sep="\t", header=T)
ragp.rlog <- as.matrix(read.table("PR1643_adjusted_rlog_90samps_15kgenes.txt"), sep="\t", header=T)

ragp.annots <- read.table("PR1643_90samples_annots.txt", sep="\t", header=T)
ragp.annots <- ragp.annots[match(colnames(ragp.raw), rownames(ragp.annots)),]
annot_samp <- data.frame(batch=ragp.annots[,1]); rownames(annot_samp) <- rownames(ragp.annots); annot_cols <- NA

batch <- as.factor(annot_samp$batch)
modcombat = model.matrix(~1, data= as.data.frame(adjusted))
dat <- as.matrix(adjusted)
combat_edata = ComBat(dat=dat, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

filt.combat2.rlog<- rlog(as.matrix(combat_edata), fitType = "local", blind = T) 
rownames(filt.combat2.rlog) <- rownames(combat_edata)


ragp.rlog.med <- t(apply(ragp.rlog,1,function(x)(x-median(x))))
scale.range <- c(1,4)
ExtractGenes(ragp.rlog, rownames(ragp.rlog), show.rownames = F)
ExtractGenes(ragp.rlog.med, rownames(ragp.rlog.med), show.rownames=F)







