
library(reshape2)
library(biomaRt)
library(DESeq2)
library(Rtsne)
library(ggplot2)
library(cluster)
library(gtools)
library(pheatmap)
library(pcaMethods)
library(rgl)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(VennDiagram)
library(limma)
library(sva)
library(edgeR)
library(SCnorm)
library(sva)

source("qPCR_data_Functions11.R")
source ("help_combat_seq.R")
source ("ComBat_seq.R")

#================================================================================================================
#Annotation - Genes and exons
#================================================================================================================
###############  For Exons #####################
raw <- read.table("PR1643-mm-exons-counts.txt",sep="\t", header = T, row.names = 1)
sums <- rowSums(raw)
raw <- raw[-which(sums==0),]

####fix gene names
genes <- rownames(raw)
ensembl <- useMart("ensembl")
ensembl <- useMart("ensembl",dataset="sscrofa_gene_ensembl", host="uswest.ensembl.org")

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","wikigene_description","wikigene_name"),values=genes,mart= ensembl)

raw2 <- raw; raw2$ensembl_gene_id <- rownames(raw2)
new.list <- G_list[-which(duplicated(G_list$ensembl_gene_id)),]
full.list <- right_join(new.list,raw2)

gene.names <- full.list$ensembl_gene_id
gene.names[which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name == full.list$wikigene_name)] <- full.list$external_gene_name[which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name == full.list$wikigene_name)]
gene.names[which(full.list$external_gene_name == "" & full.list$wikigene_name != "")] <- full.list$wikigene_name[which(full.list$external_gene_name == "" & full.list$wikigene_name != "")]
gene.names[which(full.list$external_gene_name != "" & full.list$wikigene_name == "")] <- full.list$external_gene_name[which(full.list$external_gene_name != "" & full.list$wikigene_name == "")]
gene.names[which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name != full.list$wikigene_name)] <- paste0(full.list$external_gene_name[
  which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name != full.list$wikigene_name)], "_" ,
  full.list$wikigene_name[which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name != full.list$wikigene_name)])

full.list$gene_name <- make.names(gene.names, unique = T)
full.list <- full.list[,c(1,ncol(full.list), 2:(ncol(full.list)-1))]

full.list2 <- full.list; full.list2[which(full.list=="", arr.ind = T)] <- NA
full.list2$description <- gsub(",","\\.", full.list2$description); full.list2$wikigene_description <- gsub(",","\\.", full.list2$wikigene_description)

write.csv(full.list2[,1:6], "PR1643_RNAseq_IDs_Exons_Description_rnd3.csv", row.names=F, quote=F)

full_dataset <- as.matrix(full.list2[,7:ncol(full.list2)]); rownames(full_dataset) <- full.list2$ensembl_gene_id

rownames(full_dataset) <- full.list2$gene_name
#rownames(full_dataset)[which(full.list$ensembl_gene_id != full.list$gene_name)] <- paste0(full.list$ensembl_gene_id[which(full.list$ensembl_gene_id != full.list$gene_name)],"_",full.list$gene_name[which(full.list$ensembl_gene_id != full.list$gene_name)])
write.table(full_dataset,"PR1643_raw_RNAseqData_exon_names_rnd3.txt",sep="\t",quote=F, col.names = NA)

###############  For Genes #####################
raw <- read.table("PR1643-mm-genes-counts.txt",sep="\t", header = T, row.names = 1)
sums <- rowSums(raw)
raw <- raw[-which(sums==0),]

####fix gene names
genes <- rownames(raw)
ensembl <- useMart("ensembl")
ensembl <- useMart("ensembl",dataset="sscrofa_gene_ensembl", host="uswest.ensembl.org")

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","wikigene_description","wikigene_name"),values=genes,mart= ensembl)

raw2 <- raw; raw2$ensembl_gene_id <- rownames(raw2)
new.list <- G_list[-which(duplicated(G_list$ensembl_gene_id)),]
full.list <- right_join(new.list,raw2)

gene.names <- make.names(full.list$ensembl_gene_id, unique = T)
gene.names[which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name == full.list$wikigene_name)] <- full.list$external_gene_name[which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name == full.list$wikigene_name)]
gene.names[which(full.list$external_gene_name == "" & full.list$wikigene_name != "")] <- full.list$wikigene_name[which(full.list$external_gene_name == "" & full.list$wikigene_name != "")]
gene.names[which(full.list$external_gene_name != "" & full.list$wikigene_name == "")] <- full.list$external_gene_name[which(full.list$external_gene_name != "" & full.list$wikigene_name == "")]
gene.names[which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name != full.list$wikigene_name)] <- paste0(full.list$external_gene_name[
  which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name != full.list$wikigene_name)], "_" ,
  full.list$wikigene_name[which(full.list$external_gene_name != "" & full.list$wikigene_name != "" & full.list$external_gene_name != full.list$wikigene_name)])

full.list$gene_name <- make.names(gene.names, unique=T)
full.list <- full.list[,c(1,ncol(full.list), 2:(ncol(full.list)-1))]
full.list2 <- full.list; full.list2[which(full.list=="", arr.ind = T)] <- NA
full.list2$description <- gsub(",","\\.", full.list2$description); full.list2$wikigene_description <- gsub(",","\\.", full.list2$wikigene_description)
write.csv(full.list2[,1:6], "PR1643_RNAseq_IDs_Genes_Description_rnd3.csv", row.names=F, quote=F)

full_dataset <- as.matrix(full.list2[,7:ncol(full.list2)]); rownames(full_dataset) <- full.list2$ensembl_gene_id
rownames(full_dataset) <- full.list2$gene_name
write.table(full_dataset,"PR1643_raw_RNAseqData_gene_names_rnd3.txt",sep="\t",quote=F, col.names = NA)



#================================================================================================================
#Batch correction: ComBat_seq
#================================================================================================================

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

PR1643.142.adjusted <- adjusted


#================================================================================================================
#Filtering Data
#================================================================================================================
##for every gene, look at total expression, and how many samples are non-zero

###Assess genes
samples.with.counts <- apply(PR1643.142.adjusted, 1, function(x)(sum(x!=0)))
tot.counts <- apply(PR1643.142.adjusted,1, function(x)(sum(x)))
average <- apply(PR1643.142.adjusted,1, function(x)(mean(x)))

stats <- data.frame(samples.with.counts,tot.counts,average)
ggplot(stats, aes(samples.with.counts)) + geom_histogram() + scale_x_log10() + scale_y_log10()
ggplot(stats, aes(1:142,samples.with.counts)) + geom_bar() + scale_x_log10() + scale_y_log10()

###Assess samples
samples.with.counts <- apply(PR1643.142.adjusted, 2, function(x)(sum(x!=0)))
tot.counts <- apply(PR1643.142.adjusted,2, function(x)(sum(x)))
average <- apply(PR1643.142.adjusted,2, function(x)(mean(x)))

stats <- data.frame(samples.with.counts,tot.counts,average)
ggplot(stats, aes(samples.with.counts)) + geom_histogram() #+ scale_x_log10() + scale_y_log10()
ggplot(stats, aes(tot.counts)) + geom_histogram() + scale_x_log10() + scale_y_log10()
ggplot(stats, aes(average)) + geom_histogram() + scale_x_log10() + scale_y_log10()

#First level filtering
PR1643.142.adjusted.filt <- PR1643.142.adjusted[,-which(samples.with.counts < 3000 & tot.counts < 100000)]

####make filtered sets -132 samples
write.table(PR1643.142.adjusted.filt, "PR1643-132-adjusted.txt", sep="\t", quote=F)

raw <- PR1643.142.adjusted.filt
rlog <- rlog(raw, fitType = "local", blind = T);
genes.with.counts <- apply(rlog, 1, function(x)(sum(x!=0)))
tot.counts <- apply(rlog,1, function(x)(sum(x)))
average <- apply(rlog,1, function(x)(mean(x)))
max <- apply(rlog,1, function(x)(max(x)))
stats.genes <- data.frame(genes.with.counts,tot.counts,average, max)
ggplot(stats.genes, aes(max)) + geom_histogram() + scale_y_log10()
ggplot(stats.genes, aes(1:142,samples.with.counts)) + geom_bar() + scale_x_log10() + scale_y_log10()
hist(apply(raw[which(average < 0),],1,mean))
hist(apply(raw[which(average < 0),],1,sum))

###propose to cut genes with max rld value of 0 or less

##Create log2 data
raw2 <- t(apply(raw+1,1,function(x)(log2(x))))
max.genes <- apply(raw2,1,max)
max.samples <- apply(raw2,2,max)

samples.with.counts.filt <- apply(raw, 2, function(x)(sum(x!=0)))
genes.with.counts.filt <- apply(raw, 1, function(x)(sum(x!=0)))

mat.ordered2 <- raw2[rev(order(genes.with.counts.filt)), rev(order(samples.with.counts.filt))]
mat.ordered3 <- mat.ordered2; mat.ordered3[which(mat.ordered2==0)] <- NA    ##still log2
raw.ordered2 <- raw[rev(order(genes.with.counts.filt)), rev(order(samples.with.counts.filt))]

#Filtering samples with non-zero gene counts <6000 and genes <30
filtered <- mat.ordered2[1:15000,1:90]
filtered.raw <- raw.ordered2[1:15000,1:90]
write.table(filtered.raw, "PR1643_adjusted_raw_90samps_15kgenes.txt", sep="\t", quote=F)

#================================================================================================================
#NORMALIZATION
#================================================================================================================

##pull in our 90 x 15k matrix in both raw and rlog form
#(Merged, Batch corrected, filtered matrix)
ragp.raw <- filtered.raw

#DEseq rlog normalization
####################

ragp.rlog <- rlog(ragp.raw, fitType = "local", blind = T);
rownames(ragp.rlog) <- rownames(ragp.raw)
#write.table(ragp.rlog, "PR1643_adjusted_rlog_90samps_15kgenes.txt", sep="\t", quote=F)

#Metadata for the filtered matrix
ragp.annots <- read.table("PR1643_90samples_annots.txt", sep="\t", header=T)
ragp.annots <- ragp.annots[match(colnames(ragp.raw), rownames(ragp.annots)),]
annot_samp <- data.frame(batch=ragp.annots[,1]); rownames(annot_samp) <- rownames(ragp.annots); annot_cols <- NA

#Median centering
ragp.rlog.med <- t(apply(ragp.rlog,1,function(x)(x-median(x))))

#Extract high expression genes
count.means <- apply(ragp.raw,1,mean)
highgenes <- rownames(ragp.raw)[which(count.means>10)]

#Using total counts
top200 <- rev(sort(rowSums(ragp.raw)))[1:200]
top2000 <- rev(sort(rowSums(ragp.raw)))[1:2000]

#SCnorm Normalization
####################

#Input is DESeq rlog matrix (90 samples)
ExampleSimSCData <- SingleCellExperiment::SingleCellExperiment(assays = list('counts' = ragp.rlog))
Conditions <- as.factor(annot_samp$batch)

countDeptEst <- plotCountDepth(Data = ExampleSimSCData, Conditions = Conditions,FilterCellProportion = .1, NCores=4)
str(countDeptEst)
head(countDeptEst[[1]])

ExampleSimSCData = SingleCellExperiment::counts(ExampleSimSCData)

#SCnormalization - DO not assign K
DataNorm <- SCnorm(Data = ExampleSimSCData,Conditions = Conditions,PrintProgressPlots = TRUE,FilterCellNum = 10, NCores=4, reportSF = TRUE)
DataNorm

#K=20 18 for batch47 and K=9 17 for batch95
NormalizedData.rlog <- SingleCellExperiment::normcounts(DataNorm)

row.names(NormalizedData.rlog) <- row.names(ragp.rlog)
NormalizedData.rlog[1:5,1:5]

write.table(NormalizedData.rlog,"PR1643-normalized_90samples_15kgenes.txt", sep="\t", quote=F, col.names = NA, row.names = T)



