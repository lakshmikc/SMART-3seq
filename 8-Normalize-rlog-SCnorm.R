#####Plan for Assessing RNAseq data 

library(DESeq2)
library(Rtsne)
library(ggplot2)
library(pheatmap)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(VennDiagram)
library(Rtsne)
library(ggplot2)
library(pcaMethods)
library(dplyr)


library(batchelor)
library(scater)
library(scran)

library(SCnorm)
library(sva)
library(limma)
library(Seurat)

library(preprocessCore)


prefix <- "/Users/Lakshmi/"
source(paste0(prefix,"/Dropbox (SBG)/SHR Work/NTS_analysis/qPCR_data_Functions11.R"))
`%notin%` <- Negate(`%in%`)

#############################
####Import Data and Setup####
#############################
##pull in our 90 x 15k matrix in both raw and rlog form 
#(Merged, Batch corrected, filtered matrix)
ragp.raw <- as.matrix(read.table(paste0(prefix,"Dropbox (SBG)/SPARC-Data-Acquisition/Downstream_Analysis/PIG_Data/Normalization/PR1643_adjusted_raw_90samps_15kgenes.txt"), sep="\t", header=T))

ragp.rlog <- rlog(ragp.raw, fitType = "local", blind = T); 
rownames(ragp.rlog) <- rownames(ragp.raw)
write.table(ragp.rlog, "PR1643_adjusted_rlog_90samps_15kgenes.txt", sep="\t", quote=F)

ragp.rlog <- as.matrix(read.table(paste0(prefix,"Dropbox (SBG)/SPARC-Data-Acquisition/Downstream_Analysis/PIG_Data/Normalization/PR1643_adjusted_rlog_90samps_15kgenes.txt"), sep="\t", header=T)); rownames(ragp.rlog) <- rownames(ragp.raw)

#Metadata for the filtered matrix
ragp.annots <- read.table("PR1643_90samples_annots.txt", sep="\t", header=T)
ragp.annots <- ragp.annots[match(colnames(ragp.raw), rownames(ragp.annots)),]
annot_samp <- data.frame(batch=ragp.annots[,1]); rownames(annot_samp) <- rownames(ragp.annots); annot_cols <- NA

#Median centering
ragp.rlog.med <- t(apply(ragp.rlog,1,function(x)(x-median(x))))

#Extract high expresison genes
#Using mean counts 
count.means <- apply(ragp.raw,1,mean)
highgenes <- rownames(ragp.raw)[which(count.means>10)]

#Using total counts
top200 <- rev(sort(rowSums(ragp.raw)))[1:200]
top2000 <- rev(sort(rowSums(ragp.raw)))[1:2000]


#######################################################################
#####################NORMALIZATION######################################
#######################################################################

####################
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

#DataNorm_k1 <- SCnorm(Data = ExampleSimSCData,Conditions = Conditions,PrintProgressPlots = TRUE,FilterCellNum = 10, K=1, NCores=4, reportSF = TRUE)
#DataNorm_k4 <- SCnorm(Data = ExampleSimSCData,Conditions = Conditions,PrintProgressPlots = TRUE,FilterCellNum = 10, NCores=4, K=4, reportSF = TRUE)
#DataNorm_k8 <- SCnorm(Data = ExampleSimSCData,Conditions = Conditions,PrintProgressPlots = TRUE,FilterCellNum = 10, NCores=4, K=8, reportSF = TRUE)
#DataNorm_k15 <- SCnorm(Data = ExampleSimSCData,Conditions = Conditions,PrintProgressPlots = TRUE,FilterCellNum = 10, NCores=4, K=15, reportSF = TRUE)
#NormalizedData.rlog <- SingleCellExperiment::normcounts(DataNorm)
#row.names(NormalizedData.rlog) <- row.names(ragp.rlog)
#NormalizedData.rlog.k4 <- SingleCellExperiment::normcounts(DataNorm_k4)
#row.names(NormalizedData.rlog.k4) <- row.names(ragp.rlog)
#NormalizedData.rlog.k15 <- SingleCellExperiment::normcounts(DataNorm_k15)
#row.names(NormalizedData.rlog.k15) <- row.names(ragp.rlog)


#K=20 18 for batch47 and K=9 17 for batch95
NormalizedData.rlog <- SingleCellExperiment::normcounts(DataNorm)

row.names(NormalizedData.rlog) <- row.names(ragp.rlog)
NormalizedData.rlog[1:5,1:5]

write.table(NormalizedData.rlog,"PR1643-normalized_90samples_15kgenes.txt", sep="\t", quote=F, col.names = NA, row.names = T)

#norm.rlog<- read.table("PR1643-90samples-rlog-norm.txt", sep="\t", header=T)

NData <- t(apply(NormalizedData.rlog,1,function(x)(x-median(x))))

ExtractGenes(NData, rownames(NData), show.rownames = F, main="")
hist(NData,breaks=1000, xlim=c(-5,5))

ragp.47 <- t(apply(ragp.rlog[,which(ragp.annots$batch=="Batch47")],1,function(x)(x-median(x))))
ragp.95 <- t(apply(ragp.rlog[,which(ragp.annots$batch=="Batch95")],1,function(x)(x-median(x))))

ragp.norm.47 <- t(apply(NormalizedData.rlog[,which(ragp.annots$batch=="Batch47")],1,function(x)(x-median(x))))
ragp.norm.95 <- t(apply(NormalizedData.rlog[,which(ragp.annots$batch=="Batch95")],1,function(x)(x-median(x))))

countDeptEst.SCNORM <- plotCountDepth(Data = ExampleSimSCData,NormalizedData = NormalizedData.rlog.k15,Conditions = Conditions,FilterCellProportion = .1, NCores=4)

ExtractGenes(NData,names(top200),main="Top 200 genes",show.rownames = F)
ExtractGenes(NData, toupper(gtex.genes), main="1882 Neuronal Genes from Gtex Found in Our Data", show.rownames = F)

##broad heatmap
scale.range <- c(1,4)
ExtractGenes(ragp.rlog, rownames(ragp.rlog), show.rownames = F)
ExtractGenes(ragp.rlog.med, rownames(ragp.rlog.med), show.rownames=F)

list <- c("^KCN","^SCN","^CACN")
ExtractGenes(ragp.rlog, list, exact=F)

###tsne of our samples 
ragp.tsne <- Rtsne(t(ragp.rlog), perplexity = 20)

plot(ragp.tsne$Y)

##try top 1000 genes 
ragp.tsne.top <- Rtsne(t(ragp.rlog[1:500,]), perplexity = 20)
plot(ragp.tsne.top$Y, col=ragp.annots$batch, pch=16)

ragp.tsne.top <- Rtsne(t(ragp.rlog[1:1000,]), perplexity = 10)
plot(ragp.tsne.top$Y, col=ragp.annots$batch, pch=20)


########################################
####Compare Data to Other Gene Lists####
########################################

###Compare our RNAseq data to other lists 
##Compare to Neuronal Genes from gtex --list number of genes found out of how many we're starting with 
gtex.genes <- read.table(paste0(prefix,"/Dropbox (SBG)/SPARC-manuscripts/Molecular Identity of rat ICN/Compare Neuro and Myo Affy/Gtex_BrainHighGenes.txt"),sep="\t"); gtex.genes <- as.character(gtex.genes$V1)
###some have a location and a gene name together? need to separate them out....
gtex.genes <- unlist(strsplit(gtex.genes,";"))

gtex.genes.common <- rownames(NData)[rownames(NData) %in% toupper(gtex.genes)]

ExtractGenes(NData, toupper(gtex.genes), main="1882 Neuronal Genes from Gtex Found in Our Data", show.rownames = F)

##examine our biomark list 
biomark <- t(read.table(paste0(prefix,"/Dropbox (SBG)/SPARC-Data-Acquisition/Pig 1643/PR1643/Analysis/PR1643_alldata_Raw.txt"), sep="\t"))

biomark.genes.common <- rownames(NData)[rownames(NData) %in% toupper(rownames(biomark))]
ExtractGenes(NData, toupper(rownames(biomark)))


###########################
####Compare to Hu et al####
###########################

###compare to variable genes in Hu et al.###

##need to upload 2 lists from Hu for now, their variable genes and the neuronal list I pulled out
hu.vars <- as.matrix(read.table(paste0(prefix,"Dropbox (SBG)/SPARC-Data-Acquisition/RNA Seq/2017-Hu-MolCell-scRNAseq-data/Scaled_Data_Variable_Genes_Hu.txt"),sep="\t",header=T,row.names=1))
hu.clusts <- read.table(paste0(prefix,"Dropbox (SBG)/SPARC-Data-Acquisition/RNA Seq/2017-Hu-MolCell-scRNAseq-data/Hu_Cluster_Assignments.txt"), sep="\t", header = T)
hu.neuro <- as.matrix(read.table(paste0(prefix,"Dropbox (SBG)/SPARC-Data-Acquisition/RNA Seq/2017-Hu-MolCell-scRNAseq-data/Scaled_Data_Important_Genes_Hu.txt"), sep="\t", header=T))
hu.tsne.info <- read.table(paste0(prefix,"Dropbox (SBG)/SPARC-Data-Acquisition/RNA Seq/2017-Hu-MolCell-scRNAseq-data/tSNE_coords_Group.txt"), sep="\t", header = T)
hu.color.info <- read.table(paste0(prefix,"Dropbox (SBG)/SPARC-Data-Acquisition/RNA Seq/2017-Hu-MolCell-scRNAseq-data/Hu_Cluster_Color_info.txt"), sep="\t", header = T)
#hu.genes <- read.table(paste0(prefix,"/Dropbox (SBG)/SPARC-manuscripts/Molecular Identity of rat ICN/Compare Affy to Brain Paper/ANOVA_significant_genes_AM.txt"), sep="\t", header=T); colnames(hu.genes) <- c("Gene","Pval")#neuro.genes <- neuro.genes$V1


#How many of their variable genes show up in our data? 
#Venn diagram showing overlap in variable genes and neuronal genes--could be interesting 

hu.genes.common <- rownames(ragp.rlog.med)[rownames(ragp.rlog.med) %in% toupper(rownames(hu.vars))]
ExtractGenes(ragp.rlog.med, toupper(rownames(hu.vars)), main="1223 Genes Varying in the CNS Found in Our Data")

##Comparison of variable genes in Hu vs our data 
venn.plot <- venn.diagram(x=list(CNS_Variable_Genes=hu.genes.common, Gtex_Neuronal_Genes=gtex.genes.common), 
                          filename = NULL, fill=c("blue","red"), alpha=0.2, main = paste("Cross Over Between Gtex Neuronal Genes and Variable Genes in CNS"))
grid.draw(rectGrob()); grid.draw(venn.plot)

#############################
####Scatter for Two Genes####
#############################

#plot Th vs Chat and Dbh for both sets, compare pattens--thisll be good! 
###RAGP Set 
gen1<- "CHAT"
gen2<- "TH"

dat1<-ragp.rlog[which(rownames(ragp.rlog) %in% gen1),]; 
dat2<-ragp.rlog[which(rownames(ragp.rlog) %in% gen2),]; 


dat.to.plot <- data.frame(Gene1= dat1, Gene2= dat2, Batch=ragp.annots$batch)


ggplot(dat.to.plot, aes(x=dat1,y=dat2,fill=Batch))+ geom_point(pch=21,color="black",size=5)  + 
  labs(x=paste(gen1), y= paste(gen2)) + ggtitle(paste(gen2, "vs.",gen1)) +
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40), 
                     axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = "none") #  + ylim(c(-5, 5))


####Hu et al Set
gen1<- "Dbh"
gen2<- "Th"

dat1<-hu.neuro[which(rownames(hu.neuro) %in% gen1),]; 
dat2<-hu.neuro[which(rownames(hu.neuro) %in% gen2),]; 

#dat1[which(dat1 < -5)] <- -5; dat2[which(dat2 < -5)] <- -5

dat.to.plot <- data.frame(Gene1= dat1, Gene2= dat2, Group=hu.tsne.info$Group)


ggplot(dat.to.plot, aes(x=dat1,y=dat2,fill=Group))+ geom_point(pch=21,color="black",size=5)  + 
  scale_fill_manual(values=hu.color.info$col) +labs(x=paste(gen1), y= paste(gen2)) + ggtitle(paste(gen2, "vs.",gen1)) +
  theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5, size=40), 
                     axis.text = element_text(size=25),axis.title = element_text(size=30), legend.position = "none") #  + ylim(c(-5, 5))



###tsne of just the variable genes from Hu et al in our data 
ragp.rlog.huvars <- ragp.rlog[rownames(ragp.rlog) %in% toupper(rownames(hu.vars)),]
ragp.hu.tsne <- Rtsne(t(ragp.rlog.huvars), perplexity = 20)

ragp.hu.tsne.info <- data.frame(tsne.X=ragp.hu.tsne$Y[,1], tsne.Y=ragp.hu.tsne$Y[,2], Batch=ragp.annots$batch)

ggplot(ragp.hu.tsne.info, aes(x=tsne.X, y=tsne.Y, fill=Batch)) +  geom_point(size=5, pch=21, col="black") + ggtitle("tSNE of RAGP Samples with CNS Variable Genes") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_text(hjust=0.5),
        axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank(), 
        legend.position = "none") #+ geom_text(data=label.info, aes(x=x.pos, y=y.pos,label=Group), col="black", fontface="bold")

##color for gene

gene <- "UCHL1"
rownames(ragp.rlog)[grep(paste(gene), rownames(ragp.rlog))]

gene.dat <- ragp.rlog[which(rownames(ragp.rlog) %in% gene),]

ragp.hu.tsne.gene.info <- ragp.hu.tsne.info; ragp.hu.tsne.gene.info$Gene <- (gene.dat)

ggplot(ragp.hu.tsne.gene.info, aes(x=tsne.X, y=tsne.Y, col=Gene)) +  geom_point(size=5) + scale_color_gradient(low="lightgrey",high="navy") + ggtitle(paste(gene)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank(), 
        legend.position = "right") #+ geom_text(data=label.info, aes(x=x.pos, y=y.pos,label=Group), col="black", fontface="bold")

###Extract variable genes from our data to compare 
##how to do this? Through PCA? Another method? Read paper to see how they did it, would Seurat work with this? 
##subset the Hu et al list for these genes and make a tsne to see how it compares for them 
##venn diagram of genes 

######################
####PCA of Samples####
######################

ragp.pca <- pca(t(ragp.rlog),scale="none",center=T,nPcs=90)
plot(ragp.pca)
ragp.scrs <- scores(ragp.pca)
ragp.ldgs <- loadings(ragp.pca)

#grab top from each
num.genes.to.get <- 50
num.pcs.to.examine <- 50
PC.genes <- c()

for (i in 1:num.pcs.to.examine) {
  PC.genes <- c(PC.genes,(rownames(ragp.ldgs)[order(ragp.ldgs[,i])])[1:num.genes.to.get],(rownames(ragp.ldgs)[rev(order(ragp.ldgs[,i]))])[1:num.genes.to.get])
}

PC.genes <- unique(PC.genes)

ExtractGenes(ragp.rlog.med, PC.genes, show.rownames = F, main="Genes with Most Variability Through Top 50 PCs")

##export PC genes to extract from Hu data

PC.for.Hu <- tools::toTitleCase(tolower(PC.genes))
write.table(PC.for.Hu, "PCGenes_for_Huetal.txt", sep="\t", quote=F)

####tsne of Samples with Only Variable Genes 
ragp.PC.data <- ragp.rlog[which(rownames(ragp.rlog) %in% PC.genes),]

ragp.PC.tsne <- Rtsne(t(ragp.PC.data), perplexity = 20)

ragp.PC.tsne.info <- data.frame(tsne.X=ragp.PC.tsne$Y[,1], tsne.Y=ragp.PC.tsne$Y[,2], Batch=ragp.annots$batch)

ggplot(ragp.PC.tsne.info, aes(x=tsne.X, y=tsne.Y, fill=Batch)) +  geom_point(size=5, pch=21, col="black") + ggtitle("tSNE of RAGP Samples Top Variable Genes") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_text(hjust=0.5),
        axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank(), 
        legend.position = "none") #+ geom_text(data=label.info, aes(x=x.pos, y=y.pos,label=Group), col="black", fontface="bold")



##Comparison of variable genes in Hu vs our data 
##our PC genes need to be compatible 
PC.genes1 <- 
  PC.genes1 <- PC.genes1[-grep("ENS", PC.genes1)]; PC.genes1 <- PC.genes1[-grep("LOC", PC.genes1)]
PC.genes1 <- toupper(unlist(strsplit(PC.genes,"_")))

venn.plot <- venn.diagram(x=list(CNS=toupper(rownames(hu.vars)), RAGP=PC.genes), 
                          filename = NULL, fill=c("blue","red"), alpha=0.2, main = paste("Cross Over Between Variable Genes in CNS and Variable Genes in the RAGP"))
grid.draw(rectGrob()); grid.draw(venn.plot)

length(grep("ENS", PC.genes1))
length(grep("LOC", PC.genes1))


####load PC genes in Hu data, tsne
Hu.RAGP.PC <- as.matrix(read.table("D://Dropbox (SBG)/SPARC-Data-Acquisition/RNA Seq/2017-Hu-MolCell-scRNAseq-data/Hu_et_al_RAGP_PC_genes_data.txt", sep="\t", header=T))

# hu.ragp.pc.tsne <- Rtsne(t(Hu.RAGP.PC))
# hu.ragp.tsne.info <- data.frame(tsne.X=hu.ragp.pc.tsne$Y[,1], tsne.Y=hu.ragp.pc.tsne$Y[,2], Group=hu.tsne.info$Group)
# write.table(hu.ragp.tsne.info, "Hu_tSNE_Coords_Info_RAGP_genes.txt", sep="\t", quote=F)
hu.ragp.tsne.info <- read.table("Hu_tSNE_Coords_Info_RAGP_genes.txt", sep="\t")

ggplot(hu.ragp.tsne.info, aes(x=tsne.X, y=tsne.Y, col=Group)) +  geom_point() + scale_color_manual(values=as.character(hu.cols$col)) + ggtitle("Highly Variable Genes in RAGP in CNS Dataset") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.title = element_text(hjust=0.5),
        axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank(), 
        legend.position = "none") #+ geom_text(data=label.info, aes(x=x.pos, y=y.pos,label=Group), col="black", fontface="bold")


###################
####tsne for Hu#### 
###################

hu.cols <- hu.color.info[match(levels(hu.clusts$x), hu.color.info$index.name),]


label.info <- hu.tsne.info %>% group_by(Group) %>% summarise(x.pos= median(tsne.X), y.pos=median(tsne.Y))

ggplot(hu.tsne.info, aes(x=tsne.X, y=tsne.Y, col=Group)) +  geom_point() + scale_color_manual(values=as.character(hu.cols$col)) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank(), 
        legend.position = "none") + geom_text(data=label.info, aes(x=x.pos, y=y.pos,label=Group), col="black", fontface="bold")


##color for specific gene

dat.to.use <- hu.neuro

gene <- "Hcn4"
rownames(dat.to.use)[grep(paste(gene), rownames(dat.to.use))]

gene.dat <- dat.to.use[which(rownames(dat.to.use) %in% gene),]

hu.tsne.gene.info <- hu.tsne.info; hu.tsne.gene.info$Gene <- t(gene.dat)

ggplot(hu.tsne.gene.info, aes(x=tsne.X, y=tsne.Y, col=Gene)) +  geom_point() + scale_color_gradient(low="lightgrey",high="navy") + ggtitle(paste(gene)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),axis.title.y=element_blank(),panel.background=element_blank(),
        panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank(), 
        legend.position = "right") #+ geom_text(data=label.info, aes(x=x.pos, y=y.pos,label=Group), col="black", fontface="bold")


