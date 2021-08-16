##Analyzing RNAseq Data from SHR and WKY Brainstems 
library(reshape2)
library(biomaRt)

###############  For Exons #####################
raw <- read.table("D://Dropbox (SBG)/SPARC-Data-Acquisition/Upstream_Analaysis/PIG_Data/PR1643-47/RESULTS/Multimapping-SF/PR1643-mm-exons-counts.txt",sep="\t", header = T, row.names = 1)
sums <- rowSums(raw)
raw <- raw[-which(sums==0),]

####fix gene names
genes <- rownames(raw)
ensembl <- useMart("ensembl")
#datasets<- listDatasets(ensembl)
ensembl <- useMart("ensembl",dataset="sscrofa_gene_ensembl", host="uswest.ensembl.org")

#att <- listAttributes(ensembl)[,1:2]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","wikigene_description","wikigene_name"),values=genes,mart= ensembl)

# external.name.list <- getBM(filters="ensembl_gene_id", attributes = c("ensembl_gene_id","external_gene_name"), values=genes, mart=ensembl)
# wikigene.name.list <- getBM(filters="ensembl_gene_id", attributes = c("ensembl_gene_id","wikigene_name"), values=genes, mart=ensembl)
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
raw <- read.table("D://Dropbox (SBG)/SPARC-Data-Acquisition/Upstream_Analaysis/PIG_Data/PR1643-47/RESULTS/Multimapping-SF/PR1643-mm-genes-counts.txt",sep="\t", header = T, row.names = 1)
sums <- rowSums(raw)
raw <- raw[-which(sums==0),]

####fix gene names
genes <- rownames(raw)
ensembl <- useMart("ensembl")
#datasets<- listDatasets(ensembl)
ensembl <- useMart("ensembl",dataset="sscrofa_gene_ensembl", host="uswest.ensembl.org")

#att <- listAttributes(ensembl)[,1:2]
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name","description","wikigene_description","wikigene_name"),values=genes,mart= ensembl)

# external.name.list <- getBM(filters="ensembl_gene_id", attributes = c("ensembl_gene_id","external_gene_name"), values=genes, mart=ensembl)
# wikigene.name.list <- getBM(filters="ensembl_gene_id", attributes = c("ensembl_gene_id","wikigene_name"), values=genes, mart=ensembl)
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

#rownames(full_dataset)[which(full.list$ensembl_gene_id != full.list$gene_name)] <- paste0(full.list$ensembl_gene_id[which(full.list$ensembl_gene_id != full.list$gene_name)],"_",full.list$gene_name[which(full.list$ensembl_gene_id != full.list$gene_name)])
rownames(full_dataset) <- full.list2$gene_name

write.table(full_dataset,"PR1643_raw_RNAseqData_gene_names_rnd3.txt",sep="\t",quote=F, col.names = NA)


