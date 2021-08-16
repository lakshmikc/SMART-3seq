library(Rsubread)


extension <- "sam"
fileNames <- Sys.glob(paste("*.", extension, sep=""))

savefolder  <- "/home/lck003/SPARC/PIG_Data/PR1643-97/FeatureCount"
for (fileName in fileNames) {

gene.result <- featureCounts(file=sprintf(fileName),annot.ext= "/home/lck003/SPARC/PIG_Data/Pig_genome/Sus_scrofa.Sscrofa11.1.95_clean.gtf",isGTFAnnotationFile=TRUE, GTF.featureType = "gene",readExtension3=200, ignoreDup = TRUE, isPairedEnd = FALSE,allowMultiOverlap = TRUE, nthreads = 20 )
exon.result <- featureCounts(file=sprintf(fileName),annot.ext= "/home/lck003/SPARC/PIG_Data/Pig_genome/Sus_scrofa.Sscrofa11.1.95_clean.gtf",isGTFAnnotationFile=TRUE, GTF.featureType = "exon",readExtension3=200, ignoreDup = TRUE, isPairedEnd = FALSE,allowMultiOverlap = TRUE, nthreads = 20 )
write.table(x=data.frame(gene.result$annotation,gene.result$counts), file = sprintf("%s/%s-gene-counts.txt",savefolder,fileName), quote=FALSE, sep="\t")
write.table(x=data.frame(exon.result$annotation,exon.result$counts), file = sprintf("%s/%s-exon-counts.txt",savefolder,fileName), quote=FALSE, sep="\t")

}
