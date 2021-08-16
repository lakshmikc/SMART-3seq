# raj said look through samples with least number of genes, look at genes that are present in few samples, mark and keep, average or total counts? 
#   
#   Neuronal markers......key markers (samples that do poorly for them?)
# 
# once we filter the samples we can look at genes 
# 
# one of the things he recommended.....for every single sample, make a characteristic group of stats, percentage of key genes etc 


##for every gene, look at total expression, and how many samples are non-zero

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


PR1643.142.adjusted.filt <- PR1643.142.adjusted[,-which(samples.with.counts < 3000 & tot.counts < 100000)]

samples.with.counts.filt <- apply(PR1643.142.adjusted.filt, 2, function(x)(sum(x!=0)))
tot.counts.filt <- apply(PR1643.142.adjusted.filt,2, function(x)(sum(x)))
average.filt <- apply(PR1643.142.adjusted.filt,2, function(x)(mean(x))) 

stats.filt <- data.frame(samples.with.counts.filt,tot.counts.filt,average.filt)

stats.less <- stats[which(tot.counts < 100000),]  ##remaining samples with less than 100,000 count, keep all for now 


neuro.markers <- c("UCHL1", "RBFOX3", "MAP2", "CAMK2A","CAMK2D","CAMK2B", "TUBB","TUBB1","TUBB4A","TUBB4B", "TUBB6")

neuro.genes <- PR1643.142.adjusted.filt[which(rownames(PR1643.142.adjusted.filt) %in% neuro.markers),]

neuro.genes$Gene <- rownames(neuro.genes)

melted <- melt(neuro.genes) 

ggplot(melted, aes(x=variable, y=value, group=Gene, col=Gene)) + geom_line() + facet_wrap(~Gene, scales="free_y")
ggplot(melted, aes(x=variable, y=value, group=Gene, col=Gene)) + geom_line() + ylim(0,10)


hist(colSums(neuro.genes[,1:132]))

filt <- neuro.genes[,which(colSums(neuro.genes[,1:131]) < 20)]  ##samples with less than 20 counts across neuronal markers

less.and.neuro <- stats.less[rownames(stats.less) %in% colnames(filt),]  ##samples with les than 20 counts and less than 100000 total counts 

a <- cbind(stats[rownames(stats) %in% colnames(filt),], Neuro.Sum=colSums(filt))
plot(a$Neuro.Sum, a$tot.counts, xlab="NeuroMarker.Sum",ylab="Total Counts")

####make filtered sets 

# PR1643.142.adjusted.filt <- PR1643.142.adjusted[,-which(samples.with.counts < 3000 & tot.counts < 100000)]  ##same line as above
write.table(PR1643.142.adjusted.filt, "PR1643-132-adjusted.txt", sep="\t", quote=F)

PR1643.142.adjusted.rlog.filt <- PR1643.142.adjusted.rlog[,-which(samples.with.counts < 3000 & tot.counts < 100000)]

write.table(PR1643.142.adjusted.rlog.filt, "PR1643-132-adjusted-rlog.txt", sep="\t", quote=F)








#####Filter Genes

rlog <- read.table("PR1643-132-adjusted-rlog.txt", sep="\t", header=T)
raw <- read.table("PR1643-132-adjusted.txt", sep="\t", header=T)

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

##log2 data 
raw2 <- t(apply(raw+1,1,function(x)(log2(x))))
max.genes <- apply(raw2,1,max)
max.samples <- apply(raw2,2,max)

head(rev(sort(max.genes)))
head(rev(sort(max.samples)))


ExtractGenes(raw2, rownames(raw2)[1:2000], show.rownames = F, order.by.gene = "ENSSSCG00000027279", order.by.sample = "PR1643_728.G.2")

mat.ordered <- raw2[rev(order(raw2[,which(colnames(raw2)=="PR1643_728.G.2")])), rev(order(raw2[which(rownames(raw2)=="PSD3"),]))]

write.table(mat.ordered, "PR1643_132_adjusted-log2-Ordered.txt",sep="\t", quote=F)


samples.with.counts.filt <- apply(raw, 2, function(x)(sum(x!=0)))
genes.with.counts.filt <- apply(raw, 1, function(x)(sum(x!=0)))

mat.ordered2 <- raw2[rev(order(genes.with.counts.filt)), rev(order(samples.with.counts.filt))]
write.table(mat.ordered, "PR1643_adjusted_log2-ordered_by_zeros.txt", sep="\t", quote=F)

mat.ordered3 <- mat.ordered2; mat.ordered3[which(mat.ordered2==0)] <- NA    ##still log2 


raw.ordered2 <- raw[rev(order(genes.with.counts.filt)), rev(order(samples.with.counts.filt))]
write.table(raw.ordered2, "PR1643_raw_ordered_by_zeros.txt", sep="\t", quote=F)

filtered <- mat.ordered2[1:15000,1:90]

filtered.raw <- raw.ordered2[1:15000,1:90]

write.table(filtered.raw, "PR1643_adjusted_raw_90samps_15kgenes.txt", sep="\t", quote=F)

