## RNA-seq analysis with DESeq2
# Import & pre-process ----------------------------------------------------

library("DESeq2")
# Import data from featureCounts
## Previously ran at command line something like this:
## featureCounts -a genes.gtf -o counts.txt -T 12 -t exon -g gene_id GSM*.sam
# countdata <- read.table("static_clustering.txt", header=TRUE, row.names=1)
countdata <- read.table("ExpressionTable_DESEQ.txt", header=TRUE, row.names=1)

# Create a countdata for each sample comparison. Used ceiling becuase some number are not intengers to round these numbers.
 countdata1 <- ceiling(countdata[c(1,2,3,4,5,6)])
 countdata2 <- ceiling(countdata[c(1,2,3,7,8,9)])
 countdata3 <- ceiling(countdata[c(1,2,3,10,11,12)])
 countdata4 <- ceiling(countdata[c(1,2,3,13,14,15)])
 countdata5 <- ceiling(countdata[c(1,2,3,16,17,18)])
 countdata6 <- ceiling(countdata[c(1,2,3,19,20,21)])
 #print (countdata1)


# # For  libraries
 condition1 <- factor(c("control", "control", "control", "OX","OX","OX"))
 condition2 <- factor(c("control", "control", "control", "m1_4","m1_4","m1_4"))
 condition3 <- factor(c("control", "control", "control", "m123","m123","m123"))
 condition4 <- factor(c("control", "control", "control", "m134","m134","m134"))
 condition5 <- factor(c("control", "control", "control", "m3_4","m3_4","m3_4"))
 condition6 <- factor(c("control", "control", "control", "m4","m4","m4"))


# # Convert to matrix
countdata1 <- as.matrix(countdata1)
countdata2 <- as.matrix(countdata2)
countdata3 <- as.matrix(countdata3)
countdata4 <- as.matrix(countdata4)
countdata5 <- as.matrix(countdata5)
countdata6 <- as.matrix(countdata6)

# # Analysis with DESeq2 ----------------------------------------------------


# # Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
 (coldata1 <- data.frame(row.names=colnames(countdata1), condition1))
 (coldata2 <- data.frame(row.names=colnames(countdata2), condition2))
 (coldata3 <- data.frame(row.names=colnames(countdata3), condition3))
 (coldata4 <- data.frame(row.names=colnames(countdata4), condition4))
 (coldata5 <- data.frame(row.names=colnames(countdata5), condition5))
 (coldata6 <- data.frame(row.names=colnames(countdata6), condition6))


dds1 <- DESeqDataSetFromMatrix(countData=countdata1, colData=coldata1, design=~condition1)
dds2 <- DESeqDataSetFromMatrix(countData=countdata2, colData=coldata2, design=~condition2)
dds3 <- DESeqDataSetFromMatrix(countData=countdata3, colData=coldata3, design=~condition3)
dds4 <- DESeqDataSetFromMatrix(countData=countdata4, colData=coldata4, design=~condition4)
dds5 <- DESeqDataSetFromMatrix(countData=countdata5, colData=coldata5, design=~condition5)
dds6 <- DESeqDataSetFromMatrix(countData=countdata6, colData=coldata6, design=~condition6)

dds1$condition1 <- relevel(dds1$condition1, ref="control")
dds2$condition2 <- relevel(dds2$condition2, ref="control")
dds3$condition3 <- relevel(dds3$condition3, ref="control")
dds4$condition4 <- relevel(dds4$condition4, ref="control")
dds5$condition5 <- relevel(dds5$condition5, ref="control")
dds6$condition6 <- relevel(dds6$condition6, ref="control")

# # Run the DESeq pipeline
dds1 <- DESeq(dds1)
dds2 <- DESeq(dds2)
dds3 <- DESeq(dds3)
dds4 <- DESeq(dds4)
dds5 <- DESeq(dds5)
dds6 <- DESeq(dds6)

# # Remove genes with fewer than 1 read
dds1 <- dds1[ rowSums(counts(dds1)) > 1, ]
dds2 <- dds2[ rowSums(counts(dds2)) > 1, ]
dds3 <- dds3[ rowSums(counts(dds3)) > 1, ]
dds4 <- dds4[ rowSums(counts(dds4)) > 1, ]
dds5 <- dds5[ rowSums(counts(dds5)) > 1, ]
dds6 <- dds6[ rowSums(counts(dds6)) > 1, ]

# # Create MA plots of each analysis
png('MAplot_control_OX.png')
plotMA(dds1,ylim=c(-2,2), main='DESeq2')
dev.off()

png('MAplot_control_m1_4.png')
plotMA(dds2,ylim=c(-2,2), main='DESeq2')
dev.off()

png('MAplot_control_m123.png')
plotMA(dds3,ylim=c(-2,2), main='DESeq2')
dev.off()

png('MAplot_control_m134.png')
plotMA(dds4,ylim=c(-2,2), main='DESeq2')
dev.off()

png('MAplot_control_m3_4.png')
plotMA(dds5,ylim=c(-2,2), main='DESeq2')
dev.off()

png('MAplot_control_m4.png')
plotMA(dds6,ylim=c(-2,2), main='DESeq2')
dev.off()

# # Plot dispersions
png("qc-dispersions_MAplot_control_OX.png", 1000, 1000, pointsize=20)
plotDispEsts(dds1, main="Dispersion plot")
dev.off()
png("qc-dispersions_control_m1_4.png", 1000, 1000, pointsize=20)
plotDispEsts(dds2, main="Dispersion plot")
dev.off()
png("qc-dispersions_control_m123.png", 1000, 1000, pointsize=20)
plotDispEsts(dds3, main="Dispersion plot")
dev.off()
png("qc-dispersions_control_m134.png", 1000, 1000, pointsize=20)
plotDispEsts(dds4, main="Dispersion plot")
dev.off()

png("qc-dispersions_control_m3_4.png", 1000, 1000, pointsize=20)
plotDispEsts(dds5, main="Dispersion plot")
dev.off()

png("qc-dispersions_control_m4.png", 1000, 1000, pointsize=20)
plotDispEsts(dds6, main="Dispersion plot")
dev.off()


# Get differential expression results
res1 <- results(dds1)
table(res1$padj<0.05)
res1 <- res1[order(res1$padj), ]
res2 <- results(dds2)
table(res2$padj<0.05)
res2 <- res2[order(res2$padj), ]
res3 <- results(dds3)
table(res3$padj<0.05)
res3 <- res3[order(res3$padj), ]
res4 <- results(dds4)
table(res4$padj<0.05)
res4 <- res4[order(res4$padj), ]
res5 <- results(dds5)
table(res5$padj<0.05)
res5 <- res5[order(res5$padj), ]
res6 <- results(dds6)
table(res6$padj<0.05)
res6 <- res6[order(res6$padj), ]

# Order by adjusted p-value
# Merge with normalized count data
res1data <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
res2data <- merge(as.data.frame(res2), as.data.frame(counts(dds2, normalized=TRUE)), by="row.names", sort=FALSE)
res3data <- merge(as.data.frame(res3), as.data.frame(counts(dds3, normalized=TRUE)), by="row.names", sort=FALSE)
res4data <- merge(as.data.frame(res4), as.data.frame(counts(dds4, normalized=TRUE)), by="row.names", sort=FALSE)
res5data <- merge(as.data.frame(res5), as.data.frame(counts(dds5, normalized=TRUE)), by="row.names", sort=FALSE)
res6data <- merge(as.data.frame(res6), as.data.frame(counts(dds6, normalized=TRUE)), by="row.names", sort=FALSE)

# Write results
write.csv(res1data, file="diffexpr-results_control_OX.csv")
write.csv(res2data, file="diffexpr-result_control_m1_4.csv")
write.csv(res2data, file="diffexpr-result_control_m123.csv")
write.csv(res3data, file="diffexpr-results_control_m134.csv")
write.csv(res4data, file="diffexpr-results_control_m3_4.csv")
write.csv(res5data, file="diffexpr-results_control_m4.csv")
