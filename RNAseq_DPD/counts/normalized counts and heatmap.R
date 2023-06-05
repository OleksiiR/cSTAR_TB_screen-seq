
cts <- read.delim("readcount_genename.xls")
rownames(cts) <- make.names(cts[,"gene_name"], unique = TRUE)
TB_counts_raw <- as.matrix(cts[,2:40])


#BiocManager::install("sva")
library(sva)
batches<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
groups <-c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,3,3,3,4,4,4,7,7,7,8,8,8,9,9,9,10,10,10)
TB_counts <- ComBat_seq(as.matrix(TB_counts_raw),batch = batches,group=groups)


coldata <- read.csv("TB_codes.csv", row.names = "Code")

all(rownames(coldata) %in% colnames(TB_counts))
all(rownames(coldata) == colnames(TB_counts))


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = TB_counts,
                              colData = coldata,
                              design = ~ Group)
featureData <- data.frame(gene=rownames(TB_counts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts_TB <- counts(dds, normalized=TRUE)
write.table(normalized_counts_TB, file="normalized_counts_TB.txt", sep="\t", quote=F, col.names=NA)

