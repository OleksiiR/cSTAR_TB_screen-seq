#read in dataframe
# the file "GSE115905_TB_gene_exp.diff.txt" must be downloaded from GEO database
df <- read.delim("GSE115905_TB_gene_exp.diff.txt", sep = '\t')

#remove duplicate samples
df <- df[-c(96566, 190027, 236758, 610607, 750800, 797530, 937724, 1031185, 1077916, 1218109),]
#subset significant values only
df <- subset(df, significant=="yes")


#select columns required data for processing
keep <- c("gene","sample_1","sample_2","log2.fold_change.")
df <- df[keep]
df$sample <- paste(df$sample_2, df$sample_1, sep="_")
keep1 <- c("gene","sample","log2.fold_change.")
df1 <- df[keep1]
rownames(df1) <- make.names(df1[,"gene"], unique = TRUE)


#remove the infinite numbers 
is.na(df1)<-sapply(df1, is.infinite)
df1[is.na(df1)]<-0

#convert long to wide format
library(tidyr)
library(dplyr)
library(readr)
Log2FCWTvsMut <- pivot_wider(df1,
  names_from = sample,
  names_sep = ".",
  values_from = "log2.fold_change." 
)


#convert to matrix and write table
txt <- as.matrix(Log2FCWTvsMut)
write.table(txt, file = 'WTvsMutLog2FC.txt', sep = '\t', row.names = FALSE)



