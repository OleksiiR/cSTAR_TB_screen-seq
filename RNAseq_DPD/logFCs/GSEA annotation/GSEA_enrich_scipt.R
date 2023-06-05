# this script requires the file "msigdb_v2022.1.Mm_files_to_download_locally.zip"
# to be be downloaded from GSEA website and unpacked


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("fgsea")

library("readxl")
library("dplyr")
library(fgsea)
library(ggplot2)

pthr=0.05

gene_cont_df <- read_excel("genecontributions.xlsx")



# TMR vs TM
gene_set = gene_cont_df$TMR_TM_NV
names(gene_set) <- gene_cont_df$gene
gene_set = sort(gene_set,decreasing = TRUE)

# all signatures
GO_file="msigdb_v2022.1.Mm_GMTs/msigdb.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TMR_vs_TM_STV_GSEA_all.csv",row.names = FALSE)


# hallmark signatures
GO_file="msigdb_v2022.1.Mm_GMTs/mh.all.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TMR_vs_TM_STV_GSEA_MH.csv",row.names = FALSE)


# M2 signature
GO_file="msigdb_v2022.1.Mm_GMTs/m2.all.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TMR_vs_TM_STV_GSEA_M2.csv",row.names = FALSE)

# M5 GOBP signature
GO_file="msigdb_v2022.1.Mm_GMTs/m5.go.bp.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMR_vs_TM_STV_GSEA_M5_GOBP.csv",row.names = FALSE)

# M3 signature
GO_file="msigdb_v2022.1.Mm_GMTs/m3.all.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TMR_vs_TM_STV_GSEA_M3.csv",row.names = FALSE)

# M3 TFs signature
GO_file="msigdb_v2022.1.Mm_GMTs/m3.gtrd.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TMR_vs_TM_STV_GSEA_M3TFs.csv",row.names = FALSE)
#plotting
filtRes = head(res_signf, n = 6) %>% arrange(desc(NES))
filtRes$Effect = ifelse(filtRes$NES > 0, "TB resistance up", "TB resistance down")
colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("TB resistance up", "TB resistance down"))
g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Effect, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Transcription factors", y="Normalized Enrichment Score of contribution to TB resistance",
       title="Top transcription factors modified by RocA")
#pdf("TMR_TM_M3_TFs.pdf")
g1
#dev.off() 




# TM vs M
gene_set = gene_cont_df$TM_M_NV
names(gene_set) <- gene_cont_df$gene
gene_set = sort(gene_set,decreasing = TRUE)

# M3 TFs signature
GO_file="msigdb_v2022.1.Mm_GMTs/m3.gtrd.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TM_vs_M_STV_GSEA_M3TFs.csv",row.names = FALSE)
#plotting
filtRes = head(res_signf, n = 5) %>% arrange(desc(NES))
filtRes$Significance = ifelse(filtRes$padj < 0.05, "p < 0.05", "p > 0.05")
colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("p < 0.05", "p > 0.05"))
g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Significance, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Transcription factors", y="Normalized Enrichment Score of contribution to TB resistance",
       title="Top transcription factors modified by TNF")
pdf("TM_M_M3_TFs.pdf")
g1
dev.off() 

# M5 GOBP signature
GO_file="msigdb_v2022.1.Mm_GMTs/m5.go.bp.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TM_vs_M_STV_GSEA_M5_GOBP.csv",row.names = FALSE)

# M2 signature
GO_file="msigdb_v2022.1.Mm_GMTs/m2.all.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TM_vs_M_STV_GSEA_M2.csv",row.names = FALSE)

# all signatures
GO_file="msigdb_v2022.1.Mm_GMTs/msigdb.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TM_vs_M_STV_GSEA_all.csv",row.names = FALSE)

# MH signature
GO_file="msigdb_v2022.1.Mm_GMTs/mh.all.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(pval)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((pval)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TMR_vs_TM_STV_GSEA_M3TFs.csv",row.names = FALSE)
#plotting
filtRes = head(res_signf, n = 5) %>% arrange(desc(NES))
filtRes$Effect = ifelse(filtRes$NES > 0, "TB resistance up", "TB resistance down")
colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("TB resistance up", "TB resistance down"))
g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Effect, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Biological processes", y="Normalized Enrichment Score of contribution to TB resistance",
       title="Top hallmark genesets modified by TNF")
#pdf("TM_M_MH.pdf")
g1
#dev.off() 






# TMS vs TM
gene_set = gene_cont_df$TMS_TM_v_NV
names(gene_set) <- gene_cont_df$gene
gene_set = sort(gene_set,decreasing = TRUE)


# all signatures
GO_file="msigdb_v2022.1.Mm_GMTs/msigdb.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMS_vs_TM_STV_GSEA_all.csv",row.names = FALSE)

# M3 TFs signature
GO_file="msigdb_v2022.1.Mm_GMTs/m3.gtrd.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
#library(xlsx)
#write.xlsx(df,file = "genecont_GSEA.xlsx",col.names = TRUE, row.names = FALSE, append = FALSE)
write.csv(df,"TMS_vs_TM_STV_GSEA_M3TFs.csv",row.names = FALSE)
#plotting
filtRes = head(res_signf, n = 5) %>% arrange(desc(NES))
filtRes$Significance = ifelse(filtRes$padj < 0.05, "p < 0.05", "p > 0.05")
colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("p < 0.05", "p > 0.05"))
g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Significance, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Transcription factors", y="Normalized Enrichment Score of contribution to TB resistance",
       title="Top transcription factors modified by SP600125")
pdf("TMS_TM_M3_TFs.pdf")
g1
dev.off() 





# TMI vs TM
gene_set = gene_cont_df$TMI_TM_NV
names(gene_set) <- gene_cont_df$gene
gene_set = sort(gene_set,decreasing = TRUE)


# all signatures
GO_file="msigdb_v2022.1.Mm_GMTs/msigdb.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMI_vs_TM_STV_GSEA_all.csv",row.names = FALSE)


# MH signatures
GO_file="msigdb_v2022.1.Mm_GMTs/mh.all.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMI_vs_TM_STV_GSEA_MH.csv",row.names = FALSE)
#plotting
filtRes = head(res_signf, n = 5) %>% arrange(desc(NES))
filtRes$Significance = ifelse(filtRes$padj < 0.05, "p < 0.05", "p > 0.05")
colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("p < 0.05", "p > 0.05"))
g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Significance, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Biological processes", y="Normalized Enrichment Score of contribution to TB resistance",
       title="Top hallmark genesets modified by ISRIB")
g1


# M2 signatures
GO_file="msigdb_v2022.1.Mm_GMTs/m2.all.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMI_vs_TM_STV_GSEA_M2.csv",row.names = FALSE)


# M5 GOBP signature
GO_file="msigdb_v2022.1.Mm_GMTs/m5.go.bp.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMI_vs_TM_STV_GSEA_M5_GOBP.csv",row.names = FALSE)


# M3 TFs signature
GO_file="msigdb_v2022.1.Mm_GMTs/m3.gtrd.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMI_vs_TM_STV_GSEA_M3TFs.csv",row.names = FALSE)
#plotting
filtRes = head(res_signf, n = 5) %>% arrange(desc(NES))
filtRes$Significance = ifelse(filtRes$padj < 0.05, "p < 0.05", "p > 0.05")
colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("p < 0.05", "p > 0.05"))
g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Significance, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Transcription factors", y="Normalized Enrichment Score of contribution to TB resistance",
       title="Top transcription factors modified by ISRIB")
pdf("TMI_TM_M3_TFs.pdf")
g1
dev.off() 




# TMRS vs TMR_b
gene_set = gene_cont_df$TMS_TM_v_NV - gene_cont_df$TMR_b_TM_b_NV
names(gene_set) <- gene_cont_df$gene
gene_set = sort(gene_set,decreasing = TRUE)

# all signatures
GO_file="msigdb_v2022.1.Mm_GMTs/msigdb.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMRS_vs_TM_b_STV_GSEA_all.csv",row.names = FALSE)

# M3 TFs signature
GO_file="msigdb_v2022.1.Mm_GMTs/m3.gtrd.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMRS_vs_TMR_STV_GSEA_M3TFs.csv",row.names = FALSE)
#plotting
filtRes = head(res_signf, n = 5) %>% arrange(desc(NES))
filtRes$Significance = ifelse(filtRes$padj < 0.05, "p < 0.05", "p > 0.05")
colos = setNames(c("firebrick2", "dodgerblue2"),
                 c("p < 0.05", "p > 0.05"))
g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Significance, size = size), shape=21) +
  scale_fill_manual(values = colos ) +
  scale_size_continuous(range = c(2,10)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="Transcription factors", y="Normalized Enrichment Score of contribution to TB resistance",
       title="Top transcription factors modified by JNKi in presence of RocA")
pdf("TMRS_TMR_M3_TFs.pdf")
g1
dev.off() 




# TMR_b vs TM_b
gene_set = gene_cont_df$TMR_b_TM_b_NV
names(gene_set) <- gene_cont_df$gene
gene_set = sort(gene_set,decreasing = TRUE)

# all signatures
GO_file="msigdb_v2022.1.Mm_GMTs/msigdb.v2022.1.Mm.symbols.gmt"
myGO = fgsea::gmtPathways(GO_file)
# running GSEA
res_total = fgsea(pathways=myGO, stats=gene_set) %>% as.data.frame() %>% arrange(padj)
res_signf = as.data.frame(res_total %>% filter(pval < pthr) %>% arrange((padj)))
# saving results
df <- apply(res_signf,2,as.character)
write.csv(df,"TMR_b_vs_TM_b_STV_GSEA_all.csv",row.names = FALSE)

