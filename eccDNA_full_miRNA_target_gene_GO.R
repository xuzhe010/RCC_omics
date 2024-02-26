#KCA miRNA target gene GO
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggupset)
library(stringr)
library(GOplot)

setwd("/home/xuzhe/project/ECC/KidneyCancer/")

target_gene <-openxlsx::read.xlsx("CircleSeq/miRNA-t-gene-miRTarBase.xlsx")

ego <- enrichGO(target_gene$Gene.Symbol, 
                   OrgDb=org.Hs.eg.db, 
                   ont='ALL',
                   pAdjustMethod='BH', 
                   pvalueCutoff=0.05, 
                   qvalueCutoff=0.2, 
                   keyType='SYMBOL')

#
egosimp <- simplify(ego,cutoff=0.7,by="p.adjust",
                       select_fun=min,measure="Wang")

edox <- setReadable(egosimp, 'org.Hs.eg.db', 'SYMBOL')

edox2 <- pairwise_termsim(edox)

treeplot(edox2)
