library(ggsci)
library(ggprism)
library(tidyverse)
library(ggpubr)

#eccDNA 染色体分布-------------------------------------------------------------
ref <- read_tsv("Genome/window.1MB.bed",col_names = F)
ref <- ref %>% rename(Chr=X1,Start=X2,End=X3)

getReps<-function(x){
  df<-read_tsv(stringr::str_glue("preMB.bed/{x}_preMB.bed"),col_names = F)
  df <- df %>% rename(Chr=X1,Start=X2,End=X3)
  df <- ref %>% full_join(df)
  colnames(df)[4]<-{x}
  df[is.na(df)]<-0
  df
}

clinic = read_tsv("sample2R.tsv")

fin<-Reduce(\(x,y){full_join(x,y)},lapply(clinic$sample, getReps))

makeRI.count=function(g){
  count <- fin %>% select({g}) %>% 
    rowwise() %>% 
    summarise(Value=round(mean(c_across(is.numeric)))) 
  count <- cbind(ref,count)
  count$Chr <-gsub("chr","",count$Chr)
  count
}
tumor = clinic[which(clinic$group%in%"Tumor"),]$sample
normal = clinic[which(clinic$group%in%"Normal"),]$sample
case.circle=makeRI.count(tumor)
control.circle=makeRI.count(normal)

require("RIdeogram")
data(human_karyotype, package="RIdeogram")
#data(gene_density, package="RIdeogram")
#data(Random_RNAs_500, package="RIdeogram")
ideogram(karyotype = human_karyotype, 
         overlaid = case.circle, 
         label = control.circle,
         label_type = "heatmap",
         colorset1 = c("#E7F4E9", "#a50f15"),
         colorset2 = c("#F2E5F4", "#08519c"))

convertSVG("chromosome.svg", device = "png")
