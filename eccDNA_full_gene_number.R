library(ggplot2)
library(ggprism)
library(ggpubr)
library(tidyverse)
library(openxlsx)

clinc = read_tsv("KidneyCancer/sample/sample2R.tsv")

get_full_num<-function(x){
  df<-read_tsv(stringr::str_glue("KidneyCancer/CircleSeq/raw.count/fullGene.bed/{x}_coding_full.bed"),col_names = F)
  df<-df %>% distinct(X1,X2,X3)
  df1 = data.frame(num=nrow(df))
  df1$sample=x
  df1
}

eccnumbers<-do.call("bind_rows",lapply(clinc$sample,get_full_num))

sum(eccnumbers$num)

eccnumbers = left_join(eccnumbers,clinc)

df = eccnumbers %>% left_join(clinc)

#full gene number

#table(df$group)

ggplot(df,aes(x=group,y=num,fill=group))+
  stat_boxplot(geom = "errorbar",width=0.2)+
  geom_boxplot(width=0.5,aes(fill=group),outlier.shape = NA)+
  #geom_jitter(aes(color=group))+
  theme_prism(base_size = 11,
              axis_text_angle = 45,
              border = T,
              base_line_size = 0.5,
              base_fontface = "plain")+
  theme(legend.position = 'none')+
  stat_compare_means(method = "wilcox.test",label = "p.format")+
  #scale_fill_simpsons()+
  scale_fill_manual(values = c("Tumor"="#d95f02","Normal"="#1b9e77"))+
  xlab("")+
  ylab("eccDNA full gene number")
