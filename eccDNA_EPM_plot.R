library(ggsci)
library(ggprism)
library(tidyverse)
library(ggpubr)

#################################color select###################################
library(ggsci)
pal_startrek("uniform", alpha = 1)(7)
tumor_color="#C9534F"
normal_color="#415389"


#---------------------------------EPM-------------------------------------------
#eccDNA pre Million reads
get_ecc<-function(x){
  df<-read_tsv(stringr::str_glue("KidneyCancer/FinallyData/Filter/{x}_filter.bed"),col_names = F)
  mapping<-read_tsv(stringr::str_glue("KidneyCancer/FinallyData/Total.Mapping.Reads/{x}_mapping.flag.txt"),col_names = F)
  mapping<-filter(mapping,X1!="chrM")
  df1<-data.frame(count=nrow(df),totalmapping=sum(mapping$X3))
  df1$sample<-x
  df1$eccP<-df1$count/df1$totalmapping*1e+06
  df1
}
Kca_file<-read.table("KidneyCancer/sample/sample2R.tsv",header = T)

eccnumber<-do.call("bind_rows",lapply(c(Kca_file$V1),get_ecc))

df_p_value<-eccnumbers%>%
  rstatix::wilcox_test(eccP ~ group)%>%
  rstatix::add_significance(p.col = "p")%>%
  rstatix::add_xy_position()

ggplot(eccnumbers,aes(x=group,y=eccP))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_point(aes(color=group),size=2,position = "jitter")+
  theme_prism(base_size = 12,base_line_size = 0.5)+
  add_pvalue(df_p_value,
             label = "p",
             bracket.size = 0.4)+
  scale_color_manual(values = c("tumor"="#e41a1c","normal"="#377eb8"))+
  xlab("")+
  ylab("eccDNA load (EPM)")

sample=read.xlsx("KidneyCancer/sample/clinic.xlsx")

eccnumber_tum = filter(group%in%"Tumor") %>% 
  left_join(sample,by="sample")

com=list(c("Ⅰ","Ⅱ"),
         c("Ⅰ","Ⅲ"),
         c("Ⅰ","Ⅳ"))

ggplot(eccnumbers,aes(x=stage,y=eccP))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot(width=0.5,outlier.shape = NA)+
  geom_point(aes(color=group),size=2,position = "jitter")+
  theme_prism(base_size = 12,base_line_size = 0.5)+
  stat_compare_means(comparisons = com,
                     method = "wilcox.test",
                     label="p.format")+
  scale_color_manual(values = c("Ⅰ"="#868686FF","Ⅱ"="#DF8F44FF","Ⅲ"="#374E55FF","Ⅳ"="#B24745FF"))+
  xlab("")+
  ylab("eccDNA load (EPM)")
