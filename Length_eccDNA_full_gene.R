library(ggplot2)
library(ggprism)

get_full_length<-function(x){
  df<-read_tsv(stringr::str_glue("fullGene.bed/{x}_coding_full.bed"),col_names = F)
  df1 = df %>% mutate(ecc_length = X3-X2) %>% 
    select(X16,ecc_length)
  df1$sample=x
  df1
}

clinc = read_tsv("sample2R.tsv")

clinc=clinc[-which(clinc$sample%in%epm.s),]

eccLength<-do.call("bind_rows",lapply(clinc$sample,get_full_length))

df = eccLength %>% left_join(clinc)

length.long<-df %>% filter(ecc_length > 10000 ) 

ggplot(length.long,aes(x=ecc_length,color=group,fill=group))+
    #geom_line(aes(x=length,y=Frequency))+
    geom_density(alpha = 0.3)+
    #scale_x_continuous(breaks = seq(0,1000, by = 5))+
    #geom_histogram(aes(color = group, fill = group),
    #               stat = "bin",bins=90,alpha = 0.1, position = "identity")+
    #facet_wrap(~ group)+
    scale_color_manual(values =c("Tumor"="#e41a1c","Normal"="#377eb8"))+
    scale_fill_manual(values =c("Tumor"="#e41a1c","Normal"="#377eb8"),)+
    theme_prism(base_size = 11,
                axis_text_angle = 45,
                border = T,
                base_line_size = 0.5,
                base_fontface = "plain")
