library(ggprism)
library(tidyverse)
library(ggpubr)

#----------------------------------length---------------------------------------
#length-------------------------------------------------------------------------
#input data
length_calculate <- function(x) {
  df <- read.table(stringr::str_glue("bed/{x}_filter.bed"),header = T)
  df <- df %>%
    dplyr::select(chrom,length)
  df$sample<-x
  df
}

s_file = read_tsv("sample2R.tsv")

length<-do.call("bind_rows",lapply(s_file$sample,length_calculate))

length = left_join(length,s_file,by="sample")

length = read_tsv("CircleSeq/results.count/length.tsv")

length.short<-length %>% filter(length <=2000) 

length.long<-length %>% filter(length >2000 & length<=10000) 

#-------------------short length eccDNA-------------------------------------

ggplot(length.short,aes(x=length,
                     color=group,
                     fill=group))+
  geom_density(alpha = 0.3)+
  #facet_wrap(.~ group)+
  scale_color_manual(values =c("#006699","#cc0033"))+
  scale_fill_manual(values =c("#006699","#cc0033"))+
  theme_prism(base_size = 11,
              axis_text_angle = 45,
              border = T,
              base_line_size = 0.5,
              base_fontface = "plain")+
  #theme(legend.position = 'none')+
  xlab("Length (bp)")+
  labs(title="Length density of eccDNA in KCA")

#-------------------Long length eccDNA-------------------------------------
ggplot(length.long,aes(x=length,
                       #fill=group,
                      color=group))+
  geom_density(alpha = 0.3)+
  #facet_wrap(.~ group)+
  scale_color_manual(values =c("#006699","#cc0033"))+
  #scale_fill_manual(values =c("#006699","#cc0033"))+
  theme_prism(base_size = 11,
              axis_text_angle = 45,
              border = T,
              base_line_size = 0.5,
              base_fontface = "plain")+
  #theme(legend.position = 'none')+
  xlab("Length (bp)")+
  labs(title="Length density of eccDNA in KCA")
