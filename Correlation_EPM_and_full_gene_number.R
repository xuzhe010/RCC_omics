#full gene number 与emp相关性

epm.data = read.xlsx("Kidney_circle_epms.xlsx")
epm.data = epm.data %>% left_join(df) %>% 
  drop_na(num)

library(ggplot2)
library(ggpubr)
library(ggprism)
library(ggsci)

ggplot(epm.data,aes(x=EPM,y=num))+
  geom_point(size=2) +   # draw points
  geom_smooth(method="lm", se=F,color="#8c510a",size=1) +
  stat_cor(data=epm.data1, method = "pearson")+
  theme_prism(base_size = 11,
              axis_text_angle = 45,
              border = T,
              base_line_size = 0.5,
              base_fontface = "plain")+
  labs(#subtitle="CTTN Vs T.cell", 
    y="Full gene number", 
    x="EPM", 
    title="Correlation")
