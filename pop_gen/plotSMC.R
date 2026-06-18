setwd('/Users/maria/Desktop/Kenkel_lab/Acerv_CBASS_GWAS/haplotype_panel/pop_gen/2025/morePAN/SMC++')

df=read.csv('SMC_out/all_pops_u4e09_sub8_5runs.csv')
title="sub8_5runs"

pop_colors=list(cols=c("MEX_BZ"="#AA36A5","FLKeys"="#EE962F","DR_JAM"="#C6CF18","AR_CU"="#04A2A0","PAN"="#4E66C2"))
labels=list(new=c("MEX_BZ"="BEL/MEX","FLKeys"="DRT/FLL/FLU","DR_JAM"="JAM/DOM","AR_CU"="CUR/ARU","PAN"="PAN"))

library(ggplot2)
library(scales)
# add in group=plot_num is you have multiple sims
ggplot(df,aes(x=x,y=y,color=label,group=plot_num))+geom_path(size=1.2)+
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_y_log10(labels = label_log(digits = 2))+
  scale_color_manual(values=pop_colors$cols)+
  theme_bw(base_size=20)+xlab('years (g=5, u=4e-9)')+ylab("Ne")+
  geom_vline(xintercept = 23000,linetype='dashed')+
  ggtitle(title)

ggplot(df,aes(x=x,y=y,color=label,group=plot_num))+geom_path(size=1.2)+
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_y_log10(labels = label_log(digits = 2))+
  scale_color_manual(values=pop_colors$cols)+
  xlab('years (g=5, u=4e-9)')+ylab("Ne")+
  geom_vline(xintercept = 23000,linetype='dashed')+
  facet_wrap(~label)

# average runs together
library(dplyr)
df_sum = df %>% 
  #mutate(x_int = floor(x)) %>%
  mutate(x_int = cut(x, breaks = seq(0, 3000000, by = 1000)),
         x_mid = (as.numeric(sub("\\((.+),(.+)\\]", "\\1", x_int)) +
                    as.numeric(sub("\\((.+),(.+)\\]", "\\2", x_int))) / 2) %>%
  group_by(label,x_int,x_mid) %>%
  summarise(mean=mean(y))

ggplot(df_sum,aes(x=x_mid,y=mean,color=label))+geom_path(size=1.2)+
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_y_log10(labels = label_log(digits = 2))+
  scale_color_manual(values=pop_colors$cols)+
  xlab('years (g=5, u=4e-9)')+
  ylab("Ne")+
  geom_vline(xintercept = 23000,linetype='dashed')+
  annotate("text",x=55000, y=10^5.5, label="LGM",size=6)+
  theme_bw(base_size=20)+
  theme(legend.position = 'none')

# get "contemporary" population size
df_sum[df_sum$x_mid==500,] %>% group_by(label) %>% summarise(meanNe=mean(mean,na.rm=T))
df_sum[df_sum$x_mid==500,] %>% group_by() %>% summarise(meanNe=mean(mean,na.rm=T))
