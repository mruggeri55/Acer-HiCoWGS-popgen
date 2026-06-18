setwd('/Users/maria/Desktop/Kenkel_lab/Acerv_CBASS_GWAS/haplotype_panel/pop_gen/2025/morePAN/private_alleles')

df=read.csv('private_alleles_AAD10x.csv')

library(reshape2)
df_long=melt(df,value.name='pa_count',variable.name = 'af')
df_long$af=as.numeric(gsub('AF','0',df_long$af))

library(ggplot2)
pop_colors=list(cols=c("MEX_BZ"="#AA36A5","Flkeys"="#EE962F","DR_JAM"="#C6CF18","AR_CU"="#04A2A0","PAN"="#4E66C2"))
labels=list(new=c("MEX_BZ"="BEL/MEX","Flkeys"="DRT/FLL/FLU","DR_JAM"="JAM/DOM","AR_CU"="CUR/ARU","PAN"="PAN"))

df_long$POP=factor(df_long$POP, levels=names(pop_colors$cols))

ggplot(df_long,aes(x=af,y=pa_count,group=POP,color=POP))+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_color_manual(values=pop_colors$cols)+
  theme_bw(base_size=16)+
  labs(x='allele frequency',y='# private alleles')
ggplot(df_long,aes(x=af,y=pa_count,group=POP,color=POP))+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_color_manual(values=pop_colors$cols)+
  theme_bw(base_size=16)+
  labs(x='allele frequency',y='log10(private alleles)')+
  scale_y_log10()

ggplot(df_long,aes(x=af,y=pa_count,group=POP,fill=POP))+
  geom_col()+
  scale_fill_manual(values=pop_colors$cols)+
  theme_bw(base_size=16)+
  labs(x='allele frequency',y='log10(private alleles)')+
  scale_y_log10()

#using pseudo log sclae so can visualize AR_CU fixed private allele
# also remove af of 0.1 just so consistent intervals
library(scales)
ggplot(df_long[!df_long$af==0.1,],aes(x=af,y=pa_count,group=POP,fill=POP))+
  geom_col()+
  scale_fill_manual(values=pop_colors$cols,labels = labels$new)+
  theme_bw(base_size=14)+
  labs(x='allele frequency',y='log10(private alleles)')+
  #scale_y_log10()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10),
                     breaks = c(0,1e5,1e10,1e20,1e50))+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1))

# make elluvial plot
library(ggalluvial)
df_long$af_bin <- cut(df_long$af,
                      breaks = c(0, 0.1, 0.25, 0.5, 0.75, 1),
                      include.lowest = TRUE)
df_sum <- df_long %>%
  filter(af != 0.1) %>%
  group_by(POP, af_bin) %>%
  summarise(pa_count = sum(pa_count), .groups = "drop")

ggplot(df_sum,
       aes(x = af_bin,
           stratum = POP,
           alluvium = POP,
           y = pa_count,
           fill = POP)) +
  geom_flow(alpha = 0.6) +
  geom_stratum(width = 0.5,linewidth = 0.2) +
  scale_fill_manual(values = pop_colors$cols) +
  theme_bw(base_size = 16) +
  labs(x = "Allele frequency bin",
       y = "Private alleles")+
  scale_y_continuous(trans = pseudo_log_trans(base = 10),
                     breaks = c(0,1e5,1e10,1e20,1e50))

