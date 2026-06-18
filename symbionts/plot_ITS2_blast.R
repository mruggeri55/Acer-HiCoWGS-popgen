setwd('/Users/maria/Desktop/Kenkel_lab/Acerv_CBASS_GWAS/haplotype_panel/pop_gen/2025/sym')

#df=read.delim('ITS2_A3_counts_allsamps.txt',header=F,sep=" ") # blasting A3 only
#colnames(df)=c('file','ITS2','count')

df=read.delim("ITS2blast_all_samples_Nov2025.out",header=F,sep=" ") #blasting full database
colnames(df)=c('ITS2','count','file')

library(dplyr)
library(tidyr)
df_wide <- df %>% pivot_wider(names_from = ITS2, values_from = count)
df_wide[is.na(df_wide)] <- 0

# convert counts to relative abundance
df_wide$total=rowSums(df_wide[,2:ncol(df_wide)])
hist(df_wide$total)
# one sample with only 50 reads (CU_1051)
df_wide$site=factor(sapply(strsplit(df_wide$file,'_'),FUN='[',1),levels=c('CRF','Mote','DRTO','MEX','IBBZ','JAM','DR','CU21E','AR'))
library(ggplot2)
p.depth=ggplot(df_wide,aes(x=file,y=total))+
  geom_bar(stat='identity')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off")+
  theme_bw(base_size = 18)+
  facet_grid(~site,scales = "free", switch = "x", space = "free")+
  theme(panel.spacing.x = grid:::unit(0, "lines"), 
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank(),
        axis.text.x=element_blank())+
  labs(x='',y='ITS2 depth')

# export df wide for phyloseq
write.csv(df_wide,'data/sample_otu_table.csv')

rel_abund=df_wide %>%
  mutate(across(c(2:8),.fns = ~./total))

rel_long <- rel_abund[,1:8] %>% pivot_longer(!file,names_to="ITS2",values_to="rel_abund")

rel_long$pop=factor(sapply(strsplit(rel_long$file,'_'),FUN='[',1),
                    levels=c('CRF','Mote','DRTO','MEX','IBBZ','JAM','DR','CU21E','AR'))
rel_order=rel_long[order(rel_long$pop,rel_long$rel_abund),]
sample_order=rel_order$file
labs=levels(rel_order$pop)

library(forcats)
p.comp=ggplot(rel_order, aes(fill=ITS2, y=rel_abund, x=file)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  ylab('relative abundance')+
  xlab(NULL)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + theme(axis.text.x=element_blank()) +
  facet_grid(~fct_inorder(pop), scales = "free", switch = "x", space = "free") +
  theme(panel.spacing.x = grid:::unit(0, "lines"), 
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank())

library(cowplot)
X11()
plot_grid(p.depth+theme(plot.margin = unit(c(0.2, 0.2, -0.8, 0.2), "cm")), p.comp, ncol = 1, align = "v")
dev.print(pdf,"ITS2_depth_and_comp_Nov2025.pdf")
dev.off()
