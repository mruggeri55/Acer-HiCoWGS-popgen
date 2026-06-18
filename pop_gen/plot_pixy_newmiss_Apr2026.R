setwd('/Users/maria/Desktop/Kenkel_lab/Acerv_CBASS_GWAS/haplotype_panel/pop_gen/2025/morePAN/pixy')

pi=read.delim('data/100kb_sub8_pi.txt')
TajD=read.delim('data/100kb_sub8_tajima_d.txt')
theta=read.delim('data/100kb_sub8_watterson_theta.txt')

library(ggplot2)
ggplot(pi,aes(x=pop,y=avg_pi,fill=pop))+
  geom_boxplot()
ggplot(theta,aes(x=pop,y=raw_watterson_theta,fill=pop))+
  geom_boxplot()
ggplot(TajD,aes(x=pop,y=tajima_d,fill=pop))+
  geom_boxplot()

# exluding windows with >50% missing
# also adding color scheme
pop_colors=list(cols=c("MEX_BZ"="#AA36A5","FLkeys"="#EE962F","DR_JAM"="#C6CF18","AR_CU"="#04A2A0","PAN"="#4E66C2"))
labels=list(new=c("MEX_BZ"="BEL/MEX","FLkeys"="DRT/FLL/FLU","DR_JAM"="JAM/DOM","AR_CU"="CUR/ARU","PAN"="PAN"))

pi$pop=factor(pi$pop, levels=names(pop_colors$cols))
TajD$pop=factor(TajD$pop, levels=names(pop_colors$cols))

p.pi=ggplot(pi[pi$no_sites>=0.5*100000,],aes(x=pop,y=avg_pi,fill=pop))+
  geom_boxplot()+
  scale_fill_manual(values=pop_colors$cols)+
  scale_x_discrete(labels = labels$new)+
  theme_bw(base_size=16)+
  xlab('')+
  ylab(expression(pi))+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, hjust = 1))
p.d=ggplot(TajD[TajD$no_sites>=0.5*100000,],aes(x=pop,y=tajima_d,fill=pop))+
  geom_boxplot()+
  scale_fill_manual(values=pop_colors$cols)+
  scale_x_discrete(labels = labels$new)+
  theme_bw(base_size=16)+
  xlab('')+
  ylab('Tajimas D')+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45, hjust = 1))

library(grid)
grid.draw(rbind(ggplotGrob(p.pi+theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      plot.margin = unit(c(0.2, 0.2, -0.8, 0.2), "cm"))),
                ggplotGrob(p.d+theme(plot.margin = unit(c(0.2, 0.2, -0.5, 0.2), "cm"))),
                size='first'))

# ##### plotting by chromosome
# ##### pi
# inp<-pi
# 
# # Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
# #   e.g., chr1, chr2, chr22, chrX
# chroms <- unique(inp$chromosome[!grepl("CAX",inp$chromosome)]) #remove unplaces scaffolds (CAX)
# chrOrder <- sort(chroms)
# inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)

# # Plot pi for each population found in the input file
# # Saves a copy of each plot in the working directory
# if("avg_pi" %in% colnames(inp)){
#   pops <- unique(inp$pop)
#   for (p in pops){
#     thisPop <- subset(inp, pop == p)
#     # Plot stats along all chromosomes:
#     popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color=chrOrder)) +
#       geom_point()+
#       facet_grid(. ~ chrOrder)+
#       labs(title=paste("Pi for population", p))+
#       labs(x="Position of window start", y="Pi")+
#       scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
#       theme_classic()+
#       theme(legend.position = "none")
#     ggsave(paste("piplot_", p,"100kb.png", sep=""), plot = popPlot, device = "png", dpi = 300,width=8,height=4)
#   }
# } else {
#   print("Pi not found in this file")
# }
# 
# ######## TajD
# inp <- TajD
# chroms <- unique(inp$chromosome[!grepl("CAX",inp$chromosome)]) #remove unplaces scaffolds (CAX)
# chrOrder <- sort(chroms)
# inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)
# 
# if("tajima_d" %in% colnames(inp)){
#   pops <- unique(inp$pop)
#   for (p in pops){
#     thisPop <- subset(inp, pop == p)
#     # Plot stats along all chromosomes:
#     popPlot <- ggplot(thisPop, aes(window_pos_1, tajima_d, color=chrOrder)) +
#       geom_point()+
#       facet_grid(. ~ chrOrder)+
#       labs(title=paste("Tajima's D for population", p))+
#       labs(x="Position of window start", y="tajima_d")+
#       scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
#       theme_classic()+
#       theme(legend.position = "none")
#     ggsave(paste("TajDplot_", p,"100kb.png", sep=""), plot = popPlot, device = "png", dpi = 300,width=8,height=4)
#   }
# } else {
#   print("TajD not found in this file")
# }
# 

#### summary stats ####
library(dplyr)
#### how many sites with significant tajD (>2 or <-2) 
TajD[TajD$no_sites>=0.5*100000,] %>% group_by(pop) %>% summarise(common_count=sum(tajima_d >= 2),rare_count=sum(tajima_d <= -2),n_win=n())
# get weighted average for windows (based on number of sites in each)
TajD %>% group_by(pop) %>% 
  summarise(global_D=sum(tajima_d*no_sites, na.rm=T)/sum(no_sites))
# pop    global_D
# 1 MEX_BZ   -0.256
# 2 FLkeys   -0.194
# 3 DR_JAM   -0.338
# 4 AR_CU    -0.218
# 5 PAN      -0.270
TajD[TajD$no_sites>=0.5*100000,] %>% group_by(pop) %>% 
  summarise(global_D=sum(tajima_d*no_sites, na.rm=T)/sum(no_sites))
# 1 MEX_BZ   -0.267
# 2 FLkeys   -0.202
# 3 DR_JAM   -0.353
# 4 AR_CU    -0.226
# 5 PAN      -0.277

# linear model to test for diffs in TajD
library(lmerTest)
library(emmeans)
TajD$window=paste(TajD$chromosome,TajD$window_pos_1,TajD$window_pos_2)
Tmod <- lmer(tajima_d ~ pop + (1 | chromosome/window), data = TajD[TajD$no_sites>=0.5*100000,])
summary(Tmod)
anova(Tmod)
emmeans(Tmod, pairwise ~ pop)

#### avg pi across pops (note cannot just average windows because of missingness)
pi %>% group_by(pop) %>% 
  summarise(global_pi=sum(count_diffs,na.rm = TRUE)/sum(count_comparisons,na.rm = TRUE))
# pop    global_pi
# 1 MEX_BZ   0.00240
# 2 FLkeys   0.00224
# 3 DR_JAM   0.00224
# 4 AR_CU    0.00219
# 5 PAN      0.00226
pi[pi$no_sites>=0.5*100000,] %>% group_by(pop) %>% 
  summarise(global_pi=sum(count_diffs,na.rm = TRUE)/sum(count_comparisons,na.rm = TRUE))
# 1 MEX_BZ   0.00222
# 2 FLkeys   0.00208
# 3 DR_JAM   0.00206
# 4 AR_CU    0.00202
# 5 PAN      0.00209

pi$window=paste(pi$chromosome,pi$window_pos_1,pi$window_pos_2)
mod <- lmer(avg_pi ~ pop + (1 | chromosome/window), data = pi[pi$no_sites>=0.5*100000,])
summary(mod)
anova(mod)
emmeans(mod, pairwise ~ pop)

#### extracting outliers ####
pi=read.delim('data/10kb_global_pi.txt')
#pi=read.delim('data/1kb_global_pi.txt')
#pi=read.delim('data/100kb_global_pi.txt')
window=10000
up.quant=0.9999
low.quant=0.0001

box=boxplot(pi$avg_pi[pi$no_sites>=0.5*window])
outliers=box$out
out_sub=pi[pi$avg_pi %in% outliers,]
# restrict to 0.1% or 0.01% outlier
upper_bound <- quantile(pi$avg_pi[pi$no_sites>=0.5*window], up.quant, na.rm=T)
lower_bound <- quantile(pi$avg_pi[pi$no_sites>=0.5*window], low.quant, na.rm=T)
pi_outliers <- pi[pi$no_sites>=0.5*window & !is.na(pi$avg_pi) & (pi$avg_pi <=lower_bound | pi$avg_pi >= upper_bound),]
write.csv(pi_outliers,paste('results/pi_outliers_',window/1000,'kb_',low.quant,'thresh','_miss5kb','.csv',sep=''))

library(dplyr)
pi_color <- pi %>%
  filter(!grepl("CAX|OZ035980.1",chromosome)) %>%
  filter(no_sites>=0.5*window) %>%
  mutate(
    color_group = ifelse(avg_pi <= lower_bound | avg_pi >= upper_bound, "Outside Threshold", "Within Threshold"),
  )

ggplot(pi_color, aes(window_pos_1, avg_pi,color=color_group)) +
  geom_point(size=0.8)+
  labs(x="Position of window start", y="Pi")+
  scale_color_manual(values=c("red","black"))+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(~chromosome)+
  ggtitle(paste(paste(window,'kb',sep=''),paste(low.quant*100,'%',sep=''),sep=" "))

ggplot(pi_color, aes(chromosome, avg_pi))+
  geom_boxplot()+
  geom_hline(yintercept = upper_bound,linetype='dashed')+
  geom_hline(yintercept = lower_bound,linetype='dashed')+
  labs(y="Pi")+
  theme_bw()+
  ggtitle(paste(paste(window/1000,'kb',sep=''),paste(low.quant*100,'%',sep=''),sep=" "))

p.list=list()
for (chr in unique(pi_color$chromosome)){
p=ggplot(pi_color[pi_color$chromosome==chr,], aes(window_pos_1, avg_pi,color=color_group)) +
  geom_point(size=0.8)+
  labs(x="Position of window start", y="Pi")+
  scale_color_manual(values=c("red","black"))+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle(paste(chr))
p.list[[chr]]<-p
}
p.list[["OZ035971.1"]]
p.list[["OZ035972.1"]]
p.list[["OZ035973.1"]]
p.list[["OZ035975.1"]]

######## repeat with TajD
#TajD=read.delim('data/1kb_global_tajima_d.txt')
TajD=read.delim('data/10kb_global_tajima_d.txt')
#TajD=read.delim('data/100kb_global_tajima_d.txt') 
window=10000
up.quant=0.9999
low.quant=0.0001

# restrict to 0.1% or 0.01% outlier
upper_bound <- quantile(TajD$tajima_d[TajD$no_sites>=0.5*window], up.quant, na.rm=T)
lower_bound <- quantile(TajD$tajima_d[TajD$no_sites>=0.5*window], low.quant, na.rm=T)
TajD_outliers <- TajD[TajD$no_sites>=0.5*window & !is.na(TajD$tajima_d) & (TajD$tajima_d <= lower_bound | TajD$tajima_d >= upper_bound),]
write.csv(TajD_outliers,paste('results/TajD_outliers_',window/1000,'kb_',low.quant,'thresh','_miss5kb','.csv',sep=''))

TajD_color <- TajD %>%
  filter(!grepl("CAX|OZ035980.1",chromosome)) %>%
  filter(no_sites>=0.5*window) %>%
  mutate(
    color_group = ifelse(tajima_d <= lower_bound | tajima_d >= upper_bound, "Outside Threshold", "Within Threshold"),
    color_group2 = case_when(tajima_d > 2  ~ "above", tajima_d < -2 ~ "below",TRUE          ~ "within")
  )

ggplot(TajD_color, aes(window_pos_1, tajima_d,color=color_group)) +
  geom_point(size=0.6)+
  labs(x="Position of window start", y="Tajima's D")+
  scale_color_manual(values=c("red","black"))+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(~chromosome)+
  ylim(lower_bound-1,upper_bound+1)+
  ggtitle(paste(paste(window,'kb',sep=''),paste(low.quant*100,'%',sep=''),sep=" "))

ggplot(TajD_color, aes(x=chromosome, y=tajima_d))+
  geom_boxplot(size=0.6)+
  geom_hline(yintercept = upper_bound,linetype='dashed')+
  geom_hline(yintercept = lower_bound,linetype='dashed')+
  labs(y="Tajima's D")+
  theme_bw()+
  ggtitle(paste(paste(window,'kb',sep=''),paste(low.quant*100,'%',sep=''),sep=" "))

# the following needs to be modified
p.list=list()
for (chr in unique(pi_color$chromosome)){
  p=ggplot(pi_color[pi_color$chromosome==chr,], aes(window_pos_1, avg_pi,color=color_group)) +
    geom_point(size=0.8)+
    labs(x="Position of window start", y="Pi")+
    scale_color_manual(values=c("red","black"))+
    theme_bw()+
    theme(legend.position = "none")+
    ggtitle(paste(chr))
  p.list[[chr]]<-p
}

# try plotting pi and tajimas D together for a signle chromosome
P.chr1=ggplot(pi_color[pi_color$chromosome=='OZ035972.1',], aes(x=window_pos_1, y=avg_pi,color=color_group)) +
  geom_point(size=0.8)+
  labs(x="Position of window start", y="Pi")+
  scale_color_manual(values=c("Outside Threshold"="red","Within Threshold"="black"))+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(~chromosome)
  
T.chr1=ggplot(TajD_color[TajD_color$chromosome=='OZ035972.1',], aes(x=window_pos_1, y=tajima_d,color=color_group)) +
  geom_point(size=0.6)+
  labs(x="Position of window start", y="Tajima's D")+
  scale_color_manual(values=c("Outside Threshold"="red","Within Threshold"="black"))+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(~chromosome)

library(gridExtra)
grid.arrange(P.chr1+theme(axis.title.y = element_blank(), axis.text.y = element_blank()),
             T.chr1+theme(axis.title.y = element_blank(), axis.text.y = element_blank()),ncol=1)

# plot TajD with pi colors
TajD_pi_color=TajD_color
TajD_pi_color$pi_color=pi_color$color_group
TajD_pi_color$avg_pi=pi_color$avg_pi
 
ggplot(TajD_pi_color, aes(x = window_pos_1, y = tajima_d)) +
  geom_point(
    data = TajD_pi_color[TajD_pi_color$pi_color == "Within Threshold", ],
    color = "black",
    size = 0.1
  ) +
  geom_point(
    data = TajD_pi_color[TajD_pi_color$pi_color == "Outside Threshold", ],
    color = "red",
    size = 1
  ) +
  labs(x = "Position of window start", y = "TajD") +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~chromosome) +
  ggtitle("TajD colored by pi outliers",paste(paste(window/1000,'kb',sep=''),paste(low.quant*100,'%',sep=''),sep=" "))

##### STOPPING HERE ######

# zooming in on pi outlier chr 1 identified with 10kb windows but will plot 1kb values
pi_color$out_window=pi_color$window_pos_1 >= 27680399 & pi_color$window_pos_1 <= 27690398
ggplot(pi_color[pi_color$chromosome=='OZ035966.1',], aes(x=window_pos_1, y=avg_pi,color=out_window))+
  geom_point(size=0.8)+
  labs(x="Position of window start", y="Pi")+
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(~chromosome)+
  xlim(c(27500000,27900000))+
  annotate("text", x = 27687164.5, y = 0.07, label = "OAS3")

TajD_color$out_window=TajD_color$window_pos_1 >= 27680399 & TajD_color$window_pos_1 <= 27690398
ggplot(TajD_color[TajD_color$chromosome=='OZ035966.1',], aes(x=window_pos_1, y=tajima_d,color=out_window))+
  geom_point(size=0.8)+
  labs(x="Position of window start", y="tajimas D")+
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(~chromosome)+
  xlim(c(27500000,27900000))+
  annotate("text", x = 27687164.5, y = 5.5, label = "OAS3")

#### manhattan plots
library(qqman)
TajD_pi_color$windowID=paste('window',rownames(TajD_pi_color))
# rename chromosomes
old_chr <- c(
  "OZ035966.1","OZ035967.1","OZ035968.1","OZ035969.1",
  "OZ035970.1","OZ035971.1","OZ035972.1","OZ035973.1",
  "OZ035974.1","OZ035975.1","OZ035976.1","OZ035977.1",
  "OZ035978.1","OZ035979.1"
)
new_chr <- as.character(1:14)
# apply to your dataframe column
TajD_pi_color$chr <- new_chr[match(TajD_pi_color$chr, old_chr)]
TajD_pi_color$chr <- as.numeric(TajD_pi_color$chr)

pi.highlight=TajD_pi_color$windowID[TajD_pi_color$pi_color=="Outside Threshold"]
manhattan(TajD_pi_color, chr="chr", bp="window_pos_1", snp="windowID",
          highlight=pi.highlight, p="avg_pi", logp=FALSE,
          ylab="Nucleotide Diversity (pi)",ylim = c(0, 0.055),
          genomewideline=upper_bound)
manhattan(TajD_pi_color, chr="chr", bp="window_pos_1", snp="windowID",
          highlight=pi.highlight, p="tajima_d", logp=FALSE,
          ylab="Tajima's D",ylim = c(-5,5),
          genomewideline=F, suggestiveline = F)

#### trying heatmap for OAS3
library(vcfR)
vcf <- read.vcfR("data/OAS3.vcf.gz")
gt <- extract.gt(vcf)   # matrix: variants × samples
# Replace phased | with /, just in case
gt <- gsub("\\|", "/", gt)
# Convert to numeric: 0/0=0, 0/1 or 1/0=1, 1/1=2
gt_num <- matrix(NA, nrow = nrow(gt), ncol = ncol(gt))
gt_num[gt == "0/0"] <- 0
gt_num[gt %in% c("0/1", "1/0")] <- 1
gt_num[gt == "1/1"] <- 2
# Optional: missing genotypes "./." remain NA
gt_num[gt == "./."] <- NA

mat <- t(gt_num)
rownames(mat) <- colnames(gt)          # sample names
pos <- getPOS(vcf)
colnames(mat) <- pos  

library(pheatmap)
pheatmap(
  mat,
  cluster_rows = T,   
  cluster_cols = FALSE,    # keep variant order
  color = c("white", "lightblue", "blue"), # white is no alt, lightblue is het, blue is 2 alt
  na_col = "grey",
  main = "Genotype Heatmap"
)

# filtering invariant sites
# Compute variance for each column (variant)
variant_var <- apply(mat, 2, function(x) var(x, na.rm = TRUE))

# Keep only columns where variance > 0 (i.e., polymorphic sites)
mat_var <- mat[, variant_var > 0]

pheatmap(
  mat_var,
  cluster_rows = T,   
  cluster_cols = FALSE,    # keep variant order
  color = c("white", "lightblue", "blue"), 
  na_col = "grey",
  main = "Genotype Heatmap",
  show_colnames = FALSE,
  fontsize_row = 6
)

######## plot sacsin gene
# note unnannotated in Acer genome
# based on blast to Amil sacsin, closest region appears to be OZ035970.1 14408881 to 14421074
ggplot(pi_color[pi_color$chromosome=='OZ035970.1',], aes(x=window_pos_1, y=avg_pi)) +
  geom_point(size=0.8)+
  labs(x="Position of window start", y="Pi")+
  theme_bw()+
  theme(legend.position = "none")+
  xlim(c(14300000,14600000))
# so actually kind of missing from our data
# what if we plot without filtering for window coverage
ggplot(pi[pi$chromosome=='OZ035970.1',], aes(x=window_pos_1, y=avg_pi)) +
  geom_point(size=0.8)+
  labs(x="Position of window start", y="Pi")+
  theme_bw()+
  theme(legend.position = "none")+
  xlim(c(14300000,14600000))
# then shows up but regardless is low <0.01 pi
