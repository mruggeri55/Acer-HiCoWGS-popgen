setwd('/Users/maria/Desktop/Kenkel_lab/Acerv_CBASS_GWAS/haplotype_panel/pop_gen/2025/sym')
# adapted from https://github.com/cbreusing/Provannid_host-symbiont_popgen/blob/main/Symbiont/RedundancyAnalysis.r

dat <- read.table("data/sym.freebayes.GT.GT.FORMAT",header=T,na.strings = ".")
rownames(dat) <- paste(dat$CHROM,dat$POS,sep=':')
data <- t(dat[3:37])
# fill in missing genotypes with most frequently observed
newdat <- apply(data, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
newdat[1:10,1:10]

meta=as.data.frame(rownames(newdat))
meta$site=sapply(strsplit(meta$`rownames(newdat)`,'_'),FUN='[',1)
region_df=as.data.frame(cbind('site'=c('CRF','Mote','DRTO','MEX','IBBZ','DR','JAM','AR','CU21E'),
                     'region'=c('FLkeys','FLkeys','FLkeys','MEX/BZ','MEX/BZ','DR/JAM','DR/JAM','AR/CU','AR/CU')))
library(plyr)
meta2=join(meta,region_df)
meta=meta2

# checking for missingness
nNA=as.data.frame(rowSums(is.na(data)))
nNA$sample=meta$`rownames(newdat)`
ggplot(nNA,aes(x=sample,y=rowSums(is.na(data))))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# may want to remove some samples with high missingness
nNA$miss_rate=nNA$`rowSums(is.na(data))`/21739
nNA[nNA$miss_rate>0.2,]
# 3 samples with more that 20% missing
badsamps=c('JAM_11','CU21E_1052','IBBZ_13837')

# RDA
library(vegan)
rda <- capscale(newdat ~ site, data = meta, na.action = na.exclude, add = "cailliez")
RsquareAdj(rda)
vif.cca(rda)

signif.full <- anova.cca(rda, parallel=getOption("mc.cores"))
signif.full
signif.axis <- anova.cca(rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
signif.margin <- anova.cca(rda, by="margin", parallel=getOption("mc.cores"))
signif.margin

meta$site=factor(meta$site, levels=c('CRF','Mote','DRTO','MEX','IBBZ','DR','JAM','AR','CU21E'))
color=c('CRF'="#FF7F00", 'Mote'="#FDBF6F",'DRTO'="#FFFF99",'MEX'="#E31A1C", 'IBBZ'="#FB9A99", 'DR'="#33A02C", 'JAM'="#B2DF8A", 'AR'="#1F78B4",'CU21E'="#A6CEE3")

plot(rda, type="n", scaling=3)
points(rda, pch=20, cex=2, scaling=3,col=color[meta$site])
text(rda, scaling=3, display="bp", col="black", cex=1)

# by region
rda_reg <- capscale(newdat ~ region, data = meta, na.action = na.exclude, add = "cailliez")
reg_col=c('FLkeys'="#FF7F00",'MEX/BZ'="#E31A1C",'DR/JAM'="#33A02C",'AR/CU'="#1F78B4")

plot(rda_reg, type="n", scaling=3)
points(rda_reg, pch=20, cex=2, scaling=3,col=reg_col[meta$region])
text(rda_reg, scaling=3, display="bp", col="black", cex=1)


######## try a PCA
library(tidyverse)
pca <- prcomp(newdat)
pc_scores <- as.data.frame(pca$x)
pc_scores$site=meta$site
pc_scores$region=meta$region
pc_scores$region=factor(pc_scores$region,levels=c('FLkeys','MEX/BZ','DR/JAM','AR/CU'))

pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size=3,aes(fill=site),pch=22)+
  scale_fill_manual(values=color)+
  stat_ellipse(aes(color=region))+
  scale_color_manual(values=c("#FF7F00","#E31A1C","#33A02C","#1F78B4"))+
  ggtitle('major haplotype')+
  xlab('PC1')+
  ylab('PC2')

# try excluding samples with high missingness
pca2 <- prcomp(newdat[!rownames(newdat) %in% badsamps,])
summary(pca2)
pc_scores2 <- as.data.frame(pca2$x)
pc_scores2$site=meta$site[!meta$`rownames(newdat)` %in% badsamps]
pc_scores2$region=meta$region[!meta$`rownames(newdat)` %in% badsamps]
pc_scores2$region=factor(pc_scores2$region,levels=c('FLkeys','MEX/BZ','DR/JAM','AR/CU'))


pc_scores2 %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size=3,aes(fill=site),pch=22)+
  scale_fill_manual(values=color)+
  stat_ellipse(aes(color=region))+
  scale_color_manual(values=c("#FF7F00","#E31A1C","#33A02C","#1F78B4"))+
  ggtitle('major haplotype = dominant strain')+
  xlab('PC1 (3.79%)')+
  ylab('PC2 (3.73%)')

# exclude loci with any missing
nomiss=na.omit(t(data))
# 1361 loci remaining
write.csv(nomiss,'data/sym.freebayes.GT.nomiss.csv')
write.csv(meta,'data/sym.freebayes.GT.nomiss.meta.csv')

pca3 <- prcomp(t(nomiss))
summary(pca3)
pc_scores3 <- as.data.frame(pca3$x)
pc_scores3$site=meta$site
pc_scores3$region=meta$region
pc_scores3$region=factor(pc_scores3$region,levels=c('FLkeys','MEX/BZ','DR/JAM','AR/CU'))

pc_scores3 %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size=3,aes(fill=site),pch=22)+
  scale_fill_manual(values=color)+
  stat_ellipse(aes(color=region))+
  scale_color_manual(values=c("#FF7F00","#E31A1C","#33A02C","#1F78B4"))+
  ggtitle('major haplotype = dominant strain (1,361 loci)')+
  xlab('PC1 (5%)')+
  ylab('PC2 (4.5%)')

