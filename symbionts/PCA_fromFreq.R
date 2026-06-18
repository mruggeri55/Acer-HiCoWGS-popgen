setwd('/Users/maria/Desktop/Kenkel_lab/Acerv_CBASS_GWAS/haplotype_panel/pop_gen/2025/sym')

Freq=read.table('data/sym.freebayes.RefAltCatFreq.txt',na.strings = "./.")
rownames(Freq)=paste(Freq$V1,Freq$V2,Freq$V3,sep=':')
sample_order=read.table('data/SampleOrderFreq.txt')
sample_order$site=sapply(strsplit(sample_order$V1,'_'),FUN='[',1)
region_df=as.data.frame(cbind('site'=c('CRF','Mote','DRTO','MEX','IBBZ','DR','JAM','AR','CU21E'),
                              'region'=c('FLkeys','FLkeys','FLkeys','MEX/BZ','MEX/BZ','DR/JAM','DR/JAM','AR/CU','AR/CU')))
library(plyr)
meta=join(sample_order,region_df)
sample_order=meta

sub=Freq[,c(4:40)]
# convert to decimals so numeric
sub[] <- apply(sub, c(1, 2), function(x) if(x %in% NA){NA} else {eval(parse(text = x))})

nNA=as.data.frame(colSums(is.na(sub)))
nNA$sample=sample_order$V1
ggplot(nNA,aes(x=sample,y=colSums(is.na(sub))))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
badsamps=c('JAM_11','CU21E_1052','IBBZ_13837')

library(tidyverse)
pca_mat <- sub %>% as.matrix() %>% t()
pca <- prcomp(pca_mat)
# error due to NAs in matrix

# try distance matrix to handle NAs
library(stats)
mat <- sub %>% as.matrix() %>% t()
rownames(mat)=sample_order$V1

dist <- as.matrix(dist(mat)) #euclidian distance
fit <- cmdscale(dist,eig=TRUE, k=2)
library(ggplot2)
fit_scores=as.data.frame(fit$points)
fit_scores$site=factor(sample_order$site,levels=c('CRF','Mote','DRTO','MEX','IBBZ','DR','JAM','AR','CU21E'))
ggplot(fit_scores,aes(x=V1, y=V2,color=site))+
         geom_point(size=2)

# bray curtis dissim
library(vegan)
bray_curtis_dist <- vegdist(mat, method = "bray",na.rm=T)
library(ecodist)
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
bray_curtis_pcoa_df$site=factor(sample_order$site,levels=c('CRF','Mote','DRTO','MEX','IBBZ','DR','JAM','AR','CU21E'))
ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2,color=site)) +
  geom_point(size=2) +
  labs(x = "PC1",
       y = "PC2", 
       title = "Bray-Curtis PCoA") +
  theme(title = element_text(size = 10))


#### try setting NAs to mean frequency ####
library(zoo)
new=na.aggregate(t(sub)) # fills NA in with means
rownames(new)=sample_order$V1
pca <- prcomp(new[!rownames(new) %in% badsamps,])

pc_scores <- as.data.frame(pca$x)
pc_scores$site=sample_order$site[!sample_order$V1 %in% badsamps]
#pc_scores$region=meta$region
#pc_scores$region=factor(pc_scores$region,levels=c('FLkeys','MEX/BZ','DR/JAM','AR/CU'))

sample_order$site=factor(sample_order$site, levels=c('CRF','Mote','DRTO','MEX','IBBZ','DR','JAM','AR','CU21E'))
color=c('CRF'="#FF7F00", 'Mote'="#FDBF6F",'DRTO'="#FFFF99",'MEX'="#E31A1C", 'IBBZ'="#FB9A99", 'DR'="#33A02C", 'JAM'="#B2DF8A", 'AR'="#1F78B4",'CU21E'="#A6CEE3")
pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size=3,aes(fill=site),pch=22)+
  scale_fill_manual(values=color)+
  #stat_ellipse(aes(color=region))+
  #scale_color_manual(values=c("#FF7F00","#E31A1C","#33A02C","#1F78B4"))+
  ggtitle('allele frequencies = strain abundances')+
  xlab('PC1')+
  ylab('PC2')

### let's remove loci with NAs (rows) and see how many are left
subNoNa=t(na.omit(sub))
write.csv(subNoNa,'data/sym.freebayes.allele.freqs.nomiss.csv')
# 2005 loci left
rownames(subNoNa)=sample_order$V1
pca <- prcomp(subNoNa)
summary(pca)
pc_scores <- as.data.frame(pca$x)
pc_scores$site=sample_order$site
pc_scores$region=sample_order$region
pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size=3,aes(fill=site),pch=22)+
  scale_fill_manual(values=color)+
  stat_ellipse(aes(color=region))+
  scale_color_manual(values=c("#FF7F00","#E31A1C","#33A02C","#1F78B4"))+
  ggtitle('allele frequencies = strain abundances (2,005 alleles)')+
  xlab('PC1 (5.83%)')+
  ylab('PC2 (5.55%)')

### remove bad samples first then NAs
colnames(sub)=sample_order$V1
NoBads=sub[!colnames(sub) %in% badsamps]
noBadsNoNa=na.omit(NoBads)
# 6964 loci left
noBadsNoNaMat=t(noBadsNoNa)
pca <- prcomp(noBadsNoNaMat)
pc_scores <- as.data.frame(pca$x)
pc_scores$site=sample_order$site[!sample_order$V1 %in% badsamps]
pc_scores %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size=3,aes(fill=site),pch=22)+
  scale_fill_manual(values=color)+
  #stat_ellipse(aes(color=region))+
  #scale_color_manual(values=c("#FF7F00","#E31A1C","#33A02C","#1F78B4"))+
  ggtitle('allele frequencies = strain abundances')+
  xlab('PC1')+
  ylab('PC2')
