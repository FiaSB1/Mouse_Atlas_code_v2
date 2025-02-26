setwd("/Volumes/RBMB/MouseAtlas/")
.libPaths("/Volumes/RBMB/library/")
library(limma)
library(nlme)
library(psych)
library(heatmap3)
library(gplots)
library(dplyr)
library(ggrepel)

library(nnls)
library(pheatmap)
library(ggplot2)
library(data.table)
library(ggpubr)
library(edgeR)
library(DESeq2)
library(glue, lib.loc = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library(ggfortify, lib.loc = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library(Seurat)
library(reshape2)


pred=read.csv("Predicted.csv",row.names=1,header=T)
act=read.csv("Actual.csv",row.names=1,header=T)
combined=read.csv("Combined.csv",row.names=1,header=T)

data <- reshape2::melt(as.matrix(act))
colnames(data)[3]="abundance"
data$group="actual"
data2<-melt(as.matrix(pred))
colnames(data2)[3]="abundance"
data2$group="prediction"

df=rbind(data,data2)
colnames(df)=c("sample","celltypes","abundance","group")

ggplot(df,aes(x=celltypes,y=abundance,fill=group))+geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_point(aes(y=abundance,fill=group),  size=1, shape=21
             ,position=position_jitterdodge())+
  stat_compare_means(aes(label=after_stat(p.signif)),
                    method="wilcox.test")+theme_classic()+
  theme(axis.text = element_text(size=12))

ggsave("plots/celldecon_actvpred.png",width=2500,height=1500,units="px",dpi=300)
