library(Seurat)

adata=anndata::read_h5ad("data/MouseAtlas/FibrosisAge_2207.h5ad")
s.obj=CreateSeuratObject(counts=t(as.matrix(adata$X)),meta.data=adata$obs)

s.obj=readRDS("data/MouseAtlas/FibrosisAge_2207_incomplete.rds")
#s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^mt-")
s.obj <- NormalizeData(s.obj, normalization.method = "LogNormalize", scale.factor = 10000)
s.obj <- FindVariableFeatures(s.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(s.obj), 10)

all.genes <- rownames(s.obj)

#s.obj<-readRDS("/Volumes/RBMB/ALI/data/merged_postQC.rds")

##Friday 27 May last checkpoint,next is scaledata
s.obj <- ScaleData(s.obj, features = all.genes)

s.obj <- RunPCA(s.obj, features = VariableFeatures(object = s.obj))

#DimPlot(s.obj,reduction="pca")

ElbowPlot(s.obj)

s.obj <- FindNeighbors(s.obj, dims = 1:20)
s.obj <- FindClusters(s.obj, resolution = 0.2)

s.obj <- RunUMAP(s.obj, dims = 1:20)
#DimPlot(s.obj, reduction = "umap")

s.obj=DietSeurat(s.obj,counts=T,data=T,scale.data=F,features=NULL,assays="RNA",
                        dimreducs=names(s.obj@reductions),graphs=names(s.obj@graphs),misc=T)





#cell type proportion####

proportion=data.frame(unclass(table(fibro$Disease,fibro$Diff_Exp)))
#proportion <- proportion %>% select(-Blood.Cells) #Remove Blood Cells
total=as.data.frame(rowSums(proportion))
percent=(proportion/total[,1])*100
#percent=percent[startsWith(colnames(percent),"Viral")]

percent$id=rownames(percent)

data <- reshape2::melt(percent,id="id") 
colnames(data) <- c("Group", "Celltypes", "Percentage")
data <- data.frame(data)
#data=data[startsWith(as.character(data$Celltypes),"Viral"),]
#test=data[data$Group=="AgeBA5",]
data$Group=as.factor(data$Group)
#data$Group=factor(data$Group,levels=c("Yngsham","YngBA5","YngDel","YngWuh","Agesham","AgeDel","AgeBA5","AgeWuh"))


pal <- c("cornflowerblue", "tomato", "mediumseagreen", "gold", "purple","lightgreen", 
         "lightskyblue", "lightcoral", "mediumorchid", "darkgreen", "darkslateblue",
         "firebrick", "steelblue", "yellowgreen", "deepskyblue", 
         "darkturquoise", "dodgerblue", "mediumspringgreen", "springgreen", 
         "salmon", "plum", "violet", "goldenrod", "darkorange","darkred")
pal_viral<-c("plum", "violet", "goldenrod", "darkorange","darkred")
#last_five_colors <- tail(pal, 5)
# Mute the first colors by desaturating and lightening them
#muted_colors <- colorspace::desaturate(lighten(head(pal, -5), amount = 0.3), amount = 0.3)
# Combine the muted colors with the last five original colors
#final_palette <- c(muted_colors, last_five_colors)
data$Celltypes=gsub("[.]"," ",data$Celltypes)

ggplot(data, aes(fill=Celltypes, y=Percentage, x=Group)) + 
  geom_bar(position = "stack", stat="identity", color = "black", alpha = 0.9) +
  scale_fill_manual(values = pal) +
  ylab('Proportion of Cells (%)') + xlab('')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,size = 18),plot.title = element_text(size = 24, face = "bold"),
        axis.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.text  = element_text(color = "black", size = 16, angle = 0, hjust = 1, vjust = 0, face = "bold"),
        legend.text = element_text(color = "black", size = 16, face = "bold"),
        legend.title  = element_text(color = "black", size = 16, face = "bold"))+NoLegend()
