setwd("~/work/R/data/MouseAtlas")
.libPaths("~/work/R/library")
library(limma)
library(nlme)
#library(psych)
#library(heatmap3)
#library(gplots)
library(dplyr)
library(ggrepel)
library(SeuratDisk)
#library(nnls)
library(pheatmap)
library(ggplot2)
library(data.table)
library(ggpubr)
library(edgeR)
library(DESeq2)
# library(glue, lib.loc = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
# library(ggfortify, lib.loc = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library")
library(Seurat)


#Load and preprocess Ref####
ref=readRDS("/shared/ci/kyle/FinalObject_Genotype_final.rds")
ref <- SCTransform(ref, verbose = FALSE)
ref <- RunSPCA(ref, assay = 'RNA',graph='RNA_nn')
reference <- FindNeighbors(
  object = ref,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)


#Load GSE151974 Hyperoxia########
meta=read.csv('GSE151974_cell_metadata_postfilter.csv',row.names=1)
exp=fread('GSE151974_raw_umi_matrix_postfilter.csv.gz')
exp=as.data.frame(exp)
rownames(exp)=exp[,1];exp=exp[,-1]

s.obj=CreateSeuratObject(counts=exp,meta.data=meta)
s.obj <- NormalizeData(s.obj, normalization.method = "LogNormalize", scale.factor = 10000)
s.obj <- FindVariableFeatures(s.obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(s.obj)
s.obj <- ScaleData(s.obj)
s.obj <- RunPCA(s.obj, features = VariableFeatures(object = s.obj))
ElbowPlot(s.obj)
s.obj <- FindNeighbors(s.obj, dims = 1:25)
s.obj <- FindClusters(s.obj, resolution = 0.6)
s.obj <- RunUMAP(s.obj, dims = 1:25)
#saveRDS(s.obj,"GSE151974_Hyperoxia_processed.rds")

query<-RunSPCA(s.obj, assay = 'RNA',graph='RNA_nn')
#query=RunUMAP(query, assay="RNA",dim=1:20)
#query = ScaleData(query, assay = 'RNA')
#query = RunSPCA(query, assay = 'RNA', graph = 'RNA_nn')
# query$study = query$orig.ident
# query$orig.ident = query[["barcode"]]

# query$study = query$orig.ident
# query$orig.ident = query[["barcode"]]

query2 <- SplitObject(query, split.by = "Sample")
#query2 <- lapply(X = query2, FUN = NormalizeData, verbose = FALSE)

DefaultAssay(reference)="RNA"

anchors <- list()
for (i in 1:length(query2)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = query2[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

for (i in 1:length(query2)) {
  query2[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query2[[i]],
    reference = ref, 
    refdata = list(
      celltype = "cluster_labels_res.0.6",
      newcelltype ="new_cell_type"),
    reference.reduction = "spca",
    reduction.model = NULL
  )
}

query_mapped <- merge(query2[[1]], query2[2:length(query2)],merge.dr="umap")
query_mapped=DietSeurat(query_mapped,counts=T,data=F,scale.data=F,features=NULL,assays="RNA",
                         dimreducs=names(s.obj@reductions),graphs=names(s.obj@graphs),misc=T)

##Harmonising metadata####
query_mapped$Disease="Hyperoxia"
query_mapped$platform="10X Genomics"
query_mapped$strain="C57BL/6**********"
query_mapped$sex



saveRDS(query_mapped,"GSE151974_Hyperoxia_Final.rds")
SaveH5Seurat(query_mapped,"GSE151974_Hyperoxia_Final.h5seurat")

.libPaths("/Volumes/RBMB/library/")
library(Seurat)
s.obj=readRDS("/Volumes/RBMB/MouseAtlas/GSE151974_Hyperoxia_Final.rds")


#Load GSE225552 Influenza#####
meta=as.data.frame(fread('GSE225552_Metadata_Devarajan_etal_analyses.csv.gz'))
files<-list.files(path="GSE225552_RAW", pattern="*_filtered_matrix.mtx", full.names=TRUE)
filenames=files
filenames=strsplit(filenames,"/")
filenames=sapply(filenames,`[`,2)
filenames=strsplit(filenames,"_")
filenames <- lapply(filenames, function(x) head(x, -2)) #remove the last string after strsplit
filenames <- sapply(filenames, function(parts) paste(parts, collapse = "_")) #combining all the splitted strings
 

query2=list()

for (i in 1:length(filenames)){
  
  data <- ReadMtx(
    mtx = paste0("GSE225552_RAW/",filenames[i],"_filtered_matrix.mtx"), 
    features = paste0("GSE225552_RAW/",filenames[i],"_filtered_features.tsv"),
    cells = paste0("GSE225552_RAW/",filenames[i],"_filtered_barcodes.tsv"),
  )
  dim(data)
  s.obj<-CreateSeuratObject(counts = data)
  s.obj <- NormalizeData(s.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  s.obj <- FindVariableFeatures(s.obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(s.obj)
  s.obj <- ScaleData(s.obj)
  s.obj <- RunPCA(s.obj, features = VariableFeatures(object = s.obj))
  ElbowPlot(s.obj)
  s.obj <- FindNeighbors(s.obj, dims = 1:25)
  s.obj <- FindClusters(s.obj, resolution = 0.6)
  s.obj <- RunUMAP(s.obj, dims = 1:25)
  query<-RunSPCA(s.obj, assay = 'RNA',graph='RNA_nn')
  query2[[i]]=query
  names(query2)[i]=filenames[i]
}


DefaultAssay(reference)="RNA"

anchors <- list()
for (i in 1:length(query2)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = query2[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

for (i in 1:length(query2)) {
  query2[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query2[[i]],
    reference = ref, 
    refdata = list(
      celltype = "cluster_labels_res.0.6",
      newcelltype ="new_cell_type"),
    reference.reduction = "spca",
    reduction.model = NULL
  )
}

query_mapped <- merge(query2[[1]], query2[2:length(query2)],merge.dr="umap")
query_mapped=DietSeurat(query_mapped,counts=T,data=T,scale.data=F,features=NULL,assays="RNA",
                        dimreducs=names(s.obj@reductions),graphs=names(s.obj@graphs),misc=T)
barcodes=colnames(query_mapped)
barcodes=strsplit(barcodes,"_")
barcodes=sapply(barcodes,`[`,1)
colnames(query_mapped@assays$RNA@counts)=barcodes
colnames(query_mapped@assays$RNA@data)=barcodes



s.obj=CreateSeuratObject(counts=query_mapped@assays$RNA@counts)
s.obj@meta.data=query_mapped@meta.data
s.obj=s.obj[,colnames(s.obj)%in%meta$barcode]
meta2=meta[!duplicated(meta$barcode),]
s.obj@meta.data=cbind(s.obj@meta.data,meta2)

saveRDS(s.obj,"GSE225552_Influenza_final.rds")



#Load GSE252663 KP######
meta=read.delim("GSE252663_series_matrix.txt")
files<-list.files(path="GSE252663_RAW", pattern="*_barcodes.tsv.gz", full.names=TRUE)
filenames=files
filenames=strsplit(filenames,"/")
filenames=sapply(filenames,`[`,2)
filenames=strsplit(filenames,"_")
samplenames=sapply(filenames,`[`,1)
filenames <- lapply(filenames, function(x) head(x, -1)) #remove the last string after strsplit
filenames <- sapply(filenames, function(x) paste(x, collapse = "_")) #combining all the splitted strings


query2=list()

for (i in 1:length(filenames)){
  
  data <- ReadMtx(
    mtx = paste0("GSE252663_RAW/",filenames[i],"_matrix.mtx.gz"), 
    features = paste0("GSE252663_RAW/",filenames[i],"_features.tsv.gz"),
    cells = paste0("GSE252663_RAW/",filenames[i],"_barcodes.tsv.gz"),
  )
  head(data)
  s.obj<-CreateSeuratObject(counts = data)
  s.obj <- NormalizeData(s.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  s.obj <- FindVariableFeatures(s.obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(s.obj)
  s.obj <- ScaleData(s.obj)
  s.obj <- RunPCA(s.obj, features = VariableFeatures(object = s.obj))
  ElbowPlot(s.obj)
  s.obj <- FindNeighbors(s.obj, dims = 1:25)
  s.obj <- FindClusters(s.obj, resolution = 0.6)
  s.obj <- RunUMAP(s.obj, dims = 1:25)
  query<-RunSPCA(s.obj, assay = 'RNA',graph='RNA_nn')
  query$samples=samplenames[i]
  query2[[i]]=query
  names(query2)[i]=filenames[i]
}

DefaultAssay(reference)="RNA"

anchors <- list()
for (i in 1:length(query2)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = query2[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

for (i in 1:length(query2)) {
  query2[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query2[[i]],
    reference = ref, 
    refdata = list(
      celltype = "cluster_labels_res.0.6",
      newcelltype ="new_cell_type"),
    reference.reduction = "spca",
    reduction.model = NULL
  )
}
#merging at the end because reference mapping involves splitting the dataset by samples
query_mapped <- merge(query2[[1]], query2[2:length(query2)],merge.dr="umap")
query_mapped=DietSeurat(query_mapped,counts=T,data=T,scale.data=F,features=NULL,assays="RNA",
                        dimreducs=names(s.obj@reductions),graphs=names(s.obj@graphs),misc=T)

barcodes=colnames(query_mapped)
barcodes=strsplit(barcodes,"_")
barcodes=sapply(barcodes,`[`,1)
colnames(query_mapped@assays$RNA@counts)=barcodes
colnames(query_mapped@assays$RNA@data)=barcodes


s.obj=CreateSeuratObject(counts=query_mapped@assays$RNA@counts)
s.obj@meta.data=query_mapped@meta.data

meta2=as.data.frame(t(meta));colnames(meta2)=meta2[1,];meta2=meta2[-1,]
s.obj$group=s.obj$samples
s.obj$genotype=s.obj$samples
s.obj$age="12 Weeks"

for (i in 1:length(unique(s.obj$samples))){
    s.obj$group=gsub(unique(s.obj$samples)[i],rownames(meta2[meta2[,1]==unique(s.obj$samples)[i],]),s.obj$group)
    s.obj$genotype=gsub(unique(s.obj$samples)[i],meta2[meta2[,1]==unique(s.obj$samples)[i],10],s.obj$genotype)
    s.obj$genotype=gsub("genotype: ","",s.obj$genotype)
    s.obj$genotype=gsub("Stat1-/-","KO",s.obj$genotype)
    s.obj$infection=s.obj$group
    s.obj$infection=gsub(".*KP.*","KP",s.obj$infection)
    s.obj$infection=gsub(".*Contrl.*","Control",s.obj$infection)
    
}

s.obj$sex=s.obj$group
s.obj$sex=gsub("KO.KP2|WT2.KP1","Male",s.obj$sex)
s.obj$sex=gsub("KO.Contrl2|KO.KP1|WT.Contrl1|WT6.KP2","Female",s.obj$sex)



#Load Fibrosis Age#####
meta=read.delim("FibroAge_meta_data.txt")
files<-list.files(path="FibroAge_counts", pattern="*dge.txt", full.names=TRUE)
filenames=files
filenames=strsplit(filenames,"/")
filenames=sapply(filenames,`[`,2)
filenames=strsplit(filenames,"[.]")
samplenames=sapply(filenames,`[`,1)
#filenames <- lapply(filenames, function(x) head(x, -1)) #remove the last string after strsplit
#filenames <- sapply(filenames, function(x) paste(x, collapse = "_")) #combining all the splitted strings

obj=list()

for(i in 1:length(samplenames)){
  counts=read.delim(gzfile(paste0("FibroAge_counts/",samplenames[i],".dge.txt.gz")),header=T,sep='\t');print("counts loaded")
  rownames(counts)=counts$GENE;counts=counts[,-1]
  s.obj<-CreateSeuratObject(counts = counts)
  s.obj$sample=samplenames[i]
  
  obj[[i]]=s.obj
}

#merge was done prior to ref mapping because initially thought ref mapping was not necessary
s.obj <- merge(obj[[1]], obj[2:length(obj)],merge.dr="umap")
s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^mt-")
VlnPlot(s.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0,ncol = 3)
test <- subset(s.obj, subset = nCount_RNA < 7000 & percent.mt<15)

s.obj <- NormalizeData(s.obj, normalization.method = "LogNormalize", scale.factor = 10000)
s.obj <- FindVariableFeatures(s.obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(s.obj)
s.obj <- ScaleData(s.obj)
s.obj <- RunPCA(s.obj, features = VariableFeatures(object = s.obj))
ElbowPlot(s.obj)
s.obj <- FindNeighbors(s.obj, dims = 1:25)
s.obj <- FindClusters(s.obj, resolution = 0.6)
s.obj <- RunUMAP(s.obj, dims = 1:25)
s.obj<-RunSPCA(s.obj, assay = 'RNA',graph='RNA_nn')



query <- SplitObject(s.obj, split.by = "sample")
query <- lapply(X = query, FUN = NormalizeData, verbose = FALSE)

DefaultAssay(reference)="RNA"

anchors <- list()
for (i in 1:length(query)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = query[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

for (i in 1:length(query)) {
  query[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query[[i]],
    reference = ref, 
    refdata = list(
      celltype = "cluster_labels_res.0.6",
      newcelltype ="new_cell_type"),
    reference.reduction = "spca",
    reduction.model = NULL
  )
}

query_mapped <- merge(query[[1]], query[2:length(query)],merge.dr="umap")
query_mapped=DietSeurat(query_mapped,counts=T,data=T,scale.data=F,features=NULL,assays="RNA",
                        dimreducs=names(s.obj@reductions),graphs=names(s.obj@graphs),misc=T)


for (m in colnames(meta)){
  query_mapped@meta.data[,m]=query_mapped$sample
  
  for (i in 1:length(unique(test$sample))){
  query_mapped@meta.data[,m]=gsub(unique(query_mapped$sample)[i],meta[meta[,1]==unique(query_mapped$sample)[i],m],query_mapped@meta.data[,m])
  }
}