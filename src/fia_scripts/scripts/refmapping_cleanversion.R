#Load and preprocess Reference####
ref=readRDS("FinalObject_Genotype_final.rds")
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


#Load meta data of query data
meta=as.data.frame(fread('meta_data.csv'))

#List each sample filenames of query data
files<-list.files(path="RAW", pattern="*_matrix.mtx", full.names=TRUE)
filenames=files
filenames=strsplit(filenames,"/")
filenames=sapply(filenames,`[`,2)
filenames=strsplit(filenames,"_")
filenames <- lapply(filenames, function(x) head(x, -2)) #remove the last string after strsplit
filenames <- sapply(filenames, function(parts) paste(parts, collapse = "_")) #combining all the splitted strings

#Load raw count matrix of each sample and process as a seurat object in a loop
query_obj=list()

for (i in 1:length(filenames)){
  
  data <- ReadMtx(
    mtx = paste0("RAW/",filenames[i],"_matrix.mtx"), 
    features = paste0("RAW/",filenames[i],"_features.tsv"),
    cells = paste0("RAW/",filenames[i],"_barcodes.tsv"),
  )
  dim(data)
  s.obj<-CreateSeuratObject(counts = data)
  
  #Seurat basic processing steps 
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
  query_obj[[i]]=query
  names(query_obj)[i]=filenames[i]
}


DefaultAssay(reference)="RNA"

#Find anchors of each sample 
anchors <- list()
for (i in 1:length(query_obj)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = reference,
    query = query_obj[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}


#Map the query data from the reference
for (i in 1:length(query_obj)) {
  query_obj[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query_obj[[i]],
    reference = ref, 
    refdata = list(
      celltype = "cluster_labels_res.0.6",
      newcelltype ="new_cell_type"),
    reference.reduction = "spca",
    reduction.model = NULL
  )
}


#Merge each reference mapped samples into one seurat object
query_mapped <- merge(query_obj[[1]], query_obj[2:length(query_obj)],merge.dr="umap")

#'Diet' seurat object to remove scaled data - integration only requires raw count matrix and meta data
query_mapped=DietSeurat(query_mapped,counts=T,data=T,scale.data=F,features=NULL,assays="RNA",
                        dimreducs=names(s.obj@reductions),graphs=names(s.obj@graphs),misc=T)

#Add metadata to the Seurat Object 
s.obj@meta.data=cbind(s.obj@meta.data,meta2)

saveRDS(s.obj,"mapped_query.rds")


#If an error comes up regarding cell names
# barcodes=colnames(query_mapped)
# barcodes=strsplit(barcodes,"_")
# barcodes=sapply(barcodes,`[`,1)
# colnames(query_mapped@assays$RNA@counts)=barcodes
# colnames(query_mapped@assays$RNA@data)=barcodes
# 
# 
# 
# s.obj=CreateSeuratObject(counts=query_mapped@assays$RNA@counts)
# s.obj@meta.data=query_mapped@meta.data
# s.obj=s.obj[,colnames(s.obj)%in%meta$barcode]
# meta2=meta[!duplicated(meta$barcode),]
# s.obj@meta.data=cbind(s.obj@meta.data,meta2)


