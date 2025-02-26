setwd("/Volumes/RBMB/MouseAtlas/")
.libPaths("/Volumes/RBMB/library/")

library(Matrix)
library(Seurat)
library(dplyr)
library(SeuratObject)
library(sctransform)
library(DCGL)
library(scales)
library(viper)
library(PISCES)
#BiocManager::install("glmGamPoi")
library("glmGamPoi")
library(uwot)
library(ggplot2)
library(Matrix.utils)
library(monocle)
library(mclust)


s.obj=readRDS("macrophage.rds")

s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^mt-")


s.obj <- SCTransform(s.obj, method = "glmGamPoi", variable.features.n = 25000, 
                     vars.to.regress = "percent.mt", verbose = TRUE)

# Distance matrix
s.obj <- CorDist(s.obj)

saveRDS(s.obj,"macrophage_postcordist.rds")

s.obj <- FindVariableFeatures(s.obj, selection.method = "vst", nfeatures = 5000)

s.obj <- FindVariableFeatures(s.obj, selection.method = "vst", nfeatures = 5000)

TF_list = read.table("MouseTF.txt")


new_genes = c(VariableFeatures(s.obj), TF_list$V1); length(new_genes)
table(duplicated(new_genes))

new_genes = unique(new_genes)
table(rownames(s.obj@assays$RNA@counts) %in% new_genes)

idx = rownames(s.obj@assays$RNA@counts) %in% new_genes


s.obj[["Var_counts"]] <- CreateAssayObject(counts = s.obj@assays$RNA@counts[idx, ])

#Creating new function
MetaCells = function (counts.mat, dist.mat, clust.vect, num.neighbors, 
                      subset, min.samps) 
{
  counts.mat <- as.matrix(counts.mat)
  dist.mat <- as.matrix(dist.mat)
  if (missing(clust.vect)) {
    clust.vect <- rep(1, ncol(counts.mat))
    names(clust.vect) <- colnames(counts.mat)
  }
  clust.labels <- sort(unique(clust.vect))
  meta.mats <- list()
  for (cl in clust.labels) {
    #clust.samps <- names(clust.vect)[which(clust.vect == cl)]
    clust.samps <- clust.vect == cl
    names(clust.samps) = names(clust.vect)
    if (length(clust.samps) > min.samps) {
      print(paste("Making metacell matrix for cluster ", 
                  cl, "...", sep = ""))
      clust.counts <- counts.mat[, clust.samps]
      print(dim(clust.counts))
      clust.dist <- dist.mat[names(which(clust.samps == TRUE)), names(which(clust.samps == TRUE))]
      print(dim(clust.dist))
      knn.mat <- KNN(clust.dist, k = num.neighbors)
      
      subset = ifelse(table(s.obj$Cell.Type == cl)[[2]] > min.samps, 250, (table(s.obj$Cell.Type == cl)[[2]]-20))
      print(subset)
      
      if (is.null(subset)) {
        sub.samps <- clust.samps
        subset <- length(sub.samps)
      }
      else {
        set.seed(42)
        sub.samps <- sample(names(which(clust.samps == TRUE)), subset)
        print(length(sub.samps))
      }
      imp.mat <- matrix(0L, nrow = nrow(clust.counts), 
                        ncol = subset)
      rownames(imp.mat) <- rownames(counts.mat)
      colnames(imp.mat) <- sub.samps
      for (ss in sub.samps) {
        neighbor.vect <- c(ss, rownames(knn.mat)[knn.mat[ss, 
        ]])
        ss.mat <- clust.counts[, neighbor.vect]
        imp.mat[, ss] <- rowSums(ss.mat)
      }
      imp.mat <- CPMTransform(imp.mat)
      meta.mats[[cl]] <- imp.mat
    }
  }
  return(meta.mats)
}


Clus.vec = as.factor(s.obj$new.Celltype)


idx = rownames(s.obj@assays$RNA@counts) %in% new_genes


set.seed(42)
meta.mats <- MetaCells(counts.mat = as.matrix(s.obj@assays$Var_counts@counts),
                       dist.mat = s.obj@assays$SCT@misc$dist.mat, 
                       clust.vect = Clus.vec,
                       subset = 250, min.samps = 499, #default min.samps=499
                       num.neighbors = 10)
sapply(meta.mats,dim)

#Save metaCells#######
setwd("ARACNe_input")
for (m.name in names(meta.mats)) {
  f.name <- paste('PH_', m.name, '_meta.rds', sep = '')
  saveRDS(meta.mats[[m.name]], f.name)
}

##### Convert to .tsv format #####

# Function to convert .rds tables to .tsv files (required for ARACNe)
ARACNeTable = function (dat.mat, out.file, cell.no, subset = TRUE) 
{
  dat.mat <- dat.mat[!duplicated(rownames(dat.mat)), ]
  if (subset) {
    dat.mat <- dat.mat[, sample(colnames(dat.mat), min(ncol(dat.mat), 
                                                       cell.no))]
  }
  sample.names <- colnames(dat.mat)
  gene.ids <- rownames(dat.mat)
  m <- dat.mat
  mm <- rbind(c("gene", sample.names), cbind(gene.ids, 
                                             m))
  write.table(x = mm, file = out.file, sep = "\t", quote = F, 
              row.names = F, col.names = F)
}

# Convert to .tsv
for (i in names(meta.mats)) {
  in1 = readRDS(paste("PH_", i, "_meta.rds",sep=""))
  ARACNeTable(in1, paste("PH_", i, "_meta.tsv",sep=""), cell.no = ncol(in1), subset = FALSE)
}


####Run ARACNe in Linux#########
#1. Create Directory = ARACNe_input and ARACNe_output
#2. Run Aracne script by paragraph

# ##ARACNe output######
# ct=list.dirs("ARACNe_output/",recursive=T,full.names=T) ##list full content of directory
# ct2 <- strsplit(ct, split = "/")
# ct3 <- sapply(ct2, `[`, 9)
# ct3=ct3[-1]
# 
# # Convert network .txt into .tsv
# setwd("ARACNe_output")
# for (i in ct3) {
#   anet = read.delim(paste(i, "/network.txt", sep="")); dim(anet)
#   write.table(anet, file=paste(i, "/network.tsv", sep=""), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
#   print(i)
# }
# 
# 
# # Create cell type-specific Regulon objects
# 
# for (i in ct3) {
#   S1 = readRDS(paste0("ARACNe_input/",i, ".rds"))
#   RegProcess(paste0("ARACNe_output/",i, "/network.tsv"), S1, 
#              out.dir = 'Regs/', out.name = paste(i, "_Net_", sep=""))
# }
# 
# 
# # Load regulons and compile into a list
# reg_networks=list()
# for (i in 1:length(ct3)){
#   nam = ct3[i]
#   print(nam)
#   r.net <- readRDS(paste0("Regs/",nam, "_Net_pruned.rds"))
#   reg_networks[[i]]<-r.net
# }
# names(reg_networks)=ct3
# saveRDS(reg_networks,"Regs/regulon_networks.rds")
# 
# 
# ##### Run PISCES #####
# 
# s.obj = readRDS("MouseAtlas_preARACNe.rds") 
# 
# DefaultAssay(s.obj) = "RNA"
# s.obj
# TF_list = read.table("MouseTF.txt")
# 
# idx=rownames(s.obj@assays$Var_counts)
# 
# reg_networks=readRDS("Regs/regulon_networks.rds")
# 
# 
# ## Create PISCES slot
# s.obj[["PISCES"]] <- CreateAssayObject(counts = s.obj@assays$RNA@counts[idx, ] )
# s.obj@active.assay <- "PISCES"
# dim(s.obj@assays$PISCES@counts)
# 
# # Transform by counts per million
# s.obj <- CPMTransform(s.obj)
# 
# # Create Gene expression signature
# s.obj <- GESTransform(s.obj)
# 
# # Run viper() function of individual cells 
# s.obj <- PISCESViper(data.obj=s.obj, net.list = reg_networks)
# dim(s.obj@assays$PISCES@scale.data)
# 
# # pisces.obj@assays$PISCES@scale.data was used for master regulator analysis ("PISCES_TF")
# s.obj[["PISCES_TF"]]<-CreateAssayObject(counts=s.obj@assays$PISCES@scale.data)
# 
# saveRDS(s.obj,"MouseAtlas_postPISCES.rds")
# 
# 
