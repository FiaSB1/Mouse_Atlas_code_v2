
library(monocle)
library(Seurat)



s.obj=readRDS("macrophage.rds")

total_cells <- ncol(s.obj)

subset_size <- round(total_cells/4)

random_indices <- sample(1:total_cells, size = subset_size, replace = FALSE)

s.obj2 <- s.obj[, random_indices]

pd <- new("AnnotatedDataFrame", data = data.frame(s.obj@meta.data))
gm = data.frame(geneID = (rownames(s.obj@assays$RNA@counts)), 
                gene_short_name = rownames(s.obj@assays$RNA@counts))
rownames(gm) = gm$gene_short_name
gm <- new("AnnotatedDataFrame", data = gm)
Mon1 <- newCellDataSet(as.matrix(s.obj@assays$RNA@counts),
                       phenoData = pd, featureData = gm,
                       expressionFamily=negbinomial.size())

# Run preprocessing
Mon1 <- estimateSizeFactors(Mon1)
Mon1 <- estimateDispersions(Mon1)
Mon1 <- detectGenes(Mon1, min_expr = 0.1)


##with diff test - include MS samples
diff_test_res <- differentialGeneTest(Mon1, fullModelFormulaStr = "~Group")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
Mon1 <- setOrderingFilter(Mon1, ordering_genes)

# Dim reduction
set.seed(42); Mon1 <- reduceDimension(Mon1, max_components = 2, method = 'DDRTree')


# Order cells from control
Mon1 <- orderCells(Mon1)
#Mon2<-HSMM
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Disease)[,"Control"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
Mon1 <- orderCells(Mon1, root_state = GM_state(Mon2))

saveRDS(Mon1,"monocle_macrophage.rds")

