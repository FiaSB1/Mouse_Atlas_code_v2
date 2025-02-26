import scanpy as sc
import pandas as pd
import os
os.chdir("/shared/homes/155227/work/R/data/senani")

adata=sc.read("vieira19_Bronchi_anonymised.processed.h5ad")
adata.obsm['X_umap']=adata.obsm['X_umap_hm']

deg=pd.read_csv("tT2-Asthmaonlysig.csv")

up=deg[deg["legend"]=="Up"]
down=deg[deg["legend"]=="Down"]

upgenes=up.iloc[:, 0].astype(str).tolist()
downgenes=down.iloc[:, 0].astype(str).tolist()
deggenes=deg.iloc[:, 0].astype(str).tolist()


sc.tl.score_genes(adata,upgenes, score_name="upgenes")
sc.tl.score_genes(adata,downgenes, score_name="downgenes")
sc.tl.score_genes(adata,deggenes, score_name="deggenes")


sc.pl.umap(adata,use_raw=False,color=['upgenes'],frameon=False,save="UpregulatedSignature.png")
sc.pl.umap(adata,use_raw=False,color=['downgenes'],frameon=False,save="_DownregulatedSignature.png")
sc.pl.umap(adata,use_raw=False,color=['deggenes'],frameon=False,save="_DEGSignature.png")

sc.pl.umap(adata, use_raw=False, legend_loc='on data', legend_fontsize='medium',
           legend_fontweight='normal', color=['CellType'],cmap = 'bwr', vmin=-1, vmax=1,
           save="_Dimplot.png")