kp=sc.read("/shared/homes/155227/work/R/data/MouseAtlas/GSE252663_KP_final.h5ad")
KP=kp.copy()
sc.pp.filter_cells(KP,min_genes=200)
sc.pp.filter_genes(KP,min_cells=10)
sc.pp.normalize_total(KP, target_sum=1e4)
sc.pp.log1p(KP)
sc.pp.highly_variable_genes(KP, n_top_genes=6000,min_mean=0.06)
hvgkp=KP.var_names[KP.var['highly_variable']]
len(hvgkp)


kp.obs['study']='kp'
adata=flu.concatenate(kp,hyper,join='outer',fill_value=0,index_unique=None)
study_subset = adata[adata.obs['study'] == study].copy()
sc.pp.filter_cells(study_subset,min_genes=200)
sc.pp.filter_genes(study_subset,min_cells=10)
sc.pp.normalize_total(study_subset, target_sum=1e4)
sc.pp.log1p(study_subset)
sc.pp.highly_variable_genes(study_subset, n_top_genes=6000,min_mean=0.06)
hvgadata=study_subset.var_names[study_subset.var['highly_variable']]
len(hvgadata.intersection(hvgkp))

