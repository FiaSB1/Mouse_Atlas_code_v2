# folder and files
output_dir = "results"
output_file = "output.h5ad"

import fgread
import os
import warnings
import shutil
os.chdir('/shared/ci/kyle')

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)

##in case if HPC won't let me map the entire dataset - do it per samples
#adata=sc.read("/shared/homes/155227/work/R/data/Collaborator_Xu/humanCOPD_LC_diet.h5ad")

#adata.X.data

#adata.X=sparse.csr_matrix(adata.X,dtype="float32")

#del adata.raw
#del adata.obsm
#del adata.varm
#del adata.var

#adata.obs['dataset']="My_dataset"

#samples=adata.obs.Sample.unique()

#move each samples to dataset_0002 - pipeline only wants one query at a time
#for s in samples:
#    src_path=f"/shared/homes/155227/work/R/data/Collaborator_Xu/samples/{s}.h5ad"
#    dst_path=f"/shared/homes/155227/work/R/data/HLCA/data/dataset_0002/{s}.h5ad"
#    shutil.move(src_path, dst_path)


# make sure exactly two datasets are attached
assert len(os.listdir("data/")) == 2,"ERROR: There is more than one query dataset attached to this analysis."

# find HLCA and query data folder
HLCA_emb_and_metadata = "MouseAtlas_Healthy_emb.h5ad" # HLCA_emb_and_metadata.h5ad"
query_path = "/shared/ci/Hansbro_Group/Datasets/scRNASeq/aged_young_mice/Objects/"
hlca_path = "/shared/ci/kyle/Control_ref_model_2k_seurat"

#if not os.path.exists(HLCA_emb_and_metadata):
#    HLCA_emb_and_metadata = "data/dataset_0002/HLCA_core_v1.1_emb.h5ad" # HLCA_emb_and_metadata.h5ad"
#    query_path = "data/dataset_0001"
#    hlca_path = "data/dataset_0002"

# make sure there is only one expression data file in query path
#df = fgread.ds_info(data_dir= '/shared/homes/155227/work/R/data/HLCA/data')
#expr_data = df.expressionDataFileInfos[df.path==query_path]
#assert len(expr_data.values[0]) == 1,"ERROR: There is more than one expression data file in the query data."
#test_dataset = query_path + "/" + expr_data.values[0][0].get("name")

#expr_data = fgread.ds_info().expressionDataFileInfos[df.path==query_path]
#expr_data.values[0]


#!rm -r /fastgenomics/analysis/$output_dir
#!mkdir /fastgenomics/analysis/$output_dir

import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import os
import scarches as sca
from scipy import sparse

sc.set_figure_params(figsize=(5,5))

sc.logging.print_versions()

# reference embedding, including cell/sample/subject metadata:
reference_embedding = sc.read_h5ad(HLCA_emb_and_metadata)
# path to scArches reference model
reference_model_path = hlca_path
# gene order for scArches model
adata=sc.read("/shared/ci/kyle/Control_clustering_comp.h5ad")
reference_gene_order = adata[:,adata.var['highly_variable']].var_names.tolist()
reference_gene_order=pd.DataFrame(reference_gene_order)
reference_gene_order.columns = ['gene_symbol']

#input of query data
query_data_full = sc.read_h5ad("/shared/ci/Hansbro_Group/Datasets/scRNASeq/aged_young_mice/Objects/mouse.cov.h5ad")
query_data_full.obs['study']='aged'

def subset_and_pad_adata_object(adata, ref_genes_df, min_n_genes_included=1500):
    """Function to subset a dataset to the 2000 genes used for scArches mapping.
    Missing genes are added and set to 0 counts for all cells."""
    # delete obsm and varm to enable concatenation
    print('Deleting .obsm and .varm from data if present, to enable cleaning and padding.')
    del adata.obsm
    del adata.varm
    # test if adata.var.index has gene names or ensembl names:
    # n_ids_detected = sum(adata.var.index.isin(ref_genes_df.gene_id))
    # n_symbols_detected = sum(adata.var.index.isin(ref_genes_df.gene_symbol))
    # if max(n_ids_detected, n_symbols_detected) < min_n_genes_included:
    #     # change column names to lower case
    #     adata.var.columns = adata.var.columns.str.lower()
    #     # check if gene names are in another column:
    #     if "gene_symbols" in adata.var.columns:
    #         adata.var.index = adata.var.gene_symbol
    #         n_symbols_detected = sum(adata.var.index.isin(ref_genes_df.gene_symbol))
    #     elif "gene_ids" in adata.var.columns:
    #         adata.var.index = adata.var.gene_ids
    #         n_ids_detected = sum(adata.var.index.isin(ref_genes_df.gene_id))
    #     # check if anything changed:
    #     if max(n_ids_detected, n_symbols_detected) < min_n_genes_included:
    #         raise ValueError(f"We could detect only {max(n_ids_detected, n_symbols_detected)} genes of the 2000 that we need for the mapping! The minimum overlap is {min_n_genes_included}. Contact the HLCA team for questions. Exiting")
    # else:
    #     if n_ids_detected >= n_symbols_detected:
    #         gene_type = "gene_id"
    #         print("Gene names detected: ensembl gene ids.")
    #         n_genes = n_ids_detected
    #     else:
    #         gene_type = "gene_symbol"
    #         n_genes = n_symbols_detected
    #         print("Gene names detected: ensembl gene symbols.")
    gene_type="gene_symbol"
    n_symbols_detected = sum(adata.var.index.isin(ref_genes_df.gene_symbol))
    n_genes = n_symbols_detected
    genes = adata.var.index[adata.var.index.isin(ref_genes_df[gene_type])].tolist()
    # if not all genes are included, pad:
    if n_genes > 2000:
        raise ValueError("Your gene names appear not to be unique, something must be wrong. Exiting.")
    print(f"{n_genes} genes detected out of 2000 used for mapping.")
    # Subset adata object
    adata_sub = adata[:,genes].copy()
    # Pad object with 0 genes if needed
    if n_genes < 2000:
        diff = 2000 - n_genes
        print(f'Not all genes were recovered, filling in zeros for {diff} missing genes...')
        # Genes to pad with
        genes_to_add = list(set(ref_genes_df[gene_type].values).difference(set(adata_sub.var_names))) # added list() bit here to prevent error
        df_padding = pd.DataFrame(data=np.zeros((adata_sub.shape[0],len(genes_to_add))), index=adata_sub.obs_names, columns=genes_to_add)
        adata_padding = sc.AnnData(df_padding)
        # Concatenate object
        adata_sub = anndata.concat([adata_sub, adata_padding], axis=1, join='outer', index_unique=None, merge='unique')
    # order genes:
    adata_sub = adata_sub[:,ref_genes_df[gene_type]].copy()
    return adata_sub

query_data = subset_and_pad_adata_object(query_data_full, reference_gene_order)

if (query_data.var.index == reference_gene_order.gene_symbol).all() or (
    query_data.var.index == reference_gene_order.gene_id
).all():
    print("Gene order is correct.")
else:
    print(
        "WARNING: your gene order does not match the order of the HLCA reference. Fix this before continuing!"
    )


query_data.X = sparse.csr_matrix(query_data.X)

query_data.X[:10, :10].toarray()

query_data.obs["scanvi_label"] = "unlabeled"

if not 'dataset' in query_data.obs.columns:
    print("Setting your entire query data to a single batch. For more information on how to best set your batch variable, check out the HLCA mapping instruction page.")
    query_data.obs['dataset'] = "query_batch_1"

surgery_model = sca.models.SCANVI.load_query_data(
    query_data,
    reference_model_path,
    freeze_dropout=True,
)

batch_variable = "dataset"  # the column name under which you stored your batch variable
query_batches = sorted(query_data.obs[batch_variable].unique())
print(query_batches)

surgery_epochs = 50
early_stopping_kwargs_surgery = {
    "early_stopping_monitor": "elbo_train",
    "early_stopping_patience": 10,
    "early_stopping_min_delta": 0.001,
    "plan_kwargs": {"weight_decay": 0.0},
}

surgery_model.train(max_epochs=surgery_epochs, **early_stopping_kwargs_surgery)

surgery_model.save(output_dir, overwrite=True)


query_embedding = sc.AnnData(surgery_model.get_latent_representation(query_data))
query_embedding.obs = query_data.obs

query_embedding.obs['ref_or_query'] = "query"
reference_embedding.obs['ref_or_query'] = "ref"

combined_emb = sc.concat(
    (reference_embedding, query_embedding), index_unique=None, join="outer"
)  # index_unique="_", batch_key="ref_or_query")

for cat in combined_emb.obs.columns:
    if isinstance(combined_emb.obs[cat].values, pd.Categorical):
        pass
    elif pd.api.types.is_float_dtype(combined_emb.obs[cat]):
        pass
    else:
        print(
            f"Setting obs column {cat} (not categorical neither float) to strings to prevent writing error."
        )
        combined_emb.obs[cat] = combined_emb.obs[cat].astype(str)



sc.pp.neighbors(combined_emb, n_neighbors=30)
sc.tl.umap(combined_emb)
sc.pl.umap(combined_emb, color="ref_or_query", save="_MouseAtlas_and_query_cells_from_query.png")


# cts_ordered = pd.read_csv((f"{hlca_path}/HLCA_celltypes_ordered.csv"), index_col=0).rename(
#     columns={f"Level_{lev}": f"labtransf_ann_level_{lev}" for lev in range(1, 6)}
# )

#reference_embedding.obs = reference_embedding.obs.join(cts_ordered, on="ann_finest_level")

knn_transformer = sca.utils.knn.weighted_knn_trainer(
    train_adata=reference_embedding,
    train_adata_emb="X",  # location of our joint embedding
    n_neighbors=50,
)

labels, uncert = sca.utils.knn.weighted_knn_transfer(
    query_adata=query_embedding,
    query_adata_emb="X",  # location of our embedding, query_adata.X in this case
    label_keys="Level_4",  # (start of) obs column name(s) for which to transfer labels
    knn_model=knn_transformer,
    ref_adata_obs=reference_embedding.obs,
)

uncertainty_threshold = 0.2

labels.rename(
    columns={
        f"Level_4": f"ann_level_4_transferred_label_unfiltered"
    },
    inplace=True,
)
uncert.rename(
    columns={
        "Level_4": f"ann_level_4_transfer_uncert"
    },
    inplace=True,
)

combined_emb.obs = combined_emb.obs.join(labels['ann_level_4_transferred_label_unfiltered'])
combined_emb.obs = combined_emb.obs.join(uncert['ann_level_4_transfer_uncert'])

combined_emb.obs[f"ann_level_4_transfer_uncert"] = np.float64(combined_emb.obs[f"ann_level_4_transfer_uncert"])


combined_emb.obs[f"ann_level_4_transferred_label"] = combined_emb.obs[
    f"ann_level_4_transferred_label_unfiltered"
].mask(
    combined_emb.obs[f"ann_level_4_transfer_uncert"] > uncertainty_threshold,
    "Unknown",
)

print(f"Percentage of unknown per level, with uncertainty_threshold={uncertainty_threshold}:")
print(f"Level 4: {np.round(sum(combined_emb.obs[f'ann_level_4_transferred_label'] =='Unknown')/query_data.n_obs*100,2)}%")

sc.pl.umap(
    combined_emb,
    color="ann_level_4_transfer_uncert",
    ncols=3,
    frameon=False,
    save="_MouseAtlas_and_query_lab_transf_uncert.png"
)


combined_emb.obs["ann_level_4_transferred_label"] = pd.Categorical(
    pd.Series(combined_emb.obs["ann_level_4_transferred_label"].tolist()).fillna("nan")
)
if "ann_level_4_transferred_label_colors" in combined_emb.uns.keys():
    print("removing colors for level",4)
    del combined_emb.uns["ann_level_4_transferred_label_colors"]


sc.pl.umap(query_data_full, color=["ann_level_4_transferred_label_unfiltered"],  na_color="grey", ncols=1,size=0.5, save="_MouseAtlas_and_query_transf_labels.png")

sc.pl.umap(combined_emb, color="ann_level_4_transferred_label", save="_MouseAtlas_ref_ct_annotations.png")


combined_emb.write(f"{output_dir}/{output_file}")
query=combined_emb[combined_emb.obs['ref_or_query']=='query'].copy()
query_data_full.obs['ann_level_4_transferred_label']=query.obs['ann_level_4_transferred_label'].copy()
query_data_full.obs['ann_level_4_transferred_label_unfiltered']=query.obs['ann_level_4_transferred_label_unfiltered'].copy()
query_data_full.obsm=query.obsm.copy()
query_data_full.uns=query.obsm.copy()




min_n_cells = 20
### loop through batches
for batch in query_batches:
    emb_batch = combined_emb[combined_emb.obs[batch_variable] == batch,:].copy()
    # initiate empty dataframe:
    ds_uncert_stats = pd.DataFrame(index=range(len(emb_batch.obs.ann_level_5_transferred_label_unfiltered.unique())))
    # loop through levels, and store n cells, as well as 25th, 50th and 75th percentile uncertainty for every transfered label:
    for lev in range(1,6):
        lev_cts = emb_batch.obs[f"ann_level_{lev}_transferred_label_unfiltered"].unique()
        # group by transfered label
        g = emb_batch.obs.groupby(f"ann_level_{lev}_transferred_label_unfiltered")
        # count number of cells per label
        lev_cts_n = g.agg({f"ann_level_{lev}_transferred_label_unfiltered":"count"}).loc[lev_cts,:]
        # remove all labels with fewer than min_n_cells
        lev_cts_filt = lev_cts_n.index[lev_cts_n[f"ann_level_{lev}_transferred_label_unfiltered"]>=min_n_cells].tolist()
        # now store relevant stats
        ds_uncert_stats.loc[range(0,len(lev_cts_filt)),f"Level_{lev}_label"] = lev_cts_filt
        for perc in [25,50,75]:
            ds_uncert_stats.loc[range(0,len(lev_cts_filt)),f"Level_{lev}_uncert_{perc}"] = g[f"ann_level_{lev}_transfer_uncert"].apply(lambda x: np.percentile(x, perc)).loc[lev_cts_filt].values
    # remove empty rows
    ds_uncert_stats.dropna(axis=0,how="all",inplace=True)
    # round to three decimals:
    ds_uncert_stats = round(ds_uncert_stats,3)
    ds_uncert_stats.to_csv(f"{output_dir}/summary_stats_{batch}.csv")
