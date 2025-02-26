

full.obs['Decon'] = full.obs['Level_2_5']
full.obs['Decon'] = full.obs['Decon'].astype(object)

full.obs['Decon'].loc[(full.obs['Level_2_5']== 'Megakaryocyte')] = 'Granulocyte'
full.obs['Decon'].loc[(full.obs['Level_2_5']== 'SMC population')\
	] = 'Fibroblast'
full.obs['Decon'].loc[(full.obs['Level_2_5']== 'Multiciliated_Deuterostome')\
	] = 'Multiciliated_lineage'

full.obs.Decon.value_counts()


full.X.data = full.layers['raw'].data.copy()
decon = full[full.obs['MLCA_or_query'] == 'Control']
counts = decon.layers['raw']
decon.X.data = counts.data.copy()
batches = np.unique(decon.obs['Decon'])
batches=batches[batches != 'Other cells']
batches = batches.astype(str)
cell_indices = []
for batch in batches:
	idx = np.where(decon.obs['Decon'] == batch)[0]
	np.random.shuffle(idx)
	cell_indices += idx[:200].tolist()

adata_filter = decon[cell_indices,:]

matrix = adata_filter.X.toarray()

mat = pd.DataFrame(data = matrix, index = adata_filter.obs["Decon"],\
	columns = adata_filter.var_names).T

mat.to_csv('Change to attempt NUMBER.tsv',sep = '\t', index_label = 'GeneSymbol')

