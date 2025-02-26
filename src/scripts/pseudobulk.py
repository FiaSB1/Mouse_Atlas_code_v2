import os
#os.chdir('/data/kkovacs/Python')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
#import sankey
import anndata
import pandas as pd
from pySankey.sankey import sankey
import plotly.graph_objects as go
from adpbulk import ADPBulk



full =sc.read("/shared/ci/kyle/MLCA_All_meta_done.h5ad")

def replace_obs_value(adata, old_value, new_value):
    """
    Replace a specific value in all observation variables with a new one.
    
    Parameters:
        adata (anndata.AnnData): An AnnData object.
        old_value: The value to replace.
        new_value: The value to replace with.
    """
    for variable_name in adata.obs.columns:
        adata.obs[variable_name] = adata.obs[variable_name].replace(old_value, new_value)

replace_obs_value(full, 'Multiciliated_Deuterostome', 'Multiciliated_Deuterosomal')


#full = full[~full.obs['Level_5'].isin(['Doublets'])]

#full.obs['sex']=full.obs['sex'].replace(0.0,'Female')
#full.obs['sex']=full.obs['sex'].replace(1.0,'Male')


sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(15, 15))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=12)

full.obs['Diff_Exp'] = full.obs.Level_4.copy()
full.obs['Diff_Exp'] = full.obs['Diff_Exp'].astype(object)

full.obs['Diff_Exp'].loc[full.obs['Diff_Exp'].isin(['Classical monocyte','Non-classical monocyte'])] = 'Monocytes'
full.obs['Diff_Exp'].loc[full.obs['Diff_Exp'].isin(['Migratory dendritic','pDC'])] = 'Dendritic cells'
full.obs['Diff_Exp'].loc[full.obs['Diff_Exp'].isin(['Pericyte','Smooth muscle cells','Mesothelium'])] = 'Fibroblast'
full.obs['Diff_Exp'].loc[full.obs['Diff_Exp'].isin(['IgA plasma','IgM plasma', 'Proliferating plasma'])] = 'Plasma cells'
full.obs['Diff_Exp'].loc[full.obs['Diff_Exp'].isin(['Venous endothelial','General capillary','Aerocyte capillary',
													'Arterial endothelial','Lymphatic endothelial','Progenitor endothelial',
													])] = 'Endothelial cells'
full.obs['Diff_Exp'].loc[full.obs['Diff_Exp'].isin(['Prg4+ Macrophages'])] = 'Interstitial macrophages'
full.obs['Diff_Exp'].loc[full.obs['Diff_Exp'].isin(['Neutrophils'])] = 'Granulocyte'
full.obs['Diff_Exp'].loc[full.obs['Diff_Exp'].isin(['Multiciliated_Deuterostome'])] = 'Multiciliated_Deuterosomal'



# full.obs['Disease_groups']=full.obs.Disease
# full.obs['Disease_groups'] = full.obs['Disease_groups'].astype(object)
# full.obs['MLCA_or_query'] = full.obs['MLCA_or_query'].astype(object)
# full.obs['study'] = full.obs['study'].astype(object)

# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Murid herpesvirus 4','Covid','Influenza','Post Sendai virus','Copd_Covid'])] = 'Viral'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Chlamydia','Tuberculosis'])] = 'Bacteria'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Asthma','Asthma w steroids','Severe Asthma','Severe Asthma w steroids'])] = 'Asthma'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Cancer','Cancer tumor'])] = 'Cancer'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Age - 24m'])] = 'Old mice'


# con=full[full.obs['Disease_groups']=='Control'].copy()
# dis=full[~full.obs['Disease_groups'].isin(['Control'])].copy()
# con.obs['Disease_groups']=con.obs['study']+'_'+con.obs['Disease_groups']
# con.obs['Disease_groups']=con.obs['Disease_groups'].str.replace('Mouse_','')

# full=con.concatenate(dis)

# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['TB_Control',''])] = 'Control Bacteria'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Influenza_Control','Covid_Control'])] = 'Control Viral'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['TumorKO_Control','Cancer_Control'])] = 'Control Cancer'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Murid_Asthma_Control','asthma_Control'])] = 'Control Asthma'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Copd_Control'])] = 'Control COPD'
# full.obs['Disease_groups'].loc[full.obs['Disease_groups'].isin(['Covid_Control'])] = 'Control Covid'




full.obs['Pseudo'] = full.obs['sample'] #adata.obs.Pseudo = mice samples
full.obs['Pseudo'] = full.obs['Pseudo'].astype(object)
full.obs['Pseudo'].loc[full.obs['study'] == 'Mouse_Post_Sendai'] = 'Removed'

full.obs['combined_category']=full.obs['Diff_Exp'].astype(str)+'_'+full.obs['Pseudo'].astype(str)
full.obs['combined_category'].loc[full.obs['study'] == 'Mouse_Post_Sendai'] = 'Removed'

#Pseudobulk per sample
adpb = ADPBulk(full, "Pseudo")

# perform the pseudobulking
pseudo = adpb.fit_transform().T

pseudo = pseudo.rename(columns=lambda x: x.replace('Pseudo.',''))

pseudo = pseudo.drop('Removed',axis=1)
pseudo.to_csv('/shared/homes/155227/work/R/data/MouseAtlas/Pseudo_Mouse_Atlas.tsv',sep = '\t', index_label='GeneSymbol')  
#pseudo = pd.read_csv('/data/kkovacs/Python/Pseudo_Mouse_Atlas.tsv',sep = '\t')


obs_df = pd.DataFrame({
    'Pseudo': full.obs['Pseudo'],
    'age': full.obs['age'],
    'sex': full.obs['sex'],
    'platform': full.obs['platform'],
    'Viral': full.obs['Viral'],
    'strain': full.obs['strain'],
    'Disease': full.obs['Disease'],
    'Group': full.obs['Group'],
    'study': full.obs['study'],
	'Severe Asthma':full.obs['Severe Asthma_Group'],
	'Chlamydia':full.obs['Chlamydia_Group'],
	'Oldmice': full.obs['Old mice_Group'],
	'Bleoyoung':full.obs['Bleo-young_Group'],
	'Bleoold':full.obs['Bleo-old_Group'],
	'Influenza':full.obs['Influenza_Group'],
	'Cigarette Smoke':full.obs['Cigarette Smoke_Group'],
	'Postsendai':full.obs['Post Sendai virus_Group'],
	'Age':full.obs['Age - 24m_Group'],
	'TB':full.obs['TB_Group'],
	'CopdCovid':full.obs['Copd_Covid_Group'],
	'Herpes':full.obs['Murid herpesvirus 4_Group'],
	'Asthma':full.obs['Asthma_Group'],
	'Covid':full.obs['Covid_Group'],
	'Copd':full.obs['Copd_Group'],
	'Fibrosis':full.obs['Fibrosis_Group'],
	'Mir155KO':full.obs['Mir-155 KO_Group'],
	'Severe Asthma w Steroids':full.obs['Severe Asthma w steroids_Group'],
	'Asthma w Steroids':full.obs['Asthma w steroids_Group'],
	'Cancer':full.obs['Cancer_Group'],
	'Hyperoxia':full.obs['Hyperoxia_Group'],
	'CancerTumor':full.obs['Cancer tumor_Group'],
	'Radiation':full.obs['Radiation_Group'],
	'Hypoxia':full.obs['Hypoxia_Group']
    #'Diff_expre_level': full.obs['Diff_Expre'],
    # Add other observation columns as needed
})

# Use pivot_table to reshape the DataFrame
table = pd.pivot_table(obs_df, values=['age', 'sex','platform','Viral','strain','Disease','Group','study',
									   'Severe Asthma','Chlamydia','Oldmice','Bleoyoung','Bleoold','Influenza',
									   'Cigarette Smoke','Postsendai','Age','TB','CopdCovid','Herpes',
									   'Asthma','Covid','Copd','Fibrosis','Mir155KO','Severe Asthma w Steroids',
									   'Asthma w Steroids','Cancer','Hyperoxia','CancerTumor','Radiation',
									   'Hypoxia'], index=['Pseudo'], aggfunc='first')

table = table.drop('Removed',axis=0)

table.to_csv('/shared/homes/155227/work/R/data/MouseAtlas/Meta_Pseudo.csv')







from adpbulk import ADPBulk
adpb = ADPBulk(full, "combined_category")

# perform the pseudobulking
pseudo = adpb.fit_transform().T

pseudo = pseudo.rename(columns=lambda x: x.replace('combined_category.',''))

pseudo = pseudo.drop('Removed',axis=1)
pseudo.to_csv('/shared/homes/155227/work/R/data/MouseAtlas/pseudobulk_percelltypes.tsv',sep = '\t', index_label='GeneSymbol')  



obs_df = pd.DataFrame({
	'Pseudo': full.obs['Pseudo'],
	'age': full.obs['age'],
	'sex': full.obs['sex'],
	'platform': full.obs['platform'],
	'Viral': full.obs['Viral'],
	'strain': full.obs['strain'],
	'Disease': full.obs['Disease'],
	'Group': full.obs['Group'],
	'study': full.obs['study'],
	'Diff_exp': full.obs['Diff_Exp'],
	'combined_category':full.obs['combined_category'],
	'Severe Asthma':full.obs['Severe Asthma_Group'],
	'Chlamydia':full.obs['Chlamydia_Group'],
	'Oldmice': full.obs['Old mice_Group'],
	'Bleoyoung':full.obs['Bleo-young_Group'],
	'Bleoold':full.obs['Bleo-old_Group'],
	'Influenza':full.obs['Influenza_Group'],
	'Cigarette Smoke':full.obs['Cigarette Smoke_Group'],
	'Postsendai':full.obs['Post Sendai virus_Group'],
	'Age':full.obs['Age - 24m_Group'],
	'TB':full.obs['TB_Group'],
	'CopdCovid':full.obs['Copd_Covid_Group'],
	'Herpes':full.obs['Murid herpesvirus 4_Group'],
	'Asthma':full.obs['Asthma_Group'],
	'Covid':full.obs['Covid_Group'],
	'Copd':full.obs['Copd_Group'],
	'Fibrosis':full.obs['Fibrosis_Group'],
	'Mir155KO':full.obs['Mir-155 KO_Group'],
	'Severe Asthma w Steroids':full.obs['Severe Asthma w steroids_Group'],
	'Asthma w Steroids':full.obs['Asthma w steroids_Group'],
	'Cancer':full.obs['Cancer_Group'],
	'Hyperoxia':full.obs['Hyperoxia_Group'],
	'CancerTumor':full.obs['Cancer tumor_Group'],
	'Radiation':full.obs['Radiation_Group'],
	'Hypoxia':full.obs['Hypoxia_Group']
})

table = pd.pivot_table(obs_df, values = ['age', 'sex','platform','Viral','strain','Disease','Group','study',
									   'Severe Asthma','Chlamydia','Oldmice','Bleoyoung','Bleoold','Influenza',
									   'Cigarette Smoke','Postsendai','Age','TB','CopdCovid','Herpes',
									   'Asthma','Covid','Copd','Fibrosis','Mir155KO','Severe Asthma w Steroids',
									   'Asthma w Steroids','Cancer','Hyperoxia','CancerTumor','Radiation',
									   'Hypoxia'], index = ['combined_category'], aggfunc = 'first')

table.to_csv("/shared/homes/155227/work/R/data/MouseAtlas/meta_pseudobulk_percelltypes.csv")



#full.obs['final_anno'] = full.obs['Level_5']
#dis = full[full.obs['MLCA_or_query'] == 'Disease']
con = full[full.obs['MLCA_or_query'] == 'Control']

plots_dis = pd.concat([con.obs['age'],con.obs['age_cont'],con.obs['strain'],con.obs['sex'],con.obs['study'],con.obs['platform'],con.obs['Diff_Exp'],con.obs['Pseudo']],axis=1)
cat_dtype = pd.CategoricalDtype(categories=['Young','Adult','Old'], ordered=True)
plots_dis['age_cont'] = plots_dis['age_cont'].astype(cat_dtype)
plots_dis['age_cont']=plots_dis['age_cont'].iloc['Young','Adult','Old']

plots_dis['age'] = plots_dis['age'].astype('category')

#plots_dis = plots_dis[~plots_dis['age'].isin(['5-6'])]

plots_dis.age=plots_dis.age.cat.remove_categories('5-6')
#cross_tab = plots_dis.crosstab(index=[plots_dis['combined_category'], plots_dis['age']], columns=plots_dis['Diff_Exp'])



# ax = plots_dis.groupby(['strain','Level_3']).size().unstack().\
# 	plot(kind='bar', stacked=True,figsize=(20,20))
# ax.legend(fontsize = 20, ncol =2)
# ax.legend_.set_bbox_to_anchor((1,1))
# plt.tight_layout()
# plt.xticks(fontsize = 15)
# plt.yticks(fontsize = 15)
# plt.savefig('/shared/homes/155227/work/R/plots/MouseAtlas/stackedbarplot_age.png')

meta=plots_dis.columns[0:len(plots_dis.columns)-2]
for m in meta:
	freq_df = plots_dis.groupby([m,'Diff_Exp']).size().unstack()
	pct_df = freq_df.divide(freq_df.sum(axis=1), axis=0)
	ax2 = pct_df.plot(kind="bar", stacked=True,figsize=(25,10))
	ax2.legend(fontsize = 20,ncol=2)
	ax2.legend_.set_bbox_to_anchor((1,0.8))
	plt.tight_layout()
	plt.xticks(fontsize = 24,rotation=35)
	plt.yticks(fontsize = 24) 
	plt.savefig('/shared/homes/155227/work/R/plots/MouseAtlas/stackedbarplot/stackedbarplot_'+ str(m)+ '_pct.png')
	result = plots_dis.groupby(['Pseudo','Diff_Exp',m]).size().unstack(fill_value=0)
	result=result.reset_index()
	total=plots_dis.Pseudo.value_counts()
	total=total.to_frame(0) #convert to dataframe with 0 as title because otherwise it's as a series
	total['Pseudo']=total.index.values
	merged = pd.merge(result, total, on='Pseudo', how='inner')
	#for i in range(2,len(merged.columns)-1):
	 #	merged[merged.columns[i]] = merged[merged.columns[i]] / merged[0]
	# test=pd.concat([plots_dis.Pseudo,plots_dis[m]],axis=1)
	# merged2=pd.merge(merged, test, on='Pseudo', how='inner')
	merged.to_csv('/shared/homes/155227/work/R/data/MouseAtlas/celltypeproportion/'+str(m)+'_celltypeproportion.csv')



# dis=pd.concat([full.obs['Disease_groups'],full.obs['Diff_Exp']],axis=1)

# meta=dis.columns[0:len(dis.columns)-1]
# plots_dis=dis
# for m in meta:
# 	freq_df = plots_dis.groupby([m,'Diff_Exp']).size().unstack()
# 	pct_df = freq_df.divide(freq_df.sum(axis=1), axis=0)
# 	ax2 = pct_df.plot(kind="bar", stacked=True,figsize=(30,10))
# 	ax2.legend(fontsize = 20,ncol=2)
# 	ax2.legend_.set_bbox_to_anchor((1,0.5))
# 	plt.tight_layout()
# 	plt.xticks(fontsize = 24,rotation=90)
# 	plt.yticks(fontsize = 24) 
# 	plt.savefig('/shared/homes/155227/work/R/plots/MouseAtlas/stackedbarplot_'+ str(m)+ '_pct.png')

plots_con = pd.concat([con.obs['study'],con.obs['Level_5']],axis=1)

macrophages=['Alveolar macrophages',
			 'Interstitial macrophages',
			 'Monocyte',
			 'Viral induced macrophages']
macro=full[full.obs['Level_3'].isin(macrophages)].copy()


