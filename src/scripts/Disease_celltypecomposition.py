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

all_diseases = full.obs['Disease'][full.obs['Disease'] != 'Control'].unique()
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

plots_dis = pd.DataFrame()
pct_df = pd.DataFrame()
# Iterate over all diseases
for target_disease in all_diseases:
    # Find the studies associated with the target disease
    associated_studies = full.obs['study'][full.obs['Disease'] == target_disease].unique()
    
    # Create a new observation variable 'Group' to represent control and disease groups
    full.obs[target_disease + '_Group'] = ''
    
    # Iterate over the associated studies
    for study in associated_studies:
        # Find the indices of cells belonging to the target disease and the current study
        indices_disease = (full.obs['Disease'] == target_disease) & (full.obs['study'] == study)
        # Find the indices of cells belonging to the control group and the current study
        indices_control = (full.obs['MLCA_or_query'] == 'Control') & (full.obs['study'] == study)
        full.obs.loc[indices_disease, target_disease + '_Group'] = target_disease
        full.obs.loc[indices_control, target_disease + '_Group'] = 'Control' +'_' + target_disease
    
    group = target_disease + '_Group'
    filtered_data = full.obs[(full.obs[group] != '')]
    filtered_data = filtered_data[[group,'Diff_Exp']]
    current_iteration = pd.concat([filtered_data[group], filtered_data['Diff_Exp']], axis=1)
    
    
    freq_df = current_iteration.groupby([group, 'Diff_Exp']).size().unstack()
    
    current_pct_df = freq_df.divide(freq_df.sum(axis=1), axis=0)
    current_pct_df.columns.name = None    
    
    pct_df = pd.concat([current_pct_df,pct_df])
    
    plots_dis = pd.concat([pct_df, current_iteration], ignore_index=True)


#freq_df = freq_df.sort_index(key=lambda x: x.map({v: i for i, v in enumerate(custom_order)}))
custom_order = ['Control_Severe Asthma', 'Severe Asthma', 
       'Control_Chlamydia', 'Chlamydia','Control_Old mice', 'Old mice', 
       'Control_Bleo-young', 'Bleo-young', 'Control_Bleo-old', 'Bleo-old', 
       'Control_Influenza', 'Influenza', 
       'Control_Cigarette Smoke', 'Cigarette Smoke',
       'Control_Post Sendai virus', 'Post Sendai virus', 
       'Control_TB', 'TB',
       'Control_Copd_Covid','Copd_Covid', 
       'Control_Murid herpesvirus 4', 'Murid herpesvirus 4',
       'Control_Age - 24m', 'Age - 24m',  'Control_Asthma','Asthma',
       'Control_Covid', 'Covid', 'Control_Copd', 'Copd', 'Control_Fibrosis',
       'Fibrosis', 'Control_Severe Asthma w steroids',
       'Severe Asthma w steroids', 
       'Control_Asthma w steroids', 'Asthma w steroids',
       'Control_Cancer','Cancer', 
       'Control_Hyperoxia', 'Hyperoxia', 
       'Control_Cancer tumor', 'Cancer tumor',
       'Control_Radiation', 'Radiation',
       'Control_Hypoxia', 'Hypoxia']

fig, ax2 = plt.subplots(figsize=(35, 15))  # Increase the figure size
pct_df = pct_df.sort_index(key=lambda x: x.map({v: i for i, v in enumerate(custom_order)}))
ax2 = pct_df.plot(kind="bar", stacked=True, ax=ax2)
ax2.legend(fontsize=25, ncol=2)
ax2.legend_.set_bbox_to_anchor((1, 0.8))
plt.tight_layout(pad=2.0)
ax2.set_xticklabels(ax2.get_xticklabels(), fontsize=35)
ax2.set_yticklabels(ax2.get_yticklabels(),fontsize=35)  

plt.xlabel('')
plt.ylabel('')

plt.savefig('/shared/homes/155227/work/R/plots/MouseAtlas/Disease_Cell_Type_Comp_2602.png', bbox_inches='tight')  

