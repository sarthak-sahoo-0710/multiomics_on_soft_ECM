#%%
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
from gseapy import Biomart
import os
from gseapy import barplot, dotplot
import re

# %%
names = gp.get_library_name()
print(names)

# %%
for i in names:
    print(i)
# %%
#give an input list of genes in the significant modules and specify the gene sets from the libraries available
modules_dir = "E:\\tanmay_lele_19thmay\\wgcna_23rdmay\\significant_modules_anova"
module_gene_map = pd.read_csv(os.path.join(modules_dir,"module_colors_for_data_sig.txt"), sep='\t')
#%%
degs_dir = "E:\\tanmay_lele_19thmay\\wgcna_23rdmay\\degs_overlap_module_genes"
so_sel = pd.read_csv(os.path.join(degs_dir,"so_sel_degs_p.05_with_module_info.txt"), sep="\t")
sel_anc = pd.read_csv(os.path.join(degs_dir,"sel_anc_degs_p.05_with_module_info.txt"), sep="\t")
so_st = pd.read_csv(os.path.join(degs_dir,"so_st_degs_p.05_with_module_info.txt"), sep="\t")

# %%
#list of genes in the brown module that are downregulated and show decreased accessibility in the soft selected samples
brown_down = set(pd.read_csv("E:\\tanmay_lele_19thmay\\ATAC\\bw_average_accessibility_plots\\TL_brown_down.txt", sep="\t", names=['Gene'])['Gene'].to_list())

#%%
'''module_gene_dict = {
    'black': module_gene_map[module_gene_map['Module'] == 'black']['Gene'].to_list(),
    'brown': module_gene_map[module_gene_map['Module'] == 'brown']['Gene'].to_list(),
    'blue': module_gene_map[module_gene_map['Module'] == 'blue']['Gene'].to_list(),
    'tan': module_gene_map[module_gene_map['Module'] == 'tan']['Gene'].to_list()
}'''  #dictionary of module names and their corresponding genes identified by wgcna

'''module_gene_dict = {
    'so_sel_up': list(set(so_sel[so_sel['Expression'] == "up"]["Gene symbol"])),
    'so_sel_down': list(set(so_sel[so_sel['Expression'] == "down"]["Gene symbol"])),
    'sel_anc_up': list(set(sel_anc[sel_anc['Expression'] == "up"]["Gene symbol"])),
    'sel_anc_down': list(set(sel_anc[sel_anc['Expression'] == "down"]["Gene symbol"])),
    'so_st_up': list(set(so_st[so_st['Expression'] == "up"]["Gene symbol"])),
    'so_st_down': list(set(so_st[so_st['Expression'] == "down"]["Gene symbol"]))
}''' #dictionary of genes identified by deg analysis
 
module_gene_dict = {
    'brown': list(brown_down)
}

# %%
#gene sets to be used for enrichment analysis

#gene_sets = ['GO_Biological_Process_2025','MSigDB_Hallmark_2020','KEGG_2021_Human','GO_Cellular_Component_2025','GO_Molecular_Function_2025'] 
gene_sets = ['GO_Biological_Process_2025','KEGG_2021_Human','MSigDB_Hallmark_2020']

#%%

all_enrichment_results = {}
for module_name, module_genes in module_gene_dict.items(): 
    all_enrichment_results[module_name] = {}
    print(f"Running enrichment for {module_name}...")
    enr = gp.enrichr(gene_list=module_genes, gene_sets=gene_sets, organism='Human')
    all_enrichment_results[module_name] = enr.results
#%%
#combine all the enrichment results in a dataframe 

combined_results = []
for module, df in all_enrichment_results.items():
    df = df.copy()
    df['Module'] = module #adds a column for module
    df['-log10_padj'] = -np.log10(df['Adjusted P-value'])
    combined_results.append(df)
enrichment_df = pd.concat(combined_results, ignore_index=True)
# %%
#enrichment_df.to_csv("enrichr_all_modules_genesets.txt", sep="\t", index=False, header=True)
#enrichment_df.to_csv("enrichr_on_degs_23_09_25.txt", sep="\t", index=False, header=True)

#%%
plt.rcParams['font.family'] = 'Arial'
plot_genesets = ['GO_Biological_Process_2025']
for m in ['brown']: 
#for m in ['so_sel_up','so_sel_down','sel_anc_up','sel_anc_down','so_st_up','so_st_down']:  
    print(f'Plotting for module {m}')
    try:
        df_filtered = all_enrichment_results[m]
        df_filtered = df_filtered[df_filtered['Gene_set'].isin(set(plot_genesets))]
        fig = barplot(df_filtered,
                    column="Adjusted P-value",
                    group='Gene_set',
                    size=25,
                    figsize=(3,15),
                    color={'GO_Biological_Process_2025':'teal','KEGG_2021_Human':'blue','MSigDB_Hallmark_2020':'red'},top_term=15)
        fig.set_title(f'{m}', fontsize=20)
        
        fig.set_xlabel(fig.get_xlabel(), fontsize=16, fontweight="bold", fontfamily="Arial")
        fig.set_ylabel(fig.get_ylabel(), fontsize=16, fontweight="bold", fontfamily="Arial")
        
        # Tick labels
        fig.tick_params(axis='x', labelsize=14)
        fig.tick_params(axis='y', labelsize=20)
        for label in fig.get_xticklabels() + fig.get_yticklabels():
            label.set_fontweight("bold")
            label.set_fontfamily("Arial")
        #fig.figure.savefig(f"{m}_down_enrichment_msigdb.png", dpi=500, bbox_inches='tight')
    except ValueError as e:
        print(f"Skipped module {m}: {e}")
        continue

