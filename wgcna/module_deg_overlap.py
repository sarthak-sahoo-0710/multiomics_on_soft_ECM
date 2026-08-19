#%%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
from matplotlib_venn import venn3, venn3_circles
# %%
modules_dir = "E:\\tanmay_lele_19thmay\\wgcna_23rdmay\\significant_modules_anova"
module_gene_map = pd.read_csv(os.path.join(modules_dir,"module_colors_for_data_sig.txt"), sep='\t')
black_sel = set(module_gene_map[module_gene_map['Module'] == 'black']['Gene'].to_list()) #selected module
brown_anc = set(module_gene_map[module_gene_map['Module'] == 'brown']['Gene'].to_list()) #ancestral module; look for the negative regulators of the module and remove them
blue_so_sel = set(module_gene_map[module_gene_map['Module'] == 'blue']['Gene'].to_list()) #soft selected module
tan_st = set(module_gene_map[module_gene_map['Module'] == 'tan']['Gene'].to_list()) #stiff module

#%%
gene_exp_df = pd.read_csv('D:\\sejal\\tanmay_lele\\rpkm_rna.txt', sep='\t')

soft_selected = [col for col in gene_exp_df.columns if col.startswith('soft_selected')]
stiff_selected = [col for col in gene_exp_df.columns if col.startswith('stiff_selected')]
soft_ancestral = [col for col in gene_exp_df.columns if col.startswith('soft_ancestral')]
stiff_ancestral = [col for col in gene_exp_df.columns if col.startswith('stiff_ancestral')]

selected_samples = soft_selected + stiff_selected
ancestral_samples = soft_ancestral + stiff_ancestral
cols = soft_selected + stiff_selected + soft_ancestral + stiff_ancestral

log2_gene_exp_df = gene_exp_df.copy()
log2_gene_exp_df[cols] = log2_gene_exp_df[cols].apply(lambda x: np.log2(x+1))

#%%
def find_the_degs_between_pairs(group1,group2):
    t_stats = []
    p_vals = []
    significant_genes =[]
    expression_status =[]
    for idx, row in log2_gene_exp_df.iterrows():
        group1_list = row[group1].dropna().tolist()
        group2_list = row[group2].dropna().tolist()
        t_stat, p_val = ttest_ind(group1_list,group2_list, equal_var=True, alternative='two-sided')
        #if p_val < alpha:
        significant_genes.append(f"{row['symbol']}")
        p_vals.append(p_val)
        t_stats.append(t_stat)
        mean_group1 = np.mean(group1_list)
        mean_group2 = np.mean(group2_list)
        expression_status.append("up" if mean_group1 > mean_group2 else "down")
    results_df = pd.DataFrame({'Gene symbol': significant_genes,'p_val': p_vals,'t_stat': t_stats, 'Expression':expression_status})
    results_df_filt = results_df.dropna(subset=['p_val'])
    results_df_sig = results_df_filt[results_df_filt['p_val'] < 0.05]
    results_df_mod = pd.merge(results_df_sig,module_gene_map,left_on='Gene symbol',right_on='Gene')
    return results_df_mod
# %%
sel_anc = find_the_degs_between_pairs(selected_samples,ancestral_samples)
so_st = find_the_degs_between_pairs(soft_selected+soft_ancestral,stiff_selected+stiff_ancestral)
so_sel = find_the_degs_between_pairs(soft_selected,stiff_ancestral+soft_ancestral+stiff_selected)
#%%
sel_anc.to_csv("sel_anc_degs_p.05_with_module_info.txt",sep="\t",index=False)
so_st.to_csv("so_st_degs_p.05_with_module_info.txt",sep="\t",index=False)
so_sel.to_csv("so_sel_degs_p.05_with_module_info.txt",sep="\t",index=False)
#%%
def plot_module_distribution(degs_df, filename):
    counts = degs_df.groupby(["Expression", "Module"]).size().reset_index(name="count") 
    #make a percentage column and put num of module genes divided by total upreg or downreg genes
    totals = counts.groupby("Expression")["count"].transform("sum")

    counts["percentage"] = (counts["count"] / totals) * 100

    # Pivot to wide format for stacked bar plotting
    counts_pivot = counts.pivot(index="Expression", columns="Module", values="count").fillna(0)
    module_palette = {"blue":"#0000ffff","black":"#000000ff","brown":"#a52a2aff","tan":"#d2b48cff"}
    # Plot stacked bar chart
    counts_pivot.plot(
        kind="bar", stacked=True, figsize=(3.5, 6), 
        color=[module_palette[m] for m in counts_pivot.columns], width=0.6,legend=False
    )

    plt.xlabel("Expression direction")
    plt.ylabel("Number of genes")
    #plt.title(title)
    #plt.legend(title="Module")
    plt.tight_layout()
    #plt.savefig(f"{filename}.png",dpi=500,bbox_inches='tight')
    counts.to_csv(f"{filename}.txt",sep="\t",index=False)
    plt.show()
    return counts
#%%
plot_module_distribution(sel_anc, "sel_anc_module")
plot_module_distribution(so_st, "so_st_module")
plot_module_distribution(so_sel, "so_sel_module")

# %%
#venn diagrams for paper revision
sel_anc = pd.read_csv("./degs_overlap_module_genes/sel_anc_degs_p.05_with_module_info.txt",sep="\t")
so_st = pd.read_csv("./degs_overlap_module_genes/so_st_degs_p.05_with_module_info.txt",sep="\t")
# %%
#filter both the dataframes for expression direction up 
sel_anc_up = sel_anc[sel_anc['Expression'] == 'up']
so_st_up = so_st[so_st['Expression'] == 'up']
# %%
gene_col = "Gene symbol"
set_sel_anc = set(sel_anc_up[gene_col].dropna())
set_so_st = set(so_st_up[gene_col].dropna())
set_blue_module = blue_so_sel  # Already a set
# %%

# 4. Setup Figure
fig, ax = plt.subplots(figsize=(8, 8))
custom_hex_colors = ["#3f464b", "#faa165", "#0000ffff"]
# 5. Draw Venn Diagram with Seaborn palette
venn = venn3(
    subsets=(set_sel_anc, set_so_st, set_blue_module),
    set_colors=custom_hex_colors,
    alpha=0.7,
    ax=ax,
)

# 6. Add clean, thin outer circles using Seaborn's dark gray palette tone
venn3_circles(
    subsets=(set_sel_anc, set_so_st, set_blue_module),
    linewidth=1.2,
    color="#99a3aa",  # Use the first custom color for the outer circles
    ax=ax,
)


# Styling set labels to match Seaborn's dark text color
dark_color = sns.color_palette("dark")[7]

for text in venn.subset_labels:
    if text:
        text.set_fontsize(15)
        text.set_fontweight("bold")
        text.set_color("red")

# Remove axes spines for a clean Seaborn aesthetic
sns.despine(left=True, bottom=True)
#plt.savefig("venn_sel_anc_so_st_blue_module.png", dpi=600, bbox_inches='tight')
plt.show()


# %%
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3, venn3_circles

# -------------------------------------------------------------------------
# 1. PREPARE SETS & CALCULATE BLUE MODULE OVERLAPS
# -------------------------------------------------------------------------
gene_col = "Gene symbol"
set_sel_anc = set(sel_anc_up[gene_col].dropna())
set_so_st = set(so_st_up[gene_col].dropna())
set_blue_module = blue_so_sel  # Already a set

total_blue = len(set_blue_module)

# Overlap subset calculations (for percentage summary)
blue_in_sel_anc = len(set_blue_module & set_sel_anc)
blue_in_so_st = len(set_blue_module & set_so_st)
blue_in_both = len(set_blue_module & set_sel_anc & set_so_st)
blue_in_either = len(set_blue_module & (set_sel_anc | set_so_st))

# Percentages relative to total blue module size
pct_sel_anc = (blue_in_sel_anc / total_blue) * 100 if total_blue > 0 else 0
pct_so_st = (blue_in_so_st / total_blue) * 100 if total_blue > 0 else 0
pct_both = (blue_in_both / total_blue) * 100 if total_blue > 0 else 0
pct_either = (blue_in_either / total_blue) * 100 if total_blue > 0 else 0

# -------------------------------------------------------------------------
# 2. CALCULATE EXACT 7 VENN REGIONS USING PURE PYTHON SET MATH
# -------------------------------------------------------------------------
c_100 = len(set_sel_anc - set_so_st - set_blue_module)       # Only SEL ANC
c_010 = len(set_so_st - set_sel_anc - set_blue_module)       # Only SO ST
c_001 = len(set_blue_module - set_sel_anc - set_so_st)       # Only Blue Module
c_110 = len((set_sel_anc & set_so_st) - set_blue_module)      # SEL ANC & SO ST only
c_101 = len((set_sel_anc & set_blue_module) - set_so_st)     # SEL ANC & Blue Module only
c_011 = len((set_so_st & set_blue_module) - set_sel_anc)     # SO ST & Blue Module only
c_111 = len(set_sel_anc & set_so_st & set_blue_module)       # Center (All 3)

# -------------------------------------------------------------------------
# 3. PRINT OVERLAP COUNTS & PERCENTAGES TO CONSOLE
# -------------------------------------------------------------------------
print("=" * 55)
print("1. INDIVIDUAL VENN REGION COUNTS (7 REGIONS)")
print("=" * 55)
print(f"• Only SEL ANC (100)           : {c_100}")
print(f"• Only SO ST (010)             : {c_010}")
print(f"• Only Blue Module (001)       : {c_001}")
print(f"• SEL ANC & SO ST only (110)   : {c_110}")
print(f"• SEL ANC & Blue Module (101)  : {c_101}")
print(f"• SO ST & Blue Module (011)    : {c_011}")
print(f"• Center - All 3 Sets (111)    : {c_111}")

print("\n" + "=" * 55)
print(f"2. BLUE MODULE OVERLAP SUMMARY (Total Blue = {total_blue})")
print("=" * 55)
print(f"• Overlap with SEL ANC (Up): {blue_in_sel_anc} genes ({pct_sel_anc:.2f}%)")
print(f"• Overlap with SO ST (Up):   {blue_in_so_st} genes ({pct_so_st:.2f}%)")
print(f"• In Both Up Sets:          {blue_in_both} genes ({pct_both:.2f}%)")
print(f"• Total Unique Overlap:     {blue_in_either} genes ({pct_either:.2f}%)")
print("=" * 55)

# -------------------------------------------------------------------------
# 4. GENERATE CLEAN VENN DIAGRAM (NO NUMBERS, NO LABELS)
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 8))


#custom_hex_colors = ["#4A6572", "#F9AA33", "#0000ffff"]
custom_hex_colors = ["#4A5568", "#DD9C10", "#0000FF"]
# Draw Venn Diagram without set labels or subset numbers
venn = venn3(
    subsets=(set_sel_anc, set_so_st, set_blue_module),
    set_labels=None,  # Suppresses outer set labels
    subset_label_formatter=lambda x: "",  # Suppresses numbers inside circles
    set_colors=custom_hex_colors,
    alpha=0.7,
    ax=ax,
)

# Draw boundary outlines
venn3_circles(
    subsets=(set_sel_anc, set_so_st, set_blue_module),
    linewidth=1.2,
    color="#99a3aa",
    ax=ax,
)

sns.despine(left=True, bottom=True)
plt.tight_layout()
#plt.savefig("venn_sel_anc_so_st_blue_module_clean.png", dpi=600, bbox_inches='tight')
plt.show()
