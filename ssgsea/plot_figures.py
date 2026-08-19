#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
# Create custom legend instead of colorbar
from matplotlib.patches import Patch
plt.rcParams.update({'font.size': 18, 'font.family': 'Arial'})
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

# %%
df_pid = pd.read_csv("G:/projects_part_2/tanmay_lele/ssGSEA/scores/PID_pathway_nes_scores.txt", sep="\t", index_col=0)
# get groups from column names
groups = [col.split("_Re")[0] for col in df_pid.columns]
# perform pca on the dataframe
pca = PCA(n_components=2)
pca_result = pca.fit_transform(df_pid.T)
# plot the pca result with samples colored by their condition
condition_array = ['soft_ancestral']*4 + ['soft_selected']*4 + ['stiff_ancestral']*4 + ['stiff_selected']*4
# Define color map manually
condition_list = ['soft_ancestral', 'soft_selected', 'stiff_ancestral', 'stiff_selected']
color_map = {
    'soft_selected': '#8E44AD', 
    'stiff_selected': '#27AE60', 
    'soft_ancestral': '#3498DB', 
    'stiff_ancestral': '#FF851B'}
# Map conditions to colors
colors = [color_map[cond] for cond in condition_array]
# Plot
plt.figure(figsize=(8,6))
scatter = plt.scatter(pca_result[:,0], pca_result[:,1], c=colors, s=150)
# write percentage of variance explained by each pca component on the axes
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.ylim(-22, 22)
#plt.xlim(-32, 32)
#plt.axhline(0, color='grey', lw=1)
#plt.axvline(0, color='grey', lw=1)
# Axis labels and title
plt.xlabel('PCA-1 ({:.2f}% variance)'.format(pca.explained_variance_ratio_[0]*100))
plt.ylabel('PCA-2 ({:.2f}% variance)'.format(pca.explained_variance_ratio_[1]*100))
# plot legend with custom colors
#legend_elements = [Patch(facecolor=color_map[cond], label=cond) for cond in condition_list]
#plt.legend(handles=legend_elements, loc='best', fontsize=12, frameon=False)
#plt.savefig("PCA_PID_pathways.png", dpi=500, bbox_inches='tight')
plt.show()
# %%
df_hallmark = pd.read_csv("G:/projects_part_2/tanmay_lele/ssGSEA/scores/PID_pathway_nes_scores.txt", sep="\t", index_col=0)
# get groups from column names
groups = [col.split("_Re")[0] for col in df_hallmark.columns]
# perform pca on the dataframe
pca = PCA(n_components=2)
pca_result = pca.fit_transform(df_hallmark.T)
# plot the pca result with samples colored by their condition
condition_array = ['soft_ancestral']*4 + ['soft_selected']*4 + ['stiff_ancestral']*4 + ['stiff_selected']*4
# Define color map manually
condition_list = ['soft_ancestral', 'soft_selected', 'stiff_ancestral', 'stiff_selected']
color_map = {
    'soft_selected': '#8E44AD', 
    'stiff_selected': '#27AE60', 
    'soft_ancestral': '#3498DB', 
    'stiff_ancestral': '#FF851B'}
# Map conditions to colors
colors = [color_map[cond] for cond in condition_array]
# Plot
plt.figure(figsize=(8,6))
scatter = plt.scatter(pca_result[:,0], pca_result[:,1], c=colors, s=150)
# write percentage of variance explained by each pca component on the axes
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.ylim(-22, 22)
#plt.xlim(-32, 32)
#plt.axhline(0, color='grey', lw=1)
#plt.axvline(0, color='grey', lw=1)
# Axis labels and title
plt.xlabel('PCA-1 ({:.2f}% variance)'.format(pca.explained_variance_ratio_[0]*100))
plt.ylabel('PCA-2 ({:.2f}% variance)'.format(pca.explained_variance_ratio_[1]*100))
# plot legend with custom colors
#legend_elements = [Patch(facecolor=color_map[cond], label=cond) for cond in condition_list]
#plt.legend(handles=legend_elements, loc='best', fontsize=12, frameon=False)
#plt.savefig("PCA_hallmark_pathways.png", dpi=500, bbox_inches='tight')
plt.show()

# %%
# read so_anc_st_anc_pid_pathway_ttest.txt from G:\projects_part_2\tanmay_lele\ssGSEA\pathway_t_test\PID

df_ttest_soanc_stanc = pd.read_csv("G:/projects_part_2/tanmay_lele/ssGSEA/pathway_t_test/PID/so_anc_st_anc_pid_pathway_ttest.txt", sep="\t", index_col=0)

# make a volcano plot for the ttest results
plt.figure(figsize=(8,6))
# plot all points in grey
plt.scatter(df_ttest_soanc_stanc['t_stat'], -np.log10(df_ttest_soanc_stanc['p_val']), color='grey', s=18)
# highlight points with pval < 0.05 and abs(t_stat) > 5 in red
significant = (df_ttest_soanc_stanc['p_val'] < 0.05) & ((df_ttest_soanc_stanc['t_stat']) > 5)
for pathway in df_ttest_soanc_stanc[significant].index:
    print(pathway, df_ttest_soanc_stanc.loc[pathway])
plt.scatter(df_ttest_soanc_stanc[significant]['t_stat'], -np.log10(df_ttest_soanc_stanc[significant]['p_val']), color='red', s=30)
significant = (df_ttest_soanc_stanc['p_val'] < 0.05) & ((df_ttest_soanc_stanc['t_stat']) < -5)
for pathway in df_ttest_soanc_stanc[significant].index:
    print(pathway, df_ttest_soanc_stanc.loc[pathway])
plt.scatter(df_ttest_soanc_stanc[significant]['t_stat'], -np.log10(df_ttest_soanc_stanc[significant]['p_val']), color='blue', s=30)
# add horizontal line at pval = 0.05
plt.axhline(-np.log10(0.05), color='black', linestyle='--')
# add vertical lines at t_stat = -5 and 5
plt.axvline(-5, color='black', linestyle='--')
plt.axvline(5, color='black', linestyle='--')
# add labels and title
plt.xlabel('t-statistic (Soft Ancestral vs Stiff Ancestal)')
plt.ylabel('-log10(p-value)')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(-10.1, 10.1)
#plt.ylim(0, 5)
#plt.savefig("volcano_plot_soanc_stanc_pid_pathways.png", dpi=500, bbox_inches='tight')
plt.show()
# %%
df_ttest_soanc_stanc = pd.read_csv("G:/projects_part_2/tanmay_lele/ssGSEA/pathway_t_test/PID/so_sel_so_anc_pid_pathway_ttest.txt", sep="\t", index_col=0)

# make a volcano plot for the ttest results
plt.figure(figsize=(8,6))
# plot all points in grey
plt.scatter(df_ttest_soanc_stanc['t_stat'], -np.log10(df_ttest_soanc_stanc['p_val']), color='grey', s=18)
# highlight points with pval < 0.05 and abs(t_stat) > 5 in red
significant = (df_ttest_soanc_stanc['p_val'] < 0.05) & ((df_ttest_soanc_stanc['t_stat']) > 5)
plt.scatter(df_ttest_soanc_stanc[significant]['t_stat'], -np.log10(df_ttest_soanc_stanc[significant]['p_val']), color='red', s=30)
significant = (df_ttest_soanc_stanc['p_val'] < 0.05) & ((df_ttest_soanc_stanc['t_stat']) < -5)
plt.scatter(df_ttest_soanc_stanc[significant]['t_stat'], -np.log10(df_ttest_soanc_stanc[significant]['p_val']), color='blue', s=30)
# add horizontal line at pval = 0.05
plt.axhline(-np.log10(0.05), color='black', linestyle='--')
# add vertical lines at t_stat = -5 and 5
plt.axvline(-5, color='black', linestyle='--')
plt.axvline(5, color='black', linestyle='--')
# add labels and title
plt.xlabel('t-statistic (Soft Selected vs Soft Ancestral)')
plt.ylabel('-log10(p-value)')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.xlim(-10.1, 10.1)
#plt.ylim(0, 5)
#plt.savefig("volcano_plot_sosel_soanc_pid_pathways.png", dpi=500, bbox_inches='tight')
plt.show()
#%%
df_ttest_soanc_stanc = pd.read_csv("G:/projects_part_2/tanmay_lele/ssGSEA/pathway_t_test/PID/st_sel_st_anc_pid_pathway_ttest.txt", sep="\t", index_col=0)

# make a volcano plot for the ttest results
plt.figure(figsize=(8,6))
# plot all points in grey
plt.scatter(df_ttest_soanc_stanc['t_stat'], -np.log10(df_ttest_soanc_stanc['p_val']), color='grey', s=18)
# highlight points with pval < 0.05 and abs(t_stat) > 5 in red
significant = (df_ttest_soanc_stanc['p_val'] < 0.05) & ((df_ttest_soanc_stanc['t_stat']) > 5)
plt.scatter(df_ttest_soanc_stanc[significant]['t_stat'], -np.log10(df_ttest_soanc_stanc[significant]['p_val']), color='red', s=30)
print(df_ttest_soanc_stanc[significant])
significant = (df_ttest_soanc_stanc['p_val'] < 0.05) & ((df_ttest_soanc_stanc['t_stat']) < -5)
plt.scatter(df_ttest_soanc_stanc[significant]['t_stat'], -np.log10(df_ttest_soanc_stanc[significant]['p_val']), color='blue', s=30)
print(df_ttest_soanc_stanc[significant])
# add horizontal line at pval = 0.05
plt.axhline(-np.log10(0.05), color='black', linestyle='--')
# add vertical lines at t_stat = -5 and 5
plt.axvline(-5, color='black', linestyle='--')
plt.axvline(5, color='black', linestyle='--')
# add labels and title
plt.xlabel('t-statistic (Stiff Selected vs Stiff Ancestral)')
plt.ylabel('-log10(p-value)')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(-11.5, 11.1)
#plt.ylim(0, 5)
#plt.savefig("volcano_plot_stsel_stanc_pid_pathways.png", dpi=500, bbox_inches='tight')
plt.show()

# %%
df_hallmark_t = df_hallmark.T
# make barplot with error bars for pathway "HALLMARK_E2F_TARGETS" and color according to condition
import seaborn as sns
plt.figure(figsize=(6,5))
# get groups from column names
groups = [col.split("_Re")[0] for col in df_hallmark.columns]
df_hallmark_t['Group'] = groups
color_map = {
    'soft_selected': '#8E44AD', 
    'stiff_selected': '#27AE60', 
    'soft_ancestral': '#3498DB', 
    'stiff_ancestral': '#FF851B'}
# Map conditions to colors
colors = [color_map[cond] for cond in condition_array]
# %%

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec

# Create a figure with a custom GridSpec layout
fig = plt.figure(figsize=(8, 6))  # Adjust width if needed
gs = gridspec.GridSpec(1, 2, width_ratios=[2, 3])  # 2:3 ratio

# Plot 1: Wider panel (left)
ax1 = plt.subplot(gs[1])
pathway = "PID_RXR_VDR_PATHWAY" ##################################################### add pathway name here
sns.boxplot(x=df_hallmark_t['Group'], y=df_hallmark_t[pathway],
            palette=['#FF851B', '#3498DB', '#8E44AD'],
            order=['stiff_ancestral', 'soft_ancestral', 'soft_selected'],
            ax=ax1)
sns.stripplot(x=df_hallmark_t['Group'], y=df_hallmark_t[pathway],
              color='black', size=6, jitter=True, dodge=True,
              order=['stiff_ancestral', 'soft_ancestral', 'soft_selected'],
              ax=ax1)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, fontsize=16)
ax1.set_ylim(0.20, 0.28)
ax1.set_ylabel('', fontsize=16)
ax1.set_title('Soft Adaptation', fontsize=16)

# Plot 2: Narrower panel (right)
ax2 = plt.subplot(gs[0])
sns.boxplot(x=df_hallmark_t['Group'], y=df_hallmark_t[pathway],
            palette=['#FF851B', '#27AE60'],
            order=['stiff_ancestral', 'stiff_selected'],
            ax=ax2)
sns.stripplot(x=df_hallmark_t['Group'], y=df_hallmark_t[pathway],
              color='black', size=6, jitter=True, dodge=True,
              order=['stiff_ancestral', 'stiff_selected'],
              ax=ax2)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, fontsize=16)
ax2.set_ylim(0.20, 0.28)
ax2.set_ylabel('Normalized Enrichment Score (NES)', fontsize=16)
ax2.set_title('Stiff Adaptation', fontsize=16)

# Tidy layout
plt.tight_layout()
plt.savefig("boxplot_PID_RXR_VDR_PATHWAY.png", dpi=500, bbox_inches='tight')


#%%

############## for genes

df = pd.read_csv("G:/projects_part_2/tanmay_lele/rpkm_rna.txt", sep="\t")
# remove gene_id column and make symbol as index
df = df.set_index("symbol").drop(columns=["gene_id"])
# perform log2 transformation on the dataframe
df_log2 = np.log2(df + 1)
df_t = df_log2.T

groups = [col.split("_Re")[0] for col in df_log2.columns]
df_t['Group'] = groups
color_map = {
    'Adp1kpa': '#8E44AD', 
    'Adp3kpa': '#27AE60', 
    'Ans1kpa': '#3498DB', 
    'Ans3kpa': '#FF851B'}
# Map conditions to colors
colors = [color_map[cond] for cond in df_t['Group']]

# Create a figure with a custom GridSpec layout
fig = plt.figure(figsize=(8, 6))  # Adjust width if needed
gs = gridspec.GridSpec(1, 2, width_ratios=[2, 3])  # 2:3 ratio
# remove gridlines


# Plot 1: Wider panel (left)
ax1 = plt.subplot(gs[1])
pathway = "RACGAP1" ##################################################### add gene name here
sns.boxplot(x=df_t['Group'], y=df_t[pathway],
            palette=['#FF851B', '#3498DB', '#8E44AD'],
            order=['Ans3kpa', 'Ans1kpa', 'Adp1kpa'],
            ax=ax1)
sns.stripplot(x=df_t['Group'], y=df_t[pathway],
              color='black', size=6, jitter=True, dodge=True,
              order=['Ans3kpa', 'Ans1kpa', 'Adp1kpa'],
              ax=ax1)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, fontsize=16)
ax1.set_ylabel('', fontsize=16)
ax1.set_title('Soft Adaptation', fontsize=16)

# Plot 2: Narrower panel (right)
ax2 = plt.subplot(gs[0])
sns.boxplot(x=df_t['Group'], y=df_t[pathway],
            palette=['#FF851B', '#27AE60'],
            order=['Ans3kpa', 'Adp3kpa'],
            ax=ax2)
sns.stripplot(x=df_t['Group'], y=df_t[pathway],
              color='black', size=6, jitter=True, dodge=True,
              order=['Ans3kpa', 'Adp3kpa'],
              ax=ax2)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, fontsize=16)
ax1.set_ylim(min(df_t[pathway])-0.1, max(df_t[pathway]) + 0.2)
ax2.set_ylim(min(df_t[pathway])-0.1, max(df_t[pathway]) + 0.2)
ax2.set_ylabel('log2(RPKM+1)', fontsize=16)
ax2.set_title('Stiff Adaptation', fontsize=16)

# Tidy layout
plt.tight_layout()

# compute ttest between Ans1kpa and Adp1kpa
from scipy import stats
soft_anc = list(df_t[df_t['Group'] == 'Ans1kpa'][pathway])
soft_sel = list(df_t[df_t['Group'] == 'Adp1kpa'][pathway])
stiff_anc = list(df_t[df_t['Group'] == 'Ans3kpa'][pathway])
stiff_sel = list(df_t[df_t['Group'] == 'Adp3kpa'][pathway])

print("soft_anc vs soft_sel", stats.ttest_ind(soft_anc, soft_sel))
print("stiff_anc vs stiff_sel", stats.ttest_ind(stiff_anc, stiff_sel))
print("stiff_anc vs soft_anc", stats.ttest_ind(stiff_anc, soft_anc))


plt.savefig("boxplot_gene"+pathway+".png", dpi=500, bbox_inches='tight')
# %%
import gseapy as gp
ssgsea_scores = gp.ssgsea(data=df_log2, gene_sets="G:/projects_part_2/tanmay_lele/ssGSEA/genesets/c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt", outdir=None, sample_norm_method='rank', no_plot=True, processes=4, min_size = 5, max_size=10000)
ssgsea_df = ssgsea_scores.res2d
#groups = [col.split("_Re")[0] for col in ssgsea_df.columns]

# pivot the dataframe to have groups as index and NES as values
ssgsea_df = ssgsea_df.pivot(index='Name', columns='Term', values='NES').reset_index()

ssgsea_df['Group'] = [value.split('_')[0] for value in ssgsea_df['Name']]

color_map = {
    'Adp1kpa': '#8E44AD', 
    'Adp3kpa': '#27AE60', 
    'Ans1kpa': '#3498DB', 
    'Ans3kpa': '#FF851B'}
# Map conditions to colors
colors = [color_map[cond] for cond in ssgsea_df['Group']]

# Create a figure with a custom GridSpec layout
fig = plt.figure(figsize=(8, 6))  # Adjust width if needed
gs = gridspec.GridSpec(1, 2, width_ratios=[2, 3])  # 2:3 ratio

# Plot 1: Wider panel (left)
ax1 = plt.subplot(gs[1])
pathway = "HALLMARK_HYPOXIA"
sns.boxplot(x=ssgsea_df['Group'], y=ssgsea_df[pathway],
            palette=['#FF851B', '#3498DB', '#8E44AD'],
            order=['Ans3kpa', 'Ans1kpa', 'Adp1kpa'],
            ax=ax1)
sns.stripplot(x=ssgsea_df['Group'], y=ssgsea_df[pathway],
              color='black', size=6, jitter=True, dodge=True,
              order=['Ans3kpa', 'Ans1kpa', 'Adp1kpa'],
              ax=ax1)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, fontsize=16)
ax1.set_ylabel('', fontsize=16)
ax1.set_title('Soft Adaptation', fontsize=16)

# Plot 2: Narrower panel (right)
ax2 = plt.subplot(gs[0])
sns.boxplot(x=ssgsea_df['Group'], y=ssgsea_df[pathway],
            palette=['#FF851B', '#27AE60'],
            order=['Ans3kpa', 'Adp3kpa'],
            ax=ax2)
sns.stripplot(x=ssgsea_df['Group'], y=ssgsea_df[pathway],
              color='black', size=6, jitter=True, dodge=True,
              order=['Ans3kpa', 'Adp3kpa'],
              ax=ax2)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, fontsize=16)
ax1.set_ylim(min(ssgsea_df[pathway])-0.05, max(ssgsea_df[pathway]) + 0.05)
ax2.set_ylim(min(ssgsea_df[pathway])-0.05, max(ssgsea_df[pathway]) + 0.05)
ax2.set_ylabel('Normalized Enrichment Score (NES)', fontsize=16)
ax2.set_title('Stiff Adaptation', fontsize=16)
plt.title(pathway, fontsize=16)

# Tidy layout
plt.tight_layout()
#plt.savefig("boxplot_isolation_stress.png", dpi=500, bbox_inches='tight')
#plt.savefig("boxplot_importin_signature.png", dpi=500, bbox_inches='tight')

# compute ttest between Ans1kpa and Adp1kpa
from scipy import stats
soft_anc = list(ssgsea_df[ssgsea_df['Group'] == 'Ans1kpa'][pathway])
soft_sel = list(ssgsea_df[ssgsea_df['Group'] == 'Adp1kpa'][pathway])
stiff_anc = list(ssgsea_df[ssgsea_df['Group'] == 'Ans3kpa'][pathway])
stiff_sel = list(ssgsea_df[ssgsea_df['Group'] == 'Adp3kpa'][pathway])

print("soft_anc vs soft_sel", stats.ttest_ind(soft_anc, soft_sel))
print("stiff_anc vs stiff_sel", stats.ttest_ind(stiff_anc, stiff_sel))
print("stiff_anc vs soft_anc", stats.ttest_ind(stiff_anc, soft_anc))

#%%

sns.scatterplot(x=ssgsea_df["HALLMARK_GLYCOLYSIS"], y=ssgsea_df["HALLMARK_OXIDATIVE_PHOSPHORYLATION"], hue=ssgsea_df['Group'], palette=color_map)

# %%





# --------------------------------- KEGG ssGSEA plotting ---------------------------------
import gseapy as gp
ssgsea_scores = gp.ssgsea(data=df_log2, gene_sets="G:/projects_part_2/tanmay_lele/ssGSEA/genesets/metabolic_signatures.gmt", outdir=None, sample_norm_method='rank', no_plot=True, processes=4, min_size = 5, max_size=10000)
ssgsea_df = ssgsea_scores.res2d
#groups = [col.split("_Re")[0] for col in ssgsea_df.columns]

# pivot the dataframe to have groups as index and NES as values
ssgsea_df = ssgsea_df.pivot(index='Name', columns='Term', values='NES').reset_index()
print(ssgsea_df.columns)
# make Name column as index
ssgsea_df.set_index('Name', inplace=True)
# z normalize the NES scores for each pathway
ssgsea_df = (ssgsea_df - ssgsea_df.mean()) / ssgsea_df.std()
#%%
# make a sns clustermap of the z normalized NES scores with colors for each group
import seaborn as sns
color_map = {
    'Adp1kpa': '#8E44AD', 
    'Adp3kpa': '#27AE60', 
    'Ans1kpa': '#3498DB', 
    'Ans3kpa': '#FF851B'}
groups = [value.split('_')[0] for value in ssgsea_df.index]

row_colors = [color_map[group] for group in groups]
# make sure datatype is correct
ssgsea_df = ssgsea_df.astype(float)
g = sns.clustermap(
        ssgsea_df,
        row_colors=row_colors,
        cmap="vlag",  # A diverging colormap is good for z-scores (blue-white-red)
        metric="euclidean",
        method="ward",
        figsize=(10, 16)
    )
# angle x labels 45 degrees
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=90, ha='right')
plt.show()
#%%



ssgsea_df['Group'] = [value.split('_')[0] for value in ssgsea_df['Name']]

color_map = {
    'Adp1kpa': '#8E44AD', 
    'Adp3kpa': '#27AE60', 
    'Ans1kpa': '#3498DB', 
    'Ans3kpa': '#FF851B'}
# Map conditions to colors
colors = [color_map[cond] for cond in ssgsea_df['Group']]

#%%

#%%


    # Create a figure with a custom GridSpec layout
fig = plt.figure(figsize=(8, 6))  # Adjust width if needed
gs = gridspec.GridSpec(1, 2, width_ratios=[2, 3])  # 2:3 ratio
    # Plot 1: Wider panel (left)
ax1 = plt.subplot(gs[1])
pathway = i
sns.boxplot(x=ssgsea_df['Group'], y=ssgsea_df[pathway],
            palette=['#FF851B', '#3498DB', '#8E44AD'],
            order=['Ans3kpa', 'Ans1kpa', 'Adp1kpa'],
            ax=ax1)
sns.stripplot(x=ssgsea_df['Group'], y=ssgsea_df[pathway],
                color='black', size=6, jitter=True, dodge=True,
                order=['Ans3kpa', 'Ans1kpa', 'Adp1kpa'],
                ax=ax1)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, fontsize=16)
ax1.set_ylabel('', fontsize=16)
ax1.set_title('Soft Adaptation', fontsize=16)

    # Plot 2: Narrower panel (right)
ax2 = plt.subplot(gs[0])
sns.boxplot(x=ssgsea_df['Group'], y=ssgsea_df[pathway],
                palette=['#FF851B', '#27AE60'],
                order=['Ans3kpa', 'Adp3kpa'],
                ax=ax2)
sns.stripplot(x=ssgsea_df['Group'], y=ssgsea_df[pathway],
                color='black', size=6, jitter=True, dodge=True,
                order=['Ans3kpa', 'Adp3kpa'],
                ax=ax2)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, fontsize=16)
ax1.set_ylim(min(ssgsea_df[pathway])-0.05, max(ssgsea_df[pathway]) + 0.05)
ax2.set_ylim(min(ssgsea_df[pathway])-0.05, max(ssgsea_df[pathway]) + 0.05)
ax2.set_ylabel('Normalized Enrichment Score (NES)', fontsize=16)
ax2.set_title('Stiff Adaptation', fontsize=16)
plt.title(pathway, fontsize=16)

    # Tidy layout
plt.tight_layout()
plt.show()
plt.close()
# %%
