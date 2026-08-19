import pandas as pd
import numpy as np
import gseapy as gp

# ---------------------------------------------------------
# 1. Define File Paths
# ---------------------------------------------------------
hallmark_geneset_dir = "D:\\sejal\\RNA_seq_dec24th\\genesets.gmt"
rna_data_path = "D:\\sejal\\tanmay_lele\\rpkm_rna.txt"
output_file = "MSigDb_hallmark_pathway_nes_scores.txt"

# ---------------------------------------------------------
# 2. Load and Preprocess Data
# ---------------------------------------------------------
# Read the gene expression RPKM data
gene_exp_df = pd.read_csv(rna_data_path, sep='\t')

# Group column names by condition
soft_selected = [col for col in gene_exp_df.columns if col.startswith('soft_selected')]
stiff_selected = [col for col in gene_exp_df.columns if col.startswith('stiff_selected')]
soft_ancestral = [col for col in gene_exp_df.columns if col.startswith('soft_ancestral')]
stiff_ancestral = [col for col in gene_exp_df.columns if col.startswith('stiff_ancestral')]

# Combine all relevant columns
cols = soft_selected + stiff_selected + soft_ancestral + stiff_ancestral

# Apply log2(x + 1) transformation to normalize the expression values
log2_gene_exp_df = gene_exp_df.copy()
log2_gene_exp_df[cols] = log2_gene_exp_df[cols].apply(lambda x: np.log2(x + 1))

# Set the gene symbol as the index, filtering down to just the target columns
log2_gene_exp_filt = log2_gene_exp_df.set_index('symbol')[cols]

# ---------------------------------------------------------
# 3. Run ssGSEA (Single-Sample Gene Set Enrichment Analysis)
# ---------------------------------------------------------
# Calculate enrichment scores for each sample using the hallmark gene sets
ss = gp.ssgsea(
    data=log2_gene_exp_filt,
    gene_sets=hallmark_geneset_dir,
    outdir=None,                  # Set to None since we extract the dataframe directly
    sample_norm_method='rank',    # Normalizes expression profiles by ranking
    no_plot=True,                 # Skip generating plots
    min_size=10                   # Exclude gene sets with fewer than 10 genes
)

# Pivot the results so pathways are rows, samples are columns, and values are NES
nes = ss.res2d.pivot(index='Term', columns='Name', values='NES')

# ---------------------------------------------------------
# 4. Save Results
# ---------------------------------------------------------
# Export the NES dataframe to a tab-separated text file
nes.to_csv(output_file, sep="\t")
print(f"ssGSEA complete. Results saved to: {output_file}")