library(WGCNA)
library(pheatmap)
library(ggplot2)
library(reshape2)

allowWGCNAThreads(nThreads = 10)
data <- read.delim("raw_and_filt_input_data/rpkm_rna_log_normalised.txt", row.names = 1)
gsg = goodGenes(data, verbose = 3)
View(data)
data_filt = data[, gsg]
View(data_filt)
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
powers = c(c(1:20))

#make a plot to decide the power values, recommended to choose a power with R2 value greater 0.8, if there are multiple
#values with R2 greater than 0.8, choose the one after which R2 value saturates (inflection point)
sft = pickSoftThreshold(data_filt, powerVector = powers, verbose = 3, networkType = "signed") #look again for the network type; default is unsigned, other options available are signed and hybrid signed
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"), family = "Arial")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red", family = "Arial")
abline(h = 0.90, col = "red")

png("scale_independence_plot.png", width = 3000, height = 2400, res = 300)  # 300 DPI
par(family = "Arial", font = 2)
point_colors <- ifelse(sft$fitIndices[,1] == 11, "red", "blue")
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",  # no default points yet
     main = paste("Scale independence"),
     family = "Arial",
     xaxt = "n", font = 2)  # we'll customize x-axis

# Add grid in the background
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

# Custom x-axis with all numbers visible
axis(1, 
     at = sft$fitIndices[,1], 
     labels = sft$fitIndices[,1], 
     family = "Arial", font = 2)
axis(2, font = 2, family = "Arial")

# Add scatter points as dots
points(sft$fitIndices[,1], 
       -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       pch = 19, col = point_colors)

# Add red text labels for powers near the dots
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers,
     cex = cex1,
     pos = 3,  # above points
     col = point_colors,
     family = "Arial")
y_val_power11 <- -sign(sft$fitIndices[,3])[sft$fitIndices[,1] == 11] *
  sft$fitIndices[,2][sft$fitIndices[,1] == 11]

# Draw horizontal line at that y-value
abline(h = y_val_power11, col = "red", lty = 2, lwd = 2)
dev.off()

'plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")'

png("mean_connectivity_plot.png",width = 3000, height = 2400, res = 300)
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"),
     family = "Arial",
     xaxt = "n", font = 2)

# Add grid
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

# Custom x-axis
axis(1, 
     at = sft$fitIndices[,1], 
     labels = sft$fitIndices[,1], 
     family = "Arial", font = 2)
axis(2, font = 2, family = "Arial")

# Add points
points(sft$fitIndices[,1], 
       sft$fitIndices[,5],
       pch = 19, col = point_colors)

# Add text labels
text(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     labels = powers,
     cex = cex1,
     pos = 3,  # above points
     col = point_colors,
     family = "Arial")

# Horizontal line for power 11
y_val_power11 <- sft$fitIndices[,5][sft$fitIndices[,1] == 11]
abline(h = y_val_power11, col = "red", lty = 2, lwd = 2)
dev.off()


fit = sft$fitIndices
power = sft$powerEstimate #used 11, but also check how results vary with power 15 
TOM = TOMsimilarityFromExpr(data_filt, power = 11, networkType = "signed",corType = "pearson", TOMType = "signed", nThreads = 10)
dissTOM = 1-TOM

#plot the gene tree
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

#identify modules using cutreedynamic
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 100) 
table(dynamicMods)
length(table(dynamicMods))

#convert the numeric label into colors
dynamicColors = labels2colors(dynamicMods) #instead of numbers, assign colors to modules
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(data_filt, colors = dynamicColors)
MEs = MEList$eigengenes #every module has an eigen gene profile across samples

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)

png("clustering_of_MEs.png",width = 3000, height = 2400, res = 300)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres=0.40 #check this threshold again if it works and try with other threshold values
abline(h=MEDissThres, col = "red")
dev.off()

merge = mergeCloseModules(data_filt, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors #module gene map  
mergedMEs = merge$newMEs #contains module eigengenes for the merged modules

#Plot merged module tree
png("cluster_dendrogram.png",width = 3000, height = 2400, res = 300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

group <- factor(sub("_Rep.*", "", rownames(mergedMEs)))
table(group)

#do one-way anova test on eigengenes of the merged modules (after merging the modules, we get 18 modules from 29)
anova_results <- apply(mergedMEs, 2, function(eig) {
  summary(aov(eig ~ group))
})
raw_pvals <- sapply(anova_results, function(x) x[[1]]["group", "Pr(>F)"]) #get the raw p values from the above 1 way anova test
adjusted_pvals <- p.adjust(raw_pvals, method = "fdr") #do fdr correction

#combine the p values and adjusted p values
pval_table <- data.frame(Module = names(raw_pvals),Raw_P = raw_pvals, Adjusted_P = adjusted_pvals)

#reorder the p adj values
pval_table_sorted <- pval_table[order(pval_table$Adjusted_P), ]
pval_table_sorted$Color <- gsub("^ME", "", pval_table_sorted$Module)

write.table(pval_table_sorted, file = "adjusted_pvalues_table.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

par(mar = c(10, 4, 4, 2))
png("adjusted_pvalues_plot.png", width = 1400, height = 800)

barplot(
  height = pval_table_sorted$Adjusted_P,
  names.arg = pval_table_sorted$Module,
  las = 2,                             # Rotate x-axis labels
  col = pval_table_sorted$Color,      # Use extracted colors
  main = "Adjusted P-values (FDR) for Module Eigengenes",
  ylab = "FDR (Adjusted P-value)",
  xlab = "Modules",
  cex.names = 0.8,
  ylim = c(0, max(pval_table_sorted$Adjusted_P) * 1.1)
)
dev.off()

#get significant module colors
signif_colors <- gsub("^ME", "", pval_table_sorted$Module[pval_table_sorted$Adjusted_P < 0.05])
names(mergedColors) <- colnames(data_filt)
sig_genes <- names(mergedColors)[mergedColors %in% signif_colors] #genes in the significant modules
data_sig = data_filt[,sig_genes] #filter the gene expression values of genes in significant modules
names(mergedColors)
annotation_row <- data.frame(Group = group)
rownames(annotation_row) <- rownames(data_sig)

module_colors <- mergedColors[colnames(data_sig)] #filter the gene-module object with the genes present in the significant module

annotation_col <- data.frame(Module = module_colors)
module_names <- unique(annotation_col$Module)

rownames(annotation_col) <- colnames(data_sig)
scaled_data_sig = scale(data_sig)

col_means <- colMeans(data_sig)
col_sds <- apply(data_sig, 2, sd)

annotation_colors <- list(Module = setNames(module_names, module_names))

graphics.off()
png("heatmap_of_genes_in_significant_modules.png", width = 801, height = 500, res = 300)

pheatmap(
    t(scaled_data_sig),
    annotation_col = annotation_row,
    annotation_row = annotation_col,
    annotation_colors = annotation_colors,
    show_colnames = TRUE,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    color = colorRampPalette(c("purple4", "white", "darkorange3"))(100),
    main = "Heatmap of Genes in Significant Modules",
    filename = "heatmap_of_genes_in_significant_modules.png", dpi = 300)

write.table(mergedMEs, file = "merged_module_eigengenes.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(data_filt, file = "rpkm_rna_log_normalised_good_genes.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(mergedColors, file = "merged_module_colors.txt", quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
write.table(data_sig, file = "rpkm_rna_log_normalised_sig_module_genes.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(scaled_data_sig, file = "rpkm_rna_log_&_z_normalised_sig_module_genes.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(module_colors,file = "module_colors_for_data_sig.txt",quote = FALSE,row.names = TRUE,col.names = FALSE,sep = "\t")

#check the correlation between the gene expression and the module eigengene
print(mergedMEs)
gene_module_cors <- cor(data_sig, mergedMEs, use = "p")

positive_regulators <- list()
negative_regulators <- list()

for (mod in signif_colors) {
  eigengene_name <- paste0("ME", mod)
  
  # Get genes in this module
  mod_genes <- names(mergedColors)[mergedColors == mod]
  
  # Filter gene correlations with the current module eigengene
  cor_vals <- gene_module_cors[mod_genes, eigengene_name]
  
  # Split into positive and negative regulators
  positive_regulators[[mod]] <- names(cor_vals[cor_vals > 0.5])
  negative_regulators[[mod]] <- names(cor_vals[cor_vals < -0.5])
}

for (mod in signif_colors) {
  pos_file <- paste0("positive_regulators_0.5_", mod, ".txt")
  neg_file <- paste0("negative_regulators_0.5", mod, ".txt")
  
  writeLines(positive_regulators[[mod]], con = pos_file)
  writeLines(negative_regulators[[mod]], con = neg_file)
}

signif_modules <- paste0("ME", signif_colors)
signif_MEs <- mergedMEs[, signif_modules]
signif_MEs$Sample <- rownames(mergedMEs)

#plot for num of genes in the signififcant modules
sig_gene_modules <- mergedColors[names(mergedColors) %in% sig_genes]

# Count genes per module
module_counts <- as.data.frame(table(sig_gene_modules))
colnames(module_counts) <- c("Module", "Count")

# Plot bar plot
ggplot(module_counts, aes(x = Module, y = Count, fill = Module)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Number of Genes in Significant Modules",
    x = "Module",
    y = "Gene Count"
  )

#melt the data for ggplot
melted_MEs <- melt(signif_MEs, id.vars = "Sample", variable.name = "Module", value.name = "Eigengene")

#add group information
melted_MEs$Group <- factor(sub("_Rep.*", "", melted_MEs$Sample))
melted_MEs$Module <- gsub("^ME", "", melted_MEs$Module)

write.table(
  melted_MEs,
  file = "melted_MEs_for_sig_modules.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggplot(melted_MEs, aes(x = Sample, y = Eigengene, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Module, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(title = "Module Eigengene Expression (Significant Modules)",
       y = "Module Eigengene Value",
       x = "Sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        strip.text = element_text(size = 10))