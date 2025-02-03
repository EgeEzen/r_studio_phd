#Installing the DESeq2 and other necessary libraries from Bioconductor project & others

version
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("DESeq2","RColorBrewer", "pheatmap","tidyverse","apeglm","EnhancedVolcano"))
install.packages("ggrepel")
install.packages('devtools')
install.packages('ggthemr')
devtools::install_github('Mikata-Project/ggthemr',force=TRUE)

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(apeglm)
library(EnhancedVolcano)
library(writexl)
library(ggrepel)
library(tidyr)
library(limma)
library(dplyr)
library(biomaRt)
library(ggthemr)
library(ggsignif)
library(envalysis)

#Setting the working directory

setwd("~/Desktop/PhD/2024 Second Term/RNA-seq all/")
directory <- "~/Desktop/PhD/2024 Second Term/RNA-seq all/"

# read the raw data
count_data <- read.csv("count_data_all.csv",row.names =1)
# read the metadata
metadata <- read.csv("metadata.csv",row.names=1)
metadata$Nr. <- NULL

# --------- # 
#to do healthy vs ra (not tnf stimulated samples we have to get rid of TNF stimulated)
metadata <- metadata[metadata$condition != "Stimulated", ]
count_data <- count_data[,colnames(count_data) %in% rownames(metadata)]
# -------#
#to do tnf vs non-tnf
metadata <- metadata[metadata$condition != "not_applicable",]
count_data <- count_data[,colnames(count_data) %in% rownames(metadata)]
# -------#

  
# creating the DESeq2 object
# design is the group (the base will be controlled in later stage)
# technology as a batch effect -> this only effects the differential gene expression (logfold changes, padj etc)
# for visualization batch effect remover will be used (only for visualization, not for d.g.e. analysis)
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData =  metadata,
                              design = ~ technology + group) #group or condition!
#it will give warning message, ignore it (about the base, it will automatically convert to factors)

#factor level is set here, this means that the base is now untreated (NS) condition
dds$condition <- relevel(dds$group, ref = "Healthy") #Healthy or Non_stimulated

#seeing the object
dds

#calculating the median of ratios, and adding to our object 
#(this is already done in while the deseq2 object, but you can use these values separately)
dds <- estimateSizeFactors(dds)

#extracting the normalized counts from object
nrm_counts <- counts(dds, normalized = TRUE)

# Extracting the transformed values (DESEQ2 says blind should be false to reduce noise for downstream analysis)
# DESEQ2 Manual says that bias introduced in blind=FALSE argument is not cruical
vsd <- vst(dds, blind=FALSE) 

#computing pairwise correlation values, first create a matrix from vsd object, then compute correlation
vsd_mat <- assay(vsd)
vsd_cor<- cor(vsd_mat)

#then do the correlation heatmap plot
#pheatmap data gets the correlation matrix and in the annotation part it uses the metadata 
#it matches the rownames of the metadataof to column names of the correlation matrix

#select function gets the specific column from a dataframe
pheatmap(vsd_cor, annotation_col = metadata,
         main = "Hierarchical Heatmap", filename ="clustermap.png",width= 10, height = 8, show_rownames=FALSE, show_colnames = TRUE)
#don't forget the change the name of the map (main argument)

#you can also remove the batch effects here for visualization
vsd_copy <- vsd
mat_copy <- assay(vsd_copy)
mm <- model.matrix(~group, colData(vsd_copy)) #group or condition
mat_copy <- removeBatchEffect(mat_copy, batch=vsd_copy$technology, design=mm)
assay(vsd_copy) <- mat_copy
#also for normalized counts
nrm_counts_copy <- nrm_counts
nrm_counts_copy <- removeBatchEffect(nrm_counts_copy,batch=vsd_copy$technology, design=mm)

#Doing PCA Analysis
pcaData<-plotPCA(vsd, intgroup="group",returnData=TRUE)
#or with the batch effect removed 
pcaData<-plotPCA(vsd_copy, intgroup="group",returnData=TRUE)
#when returnData is true, the graph is assigned to a variable, then you can use ggplot to modify it

percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$group <- metadata$group
pcaData$condition <- metadata$condition
pcaData$joint <- metadata$joint
pcaData$gender <- metadata$gender
pcaData$condition <- ifelse(pcaData$condition == "not_applicable", "Not Applicable",ifelse(pcaData$condition == "Non_stimulated", "Non-stimulated","Stimulated"))

#now to use ggplot

#ggthemr('pale', layout= 'clear')

p <- ggplot(pcaData, aes(PC1, PC2, color=condition, label=group )) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(title= "PCA Plot") +
  geom_label_repel(fontface = "bold", nudge_x = 1, show.legend = FALSE) +
  theme_publish()

ggsave(
  plot = p,
  filename = "pca_after_be_removed.png",
  bg = "transparent",
  width = 6, height = 6
)

#Don't use "color" in the legend, remember to change title of the plot
#ggrepel is very useful if the annotations overlap

#finally doing the DESeq analysis, dds is now the deseq2 object
dds <- DESeq(dds)

#plotting dispersion estimate graph (look further for manual)
plotDispEsts(dds)

#extracting the results
#If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
result <- results(dds, alpha = 0.05)

#viewing the results
summary(result)
result

#you can see the up and down regulation here (LFC > 0 is given, we will further limit that)
#Row names are extracted to gene_id column
result$gene_id <- rownames(result)
#extracting this results to an excel file

#write_xlsx(as.data.frame(result),"deseq2_genes.xlsx")

# getting the coefficients for the log shrinkage 
resultsNames(dds)

# log shrinkage, apeglm method will be used (2018 paper)
results_vs <- lfcShrink(dds,
                         coef="group_RA_vs_Healthy",type="apeglm")
#results to save 
results_to_save <- as.data.frame(results_vs)
results_to_save$gene <- rownames(results_to_save)
# then you can save this file
write_xlsx(results_to_save,"results_vs_healthy_vs_ra.xlsx" )

#save RDS for shiny app
results_vs_df <- as.data.frame(results_vs)
healthy_vs_ra_shiny <- list(nrm_counts_copy,mat_copy,results_vs_df, metadata)
saveRDS(healthy_vs_ra_shiny, file = "healthy_vs_ra_shiny.rds")
#------------#

# from now on continue the analysis with the shrinkage version of the results

#MA plots 
plotMA(results_vs, ylim=c(-4,4))

#plotting counts, if you want to see the change in a specific gene use this function below (change the gene argument)
#plotCounts(dds, gene=which.min(result_shr$svalue), intgroup="condition")

#only getting significants and ordering (FDR 0.05, fc 1.5 is picked here)
resSig <- results_vs[order(results_vs$padj),]
resSig <- subset(resSig, padj <= 0.05)
resSig_with_fc1 <- subset(resSig, ((log2FoldChange >= 0.585) | (log2FoldChange <= -0.585) ))

#getting the normalized counts for ONLY significant genes (then we will create a heatmap)
resSig_with_normalized_counts <- merge(as.data.frame(resSig), as.data.frame(nrm_counts), by="row.names")
resSig_with_be_removed_normalized_counts <- merge(as.data.frame(resSig), as.data.frame(nrm_counts_copy), by="row.names")

#exporting as xlsx
write_xlsx(resSig_with_normalized_counts,"results_vs_sig_with_normalized_counts.xlsx")
write_xlsx(resSig_with_be_removed_normalized_counts,"results_vs_sig_with_be_removed_normalized_counts.xlsx")

#resSig with the batch effect removed vsd counts -> this is useful for visualization
resSig_with_be_removed_vsd_counts <- merge(as.data.frame(resSig), as.data.frame(mat_copy), by="row.names")
resSig_with_fc1_with_be_removed_vsd_counts <- merge(as.data.frame(resSig_with_fc1), as.data.frame(mat_copy), by="row.names")

# draw heatmap
custom_colors <- colorRampPalette(c("blue", "white", "red"))(20)

p1 <- pheatmap(resSig_with_be_removed_vsd_counts[,7:ncol(resSig_with_be_removed_vsd_counts)], cluster_rows=TRUE, show_rownames=FALSE,scale = "row",color = custom_colors,
         cluster_cols=TRUE, annotation_col=metadata[,"condition",drop=FALSE],width=10,height=10, #group or condition
         main = "TNF vs Non-TNF\nin RA Heatmap",treeheight_row = 0) #healthy vs ra or tnf vs non tnf in ra

ggsave(
  plot = p1,
  filename = "heatmap_resSig.png",
  bg = "transparent",
  width = 10, height = 6
)

#with lfc1
p1 <- pheatmap(resSig_with_fc1_with_be_removed_vsd_counts[,7:ncol(resSig_with_fc1_with_be_removed_vsd_counts)], cluster_rows=TRUE, show_rownames=FALSE,scale = "row",color = custom_colors,
         cluster_cols=TRUE, annotation_col=metadata[,"condition",drop=FALSE],width=10,height=10,
         main = "TNF vs Non-TNF\nin RA Heatmap",treeheight_row = 0)

ggsave(
  plot = p1,
  filename = "heatmap_resSig_fc1.png",
  bg = "transparent",
  width = 10, height = 6
)


# --------- draw the volcano plot------ #

results_vs_df <- as.data.frame(results_vs)

#volcano plot drawing starts here:
#add a column to classify points as significant or not
results_vs_df$Classification <- with(results_vs_df, ifelse(padj <= 0.05 & log2FoldChange >= 0.585, "Upregulated", 
                                                                     ifelse(padj <= 0.05 & log2FoldChange <= -0.585, "Downregulated", 
                                                                            "Not Significant")))
#calculating the -log10(padj)
results_vs_df$log_adjP <- -log10(results_vs_df$padj)

#selecting points for geom text label
label_data <- subset(results_vs_df, (log2FoldChange >= 0.585 | log2FoldChange <= -0.585) & Classification != "Not Significant" )

top_25_up <- label_data %>% filter(Classification == "Upregulated") %>% arrange(desc(log2FoldChange)) %>% head(25)
top_25_down <- label_data %>% filter(Classification == "Downregulated") %>% arrange(log2FoldChange) %>% head(25)
label_data <- bind_rows(top_25_up, top_25_down)

#rownames(label_data) <- NULL
#label_data$rownames <- rownames(label_data)
#label_data <- label_data[,c(5,6,19,21,22)]


p <- ggplot(results_vs_df, aes(x = log2FoldChange, y = log_adjP)) +
  geom_point(aes(color = Classification), size = 1,alpha=0.8) + #alpha 0.5-0.8 looks "cool"
  scale_color_manual(values = c("Upregulated" = "#EE4E4E", "Downregulated" = "#2A629A", "Not Significant" = "lightgray")) +
  geom_vline(xintercept = c(-0.585,0.588), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Differentially Expressed Genes TNF vs non-TNF in RA",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted P-Value)",
    fill = "Classification"
  ) +   # Increase plot margins if necessary
  ylim(0,NA) +
  geom_label_repel(
    data = label_data,
    aes(label = rownames(label_data)),
    box.padding = unit(0.1, "lines"),
    segment.color = "#2B2A4C",
    max.overlaps= 50
  ) +
  guides(color = guide_legend(title = "Differential Expression", override.aes = list(size = 3)),alpha= FALSE) +
  theme_publish()

p

ggsave(
  plot = p,
  filename = "volcanoPlot_lf1.png",
  bg = "transparent",
  width = 10, height = 6
)

# ------ draw dotplot for a specific gene ------ #
  
gene_of_interest <- "ZFP36L1"
nrm_counts_filtered <- nrm_counts_copy[rownames(nrm_counts_copy) == gene_of_interest,,drop=FALSE]
nrm_counts_filtered <- as.data.frame(t(nrm_counts_filtered))
nrm_counts_filtered <- rownames_to_column(nrm_counts_filtered, var = "Sample") 
nrm_counts_filtered$Condition <- ifelse(grepl("Healthy",nrm_counts_filtered$Sample), "Healthy", "RA") 
#or stimulated and unstimulated
nrm_counts_filtered$Condition <- ifelse(grepl("TNF",nrm_counts_filtered$Sample), "TNF", "non-TNF") 
colnames(nrm_counts_filtered) <- c("Sample","Counts","Condition") 


p <- ggplot(nrm_counts_filtered, aes(x = Condition, y = Counts, fill = Condition)) +
  geom_violin(trim = TRUE, alpha = 0.8) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +
  labs(title = paste0("Gene Expression of\n", gene_of_interest," in Healthy vs RA"), 
       x = "Condition", y = "Normalized Counts") +
  theme_publish() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        legend.position = "none")  


checking_condition <- (results_vs_df[rownames(results_vs_df) == gene_of_interest,,drop=FALSE ])$padj
if (checking_condition <= 0.05){
  p <- p + geom_signif(comparisons = list(c("non-TNF", "TNF")), # TNF non TNF or Healthy Ra
                       annotations = "*",  # Star for significance
                       y_position = max(nrm_counts_filtered$Counts) * 1.05, 
                       tip_length = 0.02, textsize = 7)
}
  
ggsave(
  plot = p,
  filename = paste0("dotplot_of_",gene_of_interest,".png"),
  bg = "transparent",
  width = 5, height = 7
)

# --------- heatmap from selected genes ------- #

selected_genes <- c("RBPJ","AGAP2-AS1","ZFP36L1")
mat_copy_filtered <- mat_copy[rownames(mat_copy) %in% selected_genes, ,drop=FALSE]

#for the TNF vs nonTNF you have re-order the columns
metadata$condition <- factor(metadata$condition, levels = c("Non_stimulated", "Stimulated"))
ordered_samples <- rownames(metadata)[order(metadata$condition)]  # Get ordered sample names
mat_copy_filtered <- mat_copy_filtered[, ordered_samples]  # Reorder matrix columns

# draw heatmap
custom_colors <- colorRampPalette(c("blue", "white", "red"))(10)

p1 <- pheatmap(mat_copy_filtered, cluster_rows=TRUE, show_rownames=TRUE,show_colnames = FALSE,scale = "row",color = custom_colors,
               cluster_cols=FALSE, annotation_col=metadata[,"condition",drop=FALSE],width=10,height=10, #group or condition
               main = "Heatmap of Selected Genes\nin RA and Healthy",treeheight_row = 0)

ggsave(
  plot = p1,
  filename = "heatmap_vsd.png",
  bg = "transparent",
  width = 7, height = 2
)

