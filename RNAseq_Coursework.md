Gene_Expression_Coursework
================
2024-11-13

``` r
rm(list=ls())
```

``` r
data  <- load("~/NGS_Coursework/CancerAndParacancer.RData")

#set seed for reproducible results
set.seed(1234)
```

``` r
class(ExprMat)
```

    ## [1] "matrix" "array"

``` r
class(pDataMat)
```

    ## [1] "data.frame"

``` r
library(limma)
library(edgeR)
library(dplyr)

#Pseudocount of 1 to all FPKM count data so when log transformed there would be no errors
#Log2 transformation due to data being FPKM, this stabilises the variance of the data
dge_ExprMat <- log2((ExprMat + 1))

#Design matrix using paracancer as 1 and cancer as 0
#When doing DGE it will compare paracancer to cancer due to the order the the design matrix
#dge_design <- model.matrix(pDataMat)
pDataMat$Group <- factor(pDataMat$Group) #change the group to a factor as it's a categorical variable
dge_design <- model.matrix(~ Group, data = pDataMat) #create design matrix

#Makes a linear model based off of genes from series of repeats 
fit <- lmFit(dge_ExprMat, dge_design)

#Computes statistical analysis from the general linear model
#t, F, B, log of the differential expression between grouped data
fit <- eBayes(fit, trend = TRUE)
```

    ## Warning: Zero sample variances detected, have been offset away from zero

``` r
#Used this table to look for possible negative logFC values (looking for downregulated genes)
fit_test <- topTable(fit, coef=ncol(dge_design), number = Inf, sort.by = "logFC", p.value = 0.05)
#fit_test_2 <- topTable(fit, coef=ncol(dge_design), number = Inf, sort.by = "logFC", p.value = 0.1)

######UPREGULATED DEGS########

#Used this table to highlight the top 15 highest upregulated genes from our DGE
fit_sig_up <- topTable(fit, coef=ncol(dge_design), number = 15, sort.by = "logFC", p.value = 0.1)

#filter the gene names of the top 15 DEGs
top_fifteen_up_names <- rownames(fit_sig_up)

#filter the original FPKM data to contain only the data from the top 15 upregulated DEGs
top_fifteen_up <- as.data.frame(dge_ExprMat[top_fifteen_up_names,])

#filter feature data to contain only top 15 upregulated DEGs
top_fifteen_up_features <- as.data.frame(FeatureMat[top_fifteen_up_names,])
top_fifteen_up_genetype <- top_fifteen_up_features[5]
rownames(top_fifteen_up_genetype) <- top_fifteen_up_features$Symbol

#set row names of the filtered FPKM data to the GeneID
rownames(top_fifteen_up) <- top_fifteen_up_features$Symbol

#######DOWNREGULATED DEGS########

#filtered table to get top 15 downregulated genes from DGE
#fit_sig_down <- fit_test %>% filter (logFC< 0) %>% slice(1:15)
fit_sig_down_all <- fit_test[fit_test$logFC < 0, ]
fit_sig_down <- head(fit_sig_down_all, 15)


#filter the gene names of the top 15 DEGs
top_fifteen_down_names <- rownames(fit_sig_down)

#filter the original FPKM data to contain only the data from the top 15 downregulated DEGs
top_fifteen_down <- as.data.frame(dge_ExprMat[top_fifteen_down_names,])

#filter feature data to contain only top 15 upregulated DEGs
top_fifteen_down_features <- as.data.frame(FeatureMat[top_fifteen_down_names,])
top_fifteen_down_genetype <- top_fifteen_down_features[5]
rownames(top_fifteen_down_genetype) <- top_fifteen_down_features$Symbol

#set row names of the filtered FPKM data to the GeneID
rownames(top_fifteen_down) <- top_fifteen_down_features$Symbol
```

``` r
#Takes the top 15 up and down regulated genes that are differentiall expressed between cancer and paracancer
pca_genes <- bind_rows(top_fifteen_up, top_fifteen_down)

#Pseudocount of 1 and then log2 transformed to stabalise variance of FPKM data and stores it in a data frame so it can be used in ggplot
pca_genes <- log2((pca_genes+1))
T_pca_genes <- as.data.frame(t(pca_genes))

#PCA analysis using prcomp function
pca_deg <- prcomp(T_pca_genes, center = TRUE, scale. = TRUE) 
pca_deg_res <- as.data.frame(pca_deg$x)
```

``` r
#Creates the png file to store the plot in
png(filename="PCA.png")
PCA_plot <- ggplot(data = pca_deg_res, aes(x = PC1, y = PC2))+ # takes data from first 2 prinicpal components
  geom_point(aes(color = pDataMat$Group))+ #Uses cancer and paracancer groups for the data
  ggtitle("Principal component analysis of the sample data \nrelative to their gene expression")+
  scale_color_discrete(name = "Group") #Changes the legened title

#Plots the graph and then saves into the png file
plot(PCA_plot)
dev.off()                    
```

    ## png 
    ##   2

``` r
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

top_fifteen_up_zscore <- t(scale(t(top_fifteen_up))) #scale data to normalise by the z-score to allow for better comparison between samples

png("figures/heatmap_plot_up.png", width = 1500, height = 900, res = 150)
heat_plot_up <- pheatmap(top_fifteen_up_zscore, #scale = "row",
         cluster_rows = TRUE, #clustering based on rows (top 15 DEGs)
         cluster_cols = TRUE, #clustering based on columns (sample ID)
         clustering_method = "complete", #clustering method
         show_colnames = FALSE, #don't show sample ID (too many, not relevant)
         show_rownames = TRUE, #show Gene ID
         main = "Heatmap of Top 15 Upregulated Differentially Expressed Genes", #title
         annotation_row = top_fifteen_up_genetype, #group rows by genetype
         annotation_col = pDataMat, #group columns by cancerous or paracancerous
         col = brewer.pal(11, "RdYlBu"), #change colour scheme
         
         )
```

![](RNAseq_Coursework_files/figure-gfm/-%20Amy,%20Heatmap%20Plot%20for%20Upregulated%20DEGs-1.png)<!-- -->

``` r
dev.off()  
```

    ## png 
    ##   3

``` r
heat_plot_up
```

``` r
top_fifteen_down_zscore <- t(scale(t(top_fifteen_down))) #scale data to normalise by the z-score to allow for better comparison between samples

png("figures/heatmap_plot_down.png", width = 1500, height = 900, res = 150)
heat_plot_down <- pheatmap(top_fifteen_down_zscore, #scale = "row",
         cluster_rows = TRUE, #clustering based on rows (top 15 DEGs)
         cluster_cols = TRUE, #clustering based on columns (sample ID)
         clustering_method = "complete", #clustering method
         show_colnames = FALSE, #don't show sample ID (too many, not relevant)
         show_rownames = TRUE, #show Gene ID
         main = "Heatmap of Top 15 Downregulated Differentially Expressed Genes", #title
         annotation_row = top_fifteen_down_genetype, #group rows by genetype
         annotation_col = pDataMat, #group columns by cancerous or paracancerous
         col = brewer.pal(11, "RdYlBu"), #change colour scheme
         )
```

![](RNAseq_Coursework_files/figure-gfm/-%20Amy,%20Heatmap%20Plot%20for%20Downregulated%20DEGs-1.png)<!-- -->

``` r
dev.off()
```

    ## png 
    ##   3

``` r
heat_plot_down
```

``` r
library(ggplot2)
library(ggrepel)

#Add gene symbol to DEG results dataframe
#orderFeaturMat in the same way as the DEG output
FeaturMat_ordered <- FeatureMat[rownames(fit_test),]
#add a new symbol column to fit_test with the gene symbol
fit_test$Symbol <- FeaturMat_ordered[,2]

#create column for if gene is differentially expressed
fit_test$DEG <- "NO"

#Add row for upregulated genes
fit_test$DEG[fit_test$logFC > 5.749045 & fit_test$adj.P.Val<0.05] <- "UP"

#add row for downregulated genes
fit_test$DEG[fit_test$logFC < -3.236370 & fit_test$adj.P.Val<0.05] <- "DOWN"

fit_sig_up$Symbol <- rownames(top_fifteen_up)
fit_sig_down$Symbol <- rownames(top_fifteen_down)

#Create Volcano Plot
volcano_deg <- ggplot(data = fit_test, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_vline(xintercept = c(-3.236370,5.749045), col = "gray", linetype = "dashed")+
  geom_hline(yintercept = c(-log10(0.05)), col = "gray", linetype = "dashed")+
  scale_color_manual(values = c("blue", "grey", "red"), labels = c("Downregulated", "Not Significant", "Upregulated"))+
  labs ( title = "Volcano Plot Highlighting Differentially Expressed Genes",
         subtitle = "in Cancer and Paracancerous samples",
         x = "Log Fold Change",
         y = "-log10(adjusted p value)")+
  geom_point(aes(colour=DEG), alpha = 0.6, size=1)+
  ggrepel::geom_text_repel(data=fit_sig_up, aes(label=Symbol), max.overlaps =20, size =2)+
  ggrepel::geom_text_repel(data=fit_sig_down, aes(label=Symbol), max.overlaps =20, size =2)+
  theme_bw(base_size = 14)

print(volcano_deg)
```

![](RNAseq_Coursework_files/figure-gfm/-%20Amy%20-%20Volcano%20Plot%20of%20Top%2015%20Upregulated%20Genes-1.png)<!-- -->

``` r
ggsave("figures/volcano_plot.png", plot = volcano_deg, width = 8, height = 6, dpi = 300)
```

ORA-based pathway enrichment analysis

``` r
#Reactome Pathway
#install Reactome PA in Console

if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("ReactomePA")
}

#load ReactomePA
library(ReactomePA)

gene_list_up <- c(rownames(fit_sig_up)) 
gene_list_down <- c(rownames(fit_sig_down))
total_gene_list <- c(gene_list_up,gene_list_down)

#shows 
reactome_results_up <- enrichPathway(
  gene = gene_list_up, #list of top 15 upregulated DEGs
  organism = "human", #human genome
  pvalueCutoff = 0.1, #p value cutoff for results
  readable = TRUE) #include gene ID

reactome_results_down <- enrichPathway(
  gene = gene_list_down,
  organism = "human",
  pvalueCutoff = 0.1,
  readable = TRUE)
```

Barplots and map of reactome pathway for top 15 upregulated genes

``` r
#this was not used in the presentation - we compared with the ranked gene list analysis 

library(enrichplot)
up_reactome_barplot <- barplot(reactome_results_up)+ #showCategory = 10
  ggtitle("Barplot of Significant Biological Pathways Within\nTop 15 Upregulated DEGs from Reactome Database")+
  theme (
    plot.title = element_text (size = 12, face= "bold")
  )

up_reactome_emap <- emapplot(pairwise_termsim(reactome_results_up))+
  ggtitle("Map of Significant Biological Pathways Within\nTop 15 Upregulated DEGs from Reactome Database")+
  theme(
  geom_text(aes(label = description), size = 0.1),# Adjust node label size
  plot.title = element_text (face = "bold")
)

up_reactome_barplot
```

![](RNAseq_Coursework_files/figure-gfm/-%20Amy%20-%20Reactome%20Plots%20for%20Upregulated%20genes-1.png)<!-- -->

``` r
up_reactome_emap
```

![](RNAseq_Coursework_files/figure-gfm/-%20Amy%20-%20Reactome%20Plots%20for%20Upregulated%20genes-2.png)<!-- -->

``` r
ggsave("figures/up_reactome_barplot.png", plot = up_reactome_barplot, width = 8, height = 6, dpi = 300)
ggsave("figures/up_reactome_emap.png", plot = up_reactome_emap, width = 8, height = 6, dpi = 300)
```

``` r
#used to compare with the ranked gene list analysis

#download library to make plots for enrichment analysis
library(enrichplot)
down_reactome_barplot <- barplot(reactome_results_down)+
  ggtitle("Barplot of Significant Biological Pathways Within\nTop 15 Downregulated DEGs from Reactome Database")+
  theme (
    plot.title = element_text (size = 12, face= "bold")
  )

#pairwise termism is important as it is a calculation of how similar two pathways are based on the genes they share
down_reactome_emap <- emapplot(pairwise_termsim(reactome_results_down))+
  ggtitle("Map of Significant Biological Pathways Within\nTop 15 Downregulated DEGs from Reactome Database")+
  theme(
  geom_text(aes(label = description), size = 0.1),# Adjust node label size
  plot.title = element_text (face = "bold")
)

down_reactome_barplot
```

![](RNAseq_Coursework_files/figure-gfm/-%20Amy%20-%20Reactome%20plots%20for%20top%2015%20downregulated%20genes-1.png)<!-- -->

``` r
down_reactome_emap
```

![](RNAseq_Coursework_files/figure-gfm/-%20Amy%20-%20Reactome%20plots%20for%20top%2015%20downregulated%20genes-2.png)<!-- -->

``` r
ggsave("figures/down_reactome_barplot.png", plot = down_reactome_barplot, width = 8, height = 6, dpi = 300)
ggsave("figures/down_reactome_emap.png", plot = down_reactome_emap, width = 8, height = 6, dpi = 300)
```

``` r
reactome_cnet_up <- cnetplot(reactome_results_up, showCategory = 20)+
  ggtitle("Map of Upregulated Genes Related to Significant Pathways")+
  theme(plot.title = element_text (face = "bold", hjust =0.5))

reactome_cnet_down <- cnetplot(reactome_results_down, showCategory = 20)+
  ggtitle("Map of Downregulated Genes Related to Significant Pathways")+
  theme(plot.title = element_text (face = "bold", hjust =0.5))

ggsave("figures/reactome_cnet_up.png", plot = reactome_cnet_up, width = 16, height = 6, dpi = 300)
ggsave("figures/reactome_cnet_down.png", plot = reactome_cnet_down, width = 16, height = 6, dpi = 300)
```

``` r
#create ranked list of genes ranked by logFC
#create a vector of logFC values
gene_list <- c(fit_test[,1])
fit_test_names <- rownames(fit_test) #create a vector for entrez IDs
names(gene_list) <- fit_test_names #add Entrez IDs to names in list of logFC

gene_list <- sort(gene_list, decreasing = TRUE) #sort by decreasing logFC
  
head(gene_list)  #check this is highest positive     
```

    ##   221476      125     1674     6422     6423   126393 
    ## 7.437844 6.877418 6.790997 6.671820 6.283476 6.116278

``` r
tail(gene_list)  #check this is most negative     
```

    ## 101929397       810      7138    204962     11082     23581 
    ## -3.502136 -3.502591 -3.561037 -3.735724 -4.135538 -4.144859

``` r
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")
```

``` r
library(clusterProfiler)
library(enrichplot)

reactome_results <-
  gsePathway(gene_list, 
                pvalueCutoff = 0.05, seed =TRUE) #seed to ensure reproducible results
```

    ## using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (0.18% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize = minSize, : There were 19
    ## pathways for which P-values were not calculated properly due to unbalanced (positive and negative)
    ## gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA. You can try to
    ## increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)

    ## Warning in serialize(data, node$con): 'package:stats' may not be available when loading
    ## Warning in serialize(data, node$con): 'package:stats' may not be available when loading
    ## Warning in serialize(data, node$con): 'package:stats' may not be available when loading
    ## Warning in serialize(data, node$con): 'package:stats' may not be available when loading
    ## Warning in serialize(data, node$con): 'package:stats' may not be available when loading
    ## Warning in serialize(data, node$con): 'package:stats' may not be available when loading

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize = minSize, : For some of the
    ## pathways the P-values were likely overestimated. For such pathways log2err is set to NA.

    ## Warning in fgseaMultilevel(pathways = pathways, stats = stats, minSize = minSize, : For some pathways,
    ## in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.

    ## leading edge analysis...

    ## done...

``` r
reactome_results_df <- as.data.frame(reactome_results)

#ridge plot
#reactome_ridge_plot <- ridgeplot(reactome_results, showCategory = 10)

#gsea plot
png("figures/gseaplot.png", width = 1500, height = 900, res = 150)
reactome_gseaplot <- gseaplot2(reactome_results, c(1:3))
dev.off()
```

    ## png 
    ##   2

``` r
ggsave("figures/gseaplot.png", plot = reactome_gseaplot, width = 8, height = 6, dpi = 300)
```

``` r
#barplot
y <- arrange(reactome_results, desc(abs(NES)))%>% #order by descending NES values to compare between pathways
        group_by(sign(NES)) %>% #both positive and negative values for upregulated and downregulated pathways
        slice(1:5) #pathways with 5 highest NES score
```

    ## Warning: The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
    ## ℹ Please use the `.add` argument instead.
    ## ℹ The deprecated feature was likely used in the dplyr package.
    ##   Please report the issue at <https://github.com/tidyverse/dplyr/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

``` r
reactome_barplot <- ggplot(y, aes(NES, fct_reorder(Description, NES), fill=p.adjust), showCategory=10) + 
    geom_col(orientation='y') + 
    scale_fill_continuous(low='purple', high='blue', guide=guide_colorbar(reverse=TRUE)) +
    ggtitle("Barplot of 10 Pathways with Greatest Absolute Positive and Negative NES Value")+
    labs(y = "Pathways")+
    theme_minimal() + ylab(NULL)

ggsave("figures/reactome_barplot.png", plot = reactome_barplot, width = 16, height = 6, dpi = 300)
```

``` r
#library(KEGGREST)
#fit_test_names <- rownames(fit_test)
#gene_list <- c(fit_test[ ,1])
#names(gene_list) <- fit_test_names

#gene_list<- sort(gene_list, decreasing = TRUE)
#kegg_over <- gseKEGG(geneList = gene_list,
 #                    pvalueCutoff = 0.05,
  #                   organism = "hsa",
 #                    seed = TRUE)

#view(kegg_over)
```

``` r
#dotplot(kegg_over, showCategory = 10, title = "Kegg enrichment pathway dotplot")
#cnetplot(kegg_over , showCategory = 5, color.params =
#           list(FoldChange = NULL), title = "Kegg enrichment pathway cnetplot")
#heatplot(kegg_over, showCategory = 10)
#ridgeplot(kegg_over, showCategory = 10)

#gseaplot2(kegg_over, c(1,2,5))

# Convert KEGG GSEA results to a data frame
#kegg_results_df <- as.data.frame(kegg_over)
# Create a bar plot for KEGG enrichment results
#ggplot(kegg_results_df[1:10, ], aes(x = reorder(Description, NES), y = NES)) +
#  geom_bar(stat = "identity", fill = "steelblue") +
#  coord_flip() +
#  labs(title = "Top 10 KEGG Enrichment Pathways", x = "Pathway", y = "Normalized Enrichment Score (NES)") +
#  theme_minimal()
```
