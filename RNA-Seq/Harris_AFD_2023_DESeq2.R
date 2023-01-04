#Using RNA-Seq per gene read counts from TRAP samples in Harris et al
#Counts and sample metadata are found in "raw_counts_all_samples.txt" and
#"meta_all_samples.txt"
#This code will generate the data sets used in Figure 1 and Figure S1, as
#well as the volcano plot, ridgeplots, PCA plot, and dot plot.

library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(AnnotationHub)
library(ensembldb)
library(clusterProfiler)
library(org.Ce.eg.db)
library(ggridges)

#Set up working directory
setwd('/Users/nathancsharris/Documents/Harris_AFD_DESeq2/')

## Load in data
data_all <- read.table("raw_counts_all_samples.txt", header=T, row.names=1)
meta_all <- read.table("meta_all_samples.txt", header=T, row.names=1)

AFD_samples <- c('AFD_25_1','AFD_15_1','AFD_25_2','AFD_15_2','AFD_25_3','AFD_15_3')
whole_animal_samples <- c('whole_25_1','whole_15_1','whole_25_2','whole_15_2')

data_AFD <- data_all[AFD_samples]
data_whole <- data_all[whole_animal_samples]

meta_AFD <- meta_all[AFD_samples,]
meta_whole <- meta_all[whole_animal_samples,]

## Create DESeq2Dataset object for PCA
dds_all <- DESeqDataSetFromMatrix(countData = data_all, colData = meta_all, design = ~ replicate + tissue + temperature)

#get normalized counts for later plots showing all data points
dds_all <- estimateSizeFactors(dds_all)
sizeFactors(dds_all)
normalized_counts_all <- counts(dds_all, normalized=TRUE)

### Transform counts for data visualization
rld <- rlog(dds_all, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup=c("tissue","temperature"))

## Create DESeq objects for each comparison
dds_AFD <- DESeqDataSetFromMatrix(countData = data_AFD, colData = meta_AFD, design = ~ replicate + temperature)
dds_whole <- DESeqDataSetFromMatrix(countData = data_whole, colData = meta_whole, design = ~ replicate + temperature)
dds_tissue <- DESeqDataSetFromMatrix(countData = data_all, colData = meta_all, design = ~ replicate + tissue)

## Run analyses

#dds_AFD compares expression at 25 vs 15 degrees for RNA extracted from AFD ribosome IP samples
dds_AFD <- DESeq(dds_AFD)

#dds_whole compares expression at 25 vs 15 degrees for RNA from whole animal lysates
dds_whole <- DESeq(dds_whole)

#dds_tissue compares expression in all AFD IP samples vs all whole animal lysate samples, independent of temperature
dds_tissue <- DESeq(dds_tissue)

res_AFD=results(dds_AFD)
res_whole=results(dds_whole)
res_tissue=results(dds_tissue)

##Extract results table, and shrink the log2 fold changes
res_table_AFD_apeglm <- lfcShrink(dds_AFD,coef='temperature_T_25_vs_T_15',type='apeglm')
write.csv(res_table_AFD_apeglm,'res_AFD_25vs15.csv')

res_table_whole_apeglm <- lfcShrink(dds_whole,coef='temperature_T_25_vs_T_15',type='apeglm')
write.csv(res_table_whole_apeglm,'res_whole_15vs25.csv')

dds_tissue$tissue <- relevel(dds_tissue$tissue, ref = 'whole_animal')
res_table_tissue_apeglm <- lfcShrink(dds_tissue,coef='tissue_whole_animal_vs_AFD', type='apeglm')
write.csv(res_table_tissue_apeglm,'res_tissue_AFD_vs_whole_animal.csv')

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 2

#turn results into tibble
res_table_AFD_tb <- res_table_AFD_apeglm %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_table_whole_tb <- res_table_whole_apeglm %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_table_tissue_tb <- res_table_tissue_apeglm %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


#filter to only include 0.05 significant and lfc > 2
sig_AFD <- res_table_AFD_tb %>%
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sig_whole <- res_table_whole_tb %>%
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

sig_tissue <- res_table_tissue_tb %>%
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

##order by pval and export
sig_AFD_ordered <- sig_AFD[order(sig_AFD$pvalue),]
write.csv(sig_AFD_ordered,'sig_AFD_ordered.csv')

sig_whole_ordered <- sig_whole[order(sig_whole$pvalue),]
write.csv(sig_whole_ordered,'sig_whole_ordered.csv')

sig_tissue_ordered <- sig_tissue[order(sig_tissue$pvalue),]
write.csv(sig_tissue_ordered,'sig_tissue_ordered.csv')

##remove genes that are temperature-regulated in both AFD and whole animal samples
sig_AFD_specific <- subset(sig_AFD, (!(sig_AFD$gene %in% sig_whole$gene)))
sig_AFD_specific_ordered <- sig_AFD_specific[order(sig_AFD_specific$pvalue),]
write.csv(sig_AFD_specific_ordered,'sig_AFD_specific.csv')


##Volcano plot code
volc_labels <- c('F08H9.4','F23D12.3','ins-39','droe-4','dac-1','C12D8.15',
                 'gcy-29','T25B6.4','gcy-18','zig-4')

keyvals <- ifelse(
  res_table_AFD_tb$log2FoldChange < -2 & res_table_AFD_tb$padj < 0.05, 'steelblue3',
  ifelse(res_table_AFD_tb$log2FoldChange > 2 & res_table_AFD_tb$padj < 0.05, 'indianred3',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'indianred3'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'steelblue3'] <- 'low'

EnhancedVolcano(res_table_AFD_tb,
                lab = res_table_AFD_tb$gene,
                selectLab = volc_labels,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-8, 8),
                ylim = c(0,20),
                titleLabSize = .1,
                title = '',
                subtitleLabSize = .1,
                subtitle = '',
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                colCustom = keyvals,
                FCcutoff = 2,
                pointSize = .5,
                labSize = 3,
                colAlpha = 1,
                cutoffLineType = 'dashed',
                #cutoffLineCol = 'red2',
                cutoffLineWidth = 0.5,
                #hline = c(10e-20, 10e-20 * 10e-30, 10e-20 * 10e-60,10e-20 * 10e-90),
                #hlineCol = c('black', 'black', 'black', 'black'),
                #hlineType = c('longdash', 'longdash', 'dotdash', 'dotdash'),
                #hlineWidth = c(0.4, 0.4, 0.8, 0.8),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'none',
                captionLabSize = .1,
                caption = '',
                drawConnectors = TRUE)

#Gene Set Enrichment
ah <- AnnotationHub()
elegans_ens <- query(ah, c("Caenorhabditis elegans", "EnsDb"))
elegans_ens <- elegans_ens[["AH95683"]]
genes(elegans_ens, return.type = "data.frame") %>% View()
annotations_ahb <- genes(elegans_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>% 
  dplyr::filter(gene_name %in% res_table_temperature_tb$gene)

annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()

non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)

annotations_ahb <- annotations_ahb[non_duplicates_idx, ]

#use non-shrunken results tables because overshrinkage of low replicate whole
#animal samples distorts LFCs, especially for genes with near 0 expression in one condition
res_table_AFD_unshrunken <- results(dds_AFD,contrast=contrast_AFD, alpha =0.05)
res_table_AFD_unshrunken_tb <- res_table_AFD_unshrunken %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_table_whole_unshrunken <- results(dds_whole,contrast=contrast_whole, alpha =0.05)
res_table_whole_unshrunken_tb <- res_table_whole_unshrunken %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_ids_AFD <- left_join(res_table_AFD_unshrunken_tb, annotations_ahb, by=c("gene"="gene_name"))
res_ids_whole <- left_join(res_table_whole_unshrunken_tb, annotations_ahb, by=c("gene"="gene_name"))

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
all_genes_AFD <- as.character(res_ids_AFD$gene_id)
all_genes_whole <- as.character(res_ids_whole$gene_id)

## Extract significant results
sig_AFD_genes <- dplyr::filter(res_ids_AFD, padj < 0.05, abs(log2FoldChange) > 2)
sig_whole_genes <- dplyr::filter(res_ids_whole, padj < 0.05, abs(log2FoldChange) > 2)

##GSEA
## Remove any NA values (reduces the data by quite a bit)
res_entrez_AFD <- dplyr::filter(res_ids_AFD, entrezid != "NA")
res_entrez_whole <- dplyr::filter(res_ids_whole, entrezid !="NA")

## Remove any Entrez duplicates
res_entrez_AFD <- res_entrez_AFD[which(duplicated(res_entrez_AFD$entrezid) == F), ]
res_entrez_whole <- res_entrez_whole[which(duplicated(res_entrez_whole$entrezid) ==F),]

res_entrez_AFD$CELE <- paste("CELE", res_entrez_AFD$gene, sep="_")
res_entrez_whole$CELE <- paste("CELE", res_entrez_whole$gene, sep="_")

foldchanges_AFD <- res_entrez_AFD$log2FoldChange
names(foldchanges_AFD) <- res_entrez_AFD$CELE

foldchanges_whole <- res_entrez_whole$log2FoldChange
names(foldchanges_whole) <- res_entrez_whole$CELE

## Sort fold changes in decreasing order
foldchanges_AFD <- sort(foldchanges_AFD, decreasing = TRUE)
foldchanges_whole <- sort(foldchanges_whole, decreasing = TRUE)

set.seed(123456)

## GSEA using gene sets from GO terms
foldchanges_entrez_AFD <- res_entrez_AFD$log2FoldChange
foldchanges_entrez_whole <- res_entrez_whole$log2FoldChange
names(foldchanges_entrez_AFD) <- res_entrez_AFD$entrezid
names(foldchanges_entrez_whole) <- res_entrez_whole$entrezid
foldchanges_entrez_AFD <- sort(foldchanges_entrez_AFD, decreasing = TRUE)
foldchanges_entrez_whole <- sort(foldchanges_entrez_whole, decreasing = TRUE)

gseaGO_AFD <- gseGO(geneList = foldchanges_entrez_AFD, # ordered named vector of fold changes (Entrez IDs are the associated names)
                org.Ce.eg.db,
                ont = "MF",
                #keytype = "ENTREZID",
                #nPerm = 1000, # default number permutations
                minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                pvalueCutoff = 0.05, # padj cutoff value
                verbose = FALSE)
gseaGO_results_AFD <- gseaGO_AFD@result

write.csv(gseaGO_results_AFD, "gsea_GO_AFD_MF_normal_shrink.csv", quote=F)

ridgeplot(gseaGO_AFD)

gseaGO_whole <- gseGO(geneList = foldchanges_entrez_whole, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    org.Ce.eg.db,
                    ont = "MF",
                    #keytype = "ENTREZID",
                    #nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)
gseaGO_results_whole <- gseaGO_whole@result

write.csv(gseaGO_results_whole, "gsea_GO_whole_MF_no_shrink.csv", quote=F)

ridgeplot(gseaGO_whole)

## Enrichment analysis using genes from "the neuronal genome of C elegans"
TERM2GENE <-read.table("term_to_name_neuronal_genome_no_duplicates_no_TFs.csv", header=F, sep = ",")

sig_genes_gene_name_AFD <- as.character(sig_AFD$gene)
ego_neuronal_AFD <- enricher(gene = sig_genes_gene_name_AFD,
                         TERM2GENE = TERM2GENE,
                         #universe = all_genes,
                         #keyType = "ENSEMBL",
                         #OrgDb = org.Ce.eg.db, 
                         #ont = "MF", 
                         #pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05) 
#readable = TRUE)

cluster_summary_neuronal <- data.frame(ego_neuronal_AFD)
write.csv(cluster_summary_neuronal, "neuronal_genome_GO_AFD_010323.csv", quote=F)
dotplot(ego_neuronal_AFD)

sig_genes_gene_name_whole <- as.character(sig_whole$gene)
ego_neuronal_whole <- enricher(gene = sig_genes_gene_name_whole,
                               TERM2GENE = TERM2GENE,
                               #universe = all_genes,
                               #keyType = "ENSEMBL",
                               #OrgDb = org.Ce.eg.db, 
                               #ont = "MF", 
                               #pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05) 

cluster_summary_neuronal_whole <- data.frame(ego_neuronal_whole)