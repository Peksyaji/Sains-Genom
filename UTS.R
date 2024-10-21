library(GEOquery)

## Specify the GEO dataset ID
my_id <- "GSE33126"

## Download the dataset from GEO
gse <- getGEO(my_id)[[1]]

## Check the downloaded data object
gse

# Inspecting the Data
dim(gse)
pData(gse)
fData(gse)
exprs(gse)

# Normalization and Data Transformation
## Check the normalization and scales used in the data
summary(exprs(gse))

## Apply log2 transformation to the expression data
exprs(gse) <- log2(exprs(gse))

## Visualize the transformed data using a boxplot
boxplot(exprs(gse), outline = FALSE)

# Inspect clinical variables
library(dplyr)

## Extract relevant clinical information
sampleInfo <- pData(gse) %>%
  select(source_name_ch1, characteristics_ch1.1)

## Rename columns for convenience
sampleInfo <- rename(sampleInfo, group = source_name_ch1, patient = characteristics_ch1.1)
sampleInfo

# Sample clustering and PCA
library(pheatmap)

## Compute the correlation matrix and visualize it using a heatmap
corMatrix <- cor(exprs(gse), use = "c")
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix, annotation_col = sampleInfo)

## Perform PCA and plot the results
library(ggplot2)
library(ggrepel)

pca <- prcomp(t(exprs(gse)))
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y = PC2, col = group, label = paste("Patient", patient))) + 
  geom_point() + 
  geom_text_repel()

# Differential Expression Analysis
library(limma)
## Design matrix for differential expression analysis
design <- model.matrix(~0 + sampleInfo$group)
colnames(design) <- c("Normal", "Tumour")

## Perform the differential expression analysis using limma
fit <- lmFit(exprs(gse), design)
contrasts <- makeContrasts(Tumour - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

## Display the top results
topTable(fit2)

# Further Processing and Visualization of Differential Expression Results
## Annotate genes and prepare results for visualization
library(tibble)

anno <- select(fData(gse), Symbol, Entrez_Gene_ID, Chromosome, Cytoband)
fit2$genes <- anno
full_results <- topTable(fit2, number = Inf)
full_results <- rownames_to_column(full_results, "ID")

## Create a volcano plot for visualization
p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff) %>% 
  ggplot(aes(x = logFC, y = B, col = Significant)) + 
  geom_point()

# Exporting the Results
## Filter results for specific genes of interest
filter(full_results, Symbol == "SMOX")

## Filter results for genes related to TP53
filter(full_results, grepl("TP53", Symbol))

## Filter results for specific chromosomes
filter(full_results, Chromosome == 20)

## Filter results based on significance and fold change
filter(full_results, adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff)

# Heatmap of Selected Genes
## Visualize the top 20 most differentially expressed genes
topN <- 20

ids_of_interest <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)

gene_names <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Symbol)

gene_matrix <- exprs(gse)[ids_of_interest, ]
pheatmap(gene_matrix, labels_row = gene_names, scale = "row")

my_genes <- c("HIG2", "CA1","ETV4","FOXA1")
ids_of_interest <-  filter(full_results,Symbol %in% my_genes) %>% 
  pull(ID)

gene_names <-  filter(full_results,Symbol %in% my_genes) %>% 
  pull(Symbol)

# Gene Ontology Annotation
## Annotate the selected genes using Gene Ontology (GO)
library(org.Hs.eg.db)

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              columns = c("ENSEMBL", "GO"),
                              keys = my_genes,
                              keytype = "SYMBOL")
anno
