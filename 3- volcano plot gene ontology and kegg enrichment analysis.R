library(ashr)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(ggtree)
library(enrichplot)
library(tidytree)
library(tidyr)
library(ggvenn)
library(stringr)

#ggvenn plot of upregulated and downregulated data
vennDat <- tibble(Geneid=rownames(annot.results.add)) %>% 
  dplyr::mutate('upregulated genes' = annot.results.add$padj < 0.05 &
                  !is.na(annot.results.add$padj) & 
                  annot.results.add$log2FoldChange > 0) %>%
  dplyr:: mutate('downregulated genes'= annot.results.add$padj< 0.05&
                   !is.na(annot.results.add$padj) & 
                   annot.results.add$log2FoldChange < 0)
ggvenn(vennDat, set_name_size = 4)

#Shrinking the log2FoldChange estimates towards zero and annotate our results.
shrink <- lfcShrink(ddsObj.b, res = results.additive, type = "ashr")
shrink.tab <- as.data.frame(shrink) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID") %>% 
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

#volcanoplot of all results
EnhancedVolcano(shrink.tab, x='logFC', lab = shrink.tab$Symbol,
                y= 'FDR', title = 'TG Vs WT wnt1 Volcano plot', 
                subtitle = 'differentially expressed genes',
                pCutoff = 0.05,
                FCcutoff = 1) +
  ggplot2::coord_cartesian(xlim=c(-10, 8), ylim = c(0,220)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-8,8,2))


#TreePlot of enrich GO analysis

#1 enrichGO analysis for 50 genes with FDR < 0.01 and |logFC| > 2  

Seg.results <- shrink.tab %>%
  drop_na(FDR) %>%  dplyr:: filter(FDR < 0.01 & abs(logFC) > 1)

KeggSegGen.1 <- shrink.tab %>%
  drop_na(FDR) %>%  dplyr:: filter(FDR < 0.01 & abs(logFC) > 1)%>%
  pull(GeneID)


universe <- shrink.tab$GeneID

ego <- enrichGO(gene = KeggSegGen.1,
                universe = universe,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05, 
                readable = TRUE)

#generate results as dataframe
cluster_summary <- data.frame(ego)

#choose the top 50 terms according to the lowest P-value
ego.1 <-  ego %>% top_n(50, wt=-p.adjust)

egox <- setReadable(ego.1, 'org.Mm.eg.db', 'ENTREZID')

#2 treeplot
ego.tree <- pairwise_termsim(egox)
treeplot(ego.tree, showCategory= 50)

#KEGG analysis of the samples
search_kegg_organism('mouse', by='common_name')

kegggenes <- shrink.tab %>% 
  drop_na(Entrez, FDR) %>% 
  filter(FDR < 0.05) %>% 
  pull(Entrez)

kegggenes <- as.character(kegggenes)
keggRes <- enrichKEGG(gene = kegggenes, organism = 'mmu')

keggRes$Description
keggRes$ID

#Plot pathway of wnt1
logFC <- shrink.tab$logFC
names(logFC) <- shrink.tab$Entrez
pathview(gene.data = logFC, 
         pathway.id = "mmu04310", 
         species= "mmu", limit = list(gene = 1, cpd=1))

#Plot pathway of Osteoclast Differentiation
pathview(gene.data = logFC, 
         pathway.id = "mmu04380", 
         species= "mmu", limit = list(gene = 1, cpd=1))


#heatplot for mutual pathways between p0-P14
#retrieve selected pathways from prepared data frame
KeggI <- read.csv("keggI.csv")
new.keggRes <- keggRes %>% dplyr::filter(ID%in%KeggI$ID)

keggx <- setReadable(new.keggRes,'org.Mm.eg.db', 'ENTREZID')

#create a data frame for logFC to choose the top 5-10 genes
lFC <- as.data.frame(logFC)
lFC$Symbol <- shrink.tab$Symbol
lFC$abs <- abs(logFC)

# Assuming you have a data frame keggx with a "geneID" column containing comma-separated gene symbols
# Split the geneID column by comma and convert each split part into a separate data frame

geneIDs <- as.character(keggx$geneID)
data_frames <- lapply(str_split(keggx@result[,"geneID"], "/"), as.data.frame)

top_10_by_abs <- function(df) {
  mutual_genes <- intersect(df[[1]], lFC$Symbol)
  
  if (length(mutual_genes) > 0) {
    # Filter the lFC data frame to get the absolute values for mutual genes
    filtered_lFC <- lFC %>%
      dplyr::filter(Symbol %in% mutual_genes)
    
    # Arrange the filtered data frame in descending order of absolute values
    filtered_lFC <- filtered_lFC %>%
      arrange(desc(abs))
    
    # Select the top 10 genes based on absolute values
    top_10_genes <- head(filtered_lFC, n = 10)
    
    return(top_10_genes)
  } else {
    return(NULL)  # Return NULL if there are no mutual genes
  }
}

# Apply the function to each data frame
top_10_data <- lapply(data_frames, top_10_by_abs)

# Combine the results into a character vector
random.core.genes <- sapply(top_10_data, function(df) {
  if (!is.null(df)) {
    paste(df$Symbol, collapse = "/")
  } else {
    NA  # Return NA if no matching genes were found
  }
})

#replace the genes included in each pathway for Keggx by the top 5-10 with highest |logFc| genes in that pathway

keggx@result$geneID <- random.core.genes
#plot the heatplot for the selected pathways and selceted genes
heatplot(keggx, foldChange=logFC, showCategory=20)



