library(biomaRt)
library(tidyverse)
library(ensembldb)
library(ensembldb)

## set up connection to ensembl database

ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'mmusculus_gene_ensembl',
                      version = 105)

# Set the filter type and values
filterType <- "ensembl_gene_id"
# get the Ensembl IDs from our results table
filterValues <- rownames(results.additive)
# Set the list of attributes
attributeNames <- c("ensembl_gene_id",
                    "entrezgene_id",
                    "external_gene_name",
                    "description",
                    "gene_biotype",
                    "chromosome_name",
                    "start_position",
                    "end_position",
                    "strand",
                    "entrezgene_accession")

# run the query
annot.1 <- getBM(attributes=attributeNames,
                 filters = filterType,
                 values = filterValues,
                 mart = ensembl)

#check of ensembel gene id duplicates
annot.1 %>%  
  add_count(ensembl_gene_id) %>%  
  dplyr::filter(n>1) %>% 
  dplyr::count()
#check the duplicated genes info
annot.1 %>%  
  add_count(ensembl_gene_id) %>%  
  dplyr::filter(n>1) %>% 
  dplyr::select(ensembl_gene_id, external_gene_name, entrezgene_accession)


#Deduplicate the using the entrezgene_accession and the external_gene_name
dups <- annot.1 %>%  
  add_count(ensembl_gene_id) %>%  
  dplyr::filter(n>1)

fixedDuplicates <- dups %>% 
  dplyr::select(-n) %>% 
  dplyr::filter(entrezgene_accession==external_gene_name)

#check that you don't have any replicates 
fixedDuplicates %>%  
  add_count(ensembl_gene_id) %>%  
  dplyr::filter(n>1)

#Create a new annotation with all the unique entries from the full annotation plus the fixedDuplicates.
annot2 <- annot.1 %>%  
  add_count(ensembl_gene_id) %>%  
  dplyr::filter(n==1) %>% 
  dplyr::select(-n) %>% 
  bind_rows(fixedDuplicates)

nrow(annot2)

length(unique(annot.1$ensembl_gene_id))
# there is one gene remained
dups %>% 
  dplyr::filter(!ensembl_gene_id%in%annot2$ensembl_gene_id)

#choose one that is not in annot2
fixedDuplicates <- dups %>% 
  dplyr::filter(!ensembl_gene_id%in%annot2$ensembl_gene_id) %>% 
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::select(-n)

dim(fixedDuplicates)
#add the last gene to the annotation
annotUn <- bind_rows(annot2, fixedDuplicates)

#check that we have all unique ensembl gene id
nrow(annotUn)

length(unique(annot.1$ensembl_gene_id))

all(filterValues%in%annotUn$ensembl_gene_id)


#Check for duplicated Entrez IDs
annotUn %>% 
  dplyr::filter(!is.na(entrezgene_id)) %>% 
  add_count(entrezgene_id) %>% 
  dplyr::filter(n>1) %>% 
  dplyr::count(entrezgene_id)

#check our results to see if the gene are differentially expressed 
#or are important for the analysis
as.data.frame(results.additive) %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  left_join(annotUn, by = "ensembl_gene_id") %>% 
  dplyr::filter(!is.na(entrezgene_id)) %>% 
  dplyr::add_count(entrezgene_id) %>% 
  dplyr::filter(n>1) %>% 
  dplyr::select(ensembl_gene_id, entrezgene_id, padj, chromosome_name) %>% 
  dplyr::arrange(padj)

#select duplicated entrez gene id
dupEntrez <- annotUn %>% 
  add_count(entrezgene_id) %>% 
  dplyr::filter(n>1) %>% 
  dplyr::select(-n) %>% 
  dplyr::arrange(entrezgene_id, chromosome_name)
#after checking the genes, we will set the second entry to NA
dupEntrez$entrezgene_id[duplicated(dupEntrez$entrezgene_id)] <- NA

#merge annotun and fixed duplicates
annotFinal <- annotUn %>% 
  add_count(entrezgene_id) %>% 
  dplyr::filter(n==1) %>% 
  dplyr::select(-n) %>% 
  bind_rows(dupEntrez)

dim(annotFinal)

#check for duplicates again
annotFinal %>% 
  dplyr::filter(!is.na(entrezgene_id)) %>% 
  dplyr::add_count(entrezgene_id) %>% 
  dplyr::filter(n>1) 

all(filterValues%in%annotFinal$ensembl_gene_id)

#Final table
ensemblAnnot <- rownames(results.additive) %>%  
  enframe(name = NULL, value = "ensembl_gene_id")  %>%  
  left_join(annotFinal) %>%
  dplyr::select(GeneID="ensembl_gene_id", Entrez="entrezgene_id",
                Symbol="external_gene_name", Description="description",
                Biotype="gene_biotype", Chr="chromosome_name",
                Start="start_position", End="end_position",
                Strand="strand")


#extract GeneID and Symbol from the premade annotation include only our genes to a new table
symbol_annot <- ensemblAnnot %>% dplyr::filter(GeneID%in%rownames(results.additive))%>%
  dplyr::select("GeneID","Symbol")

head (symbol_annot)
dim(symbol_annot)
head(results.additive)
#bind the new table with our results 
annot.results.add <- cbind(symbol_annot, results.additive)

#in order to keep in the annotation our gene of interest
annot.add.b <- as.data.frame(results.additive) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)

