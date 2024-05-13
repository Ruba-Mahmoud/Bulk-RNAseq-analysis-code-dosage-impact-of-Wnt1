library(DESeq2)
library(tidyverse)
library(tximport)
library(dplyr)

# Reading in the sample metadata
sampleinfo <- read_csv("SampleName.csv", col_types = c("cccc"))
sampleinfo %>% arrange(Status)

#Reading in the count data
files <- str_c("Salmon_O/", sampleinfo$SampleName, "/quant.sf")
files <- set_names(files, sampleinfo$SampleName)
tx2gene <- read_tsv("tx2gene.tsv")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

#Creating the design model formula
simple.model <- as.formula(~ Status)
model.matrix(simple.model,data = sampleinfo)


# set wildtype as the intercept:
sampleinfo <- mutate(sampleinfo, Status=fct_relevel(Status,"WT"))
model.matrix(simple.model,data = sampleinfo)
# create the DESeqDataSet object
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi, colData = sampleinfo, design = simple.model)

#Filter out  unexpressed genes
keep <- rowSums(counts(ddsObj.raw))>5
ddsObj.filt <- ddsObj.raw[keep,]
# build the DESeq from the raw data, the sample meta data 
ddsObj <- DESeq(ddsObj.filt)
results.simple <- results(ddsObj, alpha = 0.5)

sum(results.simple$padj < 0.05 , na.rm = T)


#since we have batch effect in P0 additive model was considered too:
additive.model <- as.formula(~ Batch + Status)
model.matrix(additive.model, data= sampleinfo)

# create the DESeqDataSet object
ddsObj.batch <- DESeqDataSetFromTximport(txi=txi, colData = sampleinfo, design = additive.model)
#Filter out unexpressed genes
keep <- rowSums(counts(ddsObj.batch)) > 5
ddsObj.b.filt <- ddsObj.batch[keep,]
# build the DESeq object for the additive model
ddsObj.b <- DESeq(ddsObj.b.filt)
results.additive <- results(ddsObj.b, alpha = 0.05)
sum(results.additive$padj < 0.05 , na.rm = T)

#Comparing two design models
ddsObj.LRT <- DESeq(ddsObj.b, test = "LRT", reduced = simple.model)
results.additive_vs_simple <- results(ddsObj.LRT)
results.additive_vs_simple
sum(results.additive_vs_simple$padj <0.05, na.rm=T)

