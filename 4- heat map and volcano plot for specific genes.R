library(tidyverse)
library(dplyr)
library(ComplexHeatmap)
library(limma)

#heatmap showing genes of canonical wnt pathway

#we filter rows for our target genes in the symbol column
Sigsymb <- shrink.tab %>% dplyr::filter(grepl('Wnt|Dvl|Dkk|Fzd|Lrp|Axin|Ctnnb1|Notum|Sost|Serpinf1', Symbol)& FDR < 0.05)
Sigsymb.h <- Sigsymb %>% pull(Symbol)
SygGene <- Sigsymb %>% pull (GeneID)

#normalize the  results for further analysis
plotDat <- vst( ddsObj.b)[SygGene,] %>% assay()
head(plotDat)
dim(plotDat)

#annotation and options for the heatmap
z.mat <- t(scale(t(plotDat), center = T, scale = T))
myPalette <- c("royalblue", "ivory", "red")
myRamp.1 <- circlize::colorRamp2(c(-2,0,2), myPalette)
ha1 = HeatmapAnnotation(df = colData(ddsObj.b)[,c("Status", "Batch")],
                        col= list(Status = c("WT" = "palegreen", "TG" = "darkgreen"), 
                                  Batch= c("I" = "lightblue" , "II" = "darkblue")))
row.names(z.mat) <- Sigsymb.h

#plot heatmap
Heatmap(z.mat, name = "z-score", col = myRamp.1, show_row_names = T,
        split = 3, column_title = "segnificant genes heatmap",
        rect_gp = gpar(col= "white", lwd=1),
        top_annotation = ha1)

#average heatmap for important ligands
#filter the selected ligands
Sigsymb.l <- shrink.tab %>% dplyr::filter(grepl('Wnt', Symbol)& FDR < 0.05 )
#Pull the Symbol of the genes, assign the symbols and the GeneID
Sigsymb.h.l <- Sigsymb.l %>% pull(Symbol)
SygGene.l <- Sigsymb.l %>% pull (GeneID)
#prepare the plotting data
plotDat.l <- vst( ddsObj.b)[SygGene.l,] %>% assay()
z.mat.l <- t(scale(t(plotDat.l), center = T, scale = T))
row.names(z.mat.l) <- Sigsymb.h.l
#plot the heatmap of the genes
Heatmap(z.mat.l, name = "z-score", col = myRamp.1, show_row_names = T,split = 3,
        column_title = "segnificant genes heatmap", rect_gp = gpar(col= "white", lwd=1),top_annotation = ha1)

#change the column names to calculate the average of the assay
colnames(z.mat.l)[1] = "Control"
colnames(z.mat.l)[2] = "Control"
colnames(z.mat.l)[3] = "Control"

colnames(z.mat.l)[4] = "P0 Wnt1"
colnames(z.mat.l)[5] = "P0 Wnt1"
colnames(z.mat.l)[6] = "P0 Wnt1"
ag.z.mat.l <- avearrays(z.mat.l)

#plot the average heatmap of the three replicates
Heatmap(ag.z.mat.l, name = "z-score", col = myRamp.1, show_row_names = T,split = 3, column_title = "ligands coding genes heatmap",rect_gp = gpar(col= "white", lwd=1))


#repeat the process for Receptor genes
Sigsymb.r <- shrink.tab %>% dplyr::filter(grepl('Fzd|Lrp5|Lrp6', Symbol) & FDR < 0.05)
Sigsymb.h.r <- Sigsymb.r %>% pull(Symbol)

SygGene.r <- Sigsymb.r %>% pull (GeneID)
plotDat.r <- vst( ddsObj.b)[SygGene.r,] %>% assay()
z.mat.r <- t(scale(t(plotDat.r), center = T, scale = T))
row.names(z.mat.r) <- Sigsymb.h.r
Heatmap(z.mat.r, name = "z-score", col = myRamp.1, show_row_names = T,split = 3, 
        column_title = "segnificant genes heatmap",rect_gp = gpar(col= "white", lwd=1),top_annotation = ha1)

#change the column names to calculate the average of the assay
colnames(z.mat.r)[1] = "Control"
colnames(z.mat.r)[2] = "Control"
colnames(z.mat.r)[3] = "Control"
colnames(z.mat.r)[4] = "P0 Wnt1"
colnames(z.mat.r)[5] = "P0 Wnt1"
colnames(z.mat.r)[6] = "P0 Wnt1"

ag.z.mat.r <- avearrays(z.mat.r)

Heatmap(ag.z.mat.r, name = "z-score", col = myRamp.1, show_row_names = T,split = 3, column_title = "Receptors coding genes heatmap",
        rect_gp = gpar(col= "white", lwd=1))

#repeat the process for the signaling genes
Sigsymb.s <- shrink.tab %>% dplyr::filter(grepl('Dvl|Axin|Ctnnb1', Symbol)& FDR < 0.05)
Sigsymb.h.s <- Sigsymb.s %>% pull(Symbol)
SygGene.s <- Sigsymb.s %>% pull (GeneID)
plotDat.s <- vst( ddsObj.b)[SygGene.s,] %>% assay()
z.mat.s <- t(scale(t(plotDat.s), center = T, scale = T))
row.names(z.mat.s) <- Sigsymb.h.s

Heatmap(z.mat.s, name = "z-score", col = myRamp.1, show_row_names = T,split = 3, column_title = "segnificant genes heatmap",
        rect_gp = gpar(col= "white", lwd=1),top_annotation = ha1)
colnames(z.mat.s)[1] = "Control"
colnames(z.mat.s)[2] = "Control"
colnames(z.mat.s)[3] = "Control"
colnames(z.mat.s)[4] = "P0 Wnt1"
colnames(z.mat.s)[5] = "P0 Wnt1"
colnames(z.mat.s)[6] = "P0 Wnt1"
ag.z.mat.s <- avearrays(z.mat.s)
Heatmap(ag.z.mat.s, name = "z-score", col = myRamp.1, show_row_names = T,
        split = 3, column_title = "signaling protein coding genes heatmap",rect_gp = gpar(col= "white", lwd=1))

#repeat the process for inhibitors:

#repeat the process for the signaling genes
Sigsymb.i <- shrink.tab %>% dplyr::filter(grepl('Dkk1|Dkk2|Dkk3|Sost|Serpinf1|Notum', Symbol)& FDR < 0.05)
Sigsymb.h.i <- Sigsymb.i %>% pull(Symbol)
SygGene.i <- Sigsymb.i %>% pull (GeneID)
plotDat.i <- vst( ddsObj.b)[SygGene.i,] %>% assay()
z.mat.i <- t(scale(t(plotDat.i), center = T, scale = T))
row.names(z.mat.i) <- Sigsymb.h.i

Heatmap(z.mat.i, name = "z-score", col = myRamp.1, show_row_names = T,split = 3, column_title = "segnificant genes heatmap",
        rect_gp = gpar(col= "white", lwd=1),top_annotation = ha1)
colnames(z.mat.i)[1] = "Control"
colnames(z.mat.i)[2] = "Control"
colnames(z.mat.i)[3] = "Control"
colnames(z.mat.i)[4] = "P0 Wnt1"
colnames(z.mat.i)[5] = "P0 Wnt1"
colnames(z.mat.i)[6] = "P0 Wnt1"
ag.z.mat.i <- avearrays(z.mat.i)
Heatmap(ag.z.mat.i, name = "z-score", col = myRamp.1, show_row_names = T,
        split = 3, column_title = "inhibitor proteins coding genes heatmap",rect_gp = gpar(col= "white", lwd=1))

#volcano plot of specific genes

library(EnhancedVolcano)
Sigsymb <- shrink.tab %>% dplyr::filter(grepl('Wnt|Dvl|Dkk1|Dkk2|Dkk3|Fzd|Lrp5|Lrp6|Axin|Ctnnb1|Notum|Sost|Serpinf1', Symbol))
Sigsymb.l <- shrink.tab %>% dplyr::filter(grepl('Wnt', Symbol))%>% pull(Symbol)
Sigsymb.r <- shrink.tab %>% dplyr::filter(grepl('Fzd|Lrp5|Lrp6', Symbol))%>% pull(Symbol)
Sigsymb.s <- shrink.tab %>% dplyr::filter(grepl('Dvl|Axin|Ctnnb1', Symbol))%>% pull(Symbol)
sigsymb.i <- shrink.tab %>% dplyr::filter(grepl('Dkk1|Dkk2|Dkk3|Notum|Serpinf1|Sost', Symbol)) %>% pull(Symbol)

keyvals.colour <- ifelse(Sigsymb$Symbol %in% Sigsymb.l , 'gold',
                         ifelse(Sigsymb$Symbol %in% Sigsymb.r, 'coral',
                                ifelse(Sigsymb$Symbol %in% sigsymb.i, 'darkred',
                                       ifelse('royalblue'))))


keyvals.colour[is.na(keyvals.colour)] <- 'royalblue'
names(keyvals.colour)[keyvals.colour == 'coral'] <-'receptors'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'signaling protein'
names(keyvals.colour)[keyvals.colour == 'darkred'] <- 'inhibitors'
names(keyvals.colour)[keyvals.colour == 'gold'] <-  'ligands'


EnhancedVolcano(Sigsymb, x='logFC', 
                lab = Sigsymb$Symbol[which(names(keyvals.colour) %in% c('receptors', 'ligands','signaling protein','inhibitors'))], 
                labSize = 6.0, colCustom = keyvals.colour,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 6,
                pointSize = 5,
                legendIconSize = 6.0,
                drawConnectors = TRUE,
                y= 'FDR', title = 'TG Vs WT wnt1 Volcano plot', 
                subtitle = 'differentially expressed genes',
                pCutoff = 0.05, FCcutoff = 0) 

