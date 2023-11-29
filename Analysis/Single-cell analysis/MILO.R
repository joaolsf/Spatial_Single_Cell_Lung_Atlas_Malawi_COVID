
# Reference: https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#2_Visualize_the_data
# https://github.com/MarioniLab/miloR


# 1- Converting anndata to SingleCellExperiment object (SCE)
library(anndata)
library(dplyr)
library(patchwork)
adata <- anndata::read_h5ad('adata_subset2')

library(rhdf5)
library(zellkonverter)
#ad <- readH5AD('adata_subset2') #it worked and converted to SingleCell Experiment
#split the adata object by disease group to run the analysis within each group and compare later?
ad <- readH5AD('adata_COVID')

# Milo analysis for cohorts comparisons
ad <- readH5AD('adata_Malawi_Brazil_US_inner_integrated_covid')

# 2- Save SCE object
saveRDS(ad, "adata_subset2.rds")
ad <- readRDS("adata_subset2.rds")
ad

saveRDS(ad, "adata_COVID.rds")

saveRDS(ad, "adata_Malawi_Brazil_US_inner_integrated_covid.rds")
ad <- readRDS("adata_Malawi_Brazil_US_inner_integrated_covid.rds")
ad

library(SingleCellExperiment)
names(assays(ad)) <- "counts"

library(scater)
ad <- logNormCounts(ad)

#if the Error in .local(x, ...) : size factors should be positive rises
ad <- ad[, colSums(counts(ad)) > 0]
ad <- logNormCounts(ad)
ad

set.seed(1234)
ad <- runPCA(ad)

saveRDS(ad, "adata_subset2.rds")
saveRDS(ad, "adata_Malawi_Brazil_US_inner_integrated_covid.rds")
saveRDS(ad, "adata_COVID.rds")

#library(scRNAseq)

# 3- Milo set-up and make milo object
library(miloR)
library(igraph)

milo <- Milo(ad)
milo
saveRDS(milo, "milo_clinical_groups.rds")
saveRDS(milo, "milo_HIV.rds")
milo <- readRDS("./milo_clinical_groups.rds")

saveRDS(milo, "milo_cohorts.rds")
milo <- readRDS("./milo_cohorts.rds")

# 4- Add KNN graph if there is no KNN graph in the SCE object 

# Milo looks for neighbourhoods in a KNN graph to perform DA analysis. 
# This need to be stored in the graph slot of the Milo object. Here we have two options:

# A- We can add the KNN graph precomputed with sc.pp.neighbors, using the function buildFromAdjacency. 
# IN SCANPY (PYTHON) RUN THE FOLLOWING CODE:
# Save the binary connectivity matrix
# knn_adjacency = adata.obsp["connectivities"]
# Load in R and run: milo_graph <- buildFromAdjacency(knn_adjacency, k=20, is.binary=TRUE)
# graph(milo) <- miloR::graph(milo_graph)

# B- we can recompute the KNN graph using the dedicated function in miloR. USE THIS!
milo <- buildGraph(milo, k=100, d = 37, transposed = FALSE) #k = 100 because is the same k used in sc.pp.neighbors

milo <- buildGraph(milo, k=100, d = 16, transposed = FALSE) #k = 100 because is the same k used in sc.pp.neighbors

# 5- Defining representative neighbourhoods on the KNN graph

# We define the neighbourhood of a cell, the index, as the group of cells connected by an edge in the KNN graph to the index cell. For efficiency, we don’t test for DA in the neighbourhood of every cell, but we sample as indices a subset of representative cells, using a KNN sampling algorithm used by Gut et al. 2015.
# As well as d and k, for sampling we need to define a few additional parameters:
# prop: the proportion of cells to randomly sample to start with. We suggest using prop=0.1 for datasets of less than 30k cells. 
# For bigger datasets using prop=0.05 should be sufficient (and makes computation faster).
# refined: indicates whether you want to use the sampling refinement algorith, or just pick cells at random. 
# The default and recommended way to go is to use refinement. 
# The only situation in which you might consider using random instead, is if you have batch corrected your data with a graph based correction algorithm, such as BBKNN, but the results of DA testing will be suboptimal.
milo <- makeNhoods(milo, prop = 0.05, k = 100, d=37, refined = TRUE, reduced_dims = "PCA")

milo <- makeNhoods(milo, prop = 0.05, k = 100, d=16, refined = TRUE, reduced_dims = "PCA")

# Once we have defined neighbourhoods, we plot the distribution of neighbourhood sizes (i.e. how many cells form each neighbourhood) to evaluate whether the value of k used for graph building was appropriate. 
# We can check this out using the plotNhoodSizeHist function.
# As a rule of thumb we want to have an average neighbourhood size over 5 x N_samples. 
# If the mean is lower, or if the distribution is
plotNhoodSizeHist(milo)

# 6- Counting cells in neighbourhoods

# Milo leverages the variation in cell numbers between replicates for the same experimental condition to test for differential abundance. 
# Therefore we have to count how many cells from each sample are in each neighbourhood. 
# We need to use the cell metadata and specify which column contains the sample information.
colData(milo) # to check the metadata

# This adds to the Milo object a n×m matrix, where n is the number of neighbourhoods and m is the number of experimental samples. 
# Values indicate the number of cells from each sample counted in a neighbourhood. 
# This count matrix will be used for DA testing.

# Count cells per ROI
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="ROI")
# Count cells per Patient
#milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="Patient")
# Count cells per Disease group
#milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="Group")

head(nhoodCounts(milo))

# 7- Defining experimental design

# Now we are all set to test for differential abundance in neighbourhoods. 
# We implement this hypothesis testing in a generalized linear model (GLM) framework, specifically using the Negative Binomial GLM implementation in edgeR.

# We first need to think about our experimental design. 
# The design matrix should match each sample to the experimental condition of interest for DA testing. 
# In this case, we want to detect DA between disease progression, stored in the stage column of the dataset colData. 
# We also include the patient ID (Case) column in the design matrix. This represents a known technical covariate that we want to account for in DA testing.

covid_design <- data.frame(colData(milo))[,c("ROI","HIV")]

## Convert batch info from integer to factor
covid_design$HIV <- as.factor(covid_design$HIV) 
covid_design <- distinct(covid_design)

# remove duplicated ROIs
library(dplyr)
covid_design2 <- covid_design %>% distinct(ROI, .keep_all = TRUE)
rownames(covid_design2) <- covid_design2$ROI
covid_design2

# 8- Computing neighbourhood connectivity

# Milo uses an adaptation of the Spatial FDR correction introduced by cydar, where we correct p-values accounting for the amount of overlap between neighbourhoods. 
# Specifically, each hypothesis test P-value is weighted by the reciprocal of the kth nearest neighbour distance. 
# To use this statistic we first need to store the distances between nearest neighbors in the Milo object. 
# This is done by the calcNhoodDistance function (N.B. this step is the most time consuming of the analysis workflow and might take a couple of minutes for large datasets).
milo <- calcNhoodDistance(milo, d=37, reduced.dim = "PCA") #it takes a long time to run!
milo

saveRDS(milo, "milo_cohorts.rds")

saveRDS(milo, "milo_clinical_groups.rds")
milo <- readRDS("milo_clinical_groups.rds")

saveRDS(milo, "milo_HIV.rds")

# 9- DA Testing
# Now we can do the DA test, explicitly defining our experimental design. 
# In this case, we want to test for differences between experimental stages, while accounting for the variability between technical batches.
# This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates wheather there is significant differential abundance between groups. 
# The main statistics we consider here are:
# logFC: indicates the log-Fold change in cell numbers between ED and LD;
# PValue: reports P-values before FDR correction
# SpatialFDR: reports P-values corrected for multiple testing accounting for overlap between neighbourhoods

## Reorder rownames to match columns of nhoodCounts(milo)
table(covid_design2$Group)
#covid_design2 <- rownames(covid_design2$ROI)
covid_design2 <- covid_design2[colnames(nhoodCounts(milo)), , drop=FALSE]
write.csv(covid_design2, "covid_design_group.csv")

p <- colnames(nhoodCounts(milo))
write.csv(p, "samples.csv")

covid_design2 <- read.csv("covid_design_group.csv")
rownames(covid_design2) <- covid_design2$ROI
covid_design2 <- covid_design2[,-1]

contrast.1 <- c("GroupCOVID19 - GroupPneumonia") # the syntax is <VariableName><ConditionLevel> - <VariableName><ControlLevel>
contrast.2 <- c("bnh - GroupNon_Pneumonia")
contrast.3 <- c("HIVHIV_Pos - HIVHIV_Neg")

da_results1 <- testNhoods(milo, design = ~ 0 + Group, design.df = covid_design2, model.contrasts = contrast.1)
da_results2 <- testNhoods(milo, design = ~ 0 + Group, design.df = covid_design2, model.contrasts = contrast.2)
da_results3 <- testNhoods(milo, design = ~ 0 + HIV, design.df = covid_design2, model.contrasts = contrast.3)

head(da_results1)
da_results1 %>%
  arrange(SpatialFDR) %>%
  head() 

saveRDS(da_results1, "da_results1_COVID_Pneumonia.rds")
saveRDS(da_results2, "da_results2COVID_NonPneumonia.rds")
saveRDS(da_results3, "da_results3_HIV.rds")

da_results2 <- readRDS("da_results2COVID_NonPneumonia.rds")

da_results1 <- readRDS("./DA_results1/da_results1.rds")

# 9.1- Inspecting DA testing results

# We can start inspecting the results of our DA analysis from a couple of standard diagnostic plots. 
# We first inspect the distribution of uncorrected P values, to verify that the test was balanced.
ggplot(da_results1, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results2, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results3, aes(PValue)) + geom_histogram(bins=50)

# Then we visualize the test results with a volcano plot (remember that each point here represents a neighbourhood, not a cell).
ggplot(da_results1, aes(logFC, -log10(SpatialFDR))) + #or SpatialFDR
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

# To visualize DA results relating them to the embedding of single cells, we can build an abstracted graph of neighbourhoods that we can superimpose on the single-cell embedding.
# Here each node represents a neighbourhood, while edges indicate how many cells two neighbourhoods have in common. 
# Here the layout of nodes is determined by the position of the index cell in the UMAP embedding of all single-cells. 
# The neighbourhoods displaying singificant DA are colored by their log-Fold Change.

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = "X_umap", colour_by="Group", text_by = "pheno_cluster", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

umap_pl <- plotReducedDim(milo, dimred = "X_umap", colour_by="HIV", text_by = "pheno_cluster", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

umap_pl <- plotReducedDim(milo, dimred = "X_umap", colour_by="Cohort", text_by = "pheno_cluster_edited2", 
                         text_size = 3, point_size=0.5) +
  guides(fill="none")


## Plot neighbourhood graph
milo <- buildNhoodGraph(milo)

nh_graph_pl <- plotNhoodGraphDA(milo, da_results1, layout="X_umap",alpha=0.1) 

nh_graph_pl2 <- plotNhoodGraphDA(milo, da_results2, layout="X_umap",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")

umap_pl + nh_graph_pl2 +
  plot_layout(guides="collect")

# We might also be interested in visualizing wheather DA is particularly evident in certain cell types. 
# To do this, we assign a cell type label to each neighbourhood by finding the most abundant cell type within cells in each neighbourhood. 
# We can label neighbourhoods in the results data.frame using the function annotateNhoods. 
# This also saves the fraction of cells harbouring the label.
da_results1 <- annotateNhoods(milo, da_results1, coldata_col = "pheno_cluster_new")
da_results2 <- annotateNhoods(milo, da_results2, coldata_col = "pheno_cluster_new")
da_results3 <- annotateNhoods(milo, da_results3, coldata_col = "pheno_cluster_new")

da_results1 <- annotateNhoods(milo, da_results1, coldata_col = "pheno_cluster_edited2")
da_results2 <- annotateNhoods(milo, da_results2, coldata_col = "pheno_cluster_edited2")


#reorder clusters
da_results1$pheno_cluster_ordered <- factor(da_results1$pheno_cluster_edited,levels=c("CD3+ cell", "Endothelial cell", "Epithelial cell", "Alveolar Macrophage",
                                                                                      "Apoptotic Neutrophil", "Apoptotic Smooth Muscle Cell", "Apoptotic Fibroblast", "Apoptotic Epithelial cell",
                                                                                      "Smooth Muscle Cell", "CD66bLow Neutrophil", "CD66bHigh Neutrophil", "ColHigh Fibroblast", "ColLow Fibroblast", "CD8 T cell",
                                                                                      "CD4 T cell", "Interstitial Macrophage", "Classical Monocyte", "Activated Endothelial cell",
                                                                                      "SARSCoV2+ Epithelial cell", "SARSCoV2+ Alveolar Macrophage"))
library(forcats)
da_results1$pheno_cluster_ordered <- fct_rev(da_results1$pheno_cluster_ordered)


plotDAbeeswarm(da_results1, group.by = "pheno_cluster_ordered", alpha = 0.1)
#plot only neighbourhoods with spatialFDR <0.1
plotDAbeeswarm(da_results1, group.by = "pheno_cluster_new", alpha = 0.1, subset.nhoods = da_results1$SpatialFDR < 0.05) + theme_bw(base_size = 12)
plotDAbeeswarm(da_results2, group.by = "pheno_cluster_new", alpha = 0.1, subset.nhoods = da_results2$SpatialFDR < 0.05) + theme_bw(base_size = 12)
plotDAbeeswarm(da_results3, group.by = "pheno_cluster_new", alpha = 0.1, subset.nhoods = da_results3$SpatialFDR < 0.05) + theme_bw(base_size = 12)


# Cohorts: Malawi vs Brazil
plotDAbeeswarm(da_results1, group.by = "pheno_cluster_edited2", alpha = 0.1, subset.nhoods = da_results1$SpatialFDR < 0.05) + theme_bw(base_size = 12)
# Cohorts: Malawi vs US
plotDAbeeswarm(da_results2, group.by = "pheno_cluster_edited2", alpha = 0.1, subset.nhoods = da_results1$SpatialFDR < 0.05) + theme_bw(base_size = 12)


# While neighbourhoods tend to be homogeneous, we can define a threshold for celltype_fraction to exclude neighbourhoods that are a mix of cell types.
ggplot(da_results1, aes(pheno_cluster_edited_fraction)) + geom_histogram(bins=50)

da_results1$pheno_cluster_edited <- ifelse(da_results1$pheno_cluster_edited_fraction < 0.7, "Mixed", da_results1$pheno_cluster_edited)

# Now we can visualize the distribution of DA Fold Changes in different cell types after filtering
plotDAbeeswarm(da_results1, group.by = "pheno_cluster_edited")

# 10- Finding markers of DA populations
# Once you have found your neighbourhoods showindg significant DA between conditions, 
# you might want to find gene signatures specific to the cells in those neighbourhoods.
# The function findNhoodGroupMarkers runs a one-VS-all differential gene expression test to identify marker genes for a group of neighbourhoods of interest. 
# Before running this function you will need to define your neighbourhood groups depending on your biological question, that need to be stored as a NhoodGroup column in the da_results data.frame.

da_results1$NhoodGroup <- as.numeric(da_results1$SpatialFDR < 0.1 & da_results1$logFC < 0)
da_nhood_markers1 <- findNhoodGroupMarkers(milo, da_results1, subset.row = rownames(milo)[1:38])
head(da_nhood_markers1)

plotNhoodGroups(milo, da_results1, layout="X_umap") 


ggplot(da_nhood_markers1, aes(logFC_0,-log10(adj.P.Val_0 ))) + 
  geom_point() +
  geom_hline(yintercept = 2)

rownames(da_nhood_markers1) <- da_nhood_markers1$GeneID
markers <- rownames(da_nhood_markers1)[da_nhood_markers1$adj.P.Val_0 < 0.01 & da_nhood_markers1$logFC_0 > 0]
plotNhoodExpressionGroups(milo, da_results1, features=markers,
                          subset.nhoods = da_results1$NhoodGroup %in% c('0','1'), scale=FALSE,
                          grid.space = "fixed")

##############################################################################################################################################################################################################################################################################################################################################################

# PLOTTING FREQUENCY BAR PLOTS IN SEURAT WITH GGPLOT

library(Seurat)
library(scRNAseq)
library(SingleCellExperiment)
library(ggplot2)

# Plotting immune cell type proportions

library(rhdf5)
library(zellkonverter)
#ad <- readH5AD('adata_subset2') #it worked and converted to SingleCell Experiment
#split the adata object by disease group to run the analysis within each group and compare later?
ad <- readH5AD('ad_covid_immune')

#lung.immune <-ad[,ad$pheno_cluster_new %in% c('B cell', 'CD11c+ cell', 'EM CD4 T cell', 'Proliferative CD4 T cell', 'CD4 Treg cell',
 #                                           'CD8 T cell', 'EM CD8 T cell', 'Proliferative CD8 T cell' , 'CD3+ cell',  
  #                                           'SARSCoV2+ ArginaseHighVISTAHigh Activated Neutrophil', 
   #                                           'ArginaseHighVISTAHigh Activated Neutrophil', 'ArginaseLowVISTALow Activated Neutrophil', 
    #                                          'ArginaseLowVISTALow Neutrophil', 
     #                                        'Classical Monocyte', 'SARSCoV2+ Interstitial Macrophage', 'Proliferating Interstitial Macrophage', 'Interstitial Macrophage', 
      #                                       'Apoptotic SARSCoV2+  Alveolar Macrophage', 'Apoptotic Alveolar Macrophage', 'Alveolar Macrophage')]

freq.table.group <- as.data.frame(prop.table(x = table(ad$pheno_cluster_new, ad$HIV), margin = 2))



freq.table.group[freq.table.group==0] <- NA
freq.table.group<-freq.table.group[complete.cases(freq.table.group),]
write.csv(freq.table.group, "frequency_immune_groups.csv")

# to plot the reversed order of cells
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c('B cell', 'CD11c+ cell', 'EM CD4 T cell', 'Proliferative CD4 T cell', 'CD4 Treg cell',
                                                    'CD8 T cell', 'EM CD8 T cell', 'Proliferative CD8 T cell' , 'CD3+ cell',  
                                                    'SARSCoV2+ ArginaseHighVISTAHigh Activated Neutrophil', 
                                                    'ArginaseHighVISTAHigh Activated Neutrophil', 'ArginaseLowVISTALow Activated Neutrophil', 
                                                    'ArginaseLowVISTALow Neutrophil', 
                                                    'Classical Monocyte', 'SARSCoV2+ Interstitial Macrophage', 'Proliferating Interstitial Macrophage', 'Interstitial Macrophage', 
                                                    'Apoptotic SARSCoV2+  Alveolar Macrophage', 'Apoptotic Alveolar Macrophage', 'Alveolar Macrophage')))
# to plot the normal order of cells
library(forcats)
freq.table.group$Var1 <- fct_rev(freq.table.group$Var1)

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("Non_Pneumonia", "COVID-19", "Pneumonia"))))

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("HIV_Neg", "HIV_Pos"))))


celltype.colours <- setNames(c("#000000","#b5bbe3","#ff46a1","#ffe6f2", "#A64D79",
                              "#6fa8dc",'#4d7191',"#8faac2", "#d9d2e9",
                              "#336600","#8dd593","#0fcfc0","#9cded6",
                              "#b86cb9","#00B0F0","#4900EF","#4a6fe3",
                              "#FFD966","#8595e1","#A381EF"
),
c('B cell', 'CD11c+ cell', 'EM CD4 T cell', 'Proliferative CD4 T cell', 'CD4 Treg cell',
  'CD8 T cell', 'EM CD8 T cell', 'Proliferative CD8 T cell' , 'CD3+ cell',  
  'SARSCoV2+ ArginaseHighVISTAHigh Activated Neutrophil', 
  'ArginaseHighVISTAHigh Activated Neutrophil', 'ArginaseLowVISTALow Activated Neutrophil', 
  'ArginaseLowVISTALow Neutrophil', 
  'Classical Monocyte', 'SARSCoV2+ Interstitial Macrophage', 'Proliferating Interstitial Macrophage', 'Interstitial Macrophage', 
  'Apoptotic SARSCoV2+  Alveolar Macrophage', 'Apoptotic Alveolar Macrophage', 'Alveolar Macrophage'))


ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Group", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = rev(levels(freq.table.group$Var2))) + 
  scale_fill_manual(values=celltype.colours) + theme_minimal(base_size = 14)

# Plotting stromal cell type proportions

ad <- readH5AD('ad_subset2_stromal')

freq.table.group <- as.data.frame(prop.table(x = table(ad$pheno_cluster_new, ad$Group), margin = 2))
freq.table.group[freq.table.group==0] <- NA
freq.table.group<-freq.table.group[complete.cases(freq.table.group),]
write.csv(freq.table.group, "frequency_immune_groups.csv")

# to plot the reversed order of cells
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c('Activated Endothelial cell', 'Endothelial cell', 'Proliferative Endothelial cell', 'Fibroblast', 'Proliferative Fibroblast', 'Apoptotic Fibroblast',  
                                                    'Smooth Muscle cell', 'SARSCoV2+ Epithelial cell', 'Epithelial cell' , 'Proliferative Epithelial cell', 'SARSCoV2+ AT2 cell', 'AT2 cell')))
# to plot the normal order of cells
library(forcats)
freq.table.group$Var1 <- fct_rev(freq.table.group$Var1)

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("Non_Pneumonia", "COVID-19", "Pneumonia"))))

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("HIV_Neg", "HIV_Pos"))))


celltype.colours <- setNames(c("#0070C0","#e07b91","#fce5cd", "#b86cb9","#d33f6a",
                               "#FF7DA8",'#FA8000',"#b9877d", "#f0b98d",
                               "#f4cccc","#A53F02","#d6bcc0"
),
c('Activated Endothelial cell', 'Endothelial cell', 'Proliferative Endothelial cell', 'Fibroblast', 'Proliferative Fibroblast', 'Apoptotic Fibroblast',  
  'Smooth Muscle cell', 'SARSCoV2+ Epithelial cell', 'Epithelial cell' , 'Proliferative Epithelial cell', 'SARSCoV2+ AT2 cell', 'AT2 cell'))

# to plot the normal order of cells
ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Group", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = rev(levels(freq.table.group$Var2))) + 
  scale_fill_manual(values=celltype.colours) + theme_minimal(base_size = 14)

# Plotting immune cell type proportions from integrated cohorts

ad <- readH5AD('adata_Malawi_Brazil_US_inner_integrated_covid_immune')

freq.table.group <- as.data.frame(prop.table(x = table(ad$pheno_cluster_edited2, ad$Cohort), margin = 2))
freq.table.group[freq.table.group==0] <- NA
freq.table.group<-freq.table.group[complete.cases(freq.table.group),]
write.csv(freq.table.group, "frequency_immune_groups.csv")

# to plot the reversed order of cells
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c('B cell', "CD4 T cell", "CD4 Treg cell", "CD8 T cell",
                                                    'Dendritic cell',  'Mast cell',
                                                    "SARSCoV2+ NK cell",  'NK cell',
                                                    'SARSCoV2+ neutrophil', "Apoptotic neutrophil", 'Neutrophil', 
                                                    'SARSCoV2+ monocyte', "Classical monocyte", 'SARSCoV2+ IM', 'Interstitial macrophage',
                                                    'SARSCoV2+ AM', 'Apoptotic alveolar macrophage', 'Alveolar macrophage')))
# to plot the normal order of cells
library(forcats)
freq.table.group$Var1 <- fct_rev(freq.table.group$Var1)

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("Brazil", "Malawi", "US"))))

celltype.colours <- setNames(c("#000000","#ff46a1","#A64D79","#6fa8dc",
                               "#b5bbe3",'#d33f6a',"#d9d2e9", "#d9d2e9",
                               "#336600","#8dd593","#11c638","#FF7DA8",
                               "#b86cb9","#00B0F0","#4a6fe3",
                               "#FFD966","#8595e1","#A381EF"
),
c('B cell', "CD4 T cell", "CD4 Treg cell", "CD8 T cell",
  'Dendritic cell',  'Mast cell',
  "SARSCoV2+ NK cell",  'NK cell',
  'SARSCoV2+ neutrophil', "Apoptotic neutrophil", 'Neutrophil', 
  'SARSCoV2+ monocyte', "Classical monocyte", 'SARSCoV2+ IM', 'Interstitial macrophage',
  'SARSCoV2+ AM', 'Apoptotic alveolar macrophage', 'Alveolar macrophage'))


ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Group", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = rev(levels(freq.table.group$Var2))) + 
  scale_fill_manual(values=celltype.colours) + theme_minimal(base_size = 14)

# Plotting stromal cell type proportions

ad <- readH5AD('adata_Malawi_Brazil_US_inner_integrated_covid_stromal')

freq.table.group <- as.data.frame(prop.table(x = table(ad$pheno_cluster_edited2, ad$Cohort), margin = 2))
freq.table.group[freq.table.group==0] <- NA
freq.table.group<-freq.table.group[complete.cases(freq.table.group),]
write.csv(freq.table.group, "frequency_immune_groups.csv")

# to plot the reversed order of cells
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c('SARSCoV2+ epithelial cell', 'Apoptotic epithelial cell', 'Epithelial cell', 
                                                    'Activated endothelial cell', 'Endothelial cell', 'Fibroblast', 'Apoptotic fibroblast',
                                                    'SMC', 'Apoptotic SMC', "Mesenchymal")))
# to plot the normal order of cells
library(forcats)
freq.table.group$Var1 <- fct_rev(freq.table.group$Var1)

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("Brazil", "Malawi", "US"))))


celltype.colours <- setNames(c("#A53F02","#0fcfc0","#f0b98d", "#0070C0","#e07b91",
                               "#b86cb9",'#FF7DA8',"#FA8000", "#9cded6",
                               "#ffe6f2"
),
c('SARSCoV2+ epithelial cell', 'Apoptotic epithelial cell', 'Epithelial cell', 
  'Activated endothelial cell', 'Endothelial cell', 'Fibroblast', 'Apoptotic fibroblast',
  'SMC', 'Apoptotic SMC', "Mesenchymal"))

# to plot the normal order of cells
ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Group", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = rev(levels(freq.table.group$Var2))) + 
  scale_fill_manual(values=celltype.colours) + theme_minimal(base_size = 14)

##############################################################################################################################################################################################################################################################################################################################################################

# PLOTTING CIRCULAR PLOT FOR UMAP

library(plot1cell)
library(Seurat)

#library(reticulate)
#use_condaenv("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate")
#use_python("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/bin/python")

#install.packages("renv")
#renv::init()
#renv::install("reticulate")
#renv::use_python()

#py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)
#reticulate::py_install(py_pkgs)
#renv::snapshot()

#library(reticulate)
#sc <- import("scanpy")
#adata_subset2 <- sc$AnnData$read_h5ad('./adata_subset2')

#exprs <- t(adata$X)
#colnames(exprs) <- adata$obs_names$to_list()
#rownames(exprs) <- adata$var_names$to_list()
#seurat <- CreateSeuratObject(exprs)
# Add observation metadata
#seurat <- AddMetaData(seurat, adata$obs)
# Add fetaure metadata
#seurat[["RNA"]][["n_cells"]] <- adata$var["n_cells"]
# Add embedding
#embedding <- adata$obsm["X_umap"]
#rownames(embedding) <- adata$obs_names$to_list()
#colnames(embedding) <- c("umap_1", "umap_2")
#seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")

saveRDS(adata_subset2_seurat, "adata_subset2_seurat.rds")
adata_subset2_seurat <- readRDS("adata_subset2_seurat.rds")
ad <- readRDS("adata_subset2.rds")

# Add embedding to seurat object
embedding <- ad@int_colData@listData[["reducedDims"]]@listData[["X_umap"]]
rownames(embedding) <- colnames(ad)
colnames(embedding) <- c("umap_1", "umap_2")
adata_subset2_seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")

Idents(adata_subset2_seurat) <- "pheno_cluster_new"

###Prepare data for ploting - merged atlas
circ_data <- prepare_circlize_data(adata_subset2_seurat, scale = 0.75)
set.seed(1234)


celltype.colours <- setNames(c("#000000","#b5bbe3","#ff46a1","#ffe6f2", "#A64D79",
                               "#6fa8dc",'#4d7191',"#8faac2", "#d9d2e9",
                               "#336600","#8dd593","#0fcfc0","#9cded6",
                               "#b86cb9","#00B0F0","#4900EF","#4a6fe3",
                               "#FFD966","#8595e1","#A381EF", "#0070C0","#e07b91","#fce5cd", "#b86cb9","#d33f6a",
                               "#FF7DA8",'#FA8000',"#b9877d", "#f0b98d",
                               "#f4cccc","#A53F02","#d6bcc0","red"
),
c('B cell', 'CD11c+ cell', 'EM CD4 T cell', 'Proliferative CD4 T cell', 'CD4 Treg cell',
  'CD8 T cell', 'EM CD8 T cell', 'Proliferative CD8 T cell' , 'CD3+ cell',  
  'SARSCoV2+ ArginaseHighVISTAHigh Activated Neutrophil', 
  'ArginaseHighVISTAHigh Activated Neutrophil', 'ArginaseLowVISTALow Activated Neutrophil', 
  'ArginaseLowVISTALow Neutrophil', 
  'Classical Monocyte', 'SARSCoV2+ Interstitial Macrophage', 'Proliferating Interstitial Macrophage', 'Interstitial Macrophage', 
  'Apoptotic SARSCoV2+  Alveolar Macrophage', 'Apoptotic Alveolar Macrophage', 'Alveolar Macrophage', 
  'Activated Endothelial cell', 'Endothelial cell', 'Proliferative Endothelial cell', 'Fibroblast', 'Proliferative Fibroblast', 'Apoptotic Fibroblast',  
  'Smooth Muscle cell', 'SARSCoV2+ Epithelial cell', 'Epithelial cell' , 'Proliferative Epithelial cell', 'SARSCoV2+ AT2 cell', 'AT2 cell', 'RBC'))

cluster_colors <- celltype.colours
disease_colors<- c("#F1CDB1", "#CBDFB8", "#8FA9DF")
names(disease_colors) <- c("COVID-19", "Non_Pneumonia", "Pneumonia")
hiv_colors<-c("lightblue", "red")
names(hiv_colors) <- c("HIV_Neg", "HIV_Pos")

###plot and save figures
pdf('circlize_plot_lung_cleaned.pdf', width = 6, height = 6)
plot_circlize(circ_data, do.label = FALSE, pt.size = 0.1, contour.levels = c(0, 0), col.use = cluster_colors, bg.color = 'white', kde2d.n = 200, repel = TRUE, label.cex = 0,
              contour.nlevels = 0)
add_track(circ_data, group = "Group", colors = disease_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "HIV", colors = hiv_colors, track_num = 3)
dev.off()

saveRDS(circ_data, "data_circular_umap.rds")



