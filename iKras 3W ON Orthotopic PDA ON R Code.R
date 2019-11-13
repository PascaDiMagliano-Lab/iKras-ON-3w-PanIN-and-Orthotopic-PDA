# This script reproduces the single cell RNAseq analysis of mouse pancreatic samples and
# Figures of the paper Regulatory T Cell depletion causes compensatory immune suppression and accelerated pancreatic carcinogenesis. 
# The raw data was processed in line with the 
# Seurat workflow outlined on the Satija Lab website 
# (https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html) as well as in the following
# reference: 
#
# Stuart et al., Comprehensive Integration of Single-Cell Data.
# Cell, 2019;7:177.
#

#Load required packages

library(Seurat)
library(dplyr)
library(harmony)
library(tibble)
library(ggplot2)
library(stringr)

#Load Functions

MouseAutomatedClusterMarkerTable <- function(Seurat_Object){
  library(dplyr)
  library(tibble)
  library(Seurat)
  ClusterList <- list()
  Idents(object = Seurat_Object) <- "seurat_clusters"
  current.cluster.ids <- sort(as.numeric(levels(Seurat_Object@active.ident)))
  new.cluster.ids <- c()
  
  
  for(i in current.cluster.ids){
    List_Position <- i + 1
    ClusterList[[List_Position]] <- FindMarkers(object = Seurat_Object, ident.1 = i, min.pct = 0.25, only.pos = TRUE)
    Positive_Genes <- rownames(ClusterList[[List_Position]])
    Num_Positive_Genes <- length(Positive_Genes)
    
    RPS_Num <- length(grep(pattern = "^Rps", x = Positive_Genes))
    RPL_Num <- length(grep(pattern = "^Rpl", x = Positive_Genes))
    RP_Percent <- sum(RPS_Num, RPL_Num)/length(Positive_Genes)*100
    RP_Label <- paste("RP%:", RP_Percent, sep = " ")
    
    Mito_Num <- length(grep(pattern = "^mt-", x = Positive_Genes))
    Mito_Percent <- Mito_Num/length(Positive_Genes)*100
    Mito_Label <- paste("Mito%:", RP_Percent, sep = " ")
    
    ClusterCells <- WhichCells(object = Seurat_Object, idents = i)
    Cell_Barcodes <- unlist(Seurat_Object@assays$RNA@counts@Dimnames[2])
    Cell_Number <- c()
    for(k in 1:length(ClusterCells)){
      Cell_Position <- grep(pattern = ClusterCells[k],x = Cell_Barcodes, value = FALSE)
      Cell_Number <- c(Cell_Number,Cell_Position)
    }
    
    
    S_Score <- Seurat_Object@meta.data$S.Score
    G2M_Score <- Seurat_Object@meta.data$G2M.Score
    Cluster_S_Score <- S_Score[Cell_Number]
    Cluster_G2M_Score <- G2M_Score[Cell_Number]
    Avg_Cluster_S_Score <- mean(Cluster_S_Score)
    Avg_Cluster_G2M_Score <- mean(Cluster_G2M_Score)
    Cluster_S_Score_Range <- range(Cluster_S_Score)
    Cluster_G2M_Score_Range <- range(Cluster_G2M_Score)
    
    nFeature <- Seurat_Object@meta.data$nFeature_RNA
    nCount <- Seurat_Object@meta.data$nCount_RNA
    Mito <- Seurat_Object@meta.data$percent.mt
    
    Cluster_nFeature <- nFeature[Cell_Number]
    Cluster_nCount <- nCount[Cell_Number]
    Cluster_Mito <- Mito[Cell_Number]
    
    Avg_Cluster_nFeature <- as.integer(mean(Cluster_nFeature))
    Avg_Cluster_nCount <- as.integer(mean(Cluster_nCount))
    Max_Cluster_Mito <- max(Cluster_Mito)
    
    Cell_Types <- c("Epi","T Cell","Myeloid","B Cell","Fibroblast","RBC","NK", "Endo","Acinar")
    
    Epi_Markers <- c("Krt7","Krt8","Krt18","Krt19","Epcam","Cdh1")
    T_Cell_Markers <- c("Cd3e","Cd3g","Cd3d","Cd4","Il7r","Cd8a","Lef1")
    Myeloid_Markers <- c("Cd14","Itgam","Mnda","Mpeg1","Itgax")
    B_Cell_Markers <- c("Cd79a","Ms4a1","Cd19")
    Fibroblast_Markers <- c("Cdh11","Pdgfra","Pdgfrb","Acta2")
    RBC_Markers <- c("Hba-a1","Hbb-bs","Hba-a2","Hbb-bt")
    NK_Markers <- c("Ncr3","Fcgr3","Ncam1","Klrf1","Klrc1","Cd38","Klrc1")
    Endo_Markers <- c("Cdh5","Pecam1")
    Acinar_Markers <- c("Try4","Spink1","Amy2a2")
    All_Markers <- list(Epi_Markers,T_Cell_Markers,Myeloid_Markers,B_Cell_Markers,Fibroblast_Markers,RBC_Markers,NK_Markers,Endo_Markers,Acinar_Markers)
    
    Epi_Score <- 0
    T_Cell_Score <- 0
    Myeloid_Score <- 0
    B_Cell_Score <- 0
    Fibroblast_Score <- 0
    RBC_Score <- 0
    NK_Score <- 0
    Endo_Score <- 0
    Acinar_Score <- 0 
    All_Scores <- list(Epi_Score,T_Cell_Score,Myeloid_Score,B_Cell_Score,Fibroblast_Score,RBC_Score,NK_Score,Endo_Score,Acinar_Score)
    Weighted_Scores <- c()
    Score_Weights <- c(1.85,1.85,2.22,3.7,2.78,3.7,1.85,5.56,3.7) 
    
    for(h in 1:length(All_Markers)){
      Markers_to_Test<- All_Markers[[h]]
      Marker_Row <- h
      for(j in 1:length(Markers_to_Test)){
        Gene_Found <- 0
        Gene_Found <- length(grep(pattern = Markers_to_Test[j], x = Positive_Genes))
        if(Gene_Found > 0 ){
          All_Scores[[Marker_Row]] <- All_Scores[[Marker_Row]]+1
        }
      }
      Weighted_Scores[Marker_Row] <- All_Scores[[Marker_Row]]*Score_Weights[Marker_Row]
    }
    
    ClusterID <- which(Weighted_Scores >= 5.5)
    if(length(ClusterID) > 0){
      if(length(ClusterID) > 1){
        ID <- "Multiple"
      }else{
        ID <- Cell_Types[ClusterID]
      }
    }else{
      ID <- i
    } 
    if(RP_Percent > 30){
      ID <- paste("RP_",ID,sep = "")
    }
    if(Avg_Cluster_S_Score > 0.01 | Avg_Cluster_G2M_Score > 0.01){
      CellCycleID <- "Cycling"
      ID <- paste("Cycling_",ID,sep = "")
    }else{
      CellCycleID <- "N/A"
    }
    if(Avg_Cluster_nCount < 700){
      ID <- paste("G_",ID,sep = "")
    }
    new.cluster.ids <- c(new.cluster.ids,ID)
    
    Label_Row <- length(Positive_Genes) + 1
    Label_Row2 <- length(Positive_Genes) + 2
    Label_Row3 <- length(Positive_Genes) + 3
    Label_Row4 <- length(Positive_Genes) + 4
    Label_Row5 <- length(Positive_Genes) + 5
    
    Label1 <- c("Summary:",paste("Cluster",i, sep = " "), paste("ID:",ID, sep = " "),paste("Mito%:",Mito_Percent, sep = " "),paste("RP%:",RP_Percent, sep = " "))
    Label2 <- c("Immune Summary",paste("T Cell Score:",All_Scores[[2]], sep = " "),paste("Myeloid Score:",All_Scores[[3]], sep = " "),paste("B Cell Score:",All_Scores[[4]], sep = " "),
                paste("NK Score:",All_Scores[[7]], sep = " "))
    Label3 <- c(paste("Epi Score:",All_Scores[[1]], sep = " "),paste("Fib Score:",All_Scores[[5]], sep = " "),paste("Acinar Score:",All_Scores[[9]], sep = " "),
                paste("Endo Score:",All_Scores[[8]], sep = " "),paste("RBC Score:",All_Scores[[6]], sep = " "))
    Label4 <- c("Avg S Score:", Avg_Cluster_S_Score, "Avg G2M Score:", Avg_Cluster_G2M_Score, CellCycleID)
    Label5 <- c("Filter Info",paste("Avg. nGene:", Avg_Cluster_nFeature, sep = " "),paste("Avg. nCounts:", Avg_Cluster_nCount, sep = " ")
                , paste("Highest Mito:",Max_Cluster_Mito, sep = " "), paste("# Cells:",length(Cell_Number), sep = " "))
    
    ClusterList[[List_Position]][Label_Row,] <- Label1
    ClusterList[[List_Position]][Label_Row2,] <- Label2
    ClusterList[[List_Position]][Label_Row3,] <- Label3
    ClusterList[[List_Position]][Label_Row4,] <- Label4
    ClusterList[[List_Position]][Label_Row5,] <- Label5
    ClusterList[[List_Position]] <- rownames_to_column(.data = ClusterList[[List_Position]],var = "Gene")
    ClusterList[[List_Position]][Label_Row,"Gene"] <- "Summary1"
    ClusterList[[List_Position]][Label_Row2,"Gene"] <- "Summary2"
    ClusterList[[List_Position]][Label_Row3,"Gene"] <- "Summary3"
    ClusterList[[List_Position]][Label_Row4,"Gene"] <- "Summary4"
    ClusterList[[List_Position]][Label_Row5,"Gene"] <- "Summary5"
    
    
  }
  
  ClusterDataFrame <- bind_rows(ClusterList, .id = "column_label")
  ClusterDataFrame <- ClusterDataFrame[,-1]
  ClusterPackage <- list(ClusterDataFrame, new.cluster.ids)
  return(ClusterPackage)
}


#Load in raw or final objects

#Preprocess raw objects

#ON1 
#Create Seurat Object and Pre-filter Genes and Cells
ON1 <- CreateSeuratObject(counts = On1SampleRaw.data, project = "ON1", min.cells = 3, min.features = 100)

#Mouse Metadata
ON1[["Mouse_ID"]] <- "ON1"
ON1[["Genotype"]] <- "iKras_p53"
ON1[["Run_ID"]] <- "2401"
ON1[["Sequencer"]] <- "HiSEQ-4000"
ON1[["Tumor_Type"]] <- "Orthotopic"
ON1[["Time_Point"]] <- "ON"
ON1[["Background"]] <- "FVBN"


#Calculate MT %
ON1[["percent.mt"]] <- PercentageFeatureSet(object = ON1, pattern = "^mt-")

#Normalize Data
ON1 <- NormalizeData(object = ON1, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
ON1 <- FindVariableFeatures(object = ON1, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
ON1 <- CellCycleScoring(ON1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ON1$CC.Difference <- ON1$S.Score - ON1$G2M.Score

#Scale Data
all.genes <- rownames(x = ON1)
ON1 <- ScaleData(ON1, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
ON1 <- RunPCA(object = ON1, features = VariableFeatures(object = ON1))
stdev <- ON1@reductions$pca@stdev
var <- stdev^2
sum(var[1:31])/ sum(var)

#Find Neighbors + Find CLusters
ON1 <- FindNeighbors(object = ON1, dims = 1:31)
ON1 <- FindClusters(object = ON1, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot
ON1 <- RunUMAP(object = ON1, dims = 1:31)
DimPlot(object = ON1, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Run ACMT Function for Annotating
CMT_ON1 <- MouseAutomatedClusterMarkerTable(Seurat_Object = ON1)
write.csv(CMT_ON1[[1]], file = "~/Desktop/CMT_ON1.csv") 
View(CMT_ON1[[1]])

#Master FeaturePlot
FeaturePlot(object = ON1, features = c("Ptprc", "Krt18", "Krt19", "Epcam", "Cd3e", "Cd4", "Cd8a", "Il7r", "Itgam", "Itgax", "Cd14", "Cd33", "Ncam1", "Fcgr3", "Nkg7", "Klrk1", "Cd19", "Cd22", "Ms4a1", "Tnfrsf17", "Pdgfra", "Pdgfrb", "Acta2", "Cdh11"), pt.size = 0.5)

#Visualize Labelled Clusters
new.cluster.ids <- CMT_ON1[[2]]
names(x = new.cluster.ids) <- sort(as.numeric(levels(x = ON1)), decreasing = F)
ON1 <- RenameIdents(object = ON1, new.cluster.ids)
ON1[["Auto_Labelled_Clusters"]] <- ON1@active.ident
DimPlot(object = ON1, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON1@active.ident)

#Manually Label Clusters
ON1 <- RenameIdents(ON1, "1" = "G", "6" = "Myeloid", `8` = "Matrix", "10" = "Myeloid", "Cycling_13" = "Myeloid")
ON1[["Manual_Labelled_Clusters"]] <- ON1@active.ident
DimPlot(object = ON1, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Identity Switches
Idents(object = ON1) <- "seurat_clusters"
Idents(object = ON1) <- "Auto_Labelled_Clusters"
Idents(object = ON1) <- "Manual_Labelled_Clusters"

# save Seurat object
save(ON1, file = 'ON1.RData')

#Determine n Cutoffs
VlnPlot(ON1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ON1_Fil <- subset(x = ON1, subset = nCount_RNA < 55000 )
ON1_Fil <- subset(x = ON1, subset = nCount_RNA < 55000 & percent.mt < 33 )
ON1_Fil <- subset(x = ON1, subset = nCount_RNA > 700 & nCount_RNA < 55000 & percent.mt < 33)
VlnPlot(ON1_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON1_Fil@active.ident)

#Subset Seurat Object
ON1_Fil <- subset(x = ON1, subset = nCount_RNA > 700 & nCount_RNA < 55000 & percent.mt < 33, 
                  idents =  c("Epi","Myeloid","Cycling_Epi","T Cell","Fibroblast","RBC","Acinar","Endo","Cycling_RBC"))

#Normalize Data
ON1_Fil <- NormalizeData(object = ON1_Fil, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
ON1_Fil <- FindVariableFeatures(object = ON1_Fil, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
ON1_Fil <- CellCycleScoring(ON1_Fil, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ON1_Fil$CC.Difference <- ON1_Fil$S.Score - ON1_Fil$G2M.Score

#Scale Data
all.genes <- rownames(x = ON1_Fil)
ON1_Fil <- ScaleData(ON1_Fil, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
ON1_Fil <- RunPCA(object = ON1_Fil, features = VariableFeatures(object = ON1_Fil))
stdev <- ON1_Fil@reductions$pca@stdev
var <- stdev^2
sum(var[1:31])/ sum(var)

#Find Neighbors + Find CLusters
ON1_Fil <- FindNeighbors(object = ON1_Fil, dims = 1:31)
ON1_Fil <- FindClusters(object = ON1_Fil, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot
ON1_Fil <- RunUMAP(object = ON1_Fil, dims = 1:31)
DimPlot(object = ON1_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON1_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Run ACMT Function for Annotating
CMT_ON1_Fil <- MouseAutomatedClusterMarkerTable(Seurat_Object = ON1_Fil)
write.csv(CMT_ON1_Fil[[1]], file = "~/Desktop/CMT_ON1_Fil.csv") 
View(CMT_ON1_Fil[[1]])

#Master FeaturePlot
FeaturePlot(object = ON1_Fil, features = c("Ptprc", "Krt18", "Krt19", "Epcam", "Cd3e", "Cd4", "Cd8a", "Il7r", "Itgam", "Itgax", "Cd14", "Cd33", "Ncam1", "Fcgr3", "Nkg7", "Klrk1", "Cd19", "Cd22", "Ms4a1", "Tnfrsf17", "Pdgfra", "Pdgfrb", "Acta2", "Cdh11"), pt.size = 0.5)

#Visualize Labelled Clusters
Idents(object = ON1_Fil) <- "seurat_clusters"
new.cluster.ids <- CMT_ON1_Fil[[2]]
names(x = new.cluster.ids) <- sort(as.numeric(levels(x = ON1_Fil)), decreasing = F)
ON1_Fil <- RenameIdents(object = ON1_Fil, new.cluster.ids)
ON1_Fil[["Auto_Labelled_Clusters"]] <- ON1_Fil@active.ident
DimPlot(object = ON1_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON1_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON1_Fil@active.ident)

#Final Labelling
Idents(object = ON1_Fil) <- "seurat_clusters"
ON1_Fil <- RenameIdents(ON1_Fil, 
                        "0" = "Epithelial",
                        "1" = "Epithelial",
                        "2" = "Epithelial",
                        "3" = "Fibroblasts",
                        "4" = "RBC",
                        "5" = "Myeloid",
                        "6" = "RBC",
                        "7" = "T Cell",
                        "8" = "RBC",
                        "9" = "Myeloid",
                        "10" = "Myeloid",
                        "11" = "RBC",
                        "12" = "RBC",
                        "13" = "Myeloid",
                        "14" = "NK",
                        "15" = "Fibroblasts",
                        "16" = "Endothelial",
                        "17" = "Acinar")
ON1_Fil[["Final_Labelled_Clusters"]] <- ON1_Fil@active.ident
DimPlot(object = ON1_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON1_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Identity Switches
Idents(object = ON1_Fil) <- "seurat_clusters"
Idents(object = ON1_Fil) <- "Auto_Labelled_Clusters"
Idents(object = ON1_Fil) <- "Manual_Labelled_Clusters"
Idents(object = ON1_Fil) <- "Final_Labelled_Clusters"

# save Seurat object
save(ON1_Fil, file = 'ON1_Fil.RData')


#ON2
#Create Seurat Object and Pre-filter Genes and Cells
ON2 <- CreateSeuratObject(counts = On2SampleRaw.data, project = "ON2", min.cells = 3, min.features = 100)

#Mouse Metadata
ON2[["Mouse_ID"]] <- "ON2"
ON2[["Genotype"]] <- "iKras_p53"
ON2[["Run_ID"]] <- "2560"
ON2[["Sequencer"]] <- "HiSEQ-4000"
ON2[["Tumor_Type"]] <- "Orthotopic"
ON2[["Time_Point"]] <- "ON"
ON2[["Background"]] <- "FVBN"


#Calculate MT %
ON2[["percent.mt"]] <- PercentageFeatureSet(object = ON2, pattern = "^mt-")

#Normalize Data
ON2 <- NormalizeData(object = ON2, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
ON2 <- FindVariableFeatures(object = ON2, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
ON2 <- CellCycleScoring(ON2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ON2$CC.Difference <- ON2$S.Score - ON2$G2M.Score

#Scale Data
all.genes <- rownames(x = ON2)
ON2 <- ScaleData(ON2, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
ON2 <- RunPCA(object = ON2, features = VariableFeatures(object = ON2))
stdev <- ON2@reductions$pca@stdev
var <- stdev^2
sum(var[1:36])/ sum(var)

#Find Neighbors + Find CLusters
ON2 <- FindNeighbors(object = ON2, dims = 1:36)
ON2 <- FindClusters(object = ON2, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot
ON2 <- RunUMAP(object = ON2, dims = 1:36)
DimPlot(object = ON2, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Run ACMT Function for Annotating
CMT_ON2 <- MouseAutomatedClusterMarkerTable(Seurat_Object = ON2)
write.csv(CMT_ON2[[1]], file = "~/Desktop/CMT_ON2.csv") 
View(CMT_ON2[[1]])

#Master FeaturePlot
FeaturePlot(object = ON2, features = c("Ptprc", "Krt18", "Krt19", "Epcam", "Cd3e", "Cd4", "Cd8a", "Il7r", "Itgam", "Itgax", "Cd14", "Cd33", "Ncam1", "Fcgr3", "Nkg7", "Klrk1", "Cd19", "Cd22", "Ms4a1", "Tnfrsf17", "Pdgfra", "Pdgfrb", "Acta2", "Cdh11"), pt.size = 0.5)

#Visualize Labelled Clusters
new.cluster.ids <- CMT_ON2[[2]]
names(x = new.cluster.ids) <- sort(as.numeric(levels(x = ON2)), decreasing = F)
ON2 <- RenameIdents(object = ON2, new.cluster.ids)
ON2[["Auto_Labelled_Clusters"]] <- ON2@active.ident
DimPlot(object = ON2, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON2@active.ident)

#Manually Label Clusters
ON2 <- RenameIdents(ON2, "3" = "G", "6" = "Myeloid", "7" = "T Cell", "9" = "Myeloid", "11" = "Matrix", "12" = "Myeloid")
ON2[["Manual_Labelled_Clusters"]] <- ON2@active.ident
DimPlot(object = ON2, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Identity Switches
Idents(object = ON2) <- "seurat_clusters"
Idents(object = ON2) <- "Auto_Labelled_Clusters"
Idents(object = ON2) <- "Manual_Labelled_Clusters"

# save Seurat object
save(ON2, file = 'ON2.RData')

#Determine n Cutoffs
VlnPlot(ON2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ON2_Fil <- subset(x = ON2, subset = nCount_RNA < 50000 )
ON2_Fil <- subset(x = ON2, subset = nCount_RNA < 50000 & percent.mt < 60 )
ON2_Fil <- subset(x = ON2, subset = nCount_RNA > 700 & nCount_RNA < 50000 & percent.mt < 50)
VlnPlot(ON2_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON2_Fil@active.ident)

#Subset Seurat Object
ON2_Fil <- subset(x = ON2, subset = nCount_RNA > 700 & nCount_RNA < 50000 & percent.mt < 50, 
                  idents =  c("Epi","Myeloid","Cycling_T Cell","Fibroblast","Cycling_RBC","Cycling_Epi","RBC"))

#Normalize Data
ON2_Fil <- NormalizeData(object = ON2_Fil, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
ON2_Fil <- FindVariableFeatures(object = ON2_Fil, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
ON2_Fil <- CellCycleScoring(ON2_Fil, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ON2_Fil$CC.Difference <- ON2_Fil$S.Score - ON2_Fil$G2M.Score

#Scale Data
all.genes <- rownames(x = ON2_Fil)
ON2_Fil <- ScaleData(ON2_Fil, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
ON2_Fil <- RunPCA(object = ON2_Fil, features = VariableFeatures(object = ON2_Fil))
stdev <- ON2_Fil@reductions$pca@stdev
var <- stdev^2
sum(var[1:34])/ sum(var)

#Find Neighbors + Find CLusters
ON2_Fil <- FindNeighbors(object = ON2_Fil, dims = 1:34)
ON2_Fil <- FindClusters(object = ON2_Fil, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot
ON2_Fil <- RunUMAP(object = ON2_Fil, dims = 1:34)
DimPlot(object = ON2_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON2_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Run ACMT Function for Annotating
CMT_ON2_Fil <- MouseAutomatedClusterMarkerTable(Seurat_Object = ON2_Fil)
write.csv(CMT_ON2_Fil[[1]], file = "~/Desktop/CMT_ON2_Fil.csv") 
View(CMT_ON2_Fil[[1]])

#Master FeaturePlot
FeaturePlot(object = ON2_Fil, features = c("Ptprc", "Krt18", "Krt19", "Epcam", "Cd3e", "Cd4", "Cd8a", "Il7r", "Itgam", "Itgax", "Cd14", "Cd33", "Ncam1", "Fcgr3", "Nkg7", "Klrk1", "Cd19", "Cd22", "Ms4a1", "Tnfrsf17", "Pdgfra", "Pdgfrb", "Acta2", "Cdh11"), pt.size = 0.5)

#Visualize Labelled Clusters
Idents(object = ON2_Fil) <- "seurat_clusters"
new.cluster.ids <- CMT_ON2_Fil[[2]]
names(x = new.cluster.ids) <- sort(as.numeric(levels(x = ON2_Fil)), decreasing = F)
ON2_Fil <- RenameIdents(object = ON2_Fil, new.cluster.ids)
ON2_Fil[["Auto_Labelled_Clusters"]] <- ON2_Fil@active.ident
DimPlot(object = ON2_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON2_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON2_Fil@active.ident)

#Final Labelling
Idents(object = ON2_Fil) <- "seurat_clusters"
ON2_Fil <- RenameIdents(ON2_Fil, 
                        "0" = "RBC",
                        "1" = "Epithelial",
                        "2" = "RBC",
                        "3" = "RBC",
                        "4" = "Epithelial",
                        "5" = "Myeloid",
                        "6" = "Epithelial",
                        "7" = "Epithelial",
                        "8" = "T Cell",
                        "9" = "Fibroblasts",
                        "10" = "Myeloid",
                        "11" = "Myeloid",
                        "12" = "Myeloid",
                        "13" = "NK",
                        "14" = "Acinar")
ON2_Fil[["Final_Labelled_Clusters"]] <- ON2_Fil@active.ident
DimPlot(object = ON2_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON2_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Identity Switches
Idents(object = ON2_Fil) <- "seurat_clusters"
Idents(object = ON2_Fil) <- "Auto_Labelled_Clusters"
Idents(object = ON2_Fil) <- "Manual_Labelled_Clusters"
Idents(object = ON2_Fil) <- "Final_Labelled_Clusters"

# save Seurat object
save(ON2_Fil, file = 'ON2_Fil.RData')


#Analysis of final objects
ON1_Fil <- RenameCells(object = ON1_Fil, add.cell.id = "ON1")
ON2_Fil <- RenameCells(object = ON2_Fil, add.cell.id = "ON2")

#Merge Objects
Mouse_On_Merge <- merge(x = ON1_Fil, y = ON2_Fil)

#Calculate MT %
Mouse_On_Merge[["percent.mt"]] <- PercentageFeatureSet(object = Mouse_On_Merge, pattern = "^mt-")

#Normalize Data
Mouse_On_Merge <- NormalizeData(object = Mouse_On_Merge, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
Mouse_On_Merge <- FindVariableFeatures(object = Mouse_On_Merge, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
Mouse_On_Merge <- CellCycleScoring(Mouse_On_Merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Mouse_On_Merge$CC.Difference <- Mouse_On_Merge$S.Score - Mouse_On_Merge$G2M.Score

#Scale Data
all.genes <- rownames(x = Mouse_On_Merge)
Mouse_On_Merge <- ScaleData(Mouse_On_Merge, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
Mouse_On_Merge <- RunPCA(object = Mouse_On_Merge, features = VariableFeatures(object = Mouse_On_Merge))
stdev <- Mouse_On_Merge@reductions$pca@stdev
var <- stdev^2
sum(var[1:30])/ sum(var)

#Batch Correction
Mouse_On_Merge <- RunHarmony(Mouse_On_Merge, group.by.vars = c("Run_ID","Sequencer"), dims.use = 1:30, verbose = F)

#Find Neighbors + Find CLusters
Mouse_On_Merge <- FindNeighbors(object = Mouse_On_Merge, dims = 1:30)
Mouse_On_Merge <- FindNeighbors(object = Mouse_On_Merge, dims = 1:30, reduction = "harmony")
Mouse_On_Merge <- FindClusters(object = Mouse_On_Merge, resolution = 3)

#Run UMAP and get unlabelled cluster UMAP and violin plot
Mouse_On_Merge <- RunUMAP(object = Mouse_On_Merge, dims = 1:30)
Mouse_On_Merge <- RunUMAP(object = Mouse_On_Merge, dims = 1:30, reduction = "harmony")
DimPlot(object = Mouse_On_Merge, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(Mouse_On_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Identity Switches
Idents(object = Mouse_On_Merge) <- "seurat_clusters"
Idents(object = Mouse_On_Merge) <- "Labelled_Clusters"
Idents(object = Mouse_On_Merge) <- "Patient_ID"
Idents(object = Mouse_On_Merge) <- "State"
Idents(object = Mouse_On_Merge) <- "Run_ID"
Idents(object = Mouse_On_Merge) <- "Time Point"
Idents(object = Mouse_On_Merge) <- "Collection_Method"


#Run ACMT Function for Annotating
CMT_Mouse_On_Merge <- MouseAutomatedClusterMarkerTable(Seurat_Object = Mouse_On_Merge)
write.csv(CMT_Mouse_On_Merge[[1]], file = "~/Desktop/CMT_Mouse_On_Merge.csv") 
View(CMT_Mouse_On_Merge[[1]])

FeaturePlot(object = Mouse_On_Merge, features = c("Ptprc", "Krt18", "Krt19", "Epcam", "Cd3e", "Cd4", "Cd8a", "Il7r", "Itgam", "Itgax", "Cd14", "Cd33", "Ncam1", "Fcgr3", "Nkg7", 
                                                  "Klrk1", "Cd19", "Cd22", "Ms4a1", "Pdgfra", "Pdgfrb", "Acta2", "Cdh11"), pt.size = 0.5)

#Visualize Labelled Clusters
Idents(object = Mouse_On_Merge) <- "seurat_clusters"
new.cluster.ids <- CMT_Mouse_On_Merge[[2]]
names(x = new.cluster.ids) <- sort(as.numeric(levels(x = Mouse_On_Merge)), decreasing = F)
Mouse_On_Merge <- RenameIdents(object = Mouse_On_Merge, new.cluster.ids)
Mouse_On_Merge[["Auto_Labelled_Clusters"]] <- Mouse_On_Merge@active.ident
DimPlot(object = Mouse_On_Merge, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(Mouse_On_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(Mouse_On_Merge@active.ident)

#Final Labelling
Idents(object = Mouse_On_Merge) <- "seurat_clusters"
Mouse_On_Merge <- RenameIdents(Mouse_On_Merge, 
                               "0" = "RBC",
                               "1" = "Epithelial",
                               "2" = "Epithelial",
                               "3" = "RBC",
                               "4" = "Myeloid 1",
                               "5" = "RBC",
                               "6" = "Epithelial",
                               "7" = "Epithelial",
                               "8" = "Epithelial",
                               "9" = "RBC",
                               "10" = "Epithelial",
                               "11" = "Epithelial",
                               "12" = "Fibroblast",
                               "13" = "Epithelial",
                               "14" = "RBC",
                               "15" = "RBC",
                               "16" = "Treg",
                               "17" = "RBC",
                               "18" = "Myeloid 1",
                               "19" = "Myeloid 1",
                               "20" = "RBC",
                               "21" = "CD8+ T Cell",
                               "22" = "CD4+ T Cell",
                               "23" = "Fibroblast",
                               "24" = "Myeloid 2",
                               "25" = "RBC",
                               "26" = "Myeloid 1",
                               "27" = "Myeloid 3",
                               "28" = "Fibroblast",
                               "29" = "Myeloid 1",
                               "30" = "Epithelial",
                               "31" = "Myeloid 1",
                               "32" = "B Cell",
                               "33" = "Fibroblast",
                               "34" = "Acinar",
                               "35" = "Endothelial")
Mouse_On_Merge[["Final_Labelled_Clusters"]] <- Mouse_On_Merge@active.ident
DimPlot(object = Mouse_On_Merge, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(Mouse_On_Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter Out RBC
Mouse_On_Merge_Fil <- subset(x = Mouse_On_Merge, idents = c("Epithelial","Myeloid 1","Myeloid 2","Myeloid 3","Fibroblast","Treg","B Cell",
                                                            "CD4+ T Cell","CD8+ T Cell","Endothelial","Acinar"))
DimPlot(object = Mouse_On_Merge_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)

#Calculate MT %
Mouse_On_Merge_Fil[["percent.mt"]] <- PercentageFeatureSet(object = Mouse_On_Merge_Fil, pattern = "^mt-")

#Normalize Data
Mouse_On_Merge_Fil <- NormalizeData(object = Mouse_On_Merge_Fil, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
Mouse_On_Merge_Fil <- FindVariableFeatures(object = Mouse_On_Merge_Fil, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
Mouse_On_Merge_Fil <- CellCycleScoring(Mouse_On_Merge_Fil, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Mouse_On_Merge_Fil$CC.Difference <- Mouse_On_Merge_Fil$S.Score - Mouse_On_Merge_Fil$G2M.Score

#Scale Data
all.genes <- rownames(x = Mouse_On_Merge_Fil)
Mouse_On_Merge_Fil <- ScaleData(Mouse_On_Merge_Fil, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
Mouse_On_Merge_Fil <- RunPCA(object = Mouse_On_Merge_Fil, features = VariableFeatures(object = Mouse_On_Merge_Fil))
stdev <- Mouse_On_Merge_Fil@reductions$pca@stdev
var <- stdev^2
sum(var[1:32])/ sum(var)

#Batch Correction
Mouse_On_Merge_Fil <- RunHarmony(Mouse_On_Merge_Fil, group.by.vars = c("Run_ID","Sequencer"), dims.use = 1:32, verbose = F)

#Find Neighbors + Find CLusters
Mouse_On_Merge_Fil <- FindNeighbors(object = Mouse_On_Merge_Fil, dims = 1:32)
Mouse_On_Merge_Fil <- FindNeighbors(object = Mouse_On_Merge_Fil, dims = 1:32, reduction = "harmony")
Mouse_On_Merge_Fil <- FindClusters(object = Mouse_On_Merge_Fil, resolution = 5)

#Run UMAP and get unlabelled cluster UMAP and violin plot
Mouse_On_Merge_Fil <- RunUMAP(object = Mouse_On_Merge_Fil, dims = 1:32)
Mouse_On_Merge_Fil <- RunUMAP(object = Mouse_On_Merge_Fil, dims = 1:32, reduction = "harmony")
DimPlot(object = Mouse_On_Merge_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(Mouse_On_Merge_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Run ACMT Function for Annotating
CMT_Mouse_On_Merge_Fil <- MouseAutomatedClusterMarkerTable(Seurat_Object = Mouse_On_Merge_Fil)
write.csv(CMT_Mouse_On_Merge_Fil[[1]], file = "~/Desktop/CMT_Mouse_On_Merge_Fil.csv") 
View(CMT_Mouse_On_Merge[[1]])

#Final Labelling
Idents(object = Mouse_On_Merge_Fil) <- "seurat_clusters"
Mouse_On_Merge_Fil <- RenameIdents(Mouse_On_Merge_Fil, 
                                   "0" = "Epithelial",
                                   "1" = "Myeloid 1",
                                   "2" = "Epithelial",
                                   "3" = "Fibroblasts",
                                   "4" = "Epithelial",
                                   "5" = "Epithelial",
                                   "6" = "Epithelial",
                                   "7" = "Epithelial",
                                   "8" = "Tregs",
                                   "9" = "Epithelial",
                                   "10" = "Myeloid 1",
                                   "11" = "Fibroblasts",
                                   "12" = "Myeloid 1",
                                   "13" = "Epithelial",
                                   "14" = "CD4+ T Cells",
                                   "15" = "Myeloid 1",
                                   "16" = "Myeloid 1",
                                   "17" = "Epithelial",
                                   "18" = "Myeloid 2",
                                   "19" = "Epithelial",
                                   "20" = "Epithelial",
                                   "21" = "Epithelial",
                                   "22" = "Epithelial",
                                   "23" = "Epithelial",
                                   "24" = "Fibroblasts",
                                   "25" = "Myeloid 3",
                                   "26" = "Epithelial",
                                   "27" = "Myeloid 1",
                                   "28" = "CD8+ T Cells",
                                   "29" = "Epithelial",
                                   "30" = "Epithelial",
                                   "31" = "NK Cells",
                                   "32" = "Epithelial",
                                   "33" = "B Cells",
                                   "34" = "Myeloid 1",
                                   "35" = "Fibroblasts",
                                   "36" = "Epithelial",
                                   "37" = "Acinar",
                                   "38" = "Endothelial")
Mouse_On_Merge_Fil[["Final_Labeled_Clusters"]] <- Mouse_On_Merge_Fil@active.ident
DimPlot(object = Mouse_On_Merge_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(Mouse_On_Merge_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Identity Switches
Idents(object = Mouse_On_Merge_Fil) <- "seurat_clusters"
Idents(object = Mouse_On_Merge_Fil) <- "Final_Labeled_Clusters"

# save Seurat object
save(Mouse_On_Merge_Fil, file = 'Mouse_On_Merge_Fil.RData')

#ON3
#Create Seurat Object and Pre-filter Genes and Cells
ON3 <- CreateSeuratObject(counts = On3SampleRaw.data, project = "ON3", min.cells = 3, min.features = 100)

#Mouse Metadata
ON3[["Mouse_ID"]] <- "ON3"
ON3[["Genotype"]] <- "iKras"
ON3[["Run_ID"]] <- "2800"
ON3[["Sequencer"]] <- "NovaSeq-6000"
ON3[["Tumor_Type"]] <- "Early Lesions"
ON3[["Time_Point"]] <- "ON 3 WK"
ON3[["Background"]] <- "FVBN"

#Calculate MT %
ON3[["percent.mt"]] <- PercentageFeatureSet(object = ON3, pattern = "^mt-")

#Normalize Data
ON3 <- NormalizeData(object = ON3, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
ON3 <- FindVariableFeatures(object = ON3, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
ON3 <- CellCycleScoring(ON3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ON3$CC.Difference <- ON3$S.Score - ON3$G2M.Score

#Scale Data
all.genes <- rownames(x = ON3)
ON3 <- ScaleData(ON3, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
ON3 <- RunPCA(object = ON3, features = VariableFeatures(object = ON3))
stdev <- ON3@reductions$pca@stdev
var <- stdev^2
sum(var[1:31])/ sum(var)

#Find Neighbors + Find CLusters
ON3 <- FindNeighbors(object = ON3, dims = 1:31, force.recalc = T)
ON3 <- FindClusters(object = ON3, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot
ON3 <- RunUMAP(object = ON3, dims = 1:31)
DimPlot(object = ON3, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Run ACMT Function for Annotating
CMT_ON3 <- MouseAutomatedClusterMarkerTable(Seurat_Object = ON3)
write.csv(CMT_ON3[[1]], file = "~/Desktop/CMT_ON3.csv") 
View(CMT_ON3[[1]])

#Master FeaturePlot
FeaturePlot(object = ON3, features = c("Ptprc", "Krt18", "Krt19", "Epcam", "Cd3e", "Cd4", "Cd8a", "Il7r", "Itgam", "Itgax", "Cd14", "Cd33", "Ncam1", "Fcgr3", "Nkg7", "Klrk1", "Cd19", "Cd22", "Ms4a1", "Tnfrsf17", "Pdgfra", "Pdgfrb", "Acta2", "Cdh11"), pt.size = 0.5)

#Visualize Labelled Clusters
new.cluster.ids <- CMT_ON3[[2]]
names(x = new.cluster.ids) <- sort(as.numeric(levels(x = ON3)), decreasing = F)
ON3 <- RenameIdents(object = ON3, new.cluster.ids)
ON3[["Auto_Labelled_Clusters"]] <- ON3@active.ident
DimPlot(object = ON3, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON3@active.ident)

#Identity Switches
Idents(object = ON3) <- "seurat_clusters"
Idents(object = ON3) <- "Auto_Labelled_Clusters"

# save Seurat object
save(ON3, file = 'ON3.RData')

#Determine n Cutoffs
VlnPlot(ON3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ON3_Fil <- subset(x = ON3, subset = nCount_RNA < 60000 )
ON3_Fil <- subset(x = ON3, subset = nCount_RNA < 60000 & percent.mt < 15 )
ON3_Fil <- subset(x = ON3, subset = nCount_RNA > 1000 & nCount_RNA < 60000 & percent.mt < 15)
VlnPlot(ON3_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON3_Fil@active.ident)

#Normalize Data
ON3_Fil <- NormalizeData(object = ON3_Fil, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
ON3_Fil <- FindVariableFeatures(object = ON3_Fil, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
ON3_Fil <- CellCycleScoring(ON3_Fil, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ON3_Fil$CC.Difference <- ON3_Fil$S.Score - ON3_Fil$G2M.Score

#Scale Data
all.genes <- rownames(x = ON3_Fil)
ON3_Fil <- ScaleData(ON3_Fil, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
ON3_Fil <- RunPCA(object = ON3_Fil, features = VariableFeatures(object = ON3_Fil))
stdev <- ON3_Fil@reductions$pca@stdev
var <- stdev^2
sum(var[1:30])/ sum(var)

#Find Neighbors + Find CLusters
ON3_Fil <- FindNeighbors(object = ON3_Fil, dims = 1:30)
ON3_Fil <- FindClusters(object = ON3_Fil, resolution = 1.9)

#Run UMAP and get unlabelled cluster UMAP and violin plot
ON3_Fil <- RunUMAP(object = ON3_Fil, dims = 1:30)
DimPlot(object = ON3_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON3_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Run ACMT Function for Annotating
CMT_ON3_Fil <- MouseAutomatedClusterMarkerTable(Seurat_Object = ON3_Fil)
write.csv(CMT_ON3_Fil[[1]], file = "~/Desktop/CMT_ON3_Fil.csv") 
View(CMT_ON3_Fil[[1]])

#Master FeaturePlot
FeaturePlot(object = ON3_Fil, features = c("Ptprc", "Krt18", "Krt19", "Epcam", "Cd3e", "Cd4", "Cd8a", "Il7r", "Itgam", "Itgax", "Cd14", "Cd33", "Ncam1", "Fcgr3", "Nkg7", "Klrk1", "Cd19", "Cd22", "Ms4a1", "Tnfrsf17", "Pdgfra", "Pdgfrb", "Acta2", "Cdh11"), pt.size = 0.5)

#Visualize Labelled Clusters
Idents(object = ON3_Fil) <- "seurat_clusters"
new.cluster.ids <- CMT_ON3_Fil[[2]]
names(x = new.cluster.ids) <- sort(as.numeric(levels(x = ON3_Fil)), decreasing = F)
ON3_Fil <- RenameIdents(object = ON3_Fil, new.cluster.ids)
ON3_Fil[["Auto_Labelled_Clusters"]] <- ON3_Fil@active.ident
DimPlot(object = ON3_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON3_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON3_Fil@active.ident)

#Final Labelling
Idents(object = ON3_Fil) <- "seurat_clusters"
ON3_Fil <- RenameIdents(ON3_Fil, 
                        "0" = "Fibroblasts",
                        "1" = "RBC",
                        "2" = "B Cell",
                        "3" = "Fibroblasts",
                        "4" = "Treg",
                        "5" = "T Cell",
                        "6" = "Fibroblasts",
                        "7" = "T Cell",
                        "8" = "RBC",
                        "9" = "Myeloid",
                        "10" = "T Cell",
                        "11" = "NK",
                        "12" = "B Cell",
                        "13" = "Myeloid",
                        "14" = "Acinar",
                        "15" = "RBC",
                        "16" = "Fibroblasts",
                        "17" = "Fibroblasts", 
                        "18" = "T Cell", 
                        "19" = "Fibroblasts", 
                        "20" = "RBC", 
                        "21" = "Epithelial", 
                        "22" = "Myeloid", 
                        "23" = "Myeloid", 
                        "24" = "Myeloid", 
                        "25" = "B Cell")
ON3_Fil[["Manual_Labelled_Clusters"]] <- ON3_Fil@active.ident
DimPlot(object = ON3_Fil, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON3_Fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Identity Switches
Idents(object = ON3_Fil) <- "seurat_clusters"
Idents(object = ON3_Fil) <- "Auto_Labelled_Clusters"
Idents(object = ON3_Fil) <- "Manual_Labelled_Clusters"

# save Seurat object
save(ON3_Fil, file = 'ON3_Fil.RData')

#Preparation of final object ON3 for analysis

#Filter Out RBC
ON3_Final <- subset(x = ON3_Fil, idents = c("Fibroblasts","B Cell", "Treg", "T Cell","Myeloid", "NK" ,"Acinar", "Epithelial"))

#Calculate MT %
ON3_Final[["percent.mt"]] <- PercentageFeatureSet(object = ON3_Final, pattern = "^mt-")

#Normalize Data
ON3_Final <- NormalizeData(object = ON3_Final, normalization.method = "LogNormalize", scale.factor = 10000)

#Find Variable Genes
ON3_Final <- FindVariableFeatures(object = ON3_Final, selection.method = "vst", nfeatures = 2000)

#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- str_to_sentence(s.genes)
g2m.genes <- str_to_sentence(g2m.genes)
ON3_Final <- CellCycleScoring(ON3_Final, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ON3_Final$CC.Difference <- ON3_Final$S.Score - ON3_Final$G2M.Score

#Scale Data
all.genes <- rownames(x = ON3_Final)
ON3_Final <- ScaleData(ON3_Final, vars.to.regress = c("nCount_RNA","CC.Difference"), features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance
ON3_Final <- RunPCA(object = ON3_Final, features = VariableFeatures(object = ON3_Final))
stdev <- ON3_Final@reductions$pca@stdev
var <- stdev^2
sum(var[1:31])/ sum(var)

#Find Neighbors, Find CLusters
ON3_Final <- FindNeighbors(object = ON3_Final, reduction = "pca", dims = 1:31)
ON3_Final <- FindClusters(object = ON3_Final, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot
ON3_Final <- RunUMAP(ON3_Final, reduction = "pca", dims = 1:31)
DimPlot(ON3_Final, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(object = ON3_Final, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Run ACMT Function for Annotating
CMT_ON3_Final <- MouseAutomatedClusterMarkerTable(Seurat_Object = ON3_Final)
write.csv(CMT_ON3_Final[[1]], file = "~/Desktop/CMT_ON3_Final.csv") 
View(CMT_ON3_Final[[1]])

#Visualize Labelled Clusters
Idents(object = ON3_Final) <- "seurat_clusters"
new.cluster.ids <- CMT_ON3_Final[[2]]
names(x = new.cluster.ids) <- sort(as.numeric(levels(x = ON3_Final)), decreasing = F)
ON3_Final <- RenameIdents(object = ON3_Final, new.cluster.ids)
ON3_Final[["Auto_Labelled_Clusters"]] <- ON3_Final@active.ident
DimPlot(object = ON3_Final, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON3_Final, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
table(ON3_Final@active.ident)

#Final Labelling
Idents(object = ON3_Final) <- "seurat_clusters"
ON3_Final <- RenameIdents(ON3_Final, 
                          "0" = "Fibroblasts",
                          "1" = "B Cells",
                          "2" = "Fibroblasts",
                          "3" = "Tregs",
                          "4" = "CD8+ T Cells",
                          "5" = "Fibroblasts",
                          "6" = "CD4+ T Cells",
                          "7" = "Myeloid 1",
                          "8" = "NK Cells",
                          "9" = "CD4+ T Cells",
                          "10" = "B Cells",
                          "11" = "Myeloid 2",
                          "12" = "Acinar",
                          "13" = "Fibroblasts",
                          "14" = "Fibroblasts",
                          "15" = "CD4- CD8- T Cells",
                          "16" = "Myeloid 3",
                          "17" = "Myeloid 4",
                          "18" = "Epithelial",
                          "19" = "B Cells")
ON3_Final[["Final_Labeled_Clusters"]] <- ON3_Final@active.ident
DimPlot(object = ON3_Final, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(ON3_Final, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Identity Switches
Idents(object = ON3_Final) <- "seurat_clusters"
Idents(object = ON3_Fil) <- "Auto_Labelled_Clusters"
Idents(object = ON3_Final) <- "Final_Labeled_Clusters"

# save Seurat object
save(ON3_Final, file = 'ON3_Final.RData')

#Figures

#Figure 1D

#mPDA
DimPlot(object = Mouse_On_Merge_Fil, reduction = "umap", label = F, pt.size = .5, cols = c("Epithelial" = "firebrick1", "Fibroblasts" = "royalblue1", "Acinar" = "seashell2", "Myeloid 1" = "orange", 
                                                                                           "Myeloid 2" = "saddlebrown", "Myeloid 3" = "goldenrod3", "CD4+ T Cells" = "olivedrab1", 
                                                                                           "Tregs" = "darkolivegreen", "CD8+ T Cells" = "forestgreen", 
                                                                                          "NK Cells" = "magenta4", "B Cells" = "gray0", "Endothelial" ))
FeaturePlot(object = Mouse_On_Merge_Fil, features = c("Cd4","Ctla4","Foxp3","Il2ra"), pt.size = 0.5, order = T, cols = c("gainsboro","firebrick1"))

#mPanIN
DimPlot(object = ON3_Final, reduction = "umap", label = F, pt.size = .5, cols = c("Epithelial" = "firebrick1", "Fibroblasts" = "royalblue1", "Acinar" = "seashell2", "Myeloid 1" = "orange", 
                                                                                  "Myeloid 2" = "saddlebrown", "Myeloid 3" = "goldenrod3", "Myeloid 4" = "plum3", "CD4+ T Cells" = "olivedrab1", 
                                                                                  "Tregs" = "darkolivegreen", "CD8+ T Cells" = "forestgreen", "CD4- CD8- T Cells" = "lightseagreen", "NK Cells" = "magenta4", "B Cells" = "gray0"))
FeaturePlot(object = ON3_Final, features = c("Cd4","Ctla4","Foxp3","Il2ra"), pt.size = 0.5, order = T, cols = c("gainsboro","firebrick1"))


#Figure 3D

#mPDA
FeaturePlot(object = Mouse_On_Merge_Fil, features = c("Tgfb1","Tgfbr1","Tgfbr2","Tgfbr3"), pt.size = 0.5, order = T, cols = c("gainsboro","firebrick1"))

#mPanIN
FeaturePlot(object = ON3_Final, features = c("Tgfb1","Tgfbr1","Tgfbr2","Tgfbr3"), pt.size = 0.5, order = T, cols = c("gainsboro","firebrick1"))


#Figure 6G

#mPDA
FeaturePlot(object = Mouse_On_Merge_Fil, features = c("Ccr1"), pt.size = 0.5, order = T, cols = c("gainsboro","firebrick1"))

#mPanIN
FeaturePlot(object = ON3_Final, features = c("Ccr1"), pt.size = 0.5, order = T, cols = c("gainsboro","firebrick1"))

#Figure S1F
ON3_order <- rev(c("Epithelial", "Fibroblasts", "Acinar", "NK Cells", "B Cells", "CD4+ T Cells", "CD8+ T Cells", 
                   "CD4- CD8- T Cells", "Tregs", "Myeloid 1", "Myeloid 2", "Myeloid 3", "Myeloid 4")) 
ON3_Final@active.ident <- factor(ON3_Final@active.ident, levels = ON3_order)
DotPlot(ON3_Final, features = rev(c("Krt19", "Pdgfra","Try4", "Ptprc", "Nkg7", "Cd79a", "Cd3e", "Cd4", 
                                          "Cd8a", "Foxp3","Cd68", "Itgax" ,"Fcgr3","Cd14","Itgam","Adgre1","S100a8",
                                          "Ly6g","Mrc1", "H2-Eb1", "Batf3","Itgae","Clec9a", "Ly6c2")), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


#Figure S7D

#mPDA
FeaturePlot(object = Mouse_On_Merge_Fil, features = c("Il1a","Il1r1"), pt.size = 0.5, order = T, cols = c("gainsboro","firebrick1"))

#mPanIN
FeaturePlot(object = ON3_Final, features = c("Il1a","Il1r1"), pt.size = 0.5, order = T, cols = c("gainsboro","firebrick1"))

# save Seurat object
save(Mouse_On_Merge_Fil, file = 'Mouse_On_Merge_Fil.RData')
save(ON3_Final, file = 'ON3_Final.RData')