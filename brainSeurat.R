  
#1 - Load libraries
  library(dplyr)
  library(Seurat)
  library(patchwork)

#2 - Load dataset
  
  dir() # can use setwd() to find your file too. I use (dir) to find my file names.
  brain.data <- Read10X("filtered_feature_bc_matrix")
  brain <- CreateSeuratObject(counts = brain.data, project = "MBD", min.cells = 3, min.features = 200)
  head(brain) #view column names
  
  # Now that we have loaded in our data, we need to clean our data. This is what Quality Control does,
  # it will essentially remove obviously bad data points. Since scRNA-seq data is so sparse(meaning there
  # are many 0 values), we MUST run QC to make sure results are reliable.
    # Be sure to pay attention to the columns and there values throughout the project, you can do so with 
    # head(brain)
  
#3 -  QC and selecting cells for further analysis
  
  head(brain) #preview columns
  brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-") #add percent.mt column
  head(brain)
  
  # One common mistake is the pattern = "" portion. There are many different types of patterns that
  # are used within the scRNA-seq field. For example, in the Seurat dataset, they use "^MT-"
  # because it is human data.
  # Since this dataset uses mouse data, we need to use "^mt-"
  
  VlnPlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # Take a look at how sparse the data is. We have a ton of mitochondria that is not real
  # data (also known as noise). We basically need to control what data does and does not get processed, as this can 
  # cost a lot of time and money.
  # 15-20 percent can be a decent spot to cut off. Each scenario is case by case, 
  # however use the 5-10% range to ensure live cells. There can be lots of money lost with just throwing out
  # large amount of data.
 
  brainSub15 <- subset(brain, subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & percent.mt < 15) 
  # Clean up the data with these conditions
  VlnPlot(brainSub15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # Make a violin plot with the newly cleaned data
  
  
  # Notice the difference between the plots. There is less noise and the data is more focused on actual results.
  
#4 -  Normalize the Data
  
  # Before we advance in scRNA-seq, we must normalize our large amount of data. We do this because there is too much variability
  # that ranges from high values to low values when looking at gene expression. Take a look at the before and after code. 
  
  VlnPlot(brain,features="Hbb-bs") # Before Normalization, looking at one gene Hbb-bs
  brain <- NormalizeData(brainSub15, normalization.method = "LogNormalize", scale.factor = 10000) 
  VlnPlot(brain,features="Hbb-bs") # After Normalization
  
#5 - Identification of highly variable features
  
  # Now that we have cleaned and normalized our "good" data, we need to identify which are the most significant
  # in scRNAseq, by most significant, we mean most variable and most expressd.
  # This can be displayed by the bottom code.

  brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(brain), 10)
  plot1 <- VariableFeaturePlot(brain)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 
  plot2

  # Line 68 will show you the purpose of this part. We want to look for the top
  # 10 significant genes in this study, so it shows the before and after (left and right)
  # This is great when trying to get a general idea on what kinds of genes are interacting
  # or influencing the data set
  
#6 - Scaling the data
  
  # Next up, we need to scale our data down further to perform Linear dimensional
  # Reduction. This is due to a large dataset, so we need to reduce and dimensionalize 
  # the data linearly. That is; treat each point equally and condense them all
  # down together.
  
  all.genes <- rownames(brain)
  brain <- ScaleData(brain, features = all.genes)
  
#7 -Perform linear dimensional reduction
  
  brain <- RunPCA(brain, features = VariableFeatures(object = brain))
  VizDimLoadings(brain, dims = 1, reduction = "pca")
  print(x = brain[["pca"]], dims=1, nfeatures=10)
  
  # We can even look at the top genes in each Principal Component
  # This is just another step in verification checking
  
  DimHeatmap(brain, dims = 1, cells = 500, balanced = TRUE)
  
  # Also need to look at our heat maps to look for up and down regulation,
  # As well as correlation strength.There will be more correlation between
  # each gene in the corresponding PC according to the amount of Black Lines.
  # TLDR: More black lines = less correlation
  
#8 - Determine the dimensionality of the dataset
  
  # Now we need to figure out how many of these Principal Components
  # will be clustered in the next step. We want to look at the elbow plot
  # And cut off the data where it starts to run paralell to the X axis.
  # That will remove non-variable gene expression data,
  # We are looking for only the highly variable genes
  
  brain <- JackStraw(brain, num.replicate = 100)
  brain <- ScoreJackStraw(brain, dims = 1:20)
  ElbowPlot(brain,ndims = 50)
  
  # After running our elbow plot, we can clearly see that around 20-25 PC is
  # Where we want to stop our analysis at. Some would consider even doing up 
  # to 40 PCs, however it comes down to personal choice and how much time
  # you have. 

#9 - # Cluster the cells
  
  # Now we will group all of the cells together and paint them on a graph.
  # This is important because we want to group similar cell types together
  # And see what types of genes are being expressed by these cells
  # We picked 1:25 dimensions from the above step
  
  brain <- FindNeighbors(brain, dims = 1:25)
  brain <- FindClusters(brain, resolution = 0.5) 
  
#10 - Run non-linear dimensional reduction (UMAP/tSNE)
  
  # Here is the most important part of the code, our UMAP. 
  # Essentially it will portray the cells together that are related
  # and place other cells farther away based off gene expression
  
  brain <- RunUMAP(brain, dims = 1:25)
  DimPlot(brain, reduction = "umap")
  

#11 - # Finding differentially expressed features (cluster biomarkers)
  
  # This is the longest portion of code. It will also take the longest to compute.
  # While annotation of clusters is something that can be done manually, that will 
  # take more time and probably another set of code to explain in detail.
  # For now, we can ignore annotation and just focus on significant cells identified earlier
  
  neuron.markers <- FindAllMarkers(brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  neuron.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  FeaturePlot(object = brain, features = c("Hbb-bs"), pt.size=0.05)
  VlnPlot(brain,features="Hbb-bs")




















