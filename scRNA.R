install.packages("infercnv","AnnotationHub", "patchwork")
install.packages("patchwork", dependencies = TRUE, INSTALL_opts = '--no-lock')


install.packages("ensembldb", dependencies = TRUE, INSTALL_opts = '--no-lock')

install.packages("infercnv", dependencies = TRUE, INSTALL_opts = '--no-lock')

library(tidyverse)
library(Matrix)
library(DropletUtils)
library(Seurat)
library(infercnv)
library(AnnotationHub)
library(ensembldb)
library(patchwork)
library(cowplot)
library(dplyr)


#bioconductor 

LGG_Astro1 <- Read10X(data.dir = "scRNA2/LGG_Astro3/", gene.column = 1)

LGG_Astro1 <- CreateSeuratObject(counts = LGG_Astro1, min.cells = 3, min.features = 200, project="Diaz")


object_1 <- Read10X(data.dir = "scRNA/")

object_1 <- CreateSeuratObject(counts = object_1, min.cells = 3, min.features = 200, project="Diaz")


mito.features_object1 <- grep(pattern = "^MT-", x=rownames(x=object_1), value=T)
percent.mito_object1 <- Matrix::colSums(x = GetAssayData(object = object_1, slot="counts")[mito.features_object1,]) / Matrix::colSums(x = GetAssayData(object = object_1, slot = "counts"))
object_1[["percent.mito"]] <- percent.mito_object1
VlnPlot(object = object_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 


ggplot(object_1@meta.data, aes(y=`nFeature_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 2000,col="red") +
  geom_hline(yintercept = 6500,col="red")

ggplot(object_1@meta.data, aes(y=`nCount_RNA`, x=orig.ident)) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 4500,col="red") +
  geom_hline(yintercept = 22500,col="red")


mito.features_LGG_Astro1 <- grep(pattern = "^MT-", x=rownames(x=LGG_Astro1), value=T)
percent.mito_LGG_Astro1 <- Matrix::colSums(x = GetAssayData(object = LGG_Astro1, slot="counts")[mito.features_LGG_Astro1,]) / Matrix::colSums(x = GetAssayData(object = LGG_Astro1, slot = "counts"))

LGG_Astro1[["percent.mito"]] <- percent.mito_LGG_Astro1
VlnPlot(object = LGG_Astro1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01, group.by = "orig.ident") 




ggplot(LGG_Astro1@meta.data, aes(y=`nFeature_RNA`, x="orig.ident")) +
  geom_jitter(cex=0.01) +
  geom_hline(yintercept = 1500,col="red") +
  geom_hline(yintercept = 6000,col="red")

ggplot(LGG_Astro1@meta.data, aes(y=`nCount_RNA`, x="orig.ident")) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 4500,col="red") +
  geom_hline(yintercept = 23000,col="red")

ggplot(LGG_Astro1@meta.data, aes(y=`percent.mito`, x="nUMI")) +
  geom_jitter(cex=0.01)  +
  geom_hline(yintercept = 0,col="red") +
  geom_hline(yintercept = 0.2,col="red")

LGG_Astro1 <- subset(x = LGG_Astro1, subset =
                     nFeature_RNA > 1500 &
                     nFeature_RNA < 6000 &
                     nCount_RNA > 4500 &
                     nCount_RNA < 23000 &
                     percent.mito < 0.2)

LGG_Astro1 <- NormalizeData(object = LGG_Astro1, normalization.method = "LogNormalize", scale.factor = 1e4)

LGG_Astro1 <- FindVariableFeatures(object = LGG_Astro1, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(LGG_Astro1), 10)


plot1 <- VariableFeaturePlot(LGG_Astro1) + theme(legend.position = 'bottom')

plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot2


allgenesAstro1 <- rownames(LGG_Astro1)
LGG_Astro1 <- ScaleData(LGG_Astro1, features = allgenesAstro1)

LGG_Astro1 <- RunPCA(LGG_Astro1, features = VariableFeatures(object = LGG_Astro1))



print(LGG_Astro1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(LGG_Astro1, dims = 1:2, reduction = "pca")
DimPlot(LGG_Astro1, reduction = "pca")


ElbowPlot(LGG_Astro1, ndims = 45)

d <- 25

#use the first 40 PCs to generate clusters
LGG_Astro1 <- FindNeighbors(LGG_Astro1, dims = 1:40)
LGG_Astro1 <- FindClusters(LGG_Astro1, resolution = 1)
                           #, algorithm=1,dims.use = 1:10
head(Idents(LGG_Astro1), 20)




LGG_Astro1 <- RunUMAP(LGG_Astro1, reduction = "pca", dims = 1:40)

#LGG_Astro1 <- RunTSNE(object = LGG_Astro1, do.fast = TRUE)


#TSNEPlot(object = LGG_Astro1)

UMAPPlot(object = LGG_Astro1)

DimPlot(LGG_Astro1, reduction = "umap", label = T, label.size = 6)


n_cells <- FetchData(LGG_Astro1, vars = c("ident", "orig.ident")) %>% 
  dplyr::count(ident, orig.ident) %>% 
  tidyr::spread(ident, n)


DimPlot(LGG_Astro1, label = T, split.by= "Phase") 


LGG_Astro1@meta.data

metrics <- c("nUMI","nGene", "S.score", "G2M.score", "mitoRatio")

FeaturePlot(LGG_Astro1, reduction = "umap", features = ("percent.mito"), pt.size = 0.4, sort.cell=T, min.cutoff = 'q10', label = T)


cluster1.markers <- FindMarkers(object = LGG_Astro1, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))
  

cluster0.markers <- FindMarkers(object = LGG_Astro1, ident.1 = 0, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)


cluster1.markers <- FindMarkers(object = LGG_Astro1, ident.1 = 1, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)

cluster2.markers <- FindMarkers(object = LGG_Astro1, ident.1 = 2, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)

cluster3.markers <- FindMarkers(object = LGG_Astro1, ident.1 = 3, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)
cluster4.markers <- FindMarkers(object = LGG_Astro1, ident.1 = 4, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)




# find markers for every cluster compared to all remaining cells, report only the positive ones   
LGG_Astro1.markers <- FindAllMarkers(object = LGG_Astro1, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)

LGG_Astro1.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)


top10 <- LGG_Astro1.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)


DoHeatmap(object = LGG_Astro1)
?DoHeatmap


VlnPlot(object = LGG_Astro1, features = c("EGFR", "TMEM144", "GJA1", "CD74","RP11-89N17.4"))


FeaturePlot(object = LGG_Astro1, )

FeaturePlot(object = LGG_Astro1, features = c("EGFR", "VIM")) # Tumor

FeaturePlot(object = LGG_Astro1, features = "VIM") # Tumor/MES


FeaturePlot(object = LGG_Astro1, features = "PCNA") # G1/S


#see if you find proteins with a high correlation in the same cell type 



#https://github.com/yhoogstrate/scRNA-glioma-project/blob/master/scripts/scRNA_analyses_Diaz.R



#' SF12017	Astrocytoma	44	M	IDH1R132S mutant	TP53 mutation	chr1-, chr7q+, chr19q-	10x	single cell, snATAC	N
#' SF11964	Glioblastoma	64	M	IDH1R132H mutant	ATRX deletion	chr1q-, chr6-, chr14-, chr18-, chr19q-	10x	single cell, snATAC	N
#LGG_Astro3' SF11136	Astrocytoma	NA	NA	IDH1R132C mutant	ATRX stopgain	chr13-, chr17+,	10x	single cell	Y

#######################################

LGG_Astro1 <- Read10X(data.dir = "scRNA2/LGG_Astro3/", gene.column = 1)

LGG_Astro1 <- CreateSeuratObject(counts = LGG_Astro1, min.cells = 3, min.features = 200, project="Diaz")
LGG_Astro1@meta.data

# Add number of genes per UMI for each cell to metadata
LGG_Astro1$log10GenesPerUMI <- log10(LGG_Astro1$nFeature_RNA) / log10(LGG_Astro1$nCount_RNA)


# Compute percent mito ratio (^MT works for human gene names )
LGG_Astro1$mitoRatio <- PercentageFeatureSet(object = LGG_Astro1, pattern = "^MT-")
LGG_Astro1$mitoRatio <- LGG_Astro1@meta.data$mitoRatio / 100

LGG_Astro1_meta <- LGG_Astro1@meta.data

LGG_Astro1_meta$cells <- rownames(LGG_Astro1_meta)
LGG_Astro1_meta <- LGG_Astro1_meta %>% 
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Create sample column
LGG_Astro1_meta$sample <- "1"
LGG_Astro1_meta$sample[which(str_detect(LGG_Astro1_meta$cells, "1"))] <- "1"
LGG_Astro1_meta$sample[which(str_detect(LGG_Astro1_meta$cells, "2"))] <- "2"

LGG_Astro1@meta.data <- LGG_Astro1_meta

#nr of cell counts per sample
LGG_Astro1_meta %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


#nr of UMIs/transcripts per cell
LGG_Astro1_meta %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

#distribution of genes in histogram
LGG_Astro1_meta %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

#vind veel genen in weinig cellen

LGG_Astro1_meta %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

#look at mitochondrial reads as high mitochondrial reads can indicate damaged or dying cells

LGG_Astro1_meta %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
LGG_Astro1_meta %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

LGG_Astro1_meta %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

#FILTERING
LGG_Astro1 <- subset(x = LGG_Astro1, 
                          subset= (nUMI >= 500) & #nUMI
                            (nGene >= 250) & #nGene
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))


#remove genes that are expressed in less than 10 cells
counts <- GetAssayData(object = LGG_Astro1, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
LGG_Astro1 <- CreateSeuratObject(filtered_counts, meta.data = LGG_Astro1@meta.data)
LGG_Astro1_meta <- LGG_Astro1@meta.data

rm(counts, nonzero, keep_genes, filtered_counts)


# Identify the most variable genes
LGG_Astro1 <- FindVariableFeatures(LGG_Astro1, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

LGG_Astro1 <- NormalizeData(object = LGG_Astro1)

##Check if PCs separate the cell types well 

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(LGG_Astro1, 
                     vars = columns)


# Identify the most variable genes
LGG_Astro1 <- FindVariableFeatures(LGG_Astro1, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
LGG_Astro1 <- ScaleData(LGG_Astro1)

# Perform PCA
LGG_Astro1 <- RunPCA(LGG_Astro1)

# Plot the PCA colored by cell cycle phase
DimPlot(LGG_Astro1,
        reduction = "pca")

#Determine k-nearest neighbour graph
LGG_Astro1 <- FindNeighbors(LGG_Astro1, dims = 1:40)

#Determine clusters
LGG_Astro1 <- FindClusters(LGG_Astro1, resolution = 1.5)

#Run and plot UMAP
LGG_Astro1 <- RunUMAP(LGG_Astro1,
                      reduction = "pca",
                      dims = 1:40)
DimPlot(LGG_Astro1,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

#CHECK FOR CLUSTERING ARTIFACTS

n_cells <- FetchData(LGG_Astro1, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)


DimPlot(LGG_Astro1, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()

metrics <-  c("nUMI", "nGene", "mitoRatio")

FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


##colour clusters with PCs
# Adding cluster label to center of cluster on UMAP

umap_label <- FetchData(LGG_Astro1, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))


# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)


# Examine PCA results 
print(LGG_Astro1[["pca"]], dims = 1:5, nfeatures = 5)

#astrocyte markers
pdf("Astrocyte_markers.pdf", height = 7, width =7)
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("GFAP","SLC1A2", "ACSBG1", "SLC1A3","GJA1","AQP4"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
#oligodendrocyte markers
pdf("Oligo_markers.pdf", height = 7, width =7)

FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("MBP","MOG", "MAG", "TMEM144", "SOX10"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

dev.off()
#OPC
pdf("OPC_tumor_markers.pdf", height = 5, width = 7)
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("NEU4", "TNR","EPN2", "OLIG1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

dev.off()
#Tumor
pdf("Tumor_markers.pdf",height = 5, width = 7)
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("ETV1","MYC", "EGFR", "ETNPLL", "OLIG1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
#Macrophages, $ TAM 
pdf("TAM_markers.pdf", height = 7, width = 7)
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("CD163", "CD14", "ITGB2", "C1QC", "CD68", "NAAA"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()
#T-cells not found
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("TRBC2", "CD3D", "CD3G", "CD3E", "IL7R", "GZMK", "ICOS", "GZMA"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#Endothelial cells 
pdf("Endothelial_cells.pdf", height = 7, width = 7)
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("PECAM1", "EGFL7", "ID3", "FLT1", "GNG11", "MCAM"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()


FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("AURKB"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


#G1/S
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("PCNA","KIAA0101"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#G2/M
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("CCNB1","CDC20", "TMPO"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)



LGG_Astro1.markers <- FindAllMarkers(object = LGG_Astro1, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)     


top10 <- LGG_Astro1.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)


cluster0.markers <- FindMarkers(object = LGG_Astro1, ident.1 = 0, ident.2 = 10,  
                                test.use = "roc", only.pos = TRUE)

#protsignup
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("PCNA","MCM7", "NFIX", "HMGN2", "KPNA2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#protsigndown
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("PLP1","PLLP", "MBP", "MAG", "TNR"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


VlnPlot(LGG_Astro1, features = c("MBP","MOG", "MAG", "PLP1", "PLLP", "TNR", "NEFL"))


#protsignup
FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("PCNA","MCM7", "NFIX", "HMGN2", "KPNA2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


VlnPlot(LGG_Astro1, features = c("PLP1","PLLP", "MBP", "MAG"), y.max = 6)

cl01267 <- FindMarkers(LGG_Astro1, ident.1 = c(0,1,2,6,7)) 

cl01267 <- cl01267 %>% 
  rownames_to_column("gene")

cl458 <- FindMarkers(LGG_Astro1, ident.1 = c(4,5,8))

cl458 <- cl458 %>% 
  rownames_to_column("gene")

cl3 <- FindMarkers(LGG_Astro1, ident.1 = 3)

cl3 <- cl3 %>% 
  rownames_to_column("gene")

LGG_Astro1_labeled <- RenameIdents(LGG_Astro1, 
                                   "0" = "Oligodendrocytes",
                                   "1" = "Oligodendrocytes",
                                   "2" = "Oligodendrocytes",
                                   "6" = "Oligodendrocytes",
                                   "7" = "Oligodendrocytes",
                                   "4" = "Tumor",
                                   "5" = "Tumor",
                                   "8" = "Tumor",
                                   "9" = "Oligodendrocytes",
                                   "12" = "OPC/Tumor",
                                   "3" = "Astrocytes",
                                   "10" = "TAM",
                                   "11" = "Endothelial cells"
                                   )
pdf("Cell", height = 7, width = 9)
DimPlot(object = LGG_Astro1,
        reduction = "umap",
        label = T,
        pt.size = 0.5,
        repel = T,
        label.size = 5)


dev.off()
write_rds(LGG_Astro1, file = "LGG_Astro1.rds")

write_rds(LGG_Astro1_labeled, file = "LGG_Astro1_labeled.rds")

save(LGG_Astro1.markers, file = "LGG_Astro1.markers.Rdata")



# A :: LGG?, WHO grade IV, idh1: R132H ----

# A :: GSM2758472_PJ016 :: ??  ----

#' age: 49

#' gender: Female

#' location: right frontal

#' diagnosis: Glioblastoma, WHO grade IV

#' idh1 status: R132H

#' egfr status: not amplified
#' 
#' # diagnosis: Glioblastoma, WHO grade IV - idh1 status: R132H
a <- read.delim("data/GSE103224_Yuan/GSM2758471_PJ016.filtered.matrix.txt", stringsAsFactors = F,header=F) %>%
  `colnames<-`(c("ENS", "HGNC", paste0("GSM2758471_PJ016.cell.",(1:(ncol(.)-2))+2) ))





##########CHECK HIGH CORRELATING GENES IN SCRNA DATA


FeaturePlot(LGG_Astro1_labeled, 
            reduction = "umap", 
            features = c("ANXA2","HSPB1", "VIM", "ANXA1", "S100A6", "TMPO", "VAT1", "RELA", "FLNA", "LUM", "PDIA4", "SCIN", "LAP3"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


VlnPlot(LGG_Astro1_labeled, features = c("ANXA2","HSPB1", "VIM", "ANXA1", "S100A6", "TMPO", "VAT1", "RELA", "FLNA", "LUM", "PDIA4", "SCIN", "LAP3"), y.max = 5)


VlnPlot(LGG_Astro1_labeled, features = clustersprotRNA2[1:13], y.max = 3)


DotPlot(LGG_Astro1_labeled, features = clustersprotRNA3)+ RotatedAxis()


DotPlot(LGG_Astro1_labeled, features = clustersprotRNA2)+ RotatedAxis()

DotPlot(LGG_Astro1_labeled, features = clustersprotRNA1)+ RotatedAxis()

FeaturePlot(LGG_Astro1, 
            reduction = "umap", 
            features = c("EGFR", "MYC"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


#PLP1

FeaturePlot(LGG_Astro1_labeled, 
            reduction = "umap", 
            features = c("PLP1","PLLP", "MBP", "MAG"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

pdf("oligo_down_prot.pdf", width = 8, height = 8)



VlnPlot(LGG_Astro1_labeled, features = c("PLP1","PLLP", "MBP", "MAG"), ncol = 2, y.max = 6)

dev.off()

DotPlot(LGG_Astro1_labeled, features = corrnatot$gene[151:200])+ RotatedAxis()


