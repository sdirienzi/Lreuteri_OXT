library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(gridExtra)
library(tidyverse)

#can load the rds file and skip to line 118
#this file also includes the scripts to get the OXT counts for Fig S1A

#Jejunum"
Jej.data<-read.table(file = "AdultJejunum2.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultJej<-CreateSeuratObject(counts=Jej.data,project = "AdultJej",min.cells = 3, min.features = 200)
AdultJej
#An object of class Seurat 
#16503 features across 5196 samples within 1 assay 
#Active assay: RNA (16503 features, 0 variable features)

AdultJej[["percent.mt"]] <- PercentageFeatureSet(AdultJej, pattern = "^MT-") 
head(AdultJej@meta.data,10)

VlnPlot(AdultJej,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultJej,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultJej,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultJejsub<-subset(AdultJej,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultJejsub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultJejsub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultJejsub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultJejsub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #50 OXT

AdultJejsub<-NormalizeData(AdultJejsub)
AdultJejsub <- FindVariableFeatures(AdultJejsub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultJejsub), 10)
variablegenes<-VariableFeatures(AdultJejsub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultJejsub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultJejsub)
AdultJejsub <- ScaleData(AdultJejsub, features = all.genes)
AdultJejsub <- RunPCA(AdultJejsub, features = VariableFeatures(object = AdultJejsub))
print(AdultJejsub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultJejsub, reduction ="pca")
ElbowPlot(AdultJejsub)  #use 20
AdultJejsub <- FindNeighbors(AdultJejsub, dims = 1:20)
AdultJejsub <- FindClusters(AdultJejsub, resolution = 0.7)
#Maximum modularity in 10 random starts: 0.8462
#Number of communities: 17
#Elapsed time: 0 seconds

AdultJejsub <- RunUMAP(AdultJejsub, dims = 1:20)
umapplotout<-DimPlot(AdultJejsub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultJejsub, features = c("OXT"))


saveRDS(AdultJejsub, file = "AdultJejsubumap.rds")

write.table(dotout$data,file="Jejtotalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Jej.markers <- FindAllMarkers(AdultJejsub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Jej.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultJejclustermarkers.txt",sep="\t",row.names=FALSE)

Jejvillus<-subset(AdultJejsub,subset=seurat_clusters=="0" |seurat_clusters=="1" |seurat_clusters=="2" |
                    seurat_clusters=="5" |seurat_clusters=="6"|seurat_clusters=="8"|seurat_clusters=="10"
                  |seurat_clusters=="11"|seurat_clusters=="14"|seurat_clusters=="15")

Jejvilluscells<-Cells(Jejvillus)

write.table(Jejvilluscells,file="Jejvilluscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Jejvillus<-NormalizeData(Jejvillus)
Jejvillus <- FindVariableFeatures(Jejvillus, selection.method = "vst", nfeatures = 5000)

all.genesvillus <- rownames(Jejvillus)
Jejvillus <- ScaleData(Jejvillus, features = all.genesvillus)
Jejvillus <- RunPCA(Jejvillus, features = VariableFeatures(object = Jejvillus))
print(Jejvillus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Jejvillus, reduction = "pca")
ElbowPlot(Jejvillus)  #use 13
Jejvillus <- FindNeighbors(Jejvillus, dims = 1:13)
Jejvillus <- FindClusters(Jejvillus, resolution = 0.6)
#Maximum modularity in 10 random starts: 0.7939
#Number of communities: 12
#Elapsed time: 0 seconds

OXTdata<-Jejvillus@assays$RNA@data["OXT",]
OXTcount<-Jejvillus@assays$RNA@counts["OXT",]


hist(OXTcount)
hist(OXTdata)
sum(OXTcount)
length(which(OXTcount!=0)) #40 with OXT, 3468 without cells
length(Jejvilluscells)

Jejvillusnormdata<-as.data.frame(Jejvillus@assays$RNA@data)
Jejvilluscounts<-as.data.frame(Jejvillus@assays$RNA@counts)

write.table(Jejvillusnormdata,file="Jejvillusnormdata.txt",sep="\t")
write.table(Jejvilluscounts,file="Jejvillusnormcounts.txt",sep="\t")

Jejvillus <- RunUMAP(Jejvillus, dims = 1:13)

saveRDS(Jejvillus, file = "FigS3gi_AdultJejvillusumap.rds")
Jejvillus<-readRDS(file = "FigS3gi_AdultJejvillusumap.rds")

Jejumapplotoutvillus<-DimPlot(Jejvillus, reduction = "umap",label=TRUE)

pdf(file="Jejvillusumap.pdf",width=5,height=4)
Jejumapplotoutvillus
dev.off()

pdf(file="JejOXTumap.pdf",width=5,height=4)
FeaturePlot(Jejvillus, features = c("OXT"))
dev.off()

vlnclusplot<-VlnPlot(Jejvillus, features = c("APOA1","APOA4","GSTA1","FABP1","REG1A","OLFM4",
                                             "ALDOB","APOB","IGHA1","JCHAIN","RBP2","SI",
                                             "RPL21","RPS27","GSTA2","RPS29","MUC2","ITLN1",
                                             "BEST4","CFTR","CHGA","GHRL","LYZ","DEFA5","OXT"),fill.by="ident",ncol=6)

pdf(file="Hanviolinplotcluster.pdf",height=16,width=24)
vlnclusplot
dev.off()


Jejvillus.markers <- FindAllMarkers(Jejvillus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Jejvillus.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="Jejvillusclustermarkers.txt",sep="\t",row.names=FALSE)

#try with SCTransform
Jejvillus<-SCTransform(Jejvillus, method = "glmGamPoi", vars.to.regress = "percent.mt",return.only.var.genes = FALSE,min_cells=1)
OXTdata<-Jejvillus@assays$SCT@data["OXT",]
OXTcounts<-Jejvillus@assays$SCT@counts["OXT",]

write.table(OXTdata, file="SCTdata/JEJ_OXT_SCTdata.txt",sep="\t")
write.table(OXTcounts, file="SCTcounts/JEJ_OXT_SCTcounts.txt",sep="\t")



#Duodenum
Duo.data<-read.table(file = "AdultDuodenum1.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultDuo<-CreateSeuratObject(counts=Duo.data,project = "AdultDuo",min.cells = 3, min.features = 200)
AdultDuo
#An object of class Seurat 
#14648 features across 3737 samples within 1 assay 
#Active assay: RNA (14648 features, 0 variable features)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
AdultDuo[["percent.mt"]] <- PercentageFeatureSet(AdultDuo, pattern = "^MT-") 
head(AdultDuo@meta.data,10)

VlnPlot(AdultDuo,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultDuo,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultDuo,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultDuosub<-subset(AdultDuo,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultDuosub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultDuosub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultDuosub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultDuosub@assays$RNA@data["OXT",]
hist(OXTcount)
sum(OXTcount) #8 OXT

AdultDuosub<-NormalizeData(AdultDuosub)

AdultDuosub <- FindVariableFeatures(AdultDuosub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultDuosub), 10)
variablegenes<-VariableFeatures(AdultDuosub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultDuosub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultDuosub)
AdultDuosub <- ScaleData(AdultDuosub, features = all.genes)
AdultDuosub <- RunPCA(AdultDuosub, features = VariableFeatures(object = AdultDuosub))
print(AdultDuosub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultDuosub, reduction ="pca")
ElbowPlot(AdultDuosub)  #use 15
AdultDuosub <- FindNeighbors(AdultDuosub, dims = 1:15)
AdultDuosub <- FindClusters(AdultDuosub, resolution = 0.7)
#Maximum modularity in 10 random starts: 0.8282
#Number of communities: 12
#Elapsed time: 0 seconds

AdultDuosub <- RunUMAP(AdultDuosub, dims = 1:15)
umapplotout<-DimPlot(AdultDuosub, reduction = "umap",label=TRUE)

umapplotout  #total umap
saveRDS(AdultDuosub, file = "AdultDuosubumap.rds")

write.table(dotout$data,file="Duototalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Duo.markers <- FindAllMarkers(AdultDuosub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Duo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultDuoclustermarkers.txt",sep="\t",row.names=FALSE)

Duovillus<-subset(AdultDuosub,subset=seurat_clusters=="2" |seurat_clusters=="3" |seurat_clusters=="4" |
                 seurat_clusters=="5" |seurat_clusters=="6" |seurat_clusters=="10")

Duovilluscells<-Cells(Duovillus)

write.table(Duovilluscells,file="Duovilluscells.txt",sep="\n",row.names = FALSE)

#analyze just the villus
Duovillus<-NormalizeData(Duovillus)
Duovillus <- FindVariableFeatures(Duovillus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Duovillus)
Duovillus <- ScaleData(Duovillus, features = all.genesvillus)
Duovillus <- RunPCA(Duovillus, features = VariableFeatures(object = Duovillus))
print(Duovillus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Duovillus, reduction = "pca")
ElbowPlot(Duovillus)  #use 10
Duovillus <- FindNeighbors(Duovillus, dims = 1:10)
Duovillus <- FindClusters(Duovillus, resolution = 0.8)
#Maximum modularity in 10 random starts: 0.7663
#Number of communities: 8
#Elapsed time: 0 seconds

OXTdata<-Duovillus@assays$RNA@data["OXT",]
OXTcount<-Duovillus@assays$RNA@counts["OXT",]

hist(OXTcount)
hist(OXTdata)
sum(OXTcount)#50 OXT
length(which(OXTcount!=0)) #6 with OXT, 1715 without cells
#export the counts and normalized data
Duovillusnormdata<-as.data.frame(Duovillus@assays$RNA@data)
Duovilluscounts<-as.data.frame(Duovillus@assays$RNA@counts)

write.table(Duovillusnormdata,file="Duovillusnormdata.txt",sep="\t")
write.table(Duovilluscounts,file="Duovillusnormcounts.txt",sep="\t")


Duovillus <- RunUMAP(Duovillus, dims = 1:10)
Duoumapplotoutvillus<-DimPlot(Duovillus, reduction = "umap",label=TRUE)

pdf(file="Duovillusumap.pdf",width=5,height=4)
Duoumapplotoutvillus
dev.off()

pdf(file="DuoOXTumap.pdf",width=5,height=4)
FeaturePlot(Duovillus, features = c("OXT"))
dev.off()

saveRDS(Duovillus, file = "AdultDuovillusumap.rds")
Duovillus<-readRDS(file = "AdultDuovillusumap.rds")


Duovillus<-SCTransform(Duovillus, method = "glmGamPoi", vars.to.regress = "percent.mt",return.only.var.genes = FALSE,min_cells=1)
OXTdata<-Duovillus@assays$SCT@data["OXT",]
OXTcounts<-Duovillus@assays$SCT@counts["OXT",]

write.table(OXTdata, file="SCTdata/DUO_OXT_SCTdata.txt",sep="\t")
write.table(OXTcounts, file="SCTcounts/DUO_OXT_SCTcounts.txt",sep="\t")


#Ileum
Ile.data<-read.table(file = "AdultIleum2.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultIle<-CreateSeuratObject(counts=Ile.data,project = "AdultIle",min.cells = 3, min.features = 200)
AdultIle
#An object of class Seurat 
#15966 features across 3132 samples within 1 assay 
#Active assay: RNA (15966 features, 0 variable features)

AdultIle[["percent.mt"]] <- PercentageFeatureSet(AdultIle, pattern = "^MT-") 
head(AdultIle@meta.data,10)

VlnPlot(AdultIle,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultIle,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultIle,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultIlesub<-subset(AdultIle,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultIlesub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultIlesub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultIlesub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultIlesub@assays$RNA@data["OXT",]
hist(OXTcount)
sum(OXTcount) #3 OXT

AdultIlesub<-NormalizeData(AdultIlesub)
AdultIlesub <- FindVariableFeatures(AdultIlesub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultIlesub), 10)
variablegenes<-VariableFeatures(AdultIlesub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultIlesub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultIlesub)
AdultIlesub <- ScaleData(AdultIlesub, features = all.genes)
AdultIlesub <- RunPCA(AdultIlesub, features = VariableFeatures(object = AdultIlesub))
print(AdultIlesub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultIlesub, reduction ="pca")
ElbowPlot(AdultIlesub)  #use 20
AdultIlesub <- FindNeighbors(AdultIlesub, dims = 1:20)
AdultIlesub <- FindClusters(AdultIlesub, resolution = 0.7)
Maximum modularity in 10 random starts: 0.9057
Number of communities: 15
Elapsed time: 0 seconds

AdultIlesub <- RunUMAP(AdultIlesub, dims = 1:20)
umapplotout<-DimPlot(AdultIlesub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultIlesub, features = c("OXT"))


saveRDS(AdultIlesub, file = "AdultIlesubumap.rds")

write.table(dotout$data,file="Iletotalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Ile.markers <- FindAllMarkers(AdultIlesub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Ile.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultIleclustermarkers.txt",sep="\t",row.names=FALSE)

Ilevillus<-subset(AdultIlesub,subset=seurat_clusters=="0" |seurat_clusters=="2" |seurat_clusters=="5" |
                    seurat_clusters=="9" |seurat_clusters=="12")

Ilevilluscells<-Cells(Ilevillus)

write.table(Ilevilluscells,file="Ilevilluscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Ilevillus<-NormalizeData(Ilevillus)
Ilevillus <- FindVariableFeatures(Ilevillus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Ilevillus)
Ilevillus <- ScaleData(Ilevillus, features = all.genesvillus)
Ilevillus <- RunPCA(Ilevillus, features = VariableFeatures(object = Ilevillus))
print(Ilevillus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Ilevillus, reduction = "pca")
ElbowPlot(Ilevillus)  #use 10
Ilevillus <- FindNeighbors(Ilevillus, dims = 1:10)
Ilevillus <- FindClusters(Ilevillus, resolution = 0.8)
Maximum modularity in 10 random starts: 0.7817
Number of communities: 9
Elapsed time: 0 seconds

OXTdata<-Ilevillus@assays$RNA@data["OXT",]
OXTcount<-Ilevillus@assays$RNA@counts["OXT",]

hist(OXTcount)
hist(OXTdata)
sum(OXTcount)#50 OXT
length(which(OXTcount==0)) #3 with OXT, 1300 without cells

Ilevillusnormdata<-as.data.frame(Ilevillus@assays$RNA@data)
Ilevilluscounts<-as.data.frame(Ilevillus@assays$RNA@counts)

write.table(Ilevillusnormdata,file="Ilevillusnormdata.txt",sep="\t")
write.table(Ilevilluscounts,file="Ilevillusnormcounts.txt",sep="\t")

Ilevillus <- RunUMAP(Ilevillus, dims = 1:10)
Ileumapplotoutvillus<-DimPlot(Ilevillus, reduction = "umap",label=TRUE)

pdf(file="Ilevillusumap.pdf",width=5,height=4)
Ileumapplotoutvillus
dev.off()

pdf(file="IleOXTumap.pdf",width=5,height=4)
FeaturePlot(Ilevillus, features = c("OXT"))
dev.off()

saveRDS(Ilevillus, file = "AdultIlevillusumap.rds")

Ilevillus<-readRDS( file = "AdultIlevillusumap.rds")

Ilevillus<-SCTransform(Ilevillus, method = "glmGamPoi", vars.to.regress = "percent.mt",return.only.var.genes = FALSE,min_cells=1)
OXTdata<-Ilevillus@assays$SCT@data["OXT",]
OXTcounts<-Ilevillus@assays$SCT@counts["OXT",]

write.table(OXTdata, file="SCTdata/ILE_OXT_SCTdata.txt",sep="\t")
write.table(OXTcounts, file="SCTcounts/ILE_OXT_SCTcounts.txt",sep="\t")


#epityphlon
Epi.data<-read.table(file = "AdultEpityphlon1.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultEpi<-CreateSeuratObject(counts=Epi.data,project = "AdultEpi",min.cells = 3, min.features = 200)
AdultEpi
#An object of class Seurat 
#14053 features across 3799 samples within 1 assay 
#Active assay: RNA (14053 features, 0 variable features) 

AdultEpi[["percent.mt"]] <- PercentageFeatureSet(AdultEpi, pattern = "^MT-") 
head(AdultEpi@meta.data,10)

VlnPlot(AdultEpi[],features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultEpi,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultEpi,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultEpisub<-subset(AdultEpi,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultEpisub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultEpisub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultEpisub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultEpisub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #4 OXT

AdultEpisub<-NormalizeData(AdultEpisub)
AdultEpisub <- FindVariableFeatures(AdultEpisub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultEpisub), 10)
variablegenes<-VariableFeatures(AdultEpisub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultEpisub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultEpisub)
AdultEpisub <- ScaleData(AdultEpisub, features = all.genes)
AdultEpisub <- RunPCA(AdultEpisub, features = VariableFeatures(object = AdultEpisub))
print(AdultEpisub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultEpisub, reduction ="pca")
ElbowPlot(AdultEpisub)  #use 15
AdultEpisub <- FindNeighbors(AdultEpisub, dims = 1:15)
AdultEpisub <- FindClusters(AdultEpisub, resolution = 0.8)
Maximum modularity in 10 random starts: 0.8045
Number of communities: 11
Elapsed time: 0 seconds

AdultEpisub <- RunUMAP(AdultEpisub, dims = 1:15)
umapplotout<-DimPlot(AdultEpisub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultIlesub, features = c("OXT"))

saveRDS(AdultEpisub, file = "AdultEpisubumap.rds")

write.table(dotout$data,file="Iletotalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Epi.markers <- FindAllMarkers(AdultEpisub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Epi.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultEpiclustermarkers.txt",sep="\t",row.names=FALSE)

Epivillus<-subset(AdultEpisub,subset=seurat_clusters=="1" |seurat_clusters=="6")

Epivilluscells<-Cells(Epivillus)

write.table(Epivilluscells,file="Epivilluscells.txt",sep="\n",row.names = FALSE)

#analyze just the villus
Epivillus<-NormalizeData(Epivillus)
Epivillus <- FindVariableFeatures(Epivillus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Epivillus)
Epivillus <- ScaleData(Epivillus, features = all.genesvillus)
Epivillus <- RunPCA(Epivillus, features = VariableFeatures(object = Epivillus))
print(Epivillus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Epivillus, reduction = "pca")
ElbowPlot(Epivillus)  #use 10
Epivillus <- FindNeighbors(Epivillus, dims = 1:10)
Epivillus <- FindClusters(Epivillus, resolution = 0.6)
Maximum modularity in 10 random starts: 0.6995
Number of communities: 5
Elapsed time: 0 seconds

OXTdata<-Epivillus@assays$RNA@data["OXT",]
OXTcount<-Epivillus@assays$RNA@counts["OXT",]

hist(OXTcount)
hist(OXTdata)
sum(OXTcount)#50 OXT
length(which(OXTcount==0)) #0 with OXT, 967 without cells

Epivillusnormdata<-as.data.frame(Epivillus@assays$RNA@data)
Epivilluscounts<-as.data.frame(Epivillus@assays$RNA@counts)

write.table(Epivillusnormdata,file="Epivillusnormdata.txt",sep="\t")
write.table(Epivilluscounts,file="Epivillusnormcounts.txt",sep="\t")

Epivillus <- RunUMAP(Epivillus, dims = 1:10)
Epiumapplotoutvillus<-DimPlot(Epivillus, reduction = "umap",label=TRUE)

pdf(file="Epiillusumap.pdf",width=5,height=4)
Epiumapplotoutvillus
dev.off()

pdf(file="EpiOXTumap.pdf",width=5,height=4)
FeaturePlot(Epivillus, features = c("OXT"))
dev.off()

saveRDS(Epivillus, file = "AdultEpivillusumap.rds")



Epivillus<-readRDS( file = "AdultEpivillusumap.rds")

Epivillus<-SCTransform(Epivillus, method = "glmGamPoi", vars.to.regress = "percent.mt",return.only.var.genes = FALSE,min_cells=1)
OXTdata<-Epivillus@assays$SCT@data["OXT",]
OXTcounts<-Epivillus@assays$SCT@counts["OXT",]

write.table(OXTdata, file="SCTdata/EPI_OXT_SCTdata.txt",sep="\t")
write.table(OXTcounts, file="SCTcounts/EPI_OXT_SCTcounts.txt",sep="\t")


#Ascending colon
Ascolon.data<-read.table(file = "AdultAscendingColon1.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultAscol<-CreateSeuratObject(counts=Ascolon.data,project = "AdultAscol",min.cells = 3, min.features = 200)
AdultAscol
#An object of class Seurat 
#11610 features across 1459 samples within 1 assay 
#Active assay: RNA (11610 features, 0 variable features)

AdultAscol[["percent.mt"]] <- PercentageFeatureSet(AdultAscol, pattern = "^MT-") 
head(AdultAscol@meta.data,10)

VlnPlot(AdultAscol,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultAscol,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultAscol,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultAscolsub<-subset(AdultAscol,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultAscolsub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultAscolsub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultAscolsub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultAscolsub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #0 OXT

AdultAscolsub<-NormalizeData(AdultAscolsub)
AdultAscolsub <- FindVariableFeatures(AdultAscolsub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultAscolsub), 10)
variablegenes<-VariableFeatures(AdultAscolsub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultAscolsub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultAscolsub)
AdultAscolsub <- ScaleData(AdultAscolsub, features = all.genes)
AdultAscolsub <- RunPCA(AdultAscolsub, features = VariableFeatures(object = AdultAscolsub))
print(AdultAscolsub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultAscolsub, reduction ="pca")
ElbowPlot(AdultAscolsub)  #use 20
AdultAscolsub <- FindNeighbors(AdultAscolsub, dims = 1:20)
AdultAscolsub <- FindClusters(AdultAscolsub, resolution = 0.8)
Maximum modularity in 10 random starts: 0.7242
Number of communities: 10
Elapsed time: 0 seconds

AdultAscolsub <- RunUMAP(AdultAscolsub, dims = 1:20)
umapplotout<-DimPlot(AdultAscolsub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultIlesub, features = c("OXT"))


saveRDS(AdultAscolsub, file = "AdultAscolsubumap.rds")

write.table(dotout$data,file="Ascoltotalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Ascol.markers <- FindAllMarkers(AdultAscolsub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Ascol.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultAscolclustermarkers.txt",sep="\t",row.names=FALSE)

Ascolvillus<-subset(AdultAscolsub,subset=seurat_clusters=="3" |seurat_clusters=="5" )

Ascolvilluscells<-Cells(Ascolvillus)
write.table(Ascolvilluscells,file="Ascolvilluscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Ascolvillus<-NormalizeData(Ascolvillus)
Ascolvillus <- FindVariableFeatures(Ascolvillus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Ascolvillus)
Ascolvillus <- ScaleData(Ascolvillus, features = all.genesvillus)
Ascolvillus <- RunPCA(Ascolvillus, features = VariableFeatures(object = Ascolvillus))
print(Ascolvillus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Ascolvillus, reduction = "pca")
ElbowPlot(Ascolvillus)  #use 20
Ascolvillus <- FindNeighbors(Ascolvillus, dims = 1:20)
Ascolvillus <- FindClusters(Ascolvillus, resolution = 0.4)
Maximum modularity in 10 random starts: 0.6913
Number of communities: 3
Elapsed time: 0 seconds

Ascolvillusnormdata<-as.data.frame(Ascolvillus@assays$RNA@data)
Ascolvilluscounts<-as.data.frame(Ascolvillus@assays$RNA@counts)

write.table(Ascolvillusnormdata,file="Ascolvillusnormdata.txt",sep="\t")
write.table(Ascolvilluscounts,file="Ascolvillusnormcounts.txt",sep="\t")

Ascolvillus <- RunUMAP(Ascolvillus, dims = 1:20)
Ascolumapplotoutvillus<-DimPlot(Ascolvillus, reduction = "umap",label=TRUE)

Ascolvilluscells
length(Ascolvilluscells) #232 total, 0 with OXT


pdf(file="Ascolvillusumap.pdf",width=5,height=4)
Ascolumapplotoutvillus
dev.off()


saveRDS(Ascolvillus, file = "AdultAscolvillusumap.rds")

Ascolvillus<-readRDS( file = "AdultAscolvillusumap.rds")

Ascolvillus<-SCTransform(Ascolvillus, method = "glmGamPoi", vars.to.regress = "percent.mt",return.only.var.genes = FALSE,min_cells=1)
OXTdata<-Ascolvillus@assays$SCT@data["OXT",]
OXTcounts<-Ascolvillus@assays$SCT@counts["OXT",]

write.table(OXTdata, file="SCTdata/ACL_OXT_SCTdata.txt",sep="\t")
write.table(OXTcounts, file="SCTcounts/ACL_OXT_SCTcounts.txt",sep="\t")


#Transverse Colon
TC1.data<-read.table(file = "AdultTransverseColon1.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultTC1<-CreateSeuratObject(counts=TC1.data,project = "AdultTransColon1",min.cells = 3, min.features = 200)
AdultTC1
#An object of class Seurat 
#17244 features across 4427 samples within 1 assay 
#Active assay: RNA (17244 features, 0 variable features)

AdultTC1[["percent.mt"]] <- PercentageFeatureSet(AdultTC1, pattern = "^MT-") 
head(AdultTC1@meta.data,10)

VlnPlot(AdultTC1,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultTC1,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultTC1,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultTC1sub<-subset(AdultTC1,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultTC1sub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultTC1sub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultTC1sub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultTC1sub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #0 OXT

AdultTC1sub<-NormalizeData(AdultTC1sub)
AdultTC1sub <- FindVariableFeatures(AdultTC1sub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultTC1sub), 10)
variablegenes<-VariableFeatures(AdultTC1sub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultTC1sub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultTC1sub)
AdultTC1sub <- ScaleData(AdultTC1sub, features = all.genes)
AdultTC1sub <- RunPCA(AdultTC1sub, features = VariableFeatures(object = AdultTC1sub))
print(AdultTC1sub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultTC1sub, reduction ="pca")
ElbowPlot(AdultTC1sub)  #use 20
AdultTC1sub <- FindNeighbors(AdultTC1sub, dims = 1:20)
AdultTC1sub <- FindClusters(AdultTC1sub, resolution = .7)
Maximum modularity in 10 random starts: 0.8615
Number of communities: 15
Elapsed time: 0 seconds

AdultTC1sub <- RunUMAP(AdultTC1sub, dims = 1:20)
umapplotout<-DimPlot(AdultTC1sub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultTC1sub, features = c("OXT"))


saveRDS(AdultTC1sub, file = "AdultTC1subumap.rds")

write.table(dotout$data,file="TC1totalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

TC1.markers <- FindAllMarkers(AdultTC1sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TC1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultTC1clustermarkers.txt",sep="\t",row.names=FALSE)

TC1villus<-subset(AdultTC1sub,subset=seurat_clusters=="0" |seurat_clusters=="1" |seurat_clusters=="2" |
                    seurat_clusters=="3" |seurat_clusters=="5"|seurat_clusters=="6"|seurat_clusters=="7"|
                    seurat_clusters=="11"|seurat_clusters=="13")

TC1villuscells<-Cells(TC1villus)

write.table(TC1villuscells,file="TC1villuscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
TC1villus<-NormalizeData(TC1villus)
TC1villus <- FindVariableFeatures(TC1villus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(TC1villus)
TC1villus <- ScaleData(TC1villus, features = all.genesvillus)
TC1villus <- RunPCA(TC1villus, features = VariableFeatures(object = TC1villus))
print(TC1villus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(TC1villus, reduction = "pca")
ElbowPlot(TC1villus)  #use 15
TC1villus <- FindNeighbors(TC1villus, dims = 1:15)
TC1villus <- FindClusters(TC1villus, resolution = 0.4)
Maximum modularity in 10 random starts: 0.8469
Number of communities: 6
Elapsed time: 0 seconds

OXTdata<-TC1villus@assays$RNA@data["OXT",]
OXTcount<-TC1villus@assays$RNA@counts["OXT",]

hist(OXTcount)
hist(OXTdata)
sum(OXTcount)#50 OXT
length(TC1villuscells) #3192 cells; none with OXT 

TC1villusnormdata<-as.data.frame(TC1villus@assays$RNA@data)
TC1villuscounts<-as.data.frame(TC1villus@assays$RNA@counts)

write.table(TC1villusnormdata,file="TC1villusnormdata.txt",sep="\t")
write.table(TC1villuscounts,file="TC1villusnormcounts.txt",sep="\t")

TC1villus <- RunUMAP(TC1villus, dims = 1:15)
TC1umapplotoutvillus<-DimPlot(TC1villus, reduction = "umap",label=TRUE)

pdf(file="TC1villusumap.pdf",width=5,height=4)
TC1umapplotoutvillus
dev.off()

pdf(file="TC1OXTumap.pdf",width=5,height=4)
FeaturePlot(TC1villus, features = c("OXT"))
dev.off()

saveRDS(TC1villus, file = "AdultTC1villusumap.rds")


#Transverse Colon2
TC2.data<-read.table(file = "AdultTransverseColon2.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultTC2<-CreateSeuratObject(counts=TC2.data,project = "AdultTransColon1",min.cells = 3, min.features = 200)
AdultTC2
An object of class Seurat 
17866 features across 8787 samples within 1 assay 
Active assay: RNA (17866 features, 0 variable features)

AdultTC2[["percent.mt"]] <- PercentageFeatureSet(AdultTC2, pattern = "^MT-") 
head(AdultTC2@meta.data,10)

VlnPlot(AdultTC2,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultTC2,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultTC2,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultTC2sub<-subset(AdultTC2,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultTC2sub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultTC2sub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultTC2sub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultTC2sub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #0 OXT

AdultTC2sub<-NormalizeData(AdultTC2sub)
AdultTC2sub <- FindVariableFeatures(AdultTC2sub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultTC2sub), 10)
variablegenes<-VariableFeatures(AdultTC2sub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultTC2sub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultTC2sub)
AdultTC2sub <- ScaleData(AdultTC2sub, features = all.genes)
AdultTC2sub <- RunPCA(AdultTC2sub, features = VariableFeatures(object = AdultTC2sub))
print(AdultTC2sub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultTC2sub, reduction ="pca")
ElbowPlot(AdultTC2sub)  #use 20
AdultTC2sub <- FindNeighbors(AdultTC2sub, dims = 1:20)
AdultTC2sub <- FindClusters(AdultTC2sub, resolution = .8)
Maximum modularity in 10 random starts: 0.8615
Number of communities: 15
Elapsed time: 0 seconds

AdultTC2sub <- RunUMAP(AdultTC2sub, dims = 1:20)
umapplotout<-DimPlot(AdultTC2sub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultTC2sub, features = c("OXT"))


saveRDS(AdultTC2sub, file = "AdultTC2subumap.rds")

write.table(dotout$data,file="TC2totalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

TC2.markers <- FindAllMarkers(AdultTC2sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TC2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultTC2clustermarkers.txt",sep="\t",row.names=FALSE)

TC2villus<-subset(AdultTC2sub,subset=seurat_clusters=="2" |seurat_clusters=="3" |seurat_clusters=="4" |
                    seurat_clusters=="5" |seurat_clusters=="6"|seurat_clusters=="7"|
                    seurat_clusters=="9"|seurat_clusters=="11")

TC2villuscells<-Cells(TC2villus)

write.table(TC2villuscells,file="TC2villuscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
TC2villus<-NormalizeData(TC2villus)
TC2villus <- FindVariableFeatures(TC2villus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(TC2villus)
TC2villus <- ScaleData(TC2villus, features = all.genesvillus)
TC2villus <- RunPCA(TC2villus, features = VariableFeatures(object = TC2villus))
print(TC2villus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(TC2villus, reduction = "pca")
ElbowPlot(TC2villus)  #use 15
TC2villus <- FindNeighbors(TC2villus, dims = 1:15)
TC2villus <- FindClusters(TC2villus, resolution = 0.4)
Maximum modularity in 10 random starts: 0.8844
Number of communities: 10
Elapsed time: 0 seconds


length(TC2villuscells) #4698 cells; none with OXT 

TC2villusnormdata<-as.data.frame(TC2villus@assays$RNA@data)
TC2villuscounts<-as.data.frame(TC2villus@assays$RNA@counts)

write.table(TC2villusnormdata,file="TC2villusnormdata.txt",sep="\t")
write.table(TC2villuscounts,file="TC2villusnormcounts.txt",sep="\t")

TC2villus <- RunUMAP(TC2villus, dims = 1:15)
TC2umapplotoutvillus<-DimPlot(TC2villus, reduction = "umap",label=TRUE)

pdf(file="TC2villusumap.pdf",width=5,height=4)
TC2umapplotoutvillus
dev.off()

pdf(file="TC2OXTumap.pdf",width=5,height=4)
FeaturePlot(TC2villus, features = c("OXT"))
dev.off()

saveRDS(TC2villus, file = "AdultTC2villusumap.rds")


#Sig Colon
Sig.data<-read.table(file = "AdultSigmoidColon1.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultSig<-CreateSeuratObject(counts=Sig.data,project = "AdultSig",min.cells = 3, min.features = 200)
AdultSig
#An object of class Seurat 
#14511 features across 2870 samples within 1 assay 
#Active assay: RNA (14511 features, 0 variable features)

AdultSig[["percent.mt"]] <- PercentageFeatureSet(AdultSig, pattern = "^MT-") 
head(AdultSig@meta.data,10)

VlnPlot(AdultSig,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSig,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSig,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultSigsub<-subset(AdultSig,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultSigsub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSigsub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSigsub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultSigsub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #0 OXT

AdultSigsub<-NormalizeData(AdultSigsub)
AdultSigsub <- FindVariableFeatures(AdultSigsub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultSigsub), 10)
variablegenes<-VariableFeatures(AdultSigsub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultSigsub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultSigsub)
AdultSigsub <- ScaleData(AdultSigsub, features = all.genes)
AdultSigsub <- RunPCA(AdultSigsub, features = VariableFeatures(object = AdultSigsub))
print(AdultSigsub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultSigsub, reduction ="pca")
ElbowPlot(AdultSigsub)  #use 20
AdultSigsub <- FindNeighbors(AdultSigsub, dims = 1:20)
AdultSigsub <- FindClusters(AdultSigsub, resolution = .9)
Maximum modularity in 10 random starts: 0.8140
Number of communities: 14
Elapsed time: 0 seconds

AdultSigsub <- RunUMAP(AdultSigsub, dims = 1:20)
umapplotout<-DimPlot(AdultSigsub, reduction = "umap",label=TRUE)

umapplotout  #total umap


saveRDS(AdultSigsub, file = "AdultSigsubumap.rds")

write.table(dotout$data,file="Sigtotalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Sig.markers <- FindAllMarkers(AdultSigsub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultSigclustermarkers.txt",sep="\t",row.names=FALSE)

Sigvillus<-subset(AdultSigsub,subset=seurat_clusters=="0" |seurat_clusters=="2" |seurat_clusters=="3" |
                    seurat_clusters=="5" |seurat_clusters=="6")

Sigvilluscells<-Cells(Sigvillus)

write.table(Sigvilluscells,file="Sigvilluscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Sigvillus<-NormalizeData(Sigvillus)
Sigvillus <- FindVariableFeatures(Sigvillus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Sigvillus)
Sigvillus <- ScaleData(Sigvillus, features = all.genesvillus)
Sigvillus <- RunPCA(Sigvillus, features = VariableFeatures(object = Sigvillus))
print(Sigvillus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Sigvillus, reduction = "pca")
ElbowPlot(Sigvillus)  #use 15
Sigvillus <- FindNeighbors(Sigvillus, dims = 1:15)
Sigvillus <- FindClusters(Sigvillus, resolution = 0.6)
Maximum modularity in 10 random starts: 0.7003
Number of communities: 6
Elapsed time: 0 seconds


length(Sigvilluscells) #1724 cells; none with OXT 

Sigvillusnormdata<-as.data.frame(Sigvillus@assays$RNA@data)
Sigvilluscounts<-as.data.frame(Sigvillus@assays$RNA@counts)

write.table(Sigvillusnormdata,file="Sigvillusnormdata.txt",sep="\t")
write.table(Sigvilluscounts,file="Sigvillusnormcounts.txt",sep="\t")

Sigvillus <- RunUMAP(Sigvillus, dims = 1:15)
Sigumapplotoutvillus<-DimPlot(Sigvillus, reduction = "umap",label=TRUE)

pdf(file="Sigvillusumap.pdf",width=5,height=4)
Sigumapplotoutvillus
dev.off()

pdf(file="TC2OXTumap.pdf",width=5,height=4)
FeaturePlot(TC2villus, features = c("OXT"))
dev.off()

saveRDS(Sigvillus, file = "AdultSigvillusumap.rds")


#Sigmodial colon
Sig.data<-read.table(file = "AdultSigmoidColon1.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultSig<-CreateSeuratObject(counts=Sig.data,project = "AdultSig",min.cells = 3, min.features = 200)
AdultSig
#An object of class Seurat 
#14511 features across 2870 samples within 1 assay 
#Active assay: RNA (14511 features, 0 variable features)

AdultSig[["percent.mt"]] <- PercentageFeatureSet(AdultSig, pattern = "^MT-") 
head(AdultSig@meta.data,10)

VlnPlot(AdultSig,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSig,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSig,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultSigsub<-subset(AdultSig,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultSigsub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSigsub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSigsub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultSigsub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #0 OXT

AdultSigsub<-NormalizeData(AdultSigsub)
AdultSigsub <- FindVariableFeatures(AdultSigsub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultSigsub), 10)
variablegenes<-VariableFeatures(AdultSigsub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultSigsub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultSigsub)
AdultSigsub <- ScaleData(AdultSigsub, features = all.genes)
AdultSigsub <- RunPCA(AdultSigsub, features = VariableFeatures(object = AdultSigsub))
print(AdultSigsub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultSigsub, reduction ="pca")
ElbowPlot(AdultSigsub)  #use 20
AdultSigsub <- FindNeighbors(AdultSigsub, dims = 1:20)
AdultSigsub <- FindClusters(AdultSigsub, resolution = .9)
Maximum modularity in 10 random starts: 0.8140
Number of communities: 14
Elapsed time: 0 seconds

AdultSigsub <- RunUMAP(AdultSigsub, dims = 1:20)
umapplotout<-DimPlot(AdultSigsub, reduction = "umap",label=TRUE)

umapplotout  #total umap


saveRDS(AdultSigsub, file = "AdultSigsubumap.rds")

write.table(dotout$data,file="Sigtotalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Sig.markers <- FindAllMarkers(AdultSigsub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Sig.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultSigclustermarkers.txt",sep="\t",row.names=FALSE)

Sigvillus<-subset(AdultSigsub,subset=seurat_clusters=="0" |seurat_clusters=="2" |seurat_clusters=="3" |
                    seurat_clusters=="5" |seurat_clusters=="6")

Sigvilluscells<-Cells(Sigvillus)

write.table(Sigvilluscells,file="Sigvilluscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Sigvillus<-NormalizeData(Sigvillus)
Sigvillus <- FindVariableFeatures(Sigvillus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Sigvillus)
Sigvillus <- ScaleData(Sigvillus, features = all.genesvillus)
Sigvillus <- RunPCA(Sigvillus, features = VariableFeatures(object = Sigvillus))
print(Sigvillus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Sigvillus, reduction = "pca")
ElbowPlot(Sigvillus)  #use 15
Sigvillus <- FindNeighbors(Sigvillus, dims = 1:15)
Sigvillus <- FindClusters(Sigvillus, resolution = 0.6)
Maximum modularity in 10 random starts: 0.7003
Number of communities: 6
Elapsed time: 0 seconds


length(Sigvilluscells) #1724 cells; none with OXT 

Sigvillusnormdata<-as.data.frame(Sigvillus@assays$RNA@data)
Sigvilluscounts<-as.data.frame(Sigvillus@assays$RNA@counts)

write.table(Sigvillusnormdata,file="Sigvillusnormdata.txt",sep="\t")
write.table(Sigvilluscounts,file="Sigvillusnormcounts.txt",sep="\t")

Sigvillus <- RunUMAP(Sigvillus, dims = 1:15)
Sigumapplotoutvillus<-DimPlot(Sigvillus, reduction = "umap",label=TRUE)

pdf(file="Sigvillusumap.pdf",width=5,height=4)
Sigumapplotoutvillus
dev.off()

pdf(file="TC2OXTumap.pdf",width=5,height=4)
FeaturePlot(TC2villus, features = c("OXT"))
dev.off()

saveRDS(Sigvillus, file = "AdultSigvillusumap.rds")


#Recturm
Rec.data<-read.table(file = "AdultRectum1.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultRec<-CreateSeuratObject(counts=Rec.data,project = "AdultRec",min.cells = 3, min.features = 200)
AdultRec
#An object of class Seurat 
#14511 features across 2870 samples within 1 assay 
#Active assay: RNA (14511 features, 0 variable features)

AdultRec[["percent.mt"]] <- PercentageFeatureSet(AdultRec, pattern = "^MT-") 
head(AdultRec@meta.data,10)

VlnPlot(AdultRec,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultRec,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultRec,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultRecsub<-subset(AdultRec,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultRecsub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultRecsub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultRecsub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultRecsub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #0 OXT

AdultRecsub<-NormalizeData(AdultRecsub)
AdultRecsub <- FindVariableFeatures(AdultRecsub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultRecsub), 10)
variablegenes<-VariableFeatures(AdultRecsub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultRecsub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultRecsub)
AdultRecsub <- ScaleData(AdultRecsub, features = all.genes)
AdultRecsub <- RunPCA(AdultRecsub, features = VariableFeatures(object = AdultRecsub))
print(AdultRecsub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultRecsub, reduction ="pca")
ElbowPlot(AdultRecsub)  #use 20
AdultRecsub <- FindNeighbors(AdultRecsub, dims = 1:20)
AdultRecsub <- FindClusters(AdultRecsub, resolution = .8)
Maximum modularity in 10 random starts: 0.7967
Number of communities: 14
Elapsed time: 0 seconds

AdultRecsub <- RunUMAP(AdultRecsub, dims = 1:20)
umapplotout<-DimPlot(AdultRecsub, reduction = "umap",label=TRUE)

umapplotout  #total umap


saveRDS(AdultRecsub, file = "AdultRecsubumap.rds")

write.table(dotout$data,file="Rectotalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Rec.markers <- FindAllMarkers(AdultRecsub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Rec.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultRecclustermarkers.txt",sep="\t",row.names=FALSE)

Recvillus<-subset(AdultRecsub,subset=seurat_clusters=="0" |seurat_clusters=="1" |seurat_clusters=="3" |
                    seurat_clusters=="5" |seurat_clusters=="6"|seurat_clusters=="8"|seurat_clusters=="11")

Recvilluscells<-Cells(Recvillus)

write.table(Recvilluscells,file="Recvilluscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Recvillus<-NormalizeData(Recvillus)
Redvillus <- FindVariableFeatures(Recvillus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Recvillus)
Recvillus <- ScaleData(Recvillus, features = all.genesvillus)
Recvillus <- RunPCA(Recvillus, features = VariableFeatures(object = Recvillus))
print(Recvillus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Recvillus, reduction = "pca")
ElbowPlot(Recvillus)  #use 15
Recvillus <- FindNeighbors(Recvillus, dims = 1:15)
Recvillus <- FindClusters(Recvillus, resolution = 0.7)
Maximum modularity in 10 random starts: 0.7081
Number of communities: 7
Elapsed time: 0 seconds


length(Recvilluscells) #3464 cells; none with OXT 

Recvillusnormdata<-as.data.frame(Recvillus@assays$RNA@data)
Recvilluscounts<-as.data.frame(Recvillus@assays$RNA@counts)

write.table(Recvillusnormdata,file="Recvillusnormdata.txt",sep="\t")
write.table(Recvilluscounts,file="Recvillusnormcounts.txt",sep="\t")

Recvillus <- RunUMAP(Recvillus, dims = 1:15)
Recumapplotoutvillus<-DimPlot(Recvillus, reduction = "umap",label=TRUE)

pdf(file="Recvillusumap.pdf",width=5,height=4)
Recumapplotoutvillus
dev.off()

#pdf(file="TC2OXTumap.pdf",width=5,height=4)
#FeaturePlot(TC2villus, features = c("OXT"))
#dev.off()

saveRDS(Recvillus, file = "AdultRecvillusumap.rds")


#stomach1
Sto1.data<-read.table(file = "AdultStomach1.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultSto1<-CreateSeuratObject(counts=Sto1.data,project = "AdultSto1",min.cells = 3, min.features = 200)
AdultSto1
#An object of class Seurat 
#13560 features across 1809 samples within 1 assay 
#Active assay: RNA (13560 features, 0 variable features

AdultSto1[["percent.mt"]] <- PercentageFeatureSet(AdultSto1, pattern = "^MT-") 
head(AdultSto1@meta.data,10)

VlnPlot(AdultSto1,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSto1,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSto1,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultSto1sub<-subset(AdultSto1,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultSto1sub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSto1sub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSto1sub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultSto1sub@assays$RNA@data["OXT",]
hist(OXTcount)
sum(OXTcount) #0 OXT

AdultSto1sub<-NormalizeData(AdultSto1sub)
AdultSto1sub <- FindVariableFeatures(AdultSto1sub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultSto1sub), 10)
variablegenes<-VariableFeatures(AdultSto1sub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultSto1sub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultSto1sub)
AdultSto1sub <- ScaleData(AdultSto1sub, features = all.genes)
AdultSto1sub <- RunPCA(AdultSto1sub, features = VariableFeatures(object = AdultSto1sub))
print(AdultSto1sub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultSto1sub, reduction ="pca")
ElbowPlot(AdultSto1sub)  #use 20
AdultSto1sub <- FindNeighbors(AdultSto1sub, dims = 1:20)
AdultSto1sub <- FindClusters(AdultSto1sub, resolution = .8)
Maximum modularity in 10 random starts: 0.7742
Number of communities: 10
Elapsed time: 0 seconds

AdultSto1sub <- RunUMAP(AdultSto1sub, dims = 1:20)
umapplotout<-DimPlot(AdultSto1sub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultSto1sub, features = c("OXT"))


saveRDS(AdultSto1sub, file = "AdultSto1subumap.rds")

write.table(dotout$data,file="Sto1totalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Sto1.markers <- FindAllMarkers(AdultSto1sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Sto1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultSto1clustermarkers.txt",sep="\t",row.names=FALSE)

Sto1villus<-subset(AdultSto1sub,subset=seurat_clusters=="0" |seurat_clusters=="1" |seurat_clusters=="2" |
                    seurat_clusters=="3" |seurat_clusters=="6"|seurat_clusters=="7"|seurat_clusters=="8")

Sto1villuscells<-Cells(Sto1villus)

write.table(Sto1villuscells,file="Sto1villuscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Sto1villus<-NormalizeData(Sto1villus)
Sto1villus <- FindVariableFeatures(Sto1villus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Sto1villus)
Sto1villus <- ScaleData(Sto1villus, features = all.genesvillus)
Sto1villus <- RunPCA(Sto1villus, features = VariableFeatures(object = Sto1villus))
print(Sto1villus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Sto1villus, reduction = "pca")
ElbowPlot(Sto1villus)  #use 15
Sto1villus <- FindNeighbors(Sto1villus, dims = 1:15)
Sto1villus <- FindClusters(Sto1villus, resolution = 0.8)
Maximum modularity in 10 random starts: 0.7152
Number of communities: 8
Elapsed time: 0 seconds

OXTdata<-Ilevillus@assays$RNA@data["OXT",]
OXTcount<-Ilevillus@assays$RNA@counts["OXT",]

hist(OXTcount)
hist(OXTdata)
sum(OXTcount)
length(Sto1villuscells) #1498 without cells

Sto1villusnormdata<-as.data.frame(Sto1villus@assays$RNA@data)
Sto1villuscounts<-as.data.frame(Sto1villus@assays$RNA@counts)

write.table(Sto1villusnormdata,file="Sto1villusnormdata.txt",sep="\t")
write.table(Sto1villuscounts,file="Sto1villusnormcounts.txt",sep="\t")

Sto1villus <- RunUMAP(Sto1villus, dims = 1:10)
Sto1umapplotoutvillus<-DimPlot(Sto1villus, reduction = "umap",label=TRUE)

pdf(file="Sto1villusumap.pdf",width=5,height=4)
Sto1umapplotoutvillus
dev.off()

pdf(file="IleOXTumap.pdf",width=5,height=4)
FeaturePlot(Ilevillus, features = c("OXT"))
dev.off()

saveRDS(Sto1villus, file = "AdultSto1villusumap.rds")


#stomach 2

Sto2.data<-read.table(file = "AdultStomach2.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultSto2<-CreateSeuratObject(counts=Sto2.data,project = "AdultSto2",min.cells = 3, min.features = 200)
AdultSto2
#An object of class Seurat 
#14986 features across 4648 samples within 1 assay 
#Active assay: RNA (14986 features, 0 variable features)

AdultSto2[["percent.mt"]] <- PercentageFeatureSet(AdultSto2, pattern = "^MT-") 
head(AdultSto2@meta.data,10)

VlnPlot(AdultSto2,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSto2,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSto2,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultSto2sub<-subset(AdultSto2,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultSto2sub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSto2sub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSto2sub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultSto2sub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #7 OXT

AdultSto2sub<-NormalizeData(AdultSto2sub)
AdultSto2sub <- FindVariableFeatures(AdultSto2sub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultSto2sub), 10)
variablegenes<-VariableFeatures(AdultSto2sub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultSto2sub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultSto2sub)
AdultSto2sub <- ScaleData(AdultSto2sub, features = all.genes)
AdultSto2sub <- RunPCA(AdultSto2sub, features = VariableFeatures(object = AdultSto2sub))
print(AdultSto2sub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultSto2sub, reduction ="pca")
ElbowPlot(AdultSto2sub)  #use 20
AdultSto2sub <- FindNeighbors(AdultSto2sub, dims = 1:20)
AdultSto2sub <- FindClusters(AdultSto2sub, resolution = .9)
Maximum modularity in 10 random starts: 0.7840
Number of communities: 15
Elapsed time: 0 seconds

AdultSto2sub <- RunUMAP(AdultSto2sub, dims = 1:20)
umapplotout<-DimPlot(AdultSto2sub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultSto2sub, features = c("OXT"))


saveRDS(AdultSto2sub, file = "AdultSto2subumap.rds")

write.table(dotout$data,file="Sto2totalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Sto2.markers <- FindAllMarkers(AdultSto2sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Sto2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultSto2clustermarkers.txt",sep="\t",row.names=FALSE)

Sto2villus<-subset(AdultSto2sub,subset=seurat_clusters=="1" |seurat_clusters=="3" |
                    seurat_clusters=="12")


Sto2villuscells<-Cells(Sto2villus)

write.table(Sto2villuscells,file="Sto2villuscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Sto2villus<-NormalizeData(Sto2villus)
Sto2villus <- FindVariableFeatures(Sto2villus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Sto2villus)
Sto2villus <- ScaleData(Sto2villus, features = all.genesvillus)
Sto2villus <- RunPCA(Sto2villus, features = VariableFeatures(object = Sto2villus))
print(Sto2villus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Sto2villus, reduction = "pca")
ElbowPlot(Sto2villus)  #use 10
Sto2villus <- FindNeighbors(Sto2villus, dims = 1:10)
Sto2villus <- FindClusters(Sto2villus, resolution = 0.5)
Maximum modularity in 10 random starts: 0.6995
Number of communities: 4
Elapsed time: 0 seconds

OXTdata<-Sto2villus@assays$RNA@data["OXT",]
OXTcount<-Sto2villus@assays$RNA@counts["OXT",]

hist(OXTcount)
hist(OXTdata)
sum(OXTcount)#3 OXT
length(which(OXTcount==0)) #3 with OXT, 1281 total cells

Sto2villusnormdata<-as.data.frame(Sto2villus@assays$RNA@data)
Sto2villuscounts<-as.data.frame(Sto2villus@assays$RNA@counts)

write.table(Sto2villusnormdata,file="Sto2villusnormdata.txt",sep="\t")
write.table(Sto2villuscounts,file="Sto2villusnormcounts.txt",sep="\t")

Sto2villus <- RunUMAP(Sto2villus, dims = 1:10)
Sto2umapplotoutvillus<-DimPlot(Sto2villus, reduction = "umap",label=TRUE)

pdf(file="Sto2villusumap.pdf",width=5,height=4)
Sto2umapplotoutvillus
dev.off()

pdf(file="Sto2OXTumap.pdf",width=5,height=4)
FeaturePlot(Sto2villus, features = c("OXT"))
dev.off()

saveRDS(Sto2villus, file = "AdultSto2villusumap.rds")


Sto2villus<-readRDS( file = "AdultSto2villusumap.rds")

Sto2villus<-SCTransform(Sto2villus, method = "glmGamPoi", vars.to.regress = "percent.mt",return.only.var.genes = FALSE,min_cells=1)
OXTdata<-Sto2villus@assays$SCT@data["OXT",]
OXTcounts<-Sto2villus@assays$SCT@counts["OXT",]

write.table(OXTdata, file="SCTdata/PYL_OXT_SCTdata.txt",sep="\t")
write.table(OXTcounts, file="SCTcounts/PYL_OXT_SCTcounts.txt",sep="\t")


#Stomach3

Sto3.data<-read.table(file = "AdultStomach3.rmbatchdge.txt", header=TRUE,row.names=1,sep = " ", as.is=T,check.names = FALSE)
AdultSto3<-CreateSeuratObject(counts=Sto3.data,project = "AdultSto3",min.cells = 3, min.features = 200)
AdultSto3
#An object of class Seurat 
#14986 features across 4648 samples within 1 assay 
#Active assay: RNA (14986 features, 0 variable features)

AdultSto3[["percent.mt"]] <- PercentageFeatureSet(AdultSto3, pattern = "^MT-") 
head(AdultSto3@meta.data,10)

VlnPlot(AdultSto3,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSto3,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSto3,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

AdultSto3sub<-subset(AdultSto3,subset =   nCount_RNA> 50&nFeature_RNA >70 &nFeature_RNA <1800 &percent.mt <45)
VlnPlot(AdultSto3sub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(AdultSto3sub,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(AdultSto3sub,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

OXTcount<-AdultSto3sub@assays$RNA@counts["OXT",]
hist(OXTcount)
sum(OXTcount) #0 OXT

AdultSto3sub<-NormalizeData(AdultSto3sub)
AdultSto3sub <- FindVariableFeatures(AdultSto3sub, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AdultSto3sub), 10)
variablegenes<-VariableFeatures(AdultSto3sub)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(AdultSto3sub )
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(AdultSto3sub)
AdultSto3sub <- ScaleData(AdultSto3sub, features = all.genes)
AdultSto3sub <- RunPCA(AdultSto3sub, features = VariableFeatures(object = AdultSto3sub))
print(AdultSto3sub[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(AdultSto3sub, reduction ="pca")
ElbowPlot(AdultSto3sub)  #use 20
AdultSto3sub <- FindNeighbors(AdultSto3sub, dims = 1:20)
AdultSto3sub <- FindClusters(AdultSto3sub, resolution = .7)
Maximum modularity in 10 random starts: 0.8642
Number of communities: 14
Elapsed time: 0 seconds

AdultSto3sub <- RunUMAP(AdultSto3sub, dims = 1:20)
umapplotout<-DimPlot(AdultSto3sub, reduction = "umap",label=TRUE)

umapplotout  #total umap
FeaturePlot(AdultSto3sub, features = c("OXT"))


saveRDS(AdultSto3sub, file = "AdultSto3subumap.rds")

write.table(dotout$data,file="Sto3totalclustersdotplotdata.txt",sep="\t",row.names = FALSE)

Sto3.markers <- FindAllMarkers(AdultSto3sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Sto3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% write.table(file="AdultSto3clustermarkers.txt",sep="\t",row.names=FALSE)

Sto3villus<-subset(AdultSto3sub,subset=seurat_clusters=="0" |seurat_clusters=="1" |seurat_clusters=="2"|
                     seurat_clusters=="3"|seurat_clusters=="4"|seurat_clusters=="5"|seurat_clusters=="11"|
                     seurat_clusters=="12")


Sto3villuscells<-Cells(Sto3villus)

write.table(Sto3villuscells,file="Sto3villuscells.txt",sep="\n",row.names = FALSE)
#analyze just the villus
Sto3villus<-NormalizeData(Sto3villus)
Sto3villus <- FindVariableFeatures(Sto3villus, selection.method = "vst", nfeatures = 5000)
all.genesvillus <- rownames(Sto3villus)
Sto3villus <- ScaleData(Sto3villus, features = all.genesvillus)
Sto3villus <- RunPCA(Sto3villus, features = VariableFeatures(object = Sto3villus))
print(Sto3villus[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(Sto3villus, reduction = "pca")
ElbowPlot(Sto3villus)  #use 15
Sto3villus <- FindNeighbors(Sto3villus, dims = 1:15)
Sto3villus <- FindClusters(Sto3villus, resolution = 0.6)
Maximum modularity in 10 random starts: 0.8459
Number of communities: 9
Elapsed time: 0 seconds

OXTdata<-Sto3villus@assays$RNA@data["OXT",]
OXTcount<-Sto3villus@assays$RNA@counts["OXT",]

hist(OXTcount)
hist(OXTdata)
sum(OXTcount)
length(Sto3villuscells) #6161 total

Sto3villusnormdata<-as.data.frame(Sto3villus@assays$RNA@data)
Sto3villuscounts<-as.data.frame(Sto3villus@assays$RNA@counts)

write.table(Sto3villusnormdata,file="Sto3villusnormdata.txt",sep="\t")
write.table(Sto3villuscounts,file="Sto3villusnormcounts.txt",sep="\t")

Sto3villus <- RunUMAP(Sto3villus, dims = 1:15)
Sto3umapplotoutvillus<-DimPlot(Sto3villus, reduction = "umap",label=TRUE)

pdf(file="Sto3villusumap.pdf",width=5,height=4)
Sto3umapplotoutvillus
dev.off()

pdf(file="Sto3OXTumap.pdf",width=5,height=4)
FeaturePlot(Sto3villus, features = c("OXT"))
dev.off()

saveRDS(Sto3villus, file = "AdultSto3villusumap.rds")

