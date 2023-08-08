library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(gridExtra)

#how to load in data provided by Beumer et al 2020
### Can skip the beginning part using this rds file: FigS3ad_EECsubumap.rds, see line 99
EEC.data <- read.table(file = "GSE146799_EEC_atlas_raw.csv", header = TRUE, row.names=1, sep=",", as.is=T,check.names = FALSE)
EEC<-CreateSeuratObject(counts=EEC.data,project = "EEC",min.cells = 3, min.features = 200)
EEC
# An object of class Seurat 
# 20490 features across 6875 samples within 1 assay 
# Active assay: RNA (20490 features, 0 variable features)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
EEC[["percent.mt"]] <- PercentageFeatureSet(EEC, pattern = "^MT-") #I don't think this is actually mt genes
head(EEC@meta.data,10)

VlnPlot(EEC,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(EEC,feature1="nCount_RNA",feature2="percent.mt")
plot2<-FeatureScatter(EEC,feature1="nCount_RNA",feature2="nFeature_RNA")
plot1+plot2

#to add metadata
cellnames<-row.names(EEC@meta.data)
splitcellnames<-matrix(unlist(str_split(cellnames, "_",n=6)), ncol = 6,byrow = TRUE)
#unique(splitcellnames[,6])
setwithmore<-which(str_detect(splitcellnames[,6],"_")==TRUE)
setwithless<-which(str_detect(splitcellnames[,6],"_")==FALSE)

splitcellnamesless<-as.data.frame(matrix(unlist(str_split(cellnames[setwithless], "_",n=6)), ncol = 6,byrow = TRUE))
splitcellnamesmore<-as.data.frame(matrix(unlist(str_split(cellnames[setwithmore], "_",n=7)), ncol = 7,byrow = TRUE))
names(splitcellnamesless)=c("Segment","NGN3","BMP","Tph1","Plate","Well")
names(splitcellnamesmore)=c("Segment","NGN3","BMP","Tph1","Plate","Marker","Well")

splitcellnamesless$Marker="NA"
row.names(splitcellnamesless) = cellnames[setwithless]
row.names(splitcellnamesmore)=cellnames[setwithmore]

splitcells<-rbind(splitcellnamesless,splitcellnamesmore)
splitcells$allinfo<-paste(splitcells$Segment,splitcells$NGN3,splitcells$BMP,splitcells$Tph1,splitcells$Marker,sep="_")
EEC <- AddMetaData(object = EEC, metadata = splitcells)


EECsub<-subset(EEC,subset =  nFeature_RNA >1100 &nFeature_RNA <9000)
EEC$orig.ident
table(EECsub$orig.ident)
#Col  Duo  Ile   #total if not filtering
#1169 1576 2894 

#after filtering
#Col  Duo  Ile 
#693 1461 2229 
table(EECsub$NGN3)
#Ngn3- Ngn3? Ngn3+ 
#  297   420  3666 
table(EECsub$BMP)
#BMP- BMP+ 
#  2509 1874 
table(EECsub$Tph1) #fixed the thp1 typo in file
#Thp1+ Tph1- Tph1? Tph1+ 
 # 2  1671  2574   136 
table(EECsub$Marker)
#gcg- gcg+ mln- mln+   NA sst- sst+ 
 # 54  188   37  148 3823   49   84 

VlnPlot(EECsub,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
#EECsub_SI<-subset(EECsub,subset=orig.ident!="Col")
#table(EECsub_SI$orig.ident)

EECsub<-NormalizeData(EECsub)
EECsub <- FindVariableFeatures(EECsub, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(EECsub), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(EECsub)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(EECsub)
EECsub <- ScaleData(EECsub, features = all.genes)

EECsub <- RunPCA(EECsub, features = VariableFeatures(object = EECsub))
print(EECsub[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(EECsub, reduction = "pca")
ElbowPlot(EECsub)  #use 17

EECsub <- FindNeighbors(EECsub, dims = 1:17)
EECsub <- FindClusters(EECsub, resolution = .9)
#Maximum modularity in 10 random starts: 0.9146
#Number of communities: 22
#Elapsed time: 0 seconds

EECsub <- RunUMAP(EECsub, dims = 1:17)

saveRDS(EECsub, file = "FigS3ad_EECsubumap.rds")
EECsub<-readRDS(file = "FigS3ad_EECsubumap.rds")

pdf(file="Beumerumap.pdf",width=5,height=4)
DimPlot(EECsub, reduction = "umap",label=TRUE)
dev.off()

DimPlot(EECsub, reduction = "umap",group.by = "Marker")

pdf(file="Beumerumapsegment.pdf",width=5,height=4)
DimPlot(EECsub, reduction = "umap",group.by = "Segment")
dev.off()

pdf(file="BeumerumapOXT.pdf",width=5,height=4)
FeaturePlot(EECsub, features = c("OXT--chr20")) #"OXT--chr20"
dev.off()





# find markers for every cluster compared to all remaining cells, report only the positive ones
EEC.markers <- FindAllMarkers(EECsub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EEC.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) %>% write.table(file="EECclustermarkers22clusters.txt",sep="\t",row.names=FALSE)



vlnclusplot<-VlnPlot(EECsub, features = c("NEUROD1--chr2","TUBA1A--chr12","CHGA--chr14","TPH1--chr11",
                                          "GAST--chr17","GIP--chr17","GCG--chr2","NTS--chr12","PYY--chr17",
                                          "MLN--chr6","GHRL--chr3","SST--chr3","FABP2--chr4","FABP6--chr5",
                             "OAT--chr10","MTTP--chr4","RBP2--chr3","SI--chr3","APOA1--chr11","APOB--chr2",
                             "MUC1--chr1","MUC2--chr11","MUC6--chr11","LYZ--chr12",
                             "OLFM4--chr13","REG1A--chr2","MKI67--chr10",
                             "OXT--chr20"))



pdf(file="violinplotclusterOXT.pdf",height=16,width=24)
vlnclusplot
dev.off()


