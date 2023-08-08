library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(gridExtra)
library(tidyverse)
library(SeuratDisk)
library(future)
library(glmGamPoi)
options(future.globals.maxSize = 1000 * 1024^2) #or more

#script adapted from 
#https://satijalab.org/seurat/articles/integration_introduction.html
#https://satijalab.org/seurat/articles/sctransform_vignette.html
#https://satijalab.org/seurat/archive/v3.0/future_vignette.html

#data from https://www.gutcellatlas.org/

#load data from Gut Cell Atlas
#can skip to line 167 by loading the rds file for figure 3c and d.
#the beginning part of the script describes where the data for fig 1a came from
Convert("epi_raw_counts02.h5ad", dest = "h5seurat", overwrite = TRUE)
epiraw <- LoadH5Seurat("epi_raw_counts02.h5seurat")
names(epiraw@meta.data)

#subset
Adultonly<-subset(epiraw,subset =Age_group=="Adult")  #only adult samples
saveRDS(Adultonly,file="Adult_epi_raw_counts02.rds")

range(levels(factor(epiraw@meta.data$doublet_scores)))
#quality control and doublet removal
VlnPlot(Adultonly,features=c("n_genes", "nFeature_RNA","nCount_RNA","pct_counts_mt","doublet_scores"),ncol=4)
#remove doublet score of greater than 0.25 and less than 50% mitoreads and at least 500 genes with at least 3 cells
#data already look filtered for mitoreads, min number of genes; 
Adultonlysub<-subset(Adultonly,subset =   doublet_scores<0.25  &pct_counts_mt <50 & category=="Epithelial" &Age_group=="Adult")
Adultonlysub@meta.data<-droplevels(Adultonlysub@meta.data)

saveRDS(Adultonlysub,file="subAdult_epi_raw_counts02.rds")


#remove low batches
#first look at the number of cells in each batch
batchcount.list <- lapply(X = c(1:74), FUN = function(x) {
  x <- subset(Adultonlysub,subset=batch==levels(Adultonlysub@meta.data$batch)[x])
  x <- length(Cells(x))
})

unlist(batchcount.list)

#remove those with batch<50   #this means the batch has less than 50 cells
batchestorm<-which(unlist(batchcount.list)<50)
levels(Adultonlysub@meta.data$batch)[batchestorm]

Adultonlysubnolowbatch<-subset(Adultonlysub,subset= 
batch!= "Human_colon_16S8001907"& batch!="Human_colon_16S8002566" & batch!= "Human_colon_16S8002573"& batch!="Human_colon_16S8002625" &batch!="Human_colon_16S8123913" &
 batch!= "Human_colon_16S8123919" & batch!= "Human_colon_16S8159191" & batch!= "Human_colon_16S8159193" & batch!= "WTDAtest7770718" &batch!="WTDAtest7844022" &
 batch!="WTDAtest7844026" )

Adultonlysubnolowbatch@meta.data<-droplevels(Adultonlysubnolowbatch@meta.data)

#check that removal of low cell number batches was done

 batchcountcheck.list <- lapply(X = c(1:63), FUN = function(x) {
  x <- subset(Adultonlysubnolowbatch,subset=batch==levels(Adultonlysubnolowbatch@meta.data$batch)[x])
  x <- length(Cells(x))
})

saveRDS(Adultonlysubnolowbatch,file="subAdult_epi_batchgreater50_counts02.rds")

#split seurat object into batches
batched.list <- SplitObject(Adultonlysubnolowbatch, split.by = "batch")
saveRDS(batched.list,"subAdult_epi_batchgreater50_counts02_batchsplit.rds")

#batched.list<-readRDS(file="subAdult_epi_batchgreater50_counts02_batchsplit.rds") #for reloading and starting here


# needed, split up the batches here and put data back together later

#SCTransform
batched.list <- lapply(X = batched.list, FUN = function(x) {	
	x<-SCTransform(x, method = "glmGamPoi", vars.to.regress = "pct_counts_mt",return.only.var.genes = FALSE,min_cells=1)
})


#split up the batches, merge everthing back together for example:
mergedsctadult1<-readRDS(file="subAdult_epi_batchgreater50_counts02_mergedbatchsplit_1to15_SCT.rds")
mergedsctadult2<-readRDS(file="subAdult_epi_batchgreater50_counts02_mergedbatchsplit_16to30_SCT.rds")
mergedsctadult10<-merge(x=mergedsctadult1,y=mergedsctadult2)
#repeat as needed

#then split up by Region
region.list <- SplitObject(mergedsctadult10, split.by = "Region.code")

#next get out the SCTransformed data and counts 
sctdata.list <- lapply(X = region.list, FUN = function(x) {	
	x<-x@assays$SCT@data })
sctcounts.list <- lapply(X = region.list, FUN = function(x) {	
	x<-x@assays$SCT@counts })	

#pull out the OXT data
sctdata.OXTlist <- lapply(X = sctdata.list, FUN = function(x) {	
  x<-as.data.frame(x["OXT",]) })

#export data per each region, only one example shown here for counts used in Fig 1A
ACLOXTdata<-sctdata.OXTlist$ACL
write.table(ACLOXTdata,file="ACL_OXT_SCTcounts.txt",sep ="\t")


###for the jejunum umap plots in Fig3C and D

#get out the jejunum
Adultjej<-subset(Adultonlysub, subset = Region.code =="JEJ")
Adultjej@meta.data<-droplevels(Adultjej@meta.data)
levels(Adultjej@meta.data$batch)
length(levels(Adultjej@meta.data$batch))


#remove batches with low number of cells
batch.list <- lapply(X = c(1:5), FUN = function(x) {
  x <- subset(Adultjej,subset=batch==levels(Adultjej@meta.data$batch)[x])
  x <- length(Cells(x))
})

unlist(batch.list)
#[1]   99  459  499 1734   26
#just remove the last one
Adultjejbatchrm<-subset(Adultjej, subset = batch !="WTDAtest7844022")
Adultjejbatchrm@meta.data<-droplevels(Adultjejbatchrm@meta.data)
levels(Adultjejbatchrm@meta.data$batch)
saveRDS(Adultjejbatchrm,file="Adultjejbatchremoved.rds")

#split by batch
batched.list <- SplitObject(Adultjejbatchrm, split.by = "batch")

batched.list <- lapply(X = batched.list, FUN = function(x) {	
  x<-SCTransform(x, method = "glmGamPoi", vars.to.regress = "pct_counts_mt",return.only.var.genes = FALSE,min_cells=1)
})
features <- SelectIntegrationFeatures(object.list = batched.list, nfeatures = 3000)
saveRDS(features,file="jejonlySCTfeatures.rds")
batched.list <- PrepSCTIntegration(object.list = batched.list, anchor.features = features)

#get anchors and do batch integration
batched.anchors <- FindIntegrationAnchors(object.list = batched.list, normalization.method = "SCT",
    anchor.features = features)
jej.combined.sct <- IntegrateData(anchorset = batched.anchors, normalization.method = "SCT",k.weight =90)
saveRDS(jej.combined.sct,file="jejbatchedonly_integrateddata.rds")
jej.combined.sct <- RunPCA(jej.combined.sct, verbose = FALSE)

#features<-jej.combined.sct@assays$integrated@var.features  #if reloading and need to get features back
jej.combined.sct <- RunUMAP(jej.combined.sct, reduction = "pca", dims = 1:30)

ElbowPlot(jej.combined.sct)
jej.combined.sct <- FindNeighbors(jej.combined.sct, reduction = "pca", dims = 1:20)

#if not set properly
DefaultAssay(jej.combined.sct)<-"integrated"
jej.combined.sct@active.assay

jej.combined.sct <- FindClusters(jej.combined.sct, resolution = .8)
# Maximum modularity in 10 random starts: 0.7564
# Number of communities: 11
# Elapsed time: 0 seconds

saveRDS(jej.combined.sct,file="Fig3cd_Adult_jej_nolowbatch_sctintegration.rds")

#if restarting
jej.combined.sct<-readRDS(file="Fig3cd_Adult_jej_nolowbatch_sctintegration.rds")
Idents(jej.combined.sct)<-"seurat_clusters"  #if not set properly
DefaultAssay(jej.combined.sct)<-"integrated"  #this has the umap data if need to switch back to umap

umapplotoutsct<-DimPlot(jej.combined.sct, reduction = "umap",label=TRUE)



pdf(file="adultjejumap.pdf",width=5,height=4)
umapplotoutsct
dev.off()

DefaultAssay(jej.combined.sct)<-"SCT"
featureplotout<-FeaturePlot(jej.combined.sct, features = c("OXT")) #,"NTS","CHGA"))
pdf(file="adultjejfeatureplot_OXT.pdf",width=5,height=4)
featureplotout
dev.off()


Jej.markers <- FindAllMarkers(jej.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Jej.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) %>% write.table(file="fromjejsubintegration_AdultJejclustermarkersSCTmethod.txt",sep="\t",row.names=FALSE)


vlnclusplot<-VlnPlot(jej.combined.sct, features = c("APOA4","SLC46A1","MTRNR2L12","NR4A1","FABP1","RBP2",
                                       "APOC3","APOA1","RPS2","RPS8","BEST4","GUCA2A","REG1A",
                                       "GPX2","JUN","HSPA1B","OLFM4","TUBA1B","MUC2","CHGA",
                                       "BMX","SH2D6", "OXT"), fill="idents", ncol=6)

pdf(file="violinplotclusterOXT.pdf",height=16,width=24)
vlnclusplot
dev.off()



#to try to get out the OXT+ positive cells markers
#pull out the cells with OXT+ 
OXTcellsnum<-which(jej.combined.sct@assays$SCT@counts["OXT",]!=0)
OXTcells<-names(jej.combined.sct@assays$SCT@counts["OXT",OXTcellsnum])
#pull out all other cells
Negcellsnum<-which(jej.combined.sct@assays$SCT@counts["OXT",]==0)
Negcells<-names(jej.combined.sct@assays$SCT@counts["OXT",Negcellsnum])

OXTid<-c()
for (i in c(1:length(jej.combined.sct@assays$SCT@counts["OXT",]) )) {
  if (jej.combined.sct@assays$SCT@counts["OXT",i] == 0) {
    OXTstage<-"OXTNeg"
  }
  else if (jej.combined.sct@assays$SCT@counts["OXT",i] != 0) {
    OXTstage<-"OXTPos"
  }
  OXTid<-c(OXTid,OXTstage)
  
}
OXTid[332]

#name in metadata
jej.combined.sct[["OXTpos"]] <- OXTid
Idents(jej.combined.sct)<-"OXTpos"

OXTmarkers<-FindMarkers(object=jej.combined.sct, ident.1 = "OXTPos",ident.2="OXTNeg",verbose = FALSE)
write.table(OXTmarkers,file="OXTmarkerscomplete.txt",sep="\t")
