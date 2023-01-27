library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(PMCMRplus)


#Need to pick a number of cells to rarify to: use 1281
#For those with OXT+ cells, randomly select that number of cells
#load in lists of cells with data
villustotalfiles<-c("FigS1a_PYL_OXT_SCTdata.txt","FigS1a_DUO_OXT_SCTdata.txt","FigS1a_JEJ_OXT_SCTdata.txt",
                    "FigS1a_ILE_OXT_SCTdata.txt")
PYLdata<-read.delim(villustotalfiles[1],header=TRUE,check.names = FALSE)
DUOdata<-read.delim(villustotalfiles[2],header=TRUE,check.names = FALSE)
JEJdata<-read.delim(villustotalfiles[3],header=TRUE,check.names = FALSE)
ILEdata<-read.delim(villustotalfiles[4],header=TRUE,check.names = FALSE)


alldata<-list(PYLdata,DUOdata,JEJdata,ILEdata)

regionnames<-c("Pylorus","Duodenum", "Jejunum", "Ileum")
corpus<-data.frame("Region"=c(rep("Corpus",1281)),"OXTnorm"=c(rep (0,1281) ) )
gastricbody<-data.frame("Region"=c(rep("Gastric Body",1281)),"OXTnorm"=c(rep (0,1281) ) )
tracol1<-data.frame("Region"=c(rep("Mid Transverse Colon",1281)),"OXTnorm"=c(rep (0,1281) ) )
tracol2<-data.frame("Region"=c(rep("Lower Transverse Colon",1281)),"OXTnorm"=c(rep (0,1281) ) )
sigmo<-data.frame("Region"=c(rep("Sigmodial Colon",1281)),"OXTnorm"=c(rep (0,1281) ) )
rec<-data.frame("Region"=c(rep("Rectum",1281)),"OXTnorm"=c(rep (0,1281) ) )


#rare at 1281

#run with loop:
kwtestnormpvals<-c()
kwtestcountpvals<-c()
jejnormdunnvals<-data.frame("Corpus"=c(),"Gastric Body"=c(),"Pylorus"=c(),"Duodenum"=c(),"Ileum"=c(),"Mid Transverse Colon"=c(),
                            "Lower Transverse Colon"=c(), "Sigmodial Colon"=c(),"Rectum"=c())
duonormdunnvals<-data.frame("Corpus"=c(),"Gastric Body"=c(),"Pylorus"=c(),"Jejunum"=c(),"Ileum"=c(),"Mid Transverse Colon"=c(),
                            "Lower Transverse Colon"=c(), "Sigmodial Colon"=c(),"Rectum"=c())
pylnormdunnvals<-data.frame("Corpus"=c(),"Gastric Body"=c(),"Duodenum"=c(),"Jejunum"=c(),"Ileum"=c(),"Mid Transverse Colon"=c(),
                            "Lower Transverse Colon"=c(), "Sigmodial Colon"=c(),"Rectum"=c())
ilenormdunnvals<-data.frame("Corpus"=c(),"Gastric Body"=c(),"Pylorus"=c(),"Duodenum"=c(),"Jejunum"=c(),"Mid Transverse Colon"=c(),
                            "Lower Transverse Colon"=c(), "Sigmodial Colon"=c(),"Rectum"=c())

jejcountdunnvals<-jejnormdunnvals
duocountdunnvals<-duonormdunnvals
pylcountdunnvals<-pylnormdunnvals
ilecountdunnvals<-ilenormdunnvals

iterationdata<-c()



for (i in c(1:10000)) {
  totaldata<-data.frame("Region"=c(),"OXTnorm"=c())
  totalpoints<-data.frame("Region"=c(),"OXTpositive"=c(),"Iteration"=c())
  totalcount<-data.frame("Region"=c(), "OXTcount"=c())
  
for (filecount in c(1:length(villustotalfiles))) {
  villusdata<-alldata[[filecount]]
  cellnums<-sample(1:length(villusdata[,2]), 1281, replace=FALSE) #pick out 1712 cells
  villusdata<-villusdata[cellnums,2] #pick out 1281 cells
  OXTdata<-unlist(villusdata[])
  OXTpoints<-which(villusdata[]!=0)
  OXTpos<-length(OXTpoints)
  outdata<-data.frame("Region"=c(rep(regionnames[filecount],1281)), "OXTnorm"=OXTdata)
  if (OXTpos>0) {
    outpoints<-data.frame("Region"=c(rep(regionnames[filecount],OXTpos)), "OXTpositive"=OXTdata[OXTpoints],"Iteration"=i)
    totalpoints<-rbind(totalpoints,outpoints)
  }
  totaldata<-rbind(totaldata,outdata)
  countframe<-data.frame("Region"=regionnames[filecount], "OXTcount"=OXTpos)
  totalcount<-rbind(totalcount,countframe)
    }
  #put all regions together
  totaldata2<-rbind(totaldata,corpus,gastricbody,tracol1,tracol2,sigmo,rec)
  totaldata2$Region<-factor(totaldata2$Region,c("Corpus","Gastric Body","Pylorus","Duodenum","Jejunum", "Ileum","Mid Transverse Colon", "Lower Transverse Colon",
                                                "Sigmodial Colon","Rectum")) #fix the order
  iterationdata<-rbind(iterationdata,totalpoints)#for looking over all the iterations and getting the mean data
  
  #run kruskal-wallis
  kwtestnorm<-kruskal.test(OXTnorm~Region,data=totaldata2)  #norm data
  kwtestnormp<-kwtestnorm$p.value  #norm data
  kwtestnormpvals<-c(kwtestnormpvals,kwtestnormp)  #norm data
  binary<-totaldata2
  binary$OXTnorm[as.numeric(binary$OXTnorm) > 0] <- 1
  kwtestcount<-kruskal.test(OXTnorm~Region,data=binary) #count data
  kwtestcountp<-kwtestcount$p.value #count data
  kwtestcountpvals<-c(kwtestcountpvals,kwtestcountp) #count data
  #run Dunn's test
  #norm data
  dunnout<-kwAllPairsDunnTest(OXTnorm~Region,data=totaldata2,p.adjust.method = "BH")
  #pylorus
  pylvscorp<-dunnout$p.value[[2,1]]
  pylvsgb<-dunnout$p.value[[2,2]]
  pylvd<-dunnout$p.value[[3,3]]
  pylvj<-dunnout$p.value[[4,3]]
  pylvi<-dunnout$p.value[[5,3]]
  pylvmtc<-dunnout$p.value[[6,3]]
  pylvltc<-dunnout$p.value[[7,3]]
  pylsc<-dunnout$p.value[[8,3]]
  pylr<-dunnout$p.value[[9,3]]
  normdunnvals<-data.frame(pylvscorp,pylvsgb,pylvd,pylvj,pylvi,pylvmtc,pylvltc,pylsc,pylr)
  pylnormdunnvals<-rbind(pylnormdunnvals,normdunnvals)
  #jejunum
  jejvscorp<-dunnout$p.value[[4,1]]
  jejvsgb<-dunnout$p.value[[4,2]]
  jejvp<-dunnout$p.value[[4,3]]
  jejvd<-dunnout$p.value[[4,4]]
  jejvi<-dunnout$p.value[[5,5]]
  jejvmtc<-dunnout$p.value[[6,5]]
  jejvltc<-dunnout$p.value[[7,5]]
  jejsc<-dunnout$p.value[[8,5]]
  jejr<-dunnout$p.value[[9,5]]
  normdunnvals<-data.frame(jejvscorp,jejvsgb,jejvp,jejvd,jejvi,jejvmtc,jejvltc,jejsc,jejr)
  jejnormdunnvals<-rbind(jejnormdunnvals,normdunnvals)
  #duodenum
  duovscorp<-dunnout$p.value[[3,1]]
  duovsgb<-dunnout$p.value[[3,2]]
  duovp<-dunnout$p.value[[3,3]]
  duovj<-dunnout$p.value[[4,4]]
  duovi<-dunnout$p.value[[5,4]]
  duovmtc<-dunnout$p.value[[6,4]]
  duovltc<-dunnout$p.value[[7,4]]
  duosc<-dunnout$p.value[[8,4]]
  duor<-dunnout$p.value[[9,4]]
  normdunnvals<-data.frame(duovscorp,duovsgb,duovp,duovj,duovi,duovmtc,duovltc,duosc,duor)
  duonormdunnvals<-rbind(duonormdunnvals,normdunnvals)
  #ileum
  ilevscorp<-dunnout$p.value[[5,1]]
  ilevsgb<-dunnout$p.value[[5,2]]
  ilevp<-dunnout$p.value[[5,3]]
  ilevd<-dunnout$p.value[[5,4]]
  ilevj<-dunnout$p.value[[5,5]]
  ilevmtc<-dunnout$p.value[[6,6]]
  ilevltc<-dunnout$p.value[[7,6]]
  ilesc<-dunnout$p.value[[8,6]]
  iler<-dunnout$p.value[[9,6]]
  normdunnvals<-data.frame(ilevscorp,ilevsgb,ilevp,ilevd,ilevj,ilevmtc,ilevltc,ilesc,iler)
  ilenormdunnvals<-rbind(ilenormdunnvals,normdunnvals)
  
  
  #count data
  dunnout<-kwAllPairsDunnTest(OXTnorm~Region,data=binary,p.adjust.method = "BH")
  #pylorus
  pylvscorp<-dunnout$p.value[[2,1]]
  pylvsgb<-dunnout$p.value[[2,2]]
  pylvd<-dunnout$p.value[[3,3]]
  pylvj<-dunnout$p.value[[4,3]]
  pylvi<-dunnout$p.value[[5,3]]
  pylvmtc<-dunnout$p.value[[6,3]]
  pylvltc<-dunnout$p.value[[7,3]]
  pylsc<-dunnout$p.value[[8,3]]
  pylr<-dunnout$p.value[[9,3]]
  countdunnvals<-data.frame(pylvscorp,pylvsgb,pylvd,pylvj,pylvi,pylvmtc,pylvltc,pylsc,pylr)
  pylcountdunnvals<-rbind(pylcountdunnvals,countdunnvals)
  
  #duodenum
  duovscorp<-dunnout$p.value[[3,1]]
  duovsgb<-dunnout$p.value[[3,2]]
  duovp<-dunnout$p.value[[3,3]]
  duovj<-dunnout$p.value[[4,4]]
  duovi<-dunnout$p.value[[5,4]]
  duovmtc<-dunnout$p.value[[6,4]]
  duovltc<-dunnout$p.value[[7,4]]
  duosc<-dunnout$p.value[[8,4]]
  duor<-dunnout$p.value[[9,4]]
  countdunnvals<-data.frame(duovscorp,duovsgb,duovp,duovj,duovi,duovmtc,duovltc,duosc,duor)
  duocountdunnvals<-rbind(duocountdunnvals,countdunnvals)
  #jejunum
  jejvscorp<-dunnout$p.value[[4,1]]
  jejvsgb<-dunnout$p.value[[4,2]]
  jejvp<-dunnout$p.value[[4,3]]
  jejvd<-dunnout$p.value[[4,4]]
  jejvi<-dunnout$p.value[[5,5]]
  jejvmtc<-dunnout$p.value[[6,5]]
  jejvltc<-dunnout$p.value[[7,5]]
  jejsc<-dunnout$p.value[[8,5]]
  jejr<-dunnout$p.value[[9,5]]
  countdunnvals<-data.frame(jejvscorp,jejvsgb,jejvp,jejvd,jejvi,jejvmtc,jejvltc,jejsc,jejr)
  jejcountdunnvals<-rbind(jejcountdunnvals,countdunnvals)
  #ileum
  ilevscorp<-dunnout$p.value[[5,1]]
  ilevsgb<-dunnout$p.value[[5,2]]
  ilevp<-dunnout$p.value[[5,3]]
  ilevd<-dunnout$p.value[[5,4]]
  ilevj<-dunnout$p.value[[5,5]]
  ilevmtc<-dunnout$p.value[[6,6]]
  ilevltc<-dunnout$p.value[[7,6]]
  ilesc<-dunnout$p.value[[8,6]]
  iler<-dunnout$p.value[[9,6]]
  countdunnvals<-data.frame(ilevscorp,ilevsgb,ilevp,ilevd,ilevj,ilevmtc,ilevltc,ilesc,iler)
  ilecountdunnvals<-rbind(ilecountdunnvals,countdunnvals)
}



kwtestnormpvals
pylnormdunnvals
duonormdunnvals
jejnormdunnvals
ilenormdunnvals

kwtestcountpvals
pylcountdunnvals
duocountdunnvals
jejcountdunnvals
ilecountdunnvals

write.table(iterationdata,"OXTpositivealliterations.txt",row.names = FALSE,sep="\t",col.names = FALSE)
#change values to counts
iterationdata$OXTcount=1
iterationsum<-aggregate(OXTcount~Region+Iteration,FUN=sum,data=iterationdata)  #sum per iteration per region
iterationmean<-aggregate(OXTcount~Region,FUN=mean,data=iterationsum)  #mean per region
write.table(iterationmean,"OXTpositivemean.txt",row.names=FALSE,sep="\t",col.names = FALSE)


write.table(kwtestnormpvals,"TotalKWtest_normpvalsSCT.txt",row.names = FALSE,sep="\t",col.names = FALSE)
write.table(kwtestcountpvals,"TotalKWtest_countpvalsSCT.txt",row.names = FALSE,sep="\t",col.names = FALSE)
write.table(jejnormdunnvals,"JejunumDunntest_normpvalsSCT.txt",row.names = FALSE,sep="\t")
write.table(jejcountdunnvals,"JejunumDunntest_countpvalsSCT.txt",row.names = FALSE,sep="\t")
write.table(duonormdunnvals,"DuoDunntest_normpvalsSCT.txt",row.names = FALSE,sep="\t")
write.table(duocountdunnvals,"DuoDunntest_countpvalsSCT.txt",row.names = FALSE,sep="\t")
write.table(pylnormdunnvals,"PylDunntest_normpvalsSCT.txt",row.names = FALSE,sep="\t")
write.table(pylcountdunnvals,"PylDunntest_countpvalsSCT.txt",row.names = FALSE,sep="\t")
write.table(ilenormdunnvals,"IleDunntest_normpvalsSCT.txt",row.names = FALSE,sep="\t")
write.table(ilecountdunnvals,"IleDunntest_countpvalsSCT.txt",row.names = FALSE,sep="\t")



reps=10000
which(kwtestnormpvals>0.05)  #max(kwtestnormpvals) 0.005823368
which(kwtestcountpvals>0.05) #max(kwtestcountpvals) 0.005833293
length(which(jejnormdunnvals$jejvscorp>0.05))/reps #2e-04 corpus
length(which(jejnormdunnvals$jejvsgb>0.05))/reps #2e-04 gastric body
length(which(jejnormdunnvals$jejvp>0.05))/reps #0.0071 pylorus
length(which(jejnormdunnvals$jejvd>0.05))/reps #0.0554 duodenum
length(which(jejnormdunnvals$jejvi>0.05))/reps #0.0068 ileum
length(which(jejnormdunnvals$jejvmtc>0.05))/reps #2e-04 mid transverse
length(which(jejnormdunnvals$jejvltc>0.05))/reps #2e-04 lower transverse
length(which(jejnormdunnvals$jejsc>0.05))/reps #2e-04 sigmodial
length(which(jejnormdunnvals$jejr>0.05))/reps #2e-04 rectum

length(which(jejcountdunnvals$jejvscorp>0.05))/reps #2e-04 corpus
length(which(jejcountdunnvals$jejvsgb>0.05))/reps #2e-04 gastric body
length(which(jejcountdunnvals$jejvp>0.05))/reps #0.0071 pylorus
length(which(jejcountdunnvals$jejvd>0.05))/reps #0.0554 duodenum
length(which(jejcountdunnvals$jejvi>0.05))/reps #0.0068 ileum
length(which(jejcountdunnvals$jejvmtc>0.05))/reps #2e-04 mid transverse
length(which(jejcountdunnvals$jejvltc>0.05))/reps #2e-04 lower transverse
length(which(jejcountdunnvals$jejsc>0.05))/reps #2e-04 sigmodial
length(which(jejcountdunnvals$jejr>0.05))/reps #2e-04 rectum

#duodenum
length(which(duonormdunnvals$duovscorp>0.05))/reps #0.8035 corpus
length(which(duonormdunnvals$duovsgb>0.05))/reps #0.8035 gastric body
length(which(duonormdunnvals$duovp>0.05))/reps #1 pylorus
length(which(duonormdunnvals$duovj>0.05))/reps #0.0554 jejunum
length(which(duonormdunnvals$duovi>0.05))/reps #0.9998 ileum
length(which(duonormdunnvals$duovmtc>0.05))/reps #0.8035 mid transverse
length(which(duonormdunnvals$duovltc>0.05))/reps #0.8035 lower transverse
length(which(duonormdunnvals$duosc>0.05))/reps #0.8035 sigmodial
length(which(duonormdunnvals$duor>0.05))/reps #0.8035 rectum

length(which(duocountdunnvals$duovscorp>0.05))/reps #0.8035 corpus
length(which(duocountdunnvals$duovsgb>0.05))/reps #0.8035 gastric body
length(which(duocountdunnvals$duovp>0.05))/reps #1 pylorus
length(which(duocountdunnvals$duovj>0.05))/reps #0.0554 jejunum
length(which(duocountdunnvals$duovi>0.05))/reps #0.9998 ileum
length(which(duocountdunnvals$duovmtc>0.05))/reps #0.8035 mid transverse
length(which(duocountdunnvals$duovltc>0.05))/reps #0.8035 lower transverse
length(which(duocountdunnvals$duosc>0.05))/reps #0.8035 sigmodial
length(which(duocountdunnvals$duor>0.05))/reps #0.8035 rectum

#pyl
length(which(pylnormdunnvals$pylvscorp>0.05))/reps #1 corpus
length(which(pylnormdunnvals$pylvsgb>0.05))/reps #1gastric body
length(which(pylnormdunnvals$pylvd>0.05))/reps #1 duodenum
length(which(pylnormdunnvals$pylvj>0.05))/reps #0.0071 jejunum
length(which(pylnormdunnvals$pylvi>0.05))/reps #1 ileum
length(which(pylnormdunnvals$pylvmtc>0.05))/reps #1 mid transverse
length(which(pylnormdunnvals$pylvltc>0.05))/reps #1 lower transverse
length(which(pylnormdunnvals$pylsc>0.05))/reps #1 sigmodial
length(which(pylnormdunnvals$pylr>0.05))/reps #1 rectum

length(which(pylcountdunnvals$pylvscorp>0.05))/reps #1 corpus
length(which(pylcountdunnvals$pylvsgb>0.05))/reps #1 gastric body
length(which(pylcountdunnvals$pylvd>0.05))/reps #1 duodenum
length(which(pylcountdunnvals$pylvj>0.05))/reps #0.0071 jejunum
length(which(pylcountdunnvals$pylvi>0.05))/reps #1 ileum
length(which(pylcountdunnvals$pylvmtc>0.05))/reps #1 mid transverse
length(which(pylcountdunnvals$pylvltc>0.05))/reps #1 lower transverse
length(which(pylcountdunnvals$pylsc>0.05))/reps #1 sigmodial
length(which(pylcountdunnvals$pylr>0.05))/reps #1 rectum
#ile
length(which(ilenormdunnvals$ilevscorp>0.05))/reps #1 corpus
length(which(ilenormdunnvals$ilevsgb>0.05))/reps #1 gastric body
length(which(ilenormdunnvals$ilevp>0.05))/reps #1 pylorus
length(which(ilenormdunnvals$ilevd>0.05))/reps #0.9998 duodenum
length(which(ilenormdunnvals$ilevj>0.05))/reps #0.0068 jejunum
length(which(ilenormdunnvals$ilevmtc>0.05))/reps #1 mid transverse
length(which(ilenormdunnvals$ilevltc>0.05))/reps #1 lower transverse
length(which(ilenormdunnvals$ilesc>0.05))/reps #1 sigmodial
length(which(ilenormdunnvals$iler>0.05))/reps #1 rectum

length(which(ilecountdunnvals$ilevscorp>0.05))/reps #1 corpus
length(which(ilecountdunnvals$ilevsgb>0.05))/reps #1 gastric body
length(which(ilecountdunnvals$ilevp>0.05))/reps #1 pylorus
length(which(ilecountdunnvals$ilevd>0.05))/reps #0.9998 duodenum
length(which(ilecountdunnvals$ilevj>0.05))/reps #0.0068 jejunum
length(which(ilecountdunnvals$ilevmtc>0.05))/reps #1 mid transverse
length(which(ilecountdunnvals$ilevltc>0.05))/reps #1 lower transverse
length(which(ilecountdunnvals$ilesc>0.05))/reps #1 sigmodial
length(which(ilecountdunnvals$iler>0.05))/reps #1 rectum

OXTplot<-ggplot(totaldata2 ,aes(x=Region,y=OXTnorm, fill=Region))+
  geom_violin()+
#  geom_point(data=totalpoints,  aes(x=Region,y=OXTpositive,fill=Region))
  geom_jitter(data=totalpoints,  aes(x=Region,y=OXTpositive,fill=Region),pch=21,width = 0.2,height=.05)+
  ylab("OXT normalized counts")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.5)))+
  scale_fill_npg()+
  theme_classic() + theme(
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=7, colour="black"),
    axis.title.x= element_text(size=7, colour="black"),
    axis.text.x= element_text(size = 6, colour="black",angle = 45,hjust = 1),
    axis.text.y= element_text(size=6, colour="black"),
    )#legend.text=element_text(size=16))

pdf(file="Handataviolinplot.pdf",width=3,height=2.2)
OXTplot
dev.off()

