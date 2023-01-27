library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(PMCMRplus)



#Need to pick a number of cells to rarify to: use 1712
#load in lists of cells with data
#villuscellfiles<-c("Sto2villuscells.txt","Duovilluscells.txt","Ilevilluscells.txt")
villustotalfiles<-c("Fig1a_TCL_OXT_SCTdata.txt", "Fig1a_SCL_OXT_SCTdata.txt", "Fig1a_REC_OXT_SCTdata.txt", 
                    "Fig1a_JEJ_OXT_SCTdata.txt","Fig1a_ILE_OXT_SCTdata.txt","Fig1a_DUO_OXT_SCTdata.txt",
                    "Fig1a_DCL_OXT_SCTdata.txt","Fig1a_CAE_OXT_SCTdata.txt","Fig1a_APD_OXT_SCTdata.txt",
                    "Fig1a_ACL_OXT_SCTdata.txt")
TCLdata<-read.delim(villustotalfiles[1],header=TRUE,check.names = FALSE)
SCLdata<-read.delim(villustotalfiles[2],header=TRUE,check.names = FALSE)
RECdata<-read.delim(villustotalfiles[3],header=TRUE,check.names = FALSE)
JEJdata<-read.delim(villustotalfiles[4],header=TRUE,check.names = FALSE)
ILEdata<-read.delim(villustotalfiles[5],header=TRUE,check.names = FALSE)
DUOdata<-read.delim(villustotalfiles[6],header=TRUE,check.names = FALSE)
DCLdata<-read.delim(villustotalfiles[7],header=TRUE,check.names = FALSE)
CAEdata<-read.delim(villustotalfiles[8],header=TRUE,check.names = FALSE)
APDdata<-read.delim(villustotalfiles[9],header=TRUE,check.names = FALSE)
ACLdata<-read.delim(villustotalfiles[10],header=TRUE,check.names = FALSE)


alldata<-list(DUOdata,JEJdata,ILEdata,CAEdata,APDdata,ACLdata,TCLdata,DCLdata,SCLdata,RECdata)

regionnames<-c("Duodenum", "Jejunum", "Ileum","Caecum","Appendix","Ascending Colon","Transverse Colon",
               "Descending Colon","Sigmoidal Colon", "Rectum")

#rare at 1712

#run with loop:
kwtestnormpvals<-c()
kwtestcountpvals<-c()
duonormdunnvals<-data.frame("Jejunum"=c(),"Ileum"=c(),"Caecum"=c(),
                            "Appendix"=c(), "Ascending Colon"=c(),"Transverse Colon"=c(),
                            "Descending Colon"=c(), "Sigmodial Colon"=c(),"Rectum"=c())
jejnormdunnvals<-data.frame("Duodenum"=c(),"Ileum"=c(),"Caecum"=c(),
                            "Appendix"=c(), "Ascending Colon"=c(),"Transverse Colon"=c(),
                            "Descending Colon"=c(), "Sigmodial Colon"=c(),"Rectum"=c())
ilenormdunnvals<-data.frame("Duodenum"=c(),"Jejunum"=c(),"Caecum"=c(),
                            "Appendix"=c(), "Ascending Colon"=c(),"Transverse Colon"=c(),
                            "Descending Colon"=c(), "Sigmodial Colon"=c(),"Rectum"=c())

duocountdunnvals<-duonormdunnvals
jejcountdunnvals<-jejnormdunnvals
ilecountdunnvals<-ilenormdunnvals

iterationdata<-c()


for (i in c(1:10000)) {
  totaldata<-data.frame("Region"=c(),"OXTnorm"=c())
  totalpoints<-data.frame("Region"=c(),"OXTpositive"=c())
  totalcount<-data.frame("Region"=c(), "OXTcount"=c())
  
for (filecount in c(1:length(villustotalfiles))) {
  villusdata<-alldata[[filecount]]
  cellnums<-sample(1:length(villusdata[,2]), 1712, replace=FALSE) #pick out 1712 cells
  villusdata<-villusdata[cellnums,2] #pick out 1281 cells
   OXTdata<-unlist(villusdata[])
  OXTpoints<-which(villusdata[]!=0)
  OXTpos<-length(OXTpoints)
  outdata<-data.frame("Region"=c(rep(regionnames[filecount],1712)), "OXTnorm"=OXTdata)
  if (OXTpos>0) {
    outpoints<-data.frame("Region"=c(rep(regionnames[filecount],OXTpos)), "OXTpositive"=OXTdata[OXTpoints],"Iteration"=i)
    totalpoints<-rbind(totalpoints,outpoints)
  }
  totaldata<-rbind(totaldata,outdata)
  totalpoints<-rbind(totalpoints,outpoints)
  countframe<-data.frame("Region"=regionnames[filecount], "OXTcount"=OXTpos)
  totalcount<-rbind(totalcount,countframe)
    }
  #put all regions together
  totaldata2<-rbind(totaldata)
  totaldata2$Region<-factor(totaldata2$Region,regionnames) #fix the order
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
  #jejunum
  jejvd<-dunnout$p.value[[1,1]]
  jejvi<-dunnout$p.value[[2,2]]
  jejvcae<-dunnout$p.value[[3,2]]
  jejvapp<-dunnout$p.value[[4,2]]
  jejvasc<-dunnout$p.value[[5,2]]
  jejvtc<-dunnout$p.value[[6,2]]
  jejvdc<-dunnout$p.value[[7,2]]
  jejvsc<-dunnout$p.value[[8,2]]
  jejvrect<-dunnout$p.value[[9,2]]
  normdunnvals<-data.frame(jejvd,jejvi,jejvcae,jejvapp,jejvasc,jejvtc,jejvdc,jejvsc,jejvrect)
  jejnormdunnvals<-rbind(jejnormdunnvals,normdunnvals)
  #duodenum
  duovjej<-dunnout$p.value[[1,1]]
  duovi<-dunnout$p.value[[2,1]]
  duovcae<-dunnout$p.value[[3,1]]
  duovapp<-dunnout$p.value[[4,1]]
  duovasc<-dunnout$p.value[[5,1]]
  duovtc<-dunnout$p.value[[6,1]]
  duovdc<-dunnout$p.value[[7,1]]
  duovsc<-dunnout$p.value[[8,1]]
  duovrect<-dunnout$p.value[[9,1]]
  normdunnvals<-data.frame(duovjej,duovi,duovcae,duovapp,duovasc,duovtc,duovdc,duovsc,duovrect)
  duonormdunnvals<-rbind(duonormdunnvals,normdunnvals)
  #ileum
  ilevd<-dunnout$p.value[[2,1]]
  ilevjej<-dunnout$p.value[[2,2]]
  ilevcae<-dunnout$p.value[[3,3]]
  ilevapp<-dunnout$p.value[[4,3]]
  ilevasc<-dunnout$p.value[[5,3]]
  ilevtc<-dunnout$p.value[[6,3]]
  ilevdc<-dunnout$p.value[[7,3]]
  ilevsc<-dunnout$p.value[[8,3]]
  ilevrect<-dunnout$p.value[[9,3]]
  normdunnvals<-data.frame(ilevd,ilevjej,ilevcae,ilevapp,ilevasc,ilevtc,ilevdc,ilevsc,ilevrect)
  ilenormdunnvals<-rbind(ilenormdunnvals,normdunnvals)
  
  #count data
  dunnout<-kwAllPairsDunnTest(OXTnorm~Region,data=binary,p.adjust.method = "BH")
  #jejunum
  jejvd<-dunnout$p.value[[1,1]]
  jejvi<-dunnout$p.value[[2,2]]
  jejvcae<-dunnout$p.value[[3,2]]
  jejvapp<-dunnout$p.value[[4,2]]
  jejvasc<-dunnout$p.value[[5,2]]
  jejvtc<-dunnout$p.value[[6,2]]
  jejvdc<-dunnout$p.value[[7,2]]
  jejvsc<-dunnout$p.value[[8,2]]
  jejvrect<-dunnout$p.value[[9,2]]
  countdunnvals<-data.frame(jejvd,jejvi,jejvcae,jejvapp,jejvasc,jejvtc,jejvdc,jejvsc,jejvrect)
  jejcountdunnvals<-rbind(jejcountdunnvals,countdunnvals)

  #duodenum
  duovjej<-dunnout$p.value[[1,1]]
  duovi<-dunnout$p.value[[2,1]]
  duovcae<-dunnout$p.value[[3,1]]
  duovapp<-dunnout$p.value[[4,1]]
  duovasc<-dunnout$p.value[[5,1]]
  duovtc<-dunnout$p.value[[6,1]]
  duovdc<-dunnout$p.value[[7,1]]
  duovsc<-dunnout$p.value[[8,1]]
  duovrect<-dunnout$p.value[[9,1]]
  countdunnvals<-data.frame(duovjej,duovi,duovcae,duovapp,duovasc,duovtc,duovdc,duovsc,duovrect)
  duocountdunnvals<-rbind(duocountdunnvals,countdunnvals)
  
  #ileum
  ilevd<-dunnout$p.value[[2,1]]
  ilevjej<-dunnout$p.value[[2,2]]
  ilevcae<-dunnout$p.value[[3,3]]
  ilevapp<-dunnout$p.value[[4,3]]
  ilevasc<-dunnout$p.value[[5,3]]
  ilevtc<-dunnout$p.value[[6,3]]
  ilevdc<-dunnout$p.value[[7,3]]
  ilevsc<-dunnout$p.value[[8,3]]
  ilevrect<-dunnout$p.value[[9,3]]
  countdunnvals<-data.frame(ilevd,ilevjej,ilevcae,ilevapp,ilevasc,ilevtc,ilevdc,ilevsc,ilevrect)
  ilecountdunnvals<-rbind(ilecountdunnvals,countdunnvals)
  
}

kwtestnormpvals
duonormdunnvals
jejnormdunnvals
ilenormdunnvals

kwtestcountpvals
duocountdunnvals
jejcountdunnvals
ilecountdunnvals

write.table(iterationdata,"OXTpositivealliterations.txt",row.names = FALSE,sep="\t",col.names = FALSE)
#change values to counts
iterationdata$OXTcount=1
iterationsum<-aggregate(OXTcount~Region+Iteration,FUN=sum,data=iterationdata)  #sum per iteration per region
iterationmean<-aggregate(OXTcount~Region,FUN=mean,data=iterationsum)  #mean per region
write.table(iterationmean,"OXTpositivemean.txt",row.names=FALSE,sep="\t",col.names = FALSE)


write.table(kwtestnormpvals,"TotalKWtest_normpvals.txt",row.names = FALSE,sep="\t",col.names = FALSE)
write.table(kwtestcountpvals,"TotalKWtest_countpvals.txt",row.names = FALSE,sep="\t",col.names = FALSE)
write.table(jejnormdunnvals,"JejunumDunntest_normpvals.txt",row.names = FALSE,sep="\t")
write.table(jejcountdunnvals,"JejunumDunntest_countpvals.txt",row.names = FALSE,sep="\t")
write.table(duonormdunnvals,"DuoDunntest_normpvals.txt",row.names = FALSE,sep="\t")
write.table(duocountdunnvals,"DuoDunntest_countpvals.txt",row.names = FALSE,sep="\t")
write.table(ilenormdunnvals,"IleDunntest_normpvals.txt",row.names = FALSE,sep="\t")
write.table(ilecountdunnvals,"IleDunntest_countpvals.txt",row.names = FALSE,sep="\t")

which(kwtestnormpvals>0.05)  #max(kwtestnormpvals) 0
which(kwtestcountpvals>0.05) #max(kwtestcountpvals) 0

length(which(jejnormdunnvals$jejvd>0.05))/10000 #duodenum 0
length(which(jejnormdunnvals$jejvi>0.05))/10000 #ileum 0
length(which(jejnormdunnvals$jejvcae>0.05))/10000 #caecum 0
length(which(jejnormdunnvals$jejvapp>0.05))/10000 #appendix 0
length(which(jejnormdunnvals$jejvasc>0.05))/10000 #appendix 0
length(which(jejnormdunnvals$jejvtc>0.05))/10000 #transverse colon 0
length(which(jejnormdunnvals$jejvdc>0.05))/10000 #descending colon 0
length(which(jejnormdunnvals$jejvsc>0.05))/10000 #sigmodial 0
length(which(jejnormdunnvals$jejvrect>0.05))/10000 #rectum 0

length(which(jejcountdunnvals$jejvd>0.05))/10000 #duodenum 0
length(which(jejcountdunnvals$jejvi>0.05))/10000 #ileum 0
length(which(jejcountdunnvals$jejvcae>0.05))/10000 #caecum 0
length(which(jejcountdunnvals$jejvapp>0.05))/10000 #appendix 0
length(which(jejcountdunnvals$jejvasc>0.05))/10000 #asc 0
length(which(jejcountdunnvals$jejvtc>0.05))/10000 #transverse colon 0 
length(which(jejcountdunnvals$jejvdc>0.05))/10000 #descending colon 0
length(which(jejcountdunnvals$jejvsc>0.05))/10000 #sigmodial 0
length(which(jejcountdunnvals$jejvrect>0.05))/10000 #rectum 0

length(which(duonormdunnvals$duovjej>0.05))/10000 #jejunum 0
length(which(duonormdunnvals$duovi>0.05))/10000 #ileum 0.0046
length(which(duonormdunnvals$duovcae>0.05))/10000 #caecum 0
length(which(duonormdunnvals$duovapp>0.05))/10000 #appendix 0
length(which(duonormdunnvals$duovasc>0.05))/10000 #ascending colon 0
length(which(duonormdunnvals$duovtc>0.05))/10000 #transverse colon 0
length(which(duonormdunnvals$duovdc>0.05))/10000 #descending colon 0
length(which(duonormdunnvals$duovsc>0.05))/10000 #sigmodial 0
length(which(duonormdunnvals$duovrect>0.05))/10000 #rectum 0

length(which(duocountdunnvals$duovjej>0.05))/10000 #duodenum 0
length(which(duocountdunnvals$duovi>0.05))/10000 #ileum 0.0046
length(which(duocountdunnvals$duovcae>0.05))/10000 #caecum 0
length(which(duocountdunnvals$duovapp>0.05))/10000 #appendix 0
length(which(duonormdunnvals$duovasc>0.05))/10000 #ascending colon 0
length(which(duocountdunnvals$duovtc>0.05))/10000 #transverse colon 0
length(which(duocountdunnvals$duovdc>0.05))/10000 #descending colon 0
length(which(duocountdunnvals$duovsc>0.05))/10000 #sigmodial 0
length(which(duocountdunnvals$duovrect>0.05))/10000 #rectum 0

length(which(ilenormdunnvals$ilevd>0.05))/10000 #duodenum  0.0046
length(which(ilenormdunnvals$ilevjej>0.05))/10000 #jejunum 0
length(which(ilenormdunnvals$ilevcae>0.05))/10000 #caecum  0.7328
length(which(ilenormdunnvals$ilevapp>0.05))/10000 #appendix 0.7277
length(which(ilenormdunnvals$ilevasc>0.05))/10000 #ascending colon 0.7366
length(which(ilenormdunnvals$ilevtc>0.05))/10000 #transverse colon  0.6951
length(which(ilenormdunnvals$ilevdc>0.05))/10000 #descending colon 0.8229
length(which(ilenormdunnvals$ilevsc>0.05))/10000 #sigmodial 0.7841
length(which(ilenormdunnvals$ilevrect>0.05))/10000 #rectum  0.7048

length(which(ilecountdunnvals$ilevd>0.05))/10000 #duodenum  0.0046
length(which(ilecountdunnvals$ilevjej>0.05))/10000 #jejunum 0
length(which(ilecountdunnvals$ilevcae>0.05))/10000 #caecum  0.7323
length(which(ilecountdunnvals$ilevapp>0.05))/10000 #appendix 0.7278
length(which(ilenormdunnvals$ilevasc>0.05))/10000 #ascending colon 0.7366
length(which(ilecountdunnvals$ilevtc>0.05))/10000 #transverse colon  0.6951
length(which(ilecountdunnvals$ilevdc>0.05))/10000 #descending colon 0.8227
length(which(ilecountdunnvals$ilevsc>0.05))/10000 #sigmodial 0.7837
length(which(ilecountdunnvals$ilevrect>0.05))/10000 #rectum 0.7044

OXTplot<-ggplot(totaldata2 ,aes(x=Region,y=OXTnorm, fill=Region))+
  geom_violin()+
  geom_jitter(data=totalpoints,  aes(x=Region,y=OXTpositive,fill=Region),pch=21,width = 0.2)+
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

pdf(file="Elmdataviolinplot.pdf",width=3,height=2.2)
OXTplot
dev.off()

