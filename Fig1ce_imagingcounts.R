#Imaging quantification

require(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(data.table)
library(stringr)
library("Polychrome")
library(PMCMRplus)
library(ggsci)
library(lme4)
library(emmeans)
library(PMCMRplus)


humandata<-read.delim(file="Fig1c.txt")
humandata$Segment<-factor(humandata$Segment,c("Fundus","Body","Antrum","Duodenum","UpperJejunum","MidJejunum",
                                                          "LowerJejunum","UpperIleum","MidIleum","LowerIleum",
                                                          "Ascolon","Transcolon", "Descolon","Sigcolon"))
shapevalues<-c(22,23,24,25)


Human<-ggplot(humandata, aes(x=Segment,y=PercentOXT,fill=Segment)) +
  geom_boxplot(size=0.5,outlier.shape = NA,colour="black")+  #, lty=FAMI
  geom_point(aes(shape=LG),size=1,stroke=0.2 ,color="black" )+ # position=position_jitterdodge(jitter.width=0.2 )
  ylab(bquote('Percent oxytocin cells')) +
  scale_y_continuous(expand=c(0.05,.4,.4,.4))+
  scale_fill_manual(values=c("#E64B3533","#E64B357F","#E64B35CC","#E64B35FF",
                      "#4DBBD54C","#4DBBD5FF","#4DBBD599",
                      "#00A0874C","#00A087FF","#00A08799",
                      "#8491B4FF", "#91D1C2FF","#DC0000FF", "#7E6148FF"))+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.x= element_text(size=8),
    axis.title.y= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black", angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )

pdf(file="Humanimaging.pdf",width=3.5,height=2)
Human
dev.off()




#mouse

#load Nuclei counts

mousedata<-read.delim(file="Fig1e.txt",header=TRUE,sep="\t")

mousedata$Segment<-factor(mousedata$Segment,c("Stomach","ProxSI","DistalSI", "Cecum",
                                                      "ProxLI","MidLI", "DistalLI"))
shapevalues<-c(22,23,24,25)

Mouse<-ggplot(mousedata, aes(x=Segment,y=PercentOXT,fill=Segment)) +
  geom_boxplot(size=0.5,outlier.shape = NA,colour="black")+  #, lty=FAMI
  geom_point(aes(shape=Sex),size=1,stroke=0.2 ,color="black" )+ # position=position_jitterdodge(jitter.width=0.2 )
 # facet_wrap(~Sex)+
  ylab(bquote('Percent oxytocin cells')) +
  scale_y_continuous(expand=c(0.05,.4,.4,.4))+
  scale_fill_manual(values=c("#E64B35CC",
                             "#4DBBD54C",
                             "#00A08799",
                             "#8491B4FF", "#91D1C2FF","#DC0000FF", "#7E6148FF"))+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.x= element_text(size=8),
    axis.title.y= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black", angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )

pdf(file="Mouseimaging.pdf",width=3,height=2)
Mouse
dev.off()
