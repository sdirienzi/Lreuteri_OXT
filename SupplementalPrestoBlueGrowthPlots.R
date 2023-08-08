require(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(data.table)
library(stringr)
library("Polychrome")
library(PMCMRplus)
library(ggsci)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#growthcurves

growth<-as.data.table(read.delim(file="FigS2f_growthcurves.txt",header=TRUE,sep = "\t",check.names=FALSE))

#flip the data

growthlong<-melt(growth, id.vars=c("Replicate",	"Strain"),	 
                  measure.vars=c ("0","1",	"2",	"3",	"4",	"5",	"6","7", "8"),
                  variable.name="Time", value.name= "OD")

summarized<-summarySE(data=growthlong,measurevar = "OD", groupvars = c("Strain","Time"))

growthplots<-ggplot(summarized, aes(x=as.numeric(Time),y=OD,shape=as.factor(Strain))) +  
  geom_errorbar(aes(ymin=OD-sd, ymax=OD+sd,color=Strain), width=.1,alpha=0.5,) +
  #geom_ribbon(aes(ymin=OD-sd, ymax=OD+sd,fill=Strain),alpha=0.1)+
  geom_line(aes(colour=Strain)) +
  #geom_point(aes(colour=sugar),size=1) +
  #geom_text(data=averageblankedtable[name2==reactors[j]& time==max(averageblankedtable$time)], aes(x= time/60+.3,label=daynum),size=6)+  #toggle to get annotations
  #geom_text(averageblankedtable[name2==reactors[j]&&time==max(averageblankedtable$time)],aes(x=time/60,y=ODblanked,label=daynum) )+
  xlab("Hours") +
  ylab(bquote(paste('OD'['600']))) +
  scale_colour_manual(values = c("6475" = "violetred1", "B. subtilis" = "#629F87", "E. coli"="#455487")) + 
 # scale_fill_manual(values = c("trehalose" = "violetred1", "sucrose" = "#6A5ACD", "HFCS"="#FFA500")) + 
  theme_classic() +
  theme(
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=10,colour="black"),
    axis.title.y= element_text(size=10,colour="black"),
    axis.text.x= element_text(size=8,colour="black"),
    axis.text.y= element_text(size=8,colour="black")#,
    #legend.position="none"
  )

pdf(file="straingrowth.pdf",width=3.5,height=2)
growthplots
dev.off()


#PrestoBlue

Blue<-as.data.table(read.delim(file="FigS2e_PrestoBlue.txt",header=TRUE,sep="\t"))
Blue$InductionTreatment<-paste(Blue$Induction,Blue$Treatment)

Blue$InductionTreatment<-factor(Blue$InductionTreatment,c("Blank LDM4","Uninduced LDM4","Uninduced 6475",
                                                         "Induced LDM4","Induced 6475" ))
Blue$Treatment<-factor(Blue$Treatment,c("LDM4","6475"))
Blue$Induction<-factor(Blue$Induction,c("Blank", "Uninduced","Induced"))

Presto<-ggplot(Blue, aes(x=Induction,y=FoldChange,fill=Treatment,group=InductionTreatment)) +
  geom_boxplot(size=0.6,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(aes(shape=as.factor(Experiment)),size=1.5,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('FoldChange')) +
  scale_y_continuous(expand=c(0.05,.2,.2,.2))+
  scale_fill_npg()+
  scale_shape_manual(values=c("1"=22,"2"=23,"3"=24))+
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
    axis.text.x= element_text(size=7, colour = "black"), #angle=45,vjust=0.7
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )

pdf(file="prestoblue.pdf",width=2,height = 2)
Presto
dev.off()
