
require(reshape2)
library(ggplot2)
library(plyr)
library(gridExtra)
library(data.table)
library(stringr)
library("Polychrome")
library(PMCMRplus)
library(ggsci)
library(sjstats)
library(emmeans)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
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

#Fig S1C
mdt<-as.data.table(read.delim(file="FigS1c_qPCRdataD103J2NGN3.txt",header=TRUE,sep="\t",na.strings=""))
mdt$CT<-as.numeric(as.character(mdt$CT))
mdts<-as.data.table(summarySE(data=mdt[is.na(CT)==FALSE],measurevar = "CT",groupvars = c("CellType","Rep","treatmentname","Primer","Sample")),na.rm = TRUE)

#exclude those with < 2 N
mdts<-mdts[N>=2]

mdts$Sample = factor(mdts$Sample,c("Undiff","NoRTunDiff","Diff","NoRTDiff","NoDNA" ))

PGAPDH<-mdts[Primer =="GAPDH"]
names(PGAPDH)[which(names(PGAPDH)=="CT")]<-"CT.GAPDH"

mdts2<-merge(mdts,PGAPDH,by=c("CellType","Rep","Sample"),all.x=TRUE)
mdts3<-mdts2

mdts3$CTGAPDHnorm<-2^(mdts3$CT.GAPDH-mdts3$CT)

mdts3<-mdts3[Primer.x=="LGR5"|Primer.x=="SI"|Primer.x=="OXT"]
mdts3$Primer.x<-factor(mdts3$Primer.x,c("LGR5","SI","OXT"))

mdts3$FullSampleinfo<-paste(mdts3$CellType,mdts3$Sample,sep=" ")
#added extra elements to fill in the NA data as 1e-5
mdts3add<-read.delim(file="FigS1c_mdts3extras.txt",sep="\t",header=TRUE)
mdts4<-rbind(mdts3,mdts3add)

mdts4$FullSampleinfo<-factor(mdts4$FullSampleinfo,c("D103 Undiff",  "NGN3 Undiff", "D103 Diff","NGN3 Diff"))
samplecolours<-c("navyblue","slateblue1","maroon","#F0027F")
levels(mdts4$FullSampleinfo)

jitterwidth<-.75
#alphacol<-0.4
strokewidth<-0.5

PGAPDHplot<-ggplot(mdts4, aes(x=FullSampleinfo, y=CTGAPDHnorm,fill=FullSampleinfo,col=FullSampleinfo,group=FullSampleinfo) ) + 
  geom_point(size=1.5,pch=21,stroke=strokewidth, position=position_jitterdodge(jitter.width = jitterwidth-.25),color="black") + 
  geom_boxplot(size=0.6,fatten=1, outlier.color = NA,alpha=0.2,color="black") + 
  labs(y="GAPDH Normalized CN") +
  facet_wrap(~Primer.x,nrow=1,scales="fixed")+
  coord_trans(y="log10")+
  scale_fill_manual(name = "FullSampleinfo",values=samplecolours)+
 scale_x_discrete(labels=c("D103","J2NGN3","D103","J2NGN3","D103","J2NGN3","D103","J2NGN3","D103","J2NGN3","D103","J2NGN3"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4)), breaks =c(1e-5,1e-4,1e-3,1e-2,1e-1,1,7.5),labels=c("1e-5","1e-4","1e-3","0.01","0.1","1","7.5") )+
  theme_bw() + theme(
    plot.title = element_text(size = 7,hjust = 0.5),
    legend.position="none",
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=8),
    axis.title.y= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black",angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    strip.text = element_text(size=8),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background =element_blank())



pdf(file="D103NGN3OXTqPCR.pdf",width=3,height=2.25)
PGAPDHplot
dev.off()

mdts4<-as.data.table(mdts4)

SIlinear<-lmer(CTGAPDHnorm~Sample + (1|CellType),data=mdts4[Primer.x=="SI"],REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
SInull<-lmer(CTGAPDHnorm~(1|CellType),data=mdts4[Primer.x=="SI"],REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
anova(SIlinear,SInull)
# Data: mdts4[Primer.x == "SI"]
# Models:
#   SInull: CTGAPDHnorm ~ (1 | CellType)
# SIlinear: CTGAPDHnorm ~ Sample + (1 | CellType)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)   
# SInull      3 -30.157 -28.702 18.078  -36.157                        
# SIlinear    4 -34.809 -32.869 21.404  -42.809 6.6524  1   0.009902 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
performance::r2(SIlinear)
# Marginal R2: 0.447

SIemmeans<-emmeans(SIlinear,pairwise~Sample,adjust="none")
# $emmeans
# Sample   emmean     SE   df lower.CL upper.CL
# Undiff 0.000349 0.0182 5.18  -0.0459   0.0466
# Diff   0.070333 0.0182 5.18   0.0241   0.1166
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate     SE df t.ratio p.value
# Undiff - Diff    -0.07 0.0257 13  -2.722  0.0175
# 
# Degrees-of-freedom method: kenward-roger 

eff_size(SIemmeans,sigma = sigma(SIlinear),edf=13)
# contrast        effect.size    SE df lower.CL upper.CL
# (Undiff - Diff)       -1.72 0.717 13    -3.27   -0.172
# 
# sigma used for effect sizes: 0.04065 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

OXTlinear<-lmer(CTGAPDHnorm~Sample + (1|CellType),data=mdts4[Primer.x=="OXT"],REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
OXTnull<-lmer(CTGAPDHnorm~(1|CellType),data=mdts4[Primer.x=="OXT"],REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
anova(OXTlinear,OXTnull)
# Data: mdts4[Primer.x == "OXT"]
# Models:
#   OXTnull: CTGAPDHnorm ~ (1 | CellType)
# OXTlinear: CTGAPDHnorm ~ Sample + (1 | CellType)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)  
# OXTnull      3 -154.93 -153.47 80.463  -160.93                       
# OXTlinear    4 -158.91 -156.97 83.453  -166.91 5.9807  1    0.01446 *

performance::r2(OXTlinear)

# Conditional R2: 0.452
# Marginal R2: 0.40

OXTemmeans<-emmeans(OXTlinear,pairwise~Sample,adjust="none")
# $emmeans
# Sample   emmean       SE   df  lower.CL upper.CL
# Undiff 4.54e-05 0.000131 9.57 -0.000248 0.000338
# Diff   4.13e-04 0.000131 9.57  0.000120 0.000706
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast       estimate       SE   df t.ratio p.value
# Undiff - Diff -0.000368 0.000136 11.1  -2.706  0.0203
# 
# Degrees-of-freedom method: kenward-roger

eff_size(OXTemmeans,sigma = sigma(OXTlinear),edf=11.1)
# contrast        effect.size    SE   df lower.CL upper.CL
#Undiff - Diff)       -1.65 0.702 11.1    -3.19   -0.104
# 
# sigma used for effect sizes: 0.0002232 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

LGR5linear<-lmer(CTGAPDHnorm~Sample+(1|CellType),data=mdts4[Primer.x=="LGR5"], REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
LGR5null<-lmer(CTGAPDHnorm~(1|CellType),data=mdts4[Primer.x=="LGR5"], REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
anova(LGR5linear,LGR5null)
performance::r2(LGR5linear)
# Conditional R2: 0.577
# Marginal R2: 0.505

LGR5emmeans<-emmeans(LGR5linear,pairwise~Sample,adjust="none")
# $emmeans
# Sample  emmean      SE   df  lower.CL upper.CL
# Undiff 0.00345 0.00108 8.19  0.000971  0.00593
# Diff   0.00001 0.00108 8.19 -0.002469  0.00249
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate    SE   df t.ratio p.value
# Undiff - Diff  0.00344 0.001 11.1   3.438  0.0055
# 
# Degrees-of-freedom method: kenward-roger 

eff_size(LGR5emmeans,sigma = sigma(LGR5linear),edf=11.1)
# contrast        effect.size    SE   df lower.CL upper.CL
# (Undiff - Diff)        2.09 0.753 11.1    0.436     3.75
# 
# sigma used for effect sizes: 0.001644 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 
