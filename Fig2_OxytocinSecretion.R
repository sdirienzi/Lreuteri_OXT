#
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
library(sjstats)


#2B
LGLreusecret<-read.delim(file="Fig2b_AllLifeGiftsecretionLreu.txt",header=TRUE,sep="\t",check.names = FALSE)
LGtissuesize<-read.delim(file="Fig2b_LifeGiftTissuesizes.txt",header = TRUE,sep="\t")

LGLreusecmerge<-merge(LGLreusecret,LGtissuesize,by=c("LifeGift","Segment"))
LGLreusecmerge$Oxtnorm<-LGLreusecmerge$Concentration/LGLreusecmerge$Area
jitterwidth=.18
shapevalues<-c(21,22,23,24,25)
names(shapevalues)<-levels(as.factor(LGLreusecmerge$LifeGift))
LGLreusecmerge$Treatment<-factor(LGLreusecmerge$Treatment,c("LDM4","6475"))
LGLreusecmerge$SegmentTreatment=paste(LGLreusecmerge$Segment,LGLreusecmerge$Treatment,sep="_")
LGLreusecmerge$Segment<-factor(LGLreusecmerge$Segment,c("Duodenum","UpperJejunum","LowerJejunum",
                                                        "UpperIleum","LowerIleum", "AscendingColon","TransverseColon", "DescendingColon" ))
LGLreusecmerge$SegmentTreatment<-factor(LGLreusecmerge$SegmentTreatment,c("Duodenum_LDM4", "Duodenum_6475","UpperJejunum_LDM4", "UpperJejunum_6475" ,
                                        "LowerJejunum_LDM4","LowerJejunum_6475",
                                        "UpperIleum_LDM4", "UpperIleum_6475",
                                        "LowerIleum_LDM4","LowerIleum_6475",
  "AscendingColon_LDM4","AscendingColon_6475","TransverseColon_LDM4","TransverseColon_6475", 
                                        "DescendingColon_LDM4", "DescendingColon_6475"))


#pdf(file="LifeGiftLreu.pdf",width=3.6,height=2.8)
LifeGiftLreuplot<-ggplot(subset(LGLreusecmerge, LGLreusecmerge$Treatment!="4659O" & LGLreusecmerge$Treatment!="4659H"), aes(x=Segment,y=Oxtnorm,fill=Treatment,group=SegmentTreatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(aes(shape=LifeGift),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')')) +
  scale_y_continuous(expand=c(0.05,.2,.2,.2))+
  scale_fill_npg()+
  scale_shape_manual(values=shapevalues)+
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


pdf(file="LifeGiftLreusecretion.pdf",width=6,height=2)
grid.arrange(LifeGiftLreuplot)
dev.off()

humandata<-LGLreusecmerge
mmtest<-lmer(Oxtnorm~Treatment*Segment +(1|LifeGift),data=humandata, REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
mmtestnull<-lmer(Oxtnorm~(1|LifeGift),data=humandata, REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
anova(mmtest,mmtestnull)
summary(mmtest)
# Linear mixed model fit by maximum likelihood  ['lmerMod']
# Formula: Oxtnorm ~ Treatment * Segment + (1 | LifeGift)
# Data: humandata
# Control: lmerControl(optimizer = "bobyqa")
# 
# AIC      BIC   logLik deviance df.resid 
# 2343.7   2405.0  -1153.9   2307.7      204 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.5138 -0.5892 -0.0650  0.4952  4.7781 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# LifeGift (Intercept)  462.7   21.51   
# Residual             1809.9   42.54   
# Number of obs: 222, groups:  LifeGift, 5
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                           24.1088    14.6015   1.651
# Treatment6475                         64.6259    15.5345   4.160
# SegmentUpperJejunum                   -4.8743    15.5345  -0.314
# SegmentLowerJejunum                   -0.7782    15.5345  -0.050
# SegmentUpperIleum                     -7.0960    15.5345  -0.457
# SegmentLowerIleum                    -14.9818    15.5345  -0.964
# SegmentAscendingColon                -11.4809    16.5576  -0.693
# SegmentTransverseColon               -15.9065    16.5576  -0.961
# SegmentDescendingColon               -11.7237    16.5576  -0.708
# Treatment6475:SegmentUpperJejunum     29.9538    21.9691   1.363
# Treatment6475:SegmentLowerJejunum     57.9921    21.9691   2.640
# Treatment6475:SegmentUpperIleum        1.7188    21.9691   0.078
# Treatment6475:SegmentLowerIleum       -1.1548    21.9691  -0.053
# Treatment6475:SegmentAscendingColon  -47.3916    23.3018  -2.034
# Treatment6475:SegmentTransverseColon -48.1649    23.3018  -2.067
# Treatment6475:SegmentDescendingColon -55.4810    23.3018  -2.381
performance::r2.merMod(mmtest)
# Conditional R2: 0.555
# Marginal R2: 0.441

mmtest.emmeans<-emmeans(mmtest,pairwise~Treatment | Segment,adjust="BH")
# $emmeans
# Segment = Duodenum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       24.11 15.7 26.9    -8.13     56.3
# 6475       88.73 15.7 26.9    56.49    121.0
# 
# Segment = UpperJejunum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       19.23 15.7 26.9   -13.01     51.5
# 6475      113.81 15.7 26.9    81.57    146.1
# 
# Segment = LowerJejunum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       23.33 15.7 26.9    -8.91     55.6
# 6475      145.95 15.7 26.9   113.71    178.2
# 
# Segment = UpperIleum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       17.01 15.7 26.9   -15.23     49.3
# 6475       83.36 15.7 26.9    51.12    115.6
# 
# Segment = LowerIleum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4        9.13 15.7 26.9   -23.11     41.4
# 6475       72.60 15.7 26.9    40.36    104.8
# 
# Segment = AscendingColon:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       12.63 16.8 34.9   -21.49     46.7
# 6475       29.86 16.8 34.9    -4.26     64.0
# 
# Segment = TransverseColon:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4        8.20 16.8 34.9   -25.92     42.3
# 6475       24.66 16.8 34.9    -9.46     58.8
# 
# Segment = DescendingColon:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       12.39 16.8 34.9   -21.74     46.5
# 6475       21.53 16.8 34.9   -12.59     55.7
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Segment = Duodenum:
#   contrast    estimate   SE  df t.ratio p.value
# LDM4 - 6475   -64.63 16.1 233  -4.014  0.0001
# 
# Segment = UpperJejunum:
#   contrast    estimate   SE  df t.ratio p.value
# LDM4 - 6475   -94.58 16.1 233  -5.874  <.0001
# 
# Segment = LowerJejunum:
#   contrast    estimate   SE  df t.ratio p.value
# LDM4 - 6475  -122.62 16.1 233  -7.616  <.0001
# 
# Segment = UpperIleum:
#   contrast    estimate   SE  df t.ratio p.value
# LDM4 - 6475   -66.34 16.1 233  -4.121  0.0001
# 
# Segment = LowerIleum:
#   contrast    estimate   SE  df t.ratio p.value
# LDM4 - 6475   -63.47 16.1 233  -3.942  0.0001
# 
# Segment = AscendingColon:
#   contrast    estimate   SE  df t.ratio p.value
# LDM4 - 6475   -17.23 18.0 233  -0.957  0.3394
# 
# Segment = TransverseColon:
#   contrast    estimate   SE  df t.ratio p.value
# LDM4 - 6475   -16.46 18.0 233  -0.914  0.3614
# 
# Segment = DescendingColon:
#   contrast    estimate   SE  df t.ratio p.value
# LDM4 - 6475    -9.14 18.0 233  -0.508  0.6119

eff_size(mmtest.emmeans,sigma=sigma(mmtest),edf=233)
# Segment = Duodenum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (LDM4 - 6475)      -1.519 0.385 233    -2.28   -0.761
# 
# Segment = UpperJejunum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (LDM4 - 6475)      -2.223 0.392 233    -3.00   -1.450
# 
# Segment = LowerJejunum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (LDM4 - 6475)      -2.882 0.401 233    -3.67   -2.092
# 
# Segment = UpperIleum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (LDM4 - 6475)      -1.559 0.385 233    -2.32   -0.800
# 
# Segment = LowerIleum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (LDM4 - 6475)      -1.492 0.385 233    -2.25   -0.734
# 
# Segment = AscendingColon:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (LDM4 - 6475)      -0.405 0.424 233    -1.24    0.429
# 
# Segment = TransverseColon:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (LDM4 - 6475)      -0.387 0.424 233    -1.22    0.447
# 
# Segment = DescendingColon:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (LDM4 - 6475)      -0.215 0.423 233    -1.05    0.619


Fig2C

pig5<-read.delim(file="Fig2d_pig5xdata_forR.txt",header=TRUE,sep="\t")
pigmelt<-as.data.table(pig5)

#reorder
pigmelt$Treatment<-factor(pigmelt$Treatment,c("Krebs","LDM4","6475"))

#these data are 5x concentrated so right divide everything by 5

pigmelt$OXTfoldscaled<-NA
pigmelt$OXTareafoldscaled<-NA
for (i in 1:length(pigmelt$Treatment)) {
  pigmelt$OXTfoldscaled[i]<-pigmelt$OXTcorrected[i]/5
   pigmelt$OXTareafoldscaled[i]<-pigmelt$OXTfoldscaled[i]/7.07
}

pigmelt$Segment<-factor(pigmelt$Segment,c("Duodenum","Jejunum","Ileum","Ascending Colon", "Descending Colon"))
pigmelt$RegionTreatment<-paste(pigmelt$Segment,pigmelt$Treatment,sep="_")
pigmelt$RegionTreatment<-factor(pigmelt$RegionTreatment,c("Duodenum_Krebs","Duodenum_LDM4","Duodenum_6475",
                                                          "Jejunum_Krebs","Jejunum_LDM4", "Jejunum_6475",
                                                          "Ileum_Krebs","Ileum_LDM4", "Ileum_6475",
                                                          "Ascending Colon_Krebs", "Ascending Colon_LDM4","Ascending Colon_6475",
                                                          "Descending Colon_Krebs",  "Descending Colon_LDM4" , "Descending Colon_6475"))

pigplot<-ggplot(pigmelt[Segment=="Jejunum"& Treatment!="Krebs"], aes(x=Segment, y=OXTareafoldscaled,fill=Treatment,group=RegionTreatment)) +  #colour=Diet,fill=Experiment, shape=Rep
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+
  geom_point(aes(shape=as.factor(Pig)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  xlab("Segment") +
  scale_fill_npg()+
  scale_shape_manual(values=c(21,22,25))+
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')')) +
  ggtitle("Mouse small intestine")+
  coord_cartesian(
    ylim = c(0,	28))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=8),
    axis.title.x= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    plot.title=element_blank(), legend.position="none" 
  )


#LM
pigmeltnoKrebs<-pigmelt[Treatment!="Krebs"]
pigmeltnoKrebs<-droplevels(pigmeltnoKrebs)
levels(pigmeltnoKrebs$Treatment:pigmeltnoKrebs$Segment)
                               
mmtest<-lmer(OXTareafoldscaled~Treatment*Segment +(1|Pig),data=pigmeltnoKrebs, REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
mmnull<-lmer(OXTareafoldscaled~(1|Pig),data=pigmeltnoKrebs, REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
anova(mmtest,mmnull)
performance::r2(mmtest)
# R2 for Mixed Models

# Conditional R2: 0.746
# Marginal R2: 0.721

summary(mmtest)
# Linear mixed model fit by maximum likelihood  ['lmerMod']
# Formula: OXTareafoldscaled ~ Treatment * Segment + (1 | Pig)
# Data: pigmeltnoKrebs
# Control: lmerControl(optimizer = "bobyqa")
# 
# AIC      BIC   logLik deviance df.resid 
# 509.4    539.4   -242.7    485.4       78 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.7566 -0.3500  0.0870  0.2078  5.2377 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Pig      (Intercept)  1.23    1.109   
# Residual             12.30    3.507   
# Number of obs: 90, groups:  Pig, 3
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                             0.43723    1.33300   0.328
# Treatment6475                          14.69319    1.65343   8.886
# SegmentJejunum                          0.18350    1.65343   0.111
# SegmentIleum                           -0.16875    1.65343  -0.102
# SegmentAscending Colon                 -0.09915    1.65343  -0.060
# SegmentDescending Colon                -0.06545    1.65343  -0.040
# Treatment6475:SegmentJejunum            0.69406    2.33830   0.297
# Treatment6475:SegmentIleum            -10.23443    2.33830  -4.377
# Treatment6475:SegmentAscending Colon  -13.89197    2.33830  -5.941
# Treatment6475:SegmentDescending Colon -13.14373    2.33830  -5.621
# 
# Correlation of Fixed Effects:
#   (Intr) Tr6475 SgmntJ SgmntI SgmnAC SgmnDC T6475:SJ T6475:SI T6475:SAC
# Tretmnt6475 -0.620                                                               
# SegmentJjnm -0.620  0.500                                                        
# SegmentIlem -0.620  0.500  0.500                                                 
# SgmntAscndC -0.620  0.500  0.500  0.500                                          
# SgmntDscndC -0.620  0.500  0.500  0.500  0.500                                   
# Trtm6475:SJ  0.439 -0.707 -0.707 -0.354 -0.354 -0.354                            
# Trtm6475:SI  0.439 -0.707 -0.354 -0.707 -0.354 -0.354  0.500                     
# Trt6475:SAC  0.439 -0.707 -0.354 -0.354 -0.707 -0.354  0.500    0.500            
# Trt6475:SDC  0.439 -0.707 -0.354 -0.354 -0.354 -0.707  0.500    0.500    0.500   

mmtest.emmeans<-emmeans(mmtest,pairwise~Treatment | Segment,adjust="BH")

# $emmeans
# Segment = Duodenum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       0.437 1.48 38.5    -2.56     3.43
# 6475      15.130 1.48 38.5    12.13    18.13
# 
# Segment = Jejunum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       0.621 1.48 38.5    -2.38     3.62
# 6475      16.008 1.48 38.5    13.01    19.00
# 
# Segment = Ileum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       0.268 1.48 38.5    -2.73     3.26
# 6475       4.727 1.48 38.5     1.73     7.72
# 
# Segment = Ascending Colon:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       0.338 1.48 38.5    -2.66     3.33
# 6475       1.139 1.48 38.5    -1.86     4.14
# 
# Segment = Descending Colon:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       0.372 1.48 38.5    -2.62     3.37
# 6475       1.921 1.48 38.5    -1.07     4.92
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Segment = Duodenum:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - 6475  -14.693 1.75 97  -8.414  <.0001
# 
# Segment = Jejunum:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - 6475  -15.387 1.75 97  -8.812  <.0001
# 
# Segment = Ileum:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - 6475   -4.459 1.75 97  -2.553  0.0122
# 
# Segment = Ascending Colon:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - 6475   -0.801 1.75 97  -0.459  0.6474
# 
# Segment = Descending Colon:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - 6475   -1.549 1.75 97  -0.887  0.3771
# 
# Degrees-of-freedom method: kenward-roger 

eff_size(mmtest.emmeans,sigma=sigma(mmtest),edf=97)
# Segment = Duodenum:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - 6475)      -4.189 0.582 97    -5.34   -3.035
# 
# Segment = Jejunum:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - 6475)      -4.387 0.589 97    -5.56   -3.218
# 
# Segment = Ileum:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - 6475)      -1.271 0.506 97    -2.28   -0.267
# 
# Segment = Ascending Colon:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - 6475)      -0.228 0.498 97    -1.22    0.760
# 
# Segment = Descending Colon:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - 6475)      -0.442 0.499 97    -1.43    0.548
# 
# sigma used for effect sizes: 3.507 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

#Fig2D

piglet<-read.delim(file="Piglet_2.5and5xdataforR.txt",header=TRUE,sep="\t")
pigletmelt<-as.data.table(piglet)
pigletmelt$OXTfoldscaled<-NA
pigletmelt$OXTareafoldscaled<-NA

#for ascending colon 2.5 x 1 = 2.5
#duodenum, jejunum and ileum 2 x 3
for (i in 1:length(pigletmelt$Treatment)) {
 if (pigletmelt$Concentration[i]==5.0) {
    pigletmelt$OXTfoldscaled[i]<-pigletmelt$OXTcorrected[i]/5
  }
  if (pigletmelt$Intestinal.Segment[i]=="Ascending Colon") {
    #2.5
    pigletmelt$OXTareafoldscaled[i]<-pigletmelt$OXTfoldscaled[i]/2.5
  }
  else  {
    pigletmelt$OXTareafoldscaled[i]<-pigletmelt$OXTfoldscaled[i]/6
  }
}
pigletmelt$RegionTreatment<-pigletmelt$RegionTreatment<-paste(pigletmelt$Intestinal.Segment,pigletmelt$Treatment,sep="_")
pigletmelt$RegionTreatment<-factor(pigletmelt$RegionTreatment,c("Duodenum_LDM4", "Duodenum_6475",
                                                                "Jejunum_LDM4","Jejunum_6475",
                                                                "Ileum_LDM4", "Ileum_6475",
                                                                "Ascending Colon_LDM4", "Ascending Colon_6475"))

pigletmelt$Intestinal.Segment<-factor(pigletmelt$Intestinal.Segment,c("Duodenum","Jejunum","Ileum","Ascending Colon"))
pigletmelt$Treatment<-factor(pigletmelt$Treatment,c("LDM4","6475"))
pigletmelt$Treatment<-factor(pigletmelt$Treatment,c("LDM4","6475"))


shapevalues3<-c(21,22)
names(shapevalues3)<-c("Female","Male")

pigletplot<-ggplot(pigletmelt[Intestinal.Segment=="Jejunum"& Concentration==5], aes(x=Intestinal.Segment, y=OXTareafoldscaled,fill=Treatment,group=RegionTreatment)) +  #colour=Diet,fill=Experiment, shape=Rep
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+
  geom_point(aes(shape=as.factor(Sex)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
   xlab("Segment") +
  scale_fill_npg()+
  scale_shape_manual(values=shapevalues3)+
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')')) +
  ggtitle("Piglet small intestine")+
  coord_cartesian(
    ylim = c(0,	28))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=8),
    axis.title.x= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    plot.title=element_blank(), legend.position="none" 
  )

pigletdata<-pigletmelt
pigletdata<-droplevels(pigletdata)
pigletdata$Treatment<-(as.factor(pigletdata$Treatment))
pigletdata$Intestinal.Segment<-as.factor(pigletdata$Intestinal.Segment)
#note for the piglets, there are just two piglets, one male, one female
mmtest<-lmer(OXTareafoldscaled~Treatment*Intestinal.Segment +(1|Sex),data=pigletdata, REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
mmnull<-lmer(OXTareafoldscaled~(1|Sex),data=pigletdata, REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
anova(mmtest,mmnull)
performance::r2(mmtest)
summary(mmtest)
# Linear mixed model fit by maximum likelihood  ['lmerMod']
# Formula: OXTareafoldscaled ~ Treatment * Intestinal.Segment + (1 | Sex)
# Data: pigletdata
# Control: lmerControl(optimizer = "bobyqa")
# 
# AIC      BIC   logLik deviance df.resid 
# 249.9    266.3   -115.0    229.9       28 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.1945 -0.5922 -0.1977  0.5144  3.3324 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Sex      (Intercept)  9.005   3.001   
# Residual             22.178   4.709   
# Number of obs: 38, groups:  Sex, 2
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                                       1.2155     2.9932   0.406
# Treatment6475                                    20.9702     2.9785   7.041
# Intestinal.SegmentJejunum                        -0.3703     2.9785  -0.124
# Intestinal.SegmentIleum                          -0.3222     2.9785  -0.108
# Intestinal.SegmentAscending Colon                 2.6956     3.1625   0.852
# Treatment6475:Intestinal.SegmentJejunum          -5.0896     4.2122  -1.208
# Treatment6475:Intestinal.SegmentIleum            -7.5091     4.2122  -1.783
# Treatment6475:Intestinal.SegmentAscending Colon  -1.7160     4.4677  -0.384
# 
# Correlation of Fixed Effects:
#   (Intr) Tr6475 Int.SJ Int.SI In.SAC T6475:I.SJ T6475:I.SI
# Tretmnt6475 -0.498                                                  
# Intstnl.SgJ -0.498  0.500                                           
# Intstnl.SgI -0.498  0.500  0.500                                    
# Intstnl.SAC -0.471  0.471  0.471  0.471                             
# Tr6475:I.SJ  0.352 -0.707 -0.707 -0.354 -0.333                      
# Tr6475:I.SI  0.352 -0.707 -0.354 -0.707 -0.333  0.500               
# T6475:I.SAC  0.332 -0.667 -0.333 -0.333 -0.706  0.471      0.471  

mmtest.emmeans<-emmeans(mmtest,pairwise~Treatment |Intestinal.Segment)

# $emmeans
# Intestinal.Segment = Duodenum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       1.215 3.88 11.8    -7.25     9.68
# 6475      22.186 3.88 11.8    13.72    30.65
# 
# Intestinal.Segment = Jejunum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       0.845 3.88 11.8    -7.62     9.31
# 6475      16.726 3.88 11.8     8.26    25.19
# 
# Intestinal.Segment = Ileum:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       0.893 3.88 11.8    -7.57     9.36
# 6475      14.354 3.88 11.8     5.89    22.82
# 
# Intestinal.Segment = Ascending Colon:
#   Treatment emmean   SE   df lower.CL upper.CL
# LDM4       3.911 4.04 14.3    -4.74    12.57
# 6475      23.165 4.04 14.3    14.51    31.82
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Intestinal.Segment = Duodenum:
#   contrast    estimate   SE   df t.ratio p.value
# LDM4 - 6475    -21.0 3.32 44.7  -6.319  <.0001
# 
# Intestinal.Segment = Jejunum:
#   contrast    estimate   SE   df t.ratio p.value
# LDM4 - 6475    -15.9 3.32 44.7  -4.786  <.0001
# 
# Intestinal.Segment = Ileum:
#   contrast    estimate   SE   df t.ratio p.value
# LDM4 - 6475    -13.5 3.32 44.7  -4.056  0.0002
# 
# Intestinal.Segment = Ascending Colon:
#   contrast    estimate   SE   df t.ratio p.value
# LDM4 - 6475    -19.3 3.71 44.7  -5.190  <.0001
# 
# Degrees-of-freedom method: kenward-roger 

eff_size(mmtest.emmeans,sigma=sigma(mmtest),edf=44.7)
# Intestinal.Segment = Duodenum:
#   contrast      effect.size    SE   df lower.CL upper.CL
# (LDM4 - 6475)       -4.45 0.848 44.7    -6.16    -2.75
# 
# Intestinal.Segment = Jejunum:
#   contrast      effect.size    SE   df lower.CL upper.CL
# (LDM4 - 6475)       -3.37 0.790 44.7    -4.96    -1.78
# 
# Intestinal.Segment = Ileum:
#   contrast      effect.size    SE   df lower.CL upper.CL
# (LDM4 - 6475)       -2.86 0.767 44.7    -4.40    -1.31
# 
# Intestinal.Segment = Ascending Colon:
#   contrast      effect.size    SE   df lower.CL upper.CL
# (LDM4 - 6475)       -4.09 0.899 44.7    -5.90    -2.28
# 
# sigma used for effect sizes: 4.709 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 


#Fig 2E
#load in mouse
mouse<-read.delim(file="Fig2e_CombinedMouseSecretionData.txt",sep="\t",header=TRUE)
mousearea<-read.delim(file="Fig2e_mousearea.txt",sep="\t",header=TRUE)

mouse<-merge(mouse,mousearea,by="Segment",x.all=TRUE)

mouse$Treatment<-factor(mouse$Treatment,c("LDM4","Lreu"))
mouse$Segment<-factor(mouse$Segment,c("Stomach","SI","Cecum","LI","Controls"))
mouse$SegmentTreatment<-paste(mouse$Segment,mouse$Treatment)
mouse$OXTareafoldscaled<-mouse$ConcentrationSTDremovelow/mouse$Area

shapevalues3<-c(21,22)
names(shapevalues3)<-c("F","M")

#plot just the cecum
mousececum<-ggplot(subset(mouse,mouse$Segment=="Cecum"), aes(x=Segment, y=OXTareafoldscaled,fill=Treatment,group=SegmentTreatment, shape=Sex)) +  
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  #, lty=FAMI
  geom_point(size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black")+
  scale_y_continuous(expand=c(0.05,.2,.3,.3))+
  scale_shape_manual(values=shapevalues3)+
  scale_x_discrete(labels=c("Cecum"))+
  scale_fill_npg()+
  xlab("Segment") +
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')')) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.y= element_text(size=8),
    axis.title.x= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    plot.title=element_blank(), legend.position="none" 
  )

mousetest<-lm(OXTareafoldscaled~Treatment*Segment*Sex ,data=mouse)
summary(mousetest)
# Call:
#   lm(formula = OXTareafoldscaled ~ Treatment * Segment * Sex, data = mouse)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -13.3103  -2.4524  -0.4421   2.3868  30.6513 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      12.1423     2.6141   4.645 1.32e-05 ***
#   TreatmentLreu                    14.2526     3.6969   3.855 0.000232 ***
#   SegmentSI                       -10.3518     3.6969  -2.800 0.006402 ** 
#   SegmentCecum                     -2.6223     3.6969  -0.709 0.480186    
# SegmentLI                        -0.1715     3.6969  -0.046 0.963121    
# SexM                             11.8782     3.6969   3.213 0.001894 ** 
#   TreatmentLreu:SegmentSI         -14.4789     5.2282  -2.769 0.006979 ** 
#   TreatmentLreu:SegmentCecum       -5.3868     5.2282  -1.030 0.305961    
# TreatmentLreu:SegmentLI          -9.8166     5.2282  -1.878 0.064078 .  
# TreatmentLreu:SexM               -6.2859     5.2282  -1.202 0.232794    
# SegmentSI:SexM                  -12.3990     5.2282  -2.372 0.020117 *  
#   SegmentCecum:SexM               -15.3893     5.2282  -2.944 0.004246 ** 
#   SegmentLI:SexM                  -11.2042     5.2282  -2.143 0.035150 *  
#   TreatmentLreu:SegmentSI:SexM      7.7393     7.3938   1.047 0.298381    
# TreatmentLreu:SegmentCecum:SexM   4.8930     7.3938   0.662 0.510021    
# TreatmentLreu:SegmentLI:SexM      5.3656     7.3938   0.726 0.470150    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.403 on 80 degrees of freedom
# Multiple R-squared:  0.7006,	Adjusted R-squared:  0.6445 
# F-statistic: 12.48 on 15 and 80 DF,  p-value: 2.932e-15

mousececemmeans<-emmeans(mousetest,pairwise~Treatment | Segment,adjust="BH")
# $contrasts
# Segment = Stomach:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - Lreu   -11.11 2.61 80  -4.250  0.0001
# 
# Segment = SI:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - Lreu    -0.50 2.61 80  -0.191  0.8487
# 
# Segment = Cecum:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - Lreu    -8.17 2.61 80  -3.125  0.0025
# 
# Segment = LI:
#   contrast    estimate   SE df t.ratio p.value
# LDM4 - Lreu    -3.98 2.61 80  -1.521  0.1322
# 
# Results are averaged over the levels of: Sex 

eff_size(mousececemmeans, sigma=sigma(mousetest),edf=80)
# Segment = Stomach:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -1.7350 0.431 80   -2.592   -0.878
# 
# Segment = SI:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -0.0781 0.408 80   -0.891    0.734
# 
# Segment = Cecum:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -1.2758 0.421 80   -2.113   -0.439
# 
# Segment = LI:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -0.6209 0.411 80   -1.439    0.197
# 
# Results are averaged over the levels of: Sex 
# sigma used for effect sizes: 6.403 
# Confidence level used: 0.95 

pdf("PigPigletMouse.pdf", width=3,height=1.5)
grid.arrange(pigplot,pigletplot,mousececum,ncol=3)        
dev.off()

#Fig2F
luminex<-read.delim(file="Fig2f_Luminex6475",header=TRUE,sep="\t",check.names=FALSE)

luminexmelt<-reshape2::melt(luminex,id.vars=names(luminex)[c(1,2)],variable="Induction",value.name = "Oxytocin")
shapevalues2<-c(21,22,23)
names(shapevalues2)<-levels(as.factor(luminexmelt$Replicate))
luminexmelt$InductionTreatment=paste(luminexmelt$Induction,luminexmelt$Treatment,sep="_")

pdf(file="LuminexLreu.pdf",width=1.5,height=1.5)
ggplot(luminexmelt, aes(x=Induction,y=Oxytocin,fill=Treatment,group=InductionTreatment)) +
  geom_boxplot(size=1,outlier.shape = NA,alpha = 0.4,colour="black")+  #, lty=FAMI
  geom_point(aes(shape=as.factor(Replicate)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  # ylim(c(0,165))+
  coord_cartesian(
    ylim = c(0,	400))+
  scale_fill_npg()+
  scale_shape_manual(values=shapevalues2)+
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
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )
dev.off()

HIElm<-lmer(Oxytocin~Treatment*Induction+(1|Replicate),data=luminexmelt,REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
HIEnull<-lmer(Oxytocin~(1|Replicate),data=luminexmelt,REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
anova(HIElm,HIEnull)
performance::r2(HIElm)
# Conditional R2: NA
# Marginal R2: 0.836

# summary(HIElm)
HIEemmeans<-emmeans(HIElm,pairwise~Treatment|Induction,adjust="BH")
# $emmeans
# Induction = Uninduced:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4        17.2 17.4 16    -19.8     54.1
# Lreu6475    65.9 17.4 16     29.0    102.9
# 
# Induction = Induced:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4        24.7 17.4 16    -12.2     61.7
# Lreu6475   230.2 17.4 16    193.3    267.2
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Induction = Uninduced:
#   contrast        estimate   SE df t.ratio p.value
# LDM4 - Lreu6475    -48.8 24.7 18  -1.978  0.0634
# 
# Induction = Induced:
#   contrast        estimate   SE df t.ratio p.value
# LDM4 - Lreu6475   -205.5 24.7 18  -8.335  <.0001

eff_size(HIEemmeans,sigma = sigma(HIElm),edf=18)
# Induction = Uninduced:
#   contrast          effect.size    SE   df lower.CL upper.CL
# (LDM4 - Lreu6475)       -1.25 0.666 25.9    -2.62    0.118
# 
# Induction = Induced:
#   contrast          effect.size    SE   df lower.CL upper.CL
# (LDM4 - Lreu6475)       -5.27 1.083 25.9    -7.50   -3.046
# 
# sigma used for effect sizes: 38.98 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

#Supplemental Fig 2a

pigplotall<-ggplot(pigmelt[Treatment!="Krebs"], aes(x=Segment, y=OXTareafoldscaled,fill=Treatment,group=RegionTreatment)) +  #colour=Diet,fill=Experiment, shape=Rep
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+
  geom_point(aes(shape=as.factor(Pig)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.8 ),color="black" )+
  xlab("Region") +
  scale_fill_npg()+
  scale_shape_manual(values=c(21,22,25))+
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')')) +
  coord_cartesian(
    ylim = c(0,40))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.y= element_text(size=8),
    axis.title.x= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    plot.title=element_blank(), legend.position="none" 
  )

pdf(file="Supp2A.pdf",width=5,height=2)
pigplotall
dev.off()

#pig stats above  

#Supp2B
pigletplotall<-ggplot(pigletmelt, aes(x=Intestinal.Segment, y=OXTareafoldscaled,fill=Treatment,group=RegionTreatment)) +  
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+
  geom_point(aes(shape=as.factor(Sex)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.8 ),color="black" )+
   xlab("Region") +
  scale_fill_npg()+
  ggtitle("Piglet")+
  scale_shape_manual(values=c(21,22,25))+
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')')) +
  coord_cartesian(
    ylim = c(0,55))+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.y= element_text(size=8),
    axis.title.x= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    plot.title=element_blank(), legend.position="none" 
  )

pdf(file="Supp2B.pdf",width=4,height=2)
pigletplotall
dev.off()

#piglet stats above

#Supp2C
mousefemale<-ggplot(subset(mouse, mouse$Sex=="F"), aes(x=Segment, y=OXTareafoldscaled,fill=Treatment,group=SegmentTreatment)) + 
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black",shape=21)+
  scale_y_continuous(limits=c(0,75))+
  scale_x_discrete(labels=c("Stomach","Small intestine", "Cecum", "Large intestine"))+
  scale_fill_npg()+
  xlab("Segment") +
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')')) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.y= element_text(size=8),
    axis.title.x= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    plot.title=element_blank(), legend.position="none" 
  )

pdf(file="mousesecretionfemale.pdf",width=5,height=3)
mousefemale
dev.off()

#Supp2D
#plot just males
mousemale<-ggplot(subset(mouse, mouse$Sex=="M"), aes(x=Segment, y=OXTareafoldscaled,fill=Treatment,group=SegmentTreatment)) +  #colour=Diet,fill=Experiment, shape=Rep
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black",shape=21)+
  scale_y_continuous(limits=c(0,75))+
  scale_fill_npg()+
  scale_x_discrete(labels=c("Stomach","Small intestine", "Cecum", "Large intestine"))+
  xlab("Segment") +
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')')) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #axis.title.x= element_blank(), #element_text(size=1, colour="white"),
    axis.title.y= element_text(size=8),
    axis.title.x= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    plot.title=element_blank(), legend.position="none" 
  )

pdf(file="mousesecretionmale.pdf",width=5,height=3)
mousemale
dev.off()

#mousestats full
mousetest<-lm(OXTareafoldscaled~Treatment*Segment*Sex ,data=mouse)
summary(mousetest)
# Call:
#   lm(formula = OXTareafoldscaled ~ Treatment * Segment * Sex, data = mouse)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -13.3103  -2.4524  -0.4421   2.3868  30.6513 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      12.1423     2.6141   4.645 1.32e-05 ***
#   TreatmentLreu                    14.2526     3.6969   3.855 0.000232 ***
#   SegmentSI                       -10.3518     3.6969  -2.800 0.006402 ** 
#   SegmentCecum                     -2.6223     3.6969  -0.709 0.480186    
# SegmentLI                        -0.1715     3.6969  -0.046 0.963121    
# SexM                             11.8782     3.6969   3.213 0.001894 ** 
#   TreatmentLreu:SegmentSI         -14.4789     5.2282  -2.769 0.006979 ** 
#   TreatmentLreu:SegmentCecum       -5.3868     5.2282  -1.030 0.305961    
# TreatmentLreu:SegmentLI          -9.8166     5.2282  -1.878 0.064078 .  
# TreatmentLreu:SexM               -6.2859     5.2282  -1.202 0.232794    
# SegmentSI:SexM                  -12.3990     5.2282  -2.372 0.020117 *  
#   SegmentCecum:SexM               -15.3893     5.2282  -2.944 0.004246 ** 
#   SegmentLI:SexM                  -11.2042     5.2282  -2.143 0.035150 *  
#   TreatmentLreu:SegmentSI:SexM      7.7393     7.3938   1.047 0.298381    
# TreatmentLreu:SegmentCecum:SexM   4.8930     7.3938   0.662 0.510021    
# TreatmentLreu:SegmentLI:SexM      5.3656     7.3938   0.726 0.470150    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.403 on 80 degrees of freedom
# Multiple R-squared:  0.7006,	Adjusted R-squared:  0.6445 
# F-statistic: 12.48 on 15 and 80 DF,  p-value: 2.932e-15


mouseemmeans<-emmeans(mousetest,pairwise~Treatment | Segment & Sex,adjust="BH")
# $emmeans
# Segment = Stomach, Sex = F:
# Treatment emmean   SE df lower.CL upper.CL
# LDM4       12.14 2.61 80    6.940    17.34
# Lreu       26.39 2.61 80   21.193    31.60
# 
# Segment = SI, Sex = F:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4        1.79 2.61 80   -3.412     6.99
# Lreu        1.56 2.61 80   -3.638     6.77
# 
# Segment = Cecum, Sex = F:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4        9.52 2.61 80    4.318    14.72
# Lreu       18.39 2.61 80   13.184    23.59
# 
# Segment = LI, Sex = F:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4       11.97 2.61 80    6.769    17.17
# Lreu       16.41 2.61 80   11.205    21.61
# 
# Segment = Stomach, Sex = M:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4       24.02 2.61 80   18.818    29.22
# Lreu       31.99 2.61 80   26.785    37.19
# 
# Segment = SI, Sex = M:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4        1.27 2.61 80   -3.933     6.47
# Lreu        2.50 2.61 80   -2.706     7.70
# 
# Segment = Cecum, Sex = M:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4        6.01 2.61 80    0.807    11.21
# Lreu       13.48 2.61 80    8.280    18.68
# 
# Segment = LI, Sex = M:
#   Treatment emmean   SE df lower.CL upper.CL
# LDM4       12.64 2.61 80    7.443    17.85
# Lreu       16.16 2.61 80   10.958    21.36
# 
# Confidence level used: 0.95 
# 
# $contrasts
# Segment = Stomach, Sex = F:
#   contrast    estimate  SE df t.ratio p.value
# LDM4 - Lreu  -14.253 3.7 80  -3.855  0.0002
# 
# Segment = SI, Sex = F:
#   contrast    estimate  SE df t.ratio p.value
# LDM4 - Lreu    0.226 3.7 80   0.061  0.9513
# 
# Segment = Cecum, Sex = F:
#   contrast    estimate  SE df t.ratio p.value
# LDM4 - Lreu   -8.866 3.7 80  -2.398  0.0188
# 
# Segment = LI, Sex = F:
#   contrast    estimate  SE df t.ratio p.value
# LDM4 - Lreu   -4.436 3.7 80  -1.200  0.2337
# 
# Segment = Stomach, Sex = M:
#   contrast    estimate  SE df t.ratio p.value
# LDM4 - Lreu   -7.967 3.7 80  -2.155  0.0342
# 
# Segment = SI, Sex = M:
#   contrast    estimate  SE df t.ratio p.value
# LDM4 - Lreu   -1.227 3.7 80  -0.332  0.7408
# 
# Segment = Cecum, Sex = M:
#   contrast    estimate  SE df t.ratio p.value
# LDM4 - Lreu   -7.473 3.7 80  -2.021  0.0466
# 
# Segment = LI, Sex = M:
#   contrast    estimate  SE df t.ratio p.value
# LDM4 - Lreu   -3.516 3.7 80  -0.951  0.3445
eff_size(mouseemmeans, sigma=sigma(mousetest),edf=80)
# Segment = Stomach, Sex = F:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -2.2258 0.604 80    -3.43 -1.02469
# 
# Segment = SI, Sex = F:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)      0.0354 0.577 80    -1.11  1.18433
# 
# Segment = Cecum, Sex = F:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -1.3846 0.588 80    -2.55 -0.21514
# 
# Segment = LI, Sex = F:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -0.6928 0.580 80    -1.85  0.46136
# 
# Segment = Stomach, Sex = M:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -1.2442 0.586 80    -2.41 -0.07864
# 
# Segment = SI, Sex = M:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -0.1916 0.578 80    -1.34  0.95774
# 
# Segment = Cecum, Sex = M:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -1.1670 0.585 80    -2.33 -0.00351
# 
# Segment = LI, Sex = M:
#   contrast      effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu)     -0.5490 0.579 80    -1.70  0.60317
# 
# sigma used for effect sizes: 6.403 
# Confidence level used: 0.95 


#Supp 2E
luminexothers<-read.delim(file="FigS2e_Luminexothers.txt",header=TRUE,sep="\t",check.names=FALSE)
luminexothersmelt<-luminexothers
luminexothersmelt$Treatment<-factor(luminexothersmelt$Treatment,c("LDM4","Bsubtilis","Nissle"))

shapevalues2<-c(21,22,23)
names(shapevalues2)<-levels(as.factor(luminexothersmelt$Replicate))

pdf(file="LuminexOthers.pdf",width=1.6,height=2)
ggplot(subset(luminexothersmelt, luminexothersmelt$Induction=="Induced"), aes(x=Treatment,y=OxytocinCorrected,fill=Treatment)) +
  geom_boxplot(size=1,outlier.shape = NA,alpha = 0.4,colour="black")+  #, lty=FAMI
  geom_point(aes(shape=as.factor(Replicate)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.6 ),color="black" )+
#  geom_hline(yintercept = 15,linetype="dashed")+
  ylab(bquote('Oxytocin (pg/ml)')) +
  coord_cartesian(
    ylim = c(0,	200))+
  scale_fill_manual(values=c("#E64B35FF","#00A087FF", "#3C5488FF"))+
  scale_shape_manual(values=shapevalues2)+
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
    axis.text.x= element_text(size=7, colour = "black"),# angle=45,vjust=0.7),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )
dev.off()




