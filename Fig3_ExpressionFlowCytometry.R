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

data<-read.delim(file="Fig3a_qPCRdata.txt",header=TRUE,sep="\t")

mdt<-as.data.table(data)

mdt2<-mdt[is.na(SampleName)==FALSE]
mdt2$treatmentname<-paste(mdt2$Sample,mdt2$Primer,sep="_")
mdt2$CT<-gsub("Undetermined",NA, mdt2$CT)
mdt2$CT<-as.numeric(mdt2$CT)
mdts<-as.data.table(summarySE(data=mdt2,measurevar = "CT",groupvars = c("treatmentname","Primer","Samplegroup", "Batch")))


mdts$Sample = factor(mdts$Sample,c("Diff","NoRTDiff","NoDNA" ))
#need to fix sample names
mdts$Samplegroup = factor(mdts$Samplegroup,c("Diff","nortDiff","Ind","nortInd", "blank" ))


p <- ggplot(mdts[Primer =="GAPDH"], aes(x=Samplegroup, y=CT,fill=Samplegroup) ) + # data & aesthetic mapping
  geom_bar(stat="identity", position=position_dodge2(), colour="black") + # bars represent average
  labs(y="CT") +
  facet_wrap(.~Primer,ncol=3,scales="free")

ctcheck<-p +  theme(plot.title = element_text(size = 16,hjust = 0.5)) + # theme(aspect.ratio=1) +
  theme_classic() +
  theme(
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y= element_text(size=16, colour="black"),
    axis.text.x= element_text(size = 14, colour="black"),
    axis.text.y= element_text(size=16, colour="black"),#legend.text=element_text(size=16)
    legend.position="none")

ctcheck

PGAPDH<-mdts[Primer =="GAPDH"]
names(PGAPDH)[6]<-"CT.GAPDH"

mdts3<-merge(mdts,PGAPDH,by=c("Samplegroup","Batch"),all.x=TRUE)

mdts3$CTGAPDHnorm<-2^(mdts3$CT.GAPDH-mdts3$CT)

mdts3$Primer.x<-factor(mdts3$Primer.x,c("GAPDH","CHGA","OXT"))

jitterwidth<-.75
strokewidth<-0.5
inductioncolors<-c("slateblue1","#F0027F","goldenrod1")

mdts3$CTGAPDHnorm[3]=0.00001
mdts3$CTGAPDHnorm[9]=0.00001


PGAPDHplot<-ggplot(mdts3[Primer.x!="GAPDH"&(Samplegroup=="Diff"|Samplegroup=="Ind")], 
                   aes(x=Samplegroup, y=CTGAPDHnorm, fill=Samplegroup,col=Samplegroup) ) + # data & aesthetic mapping
  geom_point(size=1.5,pch=21,stroke=strokewidth, position=position_jitterdodge(jitter.width = jitterwidth-.25),
             color="black") + 
  geom_boxplot(size=0.6,fatten=1, outlier.color = NA,alpha=0.2,color="black") +  #, lty=FAMI
  labs(y="GAPDH Normalized CN") +
  coord_trans(y="log10")+
  facet_wrap(.~Primer.x,ncol=2,scales="fixed")+
  scale_fill_manual(values=inductioncolors)+
  scale_x_discrete(labels=c("Diff" = "Differentiated", 
                            "Ind"="Induced"))+
  scale_y_continuous(breaks=c(0.0001, 0.01, 0.5, 10,30), expand = expansion(mult = c(0.05, 0.3)))+
  theme_bw() + theme(
    plot.title = element_text(size = 7,hjust = 0.5),
    legend.position="none",
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=7, colour="black"),
    axis.title.x= element_text(size=7, colour="black"),
    axis.text.x= element_text(size = 6, colour="black", angle=45,vjust=0.7),
    axis.text.y= element_text(size=6, colour="black"),
    strip.text = element_text(size=7),
    panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background =element_blank())#legend.text=element_text(size=16))


pdf(file="qPCRDiffInduced.pdf",width=2,height=2.2)
PGAPDHplot
dev.off()

dataset<-mdts3[Primer.x!="GAPDH"&(Samplegroup=="Diff"|Samplegroup=="Ind")]
linear<-lm(CTGAPDHnorm~Samplegroup,data=dataset[Primer.x=="CHGA"])
summary(linear)
# Call:
#   lm(formula = CTGAPDHnorm ~ Samplegroup, data = dataset[Primer.x == 
#                                                            "CHGA"])
# 
# Residuals:
#   1        2        3        4        5        6 
# -0.01922 -0.01960  0.03882  0.57906 -0.36579 -0.21327 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     0.02124    0.20754   0.102    0.923  
# SamplegroupInd  0.96742    0.29350   3.296    0.030 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3595 on 4 degrees of freedom
# Multiple R-squared:  0.7309,	Adjusted R-squared:  0.6636 
# F-statistic: 10.86 on 1 and 4 DF,  p-value: 0.03004

emmeans(linear,pairwise~Samplegroup,adjust="none")
# $emmeans
# Samplegroup   emmean    SE df lower.CL upper.CL
# Undiff      0.000427 0.169  6   -0.414    0.415
# Diff        0.021236 0.169  6   -0.393    0.436
# Ind         0.988660 0.169  6    0.574    1.403
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast      estimate   SE df t.ratio p.value
# Undiff - Diff  -0.0208 0.24  6  -0.087  0.9336
# Undiff - Ind   -0.9882 0.24  6  -4.124  0.0062
# Diff - Ind     -0.9674 0.24  6  -4.037  0.0068

OXTlinear<-lm(CTGAPDHnorm~Samplegroup,data=dataset[Primer.x=="OXT"])
summary(OXTlinear)
# Call:
#   lm(formula = CTGAPDHnorm ~ Samplegroup, data = dataset[Primer.x == 
#                                                            "OXT"])
# 
# Residuals:
#   1          2          3          4          5          6 
# 1.081e-03 -2.371e-04 -8.434e-04 -1.057e-05  3.409e-05 -2.352e-05 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     0.0011037  0.0004018   2.747   0.0515 .
# SamplegroupInd -0.0009998  0.0005682  -1.760   0.1533  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0006959 on 4 degrees of freedom
# Multiple R-squared:  0.4363,	Adjusted R-squared:  0.2954 
# F-statistic: 3.096 on 1 and 4 DF,  p-value: 0.1533

OXTemmeans<-emmeans(OXTlinear,pairwise~Samplegroup,adjust="none")
# $emmeans
# Samplegroup   emmean       SE df  lower.CL upper.CL
# Diff        0.001104 0.000402  4 -1.18e-05  0.00222
# Ind         0.000104 0.000402  4 -1.01e-03  0.00122
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast   estimate       SE df t.ratio p.value
# Diff - Ind    0.001 0.000568  4   1.760  0.1533

eff_size(OXTemmeans,sigma = sigma(OXTlinear), edf=4)
# contrast     effect.size    SE df lower.CL upper.CL
# (Diff - Ind)        1.44 0.962  4    -1.23     4.11
# 
# sigma used for effect sizes: 0.0006959 
# Confidence level used: 0.95

CHGAlinear<-lm(CTGAPDHnorm~Samplegroup,data=dataset[Primer.x=="CHGA"])
summary(CHGAlinear)
# Call:
#   lm(formula = CTGAPDHnorm ~ Samplegroup, data = dataset[Primer.x == 
#                                                            "CHGA"])
# 
# Residuals:
#   1        2        3        4        5        6 
# -0.01922 -0.01960  0.03882  0.57906 -0.36579 -0.21327 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     0.02124    0.20754   0.102    0.923  
# SamplegroupInd  0.96742    0.29350   3.296    0.030 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3595 on 4 degrees of freedom
# Multiple R-squared:  0.7309,	Adjusted R-squared:  0.6636 
# F-statistic: 10.86 on 1 and 4 DF,  p-value: 0.03004

CHGAemmeans<-emmeans(CHGAlinear,pairwise~Samplegroup,adjust="none")
# $emmeans
# Samplegroup emmean    SE df lower.CL upper.CL
# Diff        0.0212 0.208  4   -0.555    0.597
# Ind         0.9887 0.208  4    0.412    1.565
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast   estimate    SE df t.ratio p.value
# Diff - Ind   -0.967 0.294  4  -3.296  0.0300

eff_size(CHGAemmeans,sigma=sigma(CHGAlinear), edf=4)
# contrast     effect.size   SE df lower.CL upper.CL
# (Diff - Ind)       -2.69 1.25  4    -6.17     0.79
# 
# sigma used for effect sizes: 0.3595 
# Confidence level used: 0.95

#UnDiif Diff Induced for qPCR FACS

FACSdata<-read.delim(file="Fig3b_FlowData.txt",sep="\t",header=TRUE, check.names = FALSE)

FACSdata$CellType<-factor(FACSdata$CellType,c("Diff", "Induced"))


inductioncolors<-c("slateblue1","#F0027F","goldenrod1")
FACsplot<-ggplot(FACSdata, aes(x=CellType, y=Percent,fill=CellType,col=CellType) ) + # data & aesthetic mapping
  geom_point(size=1.5,pch=21,stroke=strokewidth, position=position_jitterdodge(jitter.width = jitterwidth-.25),color="black") + 
  geom_boxplot(size=0.6,fatten=1, outlier.color = NA,alpha=0.2,color="black") +  
  labs(y="Percent stained") +
  facet_wrap(~Staining,nrow=1,scales="fixed")+
  coord_trans(y="log10")+
  scale_fill_manual(name = "FullSampleinfo",values=inductioncolors)+
  scale_y_continuous(breaks=c(0.05, 0.1, 0.5, 1,3,5), expand = expansion(mult = c(0.05, 0.3)))+
  scale_x_discrete(labels=c("Differentiated","Induced"))+
  theme_bw() + theme(
    plot.title = element_text(size = 7,hjust = 0.5),
    legend.position="none",
    axis.ticks = element_line(colour="black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.y= element_text(size=7, colour="black"),
    axis.title.x= element_text(size=7, colour="black"),
    axis.text.x= element_text(size = 6, colour="black", angle=45,vjust=0.7),
    axis.text.y= element_text(size=6, colour="black"),
    strip.text = element_text(size=7),
    panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background =element_blank())#legend.text=element_text(size=16))

pdf(file="FACSDiffInduced.pdf",width=2,height=2.2)
FACsplot
dev.off()

)

CHGAlinear<-lm(Percent~CellType,data=subset(FACSdata,FACSdata$Staining=="CHGA"))
summary(CHGAlinear)
#Call:
# lm(formula = Percent ~ CellType, data = subset(FACSdata2, FACSdata2$Staining == 
#                                                  "CHGA"))
# 
# Residuals:
#   5      6     11     12     17     18 
# 0.156  2.093 -0.142  1.063 -0.014 -3.157 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        0.224      1.137   0.197   0.8535  
# CellTypeInduced    4.573      1.608   2.843   0.0467 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.97 on 4 degrees of freedom
# Multiple R-squared:  0.6689,	Adjusted R-squared:  0.5862 
# F-statistic: 8.082 on 1 and 4 DF,  p-value: 0.04673

CHGAemmeans<-emmeans(CHGAlinear,pairwise~CellType,adjust="none")
# $emmeans
# CellType emmean   SE df lower.CL upper.CL
# Diff      0.224 1.14  4    -2.93     3.38
# Induced   4.797 1.14  4     1.64     7.95
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast       estimate   SE df t.ratio p.value
# Diff - Induced    -4.57 1.61  4  -2.843  0.0467

eff_size(CHGAemmeans,sigma = sigma(CHGAlinear), edf=4)

# contrast         effect.size   SE df lower.CL upper.CL
# (Diff - Induced)       -2.32 1.16  4    -5.54    0.893
# 
# sigma used for effect sizes: 1.97 
# Confidence level used: 0.95 

OXTlinear<-lm(Percent~CellType,data=subset(FACSdata,FACSdata$Staining=="OXT"))
summary(OXTlinear)
# Call:
#   lm(formula = Percent ~ CellType, data = subset(FACSdata2, FACSdata2$Staining == 
#                                                    "OXT"))
# 
# Residuals:
#   2         3         8         9        14        15 
# -0.006667  0.019000 -0.056667 -0.015000  0.063333 -0.004000 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)      0.15667    0.02561   6.118  0.00361 **
#   CellTypeInduced -0.10567    0.03621  -2.918  0.04334 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04435 on 4 degrees of freedom
# Multiple R-squared:  0.6804,	Adjusted R-squared:  0.6004 
# F-statistic: 8.514 on 1 and 4 DF,  p-value: 0.04334

OXTemmeans<-emmeans(OXTlinear,pairwise~CellType,adjust="none")
# $emmeans
# CellType emmean     SE df lower.CL upper.CL
# Diff      0.157 0.0256  4   0.0856    0.228
# Induced   0.051 0.0256  4  -0.0201    0.122
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast       estimate     SE df t.ratio p.value
# Diff - Induced    0.106 0.0362  4   2.918  0.0433

eff_size(OXTemmeans, sigma = sigma(OXTlinear), edf=4)
# contrast         effect.size   SE df lower.CL upper.CL
# (Diff - Induced)        2.38 1.17  4   -0.875     5.64
# 
# sigma used for effect sizes: 0.04435 
# Confidence level used: 0.95
