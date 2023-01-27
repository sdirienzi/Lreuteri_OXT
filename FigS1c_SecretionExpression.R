
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
mdt<-as.data.table(read.delim(file="qPCRdataD103J2NGN3.txt",header=TRUE,sep="\t",na.strings=""))
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
mdts3add<-read.delim(file="mdts3extras.txt",sep="\t",header=TRUE)
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
linear<-lm(CTGAPDHnorm~CellType+Sample,data=mdts4[Primer.x=="LGR5"])
summary(linear)
emmeans(linear,pairwise~CellType+Sample,adjust="none")
# $emmeans
# CellType Sample    emmean       SE df  lower.CL upper.CL
# D103     Undiff  0.004407 0.000866  9  0.002446  0.00637
# NGN3     Undiff  0.002493 0.000866  9  0.000533  0.00445
# D103     Diff    0.000967 0.000866  9 -0.000993  0.00293
# NGN3     Diff   -0.000947 0.000866  9 -0.002907  0.00101
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                  estimate      SE df t.ratio p.value
# D103 Undiff - NGN3 Undiff  0.00191 0.00100  9   1.912  0.0881
# D103 Undiff - D103 Diff    0.00344 0.00100  9   3.438  0.0074
# D103 Undiff - NGN3 Diff    0.00535 0.00141  9   3.783  0.0043
# NGN3 Undiff - D103 Diff    0.00153 0.00141  9   1.079  0.3087
# NGN3 Undiff - NGN3 Diff    0.00344 0.00100  9   3.438  0.0074
# D103 Diff - NGN3 Diff      0.00191 0.00100  9   1.912  0.0881

SIlinear<-lm(CTGAPDHnorm~CellType+Sample,data=mdts4[Primer.x=="SI"])
SIemmeans<-emmeans(linear,pairwise~CellType+Sample,adjust="none")

# $emmeans
# CellType Sample  emmean     SE df lower.CL upper.CL
# D103     Undiff  0.0112 0.0226  9 -0.04000   0.0624
# NGN3     Undiff -0.0105 0.0226  9 -0.06166   0.0407
# D103     Diff    0.0812 0.0226  9  0.02998   0.1323
# NGN3     Diff    0.0595 0.0226  9  0.00832   0.1107
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                  estimate     SE df t.ratio p.value
# D103 Undiff - NGN3 Undiff   0.0217 0.0261  9   0.829  0.4285
# D103 Undiff - D103 Diff    -0.0700 0.0261  9  -2.679  0.0252
# D103 Undiff - NGN3 Diff    -0.0483 0.0369  9  -1.308  0.2233
# NGN3 Undiff - D103 Diff    -0.0916 0.0369  9  -2.481  0.0350
# NGN3 Undiff - NGN3 Diff    -0.0700 0.0261  9  -2.679  0.0252
# D103 Diff - NGN3 Diff       0.0217 0.0261  9   0.829  0.4285
eff_size(SIemmeans,sigma = sigma(SIlinear),edf=9)
# contrast                    effect.size    SE df lower.CL upper.CL
# (D103 Undiff - NGN3 Undiff)       0.479 0.588  9   -0.852  1.80944
# (D103 Undiff - D103 Diff)        -1.547 0.683  9   -3.091 -0.00206
# (D103 Undiff - NGN3 Diff)        -1.068 0.854  9   -3.001  0.86482
# (NGN3 Undiff - D103 Diff)        -2.025 0.946  9   -4.165  0.11421
# (NGN3 Undiff - NGN3 Diff)        -1.547 0.683  9   -3.091 -0.00206
# (D103 Diff - NGN3 Diff)           0.479 0.588  9   -0.852  1.80944
# 
# sigma used for effect sizes: 0.04525 
# Confidence level used: 0.95

OXTlinear<-lm(CTGAPDHnorm~CellType+Sample,data=mdts4[Primer.x=="OXT"])
OXTemmeans<-emmeans(linear,pairwise~CellType+Sample,adjust="none")
# $emmeans
# CellType Sample    emmean       SE df  lower.CL upper.CL
# D103     Undiff  1.57e-04 0.000118  9 -1.09e-04 0.000423
# NGN3     Undiff -6.63e-05 0.000118  9 -3.32e-04 0.000200
# D103     Diff    5.25e-04 0.000118  9  2.59e-04 0.000791
# NGN3     Diff    3.01e-04 0.000118  9  3.52e-05 0.000567
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                   estimate       SE df t.ratio p.value
# D103 Undiff - NGN3 Undiff  0.000223 0.000136  9   1.644  0.1345
# D103 Undiff - D103 Diff   -0.000368 0.000136  9  -2.706  0.0242
# D103 Undiff - NGN3 Diff   -0.000144 0.000192  9  -0.751  0.4721
# NGN3 Undiff - D103 Diff   -0.000591 0.000192  9  -3.076  0.0132
# NGN3 Undiff - NGN3 Diff   -0.000368 0.000136  9  -2.706  0.0242
# D103 Diff - NGN3 Diff      0.000223 0.000136  9   1.644  0.1345

eff_size(OXTemmeans,sigma = sigma(OXTlinear),edf=9)
# contrast                    effect.size  SE df lower.CL upper.CL
# (D103 Undiff - NGN3 Undiff)          92 113  9     -164  347.934
# (D103 Undiff - D103 Diff)          -297 131  9     -594   -0.396
# (D103 Undiff - NGN3 Diff)          -205 164  9     -577  166.295
# (NGN3 Undiff - D103 Diff)          -389 182  9     -801   21.961
# (NGN3 Undiff - NGN3 Diff)          -297 131  9     -594   -0.396
# (D103 Diff - NGN3 Diff)              92 113  9     -164  347.934
# 
# sigma used for effect sizes: 0.0002353 
# Confidence level used: 0.95 

LGR5linear<-lm(CTGAPDHnorm~CellType+Sample,data=mdts4[Primer.x=="LGR5"])
LGR5emmeans<-emmeans(linear,pairwise~CellType+Sample,adjust="none")
# $emmeans
# CellType Sample  emmean     SE df lower.CL upper.CL
# D103     Undiff  0.0112 0.0226  9 -0.04000   0.0624
# NGN3     Undiff -0.0105 0.0226  9 -0.06166   0.0407
# D103     Diff    0.0812 0.0226  9  0.02998   0.1323
# NGN3     Diff    0.0595 0.0226  9  0.00832   0.1107
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                  estimate     SE df t.ratio p.value
# D103 Undiff - NGN3 Undiff   0.0217 0.0261  9   0.829  0.4285
# D103 Undiff - D103 Diff    -0.0700 0.0261  9  -2.679  0.0252
# D103 Undiff - NGN3 Diff    -0.0483 0.0369  9  -1.308  0.2233
# NGN3 Undiff - D103 Diff    -0.0916 0.0369  9  -2.481  0.0350
# NGN3 Undiff - NGN3 Diff    -0.0700 0.0261  9  -2.679  0.0252
# D103 Diff - NGN3 Diff       0.0217 0.0261  9   0.829  0.4285

eff_size(LGR5emmeans,sigma = sigma(LGR5linear),edf=9)
# contrast                    effect.size   SE df lower.CL upper.CL
# (D103 Undiff - NGN3 Undiff)        12.5 15.4  9    -22.2  47.2456
# (D103 Undiff - D103 Diff)         -40.4 17.8  9    -80.7  -0.0538
# (D103 Undiff - NGN3 Diff)         -27.9 22.3  9    -78.4  22.5810
# (NGN3 Undiff - D103 Diff)         -52.9 24.7  9   -108.7   2.9820
# (NGN3 Undiff - NGN3 Diff)         -40.4 17.8  9    -80.7  -0.0538
# (D103 Diff - NGN3 Diff)            12.5 15.4  9    -22.2  47.2456
# 
# sigma used for effect sizes: 0.001733 
# Confidence level used: 0.95 

