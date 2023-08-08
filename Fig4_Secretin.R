#Secretin plots for Figure 4


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



#4a: secretin causes release of oxytocin from NGN3 HIEs
HIEoxt<-read.delim(file="Fig4a_NGNSecretinOxytocin.txt",header=TRUE,sep="\t",check.names=FALSE)

#plot as boxplots with points and coloring
HIEFig4A<-ggplot(HIEoxt, aes(x=Treatment,y=Oxytocin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(shape=21,size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.4,.4,.4))+
  scale_fill_npg()+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
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
pdf(file="NGNsecretinoxytocin.pdf",width=1.3,height=2)
grid.arrange(HIEFig4A)
dev.off()

Fig4amodel<-lm(Oxytocin~Treatment,data=HIEoxt)
summary(Fig4amodel)
# Call:
#   lm(formula = Oxytocin ~ Treatment, data = HIEoxt)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -81.45 -28.39   0.00   0.00 109.84 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)            7.50      35.91   0.209  0.84146   
# TreatmentLreu6475    265.49      50.78   5.228  0.00196 **
#   TreatmentSecretin    149.77      50.78   2.949  0.02563 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 62.2 on 6 degrees of freedom
# Multiple R-squared:  0.8208,	Adjusted R-squared:  0.7611 
# F-statistic: 13.74 on 2 and 6 DF,  p-value: 0.005755
modelemmeans<-emmeans(Fig4amodel,pairwise~Treatment,adjust="BH")
# Treatment emmean   SE df lower.CL upper.CL
# LDM4         7.5 35.9  6    -80.4     95.4
# Lreu6475   273.0 35.9  6    185.1    360.8
# Secretin   157.3 35.9  6     69.4    245.1
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast             estimate   SE df t.ratio p.value
# LDM4 - Lreu6475          -265 50.8  6  -5.228  0.0059
# LDM4 - Secretin          -150 50.8  6  -2.949  0.0385
# Lreu6475 - Secretin       116 50.8  6   2.279  0.0629
# P value adjustment: BH method for 3 tests

eff_size(modelemmeans,sigma = sigma(Fig4amodel),edf=6)
# contrast               effect.size    SE df lower.CL upper.CL
# (LDM4 - Lreu6475)            -4.27 1.478  6   -7.886   -0.652
# (LDM4 - Secretin )           -2.41 1.072  6   -5.032    0.216
# (Lreu6475 - Secretin )        1.86 0.977  6   -0.531    4.252
# 
# sigma used for effect sizes: 62.2 
# Confidence level used: 0.95 

#4b secretin casuses release of oxtycoin from other HIEs
HIEalloxt<-read.delim(file="Fig4b_HIOSecretinOxytocin.txt",header=TRUE,sep="\t",check.names=FALSE)
HIEalloxt$SegmentTreatment=paste(HIEalloxt$Segment, HIEalloxt$Treatment)
HIEalloxt$Segment<-factor(HIEalloxt$Segment,c("Duodenum","Jejunum","Ileum", "AscendingColon"))


HIEallFig4B<-ggplot(HIEalloxt, aes(x=Segment,y=Oxytocin,fill=Treatment,group=SegmentTreatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(shape=21,size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.3,.3,.3))+
  scale_fill_manual(values=c("#E64B35FF", "#00A087FF"))+
  xlab("Enteroid type")+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
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

pdf(file="HIOsecretinoxytocin.pdf",width=2.6,height=2)
grid.arrange(HIEallFig4B)
dev.off()

Figbmodel<-lmer(Oxytocin~Treatment+(1|Enteroid),data=HIEalloxt,control=lmerControl(optimizer = "bobyqa"))
summary(Figbmodel)
Figbnullmodel<-lmer(Oxytocin~(1|Enteroid),data=HIEalloxt,control=lmerControl(optimizer = "bobyqa"))
anova(Figbmodel,Figbnullmodel)
# Data: HIEalloxt
# Models:
#   Figbnullmodel: Oxytocin ~ (1 | Enteroid)
# Figbmodel: Oxytocin ~ Treatment + (1 | Enteroid)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# Figbnullmodel    3 270.95 274.48 -132.47   264.95                         
# Figbmodel        4 219.59 224.30 -105.80   211.59 53.357  1  2.781e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
performance::r2(Figbmodel)
# Conditional R2: 0.899
# Marginal R2: 0.883

HIEemmeans<-emmeans(Figbmodel,pairwise~Treatment,adjust="BH")
HIEemmeans
# Treatment emmean   SE  df lower.CL upper.CL
# Krebs        7.5 6.89 6.6    -8.99       24
# Secretin   121.5 6.89 6.6   104.96      138
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast          estimate   SE df t.ratio p.value
# Krebs - Secretin      -114 8.03 19 -14.190  <.0001
# 
# Degrees-of-freedom method: kenward-roger

eff_size(HIEemmeans,sigma = sigma(Figbmodel),edf=19)
# contrast            effect.size   SE df lower.CL upper.CL
# (Krebs - Secretin )       -5.79 1.02 19    -7.94    -3.65
# 
# sigma used for effect sizes: 19.67 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

#4c scretin causes release of oxytocin from lifegift
LGoxt<-read.delim(file="Fig4c_LGSecretinOxytocin.txt",header=TRUE,sep="\t",check.names=FALSE)
LGoxt$Normalized<-LGoxt$Oxytocin/LGoxt$Area
shapevalues<-c(23,24,25)
names(shapevalues)<-c("LG3","LG4","LG5")

LGFig4C<-ggplot(LGoxt, aes(x=Treatment,y=Normalized,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(aes(shape=LifeGift),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')'))+
  scale_y_continuous(expand=c(0.05,.3,.3,.3))+
  xlab("Treatment")+
  scale_shape_manual(values=shapevalues)+
  scale_fill_manual(values=c("#E64B35FF",  "#4DBBD5FF","#00A087FF"))+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
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

pdf(file="~LifeGiftsecretinoxytocin.pdf",width=1.65,height=2)
grid.arrange(LGFig4C)
dev.off()

FigCmodel<-lmer(Normalized~Treatment+(1|LifeGift),data=LGoxt,REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
FigCnull<-lmer(Normalized~(1|LifeGift),data=LGoxt,REML=FALSE,control=lmerControl(optimizer = "bobyqa"))
anova(FigCmodel,FigCnull)
# Data: LGoxt
# Models:
#   FigCnull: Normalized ~ (1 | LifeGift)
# FigCmodel: Normalized ~ Treatment + (1 | LifeGift)
# npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)    
# FigCnull     3 321.80 325.69 -157.90   315.81                        
# FigCmodel    5 294.24 300.71 -142.12   284.24 31.57  2  1.396e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
performance::r2(FigCmodel)
# Conditional R2: 0.734
# Marginal R2: 0.669

summary(FigCmodel)
FigCemmeans<-emmeans(FigCmodel,pairwise~Treatment,adjust="BH")
# $emmeans
# Treatment emmean   SE   df lower.CL upper.CL
# Krebs       14.1 22.3 11.1    -34.9     63.1
# Lreu6475   167.4 22.3 11.1    118.4    216.4
# Secretin   147.9 22.3 11.1     98.9    196.9
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast            estimate   SE   df t.ratio p.value
# Krebs - Lreu6475      -153.3 21.6 26.2  -7.107  <.0001
# Krebs - Secretin      -133.8 21.6 26.2  -6.204  <.0001
# Lreu6475 - Secretin     19.5 21.6 26.2   0.903  0.3749
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: BH method for 3 tests

eff_size(FigCemmeans,sigma = sigma(FigCmodel),edf=26.2)

# contrast              effect.size    SE   df lower.CL upper.CL
# (Krebs - Lreu6475)         -3.499 0.690 26.2   -4.917    -2.08
# (Krebs - Secretin)         -3.055 0.648 26.2   -4.387    -1.72
# (Lreu6475 - Secretin)       0.444 0.496 26.2   -0.575     1.46
# 
# sigma used for effect sizes: 43.81 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

#toprow
pdf(file="~Fig4ABC.pdf",width=5.5,height=2)
grid.arrange(HIEFig4A,HIEallFig4B,LGFig4C,widths=c(1,2,1.2))
dev.off()

#Fig S4 antibody block experiment
Abblock<-read.delim(file="FigS4_HIEABblocker.txt",header=TRUE,sep="\t",check.names=FALSE)
Abblock$Treatment<-factor(Abblock$Treatment,c("LDM4","LDM4 +Ab","Lreu6475","Lreu6475 +Ab"))

HIEFig4D<-ggplot(Abblock, aes(x=Treatment,y=Oxytocin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  #, lty=FAMI
  geom_point(shape=21,size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.6,.6,.6))+
  scale_fill_manual(values=c("#E64B35FF","#F39B7FFF","#4DBBD5FF","#8491B4FF"))+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
   axis.title.x= element_blank(),
    axis.title.y= element_text(size=8),
    axis.text.x=element_blank(),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )
pdf(file="~HIEABblock_Fig4D.pdf",width=1.6,height=1.6)
grid.arrange(HIEFig4D)
dev.off()

FigS4model<-lm(Oxytocin~Treatment,data=Abblock)
summary(FigS4model)
# Call:
#   lm(formula = Oxytocin ~ Treatment, data = Abblock)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -156.850   -5.656    0.000    7.927  135.692 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             13.156     42.693   0.308    0.766    
# TreatmentLDM4 +Ab       -5.656     60.377  -0.094    0.928    
# TreatmentLreu6475      468.849     60.377   7.765 5.41e-05 ***
#   TreatmentLreu6475 +Ab  111.494     60.377   1.847    0.102    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 73.95 on 8 degrees of freedom
# Multiple R-squared:  0.9112,	Adjusted R-squared:  0.878 
# F-statistic: 27.38 on 3 and 8 DF,  p-value: 0.0001472

FigS4emmeans<-emmeans(FigS4model,pairwise~Treatment,adjust="BH")
# $emmeans
# Treatment    emmean   SE df lower.CL upper.CL
# LDM4           13.2 42.7  8    -85.3      112
# LDM4 +Ab        7.5 42.7  8    -91.0      106
# Lreu6475      482.0 42.7  8    383.6      580
# Lreu6475 +Ab  124.6 42.7  8     26.2      223
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                    estimate   SE df t.ratio p.value
# LDM4 - (LDM4 +Ab)               5.66 60.4  8   0.094  0.9277
# LDM4 - Lreu6475              -468.85 60.4  8  -7.765  0.0002
# LDM4 - (Lreu6475 +Ab)        -111.49 60.4  8  -1.847  0.1224
# (LDM4 +Ab) - Lreu6475        -474.50 60.4  8  -7.859  0.0002
# (LDM4 +Ab) - (Lreu6475 +Ab)  -117.15 60.4  8  -1.940  0.1224
# Lreu6475 - (Lreu6475 +Ab)     357.35 60.4  8   5.919  0.0007
# P value adjustment: BH method for 6 tests 

eff_size(FigDemmeans,sigma = sigma(FigDmodel),edf=8)
# contrast                      effect.size    SE df lower.CL upper.CL
# (LDM4 - (LDM4 +Ab))                0.0765 0.817  8    -1.81    1.960
# (LDM4 - Lreu6475)                 -6.3404 1.783  8   -10.45   -2.229
# (LDM4 - (Lreu6475 +Ab))           -1.5078 0.899  8    -3.58    0.566
# ((LDM4 +Ab) - Lreu6475)           -6.4169 1.800  8   -10.57   -2.266
# ((LDM4 +Ab) - (Lreu6475 +Ab))     -1.5843 0.907  8    -3.68    0.508
# (Lreu6475 - (Lreu6475 +Ab))        4.8326 1.458  8     1.47    8.195
# 
# sigma used for effect sizes: 73.95 
# Confidence level used: 0.95 

#4D 5-27 blockng experiment on NGN3s
block527data<-read.delim(file="Fig4d_527inhibitor.txt",header=TRUE,sep="\t",check.names = FALSE)

block527data$Treatment<-factor(block527data$Treatment,c("Krebs","Inhibitor","Secretin","Secretin+ Inhibitor"))

shapevalues<-c(23,24,25,22)
names(shapevalues)<-c("1","2","3","4")

HIEFig4D<-ggplot(block527data, aes(x=Treatment,y=Oxytocin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+ 
  geom_point(aes(shape=as.factor(Experiment)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.6,.6,.6))+
  scale_fill_manual(values=c("#E64B35FF","#F39B7FFF","#00A087FF","#91D1C2FF"))+
  scale_shape_manual(values=shapevalues)+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
      axis.title.x= element_blank(),
    axis.title.y= element_text(size=8),
    axis.text.x=element_blank(),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )
pdf(file="HIE527block_Fig4D.pdf",width=1.8,height=1.2)
grid.arrange(HIEFig4D)
dev.off()

FigDmodel<-lmer(Oxytocin~Treatment+(1|Experiment),data=block527data, control=lmerControl(optimizer = "bobyqa"))
FigDnull<-lmer(Oxytocin~(1|Experiment),data=block527data, control=lmerControl(optimizer = "bobyqa"))
anova(FigDmodel,FigDnull)
# Data: block527data
# Models:
#   FigDnull: Oxytocin ~ (1 | Experiment)
# FigDmodel: Oxytocin ~ Treatment + (1 | Experiment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# FigDnull     3 568.35 574.21 -281.18   562.35                         
# FigDmodel    6 435.63 447.34 -211.81   423.63 138.72  3  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
performance::r2(FigDmodel)
# Conditional R2: 0.931
# Marginal R2: 0.925

summary(FigDmodel)
FigDemmeans<-emmeans(FigDmodel,pairwise~Treatment,adjust="BH")
# Treatment           emmean   SE   df lower.CL upper.CL
# Krebs                 7.77 4.52 15.4    -1.83     17.4
# Inhibitor             7.77 4.52 15.4    -1.83     17.4
# Secretin            134.05 4.52 15.4   124.45    143.7
# Secretin+ Inhibitor  33.60 4.52 15.4    24.00     43.2
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                          estimate   SE df t.ratio p.value
# Krebs - Inhibitor                      0.0 5.65 45   0.000  1.0000
# Krebs - Secretin                    -126.3 5.65 45 -22.367  <.0001
# Krebs - (Secretin+ Inhibitor)        -25.8 5.65 45  -4.575  <.0001
# Inhibitor - Secretin                -126.3 5.65 45 -22.367  <.0001
# Inhibitor - (Secretin+ Inhibitor)    -25.8 5.65 45  -4.575  <.0001
# Secretin - (Secretin+ Inhibitor)     100.4 5.65 45  17.791  <.0001


eff_size(FigDemmeans,sigma = sigma(FigDmodel), edf=45)
# contrast                            effect.size    SE df lower.CL upper.CL
# (Krebs - Inhibitor)                        0.00 0.392 45    -0.79    0.790
# (Krebs - Secretin)                        -8.77 1.004 45   -10.80   -6.750
# (Krebs - (Secretin+ Inhibitor))           -1.79 0.435 45    -2.67   -0.918
# (Inhibitor - Secretin)                    -8.77 1.004 45   -10.80   -6.750
# (Inhibitor - (Secretin+ Inhibitor))       -1.79 0.435 45    -2.67   -0.918
# (Secretin - (Secretin+ Inhibitor))         6.98 0.834 45     5.30    8.657
# 
# sigma used for effect sizes: 14.39 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95

#4e 5-27 on HIE with Lreu
block527Lreudata<-read.delim(file="Fig4e_secretininhibitor.txt",header=TRUE,sep="\t",check.names = FALSE)

block527Lreudata$Treatment<-factor(block527Lreudata$Treatment,c("Krebs","Secretin_inhibitor","6475","6475_5-27"))

shapevalues<-c(23,24,25)
names(shapevalues)<-c("1","2","3")

HIEFig4E<-ggplot(block527Lreudata, aes(x=Treatment,y=Oxytocin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(aes(shape=as.factor(Experiment)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.6,.6,.6))+
  scale_fill_manual(values=c("#E64B35FF","#F39B7FFF","#4DBBD5FF","#8491B4FF"))+
  scale_shape_manual(values=shapevalues)+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_blank(),
    axis.title.y= element_text(size=8),
    axis.text.x=element_blank(),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )
pdf(file="block527Lreudata_FigEF.pdf",width=1.8,height=1.2)
grid.arrange(HIEFig4E)
dev.off()

FigEmodel<-lmer(Oxytocin~Treatment+(1|Experiment),data=block527Lreudata, control=lmerControl(optimizer = "bobyqa"))
FigEnull<-lmer(Oxytocin~(1|Experiment),data=block527Lreudata, control=lmerControl(optimizer = "bobyqa"))
anova(FigEmodel,FigEnull)
# Data: block527Lreudata
# Models:
#   FigFnull: Oxytocin ~ (1 | Experiment)
# FigFmodel: Oxytocin ~ Treatment + (1 | Experiment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# FigFnull     3 564.85 570.47 -279.43   558.85                         
# FigFmodel    6 496.53 507.75 -242.26   484.53 74.328  3  5.048e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
performance::r2(FigEmodel)
# Conditional R2: 0.794
# Marginal R2: 0.763

summary(FigEmodel)
FigEemmeans<-emmeans(FigEmodel,pairwise~Treatment,adjust="BH")
# $emmeans
# Treatment          emmean   SE   df lower.CL upper.CL
# Krebs                9.20 13.8 6.82 -23.7156     42.1
# Secretin_inhibitor   9.75 13.8 6.82 -23.1682     42.7
# 6475               182.43 13.8 6.82 149.5176    215.4
# 6475_5-27           32.85 13.8 6.82  -0.0673     65.8
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                         estimate   SE df t.ratio p.value
# Krebs - Secretin_inhibitor         -0.547 15.5 42  -0.035  0.9719
# Krebs - 6475                     -173.233 15.5 42 -11.197  <.0001
# Krebs - (6475_5-27)               -23.648 15.5 42  -1.528  0.1715
# Secretin_inhibitor - 6475        -172.686 15.5 42 -11.161  <.0001
# Secretin_inhibitor - (6475_5-27)  -23.101 15.5 42  -1.493  0.1715
# 6475 - (6475_5-27)                149.585 15.5 42   9.668  <.0001
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: BH method for 6 tests
eff_size(FigEemmeans, sigma = sigma(FigEmodel), edf=42)
# contrast                           effect.size    SE df lower.CL upper.CL
# (Krebs - Secretin_inhibitor)           -0.0144 0.408 42   -0.838    0.809
# (Krebs - 6475)                         -4.5710 0.645 42   -5.872   -3.270
# (Krebs - (6475_5-27))                  -0.6240 0.414 42   -1.459    0.211
# (Secretin_inhibitor - 6475)            -4.5566 0.643 42   -5.855   -3.258
# (Secretin_inhibitor - (6475_5-27))     -0.6096 0.414 42   -1.444    0.225
# (6475 - (6475_5-27))                    3.9470 0.593 42    2.749    5.145
# 
# sigma used for effect sizes: 37.9 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

#4G 5-27 on LifeGift
LGblock527data<-read.delim(file="Fig4f_LifeGift527.txt",header=TRUE,sep="\t",check.names = FALSE)
LGblock527data$Normalized<-LGblock527data$Oxytocin/LGblock527data$Area

LGblock527data$Treatment<-factor(LGblock527data$Treatment,c("Krebs","Lreu6475","Secretin","6475 + 5-27"))

LGFig4f<-ggplot(LGblock527data, aes(x=Treatment,y=Normalized,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  #, lty=FAMI
  geom_point(shape=25,size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml/'*cm^2*')'))+
  scale_y_continuous(expand=c(0.05,.42,.42,.42))+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#8491B4FF"))+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
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
pdf(file="LG527block_Fig4F.pdf",width=2,height=2)
grid.arrange(LGFig4f)
dev.off()

FigFmodel<-lm(Normalized~Treatment,data=LGblock527data)
summary(FigFmodel)
# Call:
#   lm(formula = Normalized ~ Treatment, data = LGblock527data)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -136.367   -6.105   -0.304    0.523  193.134 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)             16.70      49.76   0.336   0.7458  
# TreatmentLreu6475      162.15      70.38   2.304   0.0501 .
# TreatmentSecretin      198.59      70.38   2.822   0.0224 *
#   Treatment6475 + 5-27   -12.95      70.38  -0.184   0.8586  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 86.19 on 8 degrees of freedom
# Multiple R-squared:  0.6429,	Adjusted R-squared:  0.5089 
# F-statistic:   4.8 on 3 and 8 DF,  p-value: 0.03381


FigFemmeans<-emmeans(FigFmodel,pairwise~Treatment,adjust="BH")
# $emmeans
# Treatment   emmean   SE df lower.CL upper.CL
# Krebs        16.70 49.8  8    -98.1      131
# Lreu6475    178.86 49.8  8     64.1      294
# Secretin    215.29 49.8  8    100.5      330
# 6475 + 5-27   3.75 49.8  8   -111.0      119
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                 estimate   SE df t.ratio p.value
# Krebs - Lreu6475           -162.2 70.4  8  -2.304  0.0752
# Krebs - Secretin           -198.6 70.4  8  -2.822  0.0673
# Krebs - (6475 + 5-27)        13.0 70.4  8   0.184  0.8586
# Lreu6475 - Secretin         -36.4 70.4  8  -0.518  0.7424
# Lreu6475 - (6475 + 5-27)    175.1 70.4  8   2.488  0.0752
# Secretin - (6475 + 5-27)    211.5 70.4  8   3.006  0.0673
# 
# P value adjustment: BH method for 6 tests 

eff_size(FigGemmeans,sigma(FigGmodel),edf=8)
# contrast                   effect.size    SE df lower.CL upper.CL
# (Krebs - Lreu6475)              -1.881 0.942  8  -4.0542 0.291591
# (Krebs - Secretin)              -2.304 0.999  8  -4.6082 0.000216
# (Krebs - (6475 + 5-27))          0.150 0.817  8  -1.7346 2.035110
# (Lreu6475 - Secretin)           -0.423 0.823  8  -2.3213 1.475843
# (Lreu6475 - (6475 + 5-27))       2.032 0.962  8  -0.1858 4.248940
# (Secretin - (6475 + 5-27))       2.454 1.021  8   0.0991 4.809468
# 
# sigma used for effect sizes: 86.19 
# Confidence level used: 0.95 

#4h oxytocin and secretin coming out of the same enteroids
pairedoxtsct<-read.delim(file="Fig4gh_HIEoxytocinsecretinpaired.txt",header=TRUE,sep="\t")

shapevalues<-c(23,24)
names(shapevalues)<-c("1","2")

pairOXTFigG<-ggplot(pairedoxtsct, aes(x=Treatment,y=Oxytocin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  #, lty=FAMI
  geom_point(aes(shape=as.factor(Experiment)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.42,.42,.42))+
  scale_shape_manual(values=shapevalues)+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
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

pdf(file="HIEpairOXT_Fig4G.pdf",width=1.2,height=2)
grid.arrange(pairOXTFigG)
dev.off()

FigGmodel<-lmer(Oxytocin~Treatment+(1|Experiment),data=pairedoxtsct, control=lmerControl(optimizer = "bobyqa"))
summary(FigGmodel)
FigGnullmodel<-lmer(Oxytocin~(1|Experiment),data=pairedoxtsct, control=lmerControl(optimizer = "bobyqa"))
 anova(FigGmodel, FigGnullmodel)
# Data: pairedoxtsct
# Models:
#   FigHnullmodel: Oxytocin ~ (1 | Experiment)
# FigHmodel: Oxytocin ~ Treatment + (1 | Experiment)
# npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)    
# FigHnullmodel    3 96.654 96.892 -45.327   90.654                        
# FigHmodel        4 86.654 86.972 -39.327   78.654    12  1   0.000532 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
performance::r2(FigGmodel)
# Conditional R2: 0.824
# Marginal R2: 0.701

FigGemmeans<-emmeans(FigGmodel,pairwise~Treatment,adjust="BH")

# $emmeans
# Treatment emmean   SE   df lower.CL upper.CL
# LDM4         7.5 25.2 1.58   -134.5      149
# Lreu6475   129.4 25.2 1.58    -12.6      271
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast        estimate   SE df t.ratio p.value
# LDM4 - Lreu6475     -122 23.1  5  -5.281  0.0032
# 
# Degrees-of-freedom method: kenward-roger 

eff_size(FigGemmeans,sigma = sigma(FigGmodel), edf=5)
# contrast          effect.size   SE df lower.CL upper.CL
# (LDM4 - Lreu6475)       -3.73 1.38  5    -7.27   -0.196
# 
# sigma used for effect sizes: 32.64 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95  

#4I oxytocin and secretin coming out of the same enteroids
pairSCTFigH<-ggplot(pairedoxtsct, aes(x=Treatment,y=Secretin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+ 
  geom_point(aes(shape=as.factor(Experiment)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Secretin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.42,.42,.42))+
  scale_shape_manual(values=shapevalues)+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
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

pdf(file="HIEpairSCT_Fig4H.pdf",width=1.2,height=2)
grid.arrange(pairSCTFigH)
dev.off()

FigHmodel<-lmer(Secretin~Treatment+(1|Experiment),data=pairedoxtsct,control=lmerControl(optimizer = "bobyqa"))
summary(FigHmodel)
FigHnullmodel<-lmer(Secretin~(1|Experiment),data=pairedoxtsct,control=lmerControl(optimizer = "bobyqa"))
anova(FigHmodel, FigHnullmodel)
# Data: pairedoxtsct
# Models:
#   FigInullmodel: Secretin ~ (1 | Experiment)
# FigImodel: Secretin ~ Treatment + (1 | Experiment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# FigInullmodel    3 86.987 87.225 -40.493   80.987                         
# FigImodel        4 62.756 63.074 -27.378   54.756 26.231  1  3.029e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

performance::r2(FigHmodel)
# Conditional R2: 0.964
# Marginal R2: 0.952

FigHemmeans<-emmeans(FigHmodel,pairwise~Treatment,adjust="BH")
# $emmeans
# Treatment emmean   SE   df lower.CL upper.CL
# LDM4         7.5 5.01 1.99    -14.2     29.2
# Lreu6475    82.4 5.01 1.99     60.7    104.1
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast        estimate   SE df t.ratio p.value
# LDM4 - Lreu6475    -74.9 5.52  5 -13.579  <.0001
# 
# Degrees-of-freedom method: kenward-roger 

eff_size(FigHemmeans,sigma=sigma(FigHmodel),edf=5)
# contrast          effect.size   SE df lower.CL upper.CL
# (LDM4 - Lreu6475)        -9.6 3.12  5    -17.6    -1.59
# 
# sigma used for effect sizes: 7.803 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

#4I secretin released by Lreu from LifeGift
LifeGiftSecretin<-read.delim(file="Fig4I_LifeGiftSecretinData.txt",header=TRUE,sep="\t",check.names = FALSE)
LGtissuesize<-read.delim(file="Fig2b_LifeGiftTissuesizes.txt",header = TRUE,sep="\t")
LGSecretinmerge<-merge(LifeGiftSecretin,LGtissuesize,by=c("LifeGift","Segment"))
LGSecretinmerge$SCTnorm<-LGSecretinmerge$Secretin/LGSecretinmerge$Area

shapevalues<-c(23,24,25)
names(shapevalues)<-levels(as.factor(LGSecretinmerge$LifeGift))
LGSecretinmerge$Segment<-factor(LGSecretinmerge$Segment,c("Duodenum","UpperJejunum","LowerJejunum",
                                                        "UpperIleum","LowerIleum", "AscendingColon","TransverseColon", "DescendingColon" ))
LGSecretinmerge$SegmentTreatment=paste(LGSecretinmerge$Segment,LGSecretinmerge$Treatment,sep="_")
LGSecretinmerge$SegmentTreatment<-factor(LGSecretinmerge$SegmentTreatment,c("Duodenum_LDM4", "Duodenum_6475","UpperJejunum_LDM4", "UpperJejunum_6475" ,
                                                                          "LowerJejunum_LDM4","LowerJejunum_6475",
                                                                          "UpperIleum_LDM4", "UpperIleum_6475",
                                                                          "LowerIleum_LDM4","LowerIleum_6475",
                                                                          "AscendingColon_LDM4","AscendingColon_6475","TransverseColon_LDM4","TransverseColon_6475", 
                                                                          "DescendingColon_LDM4", "DescendingColon_6475"))

LifeGiftSCTplot<-ggplot(LGSecretinmerge, aes(x=Segment,y=SCTnorm,fill=Treatment,group=SegmentTreatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(aes(shape=LifeGift),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Secretin (pg/ml/'*cm^2*')')) +
  scale_y_continuous(expand=c(0.05,.4,.4,.4))+
  scale_fill_manual(values=c( "#4DBBD5FF","#E64B35FF"))+
  scale_shape_manual(values=shapevalues)+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.title.x= element_text(size=8),
    axis.title.y= element_text(size=8),
    axis.text.x= element_text(size=7, colour = "black"),
    axis.text.y= element_text(size=7, colour = "black"),
    panel.border = element_blank(),
    strip.background =element_blank(),
    strip.text = element_text(size=7),
    plot.title=element_blank(),
    legend.position="none"
  )


pdf(file="LifeGiftSCT.pdf",width=5,height=2)
grid.arrange(LifeGiftSCTplot)
dev.off()

pdf(file="FigGHI.pdf",width=7,height=2)
grid.arrange(pairOXTFigG,pairSCTFigH,LifeGiftSCTplot,widths=c(1,1,4))
dev.off()

FigImodel<-lmer(SCTnorm~Treatment*Segment + (1|LifeGift),data=LGSecretinmerge, control=lmerControl(optimizer = "bobyqa"))
FigInullmodel<-lmer(SCTnorm~ (1|LifeGift),data=LGSecretinmerge, control=lmerControl(optimizer = "bobyqa"))
anova(FigImodel, FigInullmodel)
# Data: LGSecretinmerge
# Models:
#   FigInullmodel: SCTnorm ~ (1 | LifeGift)
# FigImodel: SCTnorm ~ Treatment * Segment + (1 | LifeGift)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# FigInullmodel    3 1477.4 1486.3 -735.72   1471.4                         
# FigImodel       18 1341.8 1395.3 -652.91   1305.8 165.62 15  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
performance::r2(FigImodel)
# Conditional R2: 0.689
# Marginal R2: 0.614

FigIemmeans<-emmeans(FigImodel,pairwise~Treatment|Segment,adjust="BH")
# $emmeans
# Segment = Duodenum:
#   Treatment emmean   SE   df lower.CL upper.CL
# 6475       64.42 10.2 9.39    41.50     87.3
# LDM4        9.69 10.2 9.39   -13.23     32.6
# 
# Segment = UpperJejunum:
#   Treatment emmean   SE   df lower.CL upper.CL
# 6475      104.71 10.2 9.39    81.79    127.6
# LDM4       11.26 10.2 9.39   -11.66     34.2
# 
# Segment = LowerJejunum:
#   Treatment emmean   SE   df lower.CL upper.CL
# 6475       96.72 10.2 9.39    73.80    119.6
# LDM4       10.29 10.2 9.39   -12.64     33.2
# 
# Segment = UpperIleum:
#   Treatment emmean   SE   df lower.CL upper.CL
# 6475       69.20 10.2 9.39    46.28     92.1
# LDM4       11.06 10.2 9.39   -11.86     34.0
# 
# Segment = LowerIleum:
#   Treatment emmean   SE   df lower.CL upper.CL
# 6475       38.97 10.2 9.39    16.05     61.9
# LDM4       12.36 10.2 9.39   -10.56     35.3
# 
# Segment = AscendingColon:
#   Treatment emmean   SE   df lower.CL upper.CL
# 6475       13.91 10.2 9.39    -9.01     36.8
# LDM4        9.69 10.2 9.39   -13.24     32.6
# 
# Segment = TransverseColon:
#   Treatment emmean   SE   df lower.CL upper.CL
# 6475       10.61 10.2 9.39   -12.31     33.5
# LDM4        8.85 10.2 9.39   -14.07     31.8
# 
# Segment = DescendingColon:
#   Treatment emmean   SE   df lower.CL upper.CL
# 6475        9.98 10.2 9.39   -12.94     32.9
# LDM4        8.94 10.2 9.39   -13.98     31.9
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Segment = Duodenum:
#   contrast    estimate SE  df t.ratio p.value
# 6475 - LDM4    54.73 11 126   4.985  <.0001
# 
# Segment = UpperJejunum:
#   contrast    estimate SE  df t.ratio p.value
# 6475 - LDM4    93.45 11 126   8.511  <.0001
# 
# Segment = LowerJejunum:
#   contrast    estimate SE  df t.ratio p.value
# 6475 - LDM4    86.43 11 126   7.872  <.0001
# 
# Segment = UpperIleum:
#   contrast    estimate SE  df t.ratio p.value
# 6475 - LDM4    58.14 11 126   5.295  <.0001
# 
# Segment = LowerIleum:
#   contrast    estimate SE  df t.ratio p.value
# 6475 - LDM4    26.61 11 126   2.424  0.0168
# 
# Segment = AscendingColon:
#   contrast    estimate SE  df t.ratio p.value
# 6475 - LDM4     4.22 11 126   0.385  0.7012
# 
# Segment = TransverseColon:
#   contrast    estimate SE  df t.ratio p.value
# 6475 - LDM4     1.76 11 126   0.160  0.8729
# 
# Segment = DescendingColon:
#   contrast    estimate SE  df t.ratio p.value
# 6475 - LDM4     1.03 11 126   0.094  0.9252
# Degrees-of-freedom method: kenward-roger 

eff_size(FigIemmeans, sigma=sigma(FigImodel), edf=126)
# Segment = Duodenum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (6475 - LDM4)      2.3497 0.494 126    1.372    3.328
# 
# Segment = UpperJejunum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (6475 - LDM4)      4.0121 0.535 126    2.954    5.071
# 
# Segment = LowerJejunum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (6475 - LDM4)      3.7110 0.526 126    2.670    4.752
# 
# Segment = UpperIleum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (6475 - LDM4)      2.4960 0.497 126    1.513    3.479
# 
# Segment = LowerIleum:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (6475 - LDM4)      1.1425 0.477 126    0.199    2.086
# 
# Segment = AscendingColon:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (6475 - LDM4)      0.1813 0.472 126   -0.752    1.114
# 
# Segment = TransverseColon:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (6475 - LDM4)      0.0756 0.471 126   -0.857    1.009
# 
# Segment = DescendingColon:
#   contrast      effect.size    SE  df lower.CL upper.CL
# (6475 - LDM4)      0.0444 0.471 126   -0.889    0.977
# 
# sigma used for effect sizes: 23.29 
# Degrees-of-freedom method: inherited from kenward-roger when re-gridding 
# Confidence level used: 0.95 

#Secretin concentrations

secretin<-read.delim(file="FigS5A_Secretin.txt",header=TRUE,sep="\t")

shapevalues<-c(23,24,25)
names(shapevalues)<-levels(as.factor(secretin$Experiment))


SCTplot<-ggplot(secretin, aes(x=as.factor(SCT/100),y=OXT,group=SCT)) +
  geom_boxplot(size=0.5,outlier.shape = NA,alpha = 0.4,colour="black")+ 
  geom_point(aes(shape=as.factor(Experiment)),size=1.4,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black",fill="lightblue" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  xlab(bquote('Secretin (ng/ml)')) +
  scale_y_continuous(expand=c(0.05,.42,.42,.42))+
  scale_shape_manual(values=c(21,21,21))+
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
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

pdf(file="SecretinConcentration2.pdf",width=2,height=2)
SCTplot
dev.off()
