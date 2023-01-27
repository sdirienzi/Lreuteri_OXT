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
HIEoxt<-read.delim(file="Fig4a_HIEallSecretinOxytocin.txt",header=TRUE,sep="\t",check.names=FALSE)

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
pdf(file="HIEsecretinoxytocin.pdf",width=1.3,height=2)
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
HIEalloxt<-read.delim(file="Fig4b_HIEallSecretinOxytocin.txt",header=TRUE,sep="\t",check.names=FALSE)
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

pdf(file="HIEallsecretinoxytocin.pdf",width=2.6,height=2)
grid.arrange(HIEallFig4B)
dev.off()

Figbmodel<-lm(Oxytocin~Treatment+(Enteroid),data=HIEalloxt)
summary(Figbmodel)
# Call:
#   lm(formula = Oxytocin ~ Treatment + (Enteroid), data = HIEalloxt)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -31.838  -8.230  -5.060   4.056  49.207 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         15.7298     8.9784   1.752   0.0959 .  
# TreatmentSecretin  113.9540     8.0305  14.190 1.46e-11 ***
#   EnteroidD109        -0.4529    11.3569  -0.040   0.9686    
# EnteroidIL104      -23.9739    11.3569  -2.111   0.0483 *  
#   EnteroidJ11         -8.4923    11.3569  -0.748   0.4638    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 19.67 on 19 degrees of freedom
# Multiple R-squared:  0.916,	Adjusted R-squared:  0.8983 
# F-statistic:  51.8 on 4 and 19 DF,  p-value: 5.854e-10
HIEemmeans<-emmeans(Figbmodel,pairwise~Treatment+Enteroid,adjust="BH")
# $emmeans
# Treatment Enteroid emmean   SE df lower.CL upper.CL
# Krebs     ASC209    15.73 8.98 19    -3.06     34.5
# Secretin  ASC209   129.68 8.98 19   110.89    148.5
# Krebs     D109      15.28 8.98 19    -3.52     34.1
# Secretin  D109     129.23 8.98 19   110.44    148.0
# Krebs     IL104     -8.24 8.98 19   -27.04     10.5
# Secretin  IL104    105.71 8.98 19    86.92    124.5
# Krebs     J11        7.24 8.98 19   -11.55     26.0
# Secretin  J11      121.19 8.98 19   102.40    140.0
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                           estimate    SE df t.ratio p.value
# Krebs ASC209 - Secretin  ASC209    -113.954  8.03 19 -14.190  <.0001
# Krebs ASC209 - Krebs D109             0.453 11.36 19   0.040  0.9686
# Krebs ASC209 - Secretin  D109      -113.501 13.91 19  -8.160  <.0001
# Krebs ASC209 - Krebs IL104           23.974 11.36 19   2.111  0.0731
# Krebs ASC209 - Secretin  IL104      -89.980 13.91 19  -6.469  <.0001
# Krebs ASC209 - Krebs J11              8.492 11.36 19   0.748  0.5251
# Krebs ASC209 - Secretin  J11       -105.462 13.91 19  -7.582  <.0001
# Secretin  ASC209 - Krebs D109       114.407 13.91 19   8.225  <.0001
# Secretin  ASC209 - Secretin  D109     0.453 11.36 19   0.040  0.9686
# Secretin  ASC209 - Krebs IL104      137.928 13.91 19   9.916  <.0001
# Secretin  ASC209 - Secretin  IL104   23.974 11.36 19   2.111  0.0731
# Secretin  ASC209 - Krebs J11        122.446 13.91 19   8.803  <.0001
# Secretin  ASC209 - Secretin  J11      8.492 11.36 19   0.748  0.5251
# Krebs D109 - Secretin  D109        -113.954  8.03 19 -14.190  <.0001
# Krebs D109 - Krebs IL104             23.521 11.36 19   2.071  0.0731
# Krebs D109 - Secretin  IL104        -90.433 13.91 19  -6.502  <.0001
# Krebs D109 - Krebs J11                8.039 11.36 19   0.708  0.5251
# Krebs D109 - Secretin  J11         -105.915 13.91 19  -7.615  <.0001
# Secretin  D109 - Krebs IL104        137.475 13.91 19   9.884  <.0001
# Secretin  D109 - Secretin  IL104     23.521 11.36 19   2.071  0.0731
# Secretin  D109 - Krebs J11          121.993 13.91 19   8.771  <.0001
# Secretin  D109 - Secretin  J11        8.039 11.36 19   0.708  0.5251
# Krebs IL104 - Secretin  IL104      -113.954  8.03 19 -14.190  <.0001
# Krebs IL104 - Krebs J11             -15.482 11.36 19  -1.363  0.2402
# Krebs IL104 - Secretin  J11        -129.436 13.91 19  -9.306  <.0001
# Secretin  IL104 - Krebs J11          98.472 13.91 19   7.080  <.0001
# Secretin  IL104 - Secretin  J11     -15.482 11.36 19  -1.363  0.2402
# Krebs J11 - Secretin  J11          -113.954  8.03 19 -14.190  <.0001
# 
# P value adjustment: BH method for 28 tests 

eff_size(HIEemmeans,sigma = sigma(Figbmodel),edf=19)
# contrast                             effect.size    SE df lower.CL upper.CL
# (Krebs ASC209 - Secretin  ASC209)         -5.793 1.025 19  -7.9376   -3.649
# (Krebs ASC209 - Krebs D109)                0.023 0.577 19  -1.1854    1.231
# (Krebs ASC209 - Secretin  D109)           -5.770 1.173 19  -8.2254   -3.315
# (Krebs ASC209 - Krebs IL104)               1.219 0.610 19  -0.0585    2.496
# (Krebs ASC209 - Secretin  IL104)          -4.574 1.025 19  -6.7197   -2.429
# (Krebs ASC209 - Krebs J11)                 0.432 0.582 19  -0.7855    1.649
# (Krebs ASC209 - Secretin  J11)            -5.361 1.121 19  -7.7074   -3.015
# (Secretin  ASC209 - Krebs D109)            5.816 1.179 19   3.3483    8.284
# (Secretin  ASC209 - Secretin  D109)        0.023 0.577 19  -1.1854    1.231
# (Secretin  ASC209 - Krebs IL104)           7.012 1.339 19   4.2086    9.815
# (Secretin  ASC209 - Secretin  IL104)       1.219 0.610 19  -0.0585    2.496
# (Secretin  ASC209 - Krebs J11)             6.225 1.233 19   3.6446    8.805
# (Secretin  ASC209 - Secretin  J11)         0.432 0.582 19  -0.7855    1.649
# (Krebs D109 - Secretin  D109)             -5.793 1.025 19  -7.9376   -3.649
# (Krebs D109 - Krebs IL104)                 1.196 0.609 19  -0.0790    2.471
# (Krebs D109 - Secretin  IL104)            -4.597 1.028 19  -6.7484   -2.446
# (Krebs D109 - Krebs J11)                   0.409 0.581 19  -0.8076    1.625
# (Krebs D109 - Secretin  J11)              -5.384 1.124 19  -7.7365   -3.032
# (Secretin  D109 - Krebs IL104)             6.989 1.336 19   4.1922    9.785
# (Secretin  D109 - Secretin  IL104)         1.196 0.609 19  -0.0790    2.471
# (Secretin  D109 - Krebs J11)               6.202 1.230 19   3.6280    8.776
# (Secretin  D109 - Secretin  J11)           0.409 0.581 19  -0.8076    1.625
# (Krebs IL104 - Secretin  IL104)           -5.793 1.025 19  -7.9376   -3.649
# (Krebs IL104 - Krebs J11)                 -0.787 0.591 19  -2.0246    0.451
# (Krebs IL104 - Secretin  J11)             -6.580 1.280 19  -9.2600   -3.900
# (Secretin  IL104 - Krebs J11)              5.006 1.077 19   2.7523    7.260
# (Secretin  IL104 - Secretin  J11)         -0.787 0.591 19  -2.0246    0.451
# (Krebs J11 - Secretin  J11)               -5.793 1.025 19  -7.9376   -3.649
# 
# sigma used for effect sizes: 19.67 
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

#Fig 4D antibody block experiment
Abblock<-read.delim(file="Fig4d_HIEABblocker.txt",header=TRUE,sep="\t",check.names=FALSE)
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

FigDmodel<-lm(Oxytocin~Treatment,data=Abblock)
summary(FigDmodel)
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

FigDemmeans<-emmeans(FigDmodel,pairwise~Treatment,adjust="BH")
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

#4e 5-27 blockng experiment on NGN3s
block527data<-read.delim(file="Fig4e_527inhibitor.txt",header=TRUE,sep="\t",check.names = FALSE)

block527data$Treatment<-factor(block527data$Treatment,c("Krebs","Inhibitor","Secretin","Secretin+ Inhibitor"))

shapevalues<-c(23,24,25,22)
names(shapevalues)<-c("1","2","3","4")

HIEFig4E<-ggplot(block527data, aes(x=Treatment,y=Oxytocin,fill=Treatment)) +
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
pdf(file="HIE527block_Fig4E.pdf",width=1.6,height=1.6)
grid.arrange(HIEFig4E)
dev.off()

FigEmodel<-lmer(Oxytocin~Treatment+(1|Experiment),data=block527data, control=lmerControl(optimizer = "bobyqa"))
FigEnull<-lmer(Oxytocin~(1|Experiment),data=block527data, control=lmerControl(optimizer = "bobyqa"))
anova(FigEmodel,FigEnull)
performance::r2(FigEmodel)
# Conditional R2: 0.931
# Marginal R2: 0.925

summary(FigEmodel)
FigEemmeans<-emmeans(FigEmodel,pairwise~Treatment,adjust="BH")
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


eff_size(FigEemmeans,sigma = sigma(FigEmodel), edf=45)
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

#4f 5-27 on HIE with Lreu
block527Lreudata<-read.delim(file="Fig4f_secretininhibitor.txt",header=TRUE,sep="\t",check.names = FALSE)

block527Lreudata$Treatment<-factor(block527Lreudata$Treatment,c("Krebs","Secretin_inhibitor","6475","6475_5-27"))

shapevalues<-c(23,24,25)
names(shapevalues)<-c("1","2","3")

HIEFig4F<-ggplot(block527Lreudata, aes(x=Treatment,y=Oxytocin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  
  geom_point(aes(shape=as.factor(Experiment)),size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.6,.6,.6))+
  scale_fill_manual(values=c("#E64B35FF","#F39B7FFF","#00A087FF","#8491B4FF"))+
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
pdf(file="block527Lreudata_Fig4F.pdf",width=1.6,height=1.6)
grid.arrange(HIEFig4F)
dev.off()

FigFmodel<-lmer(Oxytocin~Treatment+(1|Experiment),data=block527Lreudata, control=lmerControl(optimizer = "bobyqa"))
FigFnull<-lmer(Oxytocin~(1|Experiment),data=block527Lreudata, control=lmerControl(optimizer = "bobyqa"))
anova(FigFmodel,FigFnull)
# Data: block527Lreudata
# Models:
#   FigFnull: Oxytocin ~ (1 | Experiment)
# FigFmodel: Oxytocin ~ Treatment + (1 | Experiment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# FigFnull     3 564.85 570.47 -279.43   558.85                         
# FigFmodel    6 496.53 507.75 -242.26   484.53 74.328  3  5.048e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
performance::r2(FigFmodel)
# Conditional R2: 0.794
# Marginal R2: 0.763

summary(FigFmodel)
FigFemmeans<-emmeans(FigFmodel,pairwise~Treatment,adjust="BH")
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
eff_size(FigFemmeans, sigma = sigma(FigFmodel), edf=42)
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
LGblock527data<-read.delim(file="Fig4g_LifeGift527.txt",header=TRUE,sep="\t",check.names = FALSE)
LGblock527data$Normalized<-LGblock527data$Oxytocin/LGblock527data$Area

LGblock527data$Treatment<-factor(LGblock527data$Treatment,c("Krebs","Lreu6475","Secretin","6475 + 5-27"))

LGFig4G<-ggplot(LGblock527data, aes(x=Treatment,y=Normalized,fill=Treatment)) +
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
pdf(file="LG527block_Fig4G.pdf",width=2,height=2)
grid.arrange(LGFig4G)
dev.off()

FigGmodel<-lm(Normalized~Treatment,data=LGblock527data)
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

summary(FigGmodel)
FigGemmeans<-emmeans(FigGmodel,pairwise~Treatment,adjust="BH")
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
pairedoxtsct<-read.delim(file="Figh_HIEoxytocinsecretinpaired.txt",header=TRUE,sep="\t")

pairOXTFigH<-ggplot(pairedoxtsct, aes(x=Treatment,y=Oxytocin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+  #, lty=FAMI
  geom_point(shape=21,size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Oxytocin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.42,.42,.42))+
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

pdf(file="HIEpairOXT_Fig4H.pdf",width=1.2,height=2)
grid.arrange(pairOXTFigH)
dev.off()

FigHmodel<-lm(Oxytocin~Treatment,data=pairedoxtsct)
summary(FigHmodel)
# Call:
#   lm(formula = Oxytocin ~ Treatment, data = pairedoxtsct)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -52.761  -9.259   0.000   5.172  69.107 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)           7.50      19.75   0.380  0.71718   
# TreatmentLreu6475   121.90      27.93   4.365  0.00475 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 39.5 on 6 degrees of freedom
# Multiple R-squared:  0.7605,	Adjusted R-squared:  0.7206 
# F-statistic: 19.05 on 1 and 6 DF,  p-value: 0.004746

FigHemmeans<-emmeans(FigHmodel,pairwise~Treatment,adjust="BH")

# $emmeans
# Treatment emmean   SE df lower.CL upper.CL
# LDM4         7.5 19.7  6    -40.8     55.8
# Lreu6475   129.4 19.7  6     81.1    177.7
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast        estimate   SE df t.ratio p.value
# LDM4 - Lreu6475     -122 27.9  6  -4.365  0.0047

eff_size(FigHemmeans,sigma = sigma(FigHmodel), edf=6)
# contrast          effect.size   SE df lower.CL upper.CL
# (LDM4 - Lreu6475)       -3.09 1.14  6    -5.87   -0.303
# 
# sigma used for effect sizes: 39.5 
# Confidence level used: 0.95 

#4I oxytocin and secretin coming out of the same enteroids
pairSCTFigI<-ggplot(pairedoxtsct, aes(x=Treatment,y=Secretin,fill=Treatment)) +
  geom_boxplot(size=0.8,outlier.shape = NA,alpha = 0.4,colour="black")+ 
  geom_point(shape=21,size=2,stroke=0.2, position=position_jitterdodge(jitter.width=0.2 ),color="black" )+
  ylab(bquote('Secretin (pg/ml)')) +
  scale_y_continuous(expand=c(0.05,.42,.42,.42))+
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

pdf(file="HIEpairSCT_Fig4I.pdf",width=1.2,height=2)
grid.arrange(pairSCTFigI)
dev.off()

FigImodel<-lm(Secretin~Treatment,data=pairedoxtsct)
summary(FigImodel)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -16.1162  -0.1497   0.0000   0.9122  13.0664 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          7.500      4.302   1.743    0.132    
# TreatmentLreu6475   74.925      6.084  12.316 1.75e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.604 on 6 degrees of freedom
# Multiple R-squared:  0.9619,	Adjusted R-squared:  0.9556 
# F-statistic: 151.7 on 1 and 6 DF,  p-value: 1.747e-05

FigIemmeans<-emmeans(FigImodel,pairwise~Treatment,adjust="BH")
# $emmeans
# Treatment emmean  SE df lower.CL upper.CL
# LDM4         7.5 4.3  6    -3.03       18
# Lreu6475    82.4 4.3  6    71.90       93
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast        estimate   SE df t.ratio p.value
# LDM4 - Lreu6475    -74.9 6.08  6 -12.316  <.0001
# contrast          effect.size   SE df lower.CL upper.CL
# (LDM4 - Lreu6475)       -8.71 2.61  6    -15.1    -2.32
# 
# sigma used for effect sizes: 8.604 
# Confidence level used: 0.95 

eff_size(FigIemmeans,sigma=sigma(FigImodel),edf=6)
# contrast          effect.size   SE df lower.CL upper.CL
# (LDM4 - Lreu6475)       -8.71 2.61  6    -15.1    -2.32
# 
# sigma used for effect sizes: 8.604 
# Confidence level used: 0.95 

#4J secretin released by Lreu from LifeGift
LifeGiftSecretin<-read.delim(file="Fig4j_LifeGiftSecretinData.txt",header=TRUE,sep="\t",check.names = FALSE)
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

pdf(file="FigHIJ.pdf",width=7,height=2)
grid.arrange(pairOXTFigH,pairSCTFigI,LifeGiftSCTplot,widths=c(1,1,4))
dev.off()

FigJmodel<-lmer(SCTnorm~Treatment*Segment + (1|LifeGift),data=LGSecretinmerge, control=lmerControl(optimizer = "bobyqa"))
FigJnullmodel<-lmer(SCTnorm~ (1|LifeGift),data=LGSecretinmerge, control=lmerControl(optimizer = "bobyqa"))
anova(FigJmodel, FigJnullmodel)
performance::r2(FigJmodel)
# Conditional R2: 0.689
# Marginal R2: 0.614

FigJemmeans<-emmeans(FigJmodel,pairwise~Treatment|Segment,adjust="BH")
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

eff_size(FigJemmeans, sigma=sigma(FigJmodel), edf=126)
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
