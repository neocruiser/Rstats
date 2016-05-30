# This buffer is for notes you don't want to save, and for R code.
# If you want to create an *.Rnw file, run ~/perls/knitr.pl
# then enter the file's and project's name.
setwd("/home/neo/data/Dropbox/R/heatmap")

dat <- read.table("./raw.data.protein.txt", header=T, fill=T,row.names="Accession")
dat <- read.table("./significant.95.boot.txt", header=F, fill=T,row.names=1)
colnames(dat) <- c("v1","v2","v3","d1","d2","d3")
head(dat)
dim(dat)
dat=scale(dat, center=T, scale=T)
## Version I
p = prcomp(dat, retx=T)
plot(p)
scores = p$x
loadings <- p$rotation
sd <- p$sdev
## plot 1

pdf("prot.em.final.pca.pdf",width=7,height=5)

plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2", type="n", xlim=c(min(scores[,1:2]), max(scores[,1:2])), ylim=c(min(scores[,1:2]), max(scores[,1:2])))
arrows(0,0,loadings[,1]*5,loadings[,2]*5, length=0.1,angle=20, col="red")
text(loadings[,1]*5*1.2,loadings[,2]*5*1.2, rownames(loadings), col="black", cex=0.7)
dev.off()
# plot 2

#plot(scores[,1]/sd[1], scores[,2]/sd[2], xlab="PCA 1", ylab="PCA 2", type="n")
#arrows(0,0,loadings[,1]*sd[1],loadings[,2]*sd[2], length=0.1, angle=20, col="red")
#text(loadings[,1]*sd[1]*1.2,loadings[,2]*sd[2]*1.2,rownames(loadings), col="red", cex=0.7)
# plot 3
#biplot(p)
# plot 4
biplot(scores[,1:2], loadings[,1:2],  cex=0.7)
dev.off()
## Veersion II
fit <- princomp(dat, cor=F)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
pdf("prot.em.variance2.pdf",width=7,height=5)
plot(fit,type="lines") # scree plot
dev.off()
biplot(fit,cex = .5) 

## Version III
require(FactoMineR)
result <- PCA(dat)

====================
    linear regressions Diagnostics
====================
dat <- read.table("./test", header=T, fill=T,row.names="Accession")
head(dat)
plot(dat)
library(dplyr)
require(tidyr)
dat1 <- gather(dat, "samples","count",1:6)
summary(dat1)
x=dim(dat1)[1]/2
dat1 <- data.frame(dat1,condition=gl(2,x,2*x,labels=c("cdt1","cdt2")))
head(dat1)

fit <- lm(count~condition, data = dat1)
summary(fit)
#contrasts(growth$Treatment)
anova(fit)
cutoff <- 4/((nrow(dat1)-length(fit$coefficients)-2))
pdf("prot.em.cooks.pdf",width=7,height=5)
layout(matrix(c(1,2,3,4),2,2)) # optional layout
plot(fit, which=4, cook.levels = cutoff)	## identify D values
qqPlot(fit, main = "QQ plot")	## qq plot for studentized residuals
plot(fit)
dev.off()
====================
    Biometric analyses
====================

head(dat1)
source("/home/neo/data/Dropbox/R/01funcs.R")
DF <- summarySE(dat1, measurevar="count", groupvars=c("samples","condition"))
DF
# CONFIDENCE INTERVAL OF THE MEAN OF GROWTH
require(ggplot2)
pdf("prot.em.se.pdf",width=7,height=5)
ggplot(DF, aes(x=samples, y=count)) +
  geom_errorbar(aes(ymin=count-se, ymax=count+se), width=.4) +
  geom_line(aes(group=condition)) +
  geom_point(size=4.4, fill="black",aes(shape=condition)) +
    theme_bw()
dev.off()

## two-way ANOVA
attach(dat1)
nova2way <- aov(count ~ condition, dat1)
detach(dat1)
summary(nova2way)				##type I sequential SS
drop1(nova2way, ~., test="F")				## type III marginal SS
TukeyHSD(nova2way)				## multiple comparisons test

layout(matrix(c(1,2,3,4),2,2)) # optional layout
plot(nova2way) # diagnostic plots


library(e1071)
head(dat)
skewed <- apply(dat,2, function(x) skewness(x))


## Shapiro wilk (W). Test for normality
normality <- apply(dat, 2, function(x) shapiro.test(x))

## STEM and leaf plot
neutral_stem <- apply(dat, 2, function(x) stem(x))


## Kolmogorov-Smirnov test
neutral_normality_kolSmir <- apply(dat, 2, function(x) ks.test(x, "pnorm", mean=mean(x), sd=sqrt(var(x))))




setwd("/home/neo/data/Dropbox/R/heatmap")
genre <- read.table("./test", header=T, fill=T,row.names="Accession")
#colnames(genre) <- c("H1","DM1","H2","H3","DM2","D1","H4","D2","H5","DM3","D3","H6","SW1","SW2","SW3","C")
genre=dat
head(genre)
genre <- as.matrix(genre)
## HIERARCHICAL AND BOOTSTRAP ANALYSIS
## (help) source : http://tinyurl.com/lqs76tf
require(pvclust)
require(gplots)

## load sets
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
rawdata <- genre
scaledata=scale(genre, center=T, scale=T)
#scaledata=scale(genre, center = F, scale = apply(genre, 2, sd, na.rm=T))
hra <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete")
hca <- hclust(as.dist(1-cor(scaledata, method="pearson")), method="complete")
heatmap(rawdata, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row")

## CUT THE TREE
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
mycl <- cutree(hra, h=max(hra$height)/2)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
pdf("prot.em.final.cluster.pdf",width=7,height=5)
heatmap(rawdata, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", RowSideColors=mycolhc)
dev.off()

## BOOTSTRAPING
n=10;a=0.99
#pvData <- pvclust(scale(t(rawdata)), method.dist="correlation", method.hclust="complete", nboot=n)
?pvclust::pvclust
pvData <- pvclust(scaledata, method.dist="correlation", method.hclust="complete", nboot=n, use.cor = "pairwise.complete.obs", weight=T, r=seq(0.5,0.9,by=0.1))
pvData <- pvclust(t(scaledata), method.dist="correlation", method.hclust="complete", nboot=n, use.cor = "pairwise.complete.obs",r=seq(0.5,0.9,by=0.05))
pdf("prot.em.bootstrap.pdf",width=7,height=5)
plot(pvData, hang=-1)
pvrect(pvData, alpha=a)
dev.off()

## RETRIEVE MEMBERS OF SIGNIFICANT CLUSTERS.
clsig <- unlist(pvpick(pvData, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters)
length(clsig)
write.table(clsig, "bootstrap.significat.txt",row.names=F,quote=F)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pvData$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
heatmap(rawdata, Rowv=dend_colored, Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", RowSideColors=mycolhc)

## PLOT HEATMAP WITH HEATMAP.2() FUNCTION WHICH SCALES BETTER FOR MANY ENTRIES.
pdf("prot.em.bootcluster.pdf",width=7,height=5)
x11(height=5,width =8)
heatmap.2(rawdata, Rowv=dend_colored, Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc,margins=c(8,20))

dev.off()
