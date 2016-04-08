setwd("c://Users/Neo/Downloads/")
setwd("/home/neo/data/Dropbox/R/ganglia/data")

setwd("/home/neo/data/Dropbox/R/amelie/")
dat <- read.table("bacteria.family.txt")
dat <- read.table("test.txt")
dat1=dat
head(dat1)
dat=dat1[,13:24]
dim(dat)
## Version I
p = prcomp(dat, retx=T)
plot(p)
scores = p$x
loadings <- p$rotation
sd <- p$sdev
## plot 1

scores=read.table("./gg.scores")
loadings=read.table("./gg.loadings")

pdf("GG.pca.pdf",width=7,height=5)

plot(scores[,1], scores[,2], xlab="PCA 1", ylab="PCA 2", type="n", xlim=c(min(scores[,1:2]), max(scores[,1:2])), ylim=c(min(scores[,1:2]), max(scores[,1:2])))
arrows(0,0,loadings[,1]*50,loadings[,2]*50, length=0.1,angle=20, col="red")
text(loadings[,1]*50*1.3,loadings[,2]*50*1.3, rownames(loadings), col="black", cex=0.9)
dev.off()
# plot 2
plot(scores[,1]/sd[1], scores[,2]/sd[2], xlab="PCA 1", ylab="PCA 2", type="n")
arrows(0,0,loadings[,1]*sd[1],loadings[,2]*sd[2], length=0.1, angle=20, col="red")
text(loadings[,1]*sd[1]*1.2,loadings[,2]*sd[2]*1.2,rownames(loadings), col="red", cex=0.7)
# plot 3
biplot(p)
# plot 4
biplot(scores[,1:2], loadings[,1:2],  cex=0.7)

## Veersion II
fit <- princomp(dat, cor=F)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit,type="lines") # scree plot
biplot(fit) 

## Version III
require(FactoMineR)
result <- PCA(dat)

setwd("/home/neo/data/Dropbox/R/amelie/")
dat <- read.table("bacteria.genus.txt")
head(dat)
dat <- dat[-84, -c(1,3:4,7,9,12,13:15)]
x <- apply(dat, 2, function(x) abs((x-dat[ ,7])))
head(x)
result <- PCA(x[, -7])

dat=x
