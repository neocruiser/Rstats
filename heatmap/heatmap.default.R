
?RColorBrewer
library(RColorBrewer)
display.brewer.pal(9, "YlGnBu")
palett<-brewer.pal(9,"Greens")

## Load SoftClustering Mfuzz package --- 12 clusters
setwd("/home/neo/data/Dropbox/R/heatmap")

genre <- read.table("./test", header=T, fill=T,row.names="Accession")
#colnames(genre) <- c("H1","DM1","H2","H3","DM2","D1","H4","D2","H5","DM3","D3","H6","SW1","SW2","SW3","C")
head(genre)
genre <- as.matrix(genre)
heatmap.g <- heatmap(genre, Rowv=NA, Colv=NA, col=palett, margins=c(1,40), scale="row", labCol=c(seq(1:16)))


## HIERARCHICAL AND BOOTSTRAP ANALYSIS
## (help) source : http://tinyurl.com/lqs76tf
require(pvclust)
require(gplots)

## load sets
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
## by sample
rawdata <- t(genre)
scaledata=t(scale(genre))

## by bacteria
rawdata <- genre
scaledata=scale(genre)


hra <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete")
hca <- hclust(as.dist(1-cor(scaledata, method="pearson")), method="complete")
heatmap(rawdata, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row")

## CUT THE TREE
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
mycl <- cutree(hra, h=max(hra$height)/2)
mycolhc <- sample(rainbow(256))
mycolhc <- mycolhc[as.vector(mycl)]
heatmap(rawdata, Rowv=as.dendrogram(hra), Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", RowSideColors=mycolhc)

## BOOTSTRAPING
n=200;a=0.95
pvData <- pvclust(scale(t(rawdata)), method.dist="correlation", method.hclust="ward.D2", nboot=n)
plot(pvData, hang=-1)
pvrect(pvData, alpha=a)

## RETRIEVE MEMBERS OF SIGNIFICANT CLUSTERS.
clsig <- unlist(pvpick(pvData, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R") # Import tree coloring function.
dend_colored <- dendrapply(as.dendrogram(pvData$hclust), dendroCol, keys=clsig, xPar="edgePar", bgr="black", fgr="red", pch=20)
heatmap(rawdata, Rowv=dend_colored, Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", RowSideColors=mycolhc)

## PLOT HEATMAP WITH HEATMAP.2() FUNCTION WHICH SCALES BETTER FOR MANY ENTRIES.
x11(height=5,width =8)
heatmap.2(rawdata, Rowv=dend_colored, Colv=as.dendrogram(hca), col=my.colorFct(), scale="row", trace="none", RowSideColors=mycolhc,margins=c(8,20))

dev.off()
