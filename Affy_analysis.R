##setwd and download required packages from bioconductor
setwd("")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("oligo")
library(oligo)
library(pd.mogene.2.0.st)

##read CEL files using read.celfiles func from oligo package
celFiles=list.celfiles("./", full.names=T)
celFiles
rawData = read.celfiles(celFiles)
hist(rawData, target="core")
boxplot(rawData, target="core", col=rainbow(10), las=2)
#rawexp<-exprs(rawData)

#normalize rawdata and get exp data
eset<-oligo::rma(rawData)
eset.expression<-exprs(eset)
boxplot(eset.expression, las=2, col=rainbow(10))

#clustering data
dist.mat<-dist(t(eset.expression), method="euclidean")
clusters<-hclust(dist.mat, method="ward")
plot(clusters)

#replicates
Q3InfDiff<-abs(eset.expression[,14]-eset.expression[,15])
hist(Q3InfDiff)
Q4InfDiff<-abs(eset.expression[,12]-eset.expression[,13])
hist(Q4InfDiff)
V6InfDiff<-abs(eset.expression[,10]-eset.expression[,11])
hist(V6InfDiff)
V7InfDiff<-abs(eset.expression[,8]-eset.expression[,9])
hist(V7InfDiff)
Q3MDiff<-abs(eset.expression[,6]-eset.expression[,7])
hist(Q3MDiff)
Q4MDiff<-abs(eset.expression[,4]-eset.expression[,5])
hist(Q4MDiff)
V6MDiff<-abs(eset.expression[,2]-eset.expression[,3])
hist(V6MDiff)
V7MDiff<-abs(eset.expression[,1]-eset.expression[,16])

#annotation
library(mogene20sttranscriptcluster.db)
mogene20sttranscriptcluster()
Annot <- data.frame(ENTREZID=sapply(contents(mogene20sttranscriptclusterENTREZID), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "))
Annot$ENTREZID<-as.numeric(as.character(Annot$ENTREZID))
all <- merge(Annot, eset.expression, by.x=0, by.y=0, all=T)
annot.data<-na.omit(all)

#MA plots and highlight genes
final.data<-annot.data
rownames(final.data)<-annot.data[,1]
final.data[,1]<-NULL

genes<-read.table(file="genes.txt", header=T)

library(affy)
#library(calibrate)
#y<-final.data[,c("ENTREZID","SYMBOL","LePen1_V7Int2_(MoGene-2_0-st).CEL","LePen9_V7M2_(MoGene-2_0-st).CEL")]
#MAdata<-data.frame(GeneID=y$ENTREZID,GeneSymbol=y$SYMBOL,A=rowMeans(y[,7:8]), M=y[,7]-y[,8])
y<-final.data[,c("ENTREZID","SYMBOL", ".CEL", ".CEL", ".CEL",  ".CEL")]
y$Q3_inf_avg<-rowMeans(y[,c(3,4)], na.rm=TRUE)
y$V7_inf_avg<-rowMeans(y[,c(5,6)], na.rm=TRUE)
MAdata<-data.frame(GeneID=y$ENTREZID,GeneSymbol=y$SYMBOL,A=rowMeans(y[,7:8]), M=y[,8]-y[,7])

ma.plot(MAdata[,"A"], MAdata[,"M"], cex=1, xlim=c(2,14), ylim=c(-4,4))
title("Q3_inf_avg vs V7_inf_avg")
with(MAdata[MAdata$GeneID %in% genes$ENTREZID, ], points(A, M, pch=20, col="indianred"))
with(MAdata[MAdata$GeneID %in% genes$ENTREZID & MAdata$M>2 ], textxy(A, M, labs=GeneSymbol, cex=0.8))

with(MAdata[MAdata$GeneID %in% genes$ENTREZID & MAdata$M>2, ], textxy(A, M, labs=GeneSymbol, cex=0.8))
#with(subset(MAdata, MAdata[,"M"]>4 ), points(A, M, pch=20, col="red"))
#with(subset(MAdata, MAdata[,"M"]>4 ), textxy(A, M, labs=Gene, cex=0.8))



