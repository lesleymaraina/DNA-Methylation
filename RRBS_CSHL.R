### R code from vignette source 'RRBS_CSHL.Rnw'

###################################################
### code chunk number 1: RRBS_CSHL.Rnw:26-28
###################################################
library(methylKit)
library(biomaRt)


###################################################
### code chunk number 2: RRBS_CSHL.Rnw:36-40
###################################################
file.list <- list("methylcall_test1.txt", "methylcall_test2.txt",
"methylcall_ctrl1.txt", "methylcall_ctrl2.txt")
myobj <- read(file.list, sample.id=list("test1", "test2", "ctrl1",
"ctrl2"), assembly="hg19", treatment=c(1, 1, 0, 0), context="CpG")


###################################################
### code chunk number 3: RRBS_CSHL.Rnw:47-48
###################################################
as(myobj[[2]], "GRanges")


###################################################
### code chunk number 4: RRBS_CSHL.Rnw:56-57
###################################################
getMethylationStats(myobj[[2]], plot=FALSE, both.strands=FALSE)


###################################################
### code chunk number 5: RRBS_CSHL.Rnw:63-64
###################################################
getMethylationStats(myobj[[2]], plot=TRUE, both.strands=FALSE)


###################################################
### code chunk number 6: RRBS_CSHL.Rnw:72-73
###################################################
getCoverageStats(myobj[[2]], plot=T, both.strands=F)


###################################################
### code chunk number 7: RRBS_CSHL.Rnw:80-82
###################################################
filtered.myobj <- filterByCoverage(myobj, lo.count=10, lo.perc=NULL,
hi.count=NULL, hi.perc=99.9)


###################################################
### code chunk number 8: RRBS_CSHL.Rnw:88-91
###################################################
meth <- unite(filtered.myobj, destrand=FALSE)
head(meth)
nrow(meth)


###################################################
### code chunk number 9: RRBS_CSHL.Rnw:98-99
###################################################
getCorrelation(meth, plot=T)


###################################################
### code chunk number 10: RRBS_CSHL.Rnw:110-111
###################################################
clusterSamples(meth, dist='manhattan', method='ward', plot=TRUE)


###################################################
### code chunk number 11: RRBS_CSHL.Rnw:117-119
###################################################
hc <- clusterSamples(meth, dist='manhattan', method='ward', plot=FALSE)
summary(hc)


###################################################
### code chunk number 12: RRBS_CSHL.Rnw:125-126
###################################################
PCASamples(meth)


###################################################
### code chunk number 13: RRBS_CSHL.Rnw:133-134
###################################################
PCASamples(meth, screeplot=TRUE)


###################################################
### code chunk number 14: RRBS_CSHL.Rnw:143-148
###################################################
file.list <- list("condA_sample1.txt", "condA_sample2.txt", "condB_sample1.txt",
"condB_sample2.txt", "condC_sample1.txt", "condC_sample2.txt")
clist <- read(file.list, sample.id=list("A1", "A2", "B1", "B2", "C1", "C2"),
assembly="hg19", treatment=c(2, 2, 1, 1, 0, 0), context="CpG")
newMeth <- unite(clist)


###################################################
### code chunk number 15: RRBS_CSHL.Rnw:154-155
###################################################
clusterSamples(newMeth, dist='manhattan', method='ward')


###################################################
### code chunk number 16: RRBS_CSHL.Rnw:163-167
###################################################
myDiff <- calculateDiffMeth(meth, slim=TRUE, weighted.mean=TRUE, num.cores=1)
myDiff_20p <- get.methylDiff(myDiff, difference=20, qvalue=0.01)
head(myDiff_20p)
nrow(myDiff_20p)


###################################################
### code chunk number 17: RRBS_CSHL.Rnw:172-176
###################################################
myDiff_20p.hyper <- get.methylDiff(myDiff, difference=20, qvalue=0.01, type="hyper")
nrow(myDiff_20p.hyper)
myDiff_20p.hypo <- get.methylDiff(myDiff, difference=20, qvalue=0.01, type="hypo")
nrow(myDiff_20p.hypo)


###################################################
### code chunk number 18: RRBS_CSHL.Rnw:182-184
###################################################
pooled.obj <- pool(meth, sample.ids=c("test","control"))
head(pooled.obj)


###################################################
### code chunk number 19: RRBS_CSHL.Rnw:190-192
###################################################
source("DMCplot.R")
DMCplot(pooled.obj)


###################################################
### code chunk number 20: RRBS_CSHL.Rnw:198-199
###################################################
myDiff <- calculateDiffMeth(meth, num.cores=2)


###################################################
### code chunk number 21: RRBS_CSHL.Rnw:205-209
###################################################
gene.obj <- read.transcript.features("refseq.hg19.bed.txt")
diffAnnotate.hyper <- annotate.WithGenicParts(myDiff_20p.hyper, gene.obj)
diffAnnotate.hypo <- annotate.WithGenicParts(myDiff_20p.hypo, gene.obj)
getTargetAnnotationStats(diffAnnotate.hypo, percentage=TRUE, precedence=TRUE)


###################################################
### code chunk number 22: RRBS_CSHL.Rnw:215-217
###################################################
plotTargetAnnotation(diffAnnotate.hypo, precedence=TRUE, main="Hypo DMC
distribution")


###################################################
### code chunk number 23: RRBS_CSHL.Rnw:227-230
###################################################
diffAnnotate.allRepr <- annotate.WithGenicParts(myDiff, gene.obj)
source("genomicDistBarplot.R")
genomicDistBarplot(diffAnnotate.hyper,diffAnnotate.hypo,diffAnnotate.allRepr)


###################################################
### code chunk number 24: RRBS_CSHL.Rnw:237-243
###################################################
cpg.obj <- read.feature.flank("cpgi.hg19.bed.txt", feature.flank.name=c("CpGi",
"shores"))
diffCpGann <- annotate.WithFeature.Flank(myDiff_20p.hypo, cpg.obj$CpGi, cpg.obj$shores,
feature.name = "CpGi", flank.name = "shores")
plotTargetAnnotation(diffCpGann, col=c("springgreen3", "royalblue3", "grey50"), main=
"Hypo DMC annotation")


###################################################
### code chunk number 25: RRBS_CSHL.Rnw:250-258
###################################################
tiles <- tileMethylCounts(myobj, win.size=1000, step.size=1000)
head(tiles[[2]])
meth.DMR <- unite(tiles, destrand=TRUE)
myDiff.DMR <- calculateDiffMeth(meth.DMR, slim=TRUE, weighted.mean=TRUE, num.cores=1)
myDiff.DMR_20p.hyper <- get.methylDiff(myDiff.DMR, difference=20, qvalue=0.01,
type="hyper")
myDiff.DMR_20p.hypo <- get.methylDiff(myDiff.DMR, difference=20, qvalue=0.01,
type="hypo")


###################################################
### code chunk number 26: RRBS_CSHL.Rnw:266-269
###################################################
diffAnnotate.hypo <- annotate.WithGenicParts(myDiff.DMR_20p.hypo, gene.obj)
hypoDMR_TSS <- getAssociationWithTSS(diffAnnotate.hypo)
head(hypoDMR_TSS)


###################################################
### code chunk number 27: RRBS_CSHL.Rnw:273-275
###################################################
hist(hypoDMR_TSS$dist.to.feature, col="royalblue", xlab="Distance to nearest TSS",
breaks=40)


###################################################
### code chunk number 28: RRBS_CSHL.Rnw:281-285
###################################################

diffAnnotate.hyper <- annotate.WithGenicParts(myDiff.DMR_20p.hyper, gene.obj)
hyperDMR_TSS <- getAssociationWithTSS(diffAnnotate.hyper)
hyperDMR_genes <- unique(hyperDMR_TSS[hyperDMR_TSS$dist.to.feature<2000,]$feature.name)
hypoDMR_genes <- unique(hypoDMR_TSS[hypoDMR_TSS$dist.to.feature<2000,]$feature.name)


###################################################
### code chunk number 29: RRBS_CSHL.Rnw:290-292
###################################################
fpkm <- read.table("test_ctrl_fpkm.txt", header=T, row.names=1, sep="\t")
head(fpkm)


###################################################
### code chunk number 30: RRBS_CSHL.Rnw:300-308
###################################################
hyperDMR_key <- getBM(attributes=c('hgnc_symbol', 'refseq_mrna'), filters=
'refseq_mrna', values=hyperDMR_genes, mart=useMart("ensembl", dataset=
"hsapiens_gene_ensembl"))
head(hyperDMR_key)

hypoDMR_key <- getBM(attributes=c('hgnc_symbol', 'refseq_mrna'), filters=
'refseq_mrna', values=hypoDMR_genes, mart=useMart("ensembl", dataset=
"hsapiens_gene_ensembl"))


###################################################
### code chunk number 31: RRBS_CSHL.Rnw:314-320
###################################################
hyperDMR_symbol <- hyperDMR_key$hgnc_symbol[hyperDMR_key$hgnc_symbol%in%rownames(fpkm)]

boxplot(fpkm[hyperDMR_symbol,], names=c("ctrl", "test"), ylab="Gene expression level
(FPKM)", main=paste("Hyper DMR-associated genes (n=", length(hyperDMR_symbol),") \n
Wilcoxon p=",signif(wilcox.test(fpkm[hyperDMR_symbol,1], fpkm[hyperDMR_symbol,2],
paired=T)$p.value,3),sep=""))


###################################################
### code chunk number 32: RRBS_CSHL.Rnw:325-331
###################################################
hypoDMR_symbol <- hypoDMR_key$hgnc_symbol[hypoDMR_key$hgnc_symbol%in%rownames(fpkm)]

boxplot(fpkm[hypoDMR_symbol,], names=c("ctrl", "test"), ylab="Gene expression level
(FPKM)", main=paste("Hypo DMR-associated genes (n=", length(hypoDMR_symbol),") \n
Wilcoxon p=",signif(wilcox.test(fpkm[hypoDMR_symbol,1], fpkm[hypoDMR_symbol,2],
paired=T)$p.value,3),sep=""))


###################################################
### code chunk number 33: RRBS_CSHL.Rnw:341-343
###################################################
plot(log(fpkm[hypoDMR_symbol,]), xlab="Ctrl log(FPKM)", ylab="Test log(FPKM)", pch=20,
col="grey25", main="Hypo DMR-associated genes"); abline(a=0, b=1, col="blue")


###################################################
### code chunk number 34: RRBS_CSHL.Rnw:350-354
###################################################
plot(log(fpkm[hypoDMR_symbol,]), xlab="Ctrl log(FPKM)", ylab="Test log(FPKM)", pch=20,
col="grey25", main="Hypo DMR-associated genes"); abline(a=0, b=1, col="blue")
points(log(fpkm[hypoDMR_symbol,][fpkm[hypoDMR_symbol,2]/fpkm[hypoDMR_symbol,1]>2 &
fpkm[hypoDMR_symbol,1]>5,]), pch=20, col="red")


###################################################
### code chunk number 35: RRBS_CSHL.Rnw:361-363
###################################################
rownames(fpkm[hypoDMR_symbol,][fpkm[hypoDMR_symbol,2]/fpkm[hypoDMR_symbol,1]>2 &
fpkm[hypoDMR_symbol,1]>5,])


###################################################
### code chunk number 36: RRBS_CSHL.Rnw:369-370
###################################################
bedgraph(myDiff.DMR_20p.hypo, file="hypoDMR.bedgraph", 'meth.diff')


