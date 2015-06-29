CpGid <- function(methylBase) {paste(as.character(methylBase$chr), as.numeric(methylBase$start),sep=".")}
DMCplot <- function(methylBase) {
    myDiff <- calculateDiffMeth(methylBase, slim=TRUE, weighted.mean=TRUE, num.cores=1)
    myHyper_20p <- get.methylDiff(myDiff, difference=20, qvalue=0.01, type="hyper")
    myHypo_20p <- get.methylDiff(myDiff, difference=20, qvalue=0.01, type="hypo")
    layout(matrix(c(1,1,2,2), ncol=2, byrow=TRUE), heights=c(10,4))
    par(mar=c(4, 7, 2, 5) + 0.1)
    smoothScatter(methylBase$numCs2/methylBase$coverage2*100, methylBase$numCs1/methylBase$coverage1*100, xlab="Pooled ctrl methylation (%)", ylab="Pooled test methylation (%)", main=paste("n=", nrow(methylBase), " CpGs",sep=""), col=NA, colramp=colorRampPalette(c("white", "gray90", "gray80", "gray70", "gray60", "gray50", "gray40", "gray30", "gray20", "gray10", "black")))
    lines(0:100,0:100,col="red")
    points(methylBase[CpGid(methylBase)%in%CpGid(myHyper_20p),]$numCs2/methylBase[CpGid(methylBase)%in%CpGid(myHyper_20p),]$coverage2*100, methylBase[CpGid(methylBase)%in%CpGid(myHyper_20p),]$numCs1/methylBase[CpGid(methylBase)%in%CpGid(myHyper_20p),]$coverage1*100, pch=20, col = "gold2")
    points(methylBase[CpGid(methylBase)%in%CpGid(myHypo_20p),]$numCs2/methylBase[CpGid(methylBase)%in%CpGid(myHypo_20p),]$coverage2*100, methylBase[CpGid(methylBase)%in%CpGid(myHypo_20p),]$numCs1/methylBase[CpGid(methylBase)%in%CpGid(myHypo_20p),]$coverage1*100, pch=20, col = "blue")
    par(mar=c(4, 4, 0, 0) + 0.1)
    barplot(rbind(nrow(myDiff_20p.hypo), nrow(myDiff_20p.hyper)), names="DMCs", col=c("blue", "gold2"), xlab="Number of CpGs", beside=T, horiz=T, width=0.25, xlim=c(0,max(nrow(myDiff_20p.hypo), nrow(myDiff_20p.hyper))*1.2))
    legend("topright", legend=c("Hypermethylated", "Hypomethylated"), fill=c("gold2", "blue"), bty="n")
    layout(matrix(1))
    par(mar=c(5.1,4.1,4.1,2.1))
}
