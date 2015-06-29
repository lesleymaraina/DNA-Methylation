genomicDistBarplot <- function(hyperAnn, hypoAnn, allReprAnn) {
    genomicPercent <- rbind(getTargetAnnotationStats(hyperAnn), getTargetAnnotationStats(hypoAnn), getTargetAnnotationStats(allReprAnn))
    genomicCounts <- genomicPercent*c(nrow(myDiff_20p.hyper), nrow(myDiff_20p.hypo), nrow(myDiff))/100
    barx <- barplot(genomicPercent, names=colnames(genomicPercent), beside=T, col=c("gold2", "blue", "grey50"), ylab="Percentage of CpGs", ylim=c(0,max(genomicPercent)*1.2))
    mySignificance <- NULL
    for (j in 1:ncol(genomicCounts)) {
        mySignificance_j <- c(binom.test(genomicCounts[1,j], rowSums(genomicCounts)[1], p=genomicCounts[3,j]/rowSums(genomicCounts)[3])$p.value, binom.test(genomicCounts[2,j], rowSums(genomicCounts)[2], p=genomicCounts[3,j]/rowSums(genomicCounts)[3])$p.value, binom.test(genomicCounts[3,j], rowSums(genomicCounts)[3], p=genomicCounts[3,j]/rowSums(genomicCounts)[3])$p.value)
        mySignificance <- cbind(mySignificance, mySignificance_j)
    }
    mySignificance[as.numeric(mySignificance)<0.001] <- "***"; mySignificance[as.numeric(mySignificance)<0.01] <- "**"; mySignificance[as.numeric(mySignificance)<0.05] <- "*"; mySignificance[as.numeric(mySignificance)>=0.05] <- NA
    text(barx, genomicPercent, labels=mySignificance, pos=3)
    legend("topright", legend=c("hyperDMCs", "hypoDMCs", "All Represented CpGs"),fill=c("gold2", "blue", "grey50"),bty="n", cex=0.7)
}
