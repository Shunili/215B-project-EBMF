par(mar=c(2,2,2,2))
cv.out <- PMD.cv(dna[,chrom==1],type="ordered",chrom=chrom[chrom==1],
                 nuc=nuc[chrom==1],
                 sumabsus=seq(1, sqrt(nrow(dna)), len=15))
print(cv.out)
plot(cv.out)
out <- PMD(dna[,chrom==1],type="ordered",
           sumabsu=cv.out$bestsumabsu,chrom=chrom[chrom==1],K=1,v=cv.out$v.init,
           cnames=paste("Pos",sep="",
                        nuc[chrom==1]), rnames=paste("Sample", sep=" ", 1:nrow(dna)))
print(out, verbose=TRUE)
# Which samples actually have that region of gain/loss?
par(mfrow=c(3,1))
par(mar=c(2,2,2,2))
PlotCGH(dna[which.min(out$u[,1]),chrom==1],chrom=chrom[chrom==1],
        main=paste(paste(paste("Sample ", sep="", which.min(out$u[,1])),
                         sep="; u=", round(min(out$u[,1]),3))),nuc=nuc[chrom==1])
PlotCGH(dna[88,chrom==1], chrom=chrom[chrom==1],
        main=paste("Sample 88; u=", sep="", round(out$u[88,1],3)),
        nuc=nuc[chrom==1])
PlotCGH(out$v[,1],chrom=chrom[chrom==1], main="V",nuc=nuc[chrom==1])



PlotCGH(dna[,1], chrom=chrom, main="Sample 1", nuc=nuc)
