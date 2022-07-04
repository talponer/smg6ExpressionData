######################################################################
### Some functions
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

tpm <- function(counts, len, library=NA, sFactor=FALSE) {
    x <- counts/len
    if (sFactor)
        return(colSums(x, na.rm=T))
    if (!is.na(library[1]))
        return(t(t(x)*1e6/library))
    else
        return(t(t(x)*1e6/colSums(counts, na.rm=T)))
}

rpkm <- function(counts, len, library=NA) {
    if (!is.na(library[1]))
        x <- t(t(counts)/(library[colnames(counts)]/1e6))
    else
        x <- t(t(counts)/(colSums(counts, na.rm=T)/1e6))
    return(x/len)
}

## rpkm2 <- function(counts, len, library=NA) {
##     if (!is.na(library[1])){
##         x <- counts
##         for (I in colnames(counts)){
##             x[,I] <- counts[,I]/(library[I]/1e6)
##         }
##     }else{
##         m <- colSums(counts, na.rm=T)
##         x <- counts
##         for (I in colnames(counts)){
##             x[,I] <- counts[,I]/(m[I]/1e6)
##         }
##     }
##     return(x/len)
## }

tableSubSample <- function(table, n, size, column){
    res <- 1:n
    for (I in 1:n){
        sub <- sample(1:dim(table)[1], size = size)
        x[I] <- mean(table[sub, column])
    }
    return(x)
}

maxCoverage <- function(x){
    return(max(x[,3]))
}

genomeCoverage <- function(id, geneModels, assembly, ylim, dir){
    I <- unique(geneModels[geneModels$gene == id, 'symbol'])
    pcModels <- genesModels.plt[genesModels.plt$type == 'protein_coding',]
    nmdModels <- genesModels.plt[genesModels.plt$type == 'nonsense_mediated_decay',]
    pstart <- min(genesModels.plt[genesModels.plt[,'symbol']==I, 'start'])
    pstrand <- unique(genesModels.plt[genesModels.plt[,'symbol']==I, 'strand'])
    pend <- max(genesModels.plt[genesModels.plt[,'symbol']==I, 'end'])
    dist <- (pend - pstart)*0.2
    chr <- unique(genesModels.plt[genesModels.plt[,'symbol']==I, 'chromosome'])
    fname <- paste(dir, '/', id, '_', I, '_smg6_genomePlot.pdf',
                   sep='')
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = assembly, chromosome = chr)
    gmtrack <- GeneRegionTrack(pcModels, genome = assembly,
                               chromosome = chr,
                               name = "Protein Cod.",
                               transcriptAnnotation="symbol",
                               collapseTranscripts=FALSE)
    nmdtrack <- GeneRegionTrack(nmdModels, col='blue', fill='blue',
                                genome = assembly, chromosome = chr,
                                background.title = "navyblue",
                                name = "NMD",
                                transcriptAnnotation="symbol",
                                collapseTranscripts=FALSE)
    ## Data tracks Mutant:
    if (pstrand == "+"){
        spl <- 'fwd'
    }else{
        spl <- 'rev'
    }
    smg6.ZT0.file <- paste("../bigWig/MUT.0", spl, "bw", sep='.')
    smg6.ZT0 <- DataTrack(range = smg6.ZT0.file, genome = assembly,
                          type = "l", chromosome = chr,
                          name = "smg6 ZT0", ylim=ylim,
                          background.title = "brown", type='polygon',
                          col='brown', fill='salmon')
    smg6.ZT4.file <- paste("../bigWig/MUT.4", spl, "bw", sep='.')
    smg6.ZT4 <- DataTrack(range = smg6.ZT4.file, genome = assembly,
                          type = "l", chromosome = chr,
                          name = "smg6 ZT4", ylim=ylim,
                          background.title = "brown", type='polygon',
                          col='brown', fill='salmon')
    smg6.ZT8.file <- paste("../bigWig/MUT.8", spl, "bw", sep='.')
    smg6.ZT8 <- DataTrack(range = smg6.ZT8.file, genome = assembly,
                          type = "l", chromosome = chr,
                          name = "smg6 ZT8", ylim=ylim,
                          background.title = "brown", type='polygon',
                          col='brown', fill='salmon')
    smg6.ZT12.file <- paste("../bigWig/MUT.12", spl, "bw", sep='.')
    smg6.ZT12 <- DataTrack(range = smg6.ZT12.file, genome = assembly,
                           type = "l", chromosome = chr,
                           name = "smg6 ZT12", ylim=ylim,
                           background.title = "brown", type='polygon',
                           col='brown', fill='salmon')
    smg6.ZT16.file <- paste("../bigWig/MUT.16", spl, "bw", sep='.')
    smg6.ZT16 <- DataTrack(range = smg6.ZT16.file, genome = assembly,
                           type = "l", chromosome = chr,
                           name = "smg6 ZT16", ylim=ylim,
                           background.title = "brown", type='polygon',
                           col='brown', fill='salmon')
    smg6.ZT20.file <- paste("../bigWig/MUT.20", spl, "bw", sep='.')
    smg6.ZT20 <- DataTrack(range = smg6.ZT20.file, genome = assembly,
                           type = "l", chromosome = chr,
                           name = "smg6 ZT20", ylim=ylim,
                           background.title = "brown", type='polygon',
                           col='brown', fill='salmon')
    pdf(fname, width=10, height=15)
    plotTracks(list(itrack, gtrack, smg6.ZT0, smg6.ZT4, smg6.ZT8,
                    smg6.ZT12, smg6.ZT16, smg6.ZT20, gmtrack, nmdtrack),
               from=pstart-dist, to=pend+dist, chromosome=chr)
    dev.off()
    ## Data tracks WT:
    if (pstrand == "+"){
        spl <- 'fwd'
    }else{
        spl <- 'rev'
    }
    wt.ZT0.file <- paste("../bigWig/WT.0", spl, "bw", sep='.')
    wt.ZT0 <- DataTrack(range = wt.ZT0.file, genome = assembly,
                          type = "l", chromosome = chr,
                          name = "WT ZT0", ylim=ylim,
                          background.title = "blue", type='polygon',
                          col='blue', fill='skyblue')
    wt.ZT4.file <- paste("../bigWig/WT.4", spl, "bw", sep='.')
    wt.ZT4 <- DataTrack(range = wt.ZT4.file, genome = assembly,
                          type = "l", chromosome = chr,
                          name = "WT ZT4", ylim=ylim,
                          background.title = "blue", type='polygon',
                          col='blue', fill='skyblue')
    wt.ZT8.file <- paste("../bigWig/WT.8", spl, "bw", sep='.')
    wt.ZT8 <- DataTrack(range = wt.ZT8.file, genome = assembly,
                        type = "l", chromosome = chr, name = "WT ZT8",
                        ylim=ylim, background.title = "blue",
                        type='polygon', col='blue', fill='skyblue')
    wt.ZT12.file <- paste("../bigWig/WT.12", spl, "bw", sep='.')
    wt.ZT12 <- DataTrack(range = wt.ZT12.file, genome = assembly,
                         type = "l", chromosome = chr,
                         name = "WT ZT12", ylim=ylim,
                         background.title = "blue", type='polygon',
                         col='blue', fill='skyblue')
    wt.ZT16.file <- paste("../bigWig/WT.16", spl, "bw", sep='.')
    wt.ZT16 <- DataTrack(range = wt.ZT16.file, genome = assembly,
                         type = "l", chromosome = chr,
                         name = "WT ZT16", ylim=ylim,
                         background.title = "blue", type='polygon',
                         col='blue', fill='skyblue')
    wt.ZT20.file <- paste("../bigWig/WT.20", spl, "bw", sep='.')
    wt.ZT20 <- DataTrack(range = wt.ZT20.file, genome = assembly,
                         type = "l", chromosome = chr,
                         name = "WT ZT20", ylim=ylim,
                         background.title = "blue", type='polygon',
                         col='blue', fill='skyblue')
    fname <- paste(dir, '/', id, '_', I, '_wt_genomePlot.pdf',
                   sep='')
    pdf(fname, width=10, height=15)
    plotTracks(list(itrack, gtrack, wt.ZT0, wt.ZT4, wt.ZT8,
                    wt.ZT12, wt.ZT16, wt.ZT20, gmtrack, nmdtrack),
               from=pstart-dist, to=pend+dist, chromosome=chr)
    dev.off()
}

## Plot gene coverage in a genome browser fashion
plotGeneCoverage <- function(gid, gname, sampleAnno, variable, dataDir,
                             annot=annot, normFactor=normFactor,
                             fillCol=fillCol, rylim=0, rxlim=1){
    transcriptsAnno <- annot[annot$gene_id == gid]
    ## show(transcriptsAnno)
    chrName <- as.character(seqnames(transcriptsAnno)[1])
    start.g <- start(transcriptsAnno[transcriptsAnno$type == 'gene'])
    end.g <- end(transcriptsAnno[transcriptsAnno$type == 'gene'])
    start.e.p <- start(transcriptsAnno[transcriptsAnno$type == 'exon' &
                                       transcriptsAnno$transcript_biotype == 'protein_coding'])
    end.e.p <- end(transcriptsAnno[transcriptsAnno$type == 'exon' &
                                   transcriptsAnno$transcript_biotype == 'protein_coding'])
    start.e.n <- start(transcriptsAnno[transcriptsAnno$type == 'exon' &
                                       transcriptsAnno$transcript_biotype == 'nonsense_mediated_decay'])
    end.e.n <- end(transcriptsAnno[transcriptsAnno$type == 'exon' &
                                   transcriptsAnno$transcript_biotype == 'nonsense_mediated_decay'])
    sGroups <- sort(unique(sampleAnno[,variable]))
    cSampleData <- list()
    for (S in unique(sampleAnno$SampleID)){
        cSampleData[[S]] <- read.table(paste(paste(dataDir, S, sep='/'),
                                             gid, gname, 'mpileup',
                                             sep='.'))
    }
    if (length(rxlim) == 1){
        region <- c(start.g-5000,end.g+5000)
    }else{
        region <- rxlim
    }
    if (length(rylim) == 1){
        rylim <- c(0,max(unlist(lapply(cSampleData, maxCoverage))/normFactor))
    }
    par(mar=c(0,5,1,4), mfrow=c(length(sGroups)+1, 1))
    for (G in sGroups){
        ids <- sampleAnno$SampleID[sampleAnno[,variable]==G]
        plot(cSampleData[[ids[1]]][,2],
             cSampleData[[ids[1]]][,3]/normFactor[ids[1]], type='l',
             xlim=region, col='darkgray', xlab='', ylab='Counts',
             frame.plot=F, xaxt='n',
             main=paste(variable, '-', G, sep=' '),
             ylim=rylim)
        c <- 1
        for (I in ids){
            polygon(c(0,cSampleData[[I]][,2], max(cSampleData[[I]][,2])),
                    c(0,cSampleData[[I]][,3]/normFactor[I],0), col=fillCol[c],
                    border=NA)
            c <- c+1
        }
    }
    ###
    par(mar=c(5,5,0,4))
    plot(c(min(start.e.p), max(end.e.p)), c(1,1), type='l', xlim=region,
         xlab=paste('Chromosome', chrName, sep=' '), yaxt='n', ylab='',
         frame.plot=F, ylim=c(0.2,1.15))
    rect(start.e.p, 0.9, end.e.p, 1.1, col='darkred')
    if( length(start.e.n) > 0 ){
        points(c(min(start.e.n), max(end.e.n)), c(0.7,0.7), type='l')
        rect(start.e.n, 0.6, end.e.n, 0.8, col='blue')
    }
    text(mean(min(start.e.p), max(end.e.p)), 0.3,
         labels=paste(gid, '-', gname, sep=' '))
}

### Function that plot the expression of a gene fitting an harmonic curve
plotGene <- function(x, data, title='', sub='mRNA', legend=FALSE){
    geneCounts <- plotCounts(data,
                         gene=x,
                         intgroup=c("Genotype", "ZT"),
                         returnData = TRUE)
    geneCounts$ZT <- as.integer(as.character(geneCounts$ZT))
    geneCounts$Genotype <- relevel(geneCounts$Genotype, 'MUT')
    geneCounts$count <- log2(geneCounts$count)
    mVal <- tapply(geneCounts$count,
                   paste(geneCounts$Genotype, geneCounts$ZT, sep="_"), mean,
                   na.rm=T)
    mNames <- unlist(strsplit(names(mVal), split="_"))
    fitData <- data.frame(y=c(geneCounts$count, geneCounts$count),
                          x=c(geneCounts$ZT, geneCounts$ZT+24),
                          geno = c(geneCounts$Genotype, geneCounts$Genotype))
    fit1 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)), # + cos(x*pi/24) + sin(x*pi/24),
               data=fitData[fitData$geno == 'WT',])
    fit2 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)), # + cos(x*pi/24) + sin(x*pi/24),
               data=fitData[fitData$geno == 'MUT',])
    zt <- seq(0,24,0.1)
    fData <- data.frame(x=c(zt,zt),
                        y=c(predict(fit1, newdata=list(x=zt)),
                            predict(fit2, newdata=list(x=zt))),
                        geno = factor(rep(c('WT', 'MUT'), each=length(zt)),
                                      levels=c('MUT','WT')))
    fData2 <- data.frame(x=as.numeric(mNames[seq(2,length(mNames),2)]),
                         y=mVal,
                         geno = factor(mNames[seq(1,length(mNames),2)],
                                       levels=c('MUT','WT')))
    p <- ggplot(geneCounts, aes(x=ZT, y=count, color=Genotype)) +
        ## scale_y_log10() +
        geom_point(cex = 3) +
        geom_line(data=fData2, aes(x=x, y=y, col=geno), linetype='dotted') +
        geom_line(data=fData, aes(x=x, y=y, col=geno)) +
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,6)) +
        labs(title=title,
             subtitle=sub,
             y='log2(count)',
             caption='') +
        theme_gray(14)
    if (legend){
        p+theme(legend.position=c(0.8, 0.2),
                legend.box = "horizontal",
                legend.title = element_text(size=10),
                legend.text = element_text(size=8))
    }else{
        p+theme(legend.position='none')
    }
}


## Normalisation function to adjust RNA stability values
require(MASS)
eBayesNorm <- function(exonCounts, geneCounts){
    estimate <- matrix(ncol=dim(exonCounts)[2], nrow=dim(exonCounts)[1])
    colnames(estimate) <- colnames(exonCounts)
    rownames(estimate) <- rownames(exonCounts)
    for (I in 1:dim(exonCounts)[2]){
        quant <- quantile(exonCounts[,I])
        strongPos <- which(exonCounts[,I] > quant[4])
        signal <- exonCounts[,I]
        total <- signal + geneCounts[,I] + 1
        ratio <- signal/total
        ratio <- (ratio * (length(ratio) - 1) + 0.5)/length(ratio)
        ratio.strong <- ratio[strongPos]
        ## summary(ratio.strong)
        m.norm <- fitdistr(ratio.strong, dbeta, start=list(shape1=1, shape2=10))
        alpha0 <- m.norm$estimate[1]
        beta0 <- m.norm$estimate[2]
        estimate[,I] <- (exonCounts[,I] + alpha0) / ((geneCounts[,I]-exonCounts[,I]) + alpha0 + beta0)
    }
    return(estimate)
}

zScore <- function(x){
    xm <- rowMeans(x, na.rm=T)
    xsd <- apply(x, 1, sd, na.rm=T)
    (x - xm) / xsd
}

zScoreHeat <- function(x){
    xm <- cbind(
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(1:6,13:18)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T),
        rowMeans(x[,c(7:12,19:24)], na.rm=T)
    )
    xn <- x - xm
    xsd <- apply(xn, 1, sd, na.rm=T)
    xn / xsd
}
