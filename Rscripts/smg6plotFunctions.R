### Function that plot the expression of a gene fitting an harmonic curve
plotGene.0 <- function(x, data, title='', sub='mRNA', legend=FALSE){
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


### as plotGene.0 but change the theme:
plotGene.1 <- function(x, data, title='', sub='mRNA', legend=FALSE){
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
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,4)) +
        labs(title=title,
             subtitle=sub,
             y='log2(count)',
             caption='') +
        theme_cowplot(14)
    if (legend){
        p+theme(legend.position=c(0.8, 0.2),
                legend.box = "horizontal",
                legend.title = element_text(size=10),
                legend.text = element_text(size=8))
    }else{
        p+theme(legend.position='none')
    }
}


### as plotGene.0 but change the theme:
plotGene.2 <- function(x, data, title='', sub='mRNA', legend=FALSE,
                       fit=TRUE){
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
    fit1 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
               data=fitData[fitData$geno == 'WT',])
    fit2 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
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
        geom_point(cex = 2)
    if (fit){
        p <- p +
            geom_line(data=fData, aes(x=x, y=y, col=geno), linetype='dashed')
    }
    ## geom_line(data=fData, aes(x=x, y=y, col=geno), linetype='dashed') +
    p <- p +
        geom_line(data=fData2, aes(x=x, y=y, col=geno), size=1.3) +
        geom_point(data=fData2, aes(x=x, y=y, col=geno), size=3,
                   shape=21, fill="white") +
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,4)) +
        labs(title=title,
             subtitle=sub,
             y='log2(count)',
             caption='') +
        theme_cowplot(14) +
        background_grid()
    if (legend){
        p+theme(legend.position=c(0.8, 0.2),
                legend.box = "horizontal",
                legend.title = element_text(size=10),
                legend.text = element_text(size=8))
    }else{
        p+theme(legend.position='none')
    }
}


plotGene.3 <- function(x, data, data2, title='', legend=FALSE,
                       fit=TRUE){
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
    fit1 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
               data=fitData[fitData$geno == 'WT',])
    fit2 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
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
    ###
    geneCounts2 <- plotCounts(data2,
                         gene=x,
                         intgroup=c("Genotype", "ZT"),
                         returnData = TRUE)
    geneCounts2$ZT <- as.integer(as.character(geneCounts2$ZT))
    geneCounts2$Genotype <- relevel(geneCounts2$Genotype, 'MUT')
    geneCounts2$count <- log2(geneCounts2$count)
    mVal2 <- tapply(geneCounts2$count,
                   paste(geneCounts2$Genotype, geneCounts2$ZT, sep="_"), mean,
                   na.rm=T)
    mNames2 <- unlist(strsplit(names(mVal2), split="_"))
    fitData2 <- data.frame(y=c(geneCounts2$count, geneCounts2$count),
                          x=c(geneCounts2$ZT, geneCounts2$ZT+24),
                          geno = c(geneCounts2$Genotype, geneCounts2$Genotype))
    fit1.2 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
               data=fitData2[fitData2$geno == 'WT',])
    fit2.2 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
               data=fitData2[fitData2$geno == 'MUT',])
    fData1.2 <- data.frame(x=c(zt,zt),
                         y=c(predict(fit1.2, newdata=list(x=zt)),
                             predict(fit2.2, newdata=list(x=zt))),
                         geno = factor(rep(c('WT', 'MUT'), each=length(zt)),
                                       levels=c('MUT','WT')))
    fData2.2 <- data.frame(x=as.numeric(mNames2[seq(2,length(mNames2),2)]),
                           y=mVal2,
                           geno = factor(mNames[seq(1,length(mNames2),2)],
                                         levels=c('MUT','WT')))
    
    fData3 <- rbind(fData, fData1.2)
    fData3$type <- factor(rep(c('mRNA', 'pre-mRNA'), each=dim(fData)[1]),
                           levels=c('mRNA', 'pre-mRNA'))
    fData4 <- rbind(fData2, fData2.2)
    fData4$type <- factor(rep(c('mRNA', 'pre-mRNA'), each=dim(fData2)[1]),
                            levels=c('mRNA', 'pre-mRNA'))
    geneCounts3 <- rbind(geneCounts, geneCounts2)
    geneCounts3$type <- factor(rep(c('mRNA', 'pre-mRNA'), each=dim(geneCounts)[1]),
                              levels=c('mRNA', 'pre-mRNA'))
    ###
    p <- ggplot(geneCounts3, aes(x=ZT, y=count, color=Genotype)) +
        ## scale_y_log10() +
        geom_point(cex = 2)
    if (fit){
        p <- p +
            geom_line(data=fData3, aes(x=x, y=y, col=geno), linetype='dashed')
    }
    ## geom_line(data=fData, aes(x=x, y=y, col=geno), linetype='dashed') +
    p <- p +
        geom_line(data=fData4, aes(x=x, y=y, col=geno), size=1.3) +
        geom_point(data=fData4, aes(x=x, y=y, col=geno), size=3,
                   shape=21, fill="white") +
        facet_wrap(~ type, nrow=2, scales='free_y') +
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,4)) +
        labs(title=title,
             y='log2(count)',
             col='') +
        theme_cowplot(14) +
        panel_border() +
        background_grid()
    if (legend){
        p+theme(legend.position=c(0.8, 0.2),
                legend.box = "horizontal",
                legend.title = element_text(size=10),
                legend.text = element_text(size=8))
    }else{
        p+theme(legend.position='none')
    }
}


plotGene.4 <- function(x, data, title='', sub='mRNA', legend=FALSE,
                       fit=TRUE, ymax){
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
    fit1 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
               data=fitData[fitData$geno == 'WT',])
    fit2 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
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
        geom_point(cex = 2) +
        ylim(NA, ymax)
    if (fit){
        p <- p +
            geom_line(data=fData, aes(x=x, y=y, col=geno), linetype='dashed')
    }
    ## geom_line(data=fData, aes(x=x, y=y, col=geno), linetype='dashed') +
    p <- p +
        geom_line(data=fData2, aes(x=x, y=y, col=geno), size=1.3) +
        geom_point(data=fData2, aes(x=x, y=y, col=geno), size=3,
                   shape=21, fill="white") +
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,4)) +
        labs(title=title,
             subtitle=sub,
             y='log2(count)',
             caption='') +
        theme_cowplot(14) +
        background_grid()
    if (legend){
        p+theme(legend.position=c(0.85, 0.2),
                legend.box = "horizontal",
                legend.title = element_text(size=10),
                legend.text = element_text(size=8))
    }else{
        p+theme(legend.position='none')
    }
}


plotGene.5 <- function(filename, x, data1, data2, title='', legend=FALSE,
                       fit=TRUE){
    r1 <- log2(max(counts(data1[x,], norm=T))) -
        log2(min(counts(data1[x,], norm=T)))
    r2 <- log2(max(counts(data2[x,], norm=T))) -
        log2(min(counts(data2[x,], norm=T)))
    if(is.na(r1)){
        ymax2 <- log2(max(counts(data2[x,], norm=T)))
        ymax1 <- r2
    }else if(is.na(r2)){
        ymax1 <- log2(max(counts(data1[x,], norm=T)))
        ymax2 <- r1
    }else if(r1 > r2){
        ymax1 <- log2(max(counts(data1[x,], norm=T)))
        ymax2 <- log2(min(counts(data2[x,], norm=T))) + r1
    }else{
        ymax2 <- log2(max(counts(data2[x,], norm=T)))
        ymax1 <- log2(min(counts(data1[x,], norm=T))) + r2
    }
    p1 <- plotGene.4(x=x, data=data1, title=title, sub='mRNA',
                     legend=FALSE, fit=fit, ymax=ymax1)
    p2 <- plotGene.4(x=x, data=data2, title='', sub='pre-mRNA',
                     legend=legend, fit=fit, ymax=ymax2)
    ggarrange(p1, p2, ncol=1, nrow=2) %>%
        ggexport(filename = filename, height=9, width=5)
}




plotGeneTpm <- function(x, data, anno, title='', sub='mRNA', legend=FALSE,
                       fit=TRUE, ymax){
    geneCounts <- data.frame(anno, count = data[x,])
    geneCounts$ZT <- as.integer(as.character(geneCounts$ZT))
    geneCounts$Genotype <- relevel(geneCounts$Genotype, 'MUT')
    mVal <- tapply(geneCounts$count,
                   paste(geneCounts$Genotype, geneCounts$ZT, sep="_"), mean,
                   na.rm=T)
    mNames <- unlist(strsplit(names(mVal), split="_"))
    fitData <- data.frame(y=c(geneCounts$count, geneCounts$count),
                          x=c(geneCounts$ZT, geneCounts$ZT+24),
                          geno = c(geneCounts$Genotype, geneCounts$Genotype))
    fit1 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
               data=fitData[fitData$geno == 'WT',])
    fit2 <- lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)),
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
        geom_point(cex = 2) +
        ylim(0, NA)
    if (fit){
        p <- p +
            geom_line(data=fData, aes(x=x, y=y, col=geno), linetype='dashed')
    }
    ## geom_line(data=fData, aes(x=x, y=y, col=geno), linetype='dashed') +
    p <- p +
        geom_line(data=fData2, aes(x=x, y=y, col=geno), size=1.3) +
        geom_point(data=fData2, aes(x=x, y=y, col=geno), size=3,
                   shape=21, fill="white") +
        scale_x_continuous(limits=c(0,24), breaks=seq(0,24,4)) +
        labs(title=title,
             subtitle=sub,
             y='RPKM',
             caption='') +
        theme_cowplot(14) +
        background_grid()
    if (legend){
        p+theme(legend.position=c(0.85, 0.9),
                legend.box = "horizontal",
                legend.title = element_text(size=10),
                legend.text = element_text(size=8))
    }else{
        p+theme(legend.position='none')
    }
}


plotTpm <- function(filename, x, data1, data2, anno, title='',
                    legend=FALSE, fit=TRUE){
    r1 <- max(data1[x,]) - min(data1[x,])
    r2 <- max(data2[x,]) - min(data2[x,])
    ymax <- max(c(data1[x,], data2[x,]), na.rm=T)
    p1 <- plotGeneTpm(x=x, data=data1, anno=anno, title=title,
                      sub='mRNA', legend=FALSE, fit=fit, ymax=ymax)
    p2 <- plotGeneTpm(x=x, data=data2, anno=anno, title='',
                      sub='pre-mRNA', legend=legend, fit=fit,
                      ymax=ymax)
    ggarrange(p1, p2, ncol=1, nrow=2) %>%
        ggexport(filename = filename, height=9, width=5)
}
