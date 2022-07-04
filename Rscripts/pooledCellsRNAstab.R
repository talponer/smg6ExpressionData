######################################################################
### Firs analysis on the SMG6 RNA stability data.
######################################################################

## packrat::init("./", infer.dependencies = FALSE)

library(tximport)
library(DESeq2)
library(rtracklayer) # to import GTF file
library(lattice)
library(ReactomePA)
library(clusterProfiler)
library(biomaRt)
library(org.Mm.eg.db)
library(ggplot2)
library(cowplot)
library(ggfortify)
library(plyr)
library(preprocessCore)
library(MASS)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(Gviz)
library(ggridges)
library(emojifont)
library(pheatmap)

source('smg6functions.R')


annot <- import("/mnt/data/databases/mouse/gtf/Mmusculus.GRCm38.91.gtf")
genesAnno <- annot[annot$type == 'gene']
genesSel <- genesAnno[genesAnno$gene_biotype %in% c("protein_coding")]

directory <- '../genesCounts/'

######################################################################
## Read in data at the gene level
######################################################################
## samples annotations:
sampleGeneAnno <- read.table('../samples.anno', header=T, sep="\t")
sampleGeneAnno$Genotype <- factor(sampleGeneAnno$Genotype,
                                  levels=c('+/+', 'f/+', 'f/f'))
levels(sampleGeneAnno$Genotype) <- c('ww', 'fw', 'ff')
sampleGeneAnno$Treatment <- as.factor(sampleGeneAnno$Treatment)
sampleGeneAnno$ID <- as.character(sampleGeneAnno$ID)
sampleGeneAnno$FileName <- paste(sampleGeneAnno$ID, 'gene.htseq-reverse.dat', sep='.')
str(sampleGeneAnno)

## Start with the analysis at the gene level:
ddsGeneLevel <- DESeqDataSetFromHTSeqCount(sampleTable = sampleGeneAnno,
                                           directory = directory,
                                           design = ~ -1 + Genotype : Treatment)
ddsGeneLevel <- DESeq(ddsGeneLevel)

geneCounts <- counts(ddsGeneLevel)

######################################################################
### Read data at the mRNA level (exon) and calculate pre-mRNA (intron)
### reads
######################################################################
sampleExonAnno <- sampleGeneAnno
sampleExonAnno$FileName <- paste(sampleExonAnno$ID,
                                 'exon.htseqReverse.dat', sep='.')

ddsExonLevel <- DESeqDataSetFromHTSeqCount(sampleTable = sampleExonAnno,
                                           directory = directory,
                                           design = ~ -1 + Genotype : Treatment)
ddsExonLevel <- DESeq(ddsExonLevel)

exonCounts <- counts(ddsExonLevel)
head(exonCounts)

geneExonCounts <- exonCounts
geneIntronCounts <- geneCounts - geneExonCounts
## head(geneExonCounts)
## head(geneIntronCounts)
## head(selectedExons)

geneExonCounts['ENSMUSG00000037742',]
geneCounts['ENSMUSG00000037742',]


######################################################################
### Calculate RPKMs on the mRNA and pre-mRNA reads (from total reads)
######################################################################

## geneLength <- read.table('../annotations/geneLength.txt', header=F, row.names=1)
## exonLength <- read.table('../annotations/exonLength.txt', header=F)
## geneExonLength <- tapply(exonLength[,2], exonLength[,1], sum)
## geneIntronLength <- geneLength[names(geneExonLength),1] - geneExonLength
## range(geneExonLength)

exonIntronLength <- read.table('../annotations/2021-12-07_exonAndIntronLengths.txt',
                               header=T)
geneExonLength <- exonIntronLength[,'ExonLength']
geneIntronLength <- exonIntronLength[,'IntronLength']
names(geneExonLength) <- names(geneIntronLength) <- exonIntronLength$GeneID


## select genes that have introns
geneExonLengthClean <- geneExonLength[which(geneIntronLength > 0)]
geneIntronLengthClean <- geneIntronLength[which(geneIntronLength > 0)]
geneExonCountsGood <- geneExonCounts[rownames(geneExonCounts) %in%
                                     names(geneExonLengthClean),]
geneIntronCountsGood <- geneIntronCounts[rownames(geneIntronCounts) %in%
                                         names(geneIntronLengthClean),]

## Restrict analysis to highly expressed genes that have intron counts
goodGenes <- rownames(geneExonCountsGood)[which(apply(geneExonCountsGood, 1, min) > 6 &
                                                apply(geneIntronCountsGood, 1, min) > 0)]

gExon <- geneExonCountsGood[goodGenes,]
gIntron <- geneIntronCountsGood[goodGenes,]

#########################################################################
### Evaluate RNA stability using Bayesian inference (calculated on
### strongly expressed genes)
strongGenes <- names(which(apply(gExon, 1, min) > 50))
signalRatio <- geneExonCounts[strongGenes,] / (geneCounts[strongGenes,]+1)
                                  
exonVsIntronNorm <- gExon
for (I in 1:dim(signalRatio)[2]){ 
    m.s1 <- fitdistr(signalRatio[,I], dbeta, start=list(shape1=1, shape2=10))
    alpha0.s1 <- m.s1$estimate[1]
    beta0.s1 <- m.s1$estimate[2]
    exonVsIntronNorm[,I] <- (gExon[,I] + alpha0.s1) / (gIntron[,I] + alpha0.s1 + beta0.s1)
}

pdf('figures/2021-12-01/Bayesian_norm.pdf')
par(mfrow=c(1,1))
plot(density(signalRatio[,20]), frame.plot=F)
points(density(rbeta(length(signalRatio), m.s1$estimate[1], m.s1$estimate[2])), type="l", col="red")
legend("topleft", legend=c("Prior", "Estimate"), lwd=1, col=c("black","red"), bty="n")
dev.off()

summary(exonVsIntronNorm)
exonVsIntron <- log(exonVsIntronNorm)
## exonVsIntron <- log(gExon/gIntron)
dim(exonVsIntron)



## Get the averages:
exonVsIntron.ff0.1 <- rowMeans(exonVsIntron[,1:3], na.rm=T)
exonVsIntron.ff0.2 <- rowMeans(exonVsIntron[,11:13], na.rm=T)
exonVsIntron.ff1.1 <- rowMeans(exonVsIntron[,6:8], na.rm=T)
exonVsIntron.ff1.2 <- rowMeans(exonVsIntron[,16:18], na.rm=T)
exonVsIntron.wt0 <- rowMeans(exonVsIntron[,c(4,14)], na.rm=T)
exonVsIntron.wt1 <- rowMeans(exonVsIntron[,c(9,19)], na.rm=T)
exonVsIntron.he0 <- rowMeans(exonVsIntron[,c(5,15)], na.rm=T)
exonVsIntron.he1 <- rowMeans(exonVsIntron[,c(10,20)], na.rm=T)

summary(exonVsIntron.ff0.1)
summary(exonVsIntron.ff0.2)
summary(exonVsIntron.ff1.1)
summary(exonVsIntron.ff1.2)
summary(exonVsIntron.wt0)
summary(exonVsIntron.wt1)
summary(exonVsIntron.he0)
summary(exonVsIntron.he1)

## Get the replica averages and plot them
exonVsIntron.ff0 <- (exonVsIntron.ff0.1 + exonVsIntron.ff0.2) / 2
exonVsIntron.ff1 <- (exonVsIntron.ff1.1 + exonVsIntron.ff1.2) / 2

stabData <- data.frame(stab = c(exonVsIntron.ff1, exonVsIntron.ff0,
                                exonVsIntron.wt0, exonVsIntron.wt1),
                       Sample = c(rep(c('flox/flox treated',
                                        'flox/flox untreated',
                                        'WT untreated', 'WT treated'),
                                      each = length(exonVsIntron.ff1))))

sp1 <- ggplot(stabData, aes(x=stab, fill=Sample)) +
    geom_density(alpha=0.5) +
    labs(title='Genotype effects on RNA stability',
         x='log2(mRNA / pre-mRNA)') +
    theme_gray(14)
ggsave('figures/2021-12-01/RS_densities.pdf', sp1, width=9)

sp1 <- ggplot(stabData, aes(x=Sample, y=stab, fill=Sample)) +
    geom_boxplot() +
    labs(title='Genotype effects on RNA stability',
         y='log2(mRNA / pre-mRNA)',
         x='') +
    theme_gray(14)
ggsave('figures/2021-12-01/RS_boxplots.pdf', sp1, width=9)


ecdf.ff1 <- ecdf(exonVsIntron.ff1)

pdf('figures/2021-12-01/RS_ecdf.pdf')
plot(ecdf.ff1, main='RNA stability cumulative sum',
     xlab='log2(mRNA / pre-mRNA)')
legend(x=-4, y=0.99, legend=c('flox/flox treated', 'flox/flox untreated',
                           'WT untreated', 'WT treated'), lty=1,
       col=c('black', 'red', 'blue', 'orange'), bty='n')
plot(ecdf(exonVsIntron.ff0), add=TRUE, col='red')
plot(ecdf(exonVsIntron.wt0), add=TRUE, col='blue')
plot(ecdf(exonVsIntron.wt1), add=TRUE, col='orange')
dev.off()

stabDataCS <- data.frame(
    stab = c(sort(exonVsIntron.ff1), sort(exonVsIntron.ff0),
             sort(exonVsIntron.wt0), sort(exonVsIntron.wt1)),
    cumsum = c(cumsum(abs(sort(exonVsIntron.ff1))),
             cumsum(abs(sort(exonVsIntron.ff0))),
             cumsum(abs(sort(exonVsIntron.wt0))),
             cumsum(abs(sort(exonVsIntron.wt1)))),
    Sample = c(rep(c('flox/flox treated',
                     'flox/flox untreated',
                     'WT untreated', 'WT treated'),
                   each = length(exonVsIntron.ff1)))
)

head(stabDataCS)

maxStab <- tapply(stabDataCS$cumsum, stabDataCS$Sample, max)
stabDataCS$norm <- stabDataCS$cumsum / maxStab[stabDataCS$Sample]


sp1 <- ggplot(stabDataCS, aes(x=stab, col=Sample)) +
    stat_ecdf() +
    labs(title='Genotype effects on RNA stability',
         subtitle='Cumulative sum',
         x='log2(mRNA / pre-mRNA)',
         y='Proportion',
         col='') +
    theme_cowplot(14) +
    theme(legend.position=c(0.05,0.9))
ggsave('figures/2021-12-01/RS_cumulativeSum.pdf', sp1, width=6, height=6)


### Test the distribution with a Kolmogorov-Smirnov test:

ff1vsff0.ks <- ks.test(exonVsIntron.ff1,
                       exonVsIntron.ff0)
ff1vsff0.ks

ff1vsWw1.ks <- ks.test(exonVsIntron.ff1,
                       exonVsIntron.wt1)
ff1vsWw1.ks


#######################################################################
### Calculate RNA-stability using DESeq2 (with shrinking)
allData <- cbind(gExon, gIntron)
colnames(allData) <- c(
    paste(colnames(geneExonCountsGood), 'mRNA', sep='.'),
    paste(colnames(geneExonCountsGood), 'premRNA', sep='.')
)
head(allData)

allSamplesAnno <- rbind(sampleExonAnno, sampleGeneAnno)
allSamplesAnno$ID <- c(
    paste(sampleExonAnno$ID, 'mRNA', sep='.'),
    paste(sampleGeneAnno$ID, 'premRNA', sep='.')
)
allSamplesAnno$Type <- factor(rep(c('mRNA', 'premRNA'),
                                  each=dim(sampleExonAnno)[1]),
                              levels=c('premRNA', 'mRNA'))
allSamplesAnno$Sample <- paste(
    str_replace(allSamplesAnno$Genotype, "/", ""),
    allSamplesAnno$Treatment, sep='.')
allSamplesAnno$Sample <- gsub("+", "w", allSamplesAnno$Sample,
                              fixed=TRUE)
rownames(allSamplesAnno) <- colnames(allData)
head(allSamplesAnno)

allDds <- DESeqDataSetFromMatrix(countData = allData,
                                 colData = allSamplesAnno,
                                 design = ~ Sample + Type:Sample)
keep <- rowSums(counts(allDds) >= 10) >= 3
allDds <- allDds[keep,]
### do not normalize library total counts (it is very important,
### otherwise pre-mRNA reads will be artificially increased). It is
### not a problem, later there will be a ratio between mRNA and
### pre-mRNA reads from the same library
sizeFactors(allDds) <- c(sizeFactors(ddsExonLevel),
                         sizeFactors(ddsGeneLevel))
allDds <- DESeq(allDds, fitType='local')

resultsNames(allDds)

rnaStab.ff1 <- results(allDds, name='Sampleff.1.TypemRNA')
rnaStab.ff0 <- results(allDds, name='Sampleff.0.TypemRNA')
rnaStab.ww0 <- results(allDds, name='Sampleww.0.TypemRNA')
rnaStab.ww1 <- results(allDds, name='Sampleww.1.TypemRNA')


### Plot them:
nmdGenes <- read.table('../annotations/ensembl_nmd_annotated_genes',
                       as.is=T)[,1]

rnaStabData <- data.frame(stab = c(rnaStab.ff1$log2FoldChange,
                                   rnaStab.ff0$log2FoldChange,
                                   rnaStab.ww1$log2FoldChange,
                                   rnaStab.ww0$log2FoldChange),
                          Sample = c(rep(c('flox/flox treated',
                                           'flox/flox untreated',
                                           'WT treated',
                                           'WT untreated'),
                                         each = dim(rnaStab.ff1)[1])))

geneType <- factor(rep('Other', dim(rnaStab.ff1)[1]), levels=c('Other', 'NMD'))
names(geneType) <- rownames(rnaStab.ff1)
geneType[rownames(rnaStab.ff1) %in% nmdGenes] <- 'NMD'
rnaStabData$Type <- rep(geneType, 4)
summary(rnaStabData)
head(rnaStabData)

sp1 <- ggplot(rnaStabData, aes(x=stab, fill=Sample)) +
    geom_density(alpha=0.5) +
    labs(title='Genotype effects on RNA stability',
         x='log2(mRNA / pre-mRNA)') +
    theme_gray(14)
ggsave('figures/2021-12-01/RS-DESeq2_densities.pdf', sp1, width=9)


sp1 <- ggplot(rnaStabData, aes(x=stab, fill=Sample)) +
    geom_density(alpha=0.5) +
    facet_wrap(~ Type) +
    labs(title='Genotype effects on RNA stability',
         x='log2(mRNA / pre-mRNA)') +
    theme_gray(14)
ggsave('figures/2021-12-01/RS-DESeq2_densities_types.pdf', sp1, width=14)


sp1 <- ggplot(rnaStabData, aes(x=stab, fill=Type)) +
    geom_density(alpha=0.5) +
    facet_wrap(~ Sample) +
    labs(title='Genotype effects on RNA stability',
         x='log2(mRNA / pre-mRNA)') +
    theme_gray(14)
ggsave('figures/2021-12-01/RS-DESeq2_densities_sample.pdf', sp1,
       width=12, height=10)




### This is preatty cool, I can now get genes that have different RNA
### stability between mutant and control

rnaStab.ff1vsWw0 <- results(allDds,
                            contrast=list('Sampleff.1.TypemRNA',
                                          'Sampleww.0.TypemRNA'))

rnaStab.ff1vsWw1 <- results(allDds,
                            contrast=list('Sampleff.1.TypemRNA',
                                          'Sampleww.1.TypemRNA'))

rnaStab.ff1vsFf0 <- results(allDds,
                            contrast=list('Sampleff.1.TypemRNA',
                                          'Sampleff.0.TypemRNA'))
summary(rnaStab.ff1vsFf0)
summary(rnaStab.ff1vsWw1)

summary(rnaStab.ff1vsWw1$pvalue < 0.05)

## plot the difference as well
signGenes <- rownames(rnaStab.ff1vsFf0)
annoData <- mcols(annot)
colnames(annoData)

type <- vector(length=length(signGenes))
gName <- vector(length=length(signGenes))
names(type) <- names(gName) <- signGenes
for (I in signGenes){
    eIndex <- which(annoData$gene_id == I)
    gName[I] <- annoData$gene_name[eIndex[1]]
    for (K in annoData$transcript_biotype[eIndex]){
        if (is.na(K)){
            next
        } else if (K == 'protein_coding' ){
            type[I] <- K
            next
        } else if ( K == 'nonsense_mediated_decay' ){
            type[I] <- K
            break
        } else if (K == 'retained_intron'){
            type[I] <- K
            break
        } else {
            type[I] <- 'Other'
        }
    }
}

typef <- factor(type, levels=c('nonsense_mediated_decay',
                               'retained_intron', 'Other',
                               'protein_coding'))
summary(typef)
levels(typef) <- c('NMD','Retained Intron', 'Other', 'Protein Coding')


rnaStab.ff1vsWw1$Type <- typef
data.ff1vsWw1 <- as.data.frame(rnaStab.ff1vsWw1)
class(data.ff1vsWw0)

head(data.ff1vsFf0)

sp1 <- ggplot(data.ff1vsWw1, aes(x=log2FoldChange, fill=Type)) +
    geom_density(alpha=0.5) +
    geom_vline(xintercept=0, linetype='dashed') +
    xlim(-1, 1.5) +
    labs(title='Smg6 affects RNA stability in cells',
         subtitle='Genes stratified by transcipt types',
         caption='smg6: flox/flox treated; Control: WT treated',
         fill='',
         x='\u0394 RNA stability') +
    theme_cowplot(14) +
    theme(legend.position=c(0.02, 0.9))
ggsave('figures/2021-12-01/RS-DESeq2_smg6VsWt.pdf', sp1,
       width=9, height=6)

sp1 <- ggplot(data.ff1vsWw1, aes(x=log2FoldChange, col=Type)) +
    stat_ecdf() +
    ## geom_vline(xintercept=0, linetype='dashed') +
    ## xlim(-1, 1.5) +
    labs(title='Smg6 affects RNA stability in cells',
         subtitle='Genes stratified by transcipt types',
         caption='smg6: flox/flox treated; Control: WT treated',
         col='',
         x='\u0394 RNA stability',
         y='Cumulative Sum') +
    theme_cowplot(14) +
    theme(legend.position=c(0.02, 0.9))
ggsave('figures/2021-12-01/RS-DESeq2_smg6VsWt_ecdf.pdf', sp1,
       width=6, height=6)


sp1 <- ggplot(data.ff1vsWw1, aes(x=log2FoldChange, y=Type, fill=Type)) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_density_ridges(alpha=0.5, scale=1.5) +
    xlim(-1, 1.5) +
    labs(title='Smg6 affects RNA stability in cells',
         subtitle='Genes stratified by transcipt types',
         caption='smg6: flox/flox treated; Control: WT treated',
         x='\u0394 RNA stability',
         y='') +
    ## theme_gray(14) +
    theme_cowplot(14) +
    theme(legend.position='nil')
ggsave('figures/2021-12-01/RS-DESeq2_smg6VsWt_overlay.pdf', sp1,
       width=7, height=5)

### Plot the differences flox/flox treated vs untreated
rnaStab.ff1vsFf0$Type <- typef
data.ff1vsFf0 <- as.data.frame(rnaStab.ff1vsFf0)

sp1 <- ggplot(data.ff1vsFf0, aes(x=log2FoldChange, fill=Type)) +
    geom_density(alpha=0.5) +
    geom_vline(xintercept=0, linetype='dashed') +
    xlim(-1, 1.5) +
    labs(title='Smg6 affects RNA stability in cells',
         subtitle='Genes stratified by transcipt types',
         caption='smg6: flox/flox treated; Control: flox/flox untreated',
         fill='',
         x='\u0394 RNA stability') +
    theme_cowplot(14) +
    theme(legend.position=c(0.02, 0.9))
ggsave('figures/2021-12-01/RS-DESeq2_smg6VsUntreated.pdf', sp1,
       width=9, height=6)

sp1 <- ggplot(data.ff1vsFf0, aes(x=log2FoldChange, col=Type)) +
    stat_ecdf() +
    ## geom_vline(xintercept=0, linetype='dashed') +
    ## xlim(-1, 1.5) +
    labs(title='Smg6 affects RNA stability in cells',
         subtitle='Genes stratified by transcipt types',
         caption='smg6: flox/flox treated; Control: flox/flox untreated',
         col='',
         x='\u0394 RNA stability',
         y='Cumulative Sum') +
    theme_cowplot(14) +
    theme(legend.position=c(0.02, 0.9))
ggsave('figures/2021-12-01/RS-DESeq2_smg6VsUntreated_ecdf.pdf', sp1,
       width=6, height=6)


sp1 <- ggplot(data.ff1vsFf0, aes(x=log2FoldChange, y=Type, fill=Type)) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_density_ridges(alpha=0.5, scale=1.5) +
    xlim(-1, 1.5) +
    labs(title='Smg6 affects RNA stability in cells',
         subtitle='Genes stratified by transcipt types',
         caption='smg6: flox/flox treated; Control: flox/flox untreated',
         x='\u0394 RNA stability',
         y='') +
    ## theme_gray(14) +
    theme_cowplot(14) +
    theme(legend.position='nil')
ggsave('figures/2021-12-01/RS-DESeq2_smg6VsUntreated_overlay.pdf', sp1,
       width=7, height=5)

## write.table(TukeyHSD(aov(log2FoldChange ~ Type, data = data.ff1vsFf0))$Type,
##             file = "results/2021-12-01/RNAstabilityTypesTukey.csv")
write.csv(TukeyHSD(aov(log2FoldChange ~ Type, data = data.ff1vsFf0))$Type,
            file = "results/2021-12-01/RNAstabilityTypesTukey.csv")


######################################################################
### Plot the coverage of some of them (for the moment done with data
### from liver time course, not ideal):
######################################################################
assembly <- 'GRCm38.p5'
genesModels <- import('../annotations/Mmusculus.GRCm38.91.sorted.gtf',
                      format='gtf', genome=assembly)
transcripts <- which(genesModels$type == 'transcript')
exons <- which(genesModels$type == 'exon')
cds <- which(genesModels$type == 'CDS')
table(genesModels$type)
table(genesModels$transcript_biotype)

genesModels.df <- data.frame(chromosome = paste('chr',
                                                seqnames(genesModels)[exons], sep=""),
                             start = start(genesModels)[exons],
                             end = end(genesModels)[exons],
                             strand = strand(genesModels)[exons],
                             gene = genesModels$gene_id[exons],
                             transcript = genesModels$transcript_id[exons],
                             symbol = genesModels$gene_name[exons],
                             type = genesModels$transcript_biotype[exons])
head(genesModels.df)

######################################################################
### RNA-seq data
## Default tracks:
assembly <- 'mm10'

tgenes <- c('ENSMUSG00000024462', 'ENSMUSG00000066036',
            'ENSMUSG00000029104')
gid <- which(genesModels.df$gene %in% tgenes)
genesModels.plt <- genesModels.df[gid,]

genomeCoverage(id=tgenes[1], geneModels=genesModels.plt,
               assembly=assembly, ylim=c(0,40),
               dir='figures/2021-12-01')

genomeCoverage(id=tgenes[2], geneModels=genesModels.plt,
               assembly=assembly, ylim=c(0,1000),
               dir='figures/2021-12-01')

genomeCoverage(id=tgenes[3], geneModels=genesModels.plt,
               assembly=assembly, ylim=c(0,200),
               dir='figures/2021-12-01')


######################################################################
######################################################################
### Check gene expression changes in mRNA and perfor GO analysis
######################################################################
######################################################################

ddsELevel <- DESeqDataSetFromHTSeqCount(sampleTable = sampleExonAnno,
                                           directory = directory,
                                           design = ~ -1 + Genotype : Treatment)
## keep genes with at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(ddsELevel) >= 10) >= 3
ddsELevel <- ddsELevel[keep,]
ddsELevel <- DESeq(ddsELevel)
resultsNames(ddsELevel)

mRnaRes.ff1 <- results(ddsELevel,
                       contrast=list('Genotypeff.Treatment1',
                                     'Genotypeww.Treatment1'),
                       lfcThreshold=0, alpha=0.05)

summary(mRnaRes.ff1)

mRnaRes.ff1vsFf0 <- results(ddsELevel,
                            contrast=list('Genotypeff.Treatment1',
                                          'Genotypeff.Treatment0'),
                            lfcThreshold=0, alpha=0.05)

summary(mRnaRes.ff1vsFf0)

mRnaRes.ww1vsWw0 <- results(ddsELevel,
                            contrast=list('Genotypeww.Treatment1',
                                          'Genotypeww.Treatment0'),
                            lfcThreshold=0, alpha=0.05)
summary(mRnaRes.ww1vsWw0)



######################################################################
### Pathway analysis
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gID <- rownames(mRnaRes.ff1[which(mRnaRes.ff1$padj <= 0.05 &
                                  mRnaRes.ff1$log2FoldChange > 0), ])
genesID <- getBM(attributes=c("entrezgene_id",'kegg_enzyme'),
                 mart = mart, filters='ensembl_gene_id', values=gID)

gdID <- rownames(mRnaRes.ff1[which(mRnaRes.ff1$padj <= 0.05 &
                                  mRnaRes.ff1$log2FoldChange < 0), ])
gDownID <- getBM(attributes=c("entrezgene_id",'kegg_enzyme'),
                 mart = mart, filters='ensembl_gene_id', values=gdID)


pathAnalysis <- enrichPathway(gene=genesID$entrezgene_id,
                              organism='mouse', pvalueCutoff=0.05,
                              readable=T)
dim(pathAnalysis)
head(pathAnalysis)

pdf("figures/2021-12-01/Path_analysis_genesUp.pdf", width=10, height=20)
dotplot(pathAnalysis, showCategory=62, label_format=80) +
    ggtitle('Smg6 vs Control - Pathway')
dev.off()


### kegg
keggAnalysis <- enrichKEGG(gene         = genesID$entrezgene_id,
                           organism     = 'mmu',
                           pvalueCutoff = 0.05)

dim(keggAnalysis)
head(keggAnalysis)

pdf("figures/2021-12-01/kegg_analysis_genesUp.pdf", width=6, height=4)
dotplot(keggAnalysis, showCategory=3) +
    ggtitle('Smg6 vs Control - KEGG')
dev.off()

### GO MF
goMfAnalysis <- enrichGO(gene  = genesID$entrezgene_id,
                         OrgDb = org.Mm.eg.db,
                         ont = 'MF')
dim(goMfAnalysis)
head(goMfAnalysis)

pdf("figures/2021-12-01/GO-MF_genesUp.pdf", width=8, height=10)
dotplot(goMfAnalysis, showCategory=25, label_format=50) +
    ggtitle('Smg6 vs Control - GO MF')
dev.off()

### GO BP
goBpAnalysis <- enrichGO(gene  = genesID$entrezgene_id,
                         OrgDb = org.Mm.eg.db,
                         ont = 'BP')
dim(goBpAnalysis)
head(goBpAnalysis)

intCat <- c(grep('rhythm', goBpAnalysis$Description),
            grep('cycle', goBpAnalysis$Description)
            )


pdf("figures/2021-12-01/GO-BP_genesUp.pdf", width=8, height=10)
dotplot(goBpAnalysis, showCategory=goBpAnalysis$Description[intCat],
        label_format=50) +
    ggtitle('Smg6 vs Control - GO BP')
dev.off()

### GO BPDown
goBpDAnalysis <- enrichGO(gene  = gDownID$entrezgene_id,
                          OrgDb = org.Mm.eg.db,
                          ont = 'BP')
dim(goBpDAnalysis)
head(goBpDAnalysis)

n1 <- grep('nonsense', goBpDAnalysis$Description)
n2 <- grep('RNA', goBpDAnalysis$Description)


pdf("figures/2021-12-01/GO-BP_genesDown.pdf", width=8, height=10)
dotplot(goBpDAnalysis, showCategory=goBpDAnalysis$Description[c(n1,n2)],
        label_format=50) +
    ggtitle('Smg6 vs Control - GO BP')
dev.off()

### GO CC
goCcAnalysis <- enrichGO(gene  = genesID$entrezgene_id,
                         OrgDb = org.Mm.eg.db,
                         ont = 'CC')
dim(goCcAnalysis)
head(goCcAnalysis)

pdf("figures/2021-12-01/GO-CC_genesUp.pdf", width=8, height=10)
dotplot(goCcAnalysis, showCategory=25, label_format=50) +
    ggtitle('Smg6 vs Control - GO CC')
dev.off()


######################################################################
### Chage the model to get genes that are treatment specific, genotype
### specific and treatmetn and genotype specific (these results are
### heavily skewed due to the unbalances in exprimental design):

ddsELevel2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleExonClean,
                                         directory = directory,
                                         design = ~ Genotype * Treatment)
## keep genes with at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(ddsELevel2) >= 10) >= 3
ddsELevel2 <- ddsELevel2[keep,]
ddsELevel2 <- DESeq(ddsELevel2)
resultsNames(ddsELevel2)

res.geno <- results(ddsELevel2, name='Genotype_ff_vs_ww', alpha=0.05)
summary(res.geno)

res.treat <- results(ddsELevel2, name='Treatment_1_vs_0', alpha=0.05)
summary(res.treat)

res.smg6 <- results(ddsELevel2, name='Genotypeff.Treatment1', alpha=0.05)
summary(res.smg6)


######################################################################
### Plot the differentially expressed genes

allDiff <- which(mRnaRes.ff1$padj < 0.05 | mRnaRes.ff1vsFf0$padj < 0.05)
length(allDiff)
sampleExonClean <- sampleExonAnno[which(sampleExonAnno$Genotype != 'fw'),]

exprsDiff <- counts(ddsELevel, norm=T)[allDiff, sampleExonClean[,1]]
dim(exprsDiff)
sampleOrd <- order(sampleExonClean$Genotype, sampleExonClean$Treatment)

exprsDiffZ <- zScore(exprsDiff[,sampleOrd])

head(exprsDiffZ)

heatAnno <- sampleExonClean[,3:4]
rownames(heatAnno) <- sampleExonClean[,1]

jpeg('figures/2021-12-01/differentiallyExpressedHeatmap.jpg', height=1000)
pheatmap(exprsDiffZ,
         main = "",
         cluster_cols = TRUE,
         annotation_col = heatAnno,
         show_rownames = FALSE,
         cutree_rows = 4,
         cutree_cols = 2)
dev.off()


######################################################################
### Save image
## save.image('2021-12-01_pooledCellsRNAstab.RData')
load('2021-12-01_pooledCellsRNAstab.RData')
