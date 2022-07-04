######################################################################
### Script to run some initial analyses on pre-mRNA and mRNA 

## library(packrat)
## packrat::init("./", infer.dependencies = FALSE)
library(rtracklayer) # to import GTF file
library(tximport)
library(DESeq2)
library(MetaCycle)
library(org.Mm.eg.db)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
## library(ggfortify)
## library(plyr)
library(ggrepel)
library(ggpubr)
library(pheatmap)

######################################################################
### Some functions
source('smg6functions.R')

######################################################################
### Get the annotation

annot <- import("/mnt/data/databases/mouse/gtf/Mmusculus.GRCm38.91.gtf")
genesAnno <- annot[annot$type == 'gene']
genesSel <- genesAnno[genesAnno$gene_biotype %in% c("protein_coding")]

directory <- '../genesCounts'

coreClock <- read.table('core_clock_genes.csv', header=F)

## New sequencing Run
sampleExonAnno <- read.csv('../2021-10-29_smg6Samples.anno', sep='\t')
sampleExonAnno$Genotype <- factor(sampleExonAnno$Genotype,
                                   levels=c('MUT','WT'))
sampleExonAnno$Ntime <- sampleExonAnno$ZT
sampleExonAnno$ZT <- factor(sampleExonAnno$ZT,
                             levels=c(0, 4, 8, 12, 16, 20))
sampleExonAnno$LibraryID <- as.character(sampleExonAnno$LibraryID)
sampleExonAnno$FileName <- paste(sampleExonAnno$LibraryID,
                                  'exon.htcount.dat', sep='.')
rownames(sampleExonAnno) <- sampleExonAnno$LibraryID
sampleExonAnno$Type <- factor(paste(sampleExonAnno$Genotype,
                                     sampleExonAnno$ZT, sep='.'))
sampleExonAnno$Run <- factor(1, levels=c(1,2))


summary(sampleExonAnno)

######################################################################
### Read exon reads (my script)
######################################################################
ddsExonData <- DESeqDataSetFromHTSeqCount(sampleTable = sampleExonAnno,
                                          directory = directory,
                                          design = ~ -1 + Type)
## get the differentially expressed exons:
ddsExonData <- DESeq(ddsExonData)

resultsNames(ddsExonData)

res00.mRNA <- results(ddsExonData, contrast=c("Type", "MUT.0", "WT.0"))
summary(res00.mRNA)
res04.mRNA <- results(ddsExonData, contrast=c("Type", "MUT.4", "WT.4"))
summary(res04.mRNA)
res08.mRNA <- results(ddsExonData, contrast=c("Type", "MUT.8", "WT.8"))
summary(res08.mRNA)
res12.mRNA <- results(ddsExonData, contrast=c("Type", "MUT.12", "WT.12"))
summary(res12.mRNA)
res16.mRNA <- results(ddsExonData, contrast=c("Type", "MUT.16", "WT.16"))
summary(res16.mRNA)
res20.mRNA <- results(ddsExonData, contrast=c("Type", "MUT.20", "WT.20"))
summary(res20.mRNA)

######################################################################
### Read exon reads (htseq-count)
######################################################################
sampleExonAnno2 <- sampleExonAnno
sampleExonAnno2$FileName <- paste(sampleExonAnno2$LibraryID,
                                 "exon.htseqReverse.dat", sep=".")
head(sampleExonAnno2)

ddsExonData2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleExonAnno2,
                                           directory = directory,
                                           design = ~ -1 + Type)
## get the differentially expressed exons:
ddsExonData2 <- DESeq(ddsExonData2)



######################################################################
### Read gene locus data
######################################################################
sampleGeneAnno <- sampleExonAnno
sampleGeneAnno$FileName <- paste(sampleGeneAnno$LibraryID,
                                 "gene.htseqReverse.dat", sep=".")
head(sampleGeneAnno)


ddsGeneData <- DESeqDataSetFromHTSeqCount(sampleTable = sampleGeneAnno,
                                          directory = directory,
                                          design = ~ -1 + Type)
## get the differentially expressed genes:
ddsGeneData <- DESeq(ddsGeneData)
geneSizeFactors <- sizeFactors(ddsGeneData)


######################################################################
### Calculate pre-mRNA reads (my script)
######################################################################
exonCounts <- counts(ddsExonData, norm=F)
geneCounts <- counts(ddsGeneData, norm=F)

dim(exonCounts)
dim(geneCounts)

intronCounts <- geneCounts[rownames(exonCounts),] - exonCounts

summary(intronCounts)
intronCounts[intronCounts < 0] <- 0

sum(intronCounts)
sum(exonCounts)

ddsPremRNA <- DESeqDataSetFromMatrix(intronCounts,
                                     colData=sampleGeneAnno,
                                     design = ~ -1 + Type)

## Normalize pre-mRNA reads with size factors estimated using the
## global library
sizeFactors(ddsPremRNA) <- geneSizeFactors

## Calculate differentially expressed pre-mRNA
ddsPremRNA <- DESeq(ddsPremRNA)

resultsNames(ddsPremRNA)
res00.pre <- results(ddsPremRNA, contrast=c("Type", "MUT.0", "WT.0"))
res04.pre <- results(ddsPremRNA, contrast=c("Type", "MUT.4", "WT.4"))
res08.pre <- results(ddsPremRNA, contrast=c("Type", "MUT.8", "WT.8"))
res12.pre <- results(ddsPremRNA, contrast=c("Type", "MUT.12", "WT.12"))
res16.pre <- results(ddsPremRNA, contrast=c("Type", "MUT.16", "WT.16"))
res20.pre <- results(ddsPremRNA, contrast=c("Type", "MUT.20", "WT.20"))

summary(res00.pre)
summary(res04.pre)
summary(res08.pre)
summary(res12.pre)
summary(res16.pre)
summary(res20.pre)


### In smg6, both pre-mRNA and and mRNA levels should go upfor genes
### that are NMD sensitive.
x <- rownames(res08.pre)[which(res08.pre$log2FoldChange > 0 & res08.pre$padj < 0.05)]
y <- rownames(res08.mRNA)[which(res08.mRNA$log2FoldChange > 0 & res08.mRNA$padj < 0.05)]

length(which(x %in% y))
length(y)

res08.pre[x[which(x %in% y)],]
res08.mRNA[x[which(x %in% y)],]


######################################################################
### Calculate pre-mRNA reads (htseq-count)
######################################################################
exonCounts2 <- counts(ddsExonData2, norm=F)

dim(exonCounts2)
dim(geneCounts)

intronCounts2 <- geneCounts - exonCounts2

summary(intronCounts2)
intronCounts2[intronCounts2 < 0] <- 0

sum(intronCounts2)
sum(exonCounts)

ddsPremRNA2 <- DESeqDataSetFromMatrix(intronCounts2,
                                     colData=sampleGeneAnno,
                                     design = ~ -1 + Type)

## Normalize pre-mRNA reads with size factors estimated using the
## global library
sizeFactors(ddsPremRNA2) <- geneSizeFactors

## Calculate differentially expressed pre-mRNA
ddsPremRNA2 <- DESeq(ddsPremRNA2)


######################################################################
### Plot Cry2 pre-mRNA and mRNA total reads
cry2 <- 'ENSMUSG00000068742'
cry2Data <- data.frame(rbind(sampleExonAnno, sampleExonAnno),
                       counts = log2(c(counts(ddsExonData, norm=T)[cry2,],
                                       counts(ddsPremRNA, norm=T)[cry2,])),
                       type = rep(c('mRNA', 'pre-mRNA'),
                                  each=dim(sampleExonAnno)[1])
                       )

p1 <- ggplot(cry2Data, aes(x=Genotype, y=counts, fill=Genotype)) +
    geom_boxplot() +
    geom_jitter(height=0, width=0.2) +
    facet_grid(~ type) +
    labs(title='Cry2 expression',
         y = 'Normalized counts') +
    theme_gray(14) +
    theme(legend.position='none')
ggsave('figures/2021-12-23/global_expression_cry2.pdf')


p1 <- ggplot(cry2Data, aes(x=ZT, y=counts, fill=Genotype,
                           shape=Genotype)) +
    geom_boxplot() +
    geom_jitter(height=0, width=0.05) +
    facet_grid(~ type) +
    labs(title='Cry2 expression',
         y = 'Normalized counts') +
    theme_gray(14) +
    theme(legend.position='none')
ggsave('figures/2021-12-23/global_expression_ZT_cry2.pdf', width=8)

cry2StabData <- data.frame(sampleExonAnno,
                           counts = log2(counts(ddsExonData, norm=T)[cry2,]/
                                         counts(ddsPremRNA, norm=T)[cry2,]),
                           type = rep('mRNA stability',
                                      dim(sampleExonAnno)[1])
                           )
cry2StabData$ZTn <- as.numeric(as.character(cry2StabData$ZT))
cry2StabData$Phase <- 'Light'
cry2StabData$Phase[cry2StabData$ZTn >= 12] <- 'Dark'
cry2StabData$Phase <- factor(cry2StabData$Phase, levels=c('Light', 'Dark'))

p1 <- ggplot(cry2StabData, aes(x=Genotype, y=counts, fill=Genotype)) +
    geom_boxplot() +
    geom_jitter(height=0, width=0.2) +
    labs(title='Cry2 RNA stability',
         y = 'log2(mRNA) - log2(pre-mRNA)') +
    theme_gray(14) +
    theme(legend.position='none')
ggsave('figures/2021-12-23/global_RNAstability_cry2.pdf')

p1 <- ggplot(cry2StabData, aes(x=Genotype, y=counts, fill=Genotype)) +
    geom_boxplot() +
    geom_jitter(height=0, width=0.2) +
    facet_wrap(~ Phase) +
    labs(title='Cry2 RNA stability',
         y = 'log2(mRNA) - log2(pre-mRNA)') +
    theme_gray(14) +
    theme(legend.position='none')
ggsave('figures/2021-12-23/global_RNAstability_cry2_phase.pdf')

summary(lm(counts ~ Genotype, data = subset(cry2StabData, Phase = "Dark")))


p1 <- ggplot(cry2StabData, aes(x=ZT, y=counts, fill=Genotype,
                               shape=Genotype)) +
    geom_boxplot() +
    geom_point(size=3, alpha=0.5, position=position_jitterdodge()) +
    labs(title='Cry2 RNA stability',
         y = 'log2(mRNA) - log2(pre-mRNA)') +
    theme_gray(14)
ggsave('figures/2021-12-23/global_RNAstability_ZT_cry2.pdf')


p1 <- ggplot(cry2StabData,
             aes(x=ZTn, y=counts, col=Genotype, shape=Genotype)) +
    geom_point(size=3, alpha=0.5) +
    geom_smooth(se=FALSE) +
    labs(title='Cry2 RNA stability',
         x = 'ZT',
         y = 'log2(mRNA) - log2(pre-mRNA)') +
    theme_gray(14)
ggsave('figures/2021-12-23/global_RNAstability_ZT_cry2_smooth.pdf',
       p1)



######################################################################
### Repeat the analysis for all core clock genes
cc <- coreClock[,2]
cc <- cc[which(cc %in% rownames(ddsExonData))]

ccStabData <- data.frame()
for(N in cc){
    tmpData <- data.frame(sampleExonAnno,
                          counts = log2(counts(ddsExonData, norm=T)[N,])-
                                   log2(counts(ddsPremRNA, norm=T)[N,]),
                          type = rep('mRNA stability',
                                     dim(sampleExonAnno)[1]),
                          gene = rep(N, dim(sampleExonAnno)[1])
                          )
    ccStabData <- rbind(ccStabData, tmpData)
}


head(ccStabData)
ccStabData$ZTn <- as.numeric(as.character(ccStabData$ZT))
ccStabData$Phase <- 'Light'
ccStabData$Phase[ccStabData$ZTn >= 12] <- 'Dark'
ccStabData$Phase <- factor(ccStabData$Phase, levels=c('Light', 'Dark'))
ccStabData <- ccStabData[-which(ccStabData$counts == Inf),]

p1 <- list()
for(N in unique(ccStabData$gene)){
    ngene <- coreClock[which(coreClock[,2] == N), 1]
    f1name <- paste('figures/2021-12-23/global_RNAstability_phase_',
                    ngene, '.pdf', sep='')
    p1[[N]] <- ggplot(subset(ccStabData, gene==N),
                 aes(x=Genotype, y=counts, fill=Genotype)) +
        geom_boxplot() +
        geom_jitter(height=0, width=0.2) +
        facet_wrap(~ Phase) +
        labs(title=paste(ngene, 'RNA stability'),
             y = 'log2(mRNA) - log2(pre-mRNA)') +
        theme_gray(14) +
        theme(legend.position='none')
    ##ggsave(f1name, p1)
}

pdf('figures/2021-12-23/global_RNAstability_phase_coreClock.pdf',
    width=20, height=25)
do.call(grid.arrange, p1)
dev.off()

p2 <- list()
for(N in unique(ccStabData$gene)){
    ngene <- coreClock[which(coreClock[,2] == N), 1]
    f2name <- paste('figures/2021-12-23/global_RNAstability_ZT_',
                    ngene, '.pdf', sep='')
    p2[[N]] <- ggplot(subset(ccStabData, gene==N),
                      aes(x=ZT, y=counts, fill=Genotype,
                          shape=Genotype)) +
        geom_boxplot() +
        geom_point(size=3, alpha=0.5, position=position_jitterdodge()) +
        labs(title=paste(ngene, 'RNA stability'),
             y = 'log2(mRNA) - log2(pre-mRNA)') +
        theme_gray(14)
    ## ggsave(f2name, p2)
}

pdf('figures/2021-12-23/global_RNAstability_ZT_coreClock.pdf',
    width=30, height=25)
do.call(grid.arrange, p2)
dev.off()




######################################################################
### Again but with exon reads estimated by htseq-count
cry2Data <- data.frame(rbind(sampleExonAnno, sampleExonAnno),
                       counts = log2(c(counts(ddsExonData2, norm=T)[cry2,],
                                       counts(ddsPremRNA2, norm=T)[cry2,])),
                       type = rep(c('mRNA', 'pre-mRNA'),
                                  each=dim(sampleExonAnno)[1])
                       )

p1 <- ggplot(cry2Data, aes(x=Genotype, y=counts, fill=Genotype)) +
    geom_boxplot() +
    geom_jitter(height=0, width=0.2) +
    facet_grid(~ type) +
    labs(title='Cry2 expression',
         y = 'Normalized counts') +
    theme_gray(14) +
    theme(legend.position='none')
ggsave('figures/2021-12-23/global_expression_cry2_htseqCount.pdf')

cry2StabData <- data.frame(sampleExonAnno,
                           counts = log2(counts(ddsExonData2, norm=T)[cry2,]/
                                         counts(ddsPremRNA2, norm=T)[cry2,]),
                           type = rep('mRNA stability',
                                      dim(sampleExonAnno)[1])
                           )

p1 <- ggplot(cry2StabData, aes(x=Genotype, y=counts, fill=Genotype)) +
    geom_boxplot() +
    geom_jitter(height=0, width=0.2) +
    labs(title='Cry2 RNA stability',
         y = 'log2(mRNA) - log2(pre-mRNA)') +
    theme_gray(14) +
    theme(legend.position='none')
ggsave('figures/2021-12-23/global_RNAstability_cry2_htseqCount.pdf')

p1 <- ggplot(cry2StabData, aes(x=ZT, y=counts, fill=Genotype,
                               shape=Genotype)) +
    geom_boxplot() +
    geom_point(size=3, alpha=0.5, position=position_jitterdodge()) +
    labs(title='Cry2 RNA stability',
         y = 'log2(mRNA) - log2(pre-mRNA)') +
    theme_gray(14)
ggsave('figures/2021-12-23/global_RNAstability_ZT_cry2_htseqCount.pdf')


######################################################################
### Plot all genes (as done on 2021-11-09)

geneNmaes <- genesSel$gene_name
names(geneNmaes) <- genesSel$gene_id
head(geneNmaes)

for(I in rownames(ddsExonData)){
    tname <- paste(I, geneNmaes[I], sep='_')
    pf1 <- plotGene(I, ddsExonData, title=tname, sub='mRNA')
    pf2 <- plotGene(I, ddsPremRNA, title=tname, sub='pre-mRNA', legend=T)
    fname <- paste('smg6Expression/', tname, '_mRNA-premRNA.pdf', sep="")
    ## pdf(fname, width=10, height=6)
    ggarrange(pf1, pf2, ncol = 2, nrow = 1) %>%
        ggexport(filename = fname, height=6, width=11)
    ## dev.off()
}


######################################################################
### make some test to see how RNAstability is affected by oscillatory
### expression

geneCounts <- plotCounts(ddsExonData,
                         gene='ENSMUSG00000029238',
                         intgroup=c("Genotype", "ZT"),
                         returnData = TRUE)

pmrnaCounts <- plotCounts(ddsPremRNA,
                          gene='ENSMUSG00000029238',
                          intgroup=c("Genotype", "ZT"),
                          returnData = TRUE)

counts(ddsPremRNA, norm=T)['ENSMUSG00000029238',]
counts(ddsExonData, norm=T)['ENSMUSG00000029238',]


geneCounts$RNAstab <- log2(geneCounts$count) - log2(pmrnaCounts$count)


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

summary(lm(y ~ cos(x*(2*pi/24)) + sin(x*(2*pi/24)), data=fitData))


intr <- 10
b1 <- 0.47687
b2 <- -1.12618
x <- seq(0,20,4)

cPred <- intr + b1*cos(x*(2*pi/24)) + b2*sin(x*(2*pi/24))

### pre-mRNA has the same phase as mRNA:
pPred <- (intr-2) + b1*cos(x*(2*pi/24)) + b2*sin(x*(2*pi/24))

cPred-pPred

## the RNA-stability is constant

### pre-mRNA is 2 h in advance:
pPred <- ((intr-2)*1.1) + (-b1*0.5)*cos(x*(2*pi/24)) + (b2*0.5)*sin(x*(2*pi/24))

cPred-pPred

## the RNA-stability is cyclic with a phase that coincide with the
## mRNA phase

pPred <- ((intr-2)*1.1) + (-b1*0.5)*cos(x*(2*pi/24)) + (-b2*0.5)*sin(x*(2*pi/24))

cPred-pPred


######################################################################
### Save image
## save.image('2021-12-23_liverRnaStability.RData')
load('2021-12-23_liverRnaStability.RData')
