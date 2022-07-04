######################################################################
### Script to run some initial analyses on pre-mRNA and mRNA 

library(rtracklayer) # to import GTF file
library(sva)
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
library(glmpca)


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

sampleExonAnno <- read.table('../samples.anno', header=T, sep="\t")
sampleExonAnno$Genotype <- factor(sampleExonAnno$Genotype,
                                  levels=c('+/+', 'f/+', 'f/f'))
levels(sampleExonAnno$Genotype) <- c('ww', 'fw', 'ff')
sampleExonAnno$Treatment <- as.factor(sampleExonAnno$Treatment)
sampleExonAnno$ID <- as.character(sampleExonAnno$ID)
sampleExonAnno$FileName <- paste(sampleExonAnno$ID, 'exon.coverage.dat', sep='.')
str(sampleExonAnno)

summary(sampleExonAnno)
head(sampleExonAnno)

sampleExonClean <- sampleExonAnno[which(sampleExonAnno$Genotype != 'fw'),]


######################################################################
### Read exon reads
######################################################################
ddsExonData <- DESeqDataSetFromHTSeqCount(sampleTable = sampleExonClean,
                                          directory = directory,
                                          design = ~ Genotype * Treatment)

## get the differentially expressed exons (keep exons with at least 3
## samples with a count of 10 or higher):
keep <- rowSums(counts(ddsExonData) >= 20) >= 6
ddsExonData <- ddsExonData[keep,]
ddsExonData <- DESeq(ddsExonData)

## Account for hidden variables (cell clone?)
dat  <- counts(ddsExonData, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ -1 + Genotype * Treatment, colData(ddsExonData))
mod0 <- model.matrix(~   1, colData(ddsExonData))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv

ddsExonDataHV <- ddsExonData
ddsExonDataHV$SV1 <- svseq$sv[,1]
## ddsExonDataHV$SV2 <- svseq$sv[,2]
design(ddsExonDataHV) <- ~ Genotype * Treatment + SV1
ddsExonDataHV <- DESeq(ddsExonDataHV)

resultsNames(ddsExonDataHV)

res.exon <- lfcShrink(ddsExonDataHV,
                      coef="Genotypeff.Treatment1", alpha=0.1,
                      lfcThreshold=log2(2), type="ashr")
summary(res.exon)

######################################################################
### Annotate the exons and plot them in several ways
######################################################################
annoData <- mcols(annot)
colnames(annoData)

exonToGene <- read.table(paste(directory, sampleExonAnno$FileName,
                               sep="/")[1], as.is=T)[,1:2]
rownames(exonToGene) <- exonToGene[,1]
head(exonToGene)

geneToExon <- list()
for(I in unique(exonToGene[,2])){
    geneToExon[[I]] <- exonToGene[which(exonToGene[,2] == I),1]
}

length(geneToExon)


sg.up <- rownames(res.exon[which(res.exon$padj <= 0.05 & res.exon$log2FoldChange > log2(2)),])
sg.dw <- rownames(res.exon[which(res.exon$padj <= 0.05 & res.exon$log2FoldChange < -log2(2)),])

signGenes <- unique(c(sg.up, sg.dw))
length(signGenes)


type <- vector(length=length(signGenes))
gName <- vector(length=length(signGenes))
names(type) <- names(gName) <- signGenes
for (I in signGenes){
    eIndex <- which(annoData$exon_id == I)
    gName[I] <- annoData$gene_name[eIndex[1]]
    for (K in annoData$transcript_biotype[eIndex]){
        if (K == 'protein_coding' ){
            type[I] <- K
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

type <- as.factor(type)
summary(type)

nmd.exons <- names(type[which(type == 'nonsense_mediated_decay')])
pc.exons <- names(type[which(type == 'protein_coding')])
ri.exons <- names(type[which(type == 'retained_intron')])
ot.exons <- names(type[which(type == 'Other')])

######################################################################
### Try some other visualsations to see if beter (at 0h)
tab <- data.frame(logFC = res.exon$log2FoldChange,
                  negLogPval = -log10(res.exon$padj),
                  AverageExpression = log2(res.exon$baseMean))
rownames(tab) <- res.exon@rownames
sg.df <- c(sg.up, sg.dw)
tabSig <- data.frame(tab[sg.df,], Type='Constitutive')
tabSig[sg.df[sg.df %in% nmd.exons], 'Type'] <- 'NMD'
tabSig[sg.df[sg.df %in% ri.exons], 'Type'] <- 'Retained Intron'
tabMA <- data.frame(tab, Type=rep('NC', dim(tab)[1]))
tabMA[sg.df, 'Type'] <- 'Constitutive'
tabMA[sg.df[sg.df %in% nmd.exons], 'Type'] <- 'NMD'
tabMA[sg.df[sg.df %in% ri.exons], 'Type'] <- 'Retained Intron'
pD1 <- ggplot(tabMA[which(tabMA$Type=='NC'),],
              aes(y=logFC, x=AverageExpression)) +
    geom_density_2d_filled(data=tabMA[which(tabMA$Type=='NC'),],
                           aes(y=logFC, x=AverageExpression),
                           breaks=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                           col='gray', fill='lightgray') +
    geom_density_2d(data=tabMA[sg.df,],
                    aes(y=logFC, x=AverageExpression, col=Type),
                    breaks=c(0.025, 0.5), size=2) +
    geom_point(data=tabMA[sg.df,],
               aes(y=logFC, x=AverageExpression, col=Type), alpha=0.7) +
    scale_color_manual(values=c("#999999", "darkred", "darkgreen")) +
    labs(title='Differentially expressed exon', colour='') +
    ylim(-6,6) +
    theme_classic() +
    theme(axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold"),
          axis.title = element_text(size=14, face="bold"),
          axis.line = element_line(colour='black', size = 1),
          axis.ticks = element_line(colour='black', size = 1),
          axis.ticks.length = unit(.25, "cm"),
          legend.position = c(0.15, 0.9),
          ## legend.key.size = unit(1, "cm"),
          legend.text = element_text(size=14)) +
    guides(colour = guide_legend(override.aes = list(size=2)))
ggsave('figures/2021-12-09/MAplot_KOvsCtrl_density.pdf', pD1) 



######################################################################
### Save image
## save.image('2021-12-09_pooledCellsExonAnalysis.RData')
load('2021-12-09_pooledCellsExonAnalysis.RData')
