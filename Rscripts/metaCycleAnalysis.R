######################################################################
### This script tests several plots to better represent mRNA and
### pre-mRNA expression

library(DESeq2)
library(MetaCycle)
library(org.Mm.eg.db)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

######################################################################
### Load the plotting function
source('smg6plotFunctions.R')
source('smg6functions.R')

######################################################################
### Read in the data
load('2021-12-23_liverRnaStability.RData')

######################################################################
### calculate TPM and make phase and amplitude analysis
trsLengths <- read.table('../annotations/2021-12-07_exonAndIntronLengths.txt',
                         header=T, row.names=1)
str(trsLengths)
trsLengths <- trsLengths/1000

noIntron <- rownames(trsLengths)[which(trsLengths[,2] > 0) ]

cry2 <- 'ENSMUSG00000068742'

mRNAC <- counts(ddsExonData, norm=F)
premRNAC <- counts(ddsPremRNA, norm=F)

mRNAanno <- data.frame(sampleExonAnno,
                       Counts = colSums(mRNAC),
                       Ctype = "mRNA")
premRNAanno <- data.frame(sampleExonAnno,
                          Counts = colSums(premRNAC),
                          Ctype = "pre-mRNA")

countAnno <- rbind(mRNAanno, premRNAanno)
countAnno$Type <- factor(countAnno$Type, levels = c("MUT.0", "MUT.4",
                                                    "MUT.8", "MUT.12",
                                                    "MUT.16",
                                                    "MUT.20", "WT.0",
                                                    "WT.4", "WT.8",
                                                    "WT.12", "WT.16",
                                                    "WT.20"))

countAnno <- countAnno %>%
    mutate(TR = factor(interaction(Replica, Type)))
head(countAnno)

p1 <- ggplot(countAnno, aes(x = TR, y = Counts, fill = Ctype)) +
    geom_bar(stat = "identity", position="fill") +
    coord_flip() +
    labs(x = "",
         y = "",
         fill = "") +
    theme(legend.position = "top")
ggsave("figures/2022-02-08/libraryComposition.pdf", p1, width = 6,
       height = 6)


sFactor <- tpm(mRNAC, trsLengths[rownames(mRNAC),1], sFactor=TRUE)
mRNAtpm <- tpm(mRNAC, trsLengths[rownames(mRNAC),1],
               library=sFactor)
premRNAtpm <- tpm(premRNAC, trsLengths[rownames(premRNAC),2],
                  library=sFactor)

sFactor2 <- colSums(counts(ddsExonData, norm=F))
mRNArpkm <- rpkm(mRNAC, trsLengths[rownames(mRNAC),1],
               library=sFactor2)
premRNArpkm <- rpkm(premRNAC, trsLengths[rownames(premRNAC),2],
                  library=sFactor2)

## mRNArpkm <- rpkm(mRNAC, trsLengths[rownames(mRNAC),1])
## premRNArpkm <- rpkm(premRNAC, trsLengths[rownames(premRNAC),2])

ind <- unique(c(which(is.na(rowSums(mRNArpkm))),
                which(is.na(rowSums(premRNArpkm)))))
mRNArpkm <- mRNArpkm[-ind,]
premRNArpkm <- premRNArpkm[-ind,]
clean <- rownames(mRNArpkm)[which(rownames(mRNArpkm) %in% noIntron)]
mRNArpkm <- mRNArpkm[clean,]
premRNArpkm <- premRNArpkm[clean,]


mRNArpkm[cry2,]
premRNAtpm[cry2,]

summary(colSums(premRNArpkm[, samplesMut]))
summary(colSums(premRNArpkm[, samplesWt]))
summary(sFactor2[samplesMut]/1e6)
summary(sFactor2[samplesWt]/1e6)

### It seems like rpkm is not normalising well the pre-mRNA reads
### artificially enlarging too much low depth samples

pdf('figures/2022-02-08/libraryDepthVsRpkm.pdf')
plot(sFactor2 ~ colSums(premRNArpkm), xlab = "pre-mRNA RPKM sums",
     ylab = "Gene Count Sums")
dev.off()

pdf('figures/2022-02-08/libraryDepthVsRpkm_mRNA.pdf')
plot(sFactor2 ~ colSums(mRNArpkm), xlab = "mRNA RPKM sums",
     ylab = "Gene Count Sums")
dev.off()


######################################################################
### make the phase and amplitude analysis:

### mRNA:
## Mutant
samplesMut <- sampleExonAnno$LibraryID[which(sampleExonAnno$Genotype == "MUT")]
timeMut <- as.integer(as.character(sampleExonAnno$ZT[which(sampleExonAnno$Genotype == "MUT")]))
mutantmRNArpkm <- log(mRNArpkm[, samplesMut])
mutantmRNAtpm <- log(mRNAtpm[, samplesMut])

write.table(mutantmRNArpkm,
            file='results/2022-02-08/mutantmRNArpkm.txt', col.names=F,
            row.names=T, sep="\t")

## WT:
samplesWt <- sampleExonAnno$LibraryID[which(sampleExonAnno$Genotype == "WT")]
timeWt <- as.integer(as.character(sampleExonAnno$ZT[which(sampleExonAnno$Genotype == "WT")]))
wtmRNArpkm <- log(mRNArpkm[, samplesWt])
wtmRNAtpm <- mRNAtpm[, samplesWt]

write.table(wtmRNArpkm,
            file='results/2022-02-08/wtmRNArpkm.txt', col.names=F,
            row.names=T, sep="\t")

meta2d(infile="results/2022-02-08/mutantmRNArpkm.txt",
       filestyle="txt", timepoints=timeMut,
       cycMethod=c("ARS", "JTK", "LS"), outdir="results/2022-02-08/",
       adjustPhase = "predictedPer", outIntegration="both",
       outRawData=TRUE, minper=24, maxper=24, parallelize=TRUE, nCores=20)

meta2d(infile="results/2022-02-08/wtmRNArpkm.txt",
       filestyle="txt", timepoints=timeMut,
       cycMethod=c("ARS", "JTK", "LS"), outdir="results/2022-02-08/",
       adjustPhase = "predictedPer", outIntegration="both",
       outRawData=TRUE, minper=24, maxper=24, parallelize=TRUE, nCores=20)


### premRNA:
## Mutant
mutantpremRNArpkm <- log(premRNArpkm[, samplesMut])
mutantpremRNAtpm <- premRNAtpm[, samplesMut]

write.table(mutantpremRNArpkm,
            file='results/2022-02-08/mutantPremRNArpkm.txt', col.names=F,
            row.names=T, sep="\t")

## WT:
wtpremRNArpkm <- log(premRNArpkm[, samplesWt])
wtpremRNAtpm <- premRNAtpm[, samplesWt]
write.table(wtpremRNArpkm,
            file='results/2022-02-08/wtPremRNArpkm.txt', col.names=F,
            row.names=T, sep="\t")

meta2d(infile="results/2022-02-08/mutantPremRNArpkm.txt",
       filestyle="txt", timepoints=timeMut,
       cycMethod=c("ARS", "JTK", "LS"), outdir="results/2022-02-08/",
       adjustPhase = "predictedPer", outIntegration="both",
       outRawData=TRUE, minper=24, maxper=24, parallelize=TRUE, nCores=20)

meta2d(infile="results/2022-02-08/wtPremRNArpkm.txt",
       filestyle="txt", timepoints=timeMut,
       cycMethod=c("ARS", "JTK", "LS"), outdir="results/2022-02-08/",
       adjustPhase = "predictedPer", outIntegration="both",
       outRawData=TRUE, minper=24, maxper=24, parallelize=TRUE, nCores=20)

######################################################################
### Read in the results:
phaseMrnaMut <- read.table('results/2022-02-08/meta2d_mutantmRNArpkm.txt',
                       header=T, row.names=1)
phaseMrnaWt <- read.table('results/2022-02-08/meta2d_wtmRNArpkm.txt',
                       header=T, row.names=1)

phasePreMrnaMut <- read.table('results/2022-02-08/meta2d_mutantPremRNArpkm.txt',
                       header=T, row.names=1)
phasePreMrnaWt <- read.table('results/2022-02-08/meta2d_wtPremRNArpkm.txt',
                       header=T, row.names=1)

colnames(phaseMrnaMut)

######################################################################
### Plot core clock genes:
coreClock <- read.table(file="core_clock_genes.csv", sep="\t",
                        header=FALSE)
colnames(coreClock) <- c("Name", "ID")
targetClock <- c('Clock', 'Bmal1', 'Per1', 'Per2', 'Cry1', 'Cry2',
                 'Nr1d1', 'Nr1d2', 'Rorc')


phaseMrnaMut[coreClock$ID, 1:4]


phaseData <- data.frame(rbind(coreClock, coreClock, coreClock,
                              coreClock),
                        rbind(phaseMrnaMut[coreClock$ID, c(11,14)],
                              phaseMrnaWt[coreClock$ID, c(11,14)],
                              phasePreMrnaMut[coreClock$ID, c(11,14)],
                              phasePreMrnaWt[coreClock$ID, c(11,14)]),
                        type = rep(c("mRNA", "premRNA"),
                                   each=2*dim(coreClock)[1]),
                        genotype = rep(c("smg6", "WT", "smg6", "WT"),
                                       each=dim(coreClock)[1])
                        )
head(phaseData)

phaseData$Name %in% targetClock

p1 <- ggplot(subset(phaseData, type == "mRNA" & meta2d_pvalue < 0.05),
             aes(x=meta2d_phase, y=1, col=Name)) +
    geom_segment(aes(xend = meta2d_phase, yend = 0), size = 1,
                 lineend = "butt") +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    facet_wrap(vars(genotype)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "mRNA RPKM",
         x="") +
    theme_cowplot(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
ggsave("figures/2022-02-08/phaseAnalysis_mRNA.pdf", p1, width=10,
       height=6)

p2 <- ggplot(subset(phaseData, type=="mRNA" & meta2d_pvalue < 0.05),
       aes(x=meta2d_phase, y=1, col=Name, linetype=genotype)) +
    geom_segment(aes(xend = meta2d_phase, yend = 0), size = 1,
                 lineend = "butt") +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "mRNA RPKM",
         x="") +
    theme_cowplot(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
ggsave("figures/2022-02-08/phaseAnalysis_mRNA_merged.pdf", p2, width=7,
       height=6)


ggplot(subset(phaseData, type=="premRNA" & meta2d_pvalue < 0.05),
       aes(x=meta2d_phase, y=1, col=Name, linetype=genotype)) +
    geom_segment(aes(xend = meta2d_phase, yend = 0), size = 1,
                 lineend = "butt") +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "pre-mRNA RPKM",
         x="") +
    theme_cowplot(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
## ggsave("figures/2022-02-08/phaseAnalysis_mRNA_merged.pdf", p2, width=7,
       height=6)


boundsWT <- phaseData %>%
    subset(genotype == "WT" & meta2d_pvalue < 0.05) %>%
    pivot_wider(id_cols = Name, names_from = type,
                values_from = meta2d_phase) %>%
    mutate(
        xmax = pmax(mRNA, premRNA),
        xmin = pmin(mRNA, premRNA),
        fill = mRNA >= premRNA
    )
boundsWT
boundsWT[boundsWT$Name == "Bmal1", "xmax"] <- 24
boundsWT[boundsWT$Name == "Bmal1", "xmin"] <- 21
boundsWT[boundsWT$Name == "Clock", "xmax"] <- 24
boundsWT[boundsWT$Name == "Clock", "xmin"] <- 20.2
boundsWT[boundsWT$Name == "Npas2", "xmax"] <- 26
boundsWT[boundsWT$Name == "Npas2", "xmin"] <- 22

p3 <- ggplot() +
    geom_segment(data=subset(phaseData,
                             genotype == "WT" & meta2d_pvalue < 0.05),
                 aes(x = meta2d_phase, y = 1, xend=meta2d_phase,
                     yend = 0, col=Name, linetype=type),
                 size = 1, lineend = "butt") +
    geom_rect(data = boundsWT,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                  fill = Name),
              alpha = 0.1) +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "WT RPKM",
         x="") +
    theme_half_open(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
ggsave("figures/2022-02-08/phaseAnalysis_mRNA-premRNA_WT.pdf", p3,
       width=7, height=6)

######################################################################
## Plot one against the other (only the significant):
phaseMrnaSigMut <- rownames(phaseMrnaMut)[which(phaseMrnaMut$meta2d_BH.Q < 0.05)]
phaseMrnaSigWt <- rownames(phaseMrnaWt)[which(phaseMrnaWt$meta2d_BH.Q < 0.05)]

phaseMrnaSigAll <- unique(c(phaseMrnaSigMut, phaseMrnaSigWt))
length(phaseMrnaSigAll)

phaseMrnaSigCC <- coreClock[which(coreClock[,2] %in% phaseMrnaSigAll),2]
nameCC <- coreClock[which(coreClock[,2] %in% phaseMrnaSigAll),1]

phaseMrnaData <- data.frame(smg6 = phaseMrnaMut[phaseMrnaSigAll, 'meta2d_phase'],
                        wt = phaseMrnaWt[phaseMrnaSigAll, 'meta2d_phase'],
                        smg6Amp = phaseMrnaMut[phaseMrnaSigAll, 'meta2d_AMP'],
                        wtAmp = phaseMrnaWt[phaseMrnaSigAll, 'meta2d_AMP'])
phaseMrnaData$AMP <- (phaseMrnaData$smg6Amp + phaseMrnaData$wtAmp)/2
rownames(phaseMrnaData) <- phaseMrnaSigAll

p4 <- ggplot(phaseMrnaData, aes(x=wt, y=smg6)) +
    geom_point(alpha=0.5) +
    geom_abline(slope=1, intercept=0, linetype='dashed',
                col='darkgrey') +
    geom_point(data=phaseMrnaData[phaseMrnaSigCC,],
               aes(x=wt, y=smg6, col=AMP), col='red') +
    geom_label_repel(data=phaseMrnaData[phaseMrnaSigCC,],
                     aes(label = nameCC), fill='gray',
                     color = 'black', size = 3.5, box.padding = 0.5,
                     max.overlaps = Inf) +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,4)) +
    scale_y_continuous(limits=c(0, 24), breaks=seq(0,24,4)) +
    labs(title='Phase analysis',
         subtitle='mRNA - RPKM',
         x='WT phase',
         y='smg6 phase') +
    theme_half_open(14) +
    background_grid()
ggsave('figures/2022-02-08/phase_wtVsMut.pdf', p4, width=6, height=6)


######################################################################
### Restrict the analysis to a smaller set of core clock genes:
p10 <- ggplot(subset(phaseData, Name %in% targetClock & type == "mRNA"),
             aes(x=meta2d_phase, y=1, col=Name)) +
    geom_segment(aes(xend = meta2d_phase, yend = 0), size = 1,
                 lineend = "butt") +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    facet_wrap(vars(genotype)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "mRNA RPKM",
         x="") +
    theme_cowplot(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
ggsave("figures/2022-02-08/phaseAnalysis_mRNA_subset.pdf", p10,
       width=10, height=6)


## merge them together and add the shade
boundsWT <- phaseData %>%
    subset(type=="mRNA" & Name %in% targetClock) %>%
    pivot_wider(id_cols = Name, names_from = genotype,
                values_from = meta2d_phase) %>%
    mutate(
        xmax = smg6,
        xmin = WT,
        fill = smg6 >= WT
    )
boundsWT$y <- 1
boundsWT[which(boundsWT$Name == 'Nr1d2'), 'y'] <- 0.95
p20 <- ggplot() +
    geom_segment(data=subset(phaseData,
                             type=="mRNA" & Name %in% targetClock),
                 aes(x = meta2d_phase, y = 0, xend=meta2d_phase,
                     yend = 1, col=Name, linetype=genotype),
                 size = 1, lineend = "butt") +
    geom_rect(data = boundsWT,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y,
                  fill = Name),
              alpha = 0.1) +
    geom_segment(data = subset(boundsWT, xmin + 0.1 < xmax),
                 aes(x = xmin, xend = xmax, y = y, yend = y),
                 lineend = 'square', linejoin = 'square',
                 size = 1, arrow = arrow(length = unit(0.1, "inches"))) +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "mRNA RPKM",
         x="") +
    theme_cowplot(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
ggsave("figures/2022-02-08/phaseAnalysis_mRNA_subset_merged.pdf", p20,
       width=8, height=7)


## The same for pre-mRNA
boundsWT <- phaseData %>%
    subset(type=="premRNA" & Name %in% targetClock) %>%
    pivot_wider(id_cols = Name, names_from = genotype,
                values_from = meta2d_phase) %>%
    mutate(
        xmax = smg6,
        xmin = WT,
        fill = smg6 >= WT
    )
boundsWT$y <- seq(0.6, 1, length.out=dim(boundsWT)[1])
boundsWT[which(boundsWT$Name == "Bmal1"), "xmax"] <- 24
boundsWT[which(boundsWT$Name == "Clock"), "xmax"] <- 24
boundsWT <- rbind(boundsWT, boundsWT[which(boundsWT$Name == "Clock"),])
boundsWT[10,'xmin'] <- 0
boundsWT[10,'xmax'] <- 1
p20 <- ggplot() +
    geom_segment(data=subset(phaseData,
                             type=="premRNA" & Name %in% targetClock),
                 aes(x = meta2d_phase, y = 0, xend=meta2d_phase,
                     yend = 1, col=Name, linetype=genotype),
                 size = 1, lineend = "butt") +
    geom_rect(data = boundsWT,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y,
                  fill = Name),
              alpha = 0.1) +
    geom_segment(data = subset(boundsWT, abs(xmax-xmin) > 0.1),
                 aes(x = xmin, xend = xmax, y = y, yend = y),
                 lineend = 'square', linejoin = 'square',
                 size = 1,
                 arrow = arrow(length = unit(c(0, rep(0.1,9)),
                                             "inches"))) +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "pre-mRNA RPKM",
         x="") +
    theme_cowplot(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
ggsave("figures/2022-02-08/phaseAnalysis_premRNA_subset_merged.pdf", p20,
       width=8, height=7)


######################################################################
### Split data by genotype
## Start with the WT
boundsWT <- phaseData %>%
    subset(genotype == "WT" & Name %in% targetClock) %>%
    pivot_wider(id_cols = Name, names_from = type,
                values_from = meta2d_phase) %>%
    mutate(
        xmax = mRNA,
        xmin = premRNA,
        fill = mRNA >= premRNA
    )
boundsWT$y <- seq(0.6, 1, length.out=dim(boundsWT)[1])
boundsWT[which(boundsWT$Name == "Bmal1"), "xmax"] <- 24
boundsWT[which(boundsWT$Name == "Clock"), "xmax"] <- 24
p30 <- ggplot() +
    geom_segment(data=subset(phaseData,
                             genotype == "WT" & Name %in% targetClock),
                 aes(x = meta2d_phase, y = 0, xend=meta2d_phase,
                     yend = 1, col=Name, linetype=type),
                 size = 1, lineend = "butt") +
    geom_rect(data = boundsWT,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y,
                  fill = Name),
              alpha = 0.1) +
    geom_segment(data = subset(boundsWT, abs(xmax-xmin) > 0.1),
                 aes(x = xmin, xend = xmax, y = y, yend = y),
                 lineend = 'square', linejoin = 'square',
                 size = 1,
                 arrow = arrow(length = unit(0.1, "inches"))) +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "WT RPKM",
         x="") +
    theme_cowplot(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
ggsave("figures/2022-02-08/phaseAnalysis_mRNA-premRNA_subset_WT.pdf", p30,
       width=8, height=7)

## Now smg6
boundssmg6 <- phaseData %>%
    subset(genotype == "smg6" & Name %in% targetClock) %>%
    pivot_wider(id_cols = Name, names_from = type,
                values_from = meta2d_phase) %>%
    mutate(
        xmax = mRNA,
        xmin = premRNA,
        fill = mRNA >= premRNA
    )
boundssmg6$y <- seq(0.6, 1, length.out=dim(boundssmg6)[1])
## boundssmg6[which(boundssmg6$Name == "Clock"), "xmax"] <- 24
p30 <- ggplot() +
    geom_segment(data=subset(phaseData,
                             genotype == "smg6" & Name %in% targetClock),
                 aes(x = meta2d_phase, y = 0, xend=meta2d_phase,
                     yend = 1, col=Name, linetype=type),
                 size = 1, lineend = "butt") +
    geom_rect(data = boundssmg6,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y,
                  fill = Name),
              alpha = 0.1) +
    geom_segment(data = subset(boundssmg6, abs(xmax-xmin) > 0.1),
                 aes(x = xmin, xend = xmax, y = y, yend = y),
                 lineend = 'square', linejoin = 'square',
                 size = 1,
                 arrow = arrow(length = unit(0.1, "inches"))) +
    coord_polar() +
    scale_x_continuous(limits=c(0, 24), breaks=seq(0,24,1)) +
    labs(title = "MetaCycle phase estimate",
         subtitle = "smg6 RPKM",
         x="") +
    theme_cowplot(14) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
ggsave("figures/2022-02-08/phaseAnalysis_mRNA-premRNA_subset_smg6.pdf", p30,
       width=8, height=7)



######################################################################
### Get the genes that have a periodic pre-mRNA in both WT and mutant:
colnames(phasePreMrnaWt)

highConf <- 0.1
lowConf <- 0.6

summary(phasePreMrnaWt$meta2d_BH.Q < highConf)
summary(phasePreMrnaMut$meta2d_BH.Q < highConf)
summary(phasePreMrnaWt$meta2d_BH.Q < highConf | phasePreMrnaMut$meta2d_BH.Q < highConf)

periodicGenes <- rownames(phasePreMrnaWt)[which(phasePreMrnaWt$JTK_BH.Q < highConf |
                                                phasePreMrnaMut$JTK_BH.Q < highConf)
                                          ]

wtSpec <- periodicGenes[which(phaseMrnaMut[periodicGenes, "JTK_BH.Q"] > lowConf &
                              phaseMrnaWt[periodicGenes, "JTK_BH.Q"] < highConf)]

smg6Spec <- periodicGenes[which(phaseMrnaMut[periodicGenes, "JTK_BH.Q"] < highConf &
                                phaseMrnaWt[periodicGenes, "JTK_BH.Q"] > lowConf)]

perGenesUnspec <- periodicGenes[which(! periodicGenes %in%
                                      c(wtSpec, smg6Spec))]

## Order them according to pre-mRNA phase

ordGene <- perGenesUnspec[order(phaseMrnaWt[perGenesUnspec, "meta2d_phase"])]
ordGeneWtSpec <- wtSpec[order(phaseMrnaWt[wtSpec, "meta2d_phase"])]
ordGeneMutSpec <- smg6Spec[order(phaseMrnaMut[smg6Spec, "meta2d_phase"])]


sampleMutAnno <- sampleGeneAnno[samplesMut,]
sampleWtAnno <- sampleGeneAnno[samplesWt,]
rownames(sampleMutAnno)[which(sampleMutAnno$ZT == 0)]


perGenMutZmRNA <- as.data.frame(mutantmRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(mut.mRNA.ZT0 = mean(TS110, TS112, TS130, na.rm = TRUE),
              mut.mRNA.ZT4 = mean(TS118, TS120, TS122, na.rm = TRUE),
              mut.mRNA.ZT8 = mean(TS124, TS132, TS134, na.rm = TRUE),
              mut.mRNA.ZT12 = mean(TS106, TS108, TS136, na.rm = TRUE),
              mut.mRNA.ZT16 = mean(TS114, TS116, TS126, na.rm = TRUE),
              mut.mRNA.ZT20 = mean(TS128, TS138, TS140, na.rm = TRUE)
              ) %>%
    zScore()

perGenMutZpre <- as.data.frame(mutantpremRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(mut.pre.ZT0 = mean(TS110, TS112, TS130, na.rm = TRUE),
              mut.pre.ZT4 = mean(TS118, TS120, TS122, na.rm = TRUE),
              mut.pre.ZT8 = mean(TS124, TS132, TS134, na.rm = TRUE),
              mut.pre.ZT12 = mean(TS106, TS108, TS136, na.rm = TRUE),
              mut.pre.ZT16 = mean(TS114, TS116, TS126, na.rm = TRUE),
              mut.pre.ZT20 = mean(TS128, TS138, TS140, na.rm = TRUE)
              ) %>%
    zScore()


rownames(sampleWtAnno)[which(sampleWtAnno$ZT == 0)]

perGenWtZmRNA <- as.data.frame(wtmRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(wt.mRNA.ZT0 = mean(TS109, TS111, TS129, na.rm = TRUE),
              wt.mRNA.ZT4 = mean(TS117, TS119, TS121, na.rm = TRUE),
              wt.mRNA.ZT8 = mean(TS123, TS131, TS133, na.rm = TRUE),
              wt.mRNA.ZT12 = mean(TS105, TS107, TS135, na.rm = TRUE),
              wt.mRNA.ZT16 = mean(TS113, TS115, TS125, na.rm = TRUE),
              wt.mRNA.ZT20 = mean(TS127, TS137, TS139, na.rm = TRUE)
              ) %>%
    zScore()

perGenWtZpre <- as.data.frame(wtpremRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(wt.pre.ZT0 = mean(TS109, TS111, TS129, na.rm = TRUE),
              wt.pre.ZT4 = mean(TS117, TS119, TS121, na.rm = TRUE),
              wt.pre.ZT8 = mean(TS123, TS131, TS133, na.rm = TRUE),
              wt.pre.ZT12 = mean(TS105, TS107, TS135, na.rm = TRUE),
              wt.pre.ZT16 = mean(TS113, TS115, TS125, na.rm = TRUE),
              wt.pre.ZT20 = mean(TS127, TS137, TS139, na.rm = TRUE)
              ) %>%
    zScore()


perGenesRpkm <- cbind(perGenWtZpre, perGenWtZmRNA, perGenMutZpre,
                      perGenMutZmRNA)
rownames(perGenesRpkm) <- periodicGenes

perGenesRpkmOrd <- rbind(perGenesRpkm[ordGeneWtSpec,],
                         perGenesRpkm[ordGeneMutSpec,],
                         perGenesRpkm[ordGene,])
                         

pheatmap(perGenesRpkmOrd,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(6, 12, 18),
         gaps_row = c(length(ordGeneWtSpec),
                     length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames = FALSE,
         filename = "figures/2022-02-08/metaCycle_heatmap.pdf")

col.groups <- data.frame(Type = factor(rep(c('WT pre-mRNA', 'WT mRNA',
                                           'smg6 pre-mRNA',
                                           'smg6 mRNA'), each=6),
                                     levels=c('WT pre-mRNA',
                                              'WT mRNA',
                                              'smg6 pre-mRNA',
                                              'smg6 mRNA')))
rownames(col.groups) <- colnames(perGenesRpkmOrd)
col.names <- rep(seq(0,20,4), 4)
pheatmap(perGenesRpkmOrd,
         ## color          = getCol(40),
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         gaps_col       = c(6,12,18),
         gaps_row       = c(length(ordGeneWtSpec),
                            length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames  = FALSE,
         labels_col     = col.names,
         fontsize       = 10,
         fontsize_col   = 10,
         annotation_col = col.groups,
         annotation_names_col = FALSE,
         ## cex            = 1.2,
         filename       = 'figures/2022-02-08/periodicGenes_heatmap.pdf',
         height         = 8,
         width          = 6)


######################################################################


perGenZmRNA <- as.data.frame(mRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(mut0 = mean(TS110, TS112, TS130, na.rm = TRUE),
              mut4 = mean(TS118, TS120, TS122, na.rm = TRUE),
              mut8 = mean(TS124, TS132, TS134, na.rm = TRUE),
              mut12 = mean(TS106, TS108, TS136, na.rm = TRUE),
              mut16 = mean(TS114, TS116, TS126, na.rm = TRUE),
              mut20 = mean(TS128, TS138, TS140, na.rm = TRUE),
              wt0 = mean(TS109, TS111, TS129, na.rm = TRUE),
              wt4 = mean(TS117, TS119, TS121, na.rm = TRUE),
              wt8 = mean(TS123, TS131, TS133, na.rm = TRUE),
              wt12 = mean(TS105, TS107, TS135, na.rm = TRUE),
              wt16 = mean(TS113, TS115, TS125, na.rm = TRUE),
              wt20 = mean(TS127, TS137, TS139, na.rm = TRUE)
              ) %>%
    zScore()

perGenZpre <- as.data.frame(premRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(pre.mut0 = mean(TS110, TS112, TS130, na.rm = TRUE),
              pre.mut4 = mean(TS118, TS120, TS122, na.rm = TRUE),
              pre.mut8 = mean(TS124, TS132, TS134, na.rm = TRUE),
              pre.mut12 = mean(TS106, TS108, TS136, na.rm = TRUE),
              pre.mut16 = mean(TS114, TS116, TS126, na.rm = TRUE),
              pre.mut20 = mean(TS128, TS138, TS140, na.rm = TRUE),
              pre.wt0 = mean(TS109, TS111, TS129, na.rm = TRUE),
              pre.wt4 = mean(TS117, TS119, TS121, na.rm = TRUE),
              pre.wt8 = mean(TS123, TS131, TS133, na.rm = TRUE),
              pre.wt12 = mean(TS105, TS107, TS135, na.rm = TRUE),
              pre.wt16 = mean(TS113, TS115, TS125, na.rm = TRUE),
              pre.wt20 = mean(TS127, TS137, TS139, na.rm = TRUE),
              ) %>%
    zScore()


perGenesRpkm <- cbind(perGenZpre[,7:12], perGenZmRNA[,7:12],
                      perGenZpre[,1:6], perGenZmRNA[,1:6])
rownames(perGenesRpkm) <- periodicGenes

perGenesRpkmOrd <- rbind(perGenesRpkm[ordGeneWtSpec,],
                         perGenesRpkm[ordGeneMutSpec,],
                         perGenesRpkm[ordGene,])
                         

pheatmap(perGenesRpkmOrd,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(6, 12, 18),
         gaps_row = length(ordGeneWtSpec),
         show_rownames = FALSE,
         filename = "figures/2022-02-08/metaCycle_2_heatmap.pdf",
         width = 5
         )

col.groups <- data.frame(Type = factor(rep(c('WT pre-mRNA', 'WT mRNA',
                                           'smg6 pre-mRNA',
                                           'smg6 mRNA'), each=6),
                                     levels=c('WT pre-mRNA',
                                              'WT mRNA',
                                              'smg6 pre-mRNA',
                                              'smg6 mRNA')))
rownames(col.groups) <- colnames(perGenesRpkmOrd)
col.names <- rep(seq(0,20,4), 4)


getCol <- colorRampPalette(brewer.pal(n = 11, name = 'RdYlBu'))

pheatmap(perGenesRpkmOrd,
         ## color          = getCol(40),
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         gaps_col       = c(6,12,18),
         gaps_row       = c(length(ordGeneWtSpec),
                            length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames  = FALSE,
         labels_col     = col.names,
         fontsize       = 10,
         fontsize_col   = 10,
         annotation_col = col.groups,
         annotation_names_col = FALSE,
         ## cex            = 1.2,
         filename       = 'figures/2022-02-08/wtSpecificPeriodicGenes_heatmap.pdf',
         height         = 8,
         width          = 6)


######################################################################
### Write out the genes that have lost periodicity

gene.symbols <- mapIds(org.Mm.eg.db,
                     keys = ordGeneWtSpec,
                     column=c("SYMBOL"),
                     keytype="ENSEMBL",
                     multiVals='first')
gene.descr <- mapIds(org.Mm.eg.db,
                     keys = ordGeneWtSpec,
                     column=c("GENENAME"),
                     keytype="ENSEMBL",
                     multiVals='first')

length(gene.symbols)

dataOut <- data.frame(GID = ordGeneWtSpec,
                      Name = gene.symbols,
                      Description = gene.descr)

write.csv(dataOut,
          file = 'results/2022-02-08/smg6_lost_mRNA_period.csv')


######################################################################
### Visualize them in a different way to show if global clock has
### shifted
preMrnaInd <- which(phasePreMrnaWt$meta2d_BH.Q < 0.3 &
                    phasePreMrnaMut$meta2d_BH.Q < 0.3)

mRNAind <- which(phaseMrnaWt$meta2d_BH.Q < 0.1 &
                 phaseMrnaMut$meta2d_BH.Q < 0.1)

deltaMrna <- phaseMrnaMut$meta2d_phase[mRNAind] -
    phaseMrnaWt$meta2d_phase[mRNAind]

deltaPremrna <- phasePreMrnaMut$meta2d_phase[mRNAind] -
    phasePreMrnaWt$meta2d_phase[mRNAind]


deltaData <- data.frame(delta = c(deltaMrna, deltaPremrna),
                        type = c(rep("mRNA", length(deltaMrna)),
                                 rep("pre-mRNA", length(deltaMrna)))
                        )
head(deltaData)

dp1 <- ggplot(deltaData, aes(x = delta)) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_histogram(aes(y = ..density.., fill = type),
                   col = "darkgray", alpha = 0.4,
                   position = "identity") +
    geom_density(aes(col = type)) +
    labs(x = "Phase smg6 vs. Control",
         fill = "",
         col = "") +
    theme_cowplot(14)
ggsave("figures/2022-02-08/metaCycle_phaseDelta.pdf", dp1, height = 5)


######################################################################
######################################################################
### Get the genes that have a periodic mRNA in both WT and mutant:
######################################################################
######################################################################
colnames(phasePreMrnaWt)

highConf <- 0.1
lowConf <- 0.6

t1 <- "ENSMUSG00000024827"
t2 <- "ENSMUSG00000036752"

summary(phaseMrnaWt$meta2d_BH.Q < highConf)
summary(phaseMrnaMut$meta2d_BH.Q < highConf)
summary(phaseMrnaWt$meta2d_BH.Q < highConf | phaseMrnaMut$meta2d_BH.Q < highConf)

periodicGenes <- rownames(phaseMrnaWt)[which(phaseMrnaWt$meta2d_BH.Q < highConf |
                                             phaseMrnaMut$meta2d_BH.Q < highConf)
                                          ]
length(periodicGenes)

wtSpec <- periodicGenes[which(phaseMrnaMut[periodicGenes, "JTK_BH.Q"] > lowConf &
                              phaseMrnaWt[periodicGenes, "JTK_BH.Q"] < highConf)]

wtSpecNoPre <- wtSpec[which(phasePreMrnaMut[wtSpec, "JTK_BH.Q"] > lowConf &
                            phasePreMrnaWt[wtSpec, "JTK_BH.Q"] > lowConf)]
wtSpecPre <- wtSpec[which(! wtSpec %in% wtSpecNoPre)]
length(wtSpecNoPre)
length(wtSpecPre)
                    
smg6Spec <- periodicGenes[which(phaseMrnaMut[periodicGenes, "JTK_BH.Q"] < highConf &
                                phaseMrnaWt[periodicGenes, "JTK_BH.Q"] > lowConf)]
smg6SpecNoPre <- smg6Spec[which(phasePreMrnaMut[smg6Spec, "JTK_BH.Q"] > lowConf &
                                phasePreMrnaWt[smg6Spec, "JTK_BH.Q"] > lowConf)]
smg6SpecPre <- smg6Spec[which(! smg6Spec %in% smg6SpecNoPre)]
length(smg6SpecNoPre)
length(smg6SpecPre)

which(wtSpec == t1)
which(wtSpecNoPre == t2)

which(wtSpec == t2)
which(wtSpecNoPre == t1)

perGenesUnspec <- periodicGenes[which(! periodicGenes %in%
                                      c(wtSpec, smg6Spec))]
length(perGenesUnspec)

######################################################################
### Check Bmal1 and Erb-Reva motifs around promoters in these gene
### groups
extraCols_narrowPeak <- c(signalValue="numeric", pValue="numeric",
                          qValue="numeric", peak="integer")

baml1.pwm.hits <- import("../annotations/genomic_hit_H11MO.0.B_ARNT.bed",
                         format = "BED", genome = 'mm9')
baml1.all <- import("../publicData/takahashi12/macs/BMAL1.AllCT_peaks.narrowPeak",
                    format="BED", extraCols=extraCols_narrowPeak,
                    genome='mm9')
revErb.pwm.hits <- import("../annotations/genomic_hit_MA1531.1_NR1D1.bed",
                         format = "BED", genome = 'mm9')

cry2.p16 <- import("../publicData/takahashi12/macs/GSM982752.CRY2.CT16_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
cry2.all <- import("../publicData/takahashi12/macs/CRY2.AllCT_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')

## get top 3K Cry2 peaks:
cry2.all <- cry2.all[order(cry2.all$qValue, decreasing = TRUE)[1:3000],]


## get promoter regions
all.tss <- import('../annotations/ensembl_mm9_geneTSS.bed',
                  format='BED', genome='mm9')
names(all.tss) <- all.tss$name
all.prom <- promoters(all.tss, upstream = 10000, downstream = 1000)
pg <- periodicGenes[periodicGenes %in% all.tss$name]
pgProm <- all.prom[pg,]

pgPromBmal1 <- pgProm[queryHits(findOverlaps(pgProm, baml1.pwm.hits)),]
length(pgPromBmal1)
pgPromBmal1r <- pgProm[queryHits(findOverlaps(pgProm, baml1.all)),]
length(pgPromBmal1r)
pgPromRevErb <- pgProm[queryHits(findOverlaps(pgProm, revErb.pwm.hits)),]
length(pgPromRevErb)
pgPromCry2 <- pgProm[queryHits(findOverlaps(pgProm, cry2.all)),]
length(pgPromCry2)
pgPromCry2.16 <- pgProm[queryHits(findOverlaps(pgProm, cry2.p16)),]
length(pgPromCry2.16)

motifs <- data.frame(Bmal1Mot = factor(rep('No', length(periodicGenes)),
                                       levels = c('No', 'Yes')),
                     Bmal1Exp = factor(rep('No', length(periodicGenes)),
                                       levels = c('No', 'Yes')),
                     RevErba = factor(rep('No', length(periodicGenes)),
                                      levels = c('No', 'Yes')),
                     Cry2All = factor(rep('No', length(periodicGenes)),
                                      levels = c('No', 'Yes')),
                     Cry2.16 = factor(rep('No', length(periodicGenes)),
                                      levels = c('No', 'Yes'))
                     )
rownames(motifs) <- periodicGenes
motifs[names(pgPromBmal1), "Bmal1Mot"] <- 'Yes'
motifs[names(pgPromBmal1r), "Bmal1Exp"] <- 'Yes'
motifs[names(pgPromRevErb), "RevErba"] <- 'Yes'
motifs[names(pgPromCry2), "Cry2All"] <- 'Yes'
motifs[names(pgPromCry2.16), "Cry2.16"] <- 'Yes'

head(motifs)

######################################################################
## Order them according to mRNA phase

ordGene <- perGenesUnspec[order(phaseMrnaWt[perGenesUnspec, "meta2d_phase"])]
ordGeneWtSpec <- wtSpec[order(phaseMrnaWt[wtSpec, "meta2d_phase"])]
ordGeneWtSpecPre <- wtSpecPre[order(phaseMrnaWt[wtSpecPre,
                                                "meta2d_phase"])]
ordGeneWtSpecNoPre <- wtSpecNoPre[order(phaseMrnaWt[wtSpecNoPre,
                                                    "meta2d_phase"])]
ordGeneMutSpec <- smg6Spec[order(phaseMrnaMut[smg6Spec, "meta2d_phase"])]
ordGeneMutSpecPre <- smg6SpecPre[order(phaseMrnaWt[smg6SpecPre,
                                                "meta2d_phase"])]
ordGeneMutSpecNoPre <- smg6SpecNoPre[order(phaseMrnaWt[smg6SpecNoPre,
                                                    "meta2d_phase"])]


sampleMutAnno <- sampleGeneAnno[samplesMut,]
sampleWtAnno <- sampleGeneAnno[samplesWt,]
rownames(sampleMutAnno)[which(sampleMutAnno$ZT == 0)]


perGenMutZmRNA <- as.data.frame(mutantmRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(mut.mRNA.ZT0 = mean(TS110, TS112, TS130, na.rm = TRUE),
              mut.mRNA.ZT4 = mean(TS118, TS120, TS122, na.rm = TRUE),
              mut.mRNA.ZT8 = mean(TS124, TS132, TS134, na.rm = TRUE),
              mut.mRNA.ZT12 = mean(TS106, TS108, TS136, na.rm = TRUE),
              mut.mRNA.ZT16 = mean(TS114, TS116, TS126, na.rm = TRUE),
              mut.mRNA.ZT20 = mean(TS128, TS138, TS140, na.rm = TRUE)
              ) %>%
    zScore()

perGenMutZpre <- as.data.frame(mutantpremRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(mut.pre.ZT0 = mean(TS110, TS112, TS130, na.rm = TRUE),
              mut.pre.ZT4 = mean(TS118, TS120, TS122, na.rm = TRUE),
              mut.pre.ZT8 = mean(TS124, TS132, TS134, na.rm = TRUE),
              mut.pre.ZT12 = mean(TS106, TS108, TS136, na.rm = TRUE),
              mut.pre.ZT16 = mean(TS114, TS116, TS126, na.rm = TRUE),
              mut.pre.ZT20 = mean(TS128, TS138, TS140, na.rm = TRUE)
              ) %>%
    zScore()


rownames(sampleWtAnno)[which(sampleWtAnno$ZT == 0)]

perGenWtZmRNA <- as.data.frame(wtmRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(wt.mRNA.ZT0 = mean(TS109, TS111, TS129, na.rm = TRUE),
              wt.mRNA.ZT4 = mean(TS117, TS119, TS121, na.rm = TRUE),
              wt.mRNA.ZT8 = mean(TS123, TS131, TS133, na.rm = TRUE),
              wt.mRNA.ZT12 = mean(TS105, TS107, TS135, na.rm = TRUE),
              wt.mRNA.ZT16 = mean(TS113, TS115, TS125, na.rm = TRUE),
              wt.mRNA.ZT20 = mean(TS127, TS137, TS139, na.rm = TRUE)
              ) %>%
    zScore()

perGenWtZpre <- as.data.frame(wtpremRNArpkm[periodicGenes,]) %>%
    rowwise() %>%
    summarise(wt.pre.ZT0 = mean(TS109, TS111, TS129, na.rm = TRUE),
              wt.pre.ZT4 = mean(TS117, TS119, TS121, na.rm = TRUE),
              wt.pre.ZT8 = mean(TS123, TS131, TS133, na.rm = TRUE),
              wt.pre.ZT12 = mean(TS105, TS107, TS135, na.rm = TRUE),
              wt.pre.ZT16 = mean(TS113, TS115, TS125, na.rm = TRUE),
              wt.pre.ZT20 = mean(TS127, TS137, TS139, na.rm = TRUE)
              ) %>%
    zScore()


perGenesRpkm <- cbind(perGenWtZpre, perGenWtZmRNA, perGenMutZpre,
                      perGenMutZmRNA)
rownames(perGenesRpkm) <- periodicGenes

## perGenesRpkmOrd <- rbind(perGenesRpkm[ordGeneWtSpec,],
##                          perGenesRpkm[ordGeneMutSpec,],
##                          perGenesRpkm[ordGene,])
perGenesRpkmOrd <- rbind(perGenesRpkm[ordGeneWtSpecPre,],
                         perGenesRpkm[ordGeneWtSpecNoPre,],
                         perGenesRpkm[ordGeneMutSpecPre,],
                         perGenesRpkm[ordGeneMutSpecNoPre,],
                         perGenesRpkm[ordGene,])


pheatmap(perGenesRpkmOrd,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(6, 12, 18),
         gaps_row = c(length(ordGeneWtSpecPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre) +
                      length(ordGeneMutSpecPre),
                      length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames = FALSE,
         filename = "figures/2022-02-08/mRNA_metaCycle_heatmap.pdf")

col.groups <- data.frame(Type = factor(rep(c('WT pre-mRNA', 'WT mRNA',
                                           'smg6 pre-mRNA',
                                           'smg6 mRNA'), each=6),
                                     levels=c('WT pre-mRNA',
                                              'WT mRNA',
                                              'smg6 pre-mRNA',
                                              'smg6 mRNA')))
rownames(col.groups) <- colnames(perGenesRpkmOrd)
col.names <- rep(seq(0,20,4), 4)
my_colour <- list(
    RevErba = c(Yes = "blue", No = "white"),
    Bmal1Mot = c(Yes = "green", No = "white"),
    Bmal1Exp = c(Yes = "darkgreen", No = "white")
)
pheatmap(perGenesRpkmOrd,
         ## color          = getCol(40),
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         gaps_col       = c(6,12,18),
         gaps_row       = c(length(ordGeneWtSpecPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre) +
                      length(ordGeneMutSpecPre),
                      length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames  = FALSE,
         labels_col     = col.names,
         fontsize       = 10,
         fontsize_col   = 10,
         annotation_col = col.groups,
         annotation_row = motifs,
         annotation_names_col = FALSE,
         annotation_colors = my_colour,
         ## cex            = 1.2,
         filename       = 'figures/2022-02-08/mRNA_periodicGenes_heatmap.pdf',
         height         = 8,
         width          = 6)




######################################################################
### Calculate Z-scores for mRNA together and for pre-mRNA together
######################################################################


perGenZmRNA <- zScore(data.frame(
    mut0 = apply(mRNArpkm[periodicGenes,c("TS110", "TS112", "TS130")], 1, median),
    mut4 = apply(mRNArpkm[periodicGenes,c("TS118", "TS120", "TS122")], 1, median),
    mut8 = apply(mRNArpkm[periodicGenes,c("TS124", "TS132", "TS134")], 1, median),
    mut12 = apply(mRNArpkm[periodicGenes,c("TS106", "TS108", "TS136")], 1, median),
    mut16 = apply(mRNArpkm[periodicGenes,c("TS114", "TS116", "TS126")], 1, median),
    mut20 = apply(mRNArpkm[periodicGenes,c("TS128", "TS138", "TS140")], 1, median),
    wt0 = apply(mRNArpkm[periodicGenes,c("TS109", "TS111", "TS129")], 1, median),
    wt4 = apply(mRNArpkm[periodicGenes,c("TS117", "TS119", "TS121")], 1, median),
    wt8 = apply(mRNArpkm[periodicGenes,c("TS123", "TS131", "TS133")], 1, median),
    wt12 = apply(mRNArpkm[periodicGenes,c("TS105", "TS107", "TS135")], 1, median),
    wt16 = apply(mRNArpkm[periodicGenes,c("TS113", "TS115", "TS125")], 1, median),
    wt20 = apply(mRNArpkm[periodicGenes,c("TS127", "TS137", "TS139")], 1, median)
))

perGenZpre <- zScore(data.frame(
    mut0p = apply(premRNArpkm[periodicGenes,c("TS110", "TS112", "TS130")], 1, median),
    mut4p = apply(premRNArpkm[periodicGenes,c("TS118", "TS120", "TS122")], 1, median),
    mut8p = apply(premRNArpkm[periodicGenes,c("TS124", "TS132", "TS134")], 1, median),
    mut12p = apply(premRNArpkm[periodicGenes,c("TS106", "TS108", "TS136")], 1, median),
    mut16p = apply(premRNArpkm[periodicGenes,c("TS114", "TS116", "TS126")], 1, median),
    mut20p = apply(premRNArpkm[periodicGenes,c("TS128", "TS138", "TS140")], 1, median),
    wt0p = apply(premRNArpkm[periodicGenes,c("TS109", "TS111", "TS129")], 1, median),
    wt4p = apply(premRNArpkm[periodicGenes,c("TS117", "TS119", "TS121")], 1, median),
    wt8p = apply(premRNArpkm[periodicGenes,c("TS123", "TS131", "TS133")], 1, median),
    wt12p = apply(premRNArpkm[periodicGenes,c("TS105", "TS107", "TS135")], 1, median),
    wt16p = apply(premRNArpkm[periodicGenes,c("TS113", "TS115", "TS125")], 1, median),
    wt20p = apply(premRNArpkm[periodicGenes,c("TS127", "TS137", "TS139")], 1, median)
))

perGenesRpkm <- cbind(perGenZpre[,7:12], perGenZmRNA[,7:12],
                      perGenZpre[,1:6], perGenZmRNA[,1:6])
rownames(perGenesRpkm) <- periodicGenes

perGenesRpkmOrd <- rbind(perGenesRpkm[ordGeneWtSpecPre,],
                         perGenesRpkm[ordGeneWtSpecNoPre,],
                         perGenesRpkm[ordGeneMutSpecPre,],
                         perGenesRpkm[ordGeneMutSpecNoPre,],
                         perGenesRpkm[ordGene,])

pheatmap(perGenesRpkmOrd,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(6, 12, 18),
         gaps_row = c(length(ordGeneWtSpecPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre) +
                      length(ordGeneMutSpecPre),
                      length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames = FALSE,
         filename = "figures/2022-02-08/mRNA_metaCycle_2_heatmap.pdf",
         width = 5
         )

col.groups <- data.frame(Type = factor(rep(c('WT pre-mRNA', 'WT mRNA',
                                           'smg6 pre-mRNA',
                                           'smg6 mRNA'), each=6),
                                     levels=c('WT pre-mRNA',
                                              'WT mRNA',
                                              'smg6 pre-mRNA',
                                              'smg6 mRNA')))
rownames(col.groups) <- colnames(perGenesRpkmOrd)
col.names <- rep(seq(0,20,4), 4)


getCol <- colorRampPalette(brewer.pal(n = 11, name = 'RdYlBu'))

pheatmap(perGenesRpkmOrd,
         ## color          = getCol(40),
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         gaps_col       = c(6,12,18),
         gaps_row       = c(length(ordGeneWtSpecPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre) +
                      length(ordGeneMutSpecPre),
                      length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames  = FALSE,
         labels_col     = col.names,
         fontsize       = 10,
         fontsize_col   = 10,
         annotation_col = col.groups,
         annotation_names_col = FALSE,
         ## cex            = 1.2,
         filename       = 'figures/2022-02-08/mRNA_wtSpecificPeriodicGenes_heatmap.pdf',
         height         = 8,
         width          = 6)


### The same plot but with only 3 groups:

perGenesRpkmOrd <- rbind(perGenesRpkm[ordGeneWtSpec,],
                         perGenesRpkm[ordGeneMutSpec,],
                         perGenesRpkm[ordGene,])


pheatmap(perGenesRpkmOrd,
         ## color          = getCol(40),
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         gaps_col       = c(6,12,18),
         gaps_row       = c(length(ordGeneWtSpec),
                            length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames  = FALSE,
         labels_col     = col.names,
         fontsize       = 10,
         fontsize_col   = 10,
         annotation_col = col.groups,
         annotation_names_col = FALSE,
         ## cex            = 1.2,
         filename       = 'figures/2022-02-08/mRNA_wtSpecificPeriodicGenes_3groups_heatmap.pdf',
         height         = 8,
         width          = 6)


######################################################################
### Check genes groups for phase and expression levels:
length(ordGene)
summary(perGenesRpkm[ordGeneWtSpec,])

periodicWT <- data.frame(GeneID = ordGeneWtSpec,
                          perGenesRpkm[ordGeneWtSpec,])

periodicWTL <- pivot_longer(periodicWT, cols = -c(GeneID),
                             names_to = "Samples",
                             values_to = "Zscores")
periodicWTL$type <- factor("mRNA", levels = c("mRNA", "pre-mRNA"))
periodicWTL$type[grep('p', periodicWTL$Samples)] <- "pre-mRNA"
periodicWTL$geno <- factor("WT", levels = c("WT", "MUT"))
periodicWTL$geno[grep('mut', periodicWTL$Samples)] <- "MUT"
periodicWTL

pex1 <- ggplot(subset(periodicWTL, type == "mRNA"),
               aes(y = Zscores, x = geno, fill = geno)) +
    geom_boxplot() +
    labs(title = "Genes that have lost periodicity",
         x = "",
         y = "Z scores") +
    theme(legend.position = 'none')
ggsave("figures/2022-02-08/mutantLostPer_boxplot.pdf", pex1)

### Check the Mutant specific genes
periodicMut <- data.frame(GeneID = ordGeneMutSpec,
                          perGenesRpkm[ordGeneMutSpec,])

periodicMutL <- pivot_longer(periodicMut, cols = -c(GeneID),
                             names_to = "Samples",
                             values_to = "Zscores")
periodicMutL$type <- factor("mRNA", levels = c("mRNA", "pre-mRNA"))
periodicMutL$type[grep('p', periodicMutL$Samples)] <- "pre-mRNA"
periodicMutL$geno <- factor("WT", levels = c("WT", "MUT"))
periodicMutL$geno[grep('mut', periodicMutL$Samples)] <- "MUT"
periodicMutL

pex1 <- ggplot(subset(periodicMutL, type == "mRNA"),
               aes(y = Zscores, x = geno, fill = geno)) +
    geom_boxplot() +
    labs(title = "Genes that have gained periodicity",
         x = "",
         y = "Z scores") +
    theme(legend.position = 'none')
ggsave("figures/2022-02-08/mutantGainedPer_boxplot.pdf", pex1)


### Check genes that are periodic in both
periodicAll <- data.frame(GeneID = ordGene,
                          perGenesRpkm[ordGene,])

periodicAllL <- pivot_longer(periodicAll, cols = -c(GeneID),
                             names_to = "Samples",
                             values_to = "Zscores")
periodicAllL$type <- factor("mRNA", levels = c("mRNA", "pre-mRNA"))
periodicAllL$type[grep('p', periodicAllL$Samples)] <- "pre-mRNA"
periodicAllL$geno <- factor("WT", levels = c("WT", "MUT"))
periodicAllL$geno[grep('mut', periodicAllL$Samples)] <- "MUT"
periodicAllL

pex2 <- ggplot(subset(periodicAllL, type == "mRNA"),
               aes(y = Zscores, x = geno, fill = geno)) +
    geom_boxplot() +
    labs(title = "Periodic genes in WT and Mutant",
         subtitle = "mRNA",
         x = "",
         y = "Z scores") +
    theme(legend.position = 'none')
ggsave("figures/2022-02-08/mutantWtPer_boxplot.pdf", pex2)


### Now pre-mRNA
pex2 <- ggplot(subset(periodicAllL, type == "pre-mRNA"),
               aes(y = Zscores, x = geno, fill = geno)) +
    geom_boxplot() +
    labs(title = "Periodic genes in WT and Mutant",
         subtitle = "pre-mRNA",
         x = "",
         y = "Z scores") +
    theme(legend.position = 'none')
ggsave("figures/2022-02-08/mutantWtPer_premRNA_boxplot.pdf", pex2)


## Instead of looking at all the values take only the range:
rangeDiff <- function(x){r <- range(x, na.rm = TRUE)
    return(r[2] - r[1])}


periodicAll <- data.frame(GeneID = ordGene,
                          perGenesRpkm[ordGene,])

dim(periodicAll)

wtmRNA <- grep(pattern = "wt[0-9]*$", colnames(periodicAll))
mutmRNA <- grep(pattern = "mut[0-9]*$", colnames(periodicAll))

periodicAllAmp <- data.frame(
    amplitude = c(apply(periodicAll[,wtmRNA], 1, rangeDiff),
                  apply(periodicAll[,mutmRNA], 1, rangeDiff)),
    genotype = rep(c("WT", "MUT"), each = dim(periodicAll)[1])
    )

pex3 <- ggplot(periodicAllAmp,
               aes(y = amplitude, x = genotype, fill = genotype)) +
    geom_boxplot() +
    labs(title = "Periodic genes in WT and Mutant",
         subtitle = "mRNA",
         x = "",
         y = "Z scores amplitudes") +
    theme_light(14) +
    theme(legend.position = 'none') 
ggsave("figures/2022-02-08/mutantWtPer_Zamp_boxplot.pdf", pex3)

summary(lm(amplitude ~ genotype, data = subset(periodicAllAmp, amplitude >= 0)))


wtPmRNA <- grep(pattern = "wt[0-9]*p", colnames(periodicAll))
mutPmRNA <- grep(pattern = "mut[0-9]*p", colnames(periodicAll))

periodicPreAmp <- data.frame(
    amplitude = c(apply(periodicAll[,wtPmRNA], 1, rangeDiff),
                  apply(periodicAll[,mutPmRNA], 1, rangeDiff)),
    genotype = rep(c("WT", "MUT"), each = dim(periodicAll)[1])
    )

pex4 <- ggplot(periodicPreAmp,
               aes(y = amplitude, x = genotype, fill = genotype)) +
    geom_boxplot() +
    labs(title = "Periodic genes in WT and Mutant",
         subtitle = "pre-mRNA",
         x = "",
         y = "Z scores amplitudes") +
    theme_light(14) +
    theme(legend.position = 'none') 
ggsave("figures/2022-02-08/mutantWtPer_Zamp_Pre_boxplot.pdf", pex4)

summary(lm(amplitude ~ genotype, data = subset(periodicPreAmp, amplitude >= 0)))


periodicAmpVar <- data.frame(ratio = periodicAllAmp$amplitude / periodicPreAmp$amplitude,
                             genotype = rep(c("WT", "MUT"), each = dim(periodicAll)[1])
                             )

pex4 <- ggplot(periodicAmpVar,
               aes(y = ratio, x = genotype, fill = genotype)) +
    geom_boxplot() +
    labs(title = "Periodic genes in WT and Mutant",
         subtitle = "mRNA vs pre-mRNA",
         x = "",
         y = "Z scores amplitude ratios") +
    theme_light(14) +
    theme(legend.position = 'none') 
ggsave("figures/2022-02-08/mutantWtPer_Zamp_ratio_boxplot.pdf", pex4)

### The difference is highly significant
summary(lm(ratio ~ genotype, data = subset(periodicAmpVar, ratio < 6)))



######################################################################
### Make a  circular plot with the phases
head(phaseMrnaMut[ordGene, "meta2d_phase"])

phaseData <- data.frame(phase = c(phaseMrnaMut[ordGene, "meta2d_phase"],
                                  phaseMrnaWt[ordGene, "meta2d_phase"]),
                        Amplitude = c(phaseMrnaWt[ordGene, "meta2d_AMP"],
                                      phaseMrnaMut[ordGene, "meta2d_AMP"]),
                        geno = factor(rep(c("MUT", "WT"),
                                          each = length(ordGene)))
                        )

ph1 <- ggplot(phaseData, aes(x = phase, fill = geno)) +
    geom_histogram(binwidth = 1) +
    scale_x_continuous(breaks = seq(0, 20, by = 4), limits = c(0, 24)) +
    coord_polar() +
    facet_wrap( ~ geno) +
    theme_light(14)
ggsave("figures/2022-02-08/phaseCircularPlot.pdf", ph1, height = 5,
       width = 7)

ph1 <- ggplot(subset(phaseData, Amplitude > 2),
              aes(x = phase, fill = geno)) +
    geom_histogram(binwidth = 1) +
    scale_x_continuous(breaks = seq(0, 20, by = 4), limits = c(0, 24)) +
    coord_polar() +
    facet_wrap( ~ geno) +
    theme_light(14)
ggsave("figures/2022-02-08/phaseCircularPlot_highAmplitude.pdf", ph1,
       height = 5, width = 7)

### Plot the amplitude
pa1 <- ggplot(phaseData,
              aes(y = Amplitude, x = geno, fill = geno)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(title = "Periodic genes in WT and Mutant",
         subtitle = "mRNA",
         x = "",
         y = "Expression Amplitudes") +
    theme_light(14) +
    theme(legend.position = 'none') 
ggsave("figures/2022-02-08/mutantWtPer_amp_boxplot.pdf", pa1)




######################################################################
### Calculate Z-scores globally
######################################################################

premRNArpkmMod <- premRNArpkm
colnames(premRNArpkmMod) <- paste("p", colnames(premRNArpkm), sep="")

perGenZ <- as.data.frame(cbind(mRNArpkm[periodicGenes,],
                               premRNArpkmMod[periodicGenes,])) %>%
    rowwise() %>%
    summarise(pre.wt0 = mean(pTS109, pTS111, pTS129, na.rm = TRUE),
              pre.wt4 = mean(pTS117, pTS119, pTS121, na.rm = TRUE),
              pre.wt8 = mean(pTS123, pTS131, pTS133, na.rm = TRUE),
              pre.wt12 = mean(pTS105, pTS107, pTS135, na.rm = TRUE),
              pre.wt16 = mean(pTS113, pTS115, pTS125, na.rm = TRUE),
              pre.wt20 = mean(pTS127, pTS137, pTS139, na.rm = TRUE),
              wt0 = mean(TS109, TS111, TS129, na.rm = TRUE),
              wt4 = mean(TS117, TS119, TS121, na.rm = TRUE),
              wt8 = mean(TS123, TS131, TS133, na.rm = TRUE),
              wt12 = mean(TS105, TS107, TS135, na.rm = TRUE),
              wt16 = mean(TS113, TS115, TS125, na.rm = TRUE),
              wt20 = mean(TS127, TS137, TS139, na.rm = TRUE),
              pre.mut0 = mean(pTS110, pTS112, pTS130, na.rm = TRUE),
              pre.mut4 = mean(pTS118, pTS120, pTS122, na.rm = TRUE),
              pre.mut8 = mean(pTS124, pTS132, pTS134, na.rm = TRUE),
              pre.mut12 = mean(pTS106, pTS108, pTS136, na.rm = TRUE),
              pre.mut16 = mean(pTS114, pTS116, pTS126, na.rm = TRUE),
              pre.mut20 = mean(pTS128, pTS138, pTS140, na.rm = TRUE),
              mut0 = mean(TS110, TS112, TS130, na.rm = TRUE),
              mut4 = mean(TS118, TS120, TS122, na.rm = TRUE),
              mut8 = mean(TS124, TS132, TS134, na.rm = TRUE),
              mut12 = mean(TS106, TS108, TS136, na.rm = TRUE),
              mut16 = mean(TS114, TS116, TS126, na.rm = TRUE),
              mut20 = mean(TS128, TS138, TS140, na.rm = TRUE)
              ) %>%
    zScore()

perGenesRpkm <- perGenZ
rownames(perGenesRpkm) <- periodicGenes

perGenesRpkmOrd <- rbind(perGenesRpkm[ordGeneWtSpecPre,],
                         perGenesRpkm[ordGeneWtSpecNoPre,],
                         perGenesRpkm[ordGeneMutSpecPre,],
                         perGenesRpkm[ordGeneMutSpecNoPre,],
                         perGenesRpkm[ordGene,])
                         

pheatmap(perGenesRpkmOrd,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(6, 12, 18),
         gaps_row = c(length(ordGeneWtSpecPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre) +
                      length(ordGeneMutSpecPre),
                      length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames = FALSE,
         filename = "figures/2022-02-08/mRNA_metaCycle_3_heatmap.pdf",
         width = 5
         )

col.groups <- data.frame(Type = factor(rep(c('WT pre-mRNA', 'WT mRNA',
                                           'smg6 pre-mRNA',
                                           'smg6 mRNA'), each=6),
                                     levels=c('WT pre-mRNA',
                                              'WT mRNA',
                                              'smg6 pre-mRNA',
                                              'smg6 mRNA')))
rownames(col.groups) <- colnames(perGenesRpkmOrd)
col.names <- rep(seq(0,20,4), 4)


getCol <- colorRampPalette(brewer.pal(n = 11, name = 'RdYlBu'))

pheatmap(perGenesRpkmOrd,
         ## color          = getCol(40),
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         gaps_col       = c(6,12,18),
         gaps_row       = c(length(ordGeneWtSpecPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre),
                      length(ordGeneWtSpecPre) + length(ordGeneWtSpecNoPre) +
                      length(ordGeneMutSpecPre),
                      length(ordGeneWtSpec) + length(ordGeneMutSpec)),
         show_rownames  = FALSE,
         labels_col     = col.names,
         fontsize       = 10,
         fontsize_col   = 10,
         annotation_col = col.groups,
         annotation_names_col = FALSE,
         ## cex            = 1.2,
         filename       = 'figures/2022-02-08/mRNA_wtSpecificPeriodicGenes_3_heatmap.pdf',
         height         = 10,
         width          = 6)


######################################################################
### Make a scatter plot of the phases, ordered by WT mRNA phase
######################################################################
pgPromBmal1
pgPromBmal1r
pgPromRevErb
pgPromCry2
pgPromCry2.16


yVal <- 1:length(ordGene)

phaseData <- data.frame(Rank  = rep(yVal, 2),
                        Phase = c(phaseMrnaWt[ordGene, "meta2d_phase"],
                                  phaseMrnaMut[ordGene, "meta2d_phase"]),
                        Amplitude = c(phaseMrnaWt[ordGene, "meta2d_AMP"],
                                      phaseMrnaMut[ordGene, "meta2d_AMP"]),
                        Type  = rep(c("WT", "Mut"), each = length(ordGene))
                        )

phaseDataMut <- data.frame(Rank = yVal,
                           Phase = phaseMrnaMut[ordGene, "meta2d_phase"])
rownames(phaseDataMut) <- ordGene

sp1 <- ggplot(phaseData, aes(x = Phase, y = Rank, col = Type)) +
    geom_point(alpha = 0.5) +
    geom_point(data = phaseDataMut[names(pgPromBmal1r),],
               aes(x = Phase, y = Rank), col = "darkgreen") +
    labs(title = "Bmal1 ChIP highlighted")
ggsave("figures/2022-02-08/phaseDiffWtMut_Bmal1Chip.pdf", width = 6, height = 6)

sp1 <- ggplot(phaseData, aes(x = Phase, y = Rank, col = Type)) +
    geom_point(alpha = 0.5) +
    geom_point(data = phaseDataMut[names(pgPromCry2),],
               aes(x = Phase, y = Rank), col = "darkgreen") +
    labs(title = "Cry2 highlighted",
         subtitle = "Top 3000 peaks from all ZT")
ggsave("figures/2022-02-08/phaseDiffWtMut_Cry2.pdf", width = 6, height = 6)

sp1 <- ggplot(phaseData, aes(x = Phase, y = Rank, col = Type)) +
    geom_point(alpha = 0.5) +
    geom_point(data = phaseDataMut[names(pgPromCry2.16),],
               aes(x = Phase, y = Rank), col = "darkgreen") +
    labs(title = "Cry2 highlighted",
         subtitle = "All peaks from ZT16")
ggsave("figures/2022-02-08/phaseDiffWtMut_Cry2zt16.pdf", width = 6, height = 6)

### Restrict to genes with high amplitude:
sp1 <- ggplot(subset(phaseData, Amplitude > 2),
              aes(x = Phase, y = Rank, col = Type)) +
    geom_point(alpha = 0.5) +
    geom_point(data = phaseDataMut[names(pgPromBmal1r),],
               aes(x = Phase, y = Rank), col = "darkgreen") +
    labs(title = "Bmal1 ChIP highlighted",
         subtitle = "Amplitude > 2")
ggsave("figures/2022-02-08/phaseDiffWtMut_amp_Bmal1Chip.pdf", width = 6, height = 6)

sp1 <- ggplot(subset(phaseData, Amplitude > 2),
              aes(x = Phase, y = Rank, col = Type)) +
    geom_point(alpha = 0.5) +
    geom_point(data = phaseDataMut[names(pgPromCry2),],
               aes(x = Phase, y = Rank), col = "darkgreen") +
    labs(title = "Cry2 highlighted",
         subtitle = "Top 3000 peaks all ZT, Amp > 2")
ggsave("figures/2022-02-08/phaseDiffWtMut_amp_Cry2.pdf", width = 6, height = 6)


######################################################################
### Do the same type of plot but with pre-mRNA vs mRNA
######################################################################
prePhase <- perGenesUnspec[which(phasePreMrnaWt[perGenesUnspec, "meta2d_BH.Q"] < 0.3 |
                                 phasePreMrnaMut[perGenesUnspec, "meta2d_BH.Q"] < 0.3)]


ordGene2 <- prePhase[order(phaseMrnaWt[prePhase, "meta2d_phase"])]
ordGene3 <- prePhase[order(phaseMrnaMut[prePhase, "meta2d_phase"])]

length(ordGene2)

yVal <- 1:length(ordGene2)
phaseData2 <- data.frame(Rank  = rep(yVal, 2),
                         Phase = c(phaseMrnaWt[ordGene2, "meta2d_phase"],
                                   phasePreMrnaWt[ordGene2, "meta2d_phase"]),
                         Type  = factor(rep(c("mRNA", "pre-mRNA"),
                                            each = length(ordGene2)),
                                        levels = c("mRNA", "pre-mRNA"))
                         )
head(phaseData2)

sp1 <- ggplot(phaseData2, aes(x = Phase, y = Rank, col = Type)) +
    geom_point(alpha = 0.5) +
    ggtitle("WT")
ggsave("figures/2022-02-08/phaseDiff_pre-mRNA-mRNA_WT.pdf", width = 6, height = 6)

yVal <- 1:length(ordGene3)
phaseData3 <- data.frame(Rank  = rep(yVal, 2),
                         Phase = c(phaseMrnaMut[ordGene3, "meta2d_phase"],
                                   phasePreMrnaMut[ordGene3, "meta2d_phase"]),
                         Type  = factor(rep(c("mRNA", "pre-mRNA"),
                                            each = length(ordGene3)),
                                        levels = c("mRNA", "pre-mRNA"))
                         )
head(phaseData3)

sp1 <- ggplot(phaseData3, aes(x = Phase, y = Rank, col = Type)) +
    geom_point(alpha = 0.5) +
    ggtitle("Mutant")
ggsave("figures/2022-02-08/phaseDiff_pre-mRNA-mRNA_Mutant.pdf", width = 6, height = 6)



######################################################################
## Add gene name and write gene groups to file
######################################################################
GSordGeneWtSpecPre <- mapIds(org.Mm.eg.db,
                             keys = ordGeneWtSpecPre,
                             column=c("SYMBOL"),
                             keytype="ENSEMBL",
                             multiVals='first')
GNordGeneWtSpecPre <- mapIds(org.Mm.eg.db,
                             keys = ordGeneWtSpecPre,
                             column=c("GENENAME"),
                             keytype="ENSEMBL",
                             multiVals='first')
ordGeneWtSpecPreOut <- data.frame(GeneName = GSordGeneWtSpecPre,
                                  Deascription = GNordGeneWtSpecPre,
                                  perGenesRpkm[ordGeneWtSpecPre,]
                                  )
write.csv(ordGeneWtSpecPreOut,
          file = "results/2022-02-08/WTmRNAperWTpremRNAperMUTmRNAnoPerMUTpremRNAper.csv")
###
GSordGeneWtSpecNoPre <- mapIds(org.Mm.eg.db,
                             keys = ordGeneWtSpecNoPre,
                             column=c("SYMBOL"),
                             keytype="ENSEMBL",
                             multiVals='first')
GNordGeneWtSpecNoPre <- mapIds(org.Mm.eg.db,
                             keys = ordGeneWtSpecNoPre,
                             column=c("GENENAME"),
                             keytype="ENSEMBL",
                             multiVals='first')
ordGeneWtSpecNoPreOut <- data.frame(GeneName = GSordGeneWtSpecNoPre,
                                  Deascription = GNordGeneWtSpecNoPre,
                                  perGenesRpkm[ordGeneWtSpecNoPre,]
                                  )
write.csv(ordGeneWtSpecNoPreOut,
          file = "results/2022-02-08/WTmRNAperWTpremRNAnoPerMUTmRNAnoPerMUTpremRNAnoPer.csv")
##

GSordGeneMutSpecPre <- mapIds(org.Mm.eg.db,
                             keys = ordGeneMutSpecPre,
                             column=c("SYMBOL"),
                             keytype="ENSEMBL",
                             multiVals='first')
GNordGeneMutSpecPre <- mapIds(org.Mm.eg.db,
                             keys = ordGeneMutSpecPre,
                             column=c("GENENAME"),
                             keytype="ENSEMBL",
                             multiVals='first')
ordGeneMutSpecPreOut <- data.frame(GeneName = GSordGeneMutSpecPre,
                                  Deascription = GNordGeneMutSpecPre,
                                  perGenesRpkm[ordGeneMutSpecPre,]
                                  )
write.csv(ordGeneMutSpecPreOut,
          file = "results/2022-02-08/WTmRNAnoPerWTpremRNAperMUTmRNAperMUTpremRNAper.csv")

##
GSordGeneMutSpecNoPre <- mapIds(org.Mm.eg.db,
                             keys = ordGeneMutSpecNoPre,
                             column=c("SYMBOL"),
                             keytype="ENSEMBL",
                             multiVals='first')
GNordGeneMutSpecNoPre <- mapIds(org.Mm.eg.db,
                             keys = ordGeneMutSpecNoPre,
                             column=c("GENENAME"),
                             keytype="ENSEMBL",
                             multiVals='first')
ordGeneMutSpecNoPreOut <- data.frame(GeneName = GSordGeneMutSpecNoPre,
                                  Deascription = GNordGeneMutSpecNoPre,
                                  perGenesRpkm[ordGeneMutSpecNoPre,]
                                  )
write.csv(ordGeneMutSpecNoPreOut,
          file = "results/2022-02-08/WTmRNAnoPerWTpremRNAperMUTmRNAperMUTpremRNAnoPer.csv")

##
GSordGene <- mapIds(org.Mm.eg.db,
                             keys = ordGene,
                    column=c("SYMBOL"),
                    keytype="ENSEMBL",
                    multiVals='first')
GNordGene <- mapIds(org.Mm.eg.db,
                    keys = ordGene,
                    column=c("GENENAME"),
                    keytype="ENSEMBL",
                    multiVals='first')
ordGeneOut <- data.frame(GeneName = GSordGene,
                         Deascription = GNordGene,
                         perGenesRpkm[ordGene,]
                         )
write.csv(ordGeneOut,
          file = "results/2022-02-08/WTmRNAperMUTmRNAper.csv")





######################################################################
### Save image:
## save.image("2022-02-08_metaCycleAnalysis.RData")
load("2022-02-08_metaCycleAnalysis.RData")
