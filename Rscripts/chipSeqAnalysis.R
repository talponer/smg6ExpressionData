######################################################################
######################################################################
### The chip-seq data (Takahashi 2021) was run around the clock. Macs
### was used to extract peaks from each time point. Moreover I merged
### all data for one target protein and extract all peaks.

library(rtracklayer)
library(GenomicRanges)

extraCols_narrowPeak <- c(signalValue="numeric", pValue="numeric",
                          qValue="numeric", peak="integer")


## get promoter regions
all.tss <- import('../annotations/ensembl_mm9_geneTSS.bed',
                  format='BED', genome='mm9')
names(all.tss) <- all.tss$name

all.prom <- promoters(all.tss, upstream = 10000, downstream = 1000)


## get the Cry2 peaks from MACS
cry2.p00 <- import("../publicData/takahashi12/macs/GSM982748.CRY2.CT0_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
cry2.p04 <- import("../publicData/takahashi12/macs/GSM982749.CRY2.CT4_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
cry2.p08 <- import("../publicData/takahashi12/macs/GSM982750.CRY2.CT8_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
cry2.p12 <- import("../publicData/takahashi12/macs/GSM982751.CRY2.CT12_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
cry2.p16 <- import("../publicData/takahashi12/macs/GSM982752.CRY2.CT16_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
cry2.p20 <- import("../publicData/takahashi12/macs/GSM982753.CRY2.CT20_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')

cry2.all <- import("../publicData/takahashi12/macs/CRY2.AllCT_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')


######################################################################
### Check the Cry2 peak quality comparing with BMAL1 peaks (calculated
### on the global data)
baml1.all <- import("../publicData/takahashi12/macs/BMAL1.AllCT_peaks.narrowPeak",
                    format="BED", extraCols=extraCols_narrowPeak,
                    genome='mm9')


cry2Bmal1.all <- cry2.all[queryHits(findOverlaps(cry2.all, baml1.all)),]
cry2NoBmal1.all <- cry2.all[-queryHits(findOverlaps(cry2.all, baml1.all)),]

bmal1Cry2.all <- baml1.all[queryHits(findOverlaps(baml1.all, cry2.all)),]
bmal1NoCry2.all <- baml1.all[-queryHits(findOverlaps(baml1.all, cry2.all)),]


######################################################################
### Find peaks in promoters

prom.cry2Bmal1 <- all.prom[queryHits(findOverlaps(all.prom, cry2Bmal1.all)),]
prom.cry2 <- all.prom[queryHits(findOverlaps(all.prom, cry2.all)),]


######################################################################
### Find peaks in genes that have periodic expression
wtPer <- read.table("results/2022-02-08/meta2d_wtmRNArpkm.txt",
                    header = TRUE, row.names = 1)
colnames(wtPer)
wtPerSig <- wtPer[wtPer$meta2d_pvalue < 0.05, ]
perGenes <- rownames(wtPerSig)
length(perGenes)

perGenes <- perGenes[perGenes %in% names(all.prom)]
per.prom <- all.prom[perGenes,]

perProm.cry2Bmal1 <- per.prom[queryHits(findOverlaps(per.prom, cry2Bmal1.all)),]
perProm.cry2 <- per.prom[queryHits(findOverlaps(per.prom, cry2.all)),]










pdf('figures/2019-10-16/Cry2peaks_numbers.pdf')
barplot(cry2PeaksNumber ~ times,
        main='Cry2 peaks', xlab="ZT time", ylab="Total number")
dev.off()


pdf('figures/2019-10-16/Cry2peaks_vs_BMAL1pwm_Freq.pdf')
barplot(cry2.bmal1pwm.perc ~ times,
        main='Frequency of Cry2 peaks with BMAL1 matches', xlab="ZT time", ylab="Frequency")
dev.off()

pdf('figures/2019-10-16/Cry2peaks_vs_BMAL1pwm.pdf', width=10)
par(mfrow=c(1,2))
barplot(cry2.bmal1pwm ~ times,
        main='Cry2 peaks with BMAL1 matches', xlab="ZT time", ylab="Total number")
barplot(cry2.bmal1pwm.perc ~ times,
        main='Frequency of Cry2 peaks with BMAL1 matches', xlab="ZT time", ylab="Frequency")
dev.off()

hnf4a.pwm.hits <- import("data/genomic_hit_MA0114.3_Hnf4a.bed",
                         format="BED", genome='mm9')

## now hnf4a
cry2.hnf4apwm <- c(length(cry2.p00[queryHits(findOverlaps(cry2.p00, hnf4a.pwm.hits)),]),
                   length(cry2.p04[queryHits(findOverlaps(cry2.p04, hnf4a.pwm.hits)),]),
                   length(cry2.p08[queryHits(findOverlaps(cry2.p08, hnf4a.pwm.hits)),]),
                   length(cry2.p12[queryHits(findOverlaps(cry2.p12, hnf4a.pwm.hits)),]),
                   length(cry2.p16[queryHits(findOverlaps(cry2.p16, hnf4a.pwm.hits)),]),
                   length(cry2.p20[queryHits(findOverlaps(cry2.p20, hnf4a.pwm.hits)),])
                   )
cry2.hnf4apwm.perc <- cry2.hnf4apwm/cry2PeaksNumber

barplot(cry2.hnf4apwm ~ colnames(cry2AllPeakInt), main='Cry2 peaks with HNF4a matches')
barplot(cry2.hnf4apwm.perc ~ colnames(cry2AllPeakInt), main='Frequency Cry2 peaks with HNF4a matches')

plot(cry2.bmal1pwm.perc ~ colnames(cry2AllPeakInt), main='Frequency of BMAL1 PWM matches in Cry2 peaks', type='b')

######################################################################
### Get the cyclic genes that are up-regulated in Cry2-ZT20
cry2.zt20.up <- read.csv('results/2018-11-22/cry2_transcription_zt20_phase.csv', row.names=1)
cyclicProm <- rownames(cry2.zt20.up)[which(cry2.zt20.up$cyclic)]

cry2.zt20.up.prom <- all.prom[cyclicProm,]


## get all Cry2 peaks (concatenate all peaks and collaps overlapping
## one)
all.cry2.peaks <- reduce(c(cry2.p00, cry2.p04, cry2.p08, cry2.p12, cry2.p16, cry2.p20))

## find all peaks that are near cry2 up-regulated genes
cry2.peaks.inUpProm <- all.cry2.peaks[queryHits(findOverlaps(all.cry2.peaks, cry2.zt20.up.prom)),]

## get their signalValues at each time-point
cry2PeakInt <- matrix(NA, ncol=6, nrow=length(cry2.peaks.inUpProm))
colnames(cry2PeakInt) <- c('00','04','08','12','16','20')
obj <- c('cry2.p00','cry2.p04','cry2.p08','cry2.p12','cry2.p16','cry2.p20')
for (I in 1:length(obj)){
    overlap <- findOverlaps(cry2.peaks.inUpProm, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    cry2PeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

cry2PeaksTime <- table(apply(cry2PeakInt, 1, which.max))

pNumb <- function(x) return(length(na.omit(x)))

cry2PeakIntClean <- cry2PeakInt[apply(cry2PeakInt, 1, pNumb) > 2, ]
cry2PeaksTimeClean <- table(apply(cry2PeakIntClean, 1, which.max))


## do the same for peks near all promoters
cry2.peaks.allProm <- all.cry2.peaks[queryHits(findOverlaps(all.cry2.peaks, all.prom)),]
cry2AllPeakInt <- matrix(NA, ncol=6, nrow=length(cry2.peaks.allProm))
colnames(cry2AllPeakInt) <- c('00','04','08','12','16','20')
obj <- c('cry2.p00','cry2.p04','cry2.p08','cry2.p12','cry2.p16','cry2.p20')
for (I in 1:length(obj)){
    overlap <- findOverlaps(cry2.peaks.allProm, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    cry2AllPeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

allPeaksTime <- table(apply(cry2AllPeakInt, 1, which.max))
allPeakIntClean <- cry2AllPeakInt[apply(cry2AllPeakInt, 1, pNumb) > 2, ]
allPeaksTimeClean <- table(apply(allPeakIntClean, 1, which.max))


pdf('figures/2019-10-16/Cry2peaks_All_phase_clean.pdf')
par(mfrow=c(1,1))
barplot(allPeaksTimeClean ~ colnames(cry2AllPeakInt), main='Phase all Cry2 peaks')
dev.off()

pdf('figures/2019-10-16/Cry2peaks_Upreg_phase_clean.pdf')
par(mfrow=c(1,1))
barplot(cry2PeaksTimeClean ~ colnames(cry2AllPeakInt), main='Phase all Cry2 peaks')
dev.off()

pdf('figures/2019-10-16/Cry2peaks_All_phase.pdf')
par(mfrow=c(1,1))
barplot(allPeaksTime ~ colnames(cry2AllPeakInt), main='Phase all Cry2 peaks')
dev.off()

pdf('figures/2019-10-16/Cry2peaks_Upreg_phase.pdf')
par(mfrow=c(1,1))
barplot(cry2PeaksTime ~ colnames(cry2AllPeakInt), main='Phase Cry2 peaks near Up regulated genes')
dev.off()

pdf('figures/2019-10-16/Cry2peaks_phase.pdf', height=10)
par(mfrow=c(2,1))
barplot(allPeaksTime ~ colnames(cry2AllPeakInt), main='Phase all Cry2 peaks')
barplot(cry2PeaksTime ~ colnames(cry2AllPeakInt), main='Phase Cry2 peaks near Up regulated genes')
dev.off()

## average signalValue of peaks
summary(cry2AllPeakInt)
summary(cry2PeakInt)

par(mfrow=c(1,2))
boxplot(cry2AllPeakInt)
boxplot(cry2PeakInt)



## ### Get the peaks that have signal at all time points (should be true
## ### peaks)
## cry2.peaks.allProm

## na.num <- function(x){
##     sum(is.na(x))
## }

## cry2AllPeakInt.num <- apply(cry2AllPeakInt, 1, na.num)
## cry2.peaks.allProm.clean <- cry2.peaks.allProm[which(cry2AllPeakInt.num<5)]

## ## repeat the distribution for all peaks
## cry2AllPeakInt <- matrix(NA, ncol=6, nrow=length(cry2.peaks.allProm.clean))
## colnames(cry2AllPeakInt) <- c('00','04','08','12','16','20')
## obj <- c('cry2.p00','cry2.p04','cry2.p08','cry2.p12','cry2.p16','cry2.p20')
## for (I in 1:length(obj)){
##     overlap <- findOverlaps(cry2.peaks.allProm.clean, get(obj[I]))
##     signalValues <- get(obj[I])$signalValue
##     cry2AllPeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
## }

## allPeaksTime <- table(apply(cry2AllPeakInt, 1, which.max))

## ## and for the one near up-regulated promoters 
## cry2PeakInt.num <- apply(cry2PeakInt, 1, na.num)
## cry2.peaks.inUpProm.clean <- cry2.peaks.inUpProm[which(cry2PeakInt.num<5)]
## cry2PeakInt <- matrix(NA, ncol=6, nrow=length(cry2.peaks.inUpProm.clean))
## colnames(cry2PeakInt) <- c('00','04','08','12','16','20')
## obj <- c('cry2.p00','cry2.p04','cry2.p08','cry2.p12','cry2.p16','cry2.p20')
## for (I in 1:length(obj)){
##     overlap <- findOverlaps(cry2.peaks.inUpProm.clean, get(obj[I]))
##     signalValues <- get(obj[I])$signalValue
##     cry2PeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
## }

## cry2PeaksTime <- table(apply(cry2PeakInt, 1, which.max))


## par(mfrow=c(2,1))
## barplot(allPeaksTime ~ colnames(cry2AllPeakInt), main='Phase all Cry2 peaks')
## barplot(cry2PeaksTime ~ colnames(cry2AllPeakInt), main='Phase Cry2 peaks near Up regulated genes')



######################################################################
## Do the same with BMAL1:
bmal1.p00 <- import("../publicData/takahashi12/macsPeaks/BMAL1/BMAL1_CT0_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
bmal1.p04 <- import("../publicData/takahashi12/macsPeaks/BMAL1/BMAL1_CT4_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
bmal1.p08 <- import("../publicData/takahashi12/macsPeaks/BMAL1/BMAL1_CT8_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
bmal1.p12 <- import("../publicData/takahashi12/macsPeaks/BMAL1/BMAL1_CT12_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
bmal1.p16 <- import("../publicData/takahashi12/macsPeaks/BMAL1/BMAL1_CT16_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
bmal1.p20 <- import("../publicData/takahashi12/macsPeaks/BMAL1/BMAL1_CT20_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')

## get all Bmal1 peaks (concatenate all peaks and collaps overlapping
## one)
all.bmal1.peaks <- reduce(c(bmal1.p00, bmal1.p04, bmal1.p08, bmal1.p12, bmal1.p16, bmal1.p20))

bmal1PeaksNumber <- c(length(bmal1.p00), length(bmal1.p04),
                      length(bmal1.p08), length(bmal1.p12),
                      length(bmal1.p16), length(bmal1.p20))

## find all peaks that are near bmal1 up-regulated genes
bmal1.peaks.inUpProm <- all.bmal1.peaks[queryHits(findOverlaps(all.bmal1.peaks, cry2.zt20.up.prom)),]

## get their signalValues at each time-point
bmal1PeakInt <- matrix(0, ncol=6, nrow=length(bmal1.peaks.inUpProm))
colnames(bmal1PeakInt) <- c('00','04','08','12','16','20')
obj <- c('bmal1.p00','bmal1.p04','bmal1.p08','bmal1.p12','bmal1.p16','bmal1.p20')
for (I in 1:length(obj)){
    overlap <- findOverlaps(bmal1.peaks.inUpProm, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    bmal1PeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

bmal1PeaksTime <- table(apply(bmal1PeakInt, 1, which.max))

## do the same for all peaks
bmal1.peaks.allProm <- all.bmal1.peaks[queryHits(findOverlaps(all.bmal1.peaks, all.prom)),]
bmal1AllPeakInt <- matrix(0, ncol=6, nrow=length(bmal1.peaks.allProm))
colnames(bmal1AllPeakInt) <- c('00','04','08','12','16','20')
obj <- c('bmal1.p00','bmal1.p04','bmal1.p08','bmal1.p12','bmal1.p16','bmal1.p20')
for (I in 1:length(obj)){
    overlap <- findOverlaps(bmal1.peaks.allProm, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    bmal1AllPeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

bmal1allPeaksTime <- table(apply(bmal1AllPeakInt, 1, which.max))


pdf('figures/2019-10-16/Bmal1peaksAll.pdf')
par(mfrow=c(1,1))
barplot(bmal1PeaksNumber ~ times, main='Bmal1 peaks')
dev.off()


pdf('figures/2019-10-16/Bmal1peaksAll_phase.pdf')
par(mfrow=c(1,1))
barplot(bmal1allPeaksTime ~ colnames(bmal1AllPeakInt), main='Phase all Bmal1 peaks')
dev.off()


pdf('figures/2019-10-16/Bmal1peaks_Upreg_phase.pdf')
par(mfrow=c(1,1))
barplot(c(bmal1PeaksTime,0,0) ~ colnames(bmal1AllPeakInt), main='Phase Bmal1 peaks near Up regulated genes')
dev.off()


pdf('figures/2019-10-16/Bmal1peaks_phase.pdf', height=10)
par(mfrow=c(2,1))
barplot(bmal1allPeaksTime ~ colnames(bmal1AllPeakInt), main='Phase all Bmal1 peaks')
barplot(c(bmal1PeaksTime,0,0) ~ colnames(bmal1AllPeakInt), main='Phase Bmal1 peaks near Up regulated genes')
dev.off()


## average signalValue of peaks
summary(bmal1AllPeakInt)
summary(bmal1PeakInt)

par(mfrow=c(1,2))
boxplot(log(bmal1AllPeakInt), ylim=c(3,7))
boxplot(log(bmal1PeakInt), ylim=c(3,7))


######################################################################
### Now the HNF4s
hnf4.p04 <- import("../publicData/kay18/macsPeaks/HNF4A_ZT4_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')
hnf4.p16 <- import("../publicData/kay18/macsPeaks/HNF4A_ZT16_peaks.narrowPeak",
                   format="BED", extraCols=extraCols_narrowPeak, genome='mm9')

## get all Hnf4 peaks (concatenate all peaks and collaps overlapping
## one)
all.hnf4.peaks <- reduce(c(hnf4.p04, hnf4.p16))

## find all peaks that are near hnf4 up-regulated genes
hnf4.peaks.inUpProm <- all.hnf4.peaks[queryHits(findOverlaps(all.hnf4.peaks, cry2.zt20.up.prom)),]

## get their signalValues at each time-point
hnf4PeakInt <- matrix(ncol=2, nrow=length(hnf4.peaks.inUpProm))
colnames(hnf4PeakInt) <- c('04','16')
obj <- c('hnf4.p04','hnf4.p16')
for (I in 1:length(obj)){
    overlap <- findOverlaps(hnf4.peaks.inUpProm, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    hnf4PeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

hnf4PeaksTime <- table(apply(hnf4PeakInt, 1, which.max))

## do the same for all peaks
hnf4AllPeakInt <- matrix(ncol=2, nrow=length(all.hnf4.peaks))
colnames(hnf4AllPeakInt) <- c('04','16')
obj <- c('hnf4.p04','hnf4.p16')
for (I in 1:length(obj)){
    overlap <- findOverlaps(all.hnf4.peaks, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    hnf4AllPeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

hnf4allPeaksTime <- table(apply(hnf4AllPeakInt, 1, which.max))


pdf('figures/2019-10-16/Hnf4peaksAll_phase.pdf')
barplot(hnf4allPeaksTime ~ colnames(hnf4AllPeakInt), main='Phase all Hnf4 peaks')
dev.off()

pdf('figures/2019-10-16/Hnf4peaks_Upreg_phase.pdf')
barplot(c(hnf4PeaksTime) ~ colnames(hnf4AllPeakInt), main='Phase Hnf4 peaks near Up regulated genes')
dev.off()

par(mfrow=c(2,1))
barplot(hnf4allPeaksTime ~ colnames(hnf4AllPeakInt), main='Phase all Hnf4 peaks')
barplot(c(hnf4PeaksTime) ~ colnames(hnf4AllPeakInt), main='Phase Hnf4 peaks near Up regulated genes')

## average signalValue of peaks
summary(hnf4AllPeakInt)
summary(hnf4PeakInt)

par(mfrow=c(1,2))
boxplot(log(hnf4AllPeakInt), ylim=c(3,7))
boxplot(log(hnf4PeakInt), ylim=c(3,7))



cry2Hnf4.peaks.inUpProm <- cry2.peaks.inUpProm[subjectHits(findOverlaps(hnf4.peaks.inUpProm, cry2.peaks.inUpProm))]
## of the 215 cry2 peaks near Up-promoters, 102 have a HNF4A peak that
## overlap. Their average peak width:
summary(width(ranges(cry2Hnf4.peaks.inUpProm)))


######################################################################
### Look for Cry2 peaks that have Bmap1 signal:

cry2Bmal1.peaks.inUpProm <- cry2.peaks.inUpProm[subjectHits(findOverlaps(bmal1.peaks.inUpProm, cry2.peaks.inUpProm))]
cry2Hnf4Bmal1.peaks.inUpProm <- cry2.peaks.inUpProm[subjectHits(findOverlaps(bmal1.peaks.inUpProm, cry2Hnf4.peaks.inUpProm))]

## How many genes are there:
length(unique(queryHits(findOverlaps(cry2.zt20.up.prom, cry2Hnf4.peaks.inUpProm))))
length(unique(queryHits(findOverlaps(cry2.zt20.up.prom, cry2Bmal1.peaks.inUpProm))))
length(unique(queryHits(findOverlaps(cry2.zt20.up.prom, cry2Hnf4Bmal1.peaks.inUpProm))))


cry2BmalPeakInt <- matrix(ncol=6, nrow=length(cry2Bmal1.peaks.inUpProm))
#colnames(hnf4AllPeakInt) <- c('04','16')
obj <- c('cry2.p00','cry2.p04','cry2.p08','cry2.p12','cry2.p16','cry2.p20')
for (I in 1:length(obj)){
    overlap <- findOverlaps(cry2Bmal1.peaks.inUpProm, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    cry2BmalPeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

cry2BmalPeakT <- table(apply(cry2BmalPeakInt, 1, which.max))



######################################################################
### Find Cry2 peaks that do not have Bmal1 signal

## 1. in promoters of up-regulated genes
cry2noBmal1.peaks.inUpProm <- cry2.peaks.inUpProm[-subjectHits(findOverlaps(bmal1.peaks.inUpProm, cry2.peaks.inUpProm))]

cry2noBmal1AllPeakInt <- matrix(ncol=6, nrow=length(cry2noBmal1.peaks.inUpProm))
#colnames(hnf4AllPeakInt) <- c('04','16')
obj <- c('cry2.p00','cry2.p04','cry2.p08','cry2.p12','cry2.p16','cry2.p20')
for (I in 1:length(obj)){
    overlap <- findOverlaps(cry2noBmal1.peaks.inUpProm, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    cry2noBmal1AllPeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

cry2noBmal1allPeaksTime1 <- table(apply(cry2noBmal1AllPeakInt, 1, which.max))

zt <- seq(0,20,4)

pdf('figures/2019-10-16/Cry2_phase_bmal_stratified.pdf', height=4, width=9)
par(mfrow=c(1,3))
barplot(cry2PeaksTime ~ zt, xlab="ZT", ylab="Peak numbers",
        main="Phase of Cry2 peaks in Up-regulated prom.")
barplot(cry2BmalPeakT ~ zt, xlab="ZT", ylab="Peak numbers",
        main="Phase of Cry2 peaks with BMAL1")
barplot(cry2noBmal1allPeaksTime1 ~ zt, xlab="ZT", ylab="Peak numbers",
        main="Phase of Cry2 peaks without BMAL1")
dev.off()

pdf('figures/2019-10-16/Cry2_phase_peksWithBmal.pdf', height=4)
barplot(cry2BmalPeakT ~ zt, xlab="ZT", ylab="Peak numbers",
        main="Phase of Cry2 peaks with BMAL1")
dev.off()

pdf('figures/2019-10-16/Cry2_phase_peksWithoutBmal.pdf', height=4)
barplot(cry2noBmal1allPeaksTime1 ~ zt, xlab="ZT", ylab="Peak numbers",
        main="Phase of Cry2 peaks without BMAL1")
dev.off()

## 2. in all promoters
cry2noBmal1.peaks <- all.cry2.peaks[-subjectHits(findOverlaps(all.bmal1.peaks, all.cry2.peaks))]

cry2noBmal1p <- matrix(ncol=6, nrow=length(cry2noBmal1.peaks))
#colnames(hnf4AllPeakInt) <- c('04','16')
obj <- c('cry2.p00','cry2.p04','cry2.p08','cry2.p12','cry2.p16','cry2.p20')
for (I in 1:length(obj)){
    overlap <- findOverlaps(cry2noBmal1.peaks, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    cry2noBmal1p[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}

cry2noBmal1peaksTime <- table(apply(cry2noBmal1p, 1, which.max))

par(mfrow=c(2,1))
barplot(cry2noBmal1peaksTime)
barplot(cry2noBmal1allPeaksTime1)



######################################################################
hnf4AllPeakInt <- matrix(ncol=6, nrow=length(cry2Hnf4.peaks.inUpProm))
#colnames(hnf4AllPeakInt) <- c('04','16')
obj <- c('cry2.p00','cry2.p04','cry2.p08','cry2.p12','cry2.p16','cry2.p20')
for (I in 1:length(obj)){
    overlap <- findOverlaps(cry2Hnf4.peaks.inUpProm, get(obj[I]))
    signalValues <- get(obj[I])$signalValue
    hnf4AllPeakInt[queryHits(overlap),I] <- signalValues[subjectHits(overlap)]
}


hnf4allPeaksTime1 <- table(apply(hnf4AllPeakInt, 1, which.max))

pdf('figures/2019-10-16/Cry2_phase_peksWithHnf4.pdf', height=4)
barplot(hnf4allPeaksTime1~zt, xlab="ZT", ylab="Peak numbers",
        main="Phase of Cry2 peaks with Hnf4a")
dev.off()
