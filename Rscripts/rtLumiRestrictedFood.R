## This analysis is very similar to the one performe on 2020-11-05
## with the exception that data has been masked when light goes on

## install.packages(c("behavr","ggetho","damr","scopr","sleepr","zeitgebr"))

library("zeitgebr")
library('damr')
library('behavr')
library('ggetho')
library('scopr')
library('sleepr')
library('forcats')
library('dplyr')
library('preprocessCore')
library('imputeTS')
library('zoo')
library(cowplot)

######################################################################
### Some functions:
source('rtLumiFunctions.R')

######################################################################
### Read the data (start with the activity that showed more
### differences in Georgia anlysis)

metadata <- fread("../data/rtLumi/2022-04-25/2022-04-25_metadata",
                  colClasses=c('numeric', 'character', 'character',
                               'factor', 'numeric', 'factor',
                               'factor', 'factor'))
summary(metadata)
levels(metadata[,geno])
metadata$id <- as.factor(metadata$id)
setkey(metadata, id)

## Get data from lumicycler:
dt.1 <- getWheelData('../data/rtLumi/2022-04-25/SB1_8_mouse2_WT.dat',
                     '1', type='zt')
dt.2 <- getWheelData('../data/rtLumi/2022-04-25/SB1_8_mouse3_WT.dat',
                     '2', type='zt')
dt.3 <- getWheelData('../data/rtLumi/2022-04-25/SB1_8_mouse4_MUT.dat',
                     '3', type='zt')
dt.4 <- getWheelData('../data/rtLumi/2022-04-25/SB2_8_mouse1_WT.dat',
                     '4', type='zt')
dt.5 <- getWheelData('../data/rtLumi/2022-04-25/SB2_8_mouse2_MUT.dat',
                     '5', type='zt')
dt.6 <- getWheelData('../data/rtLumi/2022-04-25/SB2_8_mouse3_MUT.dat',
                     '6', type='zt')
dt.7 <- getWheelData('../data/rtLumi/2022-04-25/SB2_8_mouse4_WT.dat',
                     '7', type='zt')
dt.8 <- getWheelData('../data/rtLumi/2022-04-25/SB3_mouse3_mUT.dat',
                     '8', type='zt')
dt.9 <- getWheelData('../data/rtLumi/2022-04-25/SB3_mouse4_WT.dat',
                     '9', type='zt')

### Exclude first days of experiment and put food shifting starting at
### the same day:
dt1 <- dt.1[Day >= 6]
dt2 <- dt.2[Day >= 6]
dt3 <- dt.3[Day >= 6]
dt4 <- dt.4[Day >= 2]
dt5 <- dt.5[Day >= 2]
dt6 <- dt.6[Day >= 2]
dt7 <- dt.7[Day >= 2]
dt8 <- dt.8[Day >= 5]
dt9 <- dt.9[Day >= 5]

dt1$Day <- dt1$Day - 5
dt2$Day <- dt2$Day - 5
dt3$Day <- dt3$Day - 5
dt4$Day <- dt4$Day - 2
dt5$Day <- dt5$Day - 2
dt6$Day <- dt6$Day - 2
dt7$Day <- dt7$Day - 2
dt8$Day <- dt8$Day - 4
dt9$Day <- dt9$Day - 4

dt1$t <- dt1$t - days(4)
dt2$t <- dt2$t - days(4)
dt3$t <- dt3$t - days(4)
#dt4$t <- dt4$t - days(1)
#dt5$t <- dt5$t - days(1)
#dt6$t <- dt6$t - days(1)
#dt7$t <- dt7$t - days(1)
dt8$t <- dt8$t - days(3)
dt9$t <- dt9$t - days(3)


scalingFactor <- c(median(dt1$Photon), median(dt2$Photon),
                   median(dt3$Photon), median(dt4$Photon),
                   median(dt5$Photon), median(dt6$Photon),
                   median(dt7$Photon), median(dt8$Photon),
                   median(dt9$Photon))
scalingFactor <- scalingFactor/min(scalingFactor)
names(scalingFactor) <- 1:9

## dt.all <- rbind(dt1, dt2, dt3, dt5, dt6, dt7, dt8, dt9)
dt.all <- rbind(dt2, dt3, dt5, dt6, dt7, dt8, dt9)
setkey(dt.all, id)

dt <- behavr(dt.all, metadata)
## dt$Photon <- dt$Photon / scalingFactor[dt$id]
dt$Photon[which(dt$Light == 1)] <- NA
dt$Photon[which(dt$Photon <= 200)] <- NA
dt


## Normalize data
photons <- data.frame(dt2 = dt[t >= days(1) & t <= days(8) & id == '2']$Photon,
                      dt3 = dt[t >= days(1) & t <= days(8) & id == '3']$Photon,
                      ## dt4 = dt[t >= days(1) & t <= days(8) & id == '4']$Photon,
                      dt5 = dt[t >= days(1) & t <= days(8) & id == '5']$Photon,
                      dt6 = dt[t >= days(1) & t <= days(8) & id == '6']$Photon,
                      dt7 = dt[t >= days(1) & t <= days(8) & id == '7']$Photon,
                      dt8 = dt[t >= days(1) & t <= days(8) & id == '8']$Photon,
                      dt9 = dt[t >= days(1) & t <= days(8) & id == '9']$Photon)
summary(photons)

photonsNorm <- round(normalize.quantiles(as.matrix(photons)))

## put normalise data in dt:
dt[t >= days(1) & t <= days(8) & id == '2']$Photon <- round(photonsNorm[,1])
dt[t >= days(1) & t <= days(8) & id == '3']$Photon <- round(photonsNorm[,2])
dt[t >= days(1) & t <= days(8) & id == '5']$Photon <- round(photonsNorm[,3])
dt[t >= days(1) & t <= days(8) & id == '6']$Photon <- round(photonsNorm[,4])
dt[t >= days(1) & t <= days(8) & id == '7']$Photon <- round(photonsNorm[,5])
dt[t >= days(1) & t <= days(8) & id == '8']$Photon <- round(photonsNorm[,6])
dt[t >= days(1) & t <= days(8) & id == '9']$Photon <- round(photonsNorm[,7])

## Make rolling mean:
dt.norm <- dt
dt.norm[t >= days(1) & t <= days(8) & id == '2']$Photon <- round(rollmean(na_kalman(photonsNorm[,1]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '3']$Photon <- round(rollmean(na_kalman(photonsNorm[,2]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '5']$Photon <- round(rollmean(na_kalman(photonsNorm[,3]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '6']$Photon <- round(rollmean(na_kalman(photonsNorm[,4]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '7']$Photon <- round(rollmean(na_kalman(photonsNorm[,5]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '8']$Photon <- round(rollmean(na_kalman(photonsNorm[,6]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '9']$Photon <- round(rollmean(na_kalman(photonsNorm[,7]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))


### Normalise activity
activity <- data.frame(dt2 = dt[t >= days(1) & t <= days(8) & id == '2']$Activity,
                      dt3 = dt[t >= days(1) & t <= days(8) & id == '3']$Activity,
                      ## dt4 = dt[t >= days(1) & t <= days(8) & id == '4']$Activity,
                      dt5 = dt[t >= days(1) & t <= days(8) & id == '5']$Activity,
                      dt6 = dt[t >= days(1) & t <= days(8) & id == '6']$Activity,
                      dt7 = dt[t >= days(1) & t <= days(8) & id == '7']$Activity,
                      dt8 = dt[t >= days(1) & t <= days(8) & id == '8']$Activity,
                      dt9 = dt[t >= days(1) & t <= days(8) & id == '9']$Activity)
summary(activity)

activityNorm <- round(normalize.quantiles(as.matrix(activity)))

## put normalise data in dt:
dt[t >= days(1) & t <= days(8) & id == '2']$Activity <- round(activityNorm[,1])
dt[t >= days(1) & t <= days(8) & id == '3']$Activity <- round(activityNorm[,2])
dt[t >= days(1) & t <= days(8) & id == '5']$Activity <- round(activityNorm[,3])
dt[t >= days(1) & t <= days(8) & id == '6']$Activity <- round(activityNorm[,4])
dt[t >= days(1) & t <= days(8) & id == '7']$Activity <- round(activityNorm[,5])
dt[t >= days(1) & t <= days(8) & id == '8']$Activity <- round(activityNorm[,6])
dt[t >= days(1) & t <= days(8) & id == '9']$Activity <- round(activityNorm[,7])

## Make rolling mean:
dt.norm <- dt.norm
dt.norm[t >= days(1) & t <= days(8) & id == '2']$Activity <- round(rollmean(na_kalman(activityNorm[,1]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '3']$Activity <- round(rollmean(na_kalman(activityNorm[,2]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '5']$Activity <- round(rollmean(na_kalman(activityNorm[,3]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '6']$Activity <- round(rollmean(na_kalman(activityNorm[,4]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '7']$Activity <- round(rollmean(na_kalman(activityNorm[,5]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '8']$Activity <- round(rollmean(na_kalman(activityNorm[,6]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))
dt.norm[t >= days(1) & t <= days(8) & id == '9']$Activity <- round(rollmean(na_kalman(activityNorm[,7]),
                                                                     k=6*60,
                                                                     fill = c("NA", "extend", "NA")))


######################################################################
## plot trends 
ggetho(dt[t >= days(1) & t <= days(8)], aes(x=t, y=id, z=Food)) +
    stat_tile_etho() +
    geom_vline(xintercept = days(3))

p1 <- ggetho(dt[t >= days(1) & t <= days(8)],
       aes(x=t, y=Photon, colour = geno)) +
    geom_smooth(methods = loess, span = 0.1) +
    ## stat_pop_etho() +
    facet_grid( sample ~ .) +
    geom_vline(xintercept = days(3)) +
    theme_light(14)
ggsave("figures/2022-04-26/samplesPhotonSignal.pdf", p1, height = 10,
       width = 7)


p2 <- ggetho(dt[t >= days(1) & t <= days(8)],
       aes(x=t, y=Photon, colour=geno)) +
    stat_ld_annotations(height=1, alpha=0.1, outline = NA) +
    ## geom_smooth(methods = loess, span = 0.01) +
    stat_pop_etho() +
    geom_vline(xintercept = days(3)) +
    theme_light(14)
ggsave("figures/2022-04-26/genotypePhotonSignal.pdf", p2, height = 7,
       width = 10)


p3 <- ggetho(dt[t >= days(1) & t <= days(8)],
       aes(x=t, y=Activity, colour=geno)) +
    stat_ld_annotations(height=1, alpha=0.1, outline = NA) +
    stat_pop_etho() +
    geom_vline(xintercept = days(3)) +
    theme_light()



p4 <- ggetho(dt[t >= days(3) & t <= days(8)],
       aes(x=t, y=Photon, colour=geno), time_wrap = hours(24)) +
    stat_pop_etho() +
    labs(title = "Daily signal after food shift") +
    theme_gray(14)
ggsave("figures/2022-04-26/dailyPhotonSignal.pdf", p4, height = 7,
       width = 7)


ggetho(dt[t >= days(5) & t <= days(8)],
       aes(x=t, y=Photon, colour=geno), time_wrap = hours(24)) +
    stat_pop_etho() +
    labs(title = "Daily signal after food shift") +
    theme_gray(14)


######################################################################
### Now with the rolling mean data

p1 <- ggetho(dt.norm[t >= days(1) & t <= days(8)],
       aes(x=t, y=Photon, colour = geno)) +
    geom_smooth(methods = loess, span = 0.1) +
    ## stat_pop_etho() +
    facet_grid( sample ~ .) +
    geom_vline(xintercept = days(3)) +
    theme_light(14)
ggsave("figures/2022-04-26/samplesPhotonSignalRM.pdf", p1, height = 10,
       width = 7)

p1 <- ggetho(dt.norm[t >= days(1) & t <= days(7)],
       aes(x=t, y=Photon, colour = geno, group = sample)) +
    geom_smooth(methods = loess, span = 0.1) +
    ## stat_pop_etho() +
    facet_grid( geno ~ .) +
    geom_vline(xintercept = days(3)) +
    theme_gray(14) +
    theme(legend.position = "none")
ggsave("figures/2022-04-26/samplesPhotonSignalRMgrouped.pdf", p1, height = 6,
       width = 7)

p1.1 <- ggetho(dt.norm[t >= days(1) & t <= days(7)],
       aes(x=t, y=Activity, colour = geno, group = sample)) +
    geom_smooth(methods = loess, span = 0.1, se = FALSE) +
    ## stat_pop_etho() +
    facet_grid( geno ~ .) +
    geom_vline(xintercept = days(3)) +
    theme_gray(14) +
    theme(legend.position = "none")
ggsave("figures/2022-04-26/samplesActivitySignalRMgrouped.pdf", p1.1,
       height = 6, width = 7)


p2 <- ggetho(dt.norm[t >= days(1) & t <= days(8)],
       aes(x=t, y=Photon, colour=geno)) +
    stat_ld_annotations(height=1, alpha=0.1, outline = NA) +
    ## geom_smooth(methods = loess, span = 0.01) +
    stat_pop_etho() +
    geom_vline(xintercept = days(3)) +
    theme_light(14)
ggsave("figures/2022-04-26/genotypePhotonSignalRM.pdf", p2, height = 7,
       width = 10)


p3 <- ggetho(dt.norm[t >= days(1) & t <= days(8)],
       aes(x=t, y=Activity, colour=geno)) +
    stat_ld_annotations(height=1, alpha=0.1, outline = NA) +
    stat_pop_etho() +
    geom_vline(xintercept = days(3)) +
    theme_light()



p4 <- ggetho(dt.norm[t >= days(3) & t <= days(7)],
       aes(x=t, y=Photon, colour=geno), time_wrap = hours(24)) +
    stat_pop_etho() +
    labs(title = "Daily signal after food shift") +
    theme_gray(14)
ggsave("figures/2022-04-26/dailyPhotonSignalRM.pdf", p4, height = 7,
       width = 7)


p5 <- ggetho(dt.norm[t >= days(5) & t <= days(7)],
       aes(x=t, y=Photon, colour=geno), time_wrap = hours(24)) +
    stat_pop_etho() +
    labs(title = "Daily signal last 2 days") +
    scale_x_continuous(breaks = seq(0, hours(24), by = hours(1))) +
    theme_gray(14)
ggsave("figures/2022-04-26/dailyPhotonSignalLateRM.pdf", p5, height = 7,
       width = 7)

p5.1 <- ggetho(dt.norm[t >= days(5) & t <= days(7)],
       aes(x=t, y=Activity, colour=geno), time_wrap = hours(24)) +
    stat_pop_etho() +
    labs(title = "Daily signal last 2 days") +
    scale_x_continuous(breaks = seq(0, hours(24), by = hours(1))) +
    theme_gray(14)
ggsave("figures/2022-04-26/dailyActivitySignalLateRM.pdf", p5.1, height = 7,
       width = 7)


dt.norm.g <- dt.norm[t >= days(5) & t <= days(7)]
dt.norm.g$geno <- metadata[dt.norm.g$id, 'geno']
dt.norm.g$t <- (dt.norm.g$t / 3600) / 24
w1 <- which(dt.norm.g$t >= 5 & dt.norm.g$t < 6)
w2 <- which(dt.norm.g$t >= 6 & dt.norm.g$t < 7)
w3 <- which(dt.norm.g$t >= 7 & dt.norm.g$t < 8)
dt.norm.g$t[w1] <- dt.norm.g$t[w1] - 5
dt.norm.g$t[w2] <- dt.norm.g$t[w2] - 6
dt.norm.g$t[w3] <- dt.norm.g$t[w3] - 7
summary(dt.norm.g$t)

per <- 24

sp <- 0.05
reslo <- loess(Photon ~ t, data = dt.norm.g[geno == "WT"], span = sp)
plot(predict(reslo)~ dt.norm.g[geno == "WT"]$t, type = "l")
reslom <- loess(Photon ~ t, data = dt.norm.g[geno == "MUT"], span = sp)
points(predict(reslom)~dt.norm.g[geno == "MUT"]$t, col = "red")

(24 * dt.norm.g[geno == "WT"]$t[which.min(predict(reslo))]) -
(24 * dt.norm.g[geno == "WT"]$t[which.min(predict(reslom))])

(24 * dt.norm.g[geno == "WT"]$t[which.max(predict(reslo))]) -
(24 * dt.norm.g[geno == "WT"]$t[which.max(predict(reslom))])

reslm <- lm(Photon ~ sin(2*pi/per*t) : geno +
                cos(2*pi/per*t) : geno +
                sin(4*pi/per*t) : geno +
                cos(4*pi/per*t) : geno + geno, data = dt.norm.g)
summary(reslm)

pMut <- predict(reslm, newdata = data.frame(t = unique(dt.norm.g$t), geno = "MUT"))
pWt <- predict(reslm, newdata = data.frame(t = unique(dt.norm.g$t), geno = "WT"))

plot(pWt ~ unique(dt.norm.g$t), col=2, type = "l")
lines(pMut ~ unique(dt.norm.g$t), col=3, lty=2)

unique(dt.norm.g$t)[which.min(pWt[1:1200])]
unique(dt.norm.g$t)[which.min(pMut[1:1200])]

unique(dt.norm.g$t)[which.max(pWt[1:1200])]
unique(dt.norm.g$t)[which.max(pMut[1:1200])]


reslm <- lm(Activity ~ sin(2*pi/per*t) : geno +
                cos(2*pi/per*t) : geno +
                sin(4*pi/per*t) : geno +
                cos(4*pi/per*t) : geno + geno, data = dt.norm.g)
summary(reslm)

pMut <- predict(reslm, newdata = data.frame(t = unique(dt.norm.g$t), geno = "MUT"))
pWt <- predict(reslm, newdata = data.frame(t = unique(dt.norm.g$t), geno = "WT"))

plot(pWt ~ unique(dt.norm.g$t), col=2, type = "l")
lines(pMut ~ unique(dt.norm.g$t), col=3, lty=2)

unique(dt.norm.g$t)[which.min(pWt[500:1200])+500]
unique(dt.norm.g$t)[which.min(pMut[500:1200])+500]

unique(dt.norm.g$t)[which.max(pWt[500:1200])+500]
unique(dt.norm.g$t)[which.max(pMut[500:1200])+500]



######################################################################
## I de-trended data, but it should not make sense here
colnames(dt1)[4] <- colnames(dt2)[4] <- colnames(dt3)[4] <- colnames(dt4)[4] <- colnames(dt5)[4] <- colnames(dt6)[4] <- colnames(dt7)[4] <- colnames(dt8)[4] <- colnames(dt9)[4] <- "Cnt"
dt.1.norm <- detrendData(dt1)
dt.2.norm <- detrendData(dt2)
dt.3.norm <- detrendData(dt3)
dt.4.norm <- detrendData(dt4)
dt.5.norm <- detrendData(dt5)
dt.6.norm <- detrendData(dt6)
dt.7.norm <- detrendData(dt7)
dt.8.norm <- detrendData(dt8)
dt.9.norm <- detrendData(dt9)

## Detranded data:
dt.all.norm <- rbind(dt.1.norm, dt.2.norm, dt.3.norm, dt.4.norm,
                     dt.5.norm, dt.6.norm, dt.7.norm, dt.8.norm,
                     dt.9.norm)
setkey(dt.all.norm, id)

dt.norm <- behavr(dt.all.norm, metadata)
dt.norm$Cnt[which(dt.norm$Light == 1)] <- NA
dt.norm$Cnt[which(dt.norm$Cnt <= 2000)] <- NA


ggetho(dt.norm[t >= days(1) & t <= days(8)],
       aes(x=t, y=Cnt, colour=geno)) +
    stat_pop_etho()

ggetho(dt.norm[t >= days(3) & t <= days(8)],
       aes(x=t, y=Cnt, colour=geno), time_wrap = hours(24)) +
    stat_pop_etho()


