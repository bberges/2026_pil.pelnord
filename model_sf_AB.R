## Run analysis, write model results

## Before: data_sf.RData and data_mf.RData (data)
## After:  config_sf.RData, results_sf.RData, 
##         config_mf.RData and results_mf.RData (model)

rm(list=ls())

library(icesTAF)

library(FLCore)
library(stockassessment)
library(FLSAM)
library(FLRef)

# taf.session()

library(methods)


mkdir("model")

dataPath               <- file.path(".", "data")
modelPath              <- file.path(".",'model')

runName <- 'pil_pelnord2026_AB_sf'

source('utilities_model_config.R')

retroFlag <- F

### ============================================================================
### Construct single fleet control object
### ============================================================================

load(file.path(dataPath,'data_pil.pelnord.AB.Rdata'))

pg <- range(pil.AB)[2]

pil.AB.tun$MA_ssb@index[,ac(2020)] <- NA
pil.AB.tun$MA_ssb <- window(pil.AB.tun$MA_ssb,end=2006)

# pil.AB.tun$MA_idx <- pil.AB.tun$MA_idx[-c(1,2),]
# pil.AB.tun <- pil.AB.tun[c("MA_idx","rec_russ")]#,'caton'


# stk <- pil.AB
# tun <- pil.AB.tun

pil.AB.tun <- pil.AB.tun[c("MA_idx","MA_ssb")]#,'caton',"rec_russ",

pil.ctrl.AB <- config_sf_pil.AB_cecaf2025(pil.AB,
                                          pil.AB.tun,
                                            pg)
#pil.ctrl@residuals <- T

### ============================================================================
### Run single fleet model then derive ref points
### ============================================================================

pil.samfit.AB     <- FLSAM(pil.AB,
                        pil.AB.tun,
                        pil.ctrl.AB,
                        return.fit=T)
# windows()
# plot(pil.samfit.AB)
# AIC(pil.samfit)

# catchplot(pil.samfit, addCI=TRUE)

#retro.samfit <- stockassessment::retro(pil.samfit,2)

# Run model
pil.ctrl.AB@residuals <- T
pil.sam.AB     <- FLSAM(pil.AB,
                     pil.AB.tun,
                     pil.ctrl.AB)# return.fit=T

#AIC(pil.sam)

pil.AB@stock.n[] <- pil.sam.AB@stock.n[,ac(range(pil.AB)["minyear"]:range(pil.AB)["maxyear"])]
pil.AB@harvest[] <- pil.sam.AB@harvest[,ac(range(pil.AB)["minyear"]:range(pil.AB)["maxyear"])]

range(pil.AB)['maxfbar'] <- 5
range(pil.sam.AB)['maxfbar'] <- 5

# dat <- subset(residuals(pil.sam.AB),fleet=="catch AB")
# print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
#         geom_point(shape=21)+
#         theme_bw()+
#         theme(legend.position = "none")+
#         scale_color_manual(values = c("red", "blue")))

### reference points

# fitting ricker and segReg
pil.spr0 <- yearMeans(spr0y(pil.AB))

pil.sr.rick.AB <- srrTMB(as.FLSR(pil.AB,
                              model=ricker),
                      spr0=pil.spr0,plim=0.1)

pil.sr.segreg.AB <- srrTMB(as.FLSR(pil.AB,
                                model=segreg),
                        spr0=pil.spr0,plim=0.1)

# use a segReg fit to infer Blim
brp.segReg <- computeFbrp(pil.AB, pil.sr.segreg.AB)

# use a Ricker fit to infer MSY ref points
brp.rick <- computeFbrp(pil.AB, pil.sr.rick.AB)

refpts.AB <- FLPars(pil=remap(refpts(brp.rick)))
refpts.AB[[1]]$SBlim <- FLPar(c(refpts(brp.segReg)['Blim','ssb']), iter=1)
refpts.AB[[1]]$SB01 <- FLPar(c(refpts(brp.rick)['Blim','ssb']), iter=1)
refpts.AB <- refpts.AB[[1]]

save(pil.AB,
     pil.AB.tun, 
     pil.ctrl.AB,
     pil.sam.AB,
     pil.samfit.AB,
     refpts.AB,
     pil.sr.rick.AB,
     file=file.path(modelPath,paste0(runName,'.RData')))
