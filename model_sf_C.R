## Run analysis, write model results

## Before: data_sf.RData and data_mf.RData (data)
## After:  config_sf.RData, results_sf.RData, 
##         config_mf.RData and results_mf.RData (model)

# .libPaths("C:/git/2025_her.27.3a47d_forecast/renv/library/windows/R-4.4/x86_64-w64-mingw32")

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

runName <- 'pil_pelnord2026_C_sf'

source('utilities_model_config.R')

retroFlag <- F

### ============================================================================
### Construct single fleet control object
### ============================================================================

load(file.path(dataPath,'data_pil.pelnord.C.Rdata'))

pg <- range(pil.C)[2]

pil.C@catch.n[,ac(2020)] <- NA
pil.C@landings.n[,ac(2020)] <- NA
pil.C@catch.wt[,ac(2020)] <- pil.C@catch.wt[,ac(2021)]
pil.C@discards.wt[,ac(2020)] <- pil.C@discards.wt[,ac(2021)]
pil.C@landings.wt[,ac(2020)] <- pil.C@landings.wt[,ac(2021)]
pil.C@stock.wt[,ac(2020)] <- pil.C@stock.wt[,ac(2021)]

#pil.C@
range(pil.C.tun$MA_ssb)['startf'] <- 0.5
range(pil.C.tun$MA_ssb)['endf'] <- 0.6

pil.C.tun$MA_ssb@index[,ac(2020)] <- NA
pil.C.tun$MA_ssb <- window(pil.C.tun$MA_ssb,end=2006)

pil.C.tun$MA_idx <- pil.C.tun$MA_idx[-c(1,2),]

pil.C.tun <- pil.C.tun[c("MA_idx","MA_ssb",'caton')]#,'caton'

# pil.C.tun$MA_ssb <- pil.C.tun$MA_ssb[,ac(1995:2000,2002:20224)]# for some ready 2001 is bolting the model
# pil.C.tun <- pil.C.tun[c("MA_idx","MA_ssb",'rec_russ')]#,"MA_ssb",'rec_russ','caton'

# pil.C.tun <- pil.C.tun[c("MA_ssb")]
# 
pil.ctrl.C <- config_sf_pil.C_cecaf2025(pil.C,
                                        pil.C.tun,
                                       pg)

pil.C@catch.n[ac(0),ac(c(1990:1992,2010))] <- NA
pil.C@catch.n[ac(5),ac(2013)] <- NA
pil.C@catch.n[ac(2),ac(2007:2008)] <- NA
# pil.C.tun$MA_ssb@index[,ac(2000:2024)] <- NA
# pil.C.tun$MA_ssb <- pil.C.tun$MA_ssb[,ac(1995:2000,2002:20224)]

#pil.ctrl@residuals <- T

### ============================================================================
### Run single fleet model then derive ref points
### ============================================================================

pil.samfit.C     <- FLSAM(pil.C,
                          pil.C.tun,
                        pil.ctrl.C,
                        return.fit=T)# 

# AIC(pil.samfit)

# catchplot(pil.samfit, addCI=TRUE)

#retro.samfit <- stockassessment::retro(pil.samfit,2)

# Run model
pil.ctrl.C@residuals <- T
pil.sam.C     <- FLSAM(pil.C,
                       pil.C.tun[c("MA_idx",'rec_russ')],
                     pil.ctrl.C)# return.fit=T

#AIC(pil.sam)

pil.C@stock.n[] <- pil.sam.C@stock.n[,ac(range(pil.C)["minyear"]:range(pil.C)["maxyear"])]
pil.C@harvest[] <- pil.sam.C@harvest[,ac(range(pil.C)["minyear"]:range(pil.C)["maxyear"])]

range(pil.C)['maxfbar'] <- 5
range(pil.sam.C)['maxfbar'] <- 5

# dat <- subset(residuals(pil.sam),fleet=="catch unique")
# print(ggplot(dat,aes(x=year,y=age,size=10*Cs(std.res),col=std.res>0))+
#         geom_point(shape=21)+
#         theme_bw()+
#         theme(legend.position = "none")+
#         scale_color_manual(values = c("red", "blue")))

### reference points

# fitting ricker and segReg
pil.spr0 <- yearMeans(spr0y(pil.C))

pil.sr.rick.C <- srrTMB(as.FLSR(pil.C,
                              model=ricker),
                      spr0=pil.spr0,plim=0.1)

pil.sr.segreg.C <- srrTMB(as.FLSR(pil.C,
                                model=segreg),
                        spr0=pil.spr0,plim=0.1)

# use a segReg fit to infer Blim
brp.segReg <- computeFbrp(pil.C, pil.sr.segreg.C)

# use a Ricker fit to infer MSY ref points
brp.rick <- computeFbrp(pil.C, pil.sr.rick.C)

refpts.C <- FLPars(pil=remap(refpts(brp.rick)))
refpts.C[[1]]$SBlim <- FLPar(c(refpts(brp.segReg)['Blim','ssb']), iter=1)
refpts.C[[1]]$SB01 <- FLPar(c(refpts(brp.rick)['Blim','ssb']), iter=1)
refpts.C <- refpts.C[[1]]

save(pil.C,
     pil.C.tun, 
     pil.ctrl.C,
     pil.sam.C,
     pil.samfit.C,
     refpts.C,
     pil.sr.rick.C,
     file=file.path(modelPath,paste0(runName,'.RData')))
