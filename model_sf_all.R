## Run analysis, write model results

## Before: data_sf.RData and data_mf.RData (data)
## After:  config_sf.RData, results_sf.RData, 
##         config_mf.RData and results_mf.RData (model)

# setwd('C:/git/2026_pil.pelnord')

# .libPaths( c( .libPaths()[2], "C:/git/2025_her.27.3a47d_forecast/renv/library/windows/R-4.4/x86_64-w64-mingw32") )  

rm(list=ls())

library(icesTAF)

library(FLCore)
library(stockassessment)
library(FLSAM)
library(FLBRP)
library(FLSRTMB)
library(FLRef)

# taf.session()

library(methods)


mkdir("model")

dataPath               <- file.path(".", "data")
modelPath              <- file.path(".",'model')

runName <- 'pil_pelnord2026_sf'

source('utilities_model_config.R')

retroFlag <- T

### ============================================================================
### Construct single fleet control object
### ============================================================================

load(file.path(dataPath,'data_pil.pelnord.all.Rdata'))

pg <- range(pil.all)[2]

pil.all@catch.n[,ac(2020)] <- NA
pil.all@landings.n[,ac(2020)] <- NA
pil.all.tun$MA_ssb@index[,ac(2020)] <- NA
pil.all.tun$MA_ssb <- window(pil.all.tun$MA_ssb,end=2006)

# plot(quantSums(pil.all.tun$MA_idx@index*pil.all@stock.wt[,ac(2007:2025)]*pil.all@mat[,ac(2007:2025)]))
# 
# plot(pil.all.tun$MA_ssb[,ac(2007:2025)])

# pil.all.tun$MA_idx <- pil.all.tun$MA_idx[-c(1,2),]

# stk <- pil.all
# tun <- pil.all.tun
# pil.all.tun <- tun

#pil.all@catch[] <- computeCatch(pil.all)

pil.all.tun <- pil.all.tun[c("MA_idx","MA_ssb","nansen")]#,'caton',"rec_russ",

pil.ctrl <- config_sf_pil.all_cecaf2025(pil.all,
                                        pil.all.tun,
                                        pg)

pil.all@catch.n[ac(0),ac(c(1990:1992,2010))] <- NA
pil.all@catch.n[ac(5),ac(2013)] <- NA
pil.all@catch.n[ac(2),ac(2007:2008)] <- NA

# temp <- new('FLSAMs')
# temp[['caton']] <- pil.sam
# temp[['base']] <- pil.sam
# 
# plot(temp)

# cohcorrplot(pil.all.tun[["MA_idx"]]@index)
# 
# cohcorrplot(pil.all@catch.n)


### ============================================================================
### Run single fleet model then derive ref points
### ============================================================================

pil.samfit     <- FLSAM(pil.all,
                     pil.all.tun,
                     pil.ctrl,
                     return.fit=T)#

# Run model
pil.ctrl@residuals <- T
pil.sam     <- FLSAM(pil.all,
                     pil.all.tun,
                     pil.ctrl)# return.fit=T

#AIC(pil.sam)

pil.all@stock.n[] <- pil.sam@stock.n[,ac(range(pil.all)["minyear"]:range(pil.all)["maxyear"])]
pil.all@harvest[] <- pil.sam@harvest[,ac(range(pil.all)["minyear"]:range(pil.all)["maxyear"])]

range(pil.all)['maxfbar'] <- 5
range(pil.sam)['maxfbar'] <- 5

# dat <- subset(residuals(pil.sam),fleet=="catch unique")
# print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
#         geom_point(shape=21)+
#         theme_bw()+
#         theme(legend.position = "none")+
#         scale_color_manual(values = c("red", "blue")))

### reference points

# fitting ricker and segReg
pil.spr0 <- yearMeans(spr0y(pil.all))

pil.sr.rick <- srrTMB(as.FLSR(pil.all,
                              model=ricker),
                      spr0=pil.spr0,plim=0.1)
pil.sr.segreg <- srrTMB(as.FLSR(pil.all,
                                model=segreg),
                        spr0=pil.spr0,plim=0.1)

# use a segReg fit to infer Blim
brp.segReg <- computeFbrp(pil.all, pil.sr.segreg)

# use a Ricker fit to infer MSY ref points
brp.rick <- computeFbrp(pil.all, pil.sr.rick,fmax = 10)

refpts <- FLPars(pil=remap(refpts(brp.rick)))
refpts[[1]]$SBlim <- FLPar(c(refpts(brp.segReg)['Blim','ssb']), iter=1)
refpts[[1]]$SB01 <- FLPar(c(refpts(brp.rick)['Blim','ssb']), iter=1)
refpts <- refpts[[1]]

save(pil.all,
     pil.all.tun, 
     pil.ctrl,
     pil.sam,
     pil.samfit,
     refpts,
     pil.sr.rick,
     file=file.path(modelPath,paste0(runName,'.RData')))

### ============================================================================
### run single fleet retro
### ============================================================================

if(retroFlag){
  # retro
  n.retro.years <- 5  # Number of years for which to run the retrospective
  pil.ctrl@residuals <- FALSE
  pil.retro <- retro(pil.all, pil.all.tun, pil.ctrl, retro=n.retro.years)
  
  save(pil.retro,
       file=file.path(modelPath,paste0(runName,'_retro.RData')))
}

### ============================================================================
### dump SAM fit
### ============================================================================

loo <- leaveout(pil.sam)

dataplot(pil.sam)

leaveout(pil.sam)

deterministicReferencepoints(pil.sam, c("0.1B0"))

catchplot(pil.sam, addCI=TRUE)
fitplot(pil.sam, fleets=4)
plot(pil.sam)
srplot(pil.sam)
plot(ypr(pil.sam))
fitplot(pil.sam, fleets=4)

setcap<-function(title="", caption=""){   
  tit.list[length(tit.list)+1]<<-paste("# Title",plotcounter,"\n")
  tit.list[length(tit.list)+1]<<-paste(title,"\n")
  tit.list[length(tit.list)+1]<<-paste("# Caption",plotcounter,"\n")
  tit.list[length(tit.list)+1]<<-paste(caption,"\n")
  plotcounter<<-plotcounter+1 
}

plotcounter<-1
tit.list<-list()

for(f in 1:pil.sam$data$noFleets){
  fitplot(pil.sam, fleets=f)
  setcap("Fit to data", "Predicted line and observed points (log scale)")
}

### ============================================================================
### dump stock recruitment
### ============================================================================

ploteq(computeFbrp(pil.all, pil.sr.rick))

refpts(pil.rp)


brp <- computeFbrp(pil.sr.rick, proxy=c("msy","f0.1"),
                   blim=0.1, type="btgt", fmax=4)

plot(res.rick)+theme(legend.position="right")

# fitting with eqSim
FIT_ricker <- eqsr_fit(pil.all,nsamp=2000, models = "Ricker")
FIT_segreg <- eqsr_fit(pil.all,nsamp=2000, models = "Segreg")

FIT_ricker$sr.det
FIT_segreg$sr.det

eqsr_plot(FIT_ricker, ggPlot=TRUE)+
  theme(legend.position="none")

eqsr_plot(FIT_segreg, ggPlot=TRUE)+
  theme(legend.position="none")
