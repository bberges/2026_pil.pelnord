
rm(list=ls())

library(icesTAF)
library(FLCore)
library(stockassessment)
library(FLSAM)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(scales)
library(kobe)
library(gridExtra)

mkdir("report")

dataSource     <- file.path(".", "bootstrap",'data')
dataPath       <- file.path(".", "data")
modelPath      <- file.path(".", "model")
reportPath     <- file.path(".", "report")

source('utilities_plots.R')

runName <- 'pil_pelnord2026_sf'

load(file.path(modelPath,paste0('pil_pelnord2026_AB_sf','.RData')))
load(file.path(modelPath,paste0('pil_pelnord2026_C_sf','.RData')))
load(file.path(modelPath,paste0(runName,'.RData')))
load(file.path(modelPath,paste0(runName,'_retro.RData')))

# load(file.path(dataPath,'data_pil.pelnord.all.Rdata'))
# load(file.path(dataPath,'data_pil.pelnord.AB.Rdata'))
# load(file.path(dataPath,'data_pil.pelnord.C.Rdata'))

#------------------------------------------
### get survey from age 0
pil.age.all <- read.csv(file.path(dataPath,'MA_idx.csv'))

ageRange <- ac(0:5)
yearRange <- an(unique(pil.age.all$year))
FLR_temp <- FLQuant(  array(0,dim = c(length(unique(pil.age.all$age)),
                                      length(min(yearRange):max(yearRange)),
                                      1,1,2,1)),
                      dimnames=list(age=ac(unique(pil.age.all$age)),
                                    year=ac(min(yearRange):max(yearRange)),
                                    unit=c("1e6"),
                                    season='all',
                                    area=c('AB','C')))

FLR_temp[,ac(unique(pil.age.all$year)),,,'AB']  <- 
  as.matrix(subset(pil.age.all,region == 'AB') %>% 
              select(-c(region)) %>%
              pivot_wider(names_from = year, values_from = idx) %>%
              select(-c(age)))

FLR_temp[,ac(unique(pil.age.all$year)),,,'C']   <- 
  as.matrix(subset(pil.age.all,region == 'C') %>% 
              select(-c(region)) %>%
              pivot_wider(names_from = year, values_from = idx) %>%
              select(-c(age)))

MA_idx <- FLIndex(name = 'MA_idx',
                  desc = 'Marocan survey index',
                  index=setPlusGroup(FLR_temp,5,na.rm=T)) # set plus group to 4 to avoid NA in the C area
range(MA_idx)[c('startf', 'endf')] <- c(0.85,0.9) # Nov-Dec
#MA_idx@index@.Data[which(is.na(MA_idx@index))] <- -1
MA_idx <- MA_idx[ageRange,] # start from age 1 as age 0 is flimsy
MA_idx@index[,ac(2020)] <- NA

# all regions combined
FLR_temp <- FLQuant(  array(0,dim = c(length(dimnames(MA_idx)$age),
                                      length(dimnames(MA_idx)$year),
                                      1,1,1,1)),
                      dimnames=list(age=dimnames(MA_idx)$age,
                                    year=dimnames(MA_idx)$year,
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))
MA_idx.all <- FLIndex(name = 'MA_idx',
                      desc = 'Marocan survey index',
                      index=FLR_temp)

MA_idx.all@index[] <- areaSums(MA_idx@index[,,,,c('AB','C')])
range(MA_idx.all)[c('startf', 'endf')] <- c(0.5,0.6)
#MA_idx.all@index@.Data[which(is.na(MA_idx.all@index))] <- -1
MA_idx.all <- MA_idx.all[ageRange,] # start from age 1 as age 0 is flimsy
MA_idx.all@index[,ac(2020)] <- NA

#------------------------------------------
### input data

prefix <- 'input'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

# natural mortality
# taf.png(paste0(prefix,"_nat_mort"))
df.m <- as.data.frame(yearMeans(pil.all@m))

p <- ggplot(df.m,aes(x=age,y=data))+
      geom_line()+
      theme_bw()+
      ylim(0,NA)+
      ylab('Total natural mortality')

ggsave(paste0(prefix,"_nat_mort.png"),
       p,
       width = 80,
       height = 80,
       units = c("mm"),
       dpi = 300)
# dev.off()

## catch figure 
# taf.png(paste0(prefix,"_prop_catches"))
prop_catches <- as.data.frame(pil.all@catch.n)
prop_catches$data[prop_catches$data == -1] <- NA

p1 <- ggplot(subset(prop_catches,year>=2000),aes(x=year,y=data))+
      theme_bw()+
      geom_bar(aes(fill = as.factor(age)),stat="identity")+#,position = "fill"
      theme_bw()+
      labs(x = 'Year',
           y = 'Proportion of Catch numbers at age',
           fill = "Age")+
      scale_y_continuous(expand = c(0, 0))

p2 <- cohcorrplot(pil.all@catch.n)

p <- grid.arrange(p1,p2)#,ncol=2

ggsave(paste0(prefix,"_catch_fig.png"),
       p,
       width = 170,
       height = 200,
       units = c("mm"),
       dpi = 300)

## index at age figure
prop_idx <- as.data.frame(MA_idx.all@index)
prop_idx$data[prop_idx$data == -1] <- NA

p1 <- ggplot(prop_idx,aes(x=year,y=data))+
      geom_bar(aes(fill = as.factor(age)),stat="identity")+
      theme_bw()+
      labs(x = 'Year',
           y = 'Proportion of idx at age',
           fill = "Age")+
      scale_y_continuous(expand = c(0, 0))

p2 <- cohcorrplot(pil.all.tun[["MA_idx"]]@index)

p <- grid.arrange(p1,p2)

ggsave(paste0(prefix,"_MA_idx_fig.png"),
       p,
       width = 170,
       height = 200,
       units = c("mm"),
       dpi = 300)


# idx proportion per stock
prop_idx <- as.data.frame(MA_idx@index)
prop_idx$data[prop_idx$data == -1] <- NA

p <- ggplot(prop_idx,aes(x=year,y=data))+
  geom_bar(aes(fill = as.factor(age)),stat="identity")+
  theme_bw()+
  labs(x = 'Year',
       y = 'Proportion of idx at age',
       fill = "Age")+
  scale_y_continuous(expand = c(0, 0))+
  facet_wrap(~area)

# p2 <- FLCore::plot(pil.all.tun[["MA_idx"]],type="internal")
# 
# p <- grid.arrange(p1,p2,ncol=2)

ggsave(paste0(prefix,"_MA_idx_stock_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)

# dev.off()

## catch SOP vs caton
df.caton.AB  <- as.data.frame(pil.AB@catch)
df.caton.AB$type <- 'caton'
df.caton.AB$region <- 'AB'
df.SOP.AB    <- as.data.frame(computeCatch(pil.AB))
df.SOP.AB$type <- 'SOP'
df.SOP.AB$region <- 'AB'

df.caton.C  <- as.data.frame(pil.C@catch)
df.caton.C$type <- 'caton'
df.caton.C$region <- 'C'
df.SOP.C    <- as.data.frame(computeCatch(pil.C))
df.SOP.C$type <- 'SOP'
df.SOP.C$region <- 'C'

df.plot <- rbind(df.caton.AB,df.SOP.AB,df.caton.C,df.SOP.C)

p1 <- ggplot(df.plot,aes(x=year,y=data*1e-6,col=type))+
  theme_bw()+
  geom_line()+
  ylim(0,NA)+
  ylab('Catch [1e6]')+
  labs(tag = "(a)")+
  theme(plot.tag = element_text(face = "bold"))+
  facet_wrap(~region)

df.caton  <- as.data.frame(pil.all@catch)
df.caton$type <- 'caton'
df.SOP    <- as.data.frame(computeCatch(pil.all))
df.SOP$type <- 'SOP'
df.plot <- rbind(df.caton,df.SOP)

p2 <- ggplot(df.plot,aes(x=year,y=data*1e-6,col=type))+
  theme_bw()+
  geom_line()+
  ylim(0,NA)+
  labs(tag = "(b)")+
  theme(plot.tag = element_text(face = "bold"))+
  ylab('Catch [1e6]')

df.plot <- as.data.frame((pil.all@catch-computeCatch(pil.all))/computeCatch(pil.all))

p3 <- ggplot(df.plot,aes(x=year,y=data*100))+
      theme_bw()+
      geom_line()+
      ylim(-50,50)+
      labs(tag = "(c)")+
      theme(plot.tag = element_text(face = "bold"))+
      ylab('Discrepency in catch SOP (%)')

p <- grid.arrange(p1,p2,p3,
                  layout_matrix = rbind(c(1,1,1,1),
                                        c(2,2,3,3)))

ggsave(paste0(prefix,"_catch_SOP_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)

# dev.off()

## weight in the stock by cohort
# taf.png(paste0(prefix,"_wstock_cohort"),
#         width = 1600, height = 1200)

stk.wt <- as.data.frame(FLCohort(window(pil.all@stock.wt,)))

p <- ggplot(subset(stk.wt,cohort > 2000),aes(x=cohort+age,y=data,col=as.factor(cohort)))+
      theme_bw()+
      geom_line()+
      geom_point()+
      guides(col="none")+
      xlab('Cohort')+
      ylab('Catch weight at age (kg)')

ggsave(paste0(prefix,"_wstock_cohort.png"),
       p,
       width = 140,
       height = 80,
       units = c("mm"),
       dpi = 300)

# dev.off()

## Internal consistency MA_idx
# taf.png(paste0(prefix,"_MA_idx_internal_consistency"))
# zoom(FLCore::plot(pil.all.tun[["MA_idx"]],type="internal"))
# dev.off()

## Catch internal consistency
# taf.png(paste0(prefix,"_catch_internal_consistency"))
# zoom(FLCore::plot(FLIndex(index=pil.all@catch.n),type="internal"))
# dev.off()

## Data map
taf.png(paste0(prefix,"_data_map"),
        width = 1600, height = 1200)
dataplot(pil.samfit)
dev.off()

## Recruitment survey
# taf.png(paste0(prefix,"_rec_survey"))
df.rec.C <- as.data.frame(pil.C.tun$rec_russ@index)
df.rec.C$region <- 'C'
df.rec.AB <- as.data.frame(pil.AB.tun$rec_russ@index)
df.rec.AB$region <- 'AB'
df.rec.all <- as.data.frame(pil.all.tun$rec_russ@index)
df.rec.all$region <- 'all'
df.rec <- rbind(df.rec.AB,df.rec.C,df.rec.all)

p <- ggplot(df.rec,aes(x=year,y=data,col=region))+
  theme_bw()+
  geom_point()+
  geom_line()+
  ylab('Recruitment index')+
  ylim(0,NA)

ggsave(paste0(prefix,"_rec_survey_fig.png"),
       p,
       width = 140,
       height = 80,
       units = c("mm"),
       dpi = 300)

## Recruitment survey
# taf.png(paste0(prefix,"_rec_survey"))
df.rec.C <- as.data.frame(pil.C.tun$rec_russ@index)
df.rec.C$region <- 'C'
df.rec.C$survey <- 'RU'
df.rec.AB <- as.data.frame(pil.AB.tun$rec_russ@index)
df.rec.AB$region <- 'AB'
df.rec.AB$survey <- 'RU'
df.rec.all <- as.data.frame(pil.all.tun$rec_russ@index)
df.rec.all$region <- 'all'
df.rec.all$survey <- 'RU'
df.rec <- rbind(df.rec.AB,df.rec.C,df.rec.all)

df.rec.C <- as.data.frame(pil.C.tun$MA_idx@index[ac(0)])
df.rec.C$region <- 'C'
df.rec.C$survey <- 'MA'
df.rec.AB <- as.data.frame(pil.AB.tun$MA_idx@index[ac(0)])
df.rec.AB$region <- 'AB'
df.rec.AB$survey <- 'MA'
df.rec.all <- as.data.frame(pil.all.tun$MA_idx@index[ac(0)])
df.rec.all$region <- 'all'
df.rec.all$survey <- 'MA'
df.rec <- rbind(df.rec,df.rec.AB,df.rec.C,df.rec.all)

p <- ggplot(df.rec,aes(x=year,y=data,col=region))+
  theme_bw()+
  geom_point()+
  geom_line()+
  ylab('Recruitment index')+
  ylim(0,NA)+
  facet_wrap(~survey,scales='free_y',ncol=1)

ggsave(paste0(prefix,"_rec_survey_comp_fig.png"),
       p,
       width = 170,
       height = 120,
       units = c("mm"),
       dpi = 300)


## SSB survey
# taf.png(paste0(prefix,"_ssb_survey"))
df.ssb.C <- as.data.frame(pil.C.tun$MA_ssb@index)
df.ssb.C$region <- 'C'
df.ssb.AB <- as.data.frame(pil.AB.tun$MA_ssb@index)
df.ssb.AB$region <- 'AB'
df.ssb.all <- as.data.frame(pil.all.tun$MA_ssb@index)
df.ssb.all$region <- 'all'
df.ssb <- rbind(df.ssb.all,df.ssb.AB,df.ssb.C)

p1 <- ggplot(df.ssb,aes(x=year,y=data,col=region))+
  theme_bw()+
  geom_point()+
  geom_line()+
  ylab('MA survey SSB index')+
  ylim(0,NA)

df.ssb.MA <- as.data.frame(pil.all.tun$MA_ssb@index)
df.ssb.MA$survey <- 'MA'
df.ssb.nansen <- as.data.frame(pil.all.tun$nansen@index)
df.ssb.nansen$survey <- 'nansen'
df.ssb <- rbind(df.ssb.nansen,df.ssb.MA)

p2 <- ggplot(df.ssb,aes(x=year,y=data,col=survey))+
  theme_bw()+
  geom_point()+
  geom_line()+
  ylab('SSB survey index')+
  ylim(0,NA)

p <- grid.arrange(p1,p2,ncol=1)

ggsave(paste0(prefix,"_ssb_survey_fig.png"),
       p,
       width = 170,
       height = 120,
       units = c("mm"),
       dpi = 300)

# dev.off()

# taf.png(paste0(prefix,"_ssb_survey_all"))

# dev.off()

setwd('../..')

################################################################################
### ============================================================================
### Single-fleet
### ============================================================================
################################################################################

prefix <- 'assessment'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

##################################################
## assessment fit
##################################################
pdf(file.path(paste0(prefix,"_fit.pdf")))
residual.diagnostics(pil.sam)
dev.off()

##################################################
## stock trajectory
##################################################
taf.png(paste0(prefix,"_stock_trajectory"),
        width = 1600, height = 1200)
zoom(plot(pil.sam))
dev.off()

##################################################
## stock trajectory AB
##################################################
taf.png(paste0(prefix,"_stock_trajectory_AB"),
        width = 1600, height = 1200)
zoom(plot(pil.sam.AB))
dev.off()

##################################################
## stock trajectory C
##################################################
taf.png(paste0(prefix,"_stock_trajectory_C"),
        width = 1600, height = 1200)
zoom(plot(pil.sam.C))
dev.off()


##################################################
## Catchability at age
##################################################
# taf.png(paste0(prefix,"_survey_catchability"),
#         width = 1600, height = 1200)
q <- catchabilities(pil.sam)
p <- ggplot(q,aes(x=age,y=value))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=lbnd, ymax=ubnd),alpha=0.2)+
  facet_wrap(~fleet,scales = 'free_y')+
  theme_bw()

ggsave(paste0(prefix,"_survey_catchability_fig.png"),
       p,
       width = 170,
       height = 100,
       units = c("mm"),
       dpi = 300)
# dev.off()

##################################################
## observation variance by data source
##################################################
# taf.png(paste0(prefix,"_observation_var_by_source"))
obv <- obs.var(pil.sam)
obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
obv$str <- factor(obv$str,levels = obv$str[order(obv$value)])
# obv$fleetAge <- paste0(obv$fleet,'-',obv$age)

p <- ggplot(subset(obv,fleet != 'caton'),aes(x=str,y=value,fill=str))+
      theme_bw()+
      geom_bar(stat="identity", alpha=0.7)+
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      ylab('Observation variance by source')+
      xlab('')+
      geom_errorbar(aes(ymin=lbnd, ymax=ubnd))

ggsave(paste0(prefix,"_observation_var_by_source.png"),
       p,
       width = 170,
       height = 100,
       units = c("mm"),
       dpi = 300)
# dev.off()

##################################################
## variance vs uncertainty for each data source
#################################################
taf.png(paste0(prefix,"_variance_vs_uncertainty"))
plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
     pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
dev.off()

##################################################
## Selectivity pattern per 5 years
#################################################
taf.png(paste0(prefix,"_selectivity"),
        width = 1600, height = 1200)
sel.pat <- merge(f(pil.sam),fbar(pil.sam),
                 by="year",suffixes=c(".f",".fbar"))
sel.pat$sel <- sel.pat$value.f/sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))
zoom(xyplot(sel ~ age|sprintf("%i's",floor((year)/5)*5),sel.pat,
       groups=year,type="l",as.table=TRUE,
       scale=list(alternating=FALSE),
       main="Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar"))
dev.off()

##################################################
## correlation matrix model params
#################################################
taf.png(paste0(prefix,"_cor_params"),
        width = 1600, height = 1200)
cor.plot(pil.sam)
dev.off()

##################################################
## catch residuals per year per age
#################################################
dat <- subset(residuals(pil.sam),fleet=="catch unique")
dat$std.res[dat$std.res == 0] <- NA
p1 <- ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue"))+
        ggtitle('Catch at age')

dat <- subset(residuals(pil.sam),fleet=="nansen")
dat$std.res[dat$std.res == 0] <- NA
p2 <- ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue"))+
        ggtitle('Nansen SSB')

dat <- subset(residuals(pil.sam),fleet=="rec_russ")
dat$std.res[dat$std.res == 0] <- NA
p3 <- ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue"))+
        ggtitle('Rec idx')

dat <- subset(residuals(pil.sam),fleet=="MA_idx")
dat$std.res[dat$std.res == 0] <- NA
p4 <- ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue"))+
        ggtitle('MA idx')

dat <- subset(residuals(pil.sam),fleet=="MA_ssb")
dat$std.res[dat$std.res == 0] <- NA
p5 <- ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
  geom_point(shape=21)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_manual(values = c("red", "blue"))+
  ggtitle('MA SSB')

p <- grid.arrange(p1,p2,p3,p4,p5)

ggsave(paste0(prefix,"_residuals.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)

##################################################
## nansen residuals per year per age
#################################################
taf.png(paste0(prefix,"_nansen_residuals"))
dat <- subset(residuals(pil.sam),fleet=="nansen")
print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue")))
dev.off()

##################################################
## MA_idx residuals per year per age
#################################################
taf.png(paste0(prefix,"_MA_idx_residuals"))
dat <- subset(residuals(pil.sam),fleet=="MA_idx")
print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue")))
dev.off()

##################################################
## IBTS0 residuals per year per age
#################################################
taf.png(paste0(prefix,"_MA_ssb_residuals"))
dat <- subset(residuals(pil.sam),fleet=="MA_ssb")
print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue")))
dev.off()

##################################################
## F at age
#################################################
# taf.png(paste0(prefix,"_fatage"))

df.plot <- as.data.frame(pil.all@harvest)

p <- ggplot(subset(df.plot,year > 2000),aes(x=year,y=data,col=as.factor(age)))+
  theme_bw()+
  guides(col=guide_legend(title="age"))+
  ylab('F')+
  geom_line()+
  geom_point()
  
ggsave(paste0(prefix,"_F_fig.png"),
       p,
       width = 80,
       height = 80,
       units = c("mm"),
       dpi = 300)
  
# dev.off()

##################################################
## F/fbar at age
#################################################
pil.temp <- pil.all
pil.temp@harvest <- sweep(pil.all@harvest,2:6,fbar(pil.all),'/')

taf.png(paste0(prefix,"_ffbaratage"))
zoom(timeseries(window(pil.temp,2000,range(pil.all)["maxyear"]),slot="harvest"))
dev.off()

##################################################
## Recruitment vs SSB
#################################################

df.rec <- rec(pil.sam)
df.rec <- df.rec %>% select(c('year','value'))
colnames(df.rec) <- c('year','rec')
df.ssb <- ssb(pil.sam)
df.ssb <- df.ssb %>% select(c('year','value'))
colnames(df.ssb) <- c('year','ssb')

df.all <- left_join(df.rec,df.ssb,by=('year'))
df.all$yearclass <- df.all$year-1
df.all$ssbR <- df.all$ssb*100/df.all$rec
df.all$log.ssbR <- log10(df.all$ssbR)
df.all$RperS <- df.all$rec/df.all$ssb

taf.png(paste0(prefix,"_SR_periods"))
ggplot(data=df.all,aes(x=ssb,y=rec),aes(x=ssb,y=rec))+
  geom_line()+
  ylab('Recruitment 0-wr (thousands)')+
  xlab('Spawning Stock Biomass (t)')+
  theme_bw()
dev.off()

taf.png(paste0(prefix,"_recruit_spawners"))
ggplot()+
  geom_smooth(data=df.all,aes(x=yearclass,y=RperS),span = 0.1,col='black', linewidth=2)+
  geom_point(data=df.all,aes(x=yearclass,y=RperS), size=2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()+
  annotation_logticks()+
  xlab('Year class')+
  ylab('Recruits per spawner')
dev.off()

##################################################
## retro stock trajectory
#################################################
taf.png(paste0(prefix,"_retrospective_stock"))
zoom(plot(pil.retro))
dev.off()

##################################################
## retro ssb
#################################################
# taf.png(paste0(prefix,"_retrospective_SSB"))
df.ssb <- ssb(pil.retro)
df.ssb$quant <- 'ssb'
df.fbar <- fbar(pil.retro)
df.fbar$quant <- 'fbar'
df.rec <- fbar(pil.retro)
df.rec$quant <- 'rec'
df.plot <- rbind(df.ssb,df.fbar,df.rec)

p <- ggplot()+
  theme_bw()+
  geom_line(data=subset(df.plot,name == '2025' & year >=2000),aes(x=year,y=value),linewidth=1)+
  geom_ribbon(data=subset(df.plot,name == '2025' & year >=2000),aes(x=year,y=value,ymin=lbnd,ymax=ubnd),col=NA,alpha=0.6,fill='grey')+
  geom_line(data=subset(df.plot,name != '2025' & year >=2000),aes(x=year,y=value,col=name),linewidth=1)+
  facet_wrap(~quant,scales='free',ncol=1)+
  theme(legend.position="top")
# 
# p1 <- ggplot()+
#   theme_bw()+
#   geom_line(data=subset(df.plot,name == '2024'),aes(x=year,y=value),linewidth=1)+
#   geom_ribbon(data=subset(df.plot,name == '2024'),aes(x=year,ymin=lbnd,ymax=ubnd),col=NA,alpha=0.6,fill='grey')+
#   geom_line(data=subset(df.plot,name != '2024'),aes(x=year,y=value,col=name),linewidth=1)+
#   ylab('SSB')
# 
# df.plot <- fbar(pil.retro)
# p2 <- ggplot()+
#   theme_bw()+
#   geom_line(data=subset(df.plot,name == '2024'),aes(x=year,y=value),linewidth=1)+
#   geom_ribbon(data=subset(df.plot,name == '2024'),aes(x=year,ymin=lbnd,ymax=ubnd),col=NA,alpha=0.6,fill='grey')+
#   geom_line(data=subset(df.plot,name != '2024'),aes(x=year,y=value,col=name),linewidth=1)+
#   ylab('Fbar')
# 
# p <- grid.arrange(p1,p2)
  
ggsave(paste0(prefix,"_retro_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)

# dev.off()

##################################################
## retro ssb
#################################################
taf.png(paste0(prefix,"_retrospective_ssb"))
df.plot <- ssb(pil.retro)
ggplot()+
  theme_bw()+
  geom_line(data=subset(df.plot,name == '2025'),aes(x=year,y=value),linewidth=1)+
  geom_ribbon(data=subset(df.plot,name == '2025'),aes(x=year,ymin=lbnd,ymax=ubnd),col=NA,alpha=0.6,fill='grey')+
  geom_line(data=subset(df.plot,name != '2025'),aes(x=year,y=value,col=name),linewidth=1)+
  ylab('SSB')
dev.off()


##################################################
## retro fbar
#################################################
taf.png(paste0(prefix,"_retrospective_fbar"))
df.plot <- fbar(pil.retro)
ggplot()+
  theme_bw()+
  geom_line(data=subset(df.plot,name == '2025'),aes(x=year,y=value),linewidth=1)+
  geom_ribbon(data=subset(df.plot,name == '2025'),aes(x=year,ymin=lbnd,ymax=ubnd),col=NA,alpha=0.6,fill='grey')+
  geom_line(data=subset(df.plot,name != '2025'),aes(x=year,y=value,col=name),linewidth=1)+
  ylab('Fbar')
dev.off()

##################################################
## retro rec
#################################################
taf.png(paste0(prefix,"_retrospective_rec"))
df.plot <- rec(pil.retro)
ggplot()+
  theme_bw()+
  geom_line(data=subset(df.plot,name == '2025'),aes(x=year,y=value),linewidth=1)+
  geom_ribbon(data=subset(df.plot,name == '2025'),aes(x=year,ymin=lbnd,ymax=ubnd),col=NA,alpha=0.6,fill='grey')+
  geom_line(data=subset(df.plot,name != '2025'),aes(x=year,y=value,col=name),linewidth=1)+
  ylab('Recruitment')
dev.off()


##################################################
## model parameters retrospective
#################################################
taf.png(paste0(prefix,"_retrospective_SAM_parameters"))
zoom(retroParams(pil.retro))
dev.off()

##################################################
## fishing selectivity retrospective
#################################################
taf.png(paste0(prefix,"_retrospective_SAM_parameters"))
zoom(retroSelectivity(pil.retro,(range(pil.all)["maxyear"]-10):range(pil.all)["maxyear"]))
dev.off()

##################################################
## catch fit
#################################################
taf.png(paste0(prefix,"_catch_fit_SOP"),
        width = 1600, height = 1200)
catchplot(pil.samfit, addCI=TRUE)
dev.off()

##################################################
## traj fit
#################################################
taf.png(paste0(prefix,"_ssb_fit"),
        width = 1600, height = 1200)
ssbplot(pil.samfit, addCI=TRUE)
dev.off()

taf.png(paste0(prefix,"_fbar_fit"),
        width = 1600, height = 1200)
fbarplot(pil.samfit, addCI=TRUE)
dev.off()

taf.png(paste0(prefix,"_rec_fit"),
        width = 1600, height = 1200)
recplot(pil.samfit, addCI=TRUE)
dev.off()

loo <- leaveout(pil.samfit,fleet=as.list(c(2:4,6)))

taf.png(paste0(prefix,"_loo_ssb_fit"),
        width = 1600, height = 1200)
ssbplot(loo, addCI=F)
dev.off()

taf.png(paste0(prefix,"_loo_fbar_fit"),
        width = 1600, height = 1200)
fbarplot(loo, addCI=F)
dev.off()

taf.png(paste0(prefix,"_loo_rec_fit"),
        width = 1600, height = 1200)
recplot(loo, addCI=F)
dev.off()

taf.png(paste0(prefix,"_MA_idx_fit"))
fitplot(pil.samfit, fleets=2)
dev.off()

taf.png(paste0(prefix,"_rec_fit"))
fitplot(pil.samfit, fleets=4)
dev.off()

taf.png(paste0(prefix,"_MA_ssb_fit"))
fitplot(pil.samfit, fleets=3)
dev.off()

taf.png(paste0(prefix,"_nansen_ssb_fit"))
fitplot(pil.samfit, fleets=5)
dev.off()

taf.png(paste0(prefix,"_catch_fit"))
fitplot(pil.samfit, fleets=1)
dev.off()

##################################################
## compare trajectories between regions
#################################################
# taf.png(paste0(prefix,"_compare_regions"))
df.ssb <- as.data.frame(ssb(pil.C))
df.ssb$quant <- 'ssb'
df.ssb$region <- 'C'
df.fbar <- as.data.frame(fbar(pil.C))
df.fbar$quant <- 'fbar'
df.fbar$region <- 'C'
df.rec <- as.data.frame(rec(pil.C))
df.rec$quant <- 'rec'
df.rec$region <- 'C'

df.all <- rbind(df.ssb,df.fbar,df.rec)

df.ssb <- as.data.frame(ssb(pil.AB))
df.ssb$quant <- 'ssb'
df.ssb$region <- 'AB'
df.fbar <- as.data.frame(fbar(pil.AB))
df.fbar$quant <- 'fbar'
df.fbar$region <- 'AB'
df.rec <- as.data.frame(rec(pil.AB))
df.rec$quant <- 'rec'
df.rec$region <- 'AB'

df.all <- rbind(df.all,df.ssb,df.fbar,df.rec)

df.ssb <- as.data.frame(ssb(pil.all))
df.ssb$quant <- 'ssb'
df.ssb$region <- 'all'
df.fbar <- as.data.frame(fbar(pil.all))
df.fbar$quant <- 'fbar'
df.fbar$region <- 'all'
df.rec <- as.data.frame(rec(pil.all))
df.rec$quant <- 'rec'
df.rec$region <- 'all'

df.all <- rbind(df.all,df.ssb,df.fbar,df.rec)

p <- ggplot(df.all,aes(x=year,y=data,col=region))+
  theme_bw()+
  geom_line()+
  facet_wrap(~quant,scales='free',ncol=1)+
  ylim(0,NA)
# dev.off()

ggsave(paste0(prefix,"_compare_regions_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)

setwd('../..')

################################################################################
### ============================================================================
### ref points
### ============================================================================
################################################################################

library(msy)

prefix <- 'refpts'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

FIT_ricker <- eqsr_fit(pil.all,nsamp=2000, models = "Ricker")
FIT_segreg <- eqsr_fit(pil.all,nsamp=2000, models = "Segreg")

FIT_ricker.AB <- eqsr_fit(pil.AB,nsamp=2000, models = "Ricker")
FIT_segreg.AB <- eqsr_fit(pil.AB,nsamp=2000, models = "Segreg")

FIT_ricker.C <- eqsr_fit(pil.C,nsamp=2000, models = "Ricker")
FIT_segreg.C <- eqsr_fit(pil.C,nsamp=2000, models = "Segreg")

# taf.png(paste0(prefix,"_ricker"))
p1 <- eqsr_plot(FIT_ricker, ggPlot=TRUE)+
  theme(legend.position="none")+
  ggtitle('Ricker')

p2 <- eqsr_plot(FIT_segreg, ggPlot=TRUE)+
  theme(legend.position="none")+
  ggtitle('Segmented regression')

p <- grid.arrange(p1,p2)

ggsave(paste0(prefix,"_fit_SRR_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)

p1 <- eqsr_plot(FIT_ricker.AB, ggPlot=TRUE)+
  theme(legend.position="none")+
  ggtitle('Ricker')

p2 <- eqsr_plot(FIT_segreg.AB, ggPlot=TRUE)+
  theme(legend.position="none")+
  ggtitle('Segmented regression')

p <- grid.arrange(p1,p2)

ggsave(paste0(prefix,"_fit_SRR_AB_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)

p1 <- eqsr_plot(FIT_ricker.C, ggPlot=TRUE)+
  theme(legend.position="none")+
  ggtitle('Ricker')

p2 <- eqsr_plot(FIT_segreg.C, ggPlot=TRUE)+
  theme(legend.position="none")+
  ggtitle('Segmented regression')

p <- grid.arrange(p1,p2)

ggsave(paste0(prefix,"_fit_SRR_C_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)


# taf.png(paste0(prefix,"_MSY"))
p <- FLRef::ploteq(FLRef::computeFbrp(pil.all, pil.sr.rick),obs=TRUE,refpts="msy")

ggsave(paste0(prefix,"_MSY_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)

# dev.off()

# taf.png(paste0(prefix,"_ssb_status"))
df.ssb <- ssb(pil.sam)
p1 <- ggplot(df.ssb,aes(x=year,y=value))+
      theme_bw()+
      geom_line()+
      geom_ribbon(aes(ymin = lbnd, ymax = ubnd),colour = NA,alpha=0.3)+
      ylim(0,NA)+
      geom_hline(yintercept = c(refpts['SBlim']),col='red')+
      geom_text(aes(1990,c(refpts['SBlim']),label = 'Blim', vjust = 1),col='red')+
      geom_hline(yintercept = c(refpts['SBMSY']),col='blue')+
      geom_text(aes(1990,c(refpts['SBMSY']),label = 'SBMSY', vjust = 1),col='blue')+
      ylab('SSB')

df.fbar <- fbar(pil.sam)
p2 <- ggplot(df.fbar,aes(x=year,y=value))+
      theme_bw()+
      geom_line()+
      geom_ribbon(aes(ymin = lbnd, ymax = ubnd),colour = NA,alpha=0.3)+
      ylim(0,NA)+
      geom_hline(yintercept = c(refpts['FMSY']),col='blue')+
      geom_text(aes(1990,c(refpts['FMSY']),label = 'FMSY', vjust = 1),col='blue')+
      ylab('fbar')

p <- grid.arrange(p1,p2)

ggsave(paste0(prefix,"_stock_status_fig.png"),
       p,
       width = 170,
       height = 120,
       units = c("mm"),
       dpi = 300)

# dev.off()

# stock status
# taf.png(paste0(prefix,"_compare_status"))
df.ssb <- as.data.frame(ssb(pil.C)/refpts.C['SBMSY'])
df.ssb$quant <- 'ssb/SBMSY'
df.ssb$region <- 'C'
df.fbar <- as.data.frame(fbar(pil.C)/refpts.C['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'C'

df.all <- rbind(df.ssb,df.fbar)

df.ssb <- as.data.frame(ssb(pil.AB)/refpts.AB['SBMSY'])
df.ssb$quant <- 'ssb/SBMSY'
df.ssb$region <- 'AB'
df.fbar <- as.data.frame(fbar(pil.AB)/refpts.AB['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'AB'

df.all <- rbind(df.all,df.ssb,df.fbar)

df.ssb <- as.data.frame(ssb(pil.all)/refpts['SBMSY'])
df.ssb$quant <- 'ssb/SBMSY'
df.ssb$region <- 'all'
df.fbar <- as.data.frame(fbar(pil.all)/refpts['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'all'

df.all <- rbind(df.all,df.ssb,df.fbar)

p <- ggplot(df.all,aes(x=year,y=data,col=region))+
  theme_bw()+
  geom_hline(yintercept = 1,col='red')+
  geom_line()+
  facet_wrap(~quant,scales='free',ncol=1)+
  ylim(0,NA)

ggsave(paste0(prefix,"_compare_status_fig.png"),
       p,
       width = 170,
       height = 120,
       units = c("mm"),
       dpi = 300)

# stock status
# taf.png(paste0(prefix,"_compare_status"))
df.ssb <- as.data.frame(ssb(pil.C)/refpts.C['SBMSY']/0.9)
df.ssb$quant <- 'ssb/SB0.1'
df.ssb$region <- 'C'
df.fbar <- as.data.frame(fbar(pil.C)/refpts.C['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'C'

df.all <- rbind(df.ssb,df.fbar)

df.ssb <- as.data.frame(ssb(pil.AB)/refpts.AB['SBMSY']/0.9)
df.ssb$quant <- 'ssb/SB0.1'
df.ssb$region <- 'AB'
df.fbar <- as.data.frame(fbar(pil.AB)/refpts.AB['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'AB'

df.all <- rbind(df.all,df.ssb,df.fbar)

df.ssb <- as.data.frame(ssb(pil.all)/refpts['SBMSY']/0.9)
df.ssb$quant <- 'ssb/SB0.1'
df.ssb$region <- 'all'
df.fbar <- as.data.frame(fbar(pil.all)/refpts['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'all'

df.all <- rbind(df.all,df.ssb,df.fbar)

p <- ggplot(df.all,aes(x=year,y=data,col=region))+
  theme_bw()+
  geom_hline(yintercept = 1,col='red')+
  geom_line()+
  facet_wrap(~quant,scales='free',ncol=1)+
  ylim(0,NA)

ggsave(paste0(prefix,"_compare_status_SB0.1_fig.png"),
       p,
       width = 170,
       height = 120,
       units = c("mm"),
       dpi = 300)

# stock status
# taf.png(paste0(prefix,"_compare_status"))
df.ssb <- as.data.frame(ssb(pil.C)/refpts.C['SBMSY']/0.3)
df.ssb$quant <- 'ssb/0.3SBMSY'
df.ssb$region <- 'C'
df.fbar <- as.data.frame(fbar(pil.C)/refpts.C['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'C'

df.all <- rbind(df.ssb,df.fbar)

df.ssb <- as.data.frame(ssb(pil.AB)/refpts.AB['SBMSY']/0.3)
df.ssb$quant <- 'ssb/0.3SBMSY'
df.ssb$region <- 'AB'
df.fbar <- as.data.frame(fbar(pil.AB)/refpts.AB['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'AB'

df.all <- rbind(df.all,df.ssb,df.fbar)

df.ssb <- as.data.frame(ssb(pil.all)/refpts['SBMSY']/0.3)
df.ssb$quant <- 'ssb/0.3SBMSY'
df.ssb$region <- 'all'
df.fbar <- as.data.frame(fbar(pil.all)/refpts['FMSY'])
df.fbar$quant <- 'fbar/FMSY'
df.fbar$region <- 'all'

df.all <- rbind(df.all,df.ssb,df.fbar)

p <- ggplot(subset(df.all,quant == 'ssb/0.3SBMSY'),aes(x=year,y=data,col=region))+
  theme_bw()+
  geom_hline(yintercept = 1,col='red')+
  geom_line()+
  facet_wrap(~quant,scales='free',ncol=1)+
  ylim(0,NA)

ggsave(paste0(prefix,"_compare_status_SBlim_fig.png"),
       p,
       width = 170,
       height = 120,
       units = c("mm"),
       dpi = 300)

# plotting kobe plots
# taf.png(paste0(prefix,"_compare_kobe"))
df.kobe.all <- as.data.frame(ssb(pil.all)/refpts['SBMSY']) %>% dplyr::rename(stock=data)
df.kobe.all$harvest <- as.data.frame(fbar(pil.all)/refpts['FMSY'])$data
df.kobe.all$region <- 'all'

df.kobe.AB <- as.data.frame(ssb(pil.AB)/refpts.AB['SBMSY']) %>% dplyr::rename(stock=data)
df.kobe.AB$harvest <- as.data.frame(fbar(pil.AB)/refpts.AB['FMSY'])$data
df.kobe.AB$region <- 'AB'

df.kobe.C <- as.data.frame(ssb(pil.C)/refpts.C['SBMSY']) %>% dplyr::rename(stock=data)
df.kobe.C$harvest <- as.data.frame(fbar(pil.AB)/refpts.C['FMSY'])$data
df.kobe.C$region <- 'C'

kobe.comb <- rbind(df.kobe.all,df.kobe.AB,df.kobe.C)

p1 <- kobePhase(subset(kobe.comb,region == 'all'))+
      theme_bw()+
      geom_text(aes(stock,harvest,label=year),
                hjust=1, vjust=-0.1,alpha=0.5)+
      geom_path(aes(stock,harvest,col=year>2020))+
      geom_point(aes(stock,harvest,col=year>2020),size=3)+
      guides(col="none")+
      ggtitle('Region all')

p2 <- kobePhase(subset(kobe.comb,region == 'AB'))+
      theme_bw()+
      geom_text(aes(stock,harvest,label=year),
                hjust=1, vjust=-0.1,alpha=0.5)+
      geom_path(aes(stock,harvest,col=year>2020))+
      geom_point(aes(stock,harvest,col=year>2020),size=3)+
      guides(col="none")+
      ggtitle('Region AB')

p3 <- kobePhase(subset(kobe.comb,region == 'C'))+
      theme_bw()+
      geom_text(aes(stock,harvest,label=year),
                hjust=1, vjust=-0.1,alpha=0.5)+
      geom_path(aes(stock,harvest,col=year>2020))+
      geom_point(aes(stock,harvest,col=year>2020),size=3)+
      guides(col="none")+
      ggtitle('Region C')

p <- grid.arrange(p1,p2,p3,layout_matrix = rbind(c(1,1,2,2),
                                                 c(NA,3,3,NA)))

ggsave(paste0(prefix,"_compare_kobe_fig.png"),
       p,
       width = 170,
       height = 140,
       units = c("mm"),
       dpi = 300)
# dev.off()

# taf.png(paste0(prefix,"_kobe_all"))
# kobePhase(subset(df.kobe.all,region == 'all'))+
#   theme_bw()+
#   geom_text(aes(stock,harvest,label=year),
#             hjust=1, vjust=-0.1,alpha=0.5)+
#   geom_path(aes(stock,harvest,col=year>2020))+
#   geom_point(aes(stock,harvest,col=year>2020),size=3)+
#   guides(col="none")+
#   ggtitle('Region all')
# dev.off()
# 
# taf.png(paste0(prefix,"_kobe_AB"))
# kobePhase(df.kobe.AB)+
#   theme_bw()+
#   geom_text(aes(stock,harvest,label=year),
#             hjust=1, vjust=-0.1,alpha=0.5)+
#   geom_path(aes(stock,harvest,col=year>2020))+
#   geom_point(aes(stock,harvest,col=year>2020),size=3)+
#   guides(col="none")+
#   ggtitle('Region AB')
# dev.off()
# 
# taf.png(paste0(prefix,"_kobe_C"))
# kobePhase(df.kobe.C)+
#   theme_bw()+
#   geom_text(aes(stock,harvest,label=year),
#             hjust=1, vjust=-0.1,alpha=0.5)+
#   geom_path(aes(stock,harvest,col=year>2020))+
#   geom_point(aes(stock,harvest,col=year>2020),size=3)+
#   guides(col="none")+
#   ggtitle('Region C')
# dev.off()

stock.status <- subset(kobe.comb,year >= (2024-5)) %>% 
                dplyr::group_by(region) %>% 
                dplyr::summarize(SSB_SBMSY=mean(stock),
                          F_FMSY=mean(harvest))

setwd('../..')

################################################################################
### ============================================================================
### stf
### ============================================================================
################################################################################

load(file.path(modelPath,paste0(runName,'_stf.RData')))

prefix <- 'stf'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

flagFirst <- T
for(idxName in names(run.fwd)){
  df.ssb <- as.data.frame(ssb(run.fwd[[idxName]]))
  df.ssb$run <- idxName
  df.ssb$quant <- 'ssb'
  
  df.fbar <- as.data.frame(fbar(run.fwd[[idxName]]))
  df.fbar$run <- idxName
  df.fbar$quant <- 'fbar'
  
  df.rec <- as.data.frame(rec(run.fwd[[idxName]]))
  df.rec$run <- idxName
  df.rec$quant <- 'rec'
  
  df.C <- as.data.frame(catch(run.fwd[[idxName]]))
  df.C$run <- idxName
  df.C$quant <- 'catch'
  
  df.comb <- rbind(df.ssb,df.fbar,df.rec,df.C)
  
  if(flagFirst){
    df.all <- df.comb
    flagFirst <- F
  }else{
    df.all <- rbind(df.all,df.comb)
  }
}

# taf.png(paste0(prefix,"_stf_ssb"))
p1 <- ggplot(subset(df.all,year >= 2020 & quant == 'ssb'),aes(x=year,y=data*1e-6,col=run,fill=run))+
  theme_bw()+
  stat_summary(geom = "line", fun = mean, linetype = "solid",linewidth=1)+
  stat_summary(geom = "ribbon", fun.min = function(z) { quantile(z,0.05) },
               fun.max = function(z) { quantile(z,0.95) }, alpha=.2,colour = NA)+
  geom_hline(yintercept = c(refpts['SBlim']*1e-6),col='red')+
  geom_text(aes(2021,c(refpts['SBlim']*1e-6),label = 'Blim', vjust = 1),col='red')+
  geom_hline(yintercept = c(refpts['SBMSY']*1e-6),col='blue')+
  geom_text(aes(2021,c(refpts['SBMSY']*1e-6),label = 'SBMSY', vjust = 1),col='blue')+
  geom_vline(xintercept = 2024,col='black',linewidth=1)+
  ylab('SSB [1e6 t]')+
  scale_x_continuous(breaks = seq(2020, 2030, by = 3),
                     minor_breaks = seq(2020, 2029, 1))+
  ylim(0,NA)+
  theme(legend.position="top")

p2 <- ggplot(subset(df.all,year >= 2020 & quant == 'fbar'),aes(x=year,y=data,col=run,fill=run))+
  theme_bw()+
  stat_summary(geom = "line", fun = mean, linetype = "solid",linewidth=1)+
  stat_summary(geom = "ribbon", fun.min = function(z) { quantile(z,0.05) },
               fun.max = function(z) { quantile(z,0.95) }, alpha=.2,colour = NA)+
  geom_hline(yintercept = c(refpts['FMSY']),col='blue')+
  geom_text(aes(2021,c(refpts['FMSY']),label = 'FMSY', vjust = 1),col='blue')+
  geom_vline(xintercept = 2024,col='black',linewidth=1)+
  ylab('Fbar (2-4)')+
  scale_x_continuous(breaks = seq(2020, 2030, by = 3),
                     minor_breaks = seq(2020, 2029, 1))+
  ylim(0,NA)+
  theme(legend.position="none")

p3 <- ggplot(subset(df.all,year >= 2020 & quant == 'catch'),aes(x=year,y=data*1e6,col=run,fill=run))+
  theme_bw()+
  stat_summary(geom = "line", fun = mean, linetype = "solid",linewidth=1)+
  stat_summary(geom = "ribbon", fun.min = function(z) { quantile(z,0.05) },
               fun.max = function(z) { quantile(z,0.95) }, alpha=.2,colour = NA)+
  geom_vline(xintercept = 2024,col='black',linewidth=1)+
  ylab('Catch [1e6 t]')+
  scale_x_continuous(breaks = seq(2020, 2030, by = 3),
                     minor_breaks = seq(2020, 2029, 1))+
  ylim(0,NA)+
  theme(legend.position="none")

df.temp <- df.all %>% group_by(quant) %>% summarize(SBlim=NA,
                                         SBMSY=NA,
                                         FMSY=NA)
df.temp[df.temp$quant == 'ssb','SBlim'] <- refpts['SBlim']
df.temp[df.temp$quant == 'ssb','SBMSY'] <- refpts['SBMSY']
df.temp[df.temp$quant == 'fbar','FMSY'] <- refpts['FMSY']
df.temp$run <- 'Csq'
df.temp <- subset(df.temp,quant != 'rec')

p <- ggplot(subset(df.all,quant != 'rec' & year >= 2020),aes(x=year,y=data,col=run,fill=run))+
  theme_bw()+
  stat_summary(geom = "line", fun = mean, linetype = "solid",linewidth=1)+
  stat_summary(geom = "ribbon", fun.min = function(z) { quantile(z,0.05) },
               fun.max = function(z) { quantile(z,0.95) }, alpha=.2,colour = NA)+
  geom_vline(xintercept = 2024,col='black',linewidth=1)+
  scale_x_continuous(breaks = seq(2020, 2030, by = 3),
                     minor_breaks = seq(2020, 2029, 1))+
  geom_hline(data=df.temp,aes(yintercept = SBlim),col='red')+
  geom_text(data=df.temp,
            aes(x=2021,y = SBlim,label = 'SBlim'), vjust = 1,col='red')+
  geom_hline(data=df.temp,aes(yintercept = SBMSY),col='blue')+
  geom_text(data=df.temp,
            aes(x=2021,y = SBMSY,label = 'SBMSY'), vjust = -1,col='blue')+
  geom_hline(data=df.temp,aes(yintercept = FMSY),col='blue')+
  geom_text(data=df.temp,
            aes(x=2021,y = FMSY,label = 'FMSY'), vjust = 1,col='blue')+
  facet_wrap(~quant,scales='free',ncol=1)+
  ylim(0,NA)+
  theme(legend.position="top")

ggsave(paste0(prefix,"_stf_fig.png"),
       p,
       width = 170,
       height = 170,
       units = c("mm"),
       dpi = 300)
# dev.off()


setwd('../..')

