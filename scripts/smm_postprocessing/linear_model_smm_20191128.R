## LOAD THE PARAMETERS
source('~/ws_idn_20190819/scripts/s0_parameters.R')


## PROCESS THE DATA
## READ THE SOIL MOISTURE MAP FOR ALL THE PHUS
smm_dir <- '~/ws_idn_20190819/data/smm_phu/all_smm/'
all.smm.ras <-  brick(paste0(smm_dir,'all_phu_smm.tif'))
all.smm.dat <-  read.csv(paste0(smm_dir,'all_phu_smm.csv'))
ts.date <- as.Date(unlist(all.smm.dat[,2]))
NAvalue(all.smm.ras) <- 0

## calculate the slope
fun_slope <- function(y) { 
if(all(is.na(y))) {
  NA
} else {
  m = lm(y ~ ts.date, na.action=na.omit); summary(m)$coefficients[2] 
}
}

## calculate the p-value
fun_pvalue <- function(y) { 
if(all(is.na(y))) {
  NA
} else {
  m = lm(y ~ ts.date, na.action=na.omit); summary(m)$coefficients[8] 
}
}
## calculate the mean
fun_mean <- function(x) calc(x, fun = mean, na.rm = T)
## calculate the standard deviation
fun_stdv <- function(x) calc(x, fun = sd, na.rm = T)

# beginCluster()
# slope <- clusterR(all.smm.ras, calc, args=list(fun=fun_slope))
# pvalue <- clusterR(all.smm.ras, calc, args=list(fun=fun_pvalue))
# mean <- clusterR(all.smm.ras, fun_mean)
# stdv <- clusterR(all.smm.ras, fun_stdv)
# endCluster()

slope <- calc(r, fun_slope)
pvalue <- calc(r,fun_pvalue)
mean <- calc(r, fun = mean, na.rm = T)
stdv <- calc(r, fun = sd, na.rm = T)

writeRaster(slope,paste0(smm_dir,"all_smm_slope.tif"),overwrite=T)
print('Completed linear regression slope map')
writeRaster(pvalue,paste0(smm_dir,"all_smm_pvalue.tif"),overwrite=T)
print('Completed linear regression signifance map')
writeRaster(mean,paste0(smm_dir,"all_smm_mean.tif"),overwrite=T)
print('Completed mean map')
writeRaster(stdv,paste0(smm_dir,"all_smm_stdv.tif"),overwrite=T)
print('Completed standard deviation map')
