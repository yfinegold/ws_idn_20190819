
### SET UP PARAMETERS FOR ANALYSIS

## COORDINATES TO ASSESS SOIL MOISTURE TIME SERIES
crd <- cbind(114.06038,-3.22315)

## DEFINE SEASONS - NOTE IT USES THE FIRST OF EACH MONTH 
dry.months <- c(6,7,8,9)
rainy.months <- c(1,2,3,4,5,10,11,12)

###############################################################
###############################################################
## LOAD THE PARAMETERS
source('~/ws_idn_20190819/scripts/s0_parameters.R')
library(ggplot2)
library(dplyr)

## READ THE SMM DATA
smm_dir <- '~/ws_idn_20190819/data/smm_phu/all_smm/'
all.smm.ras <-  brick(paste0(smm_dir,'all_phu_smm.tif'))
all.smm.dat <-  read.csv(paste0(smm_dir,'all_phu_smm.csv'))
# plot(all.smm.ras[[99]])
head(all.smm.dat[,1:3])

###############################################################
###############################################################
## FIRST ASSESS ALL THE TIME SERIES
ts.date <- unlist(all.smm.dat[,2])
x <- raster::extract(all.smm.ras,crd)
y <- rbind(x,ts.date)
z <- t(y)
colnames(z) <- c('ts','date')
z <- as.data.frame(z[z[,1]>0,])
z$date <- as.Date(z$date)
z$ts <- as.numeric(z$ts)
### explore the data
# plot a histogram of the soil moisture values for one point over time
soil.moisture <- unname(x[1,])
soil.moisture <- soil.moisture[soil.moisture>0]
hist(soil.moisture)
print(paste0('Summary statistics of the soil moisture time series for ', paste0(crd,collapse=' , ')))
print(summary(soil.moisture))
print(paste0('Standard deviation: ',round(sd(soil.moisture,  na.rm = FALSE),digits = 2)))

## PLOT THE DATA
g <- ggplot(z, aes(date, ts)) + 
  geom_line() +
  geom_point() +
  #   theme_ipsum() +
  geom_smooth(method = "lm") +
  scale_x_date(date_labels = "%Y") + 
  xlab("Date") + 
  ylab("Soil moisture value- Sentinel 1") 
# g + geom_text(aes(label = eq), data = dftext, parse = TRUE)
print(g)

m <- lm(ts ~ date, z)
print(paste0('Model: y = ',  format(unname(coef(m)[1]), digits = 2), " + ",format(unname(coef(m)[2]), digits = 2),' * x' ))
print(paste0('slope value: ',format(unname(coef(m)[2]), digits = 2)))
r2 <- format(summary(m)$r.squared, digits = 3)
print(paste0("rsquared: " ,r2))
z.mean <- format(mean(z$ts), digits = 3)
print(paste0('Mean: ', z.mean))
z.sd <- format(sd(z$ts), digits = 3)
print(paste0('Standard deviation: ', z.sd))
z.p <- format(summary(m)$coefficients[8], digits = 3)
print(paste0('P value: ', z.p))


###############################################################
###############################################################
## SECOND SEPERATE DRY AND WET SEASON
minyear <- format(as.Date(min(z$date), format="%d/%m/%Y"),"%Y")
maxyear <- format(as.Date(max(z$date), format="%d/%m/%Y"),"%Y")
dry.sub <- as.Date((paste(minyear:maxyear,rep(c(dry.months[1],dry.months[length(dry.months)]),length(minyear:maxyear)),'01',sep='-')))
dry.sub <- dry.sub[order(dry.sub)]
# str(dry.sub[1:2])

dry.sub.ts <- ts[!ts$date>=dry.sub[1] & ts$date<=dry.sub[2],]
print(dry.sub.ts)

# subset dry and wet season
## this needs to made more flexible to take any number of months of dry/wet
dry.sub.ts <- z[z$date>=dry.sub[1] & z$date<=dry.sub[2]
                |
                  z$date>=dry.sub[3] & z$date<=dry.sub[4]
                |
                  z$date>=dry.sub[5] & z$date<=dry.sub[6]
                |
                  z$date>=dry.sub[7] & z$date<=dry.sub[8]
                |
                  z$date>=dry.sub[9] & z$date<=dry.sub[10]
                #    |
                #    z$date>=dry.sub[11] & z$date<=dry.sub[12]
                ,]  
wet.sub.ts <- z[!z$date %in%dry.sub.ts$date,]

###############################################################
###############################################################
## PLOT DRY SEASON 
g <- ggplot(dry.sub.ts, aes(date, ts)) + 
  geom_line() +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_date(date_labels = "%Y") + 
  xlab("Date") + 
  ylab("Dry season soil moisture value- Sentinel 1") 
print(g)

m <- lm(ts ~ date, dry.sub.ts)
print(paste0('Model: y = ',  format(unname(coef(m)[1]), digits = 2), " + ",format(unname(coef(m)[2]), digits = 2),' * x' ))
print(paste0('slope value: ',format(unname(coef(m)[2]), digits = 2)))
r2 <- format(summary(m)$r.squared, digits = 3)
print(paste0("rsquared: " ,r2))
z.mean <- format(mean(dry.sub.ts$ts), digits = 3)
print(paste0('Mean: ', z.mean))
z.sd <- format(sd(dry.sub.ts$ts), digits = 3)
print(paste0('Standard deviation: ', z.sd))
z.p <- format(summary(m)$coefficients[8], digits = 3)
print(paste0('P value: ', z.p))

###############################################################
###############################################################
## PLOT WET SEASON 
g <- ggplot(wet.sub.ts, aes(date, ts)) + 
  geom_line() +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_date(date_labels = "%Y") + 
  xlab("Date") + 
  ylab("Wet season soil moisture value- Sentinel 1") 
print(g)

m <- lm(ts ~ date, wet.sub.ts)
print(paste0('Model: y = ',  format(unname(coef(m)[1]), digits = 2), " + ",format(unname(coef(m)[2]), digits = 2),' * x' ))
print(paste0('slope value: ',format(unname(coef(m)[2]), digits = 2)))
r2 <- format(summary(m)$r.squared, digits = 3)
print(paste0("rsquared: " ,r2))
z.mean <- format(mean(wet.sub.ts$ts), digits = 3)
print(paste0('Mean: ', z.mean))
z.sd <- format(sd(wet.sub.ts$ts), digits = 3)
print(paste0('Standard deviation: ', z.sd))
z.p <- format(summary(m)$coefficients[8], digits = 3)
print(paste0('P value: ', z.p))

