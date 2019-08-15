
df <- read.csv(paste0(gwl_dir,"data_moef_demo_20190522.csv"))
lulc <- raster(paste0(data_dir,"vegetation_change_maps/block14/lulc_block_14.tif"))

head(df)
dates_ground <- as.Date(substr(names(df)[6:ncol(df)],2,11),format = "%d.%m.%Y")

spdf <- SpatialPointsDataFrame(coords = df[,c("Longitude","Latitude")],
                               data = df,
                               proj4string = CRS("+init=epsg:4326")
                               )
plot(spdf)

map_dir <- paste0(smm_dir,"non_aceh/kalimantan_gwl_points/")

#postprocess_pysmm(map_dir)

base <- "close_SMCmap_"

system(sprintf("gdalbuildvrt -separate %s %s",
               paste0(map_dir,"stack.vrt"),
               paste0(map_dir,base,"*.tif")
))

bands <- list.files(map_dir,pattern=glob2rx(paste0(base,"*.tif")))

dates_space <- as.Date(unlist(substr(bands,nchar(base)+1,nchar(bands)-4)),format="%Y_%m_%d")

stack <- brick(paste0(map_dir,"stack.vrt"))
tt    <- data.frame(extract(x=stack,y=spdf))
lu    <- data.frame(extract(x=lulc,y=spdf))

d0 <- data.frame(cbind(spdf@data,tt,lu))
names(d0) <- c("No","JENIS.TP_TMAT","KODE_TITIK","Longitude","Latitude",
               paste0("g_",dates_ground),
               paste0("s_",dates_space),
               "lu_2016")
head(d0)
selection  <- unlist(lapply(lapply(dates_space,function(x){abs(x - dates_ground)}),min)) < 10
min_dates  <- unlist(lapply(lapply(dates_space,function(x){abs(x - dates_ground)}),which.min))
table(d0$lu_2016)
dates_space[selection]
dates_ground[min_dates[selection]]

for(match in which(selection)){
  ds <- paste0("s_",dates_space[match])
  dg <- paste0("g_",dates_ground[min_dates[match]])
  d1 <- d0[!is.na(d0[,ds]) & !is.na(d0[,dg]),c(ds,dg,"lu_2016")]
  plot(d1[,c(ds,dg)])
  points(d1[d1$lu_2016==2010,],col="blue")
  points(d1[d1$lu_2016==2007,],col="green")
  points(d1[d1$lu_2016==2006,],col="red")
  
}
