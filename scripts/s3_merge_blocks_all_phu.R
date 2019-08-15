####################  DEFINE DIRECTORIES
res_dir <- paste0(data_dir,"vegetation_change_maps/")
ope_dir <- paste0(res_dir,"operators/")

#lc <- readOGR(paste0(luc_dir,"PL_2015_2016__Planologi_KLHK__Rev_Topology.shp"))
dbf <- read.dbf(paste0(luc_dir,"PL_2015_2016__Planologi_KLHK__Rev_Topology.dbf"))
table(dbf$PL__2016,dbf$PL16_ID)

# all_phu <- readOGR(paste0(ace_dir,"all_phu_indonesia_latlon.shp"))
# all_phu$id <- row(all_phu@data)[,1]
# writeOGR(all_phu,paste0(ace_dir,"all_phu_indonesia_latlon.shp"),"all_phu_indonesia_latlon","ESRI Shapefile",overwrite_layer = T)
####################  CREATE A PSEUDO COLOR TABLE
cols <- col2rgb(c("black","beige","yellow","orange","red","darkred","palegreen","green2","forestgreen",'darkgreen'))
pct <- data.frame(cbind(c(0:9),
                        cols[1,],
                        cols[2,],
                        cols[3,]
))

write.table(pct,paste0(res_dir,"color_table.txt"),row.names = F,col.names = F,quote = F)

####################  CREATE A PSEUDO COLOR TABLE
cols <- col2rgb(c("black","white","green","orange","orange","lightblue","lightpink","darkred","purple"))
pct <- data.frame(cbind(c(0:6,21,22),
                        cols[1,],
                        cols[2,],
                        cols[3,]
))

write.table(pct,paste0(res_dir,"color_table_combine.txt"),row.names = F,col.names = F,quote = F)

####################  CREATE INDEX OF EXISTING TILES
system(sprintf("gdaltindex %s %s",
               paste0(res_dir,"index_results.shp"),
               paste0(ope_dir,"bfast_*/*2019.tif")
))

####################  MERGE AS HOMOGENEOUS BLOCKS
tiles <- readOGR(paste0(res_dir,"index_results.shp"))
tiles$fusion <- 1
u_tiles <- unionSpatialPolygons(tiles,tiles$fusion)

poly_list <- u_tiles@polygons[1][[1]]@Polygons

lp <- list()

for(i in 1:length(poly_list)){
  poly <- Polygons(list(poly_list[i][[1]]),i)
  lp <- append(lp,list(poly))
}

tiles_u <-SpatialPolygonsDataFrame(
  SpatialPolygons(lp,1:length(lp)), 
  data.frame(1:length(poly_list)), 
  match.ID = F
)
names(tiles_u) <- "fusion_id"
proj4string(tiles_u) <- proj4string(tiles)
tiles$fusion_id <- over(tiles,tiles_u)$fusion_id

table(tiles$fusion_id)

writeOGR(tiles,paste0(res_dir,"index_results.shp"),"index_results","ESRI Shapefile",overwrite_layer = T)

k <- 13

####################  LOOP THROUGH EACH BLOCK 

for(k in unique(tiles$fusion_id)[4:14]){
  
  blk_dir <- paste0(res_dir,"block",k,"/") 
  dir.create(blk_dir,showWarnings = F)
  
  to_merge <- paste0(tiles@data[tiles@data$fusion_id == k,"location"],collapse = " ")
  block    <- paste0(res_dir,"block_",k,".tif")
  date     <- paste0(blk_dir,"tmp_block_",k,"_date.tif")
  magn     <- paste0(blk_dir,"tmp_block_",k,"_magnitude.tif")
  
  ###################  MERGE TILES FOR THE BLOCK
  system(sprintf("gdal_merge.py -o %s -co COMPRESS=LZW -v %s",
                 block,
                 to_merge))
  
  ####################  EXTRACT DATE
  system(sprintf("gdal_translate -b 1 -co COMPRESS=LZW %s %s",
                 block,
                 date))
  
  ####################  EXTRACT MAGNITUDE
  system(sprintf("gdal_translate -b 2 -co COMPRESS=LZW %s %s",
                 block,
                 magn))
  
  ####################  COMPUTE  STATS FOR MAGNITUDE
  stats   <- paste0(blk_dir,"stats_block",k,".txt") 
  system(sprintf("gdalinfo -stats %s > %s",
                 magn,
                 stats
  ))
  
  s <- readLines(stats)
  maxs_b2   <- as.numeric(unlist(strsplit(s[grepl("STATISTICS_MAXIMUM",s)],"="))[2])
  mins_b2   <- as.numeric(unlist(strsplit(s[grepl("STATISTICS_MINIMUM",s)],"="))[2])
  means_b2  <- as.numeric(unlist(strsplit(s[grepl("STATISTICS_MEAN",s)],"="))[2])
  stdevs_b2 <- as.numeric(unlist(strsplit(s[grepl("STATISTICS_STDDEV",s)],"="))[2])
  
  # ####################  GET AOI LIMITS
  # system(sprintf("python %s/oft-rasterize_attr.py -v %s -i %s -o %s -a %s",
  #                scriptdir,
  #                paste0(ace_dir,"all_phu_indonesia_latlon.shp"),
  #                block,
  #                paste0(blk_dir,"aoi_block_",k,".tif"),
  #                "id"
  # ))
  
  # ####################  GET AOI LIMITS
  # system(sprintf("python %s/oft-rasterize_attr.py -v %s -i %s -o %s -a %s",
  #                scriptdir,
  #                paste0(res_dir,"index_results.shp"),
  #                block,
  #                paste0(blk_dir,"aoi_block_",k,".tif"),
  #                "fusion_id"
  # ))

  ####################  GET LCLU MAP
  system(sprintf("python %s/oft-rasterize_attr_int16.py -v %s -i %s -o %s -a %s",
                 scriptdir,
                 paste0(luc_dir,"PL_2015_2016__Planologi_KLHK__Rev_Topology.shp"),
                 block,
                 paste0(blk_dir,"lulc_block_",k,".tif"),
                 "PL16_ID"
  ))  
  
  ####################  GET LCLU MAP
  system(sprintf("python %s/oft-rasterize_attr_int16.py -v %s -i %s -o %s -a %s",
                 scriptdir,
                 paste0(ace_dir,"all_phu_indonesia_latlon.shp"),
                 block,
                 paste0(blk_dir,"phu_block_",k,".tif"),
                 "id"
  ))  
  
  
  # #############################################################
  # ### CROP AMPLITUDE TO AOI
  # system(sprintf("gdal_calc.py -A %s --A_band=2 -B %s --co=COMPRESS=LZW --overwrite --outfile=%s --calc=\"%s\"",
  #                result,
  #                paste0(aoi_dir,"aoi_block_",k,".tif"),
  #                paste0(blk_dir,"tmp",base,"_amplitude_crop.tif"),
  #                "(B>0)*A"
  # ))
  # 
  # #############################################################
  # ### CROP DATE TO AOI
  # system(sprintf("gdal_calc.py -A %s --A_band=1 -B %s --co=COMPRESS=LZW --overwrite --outfile=%s --calc=\"%s\"",
  #                result,
  #                paste0(blk_dir,"aoi.tif"),
  #                paste0(blk_dir,"tmp",base,"_date_crop.tif"),
  #                "(B>0)*A"
  # ))

  ####################  COMPUTE THRESHOLDS LAYER
  system(sprintf("gdal_calc.py -A %s -B %s -C %s --co=COMPRESS=LZW --type=Byte --overwrite --outfile=%s --calc=\"%s\"",
                 magn,
                 paste0(blk_dir,"phu_block_",k,".tif"),
                 paste0(blk_dir,"lulc_block_",k,".tif"),
                 paste0(blk_dir,"tmp_bfast_threshold_block_",k,".tif"),
                 #paste0("((B==0)+(A==0))*0+((A>0)+(A<0))*(B>0)*((B==5001)*1+((B<5001)+(B>5001))*(",
                 paste0("(B>0)*((C==5001)*1+((C>5001)+(C<5001))*(",
                        '(A<=',(maxs_b2),")*",
                        '(A>' ,(means_b2+(stdevs_b2*4)),")*9+",
                        '(A<=',(means_b2+(stdevs_b2*4)),")*",
                        '(A>' ,(means_b2+(stdevs_b2*3)),")*8+",
                        '(A<=',(means_b2+(stdevs_b2*3)),")*",
                        '(A>' ,(means_b2+(stdevs_b2*2)),")*7+",
                        '(A<=',(means_b2+(stdevs_b2*2)),")*",
                        '(A>' ,(means_b2+(stdevs_b2)),")*6+",
                        '(A<=',(means_b2+(stdevs_b2)),")*",
                        '(A>' ,(means_b2-(stdevs_b2)),")*1+",
                        '(A>=',(mins_b2),")*",
                        '(A<' ,(means_b2-(stdevs_b2*4)),")*5+",
                        '(A>=',(means_b2-(stdevs_b2*4)),")*",
                        '(A<' ,(means_b2-(stdevs_b2*3)),")*4+",
                        '(A>=',(means_b2-(stdevs_b2*3)),")*",
                        '(A<' ,(means_b2-(stdevs_b2*2)),")*3+",
                        '(A>=',(means_b2-(stdevs_b2*2)),")*",
                        '(A<' ,(means_b2-(stdevs_b2)),")*2))")
  ))
  
  
  
  ################################################################################
  ## Add pseudo color table to result
  system(sprintf("(echo %s) | oft-addpct.py %s %s",
                 paste0(res_dir,"color_table.txt"),
                 paste0(blk_dir,"tmp_bfast_threshold_block_",k,".tif"),
                 paste0(blk_dir,"tmp_bfast_threshold_block_pct_",k,".tif")
  ))
  
  ## Compress final result
  system(sprintf("gdal_translate -ot UInt16 -co COMPRESS=LZW %s %s",
                 paste0(blk_dir,"tmp_bfast_threshold_block_pct_",k,".tif"),
                 paste0(res_dir,"bfast_block_",k,"_threshold.tif")
  ))


  
  # ####################  COMBINE LULC WITH THRESHOLDS
  # system(sprintf("gdal_calc.py -A %s -B %s --co=COMPRESS=LZW --overwrite --outfile=%s --calc=\"%s\"",
  #                paste0(blk_dir,"lulc_block_",k,".tif"),
  #                paste0(res_dir,"bfast_block_",k,"_threshold.tif"),
  #                paste0(blk_dir,"tmp_bfast_combine_lulc_block_",k,".tif"),
  #                paste0("(B==1)*1 +(B>1)*(",
  #                       "(A==5001)*1+",
  #                       "(A==2010)*(B>1)*3 +",
  #                       "(A==2006)*(B>1)*4 +",
  #                       "((A==2001)+(A==2004)+(A==2005))*(B>1)*(B<6)*21 +",
  #                       "((A==2002)+(A==20041)+(A==20051))*(B>1)*(B<6)*22 +",
  #                       "((A==2001)+(A==2002)+(A==2004)+(A==2005)+(A==20041)+(A==20051))*(B>=6)*2 +",
  #                       "((A==2007)+(A==20071)+(A==2012)+(A==2014)+(A==3000)+(A>=20091))*(B>=6)*5 +",
  #                       "((A==2007)+(A==20071)+(A==2012)+(A==2014)+(A==3000)+(A>=20091))*(B<6)*6)"
  #                       )
  # ))

  # ####################  COMBINE LULC WITH THRESHOLDS
  # system(sprintf("gdal_calc.py -A %s -B %s --co=COMPRESS=LZW --overwrite --outfile=%s --calc=\"%s\"",
  #                paste0(blk_dir,"lulc_block_",k,".tif"),
  #                paste0(res_dir,"bfast_block_",k,"_threshold.tif"),
  #                paste0(blk_dir,"tmp_bfast_combine_lulc_block_",k,".tif"),
  #                paste0("((A==0)+(A==5001))*0+(A>0)*((A<5001)+(A>5001))*B")
  # ))
  # 
  # ################################################################################
  # ## Add pseudo color table to result
  # system(sprintf("(echo %s) | oft-addpct.py %s %s",
  #                paste0(res_dir,"color_table.txt"),
  #                paste0(blk_dir,"tmp_bfast_combine_lulc_block_",k,".tif"),
  #                paste0(blk_dir,"tmp_bfast_combine_lulc_block_pct_",k,".tif")
  # ))
  # 
  # ## Compress final result
  # system(sprintf("gdal_translate -ot Byte -co COMPRESS=LZW %s %s",
  #                paste0(blk_dir,"tmp_bfast_combine_lulc_block_pct_",k,".tif"),
  #                paste0(res_dir,"bfast_block_",k,"_threshold_combine.tif")
  # ))
  # 
  ## Year of disturbance
  system(sprintf("gdal_calc.py -A %s --co=COMPRESS=LZW --type=UInt16 --overwrite --outfile=%s --calc=\"%s\"",
                 date,
                 paste0(blk_dir,"bfast_block_",k,"_year.tif"),
                 "floor(A)"
  ))
  
  ## Day of disturbance
  system(sprintf("gdal_calc.py -A %s --co=COMPRESS=LZW --type=UInt16 --overwrite --outfile=%s --calc=\"%s\"",
                 date,
                 paste0(blk_dir,"bfast_block_",k,"_day.tif"),
                 "(A-floor(A))*365"
  ))
  
  system(sprintf("rm -r -f %s",
                 paste0(blk_dir,"tmp*.tif")))
}
  