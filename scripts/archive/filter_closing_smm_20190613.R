source('~/ws_idn_20190819/scripts/s0_parameters.R')

res_dir <- paste0(smm_dir,"extra/")

postprocess_pysmm <- function(tile_dir){
  
  tile_dir <- paste0(tile_dir,"/")
  base <- "SMCmap_"
  
  list <- list.files(tile_dir,pattern = glob2rx(paste0(base,"*.tif")))
  file <- list[1]
  
  
  for(file in list){
    
    date <- substr(file,nchar(base)+1,nchar(file)-4)
    
    system(sprintf("otbcli_GrayScaleMorphologicalOperation -in %s -out %s -structype %s -structype.ball.xradius %s -structype.ball.yradius %s -filter closing",
                   paste0(tile_dir,file),             
                   paste0(tile_dir,"tmp_closing.tif"),
                   "ball",
                   1,
                   1))
    
    system(sprintf("gdal_calc.py -A %s -B %s --co=COMPRESS=LZW --type=UInt16 --overwrite --outfile=%s --calc=\"%s\"",
                   paste0(tile_dir,file),
                   paste0(tile_dir,"tmp_closing.tif"),
                   paste0(tile_dir,"close_",file),
                   "(A==0)*B+A"
    ))
    
    system(sprintf("rm -r -f %s",
                   paste0(tile_dir,"tmp_closing.tif")
    ))
    
  }
  
  base <- "close_SMCmap_"
  list <- list.files(tile_dir,pattern = glob2rx(paste0(base,"*.tif")))
  
  system(sprintf("gdalbuildvrt -separate %s %s",
                 paste0(tile_dir,"ts_smm.vrt"),
                 paste0(tile_dir,list,collapse = " ")
  ))
  
  r <- brick(paste0(tile_dir,"ts_smm.vrt"))
  
  fn2 <- substr(list,nchar(base)+1,nchar(list)-4)
  
  Date <- as.Date(fn2,format='%Y_%m_%d')
  
  fun_slope <- function(y) { 
    if(all(is.na(y))) {
      NA
    } else {
      m = lm(y ~ Date); summary(m)$coefficients[2] 
    }
  }
  ## and this one is to calculate the p-value
  fun_pvalue <- function(y) { 
    if(all(is.na(y))) {
      NA
    } else {
      m = lm(y ~ Date); summary(m)$coefficients[8] 
    }
  }
  
  slope <- calc(r, fun_slope)
  pvalue <- calc(r,fun_pvalue)
  
  plot(slope)
  plot(pvalue)
  
  writeRaster(slope,paste0(tile_dir,"slope.tif"),overwrite=T)
  writeRaster(pvalue,paste0(tile_dir,"pvalue.tif"),overwrite=T)
  
}

for(tile_dir in list.dirs(res_dir)[-1]){
  postprocess_pysmm(tile_dir)}

the_folder_you_want <- "/home/waluyo/smm_maps/example/"
postprocess_pysmm(the_folder_you_want)


