list_res <- list.files(ts_dir,pattern = glob2rx("bfast*.tif"),recursive = T)

for(file in list_res){
  base <- strsplit(file,split = "/")[[1]][4]
  tile <- strsplit(file,split = "/")[[1]][1]
  file.copy(paste0(ts_dir,file),paste0(bfst_dir,username,"_tile_",tile,"_",base))
}