####################################################################################
####### Object: Google Drive to Local Drive    
####### Author:  remi.dannunzio@fao.org                               
####### Update:  2017/10/22                                    
###################################################################################

####### TIME SERIES DATA ARE GENERATED IN GOOGLE EARTH ENGINE
##      https://code.earthengine.google.com/6349290af151862c244cac3bcdc44318


###################################################################################
#### Parameters
###################################################################################

#### Root directory

####################################################################################################################
####### LOAD AUTHORIZATION KEY FOR "DRIVE" AND DOWNLOAD RESULTS
####################################################################################################################
source('~/ws_idn_20190819/scripts/s0_parameters.R')

options(echo=TRUE)
args <- commandArgs(TRUE)
print(args[1])
auth_key <- args[1]
#### Select a basename for the archives to transfer
base <- 'smcmap_'
setwd(smm_dir)

#### Initialize the DRIVE function, change the authorization key
system(sprintf("echo %s | drive init",
               auth_key))

#### Read list of files in GEDrive that contain base
system(sprintf("drive list -matches %s > %s",
               paste0(base),
               "list_down.txt"))

data_input <- basename(unlist(read.table("list_down.txt")))
data_input    
#### download
for(data in data_input){
  system(sprintf("drive pull %s",
                 data))
}

lapply(data_input,function(x){file.rename(x,paste0(smm_dir,x))})

