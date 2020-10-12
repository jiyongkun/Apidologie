genus = "Apis"
packages_needed <- c("raster", # for raster analysis
                     "dismo", # a collection of ENM/SDM tools
                     "rgeos","rgdal","sp", # spatial data analysis
                     "ENMeval", # a few new tools in ENM/SDM
                     "wallace",   # interface for Maxent and ENMeval
                     "utils", # for zip & unzip files
                     "jsonlite", # necessary for download data from GBIF
                     "rJava"
                     )
pk_to_install <- packages_needed [!( packages_needed %in% rownames(installed.packages())  )]
if(length(pk_to_install)>0 ){
  install.packages(pk_to_install)
}

library("raster")
library("dismo")
library("rgeos")
library("rgdal")
library("sp")
library("ENMeval")
library("rJava")
library("data.table")
library("devtools")
if(!require(devtools)){
    install.packages("devtools")
}
if(!require(kuenm)){
    devtools::install_github("marlonecobos/kuenm")
library("kuenm")

# download maxent.jar 3.3.3k, and place the file in the desired folder; note that, there may be a newer version of Maxent
if( !file.exists(paste0(system.file("java", package="dismo"),"/maxent.jar"))  )   {
utils::download.file(url="https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar",
                     destfile=paste0(system.file("java", package="dismo"),"/maxent.jar"),
                     mode="wb") ## wb for binary file, otherwise maxent.jar can not execute
}
# also note that both R and Java need to be the same bit (either 32 or 64) to be compatible to run

# to increase memory size of the JVM and prevent memory issues with Maxent.jar
options( java.parameters = c("-Xss2560k", "-Xmx4g") ) 

# load rasters
clim_list_current <- list.files("data/bioclim/current/",pattern=".tif$",full.names = T)#加$符号指最后的字符串
clim_current <- raster::stack(clim_list_current) 
clim_list_LGM <- list.files("data/bioclim/LGM/",pattern = ".tif$",full.names = T)
clim_LGM <- raster::stack(clim_list_LGM)
clim_list_future <- list.files("data/bioclim/future/",pattern = ".tif$",full.names = T)
clim_future <- raster::stack(clim_list_future)

#evaluate correlation between bioclimate layers
correlation=layerStats(clim_current_crop, 'pearson', na.rm=T)
corr_matrix=correlation$'pearson correlation coefficient'
write.csv(corr_matrix,file="output/bioclim_var_correlation.csv")

bio_select = c("02","04","05","06","08","15","16","19")
clim_current <- clim_current[[paste("bio",bio_select,sep = "")]]
clim_LGM <- clim_LGM[[paste("bio",bio_select,sep = "")]]
mybox <- extent(-18,   155,     -36,    71)
clim_current_crop <- crop(clim_current, mybox)
clim_LGM_crop <- crop(clim_LGM, mybox)
clim_future_crop <- crop(clim_future, mybox)

#read occurrence records of each species
csv_list <- list.files(paste0("data/occdata/",genus,"/"),pattern=".csv$",full.names = T)#加$符号指最后的字符串
for (acsv in csv_list) {
  occ_final = read.csv(file=acsv,header = T,stringsAsFactors = FALSE)
  occ_final_COPY <- occ_final
  occ_final <- subset(occ_final,(!is.na(decimalLatitude))&(!is.na(decimalLongitude)))
  species_name = as.character(occ_final$species[1])
  cat(species_name,"\n")
  coordinates(occ_final) <- ~ decimalLongitude + decimalLatitude
  myCRS1 <- CRS("+init=epsg:4326") # WGS 84
  crs(occ_final) <- myCRS1
  conditions_occ <- extract(clim_current_crop,occ_final)
  bad_records <- is.na(conditions_occ[,1] ) 
  occ_final <- occ_final[!bad_records,]
  occ_buffer <- buffer(occ_final,width=4*10^5) #unit is meter
  clim_mask <- mask(clim_current_crop, occ_buffer)
  # extract environmental conditions
  set.seed(1) 
  bg <- sampleRandom(x=clim_mask,
                     size=10000,
                     na.rm=T, #removes the 'Not Applicable' points  
                     sp=T) # return spatial points 
  
  set.seed(1) 
  # randomly select for training
  if (nrow(occ_final)>100) {
    fraction = 0.75
    selected <- sample(  1:nrow(occ_final),  nrow(occ_final)*fraction)
    occ_train <- occ_final[selected,] # this is the selection to be used for model training
    occ_test <- occ_final[-selected,] # this is the opposite of the selection which will be used for model testing
    
  } else {
    fraction = 1
    occ_train <- occ_final # this is the selection to be used for model training
    occ_test <- occ_final # this is the opposite of the selection which will be used for model testing
    
  }
  
  # extracting env conditions
  env_occ_train <- extract(clim_current_crop,occ_train)
  env_occ_test <- extract(clim_current_crop,occ_test)
  # extracting env conditions for background
  
  env_bg <- extract(clim_current_crop,bg)  
  
  #combine the conditions by row
  myPredictors <- rbind(env_occ_train,env_bg)
  
  # change matrix to dataframe
  myPredictors <- as.data.frame(myPredictors)
  
  # repeat the number 1 as many times as the number of rows in p, 
  # and repeat 0 for the rows of background points
  myResponse <- c(rep(1,nrow(env_occ_train)),
                  rep(0,nrow(env_bg))) 
  
  env <- clim_mask[[paste("bio",c("2","4","5","6","15","16","19"),sep = "")]]
  occ_coord <- occ_final@coords
  bg_coord <- bg@coords
  
  #selecte the best model parameter
  res_bl <- ENMevaluate(occ_coord, env, bg_coord, RMvalues=seq(0.5,4,0.5),method='block',
                        parallel= TRUE,
                        numCores = 4,
                        updateProgress = TRUE)
  write.csv(res_bl@results,file=paste0("output/",genus,"/",species_name,"_competition_result.csv"))
  aic_df = res_bl@results
  aic_best=aic_df[which(aic_df$delta.AICc == 0),] 
  features = as.character(aic_best$features[1])
  features = trimws(features, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  rm_value = aic_best$rm[1]
  rm_value = trimws(rm_value, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  rm_value = as.numeric(rm_value)
  cat(features,"\n")
  cat("rm:",rm_value,"\n")
  
  source("code/setPara.R")
  myparameters <- prepPara(userfeatures=features,betamultiplier=rm_value)
  mod <- maxent(x=myPredictors,
                p=myResponse,
                #path=paste0(getwd(),"/output/Apis/mod"),
                args=myparameters )
  
  
  ped_current <- predict(clim_current_crop,model=mod,progress='text') # studyArea is the clipped rasters 
  writeRaster(ped_current,format="EHdr",file=paste0("output/",genus,"/",species_name,"_current.bil"),overwrite=TRUE)
  
  ped_LGM <- predict(clim_LGM_crop,model=mod,progress='text') # studyArea is the clipped rasters 
  writeRaster(ped_LGM,format="EHdr",file=paste0("output/",genus,"/",species_name,"_LGM.bil"),overwrite=TRUE)
  
  ped_future <- predict(clim_future_crop,model=mod,progress='text') # studyArea is the clipped rasters 
  writeRaster(ped_future,format="EHdr",file=paste0("output/",genus,"/",species_name,"_future.bil"),overwrite=TRUE)
  
  sink(paste0("output/",genus,"/",species_name,"_AUC.txt"))
  mod_eval_train <- dismo::evaluate(p=env_occ_train,
                                    a=env_bg,
                                    model=mod) 
  cat("Train:","\n")
  print(mod_eval_train)
  mod_eval_test <- dismo::evaluate(p=env_occ_test,
                                    a=env_bg,
                                    model=mod)
  cat("Test:","\n")
  print(mod_eval_test)
  sink()
}
