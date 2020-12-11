# MS2Purifier
# This R script reads in the mzXML files of DDA & AIF data, outputs csv files containing the metabolic feature table with purified MS2 spectra.
# Shipei Xing, Dec 4, 2020
# Copyright @ The University of British Columbia

####################
# Note: 1. DDA and AIF mzXML files should be in the same directory
#       2. DDA and AIF files of the same sample should be named as follows:
#                         "XXXX_DDA.mzXML"     "XXXX_AIF.mzXML" (only differ by 'DDA' and 'AIF')

####################
# Libraries
library(xcms)
library(randomForest)
message('Libraries loaded.')

# Directories
file_directory <- 'E:/MS2Purifier_2020/bruker_data'       # directory containing DDA and AIF mzXML files
model_directory <- 'E:/MS2Purifier_2020/R_scripts&model'  # directory containing machine learning models

# Parameter settings
vendor <- 1     # select the vendor index: 
                # Bruker: 1; Agilent: 2; Thermo: 3     
                # due to file conversion problems, Bruker AIF mslevel=1, Agilent scan number issues
                # Bruker file conversion: DataAnalysis; Agilent & Thermo file conversion: MSConvert
analysis_time <- 30  # LC analysis time (in min), features with RT larger than this limit will be discarded
mz_tol <- 0.01       # m/z tolerance for MS2 assignment.
rt_tol <- 30         # retention time tolerance (in sec) for MS2 assignment.
mz_EIC_tol <- 0.01   # m/z tolerance window for EIC generation
rt_EIC_tol <- 15     # retention time tolerance (in sec), for EIC generation, rt window is [rt-rt_EIC_tol,rt+rt_EIC_tol]
peak_height <- T     # logical, TRUE for peak height, FALSE for peak area
multi_classification <- FALSE # logical, TRUE for conducting multi-label classification


# Functions
peak_smooth <- function(x,level=2){
  n <- level
  if(length(x) < 2*n){
    return(x)
  } else if(length(unique(x))==1){
    return(x)
  } else{
    y <- vector(length=length(x))
    for(i in 1:n){
      y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
    }
    for(i in (n+1):(length(y)-n)){
      y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
    }
    for(i in (length(y)-n+1):length(y)){
      y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
    }
    return(y)
  }
}
ideal_slope <- function(x){
  if(length(unique(x))==1){return(0)}
  if(length(unique(x)) > 1){
    max_index <- which(x %in% max(x))
    max_index <- max_index[which.min(abs(max_index-length(x)/2))]
    dx <- vector(length = length(x)) #first derivative of x, dx=(-2x[-2]-x[-1]+x[+1]+2x[+2])/10
    dx[c(1,length(x))] <- 0
    dx[2] <- (x[3] - x[1])/2
    dx[length(x)-1] <- (x[length(x)] - x[length(x)-2])/2
    for(i in 3:(length(x)-2)){dx[i] <- (-2*x[i-2]-x[i-1]+x[i+1]+2*x[i+2])/10}
    left <- 0
    right <- 0
    if(max_index!=1){
      for(i in 1:(max_index-1)){ if(dx[i]>=0){ left <- left + dx[i]} }
    }
    if(max_index!=length(x)){
      for(i in (max_index+1):(length(x)-1)){if(dx[i]<=0){ right <- right + abs(dx[i])} }
    }
    ideal_slope_ratio <- (left+right)/sum(abs(dx[-max_index]))
    return(ideal_slope_ratio)
  }
}
sharpness <- function(x){
  if(length(unique(x))==1){return(0)}
  if(length(unique(x))>1){
    max_index <- which(x %in% max(x))
    max_index <- max_index[which.min(abs(max_index-length(x)/2))]
    left <- vector(mode = 'numeric')
    right <- vector(mode = 'numeric')
    if(max_index!=1){
      for(i in 1:(max_index-1)){ left[i] <- (max(x)-x[i])/((max_index-i)*sqrt(max(x))) }
    }
    if(max_index!=length(x)){
      for(i in (max_index+1):length(x)){ right[i-max_index] <- (max(x)-x[i])/((i-max_index)*sqrt(max(x))) }
    }
    if(max_index==1){ 
      return(max(right))
    }else if(max_index==length(x)){ 
      return(max(left))
    }else{return(0.5*(max(left)+max(right)))}
  }
}

# Read in machine learning models
setwd(model_directory)
rf_binary_model <- readRDS('RF_model_binary.rds')
rf_multi_model <- readRDS('RF_model_multi.rds')
message('Machine learning models loaded.')

###### Main program ########
message('Main program starts.')
time_0 <- Sys.time()
# DDA mzXML files
setwd(file_directory)
DDAfiles <- list.files(pattern = 'DDA.mzXML')

# XCMS
xset <- xcmsSet(DDAfiles,
                method = "centWave",
                ppm=12,
                peakwidth=c(10,120),
                mzdiff = 0.01,
                snthresh = 6,
                integrate = 1,
                prefilter = c(3,100),
                noise = 100
)
if(length(DDAfiles) > 1){
  xset <- group(xset, bw = 5, minfrac = 0.5, mzwid = 0.015, minsamp = 1, max = 50)
  xset <- retcor(xset, method = "obiwarp", profStep = 1)
  xset <- group(xset, bw = 5, minfrac = 0.5, mzwid = 0.015, minsamp = 1, max = 50)
  xset <- fillPeaks(xset)
  XCMt <- data.frame(xset@groups)
  if(peak_height==T){
    xcmI <- groupval(xset, value = "maxo")
  } else {xcmI <- groupval(xset, value = "into")}
  featureTable <- as.data.frame(cbind(0,XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, XCMt$rtmax, xcmI))
}

if(length(DDAfiles) == 1){
  if(peak_height==T){featureTable <- as.data.frame(cbind(0,xset@peaks))[,c(1,2,5:7,9)]}
  if(peak_height==F){featureTable <- as.data.frame(cbind(0,xset@peaks))[,c(1,2,5:8)]}
}
colnames(featureTable)[1:5] <- c("featureID","mz", "RT", "RTmin", "RTmax")
# output DDA feature table
featureTable <- featureTable[featureTable$RT < (60*analysis_time),]
featureTable$featureID <- 1:nrow(featureTable)
#write.csv(featureTable, file = "Feature_table_DDA.csv",row.names = F)
message('DDA data processing completed.')

# add 3 columns: ms2_file, ms2mz, ms2int
featureTable$ms2_file <- NA
featureTable$ms2mz <- NA
featureTable$ms2int <- NA

# Bruker
if(vendor == 1){
  # read in DDA & AIF mzXML files
  message('Loading all the mzXML files...')
  for(m in 1:length(DDAfiles)){
    DDAfile.name <- DDAfiles[m]
    AIFfile.name <- paste0(strsplit(DDAfiles[m],'DDA.mzXML')[[1]],'AIF.mzXML')
    message(paste0('File loading: ',DDAfile.name, ' & ',AIFfile.name))
    
    # read AIF MS1 & MS2 data, DDA MS2 data
    # Bruker AIF, msLevel=1 for all spectra
    AIF_file <- readMSData(AIFfile.name, msLevel. = 1, centroided. = T)
    DDA_file <- readMSData(DDAfile.name, msLevel. = 2, centroided. = T)
    assign(paste0('AIF_file_',m), AIF_file)
    assign(paste0('DDA_file_',m), DDA_file)
    
    # AIF MS1 rt
    AIF_ms1rt <- vector(mode='numeric')

    for(i in seq(1,length(AIF_file),2)){AIF_ms1rt <- c(AIF_ms1rt, AIF_file[[i]]@rt)}
    assign(paste0('AIF_ms1rt_',m),AIF_ms1rt)

    # scan index vector
    assign(paste0('AIF_ms1scanindex_',m), seq(1,length(AIF_file),2))
    assign(paste0('AIF_ms2scanindex_',m), seq(2,length(AIF_file),2))
    
    # DDA
    DDA_rt <- vector(mode='numeric')
    DDA_premz <- vector(mode='numeric')
    DDA_peakCount <- vector(mode='numeric')
    DDA_preint <- vector(mode='numeric')

    for(i in 1:length(DDA_file)){
      DDA_rt <- c(DDA_rt,DDA_file[[i]]@rt)
      DDA_premz <- c(DDA_premz,DDA_file[[i]]@precursorMz)
      DDA_peakCount <- c(DDA_peakCount,DDA_file[[i]]@peaksCount)
      DDA_preint <- c(DDA_preint,DDA_file[[i]]@precursorIntensity)
    }
    assign(paste0('DDA_rt_',m),DDA_rt)
    assign(paste0('DDA_premz_',m),DDA_premz)
    assign(paste0('DDA_peakCount_',m),DDA_peakCount)
    assign(paste0('DDA_preint_',m),DDA_preint)
  }
  message('All mzXML files loaded.  :)')
  
  message('MS2 assignment for each feature...')
  # for each feature, find the MS2 reference file
  for(i in 1:nrow(featureTable)){
    # sort by peak int.
    fileindex <- order(featureTable[i,6:(5+length(DDAfiles))],decreasing = T)
    # target mz & rt
    t_mz <- featureTable$mz[i]
    t_rt <- featureTable$RT[i]
    for(m in 1:length(fileindex)){ # go through files in order (high to low int)
      k <- fileindex[m]
      DDA_rt <- get(paste0('DDA_rt_',k))
      DDA_premz <- get(paste0('DDA_premz_',k))
      DDA_preint <- get(paste0('DDA_preint_',k))
      # for each feature, whether there is MS2 available
      DDA_rtindex <- which(DDA_rt %in% DDA_rt[DDA_rt>(t_rt-rt_tol) & DDA_rt<(t_rt+rt_tol)])
      selected_mz <- DDA_premz[DDA_rtindex]
      index_h <- which(selected_mz %in% selected_mz[selected_mz<(t_mz+mz_tol) & selected_mz>(t_mz-mz_tol)])
      DDA_index <- DDA_rtindex[index_h]
      if(length(DDA_index)==0) next
      
      # find the MS2 scan with highest precursor intensity for each feature in DDA
      DDA_peakCount <- get(paste0('DDA_peakCount_',k))
      pre_int <- vector(mode = "numeric")
      for(p in 1:length(DDA_index)) {
        if(DDA_peakCount[DDA_index[p]]!=0){
          pre_int <- c(pre_int,DDA_preint[DDA_index[p]])
        } else{ pre_int <- c(pre_int,0) }
      }
      if(max(pre_int)==0) next 
      
      # original MS2 mz&int strings, filled in feature table
      DDA_file <- get(paste0('DDA_file_',k))
      DDA_index_ms2 <- DDA_index[which(pre_int %in% max(pre_int))[1]]
      featureTable$ms2_file[i] <- k
      featureTable$ms2mz[i] <- paste(round(DDA_file[[DDA_index_ms2]]@mz,4),collapse = ';')
      featureTable$ms2int[i] <- paste(DDA_file[[DDA_index_ms2]]@intensity,collapse = ';')
      break
    }
  }
  output_ft <- featureTable # ms2_file: index to file name
  for(i in 1:nrow(output_ft)){output_ft$ms2_file[i] <- DDAfiles[as.numeric(featureTable$ms2_file[i])]}
  message('MS2 assignment completed. :)')
  write.csv(output_ft, file = 'Feature_table_DDA_MS2assigned.csv', row.names = F)
  
  
  # feature table with MS2 assigned
  DDA_ft <- featureTable[complete.cases(featureTable),]
  # table for RF prediction
  rf_df <- data.frame(matrix(ncol = 10))
  colnames(rf_df) <- c('featureID','DDA_preMz','rt','DDA_fragMz','DDA_fragint','PPC_smoothlv2',
                       'ideal_slope_pre','sharpness_pre','ideal_slope_frag','sharpness_frag')
  rf_df_individual <- rf_df
  
  message('Generating the data matrix for machine learning prediction...This step takes some time.')
  d <- 0
  for(i in 1:nrow(DDA_ft)) { #for each metabolic feature, use AIF data to generate data matrix for RF
    # feature mz & RT
    t_mz <- DDA_ft$mz[i]
    t_rt <- DDA_ft$RT[i]
    # fragment mz vector & intensity vector
    fragMz_v <- as.numeric(strsplit(DDA_ft$ms2mz[i],';')[[1]])
    fragint_v <- as.numeric(strsplit(DDA_ft$ms2int[i],';')[[1]])
    
    # AIF
    AIF_ms1rt <- get(paste0('AIF_ms1rt_',DDA_ft$ms2_file[i]))
    AIF_file <- get(paste0('AIF_file_',DDA_ft$ms2_file[i]))
    AIF_ms1scanindex <- get(paste0('AIF_ms1scanindex_',DDA_ft$ms2_file[i]))
    AIF_ms2scanindex <- get(paste0('AIF_ms2scanindex_',DDA_ft$ms2_file[i]))
    
    # find corresponding AIF data 
    AIF_rt_index <- which(AIF_ms1rt %in% AIF_ms1rt[abs(AIF_ms1rt-t_rt) < rt_EIC_tol])
    AIF_pre_index <- AIF_ms1scanindex[AIF_rt_index]
    AIF_pre_int <- vector(mode = 'numeric',length = length(AIF_pre_index))
    for(m in 1:length(AIF_rt_index)){
      if(min(abs(AIF_file[[AIF_pre_index[m]]]@mz-t_mz)) <= mz_EIC_tol){
        mz_index <- which(abs(AIF_file[[AIF_pre_index[m]]]@mz-t_mz) %in% min(abs(AIF_file[[AIF_pre_index[m]]]@mz-t_mz)))
        AIF_pre_int[m] <- AIF_file[[AIF_pre_index[m]]]@intensity[mz_index]
      }
    }
    
    rf_df_individual$featureID <- DDA_ft$featureID[i]
    rf_df_individual$rt <- t_rt
    rf_df_individual$DDA_preMz <- t_mz
    rf_df_individual$ideal_slope_pre <- ideal_slope(AIF_pre_int)
    rf_df_individual$sharpness_pre <- sharpness(AIF_pre_int)
    
    # for each fragment ion in the MS2 scan
    AIF_frag_index <- AIF_ms2scanindex[AIF_rt_index]
    
    AIF_frag_int <- vector(mode = 'numeric',length = length(AIF_rt_index))
    for(l in 1:length(fragMz_v)) { # fragment EIC
      f_mz <- fragMz_v[l]
      for(m in 1:length(AIF_rt_index)){
        if(min(abs(AIF_file[[AIF_frag_index[m]]]@mz-f_mz)) <= mz_EIC_tol){
          mz_index <- which(abs(AIF_file[[AIF_frag_index[m]]]@mz-f_mz) %in% min(abs(AIF_file[[AIF_frag_index[m]]]@mz-f_mz)))
          AIF_frag_int[m] <- AIF_file[[AIF_frag_index[m]]]@intensity[mz_index]
        } else{AIF_frag_int[m] <- 0}
      }
      rf_df_individual$DDA_fragint <- fragint_v[l]
      rf_df_individual$DDA_fragMz <- f_mz
      
      if(length(unique(AIF_frag_int))==1){
        rf_df_individual$PPC_smoothlv2 <- 0
        rf_df_individual$ideal_slope_frag <- 0
        rf_df_individual$sharpness_frag <- 0
      } else{
        rf_df_individual$PPC_smoothlv2 <- cor(peak_smooth(AIF_pre_int,level=2),peak_smooth(AIF_frag_int,level=2),method = "pearson")
        rf_df_individual$ideal_slope_frag <- ideal_slope(AIF_frag_int)
        rf_df_individual$sharpness_frag <- sharpness(AIF_frag_int)
      }
      if(length(unique(AIF_pre_int))==1){rf_df_individual$PPC_smoothlv2 <- 0}
      rf_df <- rbind(rf_df,rf_df_individual)
    }
    if((10*i) %/% nrow(DDA_ft) > d){
      d <- ((10*i) %/% nrow(DDA_ft))
      message(paste0(10*d,'% completed.'))
    }
  }

  rf_df <- rf_df[complete.cases(rf_df),]
  
  message('Applying the machine learning model...')
  # Apply the random forest model
  pred <- predict(rf_binary_model, newdata=rf_df[,5:10])
  rf_df_pred <- rf_df
  rf_df_pred$prediction <- pred
  
  # in case that precursor EIC is missing in AIF data, keep all the fragments
  rf_df_pred$prediction[rf_df_pred$ideal_slope_frag==0 & rf_df_pred$sharpness_frag==0] <- as.factor('F')
  
  write.csv(rf_df_pred, file = 'Fragment_prediction_result.csv',row.names = F)
  
  if(multi_classification==T){
    pred_multi <- predict(rf_multi_model, newdata=rf_df[,5:10])
    rf_df_pred_multi <- rf_df
    rf_df_pred_multi$prediction <- as.character(pred_multi)
    rf_df_pred_multi$prediction[rf_df_pred_multi$ideal_slope_frag==0 &
                                  rf_df_pred_multi$sharpness_frag==0] <- as.character('Artifact_ion')
    write.csv(rf_df_pred_multi, file = 'Fragment_prediction_multi_result.csv',row.names = F)
  }
  
  # feature table with purified MS2
  message('Removing MS2 contaminant ions...')
  output_pft <- output_ft
  output_pft$purified_ms2mz <- NA
  output_pft$purified_ms2int <- NA
  for(i in 1:length(unique(rf_df_pred$featureID))){
    id <- unique(rf_df_pred$featureID)[i]
    
    p_frag_v <- rf_df_pred$DDA_fragMz[rf_df_pred$featureID == id & rf_df_pred$prediction==1]
    if(length(p_frag_v) == 0) next
    
    p_frag <- paste0(p_frag_v, collapse = ';')
    output_pft$purified_ms2mz[output_pft$featureID==id] <- p_frag
    
    p_fragint <- paste0(rf_df_pred$DDA_fragint[rf_df_pred$featureID == id & 
                                                 rf_df_pred$prediction==1], collapse = ';')
    output_pft$purified_ms2int[output_pft$featureID==id] <- p_fragint
  }
  write.csv(output_pft, file = 'Feature_table_DDA_MS2purified.csv',row.names = F)
  message('MS2 purification completed. :)')
  print(Sys.time() - time_0)
}

# Agilent
if(vendor == 2){
  # read in DDA & AIF mzXML files
  message('Loading all the mzXML files...')
  for(m in 1:length(DDAfiles)){
    DDAfile.name <- DDAfiles[m]
    AIFfile.name <- paste0(strsplit(DDAfiles[m],'DDA.mzXML')[[1]],'AIF.mzXML')
    message(paste0('File loading: ',DDAfile.name, ' & ',AIFfile.name))
    
    # read AIF MS1 & MS2 data, DDA MS2 data
    # Agilent AIF, read in all MS spectra (MS1 & MS2)
    AIF_filems1 <- readMSData(AIFfile.name, msLevel. = 1, centroided. = T)
    AIF_filems2 <- readMSData(AIFfile.name, msLevel. = 2, centroided. = T)
    DDA_file <- readMSData(DDAfile.name, msLevel. = 2, centroided. = T)
    assign(paste0('AIF_filems1_',m), AIF_filems1)
    assign(paste0('AIF_filems2_',m), AIF_filems2)
    assign(paste0('DDA_file_',m), DDA_file)
    
    # AIF MS1
    AIF_ms1rt <- vector(mode='numeric')

    for(i in 1:length(AIF_filems1)){AIF_ms1rt <- c(AIF_ms1rt, AIF_filems1[[i]]@rt)}
    assign(paste0('AIF_ms1rt_',m),AIF_ms1rt)

    
    # DDA
    DDA_rt <- vector(mode='numeric')
    DDA_premz <- vector(mode='numeric')
    DDA_peakCount <- vector(mode='numeric')
    DDA_preint <- vector(mode='numeric')
    
    for(i in 1:length(DDA_file)){
      DDA_rt <- c(DDA_rt,DDA_file[[i]]@rt)
      DDA_premz <- c(DDA_premz,DDA_file[[i]]@precursorMz)
      DDA_peakCount <- c(DDA_peakCount,DDA_file[[i]]@peaksCount)
      DDA_preint <- c(DDA_preint,DDA_file[[i]]@precursorIntensity)
    }
    assign(paste0('DDA_rt_',m),DDA_rt)
    assign(paste0('DDA_premz_',m),DDA_premz)
    assign(paste0('DDA_peakCount_',m),DDA_peakCount)
    assign(paste0('DDA_preint_',m),DDA_preint)
  }
  message('All mzXML files loaded.  :)')
  
  message('MS2 assignment for each feature...')
  # for each feature, find the MS2 reference file
  for(i in 1:nrow(featureTable)){
    # sort by peak int.
    fileindex <- order(featureTable[i,6:(5+length(DDAfiles))],decreasing = T)
    # target mz & rt
    t_mz <- featureTable$mz[i]
    t_rt <- featureTable$RT[i]
    for(m in 1:length(fileindex)){ # go through files in order (high to low int)
      k <- fileindex[m]
      DDA_rt <- get(paste0('DDA_rt_',k))
      DDA_premz <- get(paste0('DDA_premz_',k))
      DDA_preint <- get(paste0('DDA_preint_',k))
      # for each feature, whether there is MS2 available
      DDA_rtindex <- which(DDA_rt %in% DDA_rt[DDA_rt>(t_rt-rt_tol) & DDA_rt<(t_rt+rt_tol)])
      selected_mz <- DDA_premz[DDA_rtindex]
      index_h <- which(selected_mz %in% selected_mz[selected_mz<(t_mz+mz_tol) & selected_mz>(t_mz-mz_tol)])
      DDA_index <- DDA_rtindex[index_h]
      if(length(DDA_index)==0) next
      
      # find the MS2 scan with highest precursor intensity for each feature in DDA
      DDA_peakCount <- get(paste0('DDA_peakCount_',k))
      pre_int <- vector(mode = "numeric")
      for(p in 1:length(DDA_index)) {
        if(DDA_peakCount[DDA_index[p]]!=0){
          pre_int <- c(pre_int,DDA_preint[DDA_index[p]])
        } else{ pre_int <- c(pre_int,0) }
      }
      if(max(pre_int)==0) next 
      
      # original MS2 mz&int strings, filled in feature table
      DDA_file <- get(paste0('DDA_file_',k))
      DDA_index_ms2 <- DDA_index[which(pre_int %in% max(pre_int))[1]]
      featureTable$ms2_file[i] <- k
      featureTable$ms2mz[i] <- paste(round(DDA_file[[DDA_index_ms2]]@mz,4),collapse = ';')
      featureTable$ms2int[i] <- paste(DDA_file[[DDA_index_ms2]]@intensity,collapse = ';')
      break
    }
  }
  output_ft <- featureTable # ms2_file: index to file name
  for(i in 1:nrow(output_ft)){output_ft$ms2_file[i] <- DDAfiles[as.numeric(featureTable$ms2_file[i])]}
  message('MS2 assignment completed. :)')
  write.csv(output_ft, file = 'Feature_table_DDA_MS2assigned.csv', row.names = F)
  
  
  # feature table with MS2 assigned
  DDA_ft <- featureTable[complete.cases(featureTable),]
  # table for RF prediction
  rf_df <- data.frame(matrix(ncol = 10))
  colnames(rf_df) <- c('featureID','DDA_preMz','rt','DDA_fragMz','DDA_fragint','PPC_smoothlv2',
                       'ideal_slope_pre','sharpness_pre','ideal_slope_frag','sharpness_frag')
  rf_df_individual <- rf_df
  
  message('Generating the data matrix for machine learning prediction...This step takes some time.')
  d <- 0
  for(i in 1:nrow(DDA_ft)) { #for each metabolic feature, use AIF data to generate data matrix for RF
    # feature mz & RT
    t_mz <- DDA_ft$mz[i]
    t_rt <- DDA_ft$RT[i]
    # fragment mz vector & intensity vector
    fragMz_v <- as.numeric(strsplit(DDA_ft$ms2mz[i],';')[[1]])
    fragint_v <- as.numeric(strsplit(DDA_ft$ms2int[i],';')[[1]])
    
    # AIF
    AIF_ms1rt <- get(paste0('AIF_ms1rt_',DDA_ft$ms2_file[i]))
    AIF_filems1 <- get(paste0('AIF_filems1_',DDA_ft$ms2_file[i]))
    AIF_filems2 <- get(paste0('AIF_filems2_',DDA_ft$ms2_file[i]))
    
    
    # find corresponding AIF data 
    AIF_rt_index <- which(AIF_ms1rt %in% AIF_ms1rt[abs(AIF_ms1rt-t_rt) < rt_EIC_tol])
    AIF_rt_index <- AIF_rt_index[AIF_rt_index <= length(AIF_filems2)]
    AIF_pre_int <- vector(mode = 'numeric',length = length(AIF_rt_index))
    for(m in 1:length(AIF_rt_index)){
      if(min(abs(AIF_filems1[[AIF_rt_index[m]]]@mz-t_mz)) <= mz_EIC_tol){
        mz_index <- which(abs(AIF_filems1[[AIF_rt_index[m]]]@mz-t_mz) %in% min(abs(AIF_filems1[[AIF_rt_index[m]]]@mz-t_mz)))
        AIF_pre_int[m] <- AIF_filems1[[AIF_rt_index[m]]]@intensity[mz_index]
      }
    }
    
    rf_df_individual$featureID <- DDA_ft$featureID[i]
    rf_df_individual$rt <- t_rt
    rf_df_individual$DDA_preMz <- t_mz
    rf_df_individual$ideal_slope_pre <- ideal_slope(AIF_pre_int)
    rf_df_individual$sharpness_pre <- sharpness(AIF_pre_int)
    
    
    AIF_frag_int <- vector(mode = 'numeric',length = length(AIF_rt_index))
    for(l in 1:length(fragMz_v)) { # fragment EIC
      f_mz <- fragMz_v[l]
      for(m in 1:length(AIF_rt_index)){
        if(min(abs(AIF_filems2[[AIF_rt_index[m]]]@mz-f_mz)) <= mz_EIC_tol){
          mz_index <- which(abs(AIF_filems2[[AIF_rt_index[m]]]@mz-f_mz) %in% min(abs(AIF_filems2[[AIF_rt_index[m]]]@mz-f_mz)))
          AIF_frag_int[m] <- AIF_filems2[[AIF_rt_index[m]]]@intensity[mz_index]
        } else{AIF_frag_int[m] <- 0}
      }
      rf_df_individual$DDA_fragint <- fragint_v[l]
      rf_df_individual$DDA_fragMz <- f_mz
      
      if(length(unique(AIF_frag_int))==1){
        rf_df_individual$PPC_smoothlv2 <- 0
        rf_df_individual$ideal_slope_frag <- 0
        rf_df_individual$sharpness_frag <- 0
      } else{
        rf_df_individual$PPC_smoothlv2 <- cor(peak_smooth(AIF_pre_int,level=2),peak_smooth(AIF_frag_int,level=2),method = "pearson")
        rf_df_individual$ideal_slope_frag <- ideal_slope(AIF_frag_int)
        rf_df_individual$sharpness_frag <- sharpness(AIF_frag_int)
      }
      if(length(unique(AIF_pre_int))==1){rf_df_individual$PPC_smoothlv2 <- 0}
      rf_df <- rbind(rf_df,rf_df_individual)
    }
    if((10*i) %/% nrow(DDA_ft) > d){
      d <- ((10*i) %/% nrow(DDA_ft))
      message(paste0(10*d,'% completed.'))
    }
  }
  
  rf_df <- rf_df[complete.cases(rf_df),]
  
  message('Applying the machine learning model...')
  # Apply the random forest model
  pred <- predict(rf_binary_model, newdata=rf_df[,5:10])
  rf_df_pred <- rf_df
  rf_df_pred$prediction <- pred
  
  # in case that precursor EIC is missing in AIF data, keep all the fragments
  rf_df_pred$prediction[rf_df_pred$ideal_slope_frag==0 & rf_df_pred$sharpness_frag==0] <- as.factor('F')
  
  write.csv(rf_df_pred, file = 'Fragment_prediction_result.csv',row.names = F)
  
  if(multi_classification==T){
    pred_multi <- predict(rf_multi_model, newdata=rf_df[,5:10])
    rf_df_pred_multi <- rf_df
    rf_df_pred_multi$prediction <- as.character(pred_multi)
    rf_df_pred_multi$prediction[rf_df_pred_multi$ideal_slope_frag==0 &
                                  rf_df_pred_multi$sharpness_frag==0] <- as.character('Artifact_ion')
    write.csv(rf_df_pred_multi, file = 'Fragment_prediction_multi_result.csv',row.names = F)
  }
  
  # feature table with purified MS2
  message('Removing MS2 contaminant ions...')
  output_pft <- output_ft
  output_pft$purified_ms2mz <- NA
  output_pft$purified_ms2int <- NA
  for(i in 1:length(unique(rf_df_pred$featureID))){
    id <- unique(rf_df_pred$featureID)[i]
    
    p_frag_v <- rf_df_pred$DDA_fragMz[rf_df_pred$featureID == id & rf_df_pred$prediction==1]
    if(length(p_frag_v) == 0) next
    
    p_frag <- paste0(p_frag_v, collapse = ';')
    output_pft$purified_ms2mz[output_pft$featureID==id] <- p_frag
    
    p_fragint <- paste0(rf_df_pred$DDA_fragint[rf_df_pred$featureID == id & 
                                                 rf_df_pred$prediction==1], collapse = ';')
    output_pft$purified_ms2int[output_pft$featureID==id] <- p_fragint
  }
  write.csv(output_pft, file = 'Feature_table_DDA_MS2purified.csv',row.names = F)
  message('MS2 purification completed. :)')
  print(Sys.time() - time_0)
}

# Thermo
if(vendor == 3){
  # read in DDA & AIF mzXML files
  message('Loading all the mzXML files...')
  for(m in 1:length(DDAfiles)){
    DDAfile.name <- DDAfiles[m]
    AIFfile.name <- paste0(strsplit(DDAfiles[m],'DDA.mzXML')[[1]],'AIF.mzXML')
    message(paste0('File loading: ',DDAfile.name, ' & ',AIFfile.name))
    
    # read AIF MS1 & MS2 data, DDA MS2 data
    # Thermo AIF, read in all MS spectra (MS1 & MS2)
    AIF_filems1 <- readMSData(AIFfile.name, msLevel. = 1, centroided. = T)
    AIF_filems2 <- readMSData(AIFfile.name, msLevel. = 2, centroided. = T)
    DDA_file <- readMSData(DDAfile.name, msLevel. = 2, centroided. = T)
    assign(paste0('AIF_filems1_',m), AIF_filems1)
    assign(paste0('AIF_filems2_',m), AIF_filems2)
    assign(paste0('DDA_file_',m), DDA_file)
    
    # AIF MS1
    AIF_ms1rt <- vector(mode='numeric')

    for(i in 1:length(AIF_filems1)){AIF_ms1rt <- c(AIF_ms1rt, AIF_filems1[[i]]@rt)}
    assign(paste0('AIF_ms1rt_',m),AIF_ms1rt)

    
    # DDA
    DDA_rt <- vector(mode='numeric')
    DDA_premz <- vector(mode='numeric')
    DDA_peakCount <- vector(mode='numeric')
    DDA_preint <- vector(mode='numeric')
    
    for(i in 1:length(DDA_file)){
      DDA_rt <- c(DDA_rt,DDA_file[[i]]@rt)
      DDA_premz <- c(DDA_premz,DDA_file[[i]]@precursorMz)
      DDA_peakCount <- c(DDA_peakCount,DDA_file[[i]]@peaksCount)
      DDA_preint <- c(DDA_preint,DDA_file[[i]]@precursorIntensity)
    }
    assign(paste0('DDA_rt_',m),DDA_rt)
    assign(paste0('DDA_premz_',m),DDA_premz)
    assign(paste0('DDA_peakCount_',m),DDA_peakCount)
    assign(paste0('DDA_preint_',m),DDA_preint)
  }
  message('All mzXML files loaded.  :)')
  
  message('MS2 assignment for each feature...')
  # for each feature, find the MS2 reference file
  for(i in 1:nrow(featureTable)){
    # sort by peak int.
    fileindex <- order(featureTable[i,6:(5+length(DDAfiles))],decreasing = T)
    # target mz & rt
    t_mz <- featureTable$mz[i]
    t_rt <- featureTable$RT[i]
    for(m in 1:length(fileindex)){ # go through files in order (high to low int)
      k <- fileindex[m]
      DDA_rt <- get(paste0('DDA_rt_',k))
      DDA_premz <- get(paste0('DDA_premz_',k))
      DDA_preint <- get(paste0('DDA_preint_',k))
      # for each feature, whether there is MS2 available
      DDA_rtindex <- which(DDA_rt %in% DDA_rt[DDA_rt>(t_rt-rt_tol) & DDA_rt<(t_rt+rt_tol)])
      selected_mz <- DDA_premz[DDA_rtindex]
      index_h <- which(selected_mz %in% selected_mz[selected_mz<(t_mz+mz_tol) & selected_mz>(t_mz-mz_tol)])
      DDA_index <- DDA_rtindex[index_h]
      if(length(DDA_index)==0) next
      
      # find the MS2 scan with highest precursor intensity for each feature in DDA
      DDA_peakCount <- get(paste0('DDA_peakCount_',k))
      pre_int <- vector(mode = "numeric")
      for(p in 1:length(DDA_index)) {
        if(DDA_peakCount[DDA_index[p]]!=0){
          pre_int <- c(pre_int,DDA_preint[DDA_index[p]])
        } else{ pre_int <- c(pre_int,0) }
      }
      if(max(pre_int)==0) next 
      
      # original MS2 mz&int strings, filled in feature table
      DDA_file <- get(paste0('DDA_file_',k))
      DDA_index_ms2 <- DDA_index[which(pre_int %in% max(pre_int))[1]]
      featureTable$ms2_file[i] <- k
      featureTable$ms2mz[i] <- paste(round(DDA_file[[DDA_index_ms2]]@mz,4),collapse = ';')
      featureTable$ms2int[i] <- paste(DDA_file[[DDA_index_ms2]]@intensity,collapse = ';')
      break
    }
  }
  output_ft <- featureTable # ms2_file: index to file name
  for(i in 1:nrow(output_ft)){output_ft$ms2_file[i] <- DDAfiles[as.numeric(featureTable$ms2_file[i])]}
  message('MS2 assignment completed. :)')
  write.csv(output_ft, file = 'Feature_table_DDA_MS2assigned.csv', row.names = F)
  
  
  # feature table with MS2 assigned
  DDA_ft <- featureTable[complete.cases(featureTable),]
  # table for RF prediction
  rf_df <- data.frame(matrix(ncol = 10))
  colnames(rf_df) <- c('featureID','DDA_preMz','rt','DDA_fragMz','DDA_fragint','PPC_smoothlv2',
                       'ideal_slope_pre','sharpness_pre','ideal_slope_frag','sharpness_frag')
  rf_df_individual <- rf_df
  
  message('Generating the data matrix for machine learning prediction...This step takes some time.')
  d <- 0
  for(i in 1:nrow(DDA_ft)) { #for each metabolic feature, use AIF data to generate data matrix for RF
    # feature mz & RT
    t_mz <- DDA_ft$mz[i]
    t_rt <- DDA_ft$RT[i]
    # fragment mz vector & intensity vector
    fragMz_v <- as.numeric(strsplit(DDA_ft$ms2mz[i],';')[[1]])
    fragint_v <- as.numeric(strsplit(DDA_ft$ms2int[i],';')[[1]])
    
    # AIF
    AIF_ms1rt <- get(paste0('AIF_ms1rt_',DDA_ft$ms2_file[i]))
    AIF_filems1 <- get(paste0('AIF_filems1_',DDA_ft$ms2_file[i]))
    AIF_filems2 <- get(paste0('AIF_filems2_',DDA_ft$ms2_file[i]))

    
    # find corresponding AIF data 
    AIF_rt_index <- which(AIF_ms1rt %in% AIF_ms1rt[abs(AIF_ms1rt-t_rt) < rt_EIC_tol])
    AIF_rt_index <- AIF_rt_index[AIF_rt_index <= length(AIF_filems2)]
    AIF_pre_int <- vector(mode = 'numeric',length = length(AIF_rt_index))
    for(m in 1:length(AIF_rt_index)){
      if(min(abs(AIF_filems1[[AIF_rt_index[m]]]@mz-t_mz)) <= mz_EIC_tol){
        mz_index <- which(abs(AIF_filems1[[AIF_rt_index[m]]]@mz-t_mz) %in% min(abs(AIF_filems1[[AIF_rt_index[m]]]@mz-t_mz)))
        AIF_pre_int[m] <- AIF_filems1[[AIF_rt_index[m]]]@intensity[mz_index]
      }
    }
    
    rf_df_individual$featureID <- DDA_ft$featureID[i]
    rf_df_individual$rt <- t_rt
    rf_df_individual$DDA_preMz <- t_mz
    rf_df_individual$ideal_slope_pre <- ideal_slope(AIF_pre_int)
    rf_df_individual$sharpness_pre <- sharpness(AIF_pre_int)
    

    AIF_frag_int <- vector(mode = 'numeric',length = length(AIF_rt_index))
    for(l in 1:length(fragMz_v)) { # fragment EIC
      f_mz <- fragMz_v[l]
      for(m in 1:length(AIF_rt_index)){
        if(min(abs(AIF_filems2[[AIF_rt_index[m]]]@mz-f_mz)) <= mz_EIC_tol){
          mz_index <- which(abs(AIF_filems2[[AIF_rt_index[m]]]@mz-f_mz) %in% min(abs(AIF_filems2[[AIF_rt_index[m]]]@mz-f_mz)))
          AIF_frag_int[m] <- AIF_filems2[[AIF_rt_index[m]]]@intensity[mz_index]
        } else{AIF_frag_int[m] <- 0}
      }
      rf_df_individual$DDA_fragint <- fragint_v[l]
      rf_df_individual$DDA_fragMz <- f_mz
      
      if(length(unique(AIF_frag_int))==1){
        rf_df_individual$PPC_smoothlv2 <- 0
        rf_df_individual$ideal_slope_frag <- 0
        rf_df_individual$sharpness_frag <- 0
      } else{
        rf_df_individual$PPC_smoothlv2 <- cor(peak_smooth(AIF_pre_int,level=2),peak_smooth(AIF_frag_int,level=2),method = "pearson")
        rf_df_individual$ideal_slope_frag <- ideal_slope(AIF_frag_int)
        rf_df_individual$sharpness_frag <- sharpness(AIF_frag_int)
      }
      if(length(unique(AIF_pre_int))==1){rf_df_individual$PPC_smoothlv2 <- 0}
      rf_df <- rbind(rf_df,rf_df_individual)
    }
    if((10*i) %/% nrow(DDA_ft) > d){
      d <- ((10*i) %/% nrow(DDA_ft))
      message(paste0(10*d,'% completed.'))
    }
  }
  
  rf_df <- rf_df[complete.cases(rf_df),]
  
  message('Applying the machine learning model...')
  # Apply the random forest model
  pred <- predict(rf_binary_model, newdata=rf_df[,5:10])
  rf_df_pred <- rf_df
  rf_df_pred$prediction <- pred
  
  # in case that precursor EIC is missing in AIF data, keep all the fragments
  rf_df_pred$prediction[rf_df_pred$ideal_slope_frag==0 & rf_df_pred$sharpness_frag==0] <- as.factor('F')
  
  write.csv(rf_df_pred, file = 'Fragment_prediction_result.csv',row.names = F)
  
  if(multi_classification==T){
    pred_multi <- predict(rf_multi_model, newdata=rf_df[,5:10])
    rf_df_pred_multi <- rf_df
    rf_df_pred_multi$prediction <- as.character(pred_multi)
    rf_df_pred_multi$prediction[rf_df_pred_multi$ideal_slope_frag==0 &
                                  rf_df_pred_multi$sharpness_frag==0] <- as.character('Artifact_ion')
    write.csv(rf_df_pred_multi, file = 'Fragment_prediction_multi_result.csv',row.names = F)
  }
  
  # feature table with purified MS2
  message('Removing MS2 contaminant ions...')
  output_pft <- output_ft
  output_pft$purified_ms2mz <- NA
  output_pft$purified_ms2int <- NA
  for(i in 1:length(unique(rf_df_pred$featureID))){
    id <- unique(rf_df_pred$featureID)[i]
    
    p_frag_v <- rf_df_pred$DDA_fragMz[rf_df_pred$featureID == id & rf_df_pred$prediction==1]
    if(length(p_frag_v) == 0) next
    
    p_frag <- paste0(p_frag_v, collapse = ';')
    output_pft$purified_ms2mz[output_pft$featureID==id] <- p_frag
    
    p_fragint <- paste0(rf_df_pred$DDA_fragint[rf_df_pred$featureID == id & 
                                                 rf_df_pred$prediction==1], collapse = ';')
    output_pft$purified_ms2int[output_pft$featureID==id] <- p_fragint
  }
  write.csv(output_pft, file = 'Feature_table_DDA_MS2purified.csv',row.names = F)
  message('MS2 purification completed. :)')
  print(Sys.time() - time_0)
}

