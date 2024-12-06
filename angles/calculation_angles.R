require(dplyr)
require(stringr)
require(ggplot2)

rm(list=ls())

current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
dirn <- paste0(current_dir, "/")

################################################ BODY ANGLE #########################################
## Calculate the body angle alpha (between point 1 [snout tip] and point 10 [insertion pelvic fin])
#In the publication, point 10 has been renamed point 4 (see Figure 1)
files_2Dcoords <- list.files(dirn, pattern="2Dcoords_lateral")

body_angle <- data.frame(matrix(ncol = 3, nrow = 0))

for (i in 1:length(files_2Dcoords)) {
  dataset_2Dcoords <- read.csv(paste0(dirn,files_2Dcoords[i]), header=T)
  #Get segment length opposed to alpha
  dOpp <- dataset_2Dcoords$pt10_cam2_Y-dataset_2Dcoords$pt1_cam2_Y
  #Get segment length adjacent to alpha
  dAdj <- dataset_2Dcoords$pt10_cam2_X-dataset_2Dcoords$pt1_cam2_X
  #Get alpha, in degrees (degrees = radians*180/pi)
  alpha <- atan(dOpp/dAdj)*180/pi
  #Save alpha values with details on specimens and speed
  details_video <- str_split(files_2Dcoords[i],"_")
  temp <- as.data.frame(cbind(details_video[[1]][3], 
                              substr(details_video[[1]][4],1,nchar(details_video[[1]][4])-2), 
                              alpha))
  body_angle <- rbind(body_angle, temp)
}

colnames(body_angle) <- c("specimen", "speedBL", "alpha")

write.table(body_angle, paste0(dirn, "body_angle.txt"), sep=",", dec= ".", row.names=F, col.names=T)


################################################ ANGLE OF ATTACK #######################################

##Calculate the angle of attack (how pectoral fins are used): 
#3D angle between pectoral-to-pelvic fin insertion (body baseline) 
#and pect fin insertion-pect fin tip (pect extension);
#with pect fin insertion=pt6, pelvic fin insertion=pt10, pect fin tip=pt8
#In the publication, point 6 has been renamed as point 2, 
# point 8 has been renamed as point 3, and point 10 has been renamed as point 4 (see Figure 1)

files_3Dcoords <- list.files(dirn, pattern="3Dcoords")

angle_attack <- data.frame(matrix(ncol = 3, nrow = 0))

for (i in 1:length(files_3Dcoords)) {
  dataset_3Dcoords <- read.csv(paste0(dirn,files_3Dcoords[i]), header=T)
  
  #Calculate vectors between points 6 and 10
  body_baseline <- cbind(dataset_3Dcoords$pt10_X-dataset_3Dcoords$pt6_X,
                         dataset_3Dcoords$pt10_Y-dataset_3Dcoords$pt6_Y,
                         dataset_3Dcoords$pt10_Z-dataset_3Dcoords$pt6_Z)
  #and 6 and 8
  pect_extension <- cbind(dataset_3Dcoords$pt8_X-dataset_3Dcoords$pt6_X,
                          dataset_3Dcoords$pt8_Y-dataset_3Dcoords$pt6_Y,
                          dataset_3Dcoords$pt8_Z-dataset_3Dcoords$pt6_Z)
  #Calculate the dot product of the two vectors
  baseline.extension <- rowSums(body_baseline*pect_extension)
  #Calculate the magnitude of vectors (absolute value)
  magnitude_baseline <- sqrt(rowSums(body_baseline^2))
  magnitude_extension <- sqrt(rowSums(pect_extension^2))
  
  #Get the angle of attack
  angle <- (acos(baseline.extension/(magnitude_baseline*magnitude_extension)))*180/pi
  
  #Add to dataset
  details_video <- str_split(files_3Dcoords[i],"_")
  temp <- as.data.frame(cbind(details_video[[1]][2], 
                              substr(details_video[[1]][3],1,nchar(details_video[[1]][3])-2), 
                              angle))
  angle_attack <- rbind(angle_attack, temp)
}

colnames(angle_attack) <- c("specimen", "speedBL", "angle_attack")

write.table(angle_attack, paste0(dirn, "angle_attack.txt"), sep=",", dec= ".", row.names=F, col.names=T)
