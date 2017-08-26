#' To calculate the c3 Comobidity Index (targeted at cancer patients)
#' @author Miao Cai <email: miao.cai@slu.edu>
#' @description This file aims to calculate the c3 Comobidity Index
#' @param data Your data file in which Elixhauser Comorbidity Index is to be calculated
#' @param comorbidity A vector of all comorbidity variables
#' @return data_file: a new data.frame named "data_file". This data frame contains a new variable "c3": the c3 Index was developed by Sarfati et al. in 2014. This is a comorbidity index targeted at cancer patients.
#' @export

# The c3 index refers to:
#　Sarfati, D., et al., Cancer-specific administrative data-based comorbidity indices provided valid alternative to Charlson and National Cancer Institute Indices. J Clin Epidemiol, 2014. 67(5): p. 586-95

c3 <- function(data, comorbidity) {
  ### PART A: SUBSETTING DATA---------------------------
  suppressMessages(library(dplyr))
  data_comorbidity <- subset(data_file, select = comorbidity)# subset the whole dataset into "data_comorbidity": subset data that only contains comorbidities
  data_comorbidity %>% mutate_if(is.factor, as.character) -> data_comorbidity # this converts factor values into characters
  data_comorbidity[is.na(data_comorbidity)] <- 0 #imputate missing values with zeros
  dim_comorbidity <- dim(data_comorbidity)# save the dimensionality of comorbidity subset data
  unlisted_data <- unlist(data_comorbidity)# unlist the comorbidity subset data(converts the data into a vector. Vectorization speeds up the functions)
  
  ### PART B: FILTERING COMORBIDITIES---------------------------
  ##1 Angina: Angina -- weight: 0.51
  bidata_com <- (substr(unlisted_data, 1, 3) == "I20")  # regular expression to find out cases with "Angina" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Angina <- rowSums(bidata_com)
  data_file$Angina[data_file$Angina >= 1] <- 1
  
  #2 inflam_bd: Inflammatory bowel disease -- Weight: 0.52 
  bidata_com <- (substr(unlisted_data,1,5) %in% c("K52.2","K52.8","K52.9")) | (substr(unlisted_data, 1, 3) %in% c("K50","K51")) # regular expression to find out cases with "Inflammatory bowel disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$inflam_bd <- rowSums(bidata_com)
  data_file$inflam_bd[data_file$inflam_bd >= 1] <- 1
  
  #3 Anemia: Anemia -- Weight: 0.59
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("D50","D51","D52","D53")) # regular expression to find out cases with "Anemia" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Anemia <- rowSums(bidata_com)
  data_file$Anemia[data_file$Anemia >= 1] <- 1
  
  #4 Metabolic: Metabolic conditions -- Weight: 0.61
  bidata_com <- (substr(unlisted_data,1,5) %in% c("E79.1","E79.8","E79.9")) | (substr(unlisted_data, 1, 3) %in% c("E70","E71","E72","E74","E75","E76","E77","E78","E80","E83","E85","E88")) # regular expression to find out cases with "Metabolic conditions" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Metabolic <- rowSums(bidata_com)
  data_file$Metabolic[data_file$Metabolic >= 1] <- 1
  
  #5 Othercard: Other cardiac conditions -- Weight: 0.62
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I24.8","I24.9","I25.0","I25.1","I25.3","I25.4","I25.6","I25.8","I25.9","I31.0","I31.1","I42.1","I42.2","I42.4")) # regular expression to find out cases with "Other cardiac conditions" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Othercard <- rowSums(bidata_com)
  data_file$Othercard[data_file$Othercard >= 1] <- 1
  
  #6 eye: (Major) Eye conditions affecting vision -- Weight: 0.63
  bidata_com <- (substr(unlisted_data,1,5) %in% c("H18.1","H18.4","H18.5","H18.6","H20.1","H21.2","H30.1","H31.1","H31.2","H31.3","H31.4","H33.0","H33.2","H33.3","H33.4","H33.5","H53.0","H53.1","H53.2","H53.3","H53.4","H53.6","H53.8","H53.9")) | (substr(unlisted_data, 1, 3) %in% c("H16","H34","H35","H43","H46","H47","H49","H50","H51","H54","Q12","Q13","Q14","Q15")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$eye <- rowSums(bidata_com)
  data_file$eye[data_file$eye >= 1] <- 1
  
  #7 Hypertension: Hypertension -- Weight: 0.72
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I11.9","I12.9","I13.9")) | (substr(unlisted_data, 1, 3) %in% c("I10")) # regular expression to find out cases with "Hypertension" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Hypertension <- rowSums(bidata_com)
  data_file$Hypertension[data_file$Hypertension >= 1] <- 1
  
  #8 Coagulopathies: Coagulopathies and other blood disorders -- Weight: 0.75
  bidata_com <- (substr(unlisted_data,1,5) %in% c("D59.0","D59.1","D59.2","D59.3","D59.4","D59.8","D59.9","D68.0","D68.1","D68.2","D68.8","D68.9","D69.1","D69.2","D69.3","D69.4","D69.6","D69.8","D69.9","D75.0","D75.2","D75.8","D75.9")) | (substr(unlisted_data, 1, 3) %in% c("D55","D56","D57","D58","D60","D61","D64","D64","D66","D67","D70","D71","D72","D74")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Coagulopathies <- rowSums(bidata_com)
  data_file$Coagulopathies[data_file$Coagulopathies >= 1] <- 1
  
  #9 arrhythmia: Cardiac arrhythmia -- Weight: 0.77
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I44.1","I44.2","I44.3","I45.6","I45.9","T82.1","Z45.0","Z95.0")) | (substr(unlisted_data, 1, 3) %in% c("I47","I48","I49")) # regular expression to find out cases with "Cardiac arrhythmia" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$arrhythmia <- rowSums(bidata_com)
  data_file$arrhythmia[data_file$arrhythmia >= 1] <- 1
  
  #10 Obesity: Obesity -- Weight: 0.83
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("E66")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Obesity <- rowSums(bidata_com)
  data_file$Obesity[data_file$Obesity >= 1] <- 1
  
  #11 Diabetes: Diabetes with complications -- Weight: 0.88
  bidata_com <- (substr(unlisted_data,1,5) %in% c(paste("E10.",2:8,sep = ""),paste("E11.",2:8,sep = ""),paste("E12.",2:8,sep = ""),paste("E13.",2:8,sep = ""),paste("E14.",2:8,sep = ""))) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Diabetes <- rowSums(bidata_com)
  data_file$Diabetes[data_file$Diabetes >= 1] <- 1
  
  #12 Liver: Liver disease -- Weight: 0.92
  bidata_com <- (substr(unlisted_data,1,5) %in% c("K71.1","K71.3","K71.4","K71.5","K71.7","K72.1","K72.9","K76.0","K76.2","K76.3","K76.4","K76.5","K76.6","K76.7","K76.8","K76.9","I86.4","I98.2","Z94.4")) | (substr(unlisted_data, 1, 3) %in% c("K70","K73","K74","I85")) # regular expression to find out cases with "Liver disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Liver <- rowSums(bidata_com)
  data_file$Liver[data_file$Liver >= 1] <- 1
  
  #13 MI: Previous MI -- Weight: 0.93
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I24.1","I25.2")) | (substr(unlisted_data, 1, 3) %in% c("I21","I22","I23")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$MI <- rowSums(bidata_com)
  data_file$MI[data_file$MI >= 1] <- 1
  
  #14 PVD: Peripheral vascular disease -- Weight: 0.98
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I73.1","I73.8","I73.9","I77.1","K55.1","K55.2","K55.8","K55.9")) | (substr(unlisted_data, 1, 3) %in% c("I70","I71","I72")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$PVD <- rowSums(bidata_com)
  data_file$PVD[data_file$PVD >= 1] <- 1
  
  #15 Cerebrovascular: Cerebrovascular disease -- Weight: 1.09
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("I60","I61","I62","I63","I64","I65","I66","I67","I69","G45","G46")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Cerebrovascular <- rowSums(bidata_com)
  data_file$Cerebrovascular[data_file$Cerebrovascular >= 1] <- 1
  
  #16 COPD: COPD and asthma -- Weight: 1.09
  bidata_com <- (substr(unlisted_data,1,5) %in% c("J68.4","J70.1","J70.3","J96.1","J98.0","J98.2","J98.3","J98.4")) | (substr(unlisted_data, 1, 3) %in% c("E84",paste("J",40:47,sep = ""),paste("J",60:67,sep = ""),"J84")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$COPD <- rowSums(bidata_com)
  data_file$COPD[data_file$COPD >= 1] <- 1
  
  #17 valve: Cardiac valve disorders -- Weight: 1.10
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I09.1","I09.8","T82.0","Q23.0","Q23.1","Q23.2","Q23.3","Q23.8","Q23.9","Z95.2","Z95.3","Z95.4")) | (substr(unlisted_data, 1, 3) %in% c("I05","I06","I07","I08","I34","I35","I36","I37","I38")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$valve <- rowSums(bidata_com)
  data_file$valve[data_file$valve >= 1] <- 1
  
  #18 CHF: Congestive heart failure -- Weight: 1.26
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I09.9","I11.0","I13.0","I13.2","I25.5","I42.0","I42.5","I42.6","I42.7","I42.8","I42.9")) | (substr(unlisted_data, 1, 3) %in% c("I43","I50")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$CHF <- rowSums(bidata_com)
  data_file$CHF[data_file$CHF >= 1] <- 1
  
  #19 Renal: Renal disease -- Weight: 1.38
  bidata_com <- (substr(unlisted_data,1,5) %in% c(paste("N03.",2:9,sep=""),paste("N04.",2:9,sep=""),paste("N05.",2:9,sep=""),"N25.0","N25.8","N25.9","I12.0","I13.1","Z94.0","Z99.2")) | (substr(unlisted_data, 1, 3) %in% c("N11","N18","N19","Z49")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Renal <- rowSums(bidata_com)
  data_file$Renal[data_file$Renal >= 1] <- 1
  
  #20 Alcohol: Alcohol abuse -- Weight: 1.08
  bidata_com <- (substr(unlisted_data,1,5) %in% c(paste("F10.",1:9,sep=""),"Z50.2","Z71.4")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Alcohol <- rowSums(bidata_com)
  data_file$Alcohol[data_file$Alcohol >= 1] <- 1
  
  #21 Anxiety: Anxiety/behavioral disorders -- Weight: 0.57
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("F40","F41","F42","F44","F45","F48","F50","F55","F59","F60","F61","F63","F64","F65","F66","F68","F69")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Anxiety <- rowSums(bidata_com)
  data_file$Anxiety[data_file$Anxiety >= 1] <- 1
  
  #22 Connective: Connective tissue disorders -- Weight: 0.51
  bidata_com <- (substr(unlisted_data,1,5) %in% c("M12.0","M12.3","N25.9","I12.0",paste("M35.",0:6,sep=""),"M35.8","M35.9")) | (substr(unlisted_data, 1, 3) %in% c("L93","M05","M06","M08",paste("M3",0:4,sep=""))) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Connective <- rowSums(bidata_com)
  data_file$Connective[data_file$Connective >= 1] <- 1
  
  #23 Dementia: Dementia -- Weight: 1.35
  bidata_com <- (substr(unlisted_data,1,5) %in% c("F02.0","F02.1","F02.2","F02.3","F05.1","G31.0","G31.1")) | (substr(unlisted_data, 1, 3) %in% c("F00","F01","F03","G30")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Dementia <- rowSums(bidata_com)
  data_file$Dementia[data_file$Dementia >= 1] <- 1
  
  #24 DiabetesOUT: Diabetes without complications -- Weight: -0.03
  bidata_com <- (substr(unlisted_data,1,5) %in% c("E10.0","E10.1","E10.9","E11.0","E11.1","E11.9","E12.0","E12.1","E12.9","E13.0","E13.1","E13.9","E14.0","E14.1","E14.9")) # regular expression to find out cases with "Diabetes without complications" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$DiabetesOUT <- rowSums(bidata_com)
  data_file$DiabetesOUT[data_file$DiabetesOUT >= 1] <- 1
  
  #25 Endocrine: Endocrine disorders -- Weight: 0.77
  bidata_com <- (substr(unlisted_data,1,5) %in% c("E06.2","E06.3","E06.5","E16.3","E16.4","E16.8","E16.9","E21.0","E21.2","E21.3","E21.4","E21.5","E23.0","E23.2","E23.3","E23.6","E23.7","E24.0","E24.1","E24.3","E24.4","E24.8","E24.9","E34.5","E34.8","E34.9")) | (substr(unlisted_data, 1, 3) %in% c("E01","E02","E03","E05","E07","E20","E22","E25","E26","E27","E31","E32")) # regular expression to find out cases with "Endocrine disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Endocrine <- rowSums(bidata_com)
  data_file$Endocrine[data_file$Endocrine >= 1] <- 1
  
  #26 uppergi: Upper GI disorders -- Weight: 0.11
  bidata_com <- (substr(unlisted_data,1,5) %in% c("K22.0","K22.1","K22.4","K22.5","K22.8","K22.9","K31.1","K31.2","K31.4","K31.6")) | (substr(unlisted_data, 1, 3) %in% c("K25","K26","K27","K28")) # regular expression to find out cases with "Upper GI disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$uppergi <- rowSums(bidata_com)
  data_file$uppergi[data_file$uppergi >= 1] <- 1
  
  #27 InnerEar: Inner ear disorders -- Weight: 0.54
  bidata_com <- (substr(unlisted_data,1,5) %in% c("H91.1","H91.0","H91.3","H91.8","H91.9","H93.0","H93.1","H93.2","H93.3")) | (substr(unlisted_data, 1, 3) %in% c("H80","H81","H83","H90")) # regular expression to find out cases with "Inner ear disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$InnerEar <- rowSums(bidata_com)
  data_file$InnerEar[data_file$InnerEar >= 1] <- 1
  
  #28 Intestinal: Intestinal disorders -- Weight: 0.11
  bidata_com <- (substr(unlisted_data,1,5) %in% c("K59.2","K59.3")) | (substr(unlisted_data, 1, 3) %in% c("K57","K90")) # regular expression to find out cases with "Intestinal disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Intestinal <- rowSums(bidata_com)
  data_file$Intestinal[data_file$Intestinal >= 1] <- 1
  
  #29 spinal: Joint and spinal disorders -- Weight: 0.69
  bidata_com <- (substr(unlisted_data,1,5) %in% c("M15.0","M15.1","M15.2","M15.4","M15.8","M15.9","M40.0","M40.2","M40.3","M40.4","M40.5","M46.0","M46.1","M46.2","M48.0","M48.1","M48.2","M48.5","M48.8","M48.9","G95.0","G95.1")) | (substr(unlisted_data, 1, 3) %in% c("M41","M42","M43","M45","M47")) # regular expression to find out cases with "Joint and spinal disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$spinal <- rowSums(bidata_com)
  data_file$spinal[data_file$spinal >= 1] <- 1
  
  #30 psychiatric: Major psychiatric disorders -- Weight: 0.79
  bidata_com <- (substr(unlisted_data,1,5) %in% c("F30.2","F32.1","F32.2","F32.3","F32.8","F32.9")) | (substr(unlisted_data, 1, 3) %in% c("F20","F22","F25","F28","F29","F31","F33","F39")) # regular expression to find out cases with "Major psychiatric disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$psychiatric <- rowSums(bidata_com)
  data_file$psychiatric[data_file$psychiatric >= 1] <- 1
  
  #31 Nutritional: Nutritional disorders -- Weight: 1.16
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("E40","E41","E42","E43","E44","E45","E46","E50","E51","E52","E53","E54","E55","E56","E58","E59","E60","E61","E63","E64")) # regular expression to find out cases with "Nutritional disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Nutritional <- rowSums(bidata_com)
  data_file$Nutritional[data_file$Nutritional >= 1] <- 1
  
  #32 Neurological: Neurological disorders excluding epilepsy -- Weight: 1.06
  bidata_com <- (substr(unlisted_data,1,5) %in% c("G11.0","G11.1","G11.2","G11.3","G11.8","G11.9","G25.5","G31.2","G31.8","G31.9","G93.4","R47.0")) | (substr(unlisted_data, 1, 3) %in% c("G10","G12","G13","G20","G21","G23","G35","G36","G37","G90")) # regular expression to find out cases with "Neurological disorders excluding epilepsy" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Neurological <- rowSums(bidata_com)
  data_file$Neurological[data_file$Neurological >= 1] <- 1
  
  #33 Osteoporosis: Osteoporosis and bone disorder -- Weight: 0.49
  bidata_com <- (substr(unlisted_data,1,5) %in% c("M81.0","M81.1","M81.5","M81.8","M81.9","M83.1","M83.2","M83.3","M83.4","M83.5","M83.8","M83.9","M86.3","M86.4","M86.5","M86.6")) | (substr(unlisted_data, 1, 3) %in% c("M80","M85","M88")) # regular expression to find out cases with "Osteoporosis and bone disorder" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Osteoporosis <- rowSums(bidata_com)
  data_file$Osteoporosis[data_file$Osteoporosis >= 1] <- 1
  
  #34 quadriplegia: Hemi/para/quadriplegia -- Weight: 1.03
  bidata_com <- (substr(unlisted_data,1,5) %in% c("G04.1","G11.4","G80.0","G80.1","G80.2","G83.0","G83.1","G83.2","G83.3","G83.4","G83.9")) | (substr(unlisted_data, 1, 3) %in% c("G81","G82")) # regular expression to find out cases with "Hemi/para/quadriplegia" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$quadriplegia <- rowSums(bidata_com)
  data_file$quadriplegia[data_file$quadriplegia >= 1] <- 1
  
  #35 PeripheralNerve: Peripheral nerve/muscular disorder -- Weight: 1.20
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I73.1","I73.8","I73.9","I77.1","K55.1","K55.2","K55.8","K55.9")) | (substr(unlisted_data, 1, 3) %in% c("I70","I71","I72")) # regular expression to find out cases with "Peripheral nerve/muscular disorder" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$PeripheralNerve <- rowSums(bidata_com)
  data_file$PeripheralNerve[data_file$PeripheralNerve >= 1] <- 1
  
  #36 Pulmonarycir: Pulmonary circulation disorder -- Weight: 0.95
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I28.0","I28.1","I28.8","I28.9")) | (substr(unlisted_data, 1, 3) %in% c("I26","I27")) # regular expression to find out cases with "Pulmonary circulation disorder" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Pulmonarycir <- rowSums(bidata_com)
  data_file$Pulmonarycir[data_file$Pulmonarycir >= 1] <- 1
  
  #37 Sleep: Sleep disorder -- Weight: 1.41
  bidata_com <- (substr(unlisted_data,1,5) %in% c("G47.0","G47.1","G47.2","G47.3")) | (substr(unlisted_data, 1, 3) %in% c("F51")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Sleep <- rowSums(bidata_com)
  data_file$Sleep[data_file$Sleep >= 1] <- 1
  
  #38 Urinary: (Chronic)Urinary tract disorders -- Weight: 0.12
  bidata_com <- (substr(unlisted_data,1,5) %in% c("N30.1","N30.2","N31")) | (substr(unlisted_data, 1, 3) %in% c("N32","N35","N36")) # regular expression to find out cases with "(Chronic)Urinary tract disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Urinary <- rowSums(bidata_com)
  data_file$Urinary[data_file$Urinary >= 1] <- 1
  
  #39 Venous: Venous insufficiency -- Weight: 0.70
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I83.0","I83.2","I87.2")) # regular expression to find out cases with "Venous insufficiency" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Venous <- rowSums(bidata_com)
  data_file$Venous[data_file$Venous >= 1] <- 1
  
  #40 OtherMAL: Other malignacy -- Weight: 0.17
  bidata_com <- (substr(unlisted_data,1,5) %in% c("M81.0","M81.1","M81.5","M81.8","M81.9","M83.1","M83.2","M83.3","M83.4","M83.5","M83.8","M83.9","M86.3","M86.4","M86.5","M86.6")) | (substr(unlisted_data, 1, 3) %in% c(paste("C0",0:9,sep=""),paste("C1",0:9,sep=""),"C20","C21","C23","C24","C25","C26","C30","C31","C32","C33","C37","C38","C39","C43",paste("C4",5:9,sep=""),paste("C5",0:8,sep=""),paste("C6",0:9,sep=""),"C70",paste("C7",2:5,sep=""),paste("C8",1:5,sep=""),"C88",paste("C9",0:5,sep=""))) # regular expression to find out cases with "Other malignacy" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$OtherMAL <- rowSums(bidata_com)
  data_file$OtherMAL[data_file$OtherMAL >= 1] <- 1
  
  #41 Chronichepatitis: Chronic hepatitis -- Weight: 0.39
  bidata_com <- (substr(unlisted_data,1,5) %in% c("B94.2","Z22.5")) | (substr(unlisted_data, 1, 3) %in% c("B18")) # regular expression to find out cases with "Chronic hepatitis" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Chronichepatitis <- rowSums(bidata_com)
  data_file$Chronichepatitis[data_file$Chronichepatitis >= 1] <- 1
  
  #42 Epilepsy: Epilepsy -- Weight: 1.04
  bidata_com <- (substr(unlisted_data,1,5) %in% c("G40.0","G40.1","G40.2","G40.3","G40.4","G40.6","G40.7","G40.8","G40.9")) | (substr(unlisted_data, 1, 3) %in% c("G41")) # regular expression to find out cases with "Epilepsy" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Epilepsy <- rowSums(bidata_com)
  data_file$Epilepsy[data_file$Epilepsy >= 1] <- 1
  
  ### PART B -- Weighted sum for C3 index---------------------------
  data_file$c3 <- 0.51*data_file$Angina + 0.52*data_file$inflam_bd + 0.59*data_file$Anemia + 0.61*data_file$Metabolic + 0.62*data_file$Othercard + 0.63*data_file$eye + 0.72*data_file$Hypertension + 0.75*data_file$Coagulopathies + 0.77*data_file$arrhythmia + 0.83*data_file$Obesity + 0.88*data_file$Diabetes + 0.92*data_file$Liver + 0.93*data_file$MI + 0.98*data_file$PVD + 1.09*data_file$Cerebrovascular + 1.09*data_file$COPD + 1.10*data_file$valve + 1.26*data_file$CHF + 1.38*data_file$Renal + 1.08*data_file$Alcohol + 0.57*data_file$Anxiety + 0.51*data_file$Connective + 1.35*data_file$Dementia - 0.03*data_file$DiabetesOUT + 0.77*data_file$Endocrine + 0.11*data_file$uppergi + 0.54*data_file$InnerEar + 0.11*data_file$Intestinal + 0.69*data_file$spinal + 0.79*data_file$psychiatric + 1.16*data_file$Nutritional + 1.06*data_file$Neurological + 0.49*data_file$Osteoporosis + 1.03*data_file$quadriplegia + 1.20*data_file$PeripheralNerve + 0.95*data_file$Pulmonarycir + 1.41*data_file$Sleep + 0.12*data_file$Urinary + 0.70*data_file$Venous + 0.17*data_file$OtherMAL + 0.39*data_file$Chronichepatitis + 1.04*data_file$Epilepsy
  
  # All previous ICD-10 coding algorithms were refered to:
  # Sarfati, D., Developing new comorbidity indices for cancer populations using administrative data. 2013, University of Otago: Dunedin. <Dr. Sarfati's doctoral dissertation>
  
  return(data_file)
  
}