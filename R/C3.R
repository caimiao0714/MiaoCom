#' To calculate the C3 Comobidity Index (targeted at cancer patients)
#' @author Miao Cai <email: miao.cai@slu.edu>
#' @description This file aims to calculate the C3 Comobidity Index
#' @param data Your data file in which Elixhauser Comorbidity Index is to be calculated
#' @param comorbidity A vector of all comorbidity variables
#' @return data_file: a new data.frame named "data_file". This data frame contains a new variable "c3": c3 Index developed by Sarfati et al. in 2014. This is a comorbidity index targeted at cancer patients.
#' @export

# The C3 index refers to:
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
  ##1 Angina: Angina -- 0.51
  bidata_com <- (substr(unlisted_data, 1, 3) == "I20")  # regular expression to find out cases with "Angina" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Angina <- rowSums(bidata_com)
  data_file$Angina[data_file$Angina >= 1] <- 1
  
  #2 inflam_bd: Inflammatory bowel disease -- 0.52 
  bidata_com <- (substr(unlisted_data,1,5) %in% c("K52.2","K52.8","K52.9")) | (substr(unlisted_data, 1, 3) %in% c("K50","K51")) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$inflam_bd <- rowSums(bidata_com)
  data_file$inflam_bd[data_file$inflam_bd >= 1] <- 1
  
  #3 Anemia: Anemia -- 0.59
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$inflam_bd <- rowSums(bidata_com)
  data_file$inflam_bd[data_file$inflam_bd >= 1] <- 1
  
  #4 Metabolic: Metabolic conditions -- 0.61
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Metabolic <- rowSums(bidata_com)
  data_file$Metabolic[data_file$Metabolic >= 1] <- 1
  
  #5 Othercard: Other cardiac conditions -- 0.62
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Othercard <- rowSums(bidata_com)
  data_file$Othercard[data_file$Othercard >= 1] <- 1
  
  #6 eye: (Major) Eye conditions affecting vision -- 0.63
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$eye <- rowSums(bidata_com)
  data_file$eye[data_file$eye >= 1] <- 1
  
  #7 Hypertension: Hypertension -- 0.72
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Hypertension <- rowSums(bidata_com)
  data_file$Hypertension[data_file$Hypertension >= 1] <- 1
  
  #8 Coagulopathies: Coagulopathies and other blood disorders -- 0.75
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Coagulopathies <- rowSums(bidata_com)
  data_file$Coagulopathies[data_file$Coagulopathies >= 1] <- 1
  
  #9 arrhythmia: Cardiac arrhythmia -- 0.77
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$arrhythmia <- rowSums(bidata_com)
  data_file$arrhythmia[data_file$arrhythmia >= 1] <- 1
  
  #10 Obesity: Obesity -- 0.83
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Obesity <- rowSums(bidata_com)
  data_file$Obesity[data_file$Obesity >= 1] <- 1
  
  #11 Diabetes: Diabetes with complications -- 0.88
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Diabetes <- rowSums(bidata_com)
  data_file$Diabetes[data_file$Diabetes >= 1] <- 1
  
  #12 Liver: Liver disease -- 0.92
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Liver <- rowSums(bidata_com)
  data_file$Liver[data_file$Liver >= 1] <- 1
  
  #13 MI: Previous MI -- 0.93
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$MI <- rowSums(bidata_com)
  data_file$MI[data_file$MI >= 1] <- 1
  
  #14 PVD: Peripheral vascular disease -- 0.98
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$PVD <- rowSums(bidata_com)
  data_file$PVD[data_file$PVD >= 1] <- 1
  
  #15 Cerebrovascular: Cerebrovascular disease -- 1.09
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Cerebrovascular <- rowSums(bidata_com)
  data_file$Cerebrovascular[data_file$Cerebrovascular >= 1] <- 1
  
  #16 COPD: COPD and asthma -- 1.09
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$COPD <- rowSums(bidata_com)
  data_file$COPD[data_file$COPD >= 1] <- 1
  
  #17 valve: Cardiac valve disorders -- 1.10
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$valve <- rowSums(bidata_com)
  data_file$valve[data_file$valve >= 1] <- 1
  
  #18 CHF: Congestive heart failure -- 1.26
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$CHF <- rowSums(bidata_com)
  data_file$CHF[data_file$CHF >= 1] <- 1
  
  #19 Renal: Renal disease -- 1.38
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Renal <- rowSums(bidata_com)
  data_file$Renal[data_file$Renal >= 1] <- 1
  
  #20 Alcohol: Alcohol abuse -- 1.08
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Alcohol <- rowSums(bidata_com)
  data_file$Alcohol[data_file$Alcohol >= 1] <- 1
  
  #21 Anxiety: Anxiety/behavioral disorders -- 0.57
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Anxiety <- rowSums(bidata_com)
  data_file$Anxiety[data_file$Anxiety >= 1] <- 1
  
  #22 Connective: Connective tissue disorders -- 0.51
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Connective <- rowSums(bidata_com)
  data_file$Connective[data_file$Connective >= 1] <- 1
  
  #23 Dementia: Dementia -- 1.35
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Dementia <- rowSums(bidata_com)
  data_file$Dementia[data_file$Dementia >= 1] <- 1
  
  #24 DiabetesOUT: Diabetes without complications -- -0.03
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$DiabetesOUT <- rowSums(bidata_com)
  data_file$DiabetesOUT[data_file$DiabetesOUT >= 1] <- 1
  
  #25 Endocrine: Endocrine disorders -- 0.77
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Endocrine <- rowSums(bidata_com)
  data_file$Endocrine[data_file$Endocrine >= 1] <- 1
  
  #26 uppergi: Upper GI disorders -- 0.11
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$uppergi <- rowSums(bidata_com)
  data_file$uppergi[data_file$uppergi >= 1] <- 1
  
  #27 InnerEar: Inner ear disorders -- 0.54
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$InnerEar <- rowSums(bidata_com)
  data_file$InnerEar[data_file$InnerEar >= 1] <- 1
  
  #28 Intestinal: Intestinal disorders -- 0.11
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Intestinal <- rowSums(bidata_com)
  data_file$Intestinal[data_file$Intestinal >= 1] <- 1
  
  #29 spinal: Joint and spinal disorders -- 0.69
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$spinal <- rowSums(bidata_com)
  data_file$spinal[data_file$spinal >= 1] <- 1
  
  #30 psychiatric: Major psychiatric disorders -- 0.79
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$psychiatric <- rowSums(bidata_com)
  data_file$psychiatric[data_file$psychiatric >= 1] <- 1
  
  #31 Nutritional: Nutritional disorders -- 1.16
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Nutritional <- rowSums(bidata_com)
  data_file$Nutritional[data_file$Nutritional >= 1] <- 1
  
  #32 Neurological: Neurological disorders excluding epilepsy -- 1.06
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Neurological <- rowSums(bidata_com)
  data_file$Neurological[data_file$Neurological >= 1] <- 1
  
  #33 Osteoporosis: Osteoporosis and bone disorder -- 0.49
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Osteoporosis <- rowSums(bidata_com)
  data_file$Osteoporosis[data_file$Osteoporosis >= 1] <- 1
  
  #34 quadriplegia: Hemi/para/quadriplegia -- 1.03
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$quadriplegia <- rowSums(bidata_com)
  data_file$quadriplegia[data_file$quadriplegia >= 1] <- 1
  
  #35 PeripheralNerve: Peripheral nerve/muscular disorder -- 1.20
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$PeripheralNerve <- rowSums(bidata_com)
  data_file$PeripheralNerve[data_file$PeripheralNerve >= 1] <- 1
  
  #36 Pulmonarycir: Pulmonary circulation disorder -- 0.95
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Pulmonarycir <- rowSums(bidata_com)
  data_file$Pulmonarycir[data_file$Pulmonarycir >= 1] <- 1
  
  #37 Sleep: Sleep disorder -- 1.41
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Sleep <- rowSums(bidata_com)
  data_file$Sleep[data_file$Sleep >= 1] <- 1
  
  #38 Urinary: (Chronic)Urinary tract disorders -- 0.12
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Urinary <- rowSums(bidata_com)
  data_file$Urinary[data_file$Urinary >= 1] <- 1
  
  #39 Venous: Venous insufficiency -- 0.70
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Venous <- rowSums(bidata_com)
  data_file$Venous[data_file$Venous >= 1] <- 1
  
  #40 OtherMAL: Other malignacy -- 0.17
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$OtherMAL <- rowSums(bidata_com)
  data_file$OtherMAL[data_file$OtherMAL >= 1] <- 1
  
  #41 Chronichepatitis: Chronic hepatitis -- 0.39
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Chronichepatitis <- rowSums(bidata_com)
  data_file$Chronichepatitis[data_file$Chronichepatitis >= 1] <- 1
  
  #42 Epilepsy: Epilepsy -- 1.04
  bidata_com <- (substr(unlisted_data,1,5) %in% c()) | (substr(unlisted_data, 1, 3) %in% c()) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Epilepsy <- rowSums(bidata_com)
  data_file$Epilepsy[data_file$Epilepsy >= 1] <- 1
  
  ### PART B -- Weighted sum for C3 index---------------------------
  data_file$C3 <- 0.51*data_file$Angina + 0.52*data_file$inflam_bd + 0.59*data_file$Anemia + 0.61*data_file$Metabolic + 0.62*data_file$Othercard + 0.63*data_file$eye + 0.72*data_file$Hypertension + 0.75*data_file$Coagulopathies + 0.77*data_file$arrhythmia + 0.83*data_file$Obesity + 0.88*data_file$Diabetes + 0.92*data_file$Liver + 0.93*data_file$MI + 0.98*data_file$PVD + 1.09*data_file$Cerebrovascular + 1.09*data_file$COPD + 1.10*data_file$valve + 1.26*data_file$CHF + 1.38*data_file$Renal + 1.08*data_file$Alcohol + 0.57*data_file$Anxiety + 0.51*data_file$Connective + 1.35*data_file$Dementia - 0.03*data_file$DiabetesOUT + 0.77*data_file$Endocrine + 0.11*data_file$uppergi + 0.54*data_file$InnerEar + 0.11*data_file$Intestinal + 0.69*data_file$spinal + 0.79*data_file$psychiatric + 1.16*data_file$Nutritional + 1.06*data_file$Neurological + 0.49*data_file$Osteoporosis + 1.03*data_file$quadriplegia + 1.20*data_file$PeripheralNerve + 0.95*data_file$Pulmonarycir + 1.41*data_file$Sleep + 0.12*data_file$Urinary + 0.70*data_file$Venous + 0.17*data_file$OtherMAL + 0.39*data_file$Chronichepatitis + 1.04*data_file$Epilepsy
  
  
  
  
  
}