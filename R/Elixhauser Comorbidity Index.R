#' To calculate the Elixhauser Comorbidity Index
#'
#' @author Miao Cai <email: miao.cai@slu.edu>
#' @description This file aims to calculate the Elixhauser Comorbidity Index
#' @param data Your data file in which Elixhauser Comorbidity Index is to be calculated
#' @param comorbidity A vector of all comorbidity variables
#' @return data: a new data.frame named "data". This data frame contains a new variable "Elix_Index": The Elixhauser Comorbidity Index, developed by Anne Elixhauser in 1998
#' @references Elixhauser, A., Steiner, C., Harris, D. R., & Coffey, R. M. (1998). Comorbidity measures for use with administrative data. Medical care, 36(1), 8-27.
#' @import dplyr
#' @export
eci <- function(data, comorbidity) {

### PART A: SUBSETTING DATA---------------------------

  data_comorbidity <- subset(data, select = comorbidity)# subset the whole dataset into "data_comorbidity": subset data that only contains comorbidities
  data_comorbidity %>% mutate_if(is.factor, as.character) -> data_comorbidity # this converts factor values into characters
  data_comorbidity[is.na(data_comorbidity)] <- 0 #imputate missing values with zeros
  dim_comorbidity <- dim(data_comorbidity)# save the dimensionality of comorbidity subset data
  unlisted_data <- unlist(data_comorbidity)# unlist the comorbidity subset data(converts the data into a vector. Vectorization speeds up the functions)

  ### PART B: FILTERING COMORBIDITIES---------------------------
  ##1:CHF_Elix, Congestive heart failure
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I09.9", "I11.0", "I13.0", "I13.2", "I25.5", "I42.0", "I42.5", "I42.6", "I42.7", "I42.8", "I42.9", "P29.0")) | (substr(unlisted_data, 1, 3) %in% c("I43", "I50")) # regular expression to find out cases with "Congestive heart failure" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$CHF_Elix <- rowSums(bidata_com)
  data$CHF_Elix[data$CHF_Elix >= 1] <- 1

  #2: CA_Elix, Cardiac arrhythmias
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I44.1", "I44.2", "I44.3", "I45.6", "I45.9", "R00.0", "R00.1", "R00.8", "T82.1", "Z45.0", "Z95.0")) | (substr(unlisted_data, 1, 3) %in% c("I47", "I48", "I49")) # regular expression to find out cases with "Cardiac arrhythmias" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$CA_Elix <- rowSums(bidata_com)
  data$CA_Elix[data$CA_Elix >= 1] <- 1

  #3:VD_Elix, Valvular disease
  bidata_com <- (substr(unlisted_data,1,5) %in% c("A52.0", "I09.1", "I09.8", "Q23.0", "Q23.1", "Q23.2", "Q23.3", "Z95.2", "Z95.3", "Z95.4")) | (substr(unlisted_data, 1, 3) %in% c("I05", "I06", "I07", "I08", "I34", "I35", "I36", "I37", "I38", "I39")) # regular expression to find out cases with "Valvular disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$VD_Elix <- rowSums(bidata_com)
  data$VD_Elix[data$VD_Elix >= 1] <- 1

  #4.PCD_Elix, Pulmonary circulation disorders
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I28.0", "I28.8", "I28.9")) | (substr(unlisted_data, 1, 3) %in% c("I26", "I27")) # regular expression to find out cases with "Pulmonary circulation disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$PCD_Elix <- rowSums(bidata_com)
  data$PCD_Elix[data$PCD_Elix >= 1] <- 1

  #5:PVD_Elix, Peripheral vascular disorders
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I73.1", "I73.8", "I73.9", "I77.1", "I79.0", "I79.2", "K55.1", "K55.8", "K55.9", "Z95.8", "Z95.9")) | (substr(unlisted_data, 1, 3) %in% c("I70", "I71")) # regular expression to find out cases with "Valvular disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$PVD_Elix <- rowSums(bidata_com)
  data$PVD_Elix[data$PVD_Elix >= 1] <- 1

  #6: Hypertensionun_Elix, Hypertension, uncomplicated
  bidata_com <- (substr(unlisted_data, 1, 3) == "I10") # regular expression to find out cases with "Hypertension, uncomplicated" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Hypertensionun_Elix <- rowSums(bidata_com)
  data$Hypertensionun_Elix[data$Hypertensionun_Elix >= 1] <- 1

  #7: Hypertension_Elix, Hypertension, complicated
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("I11", "I12", "I13", "I15")) # regular expression to find out cases with "Hypertension, complicated" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Hypertension_Elix <- rowSums(bidata_com)
  data$Hypertension_Elix[data$Hypertension_Elix >= 1] <- 1

  #8: Paralysis_Elix
  bidata_com <- (substr(unlisted_data,1,5) %in% c("G04.1", "G11.4", "G80.1", "G80.2", "G83.9", "G83.0", "G83.1", "G83.2", "G83.3", "G83.4")) | (substr(unlisted_data, 1, 3) %in% c("G81", "G82")) # regular expression to find out cases with "Paralysis" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Paralysis_Elix <- rowSums(bidata_com)
  data$Paralysis_Elix[data$Paralysis_Elix >= 1] <- 1

  #9: OND_Elix, Other neurological disorders
  bidata_com <- (substr(unlisted_data,1,5) %in% c("G25.4", "G25.5", "G31.2", "G31.8", "G31.9", "G93.1", "G93.4", "R47.0")) | (substr(unlisted_data, 1, 3) %in% c("G10", "G20", "G35", "G41", "G11", "G21", "G36", "R56", "G12", "G22", "G37", "G13", "G32", "G40")) # regular expression to find out cases with "Other neurological disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$OND_Elix <- rowSums(bidata_com)
  data$OND_Elix[data$OND_Elix >= 1] <- 1

  #10: CPD_Elix, Chronic pulmonary disease
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I27.8", "I27.9", "J68.4", "J70.1", "J70.3")) | (substr(unlisted_data, 1, 3) %in% c("J40", "J60", "J41", "J61", "J42", "J62", "J43", "J63", "J44", "J64", "J45", "J65", "J46", "J66", "J47", "J67")) # regular expression to find out cases with "Chronic pulmonary disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$CPD_Elix <- rowSums(bidata_com)
  data$CPD_Elix[data$CPD_Elix >= 1] <- 1

  #11:DMun_Elix, Diabetes, uncomplicated
  bidata_com <- (substr(unlisted_data,1,5) %in% c("E10.0", "E10.1", "E10.9", "E11.0", "E11.1", "E11.9", "E12.0", "E12.1", "E12.9", "E13.0", "E13.1", "E13.9", "E14.0", "E14.1", "E14.9")) # regular expression to find out cases with "Diabetes, uncomplicated" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$DMun_Elix <- rowSums(bidata_com)
  data$DMun_Elix[data$DMun_Elix >= 1] <- 1

  #12: DM_Elix, Diabetes, complicated
  bidata_com <- (substr(unlisted_data,1,5) %in% c("E10.2", "E11.2", "E12.2", "E13.2", "E14.2", "E10.3", "E11.3", "E12.3", "E13.3", "E14.3", "E10.4", "E11.4", "E12.4", "E13.4", "E14.4", "E10.5", "E11.5", "E12.5", "E13.5", "E14.5", "E10.6", "E11.6", "E12.6", "E13.6", "E14.6", "E10.7", "E11.7", "E12.7", "E13.7", "E14.7", "E10.8", "E11.8", "E12.8", "E13.8", "E14.8")) # regular expression to find out cases with "Diabetes, complicated" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$DM_Elix <- rowSums(bidata_com)
  data$DM_Elix[data$DM_Elix >= 1] <- 1

  #13: Hypothyroidism_Elix
  bidata_com <- (substr(unlisted_data,1,5) == "E89.0") | (substr(unlisted_data, 1, 3) %in% c("E00", "E01", "E02", "E03")) # regular expression to find out cases with "Hypothyroidism_Elix" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Hypothyroidism_Elix <- rowSums(bidata_com)
  data$Hypothyroidism_Elix[data$Hypothyroidism_Elix >= 1] <- 1

  #14: Renal_Elix, Renal failure
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I12.0", "I13.1", "N25.0", "Z49.0", "Z49.1", "Z49.2", "Z94.0", "Z99.2")) | (substr(unlisted_data, 1, 3) %in% c("N18", "N19")) # regular expression to find out cases with "Renal failure" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Renal_Elix <- rowSums(bidata_com)
  data$Renal_Elix[data$Renal_Elix >= 1] <- 1

  #15: LD_Elix, Liver disease
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I86.4", "I98.2", "K71.1", "K71.3", "K71.4", "K71.5", "K71.7", "K76.0", "Z94.4", "k76.2", "k76.3", "k76.4", "k76.5", "k76.6", "k76.7", "k76.8", "k76.9")) | (substr(unlisted_data, 1, 3) %in% c("B18", "I85", "K70", "K72", "K73", "K74")) # regular expression to find out cases with "Liver disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$LD_Elix <- rowSums(bidata_com)
  data$LD_Elix[data$LD_Elix >= 1] <- 1

  #16: PUD_Elix, Peptic ulcer disease excluding bleeding
  bidata_com <- (substr(unlisted_data,1,5) %in% c("K25.7", "K25.9", "K26.7", "K26.9", "K27.7", "K27.9", "K28.7", "K28.9")) # regular expression to find out cases with "Peptic ulcer disease excluding bleeding" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$PUD_Elix <- rowSums(bidata_com)
  data$PUD_Elix[data$PUD_Elix >= 1] <- 1

  #17: AIDS_Elix,  AIDS/HIV
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("B20", "B21", "B22", "B24")) # regular expression to find out cases with "AIDS/HIV" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$AIDS_Elix <- rowSums(bidata_com)
  data$AIDS_Elix[data$AIDS_Elix >= 1] <- 1

  #18: Lymphoma_Elix
  bidata_com <- (substr(unlisted_data,1,5) %in% c("C90.0", "C90.2")) | (substr(unlisted_data, 1, 3) %in% c("C81", "C82", "C83", "C84", "C85", "C88", "C96")) # regular expression to find out cases with "Lymphoma" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Lymphoma_Elix <- rowSums(bidata_com)
  data$Lymphoma_Elix[data$Lymphoma_Elix >= 1] <- 1

  #19: Metastatic_Elix, Metastatic cancer
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("C77", "C78", "C79", "C80")) # regular expression to find out cases with "Metastatic cancer" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Metastatic_Elix <- rowSums(bidata_com)
  data$Metastatic_Elix[data$Metastatic_Elix >= 1] <- 1

  #20: Solid_Elix, Solid tumor without metastasis，
  bidata_com <- (substr(unlisted_data, 1, 3) %in% c("C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C30", "C31", "C32", "C33", "C34", "C37", "C38", "C39", "C40", "C41", "C43", "C45", "C46", "C47", "C48", "C49", "C50", "C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58", "C60", "C61", "C62", "C63", "C64", "C65", "C66", "C67", "C68", "C69", "C70", "C71", "C72", "C73", "C74", "C75", "C76", "C97")) # regular expression to find out cases with "Solid_Elix" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Solid_Elix <- rowSums(bidata_com)
  data$Solid_Elix[data$Solid_Elix >= 1] <- 1

  #21: Rheumatoid_Elix, Rheumatoid arthritis/collagen vascular diseases
  bidata_com <- (substr(unlisted_data,1,5) %in% c("L94.0", "L94.1", "L94.3", "M46.1", "M46.8", "M46.9", "M12.0", "M12.3", "M31.0", "M31.1", "M31.2", "M31.3")) | (substr(unlisted_data, 1, 3) %in% c("M32", "M33", "M34", "M35", "M45", "M05", "M06", "M08", "M30")) # regular expression to find out cases with "Rheumatoid arthritis/collagen vascular diseases" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Rheumatoid_Elix <- rowSums(bidata_com)
  data$Rheumatoid_Elix[data$Rheumatoid_Elix >= 1] <- 1

  #22: Coagulopathy_Elix
  bidata_com <- (substr(unlisted_data,1,5) %in% c("D69.1", "D69.3", "D69.4", "D69.5", "D69.6")) | (substr(unlisted_data, 1, 3) %in% c("D65", "D66", "D67", "D68")) # regular expression to find out cases with "Coagulopathy" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Coagulopathy_Elix <- rowSums(bidata_com)
  data$Coagulopathy_Elix[data$Coagulopathy_Elix >= 1] <- 1

  #23: Obesity_Elix
  bidata_com <- (substr(unlisted_data, 1, 3) == "E66") # regular expression to find out cases with "Obesity" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Obesity_Elix <- rowSums(bidata_com)
  data$Obesity_Elix[data$Obesity_Elix >= 1] <- 1

  #24: Weight_Elix, Weight loss
  bidata_com <- (substr(unlisted_data,1,5) == "R63.4") | (substr(unlisted_data, 1, 3) %in% c("E40", "E41", "E42", "E43", "E44", "E45", "E46", "R64")) # regular expression to find out cases with "Valvular disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Weight_Elix <- rowSums(bidata_com)
  data$Weight_Elix[data$Weight_Elix >= 1] <- 1

  #25: Fluid_Elix, Fluid and electrolyte disorders
  bidata_com <- (substr(unlisted_data,1,5) == "E22.2") | (substr(unlisted_data, 1, 3) %in% c("E86", "E87")) # regular expression to find out cases with "Valvular disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Fluid_Elix <- rowSums(bidata_com)
  data$Fluid_Elix[data$Fluid_Elix >= 1] <- 1

  #26: Blood_Elix, Blood loss anemia
  bidata_com <- (substr(unlisted_data,1,5) == "D50.0") # regular expression to find out cases with "Valvular disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Blood_Elix <- rowSums(bidata_com)
  data$Blood_Elix[data$Blood_Elix >= 1] <- 1

  #27: Deficiency_Elix, Deficiency anemia
  bidata_com <- (substr(unlisted_data,1,5) %in% c("D50.8", "D50.9")) | (substr(unlisted_data, 1, 3) %in% c("D51", "D52", "D53")) # regular expression to find out cases with "Valvular disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Deficiency_Elix <- rowSums(bidata_com)
  data$Deficiency_Elix[data$Deficiency_Elix >= 1] <- 1

  #28: Alcohol_Elix, Alcohol abuse
  bidata_com <- (substr(unlisted_data,1,5) %in% c("G62.1", "I42.6", "K29.2", "K70.0", "K70.3", "K70.9", "Z50.2", "Z71.4", "Z72.1")) | (substr(unlisted_data, 1, 3) %in% c("F10", "E52", "T51")) # regular expression to find out cases with "Alcohol abuse" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Alcohol_Elix <- rowSums(bidata_com)
  data$Alcohol_Elix[data$Alcohol_Elix >= 1] <- 1

  #29: Drug_Elix, Drug abuse
  bidata_com <- (substr(unlisted_data,1,5) %in% c("Z71.5", "Z72.2")) | (substr(unlisted_data, 1, 3) %in% c("F11", "F12", "F13", "F14", "F15", "F16", "F18", "F19")) # regular expression to find out cases with "Drug abuse" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Drug_Elix <- rowSums(bidata_com)
  data$Drug_Elix[data$Drug_Elix >= 1] <- 1

  #30: Psychoses_Elix
  bidata_com <- (substr(unlisted_data,1,5) %in% c("F30.2", "F31.2", "F31.5")) | (substr(unlisted_data, 1, 3) %in% c("F20", "F22", "F23", "F24", "F25", "F28", "F29")) # regular expression to find out cases with "Psychoses" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Psychoses_Elix <- rowSums(bidata_com)
  data$Psychoses_Elix[data$Psychoses_Elix >= 1] <- 1

  #31: Depression_Elix
  bidata_com <- (substr(unlisted_data,1,5) %in% c("F20.4", "F31.3", "F31.4", "F31.5", "F34.1", "F41.2", "F43.2")) | (substr(unlisted_data, 1, 3) %in% c("F32", "F33")) # regular expression to find out cases with "Depression" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Depression_Elix <- rowSums(bidata_com)
  data$Depression_Elix[data$Depression_Elix >= 1] <- 1

  # About age: age is not included in calculating the Elixhauser Comorbidity Index.

  ###PART D---WEIGHTED SUM---------------------------
  # Refer to:
  # Elixhauser, A., Steiner, C., Harris, D. R., & Coffey, R. M. (1998). Comorbidity measures for use with administrative data. Medical care, 36(1), 8-27.

  data$Elix_Index <- 7 * data$CHF_Elix + 5 * data$CA_Elix - 1 * data$VD_Elix + 4 * data$PCD_Elix + 2 * data$PVD_Elix + 0 * data$Hypertensionun_Elix + 0 * data$Hypertension_Elix + 7 * data$Paralysis_Elix + 6 * data$OND_Elix + 3 * data$CPD_Elix + 0 * data$DMun_Elix + 0 * data$DM_Elix + 0 * data$Hypothyroidism_Elix + 5 * data$Renal_Elix + 11 * data$LD_Elix + 0 * data$PUD_Elix + 0 * data$AIDS_Elix + 9 * data$Lymphoma_Elix + 12 * data$Metastatic_Elix + 4 * data$Solid_Elix + 0 * data$Rheumatoid_Elix + 3 * data$Coagulopathy_Elix - 4 * data$Obesity_Elix + 6 * data$Weight_Elix + 5 * data$Fluid_Elix - 2 * data$Blood_Elix - 2 * data$Deficiency_Elix + 0 * data$Alcohol_Elix - 7 * data$Drug_Elix + 0 * data$Psychoses_Elix - 3 * data$Depression_Elix

  # All previous ICD-10 coding algorithms were refered to:
  # Quan, H., Sundararajan, V., Halfon, P., Fong, A., Burnand, B., Luthi, J. C., ... & Ghali, W. A. (2005). Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data. Medical care, 1130-1139.

  return(data)
  }
