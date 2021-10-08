#' To calculate the c3 Comobidity Index (targeted at cancer patients)
#'
#' @author Miao Cai <email: miao.cai@outlook.com>
#' @description This file aims to calculate the c3 Comobidity Index.
#' @note This C3 index is used specifically for cancer patients.
#' @param data Your data file in which Elixhauser Comorbidity Index is to be calculated
#' @param comorbidity A vector of all comorbidity variables
#' @return data: a new data.frame named "data". This data frame contains a new variable "c3": the c3 Index was developed by Sarfati et al. in 2014. This is a comorbidity index targeted at cancer patients.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate_if
#' @references Sarfati, D., et al., Cancer-specific administrative data-based comorbidity indices provided valid alternative to Charlson and National Cancer Institute Indices. J Clin Epidemiol, 2014. 67(5): p. 586-95
#' @references Sarfati, D., Developing new comorbidity indices for cancer populations using administrative data. 2013, University of Otago: Dunedin. <Dr. Sarfati's doctoral dissertation>
#' @export c3
c3 <- function(data, comorbidity) {
  ### PART A: SUBSETTING DATA---------------------------

  data_comorbidity <- subset(data, select = comorbidity)# subset the whole dataset into "data_comorbidity": subset data that only contains comorbidities
  data_comorbidity %>% dplyr::mutate_if(is.factor, as.character) -> data_comorbidity # this converts factor values into characters
  data_comorbidity[is.na(data_comorbidity)] <- 'NA' #imputate missing values with zeros
  dim_comorbidity <- dim(data_comorbidity)# save the dimensionality of comorbidity subset data
  unlisted_data <- unlist(data_comorbidity)# unlist the comorbidity subset data(converts the data into a vector. Vectorization speeds up the functions)
  unlisted_data <- toupper(unlisted_data)

  ### PART B: FILTERING COMORBIDITIES---------------------------
  ##1 Angina: Angina -- weight: 0.51
  bidata_com <- grepl('I20', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Angina" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Angina <- rowSums(bidata_com)
  data$Angina[data$Angina >= 1] <- 1

  #2 inflam_bd: Inflammatory bowel disease -- Weight: 0.52
  bidata_com <- grepl('K52\\.2|K52\\.8|K52\\.9|K50|K51', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Inflammatory bowel disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$inflam_bd <- rowSums(bidata_com)
  data$inflam_bd[data$inflam_bd >= 1] <- 1

  #3 Anemia: Anemia -- Weight: 0.59
  bidata_com <- grepl('C50|C51|C52|C53', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Anemia" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Anemia <- rowSums(bidata_com)
  data$Anemia[data$Anemia >= 1] <- 1

  #4 Metabolic: Metabolic conditions -- Weight: 0.61
  bidata_com <- grepl('E79\\.1|E79\\.8|E79\\.9|E70|E71|E72|E74|E75|E76|E77|E78|E80|E83|E85|E88', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Metabolic conditions" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Metabolic <- rowSums(bidata_com)
  data$Metabolic[data$Metabolic >= 1] <- 1

  #5 Othercard: Other cardiac conditions -- Weight: 0.62
  bidata_com <- grepl('I24\\.8|I24\\.9|I25\\.0|I25\\.1|I25\\.3|I25\\.4|I25\\.6|I25\\.8|I25\\.9|I31\\.0|I31\\.1|I42\\.1|I42\\.2|I42\\.4', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Other cardiac conditions" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Othercard <- rowSums(bidata_com)
  data$Othercard[data$Othercard >= 1] <- 1

  #6 eye: (Major) Eye conditions affecting vision -- Weight: 0.63
  bidata_com <- grepl('H18\\.1|H18\\.4|H18\\.5|H18\\.6|H20\\.1|H21\\.2|H30\\.1|H31\\.1|H31\\.2|H31\\.3|H31\\.4|H33\\.0|H33\\.2|H33\\.3|H33\\.4|H33\\.5|H53\\.0|H53\\.1|H53\\.2|H53\\.3|H53\\.4|H53\\.6|H53\\.8|H53\\.9|H16|H34|H35|H43|H46|H47|H49|H50|H51|H54|Q12|Q13|Q14|Q15', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$eye <- rowSums(bidata_com)
  data$eye[data$eye >= 1] <- 1

  #7 Hypertension: Hypertension -- Weight: 0.72
  bidata_com <- grepl('I11\\.9|I12\\.9|I13\\.9|I10', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Hypertension" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Hypertension <- rowSums(bidata_com)
  data$Hypertension[data$Hypertension >= 1] <- 1

  #8 Coagulopathies: Coagulopathies and other blood disorders -- Weight: 0.75
  bidata_com <- grepl('D59\\.0|D59\\.1|D59\\.2|D59\\.3|D59\\.4|D59\\.8|D59\\.9|D68\\.0|D68\\.1|D68\\.2|D68\\.8|D68\\.9|D69\\.1|D69\\.2|D69\\.3|D69\\.4|D69\\.6|D69\\.8|D69\\.9|D75\\.0|D75\\.2|D75\\.8|D75\\.9|D55|D56|D57|D58|D60|D61|D64|D64|D66|D67|D70|D71|D72|D74', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Coagulopathies <- rowSums(bidata_com)
  data$Coagulopathies[data$Coagulopathies >= 1] <- 1

  #9 arrhythmia: Cardiac arrhythmia -- Weight: 0.77
  bidata_com <- grepl('I44\\.1|I44\\.2|I44\\.3|I45\\.6|I45\\.9|T82\\.1|Z45\\.0|Z95\\.0|I47|I48|I49', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Cardiac arrhythmia" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$arrhythmia <- rowSums(bidata_com)
  data$arrhythmia[data$arrhythmia >= 1] <- 1

  #10 Obesity: Obesity -- Weight: 0.83
  bidata_com <- grepl('E66', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Obesity <- rowSums(bidata_com)
  data$Obesity[data$Obesity >= 1] <- 1

  #11 Diabetes: Diabetes with complications -- Weight: 0.88
  bidata_com <- grepl('E10\\.2|E10\\.3|E10\\.4|E10\\.5|E10\\.6|E10\\.7|E10\\.8|E11\\.2|E11\\.3|E11\\.4|E11\\.5|E11\\.6|E11\\.7|E11\\.8|E12\\.2|E12\\.3|E12\\.4|E12\\.5|E12\\.6|E12\\.7|E12\\.8|E13\\.2|E13\\.3|E13\\.4|E13\\.5|E13\\.6|E13\\.7|E13\\.8|E14\\.2|E14\\.3|E14\\.4|E14\\.5|E14\\.6|E14\\.7|E14\\.8', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Diabetes <- rowSums(bidata_com)
  data$Diabetes[data$Diabetes >= 1] <- 1

  #12 Liver: Liver disease -- Weight: 0.92
  bidata_com <- grepl('K71\\.1|K71\\.3|K71\\.4|K71\\.5|K71\\.7|K72\\.1|K72\\.9|K76\\.0|K76\\.2|K76\\.3|K76\\.4|K76\\.5|K76\\.6|K76\\.7|K76\\.8|K76\\.9|I86\\.4|I98\\.2|Z94\\.4|K70|K73|K74|I85', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Liver disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Liver <- rowSums(bidata_com)
  data$Liver[data$Liver >= 1] <- 1

  #13 MI: Previous MI -- Weight: 0\\.93
  bidata_com <- grepl('I24\\.1|I25\\.2|I21|I22|I23', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$MI <- rowSums(bidata_com)
  data$MI[data$MI >= 1] <- 1

  #14 PVD: Peripheral vascular disease -- Weight: 0.98
  bidata_com <- grepl('I73\\.1|I73\\.8|I73\\.9|I77\\.1|K55\\.1|K55\\.2|K55\\.8|K55\\.9|I70|I71|I72', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$PVD <- rowSums(bidata_com)
  data$PVD[data$PVD >= 1] <- 1

  #15 Cerebrovascular: Cerebrovascular disease -- Weight: 1.09
  bidata_com <- grepl('I60|I61|I62|I63|I64|I65|I66|I67|I69|G45|G46', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Cerebrovascular <- rowSums(bidata_com)
  data$Cerebrovascular[data$Cerebrovascular >= 1] <- 1

  #16 COPD: COPD and asthma -- Weight: 1.09
  bidata_com <- grepl('J68\\.4|J70\\.1|J70\\.3|J96\\.1|J98\\.0|J98\\.2|J98\\.3|J98\\.4|E84|J40|J41|J42|J43|J44|J45|J46|J47|J60|J61|J62|J63|J64|J65|J66|J67|J84', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$COPD <- rowSums(bidata_com)
  data$COPD[data$COPD >= 1] <- 1

  #17 valve: Cardiac valve disorders -- Weight: 1.10
  bidata_com <- grepl('I09\\.1|I09\\.8|T82\\.0|Q23\\.0|Q23\\.1|Q23\\.2|Q23\\.3|Q23\\.8|Q23\\.9|Z95\\.2|Z95\\.3|Z95\\.4|I05|I06|I07|I08|I34|I35|I36|I37|I38', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$valve <- rowSums(bidata_com)
  data$valve[data$valve >= 1] <- 1

  #18 CHF: Congestive heart failure -- Weight: 1.26
  bidata_com <- grepl('I09\\.9|I11\\.0|I13\\.0|I13\\.2|I25\\.5|I42\\.0|I42\\.5|I42\\.6|I42\\.7|I42\\.8|I42\\.9|I43|I50', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$CHF <- rowSums(bidata_com)
  data$CHF[data$CHF >= 1] <- 1

  #19 Renal: Renal disease -- Weight: 1.38
  bidata_com <- grepl('N03\\.2|N03\\.3|N03\\.4|N03\\.5|N03\\.6|N03\\.7|N03\\.8|N03\\.9|N04\\.2|N04\\.3|N04\\.4|N04\\.5|N04\\.6|N04\\.7|N04\\.8|N04\\.9|N05\\.2|N05\\.3|N05\\.4|N05\\.5|N05\\.6|N05\\.7|N05\\.8|N05\\.9|N25\\.0|N25\\.8|N25\\.9|I12\\.0|I13\\.1|Z94\\.0|Z99\\.2|N11|N18|N19|Z49', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Renal <- rowSums(bidata_com)
  data$Renal[data$Renal >= 1] <- 1

  #20 Alcohol: Alcohol abuse -- Weight: 1.08
  bidata_com <- grepl('F10\\.1|F10\\.2|F10\\.3|F10\\.4|F10\\.5|F10\\.6|F10\\.7|F10\\.8|F10\\.9|Z50\\.2|Z71\\.4', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Alcohol <- rowSums(bidata_com)
  data$Alcohol[data$Alcohol >= 1] <- 1

  #21 Anxiety: Anxiety/behavioral disorders -- Weight: 0.57
  bidata_com <- grepl('F40|F41|F42|F44|F45|F48|F50|F55|F59|F60|F61|F63|F64|F65|F66|F68|F69', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Anxiety <- rowSums(bidata_com)
  data$Anxiety[data$Anxiety >= 1] <- 1

  #22 Connective: Connective tissue disorders -- Weight: 0.51
  bidata_com <- grepl('M12\\.0|M12\\.3|N25\\.9|I12\\.0|M35\\.0|M35\\.1|M35\\.2|M35\\.3|M35\\.4|M35\\.5|M35\\.6|M35\\.8|M35\\.9|L93|M05|M06|M08|M30|M31|M32|M33|M34', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Connective <- rowSums(bidata_com)
  data$Connective[data$Connective >= 1] <- 1

  #23 Dementia: Dementia -- Weight: 1.35
  bidata_com <- grepl('F02\\.0|F02\\.1|F02\\.2|F02\\.3|F05\\.1|G31\\.0|G31\\.1|F00|F01|F03|G30', unlisted_data, perl = TRUE, ignore.case = TRUE)
   # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Dementia <- rowSums(bidata_com)
  data$Dementia[data$Dementia >= 1] <- 1

  #24 DiabetesOUT: Diabetes without complications -- Weight: -0.03
  bidata_com <- grepl('E10\\.0|E10\\.1|E10\\.9|E11\\.0|E11\\.1|E11\\.9|E12\\.0|E12\\.1|E12\\.9|E13\\.0|E13\\.1|E13\\.9|E14\\.0|E14\\.1|E14\\.9', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Diabetes without complications" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$DiabetesOUT <- rowSums(bidata_com)
  data$DiabetesOUT[data$DiabetesOUT >= 1] <- 1

  #25 Endocrine: Endocrine disorders -- Weight: 0.77
  bidata_com <- grepl('E06\\.2|E06\\.3|E06\\.5|E16\\.3|E16\\.4|E16\\.8|E16\\.9|E21\\.0|E21\\.2|E21\\.3|E21\\.4|E21\\.5|E23\\.0|E23\\.2|E23\\.3|E23\\.6|E23\\.7|E24\\.0|E24\\.1|E24\\.3|E24\\.4|E24\\.8|E24\\.9|E34\\.5|E34\\.8|E34\\.9|E01|E02|E03|E05|E07|E20|E22|E25|E26|E27|E31|E32', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Endocrine disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Endocrine <- rowSums(bidata_com)
  data$Endocrine[data$Endocrine >= 1] <- 1

  #26 uppergi: Upper GI disorders -- Weight: 0.11
  bidata_com <- grepl('K22\\.0|K22\\.1|K22\\.4|K22\\.5|K22\\.8|K22\\.9|K31\\.1|K31\\.2|K31\\.4|K31\\.6|K25|K26|K27|K28', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Upper GI disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$uppergi <- rowSums(bidata_com)
  data$uppergi[data$uppergi >= 1] <- 1

  #27 InnerEar: Inner ear disorders -- Weight: 0.54
  bidata_com <- grepl('H91\\.1|H91\\.0|H91\\.3|H91\\.8|H91\\.9|H93\\.0|H93\\.1|H93\\.2|H93\\.3|H80|H81|H83|H90', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Inner ear disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$InnerEar <- rowSums(bidata_com)
  data$InnerEar[data$InnerEar >= 1] <- 1

  #28 Intestinal: Intestinal disorders -- Weight: 0.11
  bidata_com <- grepl('K59\\.2|K59\\.3|K57|K90', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Intestinal disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Intestinal <- rowSums(bidata_com)
  data$Intestinal[data$Intestinal >= 1] <- 1

  #29 spinal: Joint and spinal disorders -- Weight: 0.69
  bidata_com <- grepl('M15\\.0|M15\\.1|M15\\.2|M15\\.4|M15\\.8|M15\\.9|M40\\.0|M40\\.2|M40\\.3|M40\\.4|M40\\.5|M46\\.0|M46\\.1|M46\\.2|M48\\.0|M48\\.1|M48\\.2|M48\\.5|M48\\.8|M48\\.9|G95\\.0|G95\\.1|M41|M42|M43|M45|M47', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Joint and spinal disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$spinal <- rowSums(bidata_com)
  data$spinal[data$spinal >= 1] <- 1

  #30 psychiatric: Major psychiatric disorders -- Weight: 0.79
  bidata_com <- grepl('F30\\.2|F32\\.1|F32\\.2|F32\\.3|F32\\.8|F32\\.9|F20|F22|F25|F28|F29|F31|F33|F39', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Major psychiatric disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$psychiatric <- rowSums(bidata_com)
  data$psychiatric[data$psychiatric >= 1] <- 1

  #31 Nutritional: Nutritional disorders -- Weight: 1.16
  bidata_com <- grepl('E40|E41|E42|E43|E44|E45|E46|E50|E51|E52|E53|E54|E55|E56|E58|E59|E60|E61|E63|E64', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Nutritional disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Nutritional <- rowSums(bidata_com)
  data$Nutritional[data$Nutritional >= 1] <- 1

  #32 Neurological: Neurological disorders excluding epilepsy -- Weight: 1.06
  bidata_com <- grepl('G11\\.0|G11\\.1|G11\\.2|G11\\.3|G11\\.8|G11\\.9|G25\\.5|G31\\.2|G31\\.8|G31\\.9|G93\\.4|R47\\.0|G10|G12|G13|G20|G21|G23|G35|G36|G37|G90', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Neurological disorders excluding epilepsy" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Neurological <- rowSums(bidata_com)
  data$Neurological[data$Neurological >= 1] <- 1

  #33 Osteoporosis: Osteoporosis and bone disorder -- Weight: 0.49
  bidata_com <- grepl('M81\\.0|M81\\.1|M81\\.5|M81\\.8|M81\\.9|M83\\.1|M83\\.2|M83\\.3|M83\\.4|M83\\.5|M83\\.8|M83\\.9|M86\\.3|M86\\.4|M86\\.5|M86\\.6|M80|M85|M88', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Osteoporosis and bone disorder" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Osteoporosis <- rowSums(bidata_com)
  data$Osteoporosis[data$Osteoporosis >= 1] <- 1

  #34 quadriplegia: Hemi/para/quadriplegia -- Weight: 1.03
  bidata_com <- grepl('G04\\.1|G11\\.4|G80\\.0|G80\\.1|G80\\.2|G83\\.0|G83\\.1|G83\\.2|G83\\.3|G83\\.4|G83\\.9|G81|G82', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Hemi/para/quadriplegia" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$quadriplegia <- rowSums(bidata_com)
  data$quadriplegia[data$quadriplegia >= 1] <- 1

  #35 PeripheralNerve: Peripheral nerve/muscular disorder -- Weight: 1.20
  bidata_com <- grepl('I73\\.1|I73\\.8|I73\\.9|I77\\.1|K55\\.1|K55\\.2|K55\\.8|K55\\.9|I70|I71|I72', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Peripheral nerve/muscular disorder" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$PeripheralNerve <- rowSums(bidata_com)
  data$PeripheralNerve[data$PeripheralNerve >= 1] <- 1

  #36 Pulmonarycir: Pulmonary circulation disorder -- Weight: 0.95
  bidata_com <- grepl('I28\\.0|I28\\.1|I28\\.8|I28\\.9|I26|I27', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Pulmonary circulation disorder" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Pulmonarycir <- rowSums(bidata_com)
  data$Pulmonarycir[data$Pulmonarycir >= 1] <- 1

  #37 Sleep: Sleep disorder -- Weight: 1.41
  bidata_com <- grepl('G47\\.0|G47\\.1|G47\\.2|G47\\.3|F51', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Sleep <- rowSums(bidata_com)
  data$Sleep[data$Sleep >= 1] <- 1

  #38 Urinary: (Chronic)Urinary tract disorders -- Weight: 0.12
  bidata_com <- grepl('N30\\.1|N30\\.2|N31|N32|N35|N36', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "(Chronic)Urinary tract disorders" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Urinary <- rowSums(bidata_com)
  data$Urinary[data$Urinary >= 1] <- 1

  #39 Venous: Venous insufficiency -- Weight: 0.70
  bidata_com <- grepl('I83\\.0|I83\\.2|I87\\.2', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Venous insufficiency" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Venous <- rowSums(bidata_com)
  data$Venous[data$Venous >= 1] <- 1

  #40 OtherMAL: Other malignacy -- Weight: 0.17
  bidata_com <- grepl('M81\\.0|M81\\.1|M81\\.5|M81\\.8|M81\\.9|M83\\.1|M83\\.2|M83\\.3|M83\\.4|M83\\.5|M83\\.8|M83\\.9|M86\\.3|M86\\.4|M86\\.5|M86\\.6|C00|C01|C02|C03|C04|C05|C06|C07|C08|C09|C10|C11|C12|C13|C14|C15|C16|C17|C18|C19|C20|C21|C23|C24|C25|C26|C30|C31|C32|C33|C37|C38|C39|C43|C45|C46|C47|C48|C49|C50|C51|C52|C53|C54|C55|C56|C57|C58|C60|C61|C62|C63|C64|C65|C66|C67|C68|C69|C70|C72|C73|C74|C75|C81|C82|C83|C84|C85|C88|C90|C91|C92|C93|C94|C95', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Other malignacy" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$OtherMAL <- rowSums(bidata_com)
  data$OtherMAL[data$OtherMAL >= 1] <- 1

  #41 Chronichepatitis: Chronic hepatitis -- Weight: 0.39
  bidata_com <- grepl('B94\\.2|Z22\\.5|B18', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Chronic hepatitis" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Chronichepatitis <- rowSums(bidata_com)
  data$Chronichepatitis[data$Chronichepatitis >= 1] <- 1

  #42 Epilepsy: Epilepsy -- Weight: 1.04
  bidata_com <- grepl('G40\\.0|G40\\.1|G40\\.2|G40\\.3|G40\\.4|G40\\.6|G40\\.7|G40\\.8|G40\\.9|G41', unlisted_data, perl = TRUE, ignore.case = TRUE) # regular expression to find out cases with "Epilepsy" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data$Epilepsy <- rowSums(bidata_com)
  data$Epilepsy[data$Epilepsy >= 1] <- 1

  ### PART B -- Weighted sum for C3 index---------------------------
  data$c3 <- 0.51*data$Angina + 0.52*data$inflam_bd + 0.59*data$Anemia + 0.61*data$Metabolic + 0.62*data$Othercard + 0.63*data$eye + 0.72*data$Hypertension + 0.75*data$Coagulopathies + 0.77*data$arrhythmia + 0.83*data$Obesity + 0.88*data$Diabetes + 0.92*data$Liver + 0.93*data$MI + 0.98*data$PVD + 1.09*data$Cerebrovascular + 1.09*data$COPD + 1.10*data$valve + 1.26*data$CHF + 1.38*data$Renal + 1.08*data$Alcohol + 0.57*data$Anxiety + 0.51*data$Connective + 1.35*data$Dementia - 0.03*data$DiabetesOUT + 0.77*data$Endocrine + 0.11*data$uppergi + 0.54*data$InnerEar + 0.11*data$Intestinal + 0.69*data$spinal + 0.79*data$psychiatric + 1.16*data$Nutritional + 1.06*data$Neurological + 0.49*data$Osteoporosis + 1.03*data$quadriplegia + 1.20*data$PeripheralNerve + 0.95*data$Pulmonarycir + 1.41*data$Sleep + 0.12*data$Urinary + 0.70*data$Venous + 0.17*data$OtherMAL + 0.39*data$Chronichepatitis + 1.04*data$Epilepsy

  # All previous ICD-10 coding algorithms were refered to:
  # Sarfati, D., Developing new comorbidity indices for cancer populations using administrative data. 2013, University of Otago: Dunedin. <Dr. Sarfati's doctoral dissertation>

  return(data)

}
