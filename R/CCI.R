# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' To calculate the Charlson Comorbidity index(1985 orginal version and the 2011 Quan version)
#' @author Miao Cai <email: miaocai@slu.edu>
#' @description This file aims to calculate the Charlson Comorbidity index(1985 orginal version and the 2011 Quan version)
#' @param data_file Your data file in which Charlson Comorbidity index is to be calculated
#' @param x A vector of all comorbidity variables
#' @param age The name of the age variable
#' @return cci_1987: The Charlson Comorbidity index, developed by Mary E. Charlson in 1987
#' @return cci_2011: The Charlson Comorbidity index, updated by Hude Quan in 2011
#' @export

cci <- function(data_file, x, age) {
  start.time <- Sys.time()

  ### PART A: SUBSETTING DATA
  suppressMessages(library(dplyr))
  data_comorbidity <- subset(data_file, select = x)# subset the whole dataset into "data_comorbidity": subset data that only contains comorbidities
  data_comorbidity %>% mutate_if(is.factor, as.character) -> data_comorbidity # this converts factor values into characters
  data_comorbidity[is.na(data_comorbidity)] <- 0 #imputate missing values with zeros
  dim_comorbidity <- dim(data_comorbidity)# save the dimensionality of comorbidity subset data
  unlisted_data <- unlist(data_comorbidity)# unlist the comorbidity subset data(converts the data into a vector. Vectorization speeds up the functions)

### PART B: FILTERING COMORBIDITIES
  ##1:MI, Myocardial Infarction
  bidata_com <- (substr(unlisted_data,1,5) == "I25.2") | (substr(unlisted_data, 1, 4) %in% c("I21.","I22.")) # regular expression to find out cases with "Myocardial Infarction" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$MI <- rowSums(bidata_com)
  data_file$MI[data_file$MI >= 1] <- 1

  #2:CCF, Congestive Cardiac Failure
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I09.9", "I11.0", "I13.0", "I13.2", "I25.5", "I42.0", "I42.5", "I42.6", "I42.7", "I42.8", "I42.9", "P29.0")) | (substr(unlisted_data, 1, 4) %in% c("I50.","I43.")) # regular expression to find out cases with "Congestive Cardiac Failure" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$CCF <- rowSums(bidata_com)
  data_file$CCF[data_file$CCF >= 1] <- 1

  #3:PVD, Peripheral Vascular Disease
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I73.1", "I73.8", "I73.9", "I77.1", "I79.0", "I79.2", "K55.1", "K55.8", "K55.9", "Z95.8", "Z95.9")) | (substr(unlisted_data, 1, 4) %in% c("I70.", "I71.")) # regular expression to find out cases with "Peripheral Vascular Disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$PVD <- rowSums(bidata_com)
  data_file$PVD[data_file$PVD >= 1] <- 1

  #4:CD, Cerebrovascular Disease
  bidata_com <- (substr(unlisted_data,1,5) == "H34.0") | (substr(unlisted_data, 1, 4) %in% c("I60.", "I61.", "I62.", "I63.", "I64.", "I65.", "I66.", "I67.", "I68.", "I69")) # regular expression to find out cases with "Cerebrovascular Disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$CD <- rowSums(bidata_com)
  data_file$CD[data_file$CD >= 1] <- 1

  #5:Dementia
  bidata_com <- (substr(unlisted_data,1,5) %in% c("F05.1", "G31.1")) | (substr(unlisted_data, 1, 4) %in% c("F00.", "F01.", "F02.", "F03.", "G30.")) # regular expression to find out cases with "Dementia" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Dementia <- rowSums(bidata_com)
  data_file$Dementia[data_file$Dementia >= 1] <- 1

  #6:COPD, chronic obstructive pulmonary disease
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I27.8", "I27.9", "J68.4", "J70.1", "J70.3")) | (substr(unlisted_data, 1, 4) %in% c("J40.", "J41.", "J42.", "J43.", "J44.", "J45.", "J46.", "J47.", "J60.", "J61.", "J62.", "J63.", "J64.", "J65.", "J66.", "J67.")) # regular expression to find out cases with "COPD" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$COPD <- rowSums(bidata_com)
  data_file$COPD[data_file$COPD >= 1] <- 1

  #7:RD - Connective tissue disease
  bidata_com <- (substr(unlisted_data,1,5) %in% c("M31.5", "M35.1", "M35.3", "M36.0")) | (substr(unlisted_data, 1, 4) %in% c("M32.", "M33.", "M34.", "M05.", "M06.")) # regular expression to find out cases with "RD" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$RD <- rowSums(bidata_com)
  data_file$RD[data_file$RD >= 1] <- 1

  #8:PUD, Ulcers
  bidata_com <- substr(unlisted_data, 1, 4) %in% c("K25.", "K26.", "K27.", "K28.") # regular expression to find out cases with "PUD" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$PUD <- rowSums(bidata_com)
  data_file$PUD[data_file$PUD >= 1] <- 1

  #9:MLD, Mild liver disease
  bidata_com <- (substr(unlisted_data,1,5) %in% c("K70.0", "K70.1", "K70.2", "K70.3", "K70.9", "K71.3", "K71.4", "K71.5", "K71.7", "K76.0", "K76.2", "K76.3", "K76.4", "K76.8", "K76.9", "Z94.4")) | (substr(unlisted_data, 1, 4) %in% c("B18.", "K73.", "K74.")) # regular expression to find out cases with "MLD" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$MLD <- rowSums(bidata_com)
  data_file$MLD[data_file$MLD >= 1] <- 1

  #10:DMandEOD (with end-organ damage)
  bidata_com <- substr(unlisted_data, 1, 5) %in% c("E10.2", "E10.3", "E10.4", "E10.5", "E10.7", "E11.2", "E11.3", "E11.4", "E11.5", "E11.7", "E12.2", "E12.3", "E12.4", "E12.5", "E12.7", "E13.2", "E13.3", "E13.4", "E13.5", "E13.7", "E14.2", "E14.3", "E14.4", "E14.5", "E14.7") # regular expression to find out cases with "DMandEOD (with end-organ damage)" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$DMandEOD <- rowSums(bidata_com)
  data_file$DMandEOD[data_file$DMandEOD >= 1] <- 1

  #11:DMnotEOD (without end-organ damage)
  bidata_com <- substr(unlisted_data, 1, 5) %in% c("E10.0", "E10.1", "E10.6", "E10.8", "E10.9", "E11.0", "E11.1", "E11.6", "E11.8", "E11.9", "E12.0", "E12.1", "E12.6", "E12.8", "E12.9", "E13.0", "E13.1", "E13.6", "E13.8", "E13.9", "E14.0", "E14.1", "E14.6", "E14.8", "E14.9") # regular expression to find out cases with "DMnotEOD (without end-organ damage)" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$DMnotEOD <- rowSums(bidata_com)
  data_file$DMnotEOD[data_file$DMnotEOD >= 1] <- 1

  #12: Hemiplegia
  bidata_com <- (substr(unlisted_data,1,5) %in% c("G04.1", "G11.4", "G80.1", "G80.2", "G83.0", "G83.1", "G83.2", "G83.4", "G83.9")) | (substr(unlisted_data, 1, 4) %in% c("G81.", "G82.")) # regular expression to find out cases with "Hemiplegia" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Hemiplegia <- rowSums(bidata_com)
  data_file$Hemiplegia[data_file$Hemiplegia >= 1] <- 1

  #13:MSCKD, Moderate to Severe Chronic Kidney Disease
  bidata_com <- (substr(unlisted_data,1,5) %in% c("I12.0", "I13.1", "N03.2", "N03.3", "N03.4", "N03.5", "N03.6", "N03.7", "N05.2", "N05.3", "N05.4", "N05.5", "N05.6", "N05.7", "N25.0", "Z49.0", "Z49.1", "Z49.2", "Z94.0", "Z99.2")) | (substr(unlisted_data, 1, 4) %in% c("N18.", "N19.")) # regular expression to find out cases with "MSCKD" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$MSCKD <- rowSums(bidata_com)
  data_file$MSCKD[data_file$MSCKD >= 1] <- 1

  #14:Malignancy, Maligant Tumor
  bidata_com <- substr(unlisted_data, 1, 4) %in% c("C00.", "C01.", "C02.", "C03.", "C04.", "C05.", "C06.", "C07.", "C08.", "C09.", "C10.", "C11.", "C12.", "C13.", "C14.", "C15.", "C16.", "C17.", "C18.", "C19.", "C20.", "C21.", "C22.", "C23.", "C24.", "C25.", "C26.", "C30.", "C31.", "C32.", "C33.", "C34.", "C37.", "C38.", "C39.", "C40.", "C41.", "C43.", "C45.", "C46.", "C47.", "C48.", "C49.", "C50.", "C51.", "C52.", "C53.", "C54.", "C55.", "C56.", "C57.", "C58.", "C60.", "C61.", "C62.", "C63.", "C64.", "C65.", "C66.", "C67.", "C68.", "C69.", "C70.", "C71.", "C72.", "C73.", "C74.", "C75.", "C76.", "C81.", "C82.", "C83.", "C84.", "C85.", "C88.", "C90.", "C91.", "C92.", "C93.", "C94.", "C95.", "C96.", "C97") # regular expression to find out cases with "Maligant Tumor" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$Malignancy <- rowSums(bidata_com)
  data_file$Malignancy[data_file$Malignancy >= 1] <- 1

  #15:MSLD, Moderate–severe liver disease
  bidata_com <- substr(unlisted_data, 1, 5) %in% c("I85.0", "I85.9", "I86.4", "I98.2", "K70.4", "K71.1", "K72.1", "K72.9", "K76.5", "K76.6", "K76.7") # regular expression to find out cases with "Moderate–severe liver disease" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$MSLD <- rowSums(bidata_com)
  data_file$MSLD[data_file$MSLD >= 1] <- 1

  #16:MST, Metastatic solid tumour
  bidata_com <- substr(unlisted_data, 1, 4) %in% c("C77.", "C78.", "C79.", "C80.") # regular expression to find out cases with "Metastatic solid tumour" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$MST <- rowSums(bidata_com)
  data_file$MST[data_file$MST >= 1] <- 1

  #17: AIDS
  bidata_com <- substr(unlisted_data, 1, 4) %in% c("B20.", "B21.", "B22.", "B24.") # regular expression to find out cases with "AIDS" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert the comorbidity data vector back into data.frame
  data_file$AIDS <- rowSums(bidata_com)
  data_file$AIDS[data_file$AIDS >= 1] <- 1

  ###PART C---AGE
  # About age: There are a variety of calculating weights for patient ages. In the orginal 1987 paper, it wrote
  # "Using this approach, a patient 40 yr of age would be assumed to have no risk of comorbid death attributable to age and a patient with a comorbidity index score of 0 would have no risk attributable to pre-existing comorbid disease. Each decade of  age over 40 would add 1 point to risk (i.e. 50 yr, 1; 60 yr, 2; 70 yr 3; etc.) and the "age points" would be added to the score from the  comorbidity index (i.e. 0, 1, 2, 3, etc.)."
  # Therefore, we use this approach for my algorithm
  data_file$age_group <- 0
  data_file$age_group[data_file$age <= 59 & data_file$age >= 50] <- 1
  data_file$age_group[data_file$age <= 69 & data_file$age >= 60] <- 2
  data_file$age_group[data_file$age <= 79 & data_file$age >= 70] <- 3
  data_file$age_group[data_file$age >= 80] <- 4

  ###PART D--WEIGHTED SUM
  #1:CCI_1987, by Mary E. Charlson in 1987
  # Reference: Charlson, M. E., Pompei, P., Ales, K. L., & MacKenzie, C. R. (1987). A new method of classifying prognostic comorbidity in longitudinal studies: development and validation. Journal of chronic diseases, 40(5), 373-383.

  data_file <- within(data_file, {
    CCI_1987 <-
      MI + CCF + PVD + CD + Dementia + COPD + RD + PUD + MLD + DMnotEOD + DMandEOD *
      2 + Hemiplegia * 2 + MSCKD * 2 + Malignancy * 2 + MSLD * 3 + MST * 6 + AIDS *
      6 + age_group
  })

  #2:CCI_2011, weights were updated by Hude Quan et al.
  #Reference:
  data_file$CCI_2011 <- data_file$CCF * 2 + data_file$Dementia * 2 + data_file$COPD + data_file$RD + data_file$MLD *
    2 + data_file$DMandEOD + data_file$Hemiplegia * 2 + data_file$MSCKD + data_file$Malignancy *
    2 + data_file$MSLD * 4 + data_file$MST * 6 + data_file$AIDS * 4 + data_file$age_group

  # All previous ICD-10 coding algorithms were refered to:
  # Quan, H., Sundararajan, V., Halfon, P., Fong, A., Burnand, B., Luthi, J. C., ... & Ghali, W. A. (2005). Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data. Medical care, 1130-1139.

  end.time <- Sys.time()

  print("Total time cost of this function is: ")
  return(print(end.time - start.time))
  }
