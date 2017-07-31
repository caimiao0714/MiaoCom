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
#' @param x A vector of all comorbidity variables
#' @param data_file Your data file in which Charlson Comorbidity index is to be calculated
#' @return cci_1987: The Charlson Comorbidity index, developed by Mary E. Charlson in 1987
#' @return cci_2011: The Charlson Comorbidity index, updated by Hude Quan in 2011
#' @export

cci <- function(x, data_file) {
  ## subsetting data
  library(dplyr)
  data_comorbidity <- subset(data_file, select = x)# subset the whole dataset into "data_comorbidity": subset data that only contains comorbidities
  data_comorbidity %>% mutate_if(is.factor, as.character) -> data_comorbidity # this converts factor values into characters
  data_comorbidity[is.na(data_comorbidity)] <- 0 #imputate missing values with zeros
  dim_comorbidity <- dim(data_comorbidity)# save the dimensionality of comorbidity subset data
  unlisted_data <- unlist(data_comorbidity)# unlist the comorbidity subset data(converts the data into a vector. Vectorization speeds up the functions)

  ##1:MI, myocardial infarction
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  MIchk <- function(x){
    grepl("I25\\.2|I21\\.|I22\\.",x)
  }
  test1 <- sapply(test,MIchk)
  J00_J47$MI <- rowSums(test1)
  J00_J47$MI[J00_J47$MI >=1] <-1

  bidata_com <- (substr(unlisted_data,1,5) == "I25.2") | (substr(unlisted_data, 1, 4) %in% c("I21.","I22.")) # regular expression to find out cases with "myocardial infarction" comorbidity ICD-10 codes
  dim(bidata_com) <- dim_comorbidity # convert
  data_file$MI <- rowSums(bidata_com)


  J00_J47$CHF_Elix <- rowSums(res)
  J00_J47$CHF_Elix[J00_J47$CHF_Elix >=1] <-1


  #2:CCF, congestive cardiac failure
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  CCFchk <- function(x){
    grepl("I09\\.9|I11\\.0|I13\\.0|I13\\.2|I25\\.5|I42\\.0|I42\\.5|I42\\.6|I42\\.7|I42\\.8|I42\\.9|P29\\.0|I50\\.|I43\\.",x)
  }
  test1 <- sapply(test,CCFchk)
  J00_J47$CCF <- rowSums(test1)
  J00_J47$CCF[J00_J47$CCF >=1] <-1

  # 合并症筛选3-周围性血管疾病PVD, peripheral vascular disease
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  PVDchk <- function(x){
    grepl("I73\\.1|I73\\.8|I73\\.9|I77\\.1|I79\\.0|I79\\.2|K55\\.1|K55\\.8|K55\\.9|Z95\\.8|Z95\\.9|I70\\.|I71",x)
  }
  test1 <- sapply(test,PVDchk)
  J00_J47$PVD <- rowSums(test1)
  J00_J47$PVD[J00_J47$PVD >=1] <-1

  # 合并症筛选4-脑血管疾病 CD, Cerebrovascular disease
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  CDchk <- function(x){
    grepl("H34\\.0|G45|G46|I60|I61|I62|I63|I64|I65|I66|I67|I68|I69",x)
  }
  test1 <- sapply(test,CDchk)
  J00_J47$CD <- rowSums(test1)
  J00_J47$CD[J00_J47$CD >=1] <-1

  # 合并症筛选5-痴呆Dementia
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  Dementiachk <- function(x){
    grepl("F05\\.1|G31\\.1|F00|F01|F02|F03|G30",x)
  }
  test1 <- sapply(test,Dementiachk)
  J00_J47$Dementia <- rowSums(test1)
  J00_J47$Dementia[J00_J47$Dementia >=1] <-1

  # 合并症筛选6 - 慢性阻塞性肺病COPD, chronic obstructive pulmonary disease
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  COPDchk <- function(x){
    grepl("I27\\.8|I27\\.9|J68\\.4|J70\\.1|J70\\.3|J40|J41|J42|J43|J44|J45|J46|J47|J60|J61|J62|J63|J64|J65|J66|J67",x)
  }
  test1 <- sapply(test,COPDchk)
  J00_J47$COPD <- rowSums(test1)
  J00_J47$COPD[J00_J47$COPD >=1] <-1

  # 合并症筛选7 -结缔组织病 RD - Connective tissue disease
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  RDchk <- function(x){
    grepl("M31\\.5|M35\\.1|M35\\.3|M36\\.0|M32|M33|M34|M05|M06",x)
  }
  test1 <- sapply(test,RDchk)
  J00_J47$RD <- rowSums(test1)
  J00_J47$RD[J00_J47$RD >=1] <-1

  # 合并症筛选8 –溃疡PUD, Ulcers
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  PUDchk <- function(x){
    grepl("K25|K26|K27|K28",x)
  }
  test1 <- sapply(test,PUDchk)
  J00_J47$PUD <- rowSums(test1)
  J00_J47$PUD[J00_J47$PUD >=1] <-1

  # 合并症筛选9-轻微的肝脏疾病MLD, Mild liver disease
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  MLDchk <- function(x){
    grepl("K70\\.0|K70\\.1|K70\\.2|K70\\.3|K70\\.9|K71\\.3|K71\\.4|K71\\.5|K71\\.7|K76\\.0|K76\\.2|K76\\.3|K76\\.4|K76\\.8|K76\\.9|Z94\\.4|B18|K73|K74",x)
  }
  test1 <- sapply(test,MLDchk)
  J00_J47$MLD <- rowSums(test1)
  J00_J47$MLD[J00_J47$MLD >=1] <-1

  # 合并症筛选10 - 糖尿病（有终末器官损害）DMandEOD (with end-organ damage)
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  DMandEODchk <- function(x){
    grepl("E10\\.2|E10\\.3|E10\\.4|E10\\.5|E10\\.7|E11\\.2|E11\\.3|E11\\.4|E11\\.5|E11\\.7|E12\\.2|E12\\.3|E12\\.4|E12\\.5|E12\\.7|E13\\.2|E13\\.3|E13\\.4|E13\\.5|E13\\.7|E14\\.2|E14\\.3|E14\\.4|E14\\.5|E14\\.7",x)
  }
  test1 <- sapply(test,DMandEODchk)
  J00_J47$DMandEOD <- rowSums(test1)
  J00_J47$DMandEOD[J00_J47$DMandEOD >=1] <-1

  # 合并症筛选 11 -糖尿病（没有终末器官损害）DMnotEOD (without end-organ damage)
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  DMnotEODchk <- function(x){
    grepl("E10\\.0|E10\\.1|E10\\.6|E10\\.8|E10\\.9|E11\\.0|E11\\.1|E11\\.6|E11\\.8|E11\\.9|E12\\.0|E12\\.1|E12\\.6|E12\\.8|E12\\.9|E13\\.0|E13\\.1|E13\\.6|E13\\.8|E13\\.9|E14\\.0|E14\\.1|E14\\.6|E14\\.8|E14\\.9",x)
  }
  test1 <- sapply(test,DMnotEODchk)
  J00_J47$DMnotEOD <- rowSums(test1)
  J00_J47$DMnotEOD[J00_J47$DMnotEOD >=1] <-1

  # 合并症筛选 12 - 偏瘫Hemiplegia
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  Hemiplegiachk <- function(x){
    grepl("G04\\.1|G11\\.4|G80\\.1|G80\\.2|G83\\.0|G83\\.1|G83\\.2|G83\\.4|G83\\.9|G81|G82",x)
  }
  test1 <- sapply(test,Hemiplegiachk)
  J00_J47$Hemiplegia <- rowSums(test1)
  J00_J47$Hemiplegia[J00_J47$Hemiplegia >=1] <-1

  # 合并症筛选 13 - 中度到重度的慢性肾脏疾病MSCKD, Moderate to Severe Chronic Kidney Disease
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  MSCKDchk <- function(x){
    grepl("I12\\.0|I13\\.1|N03\\.2|N03\\.3|N03\\.4|N03\\.5|N03\\.6|N03\\.7|N05\\.2|N05\\.3|N05\\.4|N05\\.5|N05\\.6|N05\\.7|N25\\.0|Z49\\.0|Z49\\.1|Z49\\.2|Z94\\.0|Z99\\.2|N18|N19",x)
  }
  test1 <- sapply(test,MSCKDchk)
  J00_J47$MSCKD <- rowSums(test1)
  J00_J47$MSCKD[J00_J47$MSCKD >=1] <-1

  # 合并症筛选 14 - 恶性（肿瘤等）Malignancy
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  Malignancychk <- function(x){
    grepl("C00|C01|C02|C03|C04|C05|C06|C07|C08|C09|C10|C11|C12|C13|C14|C15|C16|C17|C18|C19|C20|C21|C22|C23|C24|C25|C26|C30|C31|C32|C33|C34|C37|C38|C39|C40|C41|C43|C45|C46|C47|C48|C49|C50|C51|C52|C53|C54|C55|C56|C57|C58|C60|C61|C62|C63|C64|C65|C66|C67|C68|C69|C70|C71|C72|C73|C74|C75|C76|C81|C82|C83|C84|C85|C88|C90|C91|C92|C93|C94|C95|C96|C97",x)
  }
  test1 <- sapply(test,Malignancychk)
  J00_J47$Malignancy <- rowSums(test1)
  J00_J47$Malignancy[J00_J47$Malignancy >=1] <-1

  # 合并症筛选 15 - 中等到重度肝脏疾病MSLD, Moderate–severe liver disease
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  MSLDchk <- function(x){
    grepl("I85\\.0|I85\\.9|I86\\.4|I98\\.2|K70\\.4|K71\\.1|K72\\.1|K72\\.9|K76\\.5|K76\\.6|K76\\.7",x)
  }
  test1 <- sapply(test,MSLDchk)
  J00_J47$MSLD <- rowSums(test1)
  J00_J47$MSLD[J00_J47$MSLD >=1] <-1

  # 合并症筛选 16 - 转移固体肿瘤MST, Metastatic solid tumour
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  MSTchk <- function(x){
    grepl("C77|C78|C79|C80",x)
  }
  test1 <- sapply(test,MSTchk)
  J00_J47$MST <- rowSums(test1)
  J00_J47$MST[J00_J47$MST >=1] <-1

  # 合并症筛选 17 - 艾滋病AIDS
  test <- subset(J00_J47, select = c("p324", "p327", "p3291", "p3294", "p3297", "p3281", "p3284", "p3287", "p3271", "p3274"))
  AIDSchk <- function(x){
    grepl("B20|B21|B22|B24",x)
  }
  test1 <- sapply(test,AIDSchk)
  J00_J47$AIDS <- rowSums(test1)
  J00_J47$AIDS[J00_J47$AIDS >=1] <-1

  #PART B---AGE
  J00_J47$年龄G[J00_J47$p7 < 50] <- 0
  J00_J47$年龄G[J00_J47$p7 <= 59 & J00_J47$p7 >= 50] <- 1
  J00_J47$年龄G[J00_J47$p7 <= 69 & J00_J47$p7 >= 60] <- 2
  J00_J47$年龄G[J00_J47$p7 >= 70] <- 3

  # 原始CCI的计算
  J00_J47 <- within(J00_J47, {
    CCI_1987 <-
      MI + CCF + PVD + CD + Dementia + COPD + RD + PUD + MLD + DMnotEOD + DMandEOD *
      2 + Hemiplegia * 2 + MSCKD * 2 + Malignancy * 2 + MSLD * 3 + MST * 6 + AIDS *
      6 + 年龄G
  })

  #2011 Quan 等人调整的CCI
  J00_J47$CCI_2011 <-
    J00_J47$CCF * 2 + J00_J47$Dementia * 2 + J00_J47$COPD + J00_J47$RD + J00_J47$MLD *
    2 + J00_J47$DMandEOD + J00_J47$Hemiplegia * 2 + J00_J47$MSCKD + J00_J47$Malignancy *
    2 + J00_J47$MSLD * 4 + J00_J47$MST * 6 + J00_J47$AIDS * 4 + J00_J47$年龄G
  end.time <- Sys.time()
  end.time - start.time
  table(J00_J47$CCI_2011,exclude = NULL)
  table(J00_J47$CCI_1987,exclude = NULL)
}
