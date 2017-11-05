# MiaoCom
[![Build status](https://ci.appveyor.com/api/projects/status/00t1917ftctbv0p7/branch/miao?svg=true)](https://ci.appveyor.com/project/caimiao0714/miaocom/branch/miao)

## Overview
A R package for calculating comorbidity indexes (i.e. [Charlson Comorbidity Index](https://en.wikipedia.org/wiki/Comorbidity#Charlson_index), [Elixhauser Comorbidity Index](https://en.wikipedia.org/wiki/Comorbidity#Elixhauser_comorbidity_measure), and [C3 Index](https://www.ncbi.nlm.nih.gov/pubmed/24582212)) for epidemiologists.


## Installation
    install.packages("devtools")
    devtools::install_github("caimiao0714/MiaoCom@miao") 
    # The master branch is just a initialized branch. The miao branch is the branch that I really work on.



## Overview of the functions
    library(MiaoCom)


*  `cci()` is a function that calculates [Charlson Comorbidity Index](https://en.wikipedia.org/wiki/Comorbidity#Charlson_index) for patients. The ICD-10 coding algorithms for defining comorbidities refers to [Hude Quan's paper](http://www.jstor.org/stable/3768193?seq=1#page_scan_tab_contents) in 2005. This function returns two Charlson Comorbidity Indexes:
    1. [the 1987 version](http://www.sciencedirect.com/science/article/pii/0021968187901718) developed by Mary E. Charlson; 
    2. [the 2011 updated version](https://academic.oup.com/aje/article/173/6/676/182985/Updating-and-Validating-the-Charlson-Comorbidity) by Hude Quan. 

* `eci()` is a function that calculates [Elixhauser Comorbidity Index](https://en.wikipedia.org/wiki/Comorbidity#Elixhauser_comorbidity_measure) for patients. The ICD-10 coding algorithms for defining comorbidities also refers to [Hude Quan's paper](http://www.jstor.org/stable/3768193?seq=1#page_scan_tab_contents) in 2005.
*  `c3()` is a function that calculates [C3 Index](https://www.ncbi.nlm.nih.gov/pubmed/24582212) specifically for cancer patients



## Usage of functions
    new_data <- cci(data, comorbidity, age)
    new_data <- eci(data, comorbidity, age)
    new_data <-  c3(data, comorbidity, age)
* `data` is a data.frame  or alike objects from which you want to calculate comorbidity indexed from.
* `comorbidity` is a vector of the variable names of patients' comorbidity ICD-10 codings in the data. For example, `c("comorbidity1", "comorbidity2","comorbidity3","comorbidity4","comorbidity5")`.
* `age` is the variable name of the age of the patients in your data, the name should **NOT** be surrounded by single or double quotes (i.e. `'`, `"`).



## A simple example
> Step 1:　Construct a demo dataframe.

    demo_data <- data.frame(
       comorbidityICD1 = c("I25.105", "I50.907", "I25.903", "I50.907", "I50.907", "I25.903",  
                           "I50.907", "I50.910", "I25.105", "I25.903"),
       comorbidityICD2 = c("I25.210", "I25.903", "I50.908", "J44.103", "I10.003", "I50.907", 
                           "K21.001", "J44.101", "I25.203", "I50.907"),
       comorbidityICD3 = c("I51.709", "K27.906", "E11.732", "I10.005", "E11.901", "E78.501", 
                           "Z98.8123", "E11.901", "I50.908", "E11.901"),
       comorbidityICD4 = c("I49.904", "I10.005", "J44.101", "K73.901", "I63.902", "J98.402", 
                           "B18.106", "E11.221", "E11.423", "E11.2211"),
       comorbidityICD5 = c("I45.101", "E11.901", "J96.903", "K70.301", "K76.811", "K72.901", 
                           "K74.602", "I15.002", "J44.003", "N18.916"),
       patient_age = c(72, 66, 81, 86, 75, 33, 63, 70, 71, 47)
    )


> Step 2: Apply `cci()` function to the demo data
    
    newdata1 <- cci(data = demo_data, comorbidity = c("comorbidityICD1", "comorbidityICD2", 
                    "comorbidityICD3", "comorbidityICD4", "comorbidityICD5"), age = patient_age)

The generated `newdata1` data.frame includes multiple new columns. The `CCI_1987` and `CCI_2011` columns deserve your special attention.

1. `CCI_1987` is  [the 1987 version CCI](http://www.sciencedirect.com/science/article/pii/0021968187901718) developed by Mary E. Charlson;
2. `CCI_2011` is [the 2011 updated version CCI](https://academic.oup.com/aje/article/173/6/676/182985/Updating-and-Validating-the-Charlson-Comorbidity) by Hude Quan.


> Step 3: Apply `eci()` function to the demo data
    
    newdata2 <- eci(data = demo_data, comorbidity = c("comorbidityICD1", "comorbidityICD2", 
                    "comorbidityICD3", "comorbidityICD4", "comorbidityICD5"))
The generated `newdata2` data.frame includes multiple new columns. The `CCI_1987` is your desired [Elixhauser Comorbidity Index](https://en.wikipedia.org/wiki/Comorbidity#Elixhauser_comorbidity_measure).


> Step 4: Apply `c3()` function to the demo data
    
    newdata3 <- c3(data = demo_data, comorbidity = c("comorbidityICD1", "comorbidityICD2", 
                    "comorbidityICD3", "comorbidityICD4", "comorbidityICD5"))
The generated `newdata3` data.frame includes multiple new columns. The `C3` is your desired [C3 Index](https://www.ncbi.nlm.nih.gov/pubmed/24582212). Please note that the C3 index is developed specifically for cancer patients. This is just a demonstration of the `c3()` function.

Note that CCI includes patient age as a variable for calculating the index. However, patient age is not included in calculating Elixhauser Comorbidity Index and C3 Index. 
