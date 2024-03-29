% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cci.R
\name{cci}
\alias{cci}
\title{cci()}
\usage{
cci(data, comorbidity, age)
}
\arguments{
\item{data}{Your data file in which Charlson Comorbidity index is to be calculated}

\item{comorbidity}{A vector of all comorbidity variables}

\item{age}{The name of the age variable}
}
\value{
data: a new data.frame named "data". This data frame contains two new variables: "cci_1987" & "cci_2011"
}
\description{
This file aims to calculate the Charlson Comorbidity index(1985 orginal version and the 2011 Quan version)
}
\details{
To calculate the Charlson Comorbidity index ("cci_1987" 1987 orginal version and cci_2011 the 2011 Quan version)
"cci_1987": The Charlson Comorbidity index, developed by Mary E. Charlson in 1987.
"cci_2011": The Charlson Comorbidity index, updated by Hude Quan in 2011.
}
\references{
Charlson, M. E., Pompei, P., Ales, K. L., & MacKenzie, C. R. (1987). A new method of classifying prognostic comorbidity in longitudinal studies: development and validation. Journal of chronic diseases, 40(5), 373-383.

Quan, H., Sundararajan, V., Halfon, P., Fong, A., Burnand, B., Luthi, J. C., ... & Ghali, W. A. (2005). Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data. Medical care, 1130-1139.

Quan, H., Li, B., Couris, C. M., Fushimi, K., Graham, P., Hider, P., ... & Sundararajan, V. (2011). Updating and validating the Charlson comorbidity index and score for risk adjustment in hospital discharge abstracts using data from 6 countries. American journal of epidemiology, 173(6), 676-682.
}
\author{
Miao Cai <email: miao.cai@outlook.com>
}
