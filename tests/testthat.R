library(testthat)
library(MiaoCom)

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


test_that("a test for the CCI function: ", {
  expect_equal(MiaoCom::cci(data = demo_data, comorbidity = c("comorbidityICD1", "comorbidityICD2",
                                                     "comorbidityICD3", "comorbidityICD4", "comorbidityICD5"), age = "patient_age")[,"CCI_2011"], c(3, 4, 8, 9, 7, 6, 6, 7, 7, 4))
})
