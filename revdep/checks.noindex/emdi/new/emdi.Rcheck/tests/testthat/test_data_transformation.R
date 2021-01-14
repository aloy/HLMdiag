# Test the data transformation functions


# The data that is used for testing is the data from the sae package. 
load("EBP/incomedata.RData")
load("EBP/incomedata_woTeruel.RData")
load("EBP/Xoutsamp_AuxVar.RData")


# Test if return is a data.frame
test_that("Test that transformed_data is data_frame", {
  expect_equal(class(data_transformation(income~educ1, smp_data=incomedata, 
                                         lambda=NULL,
                            transformation="no")$transformed_data),
            "data.frame")
  expect_equal(class(std_data_transformation(income~educ1, incomedata, lambda=NULL,
                                transformation="no")),
            "data.frame")
  expect_equal(class(data_transformation(income~educ1, smp_data=incomedata, lambda=NULL,
                               transformation="log")$transformed_data),
               "data.frame")
  expect_equal(class(std_data_transformation(income~educ1, incomedata, lambda=NULL,
                                   transformation="log")),
               "data.frame")
  expect_equal(class(data_transformation(income~educ1, smp_data=incomedata, lambda=1,
                               transformation="box.cox")$transformed_data),
               "data.frame")
  expect_equal(class(std_data_transformation(income~educ1, incomedata, lambda=1,
                                   transformation="box.cox")),
               "data.frame")

})


# Test if the data frame returns 9 variables or rather a variable y

# Load benchmark transformed data or rather y vector and lambda
data_bc <- read.csv2("EBP/data_bc.csv", sep = ",", stringsAsFactors = TRUE)
data_bc_std <- read.csv2("EBP/data_bc_std.csv", sep = ",", stringsAsFactors = TRUE)
data_log <- read.csv2("EBP/data_log.csv", sep = ",", stringsAsFactors = TRUE)
ebp_optpar_bc <- read.csv2("EBP/ebp_optpar_bc.csv", sep = ",", stringsAsFactors = TRUE)

test_that("Test if data_transformation returns correctly 
          transformed y and a correct shift paramter",{

  lambda <- as.numeric(as.character(ebp_optpar_bc[,"Optpar"]))
  
  # Compare benchmark data with transformed y-vector and shift parameter
  # by data_transformation
  
  # Box-Cox transformation
  transformed_y_bc <-  data_transformation(income~educ1, 
                                    smp_data=incomedata, 
                                    lambda=lambda,
                                    transformation="box.cox")$transformed_data$income
  transformed_y_bc <- as.vector(transformed_y_bc)
  
  shift_bc <-  data_transformation(income~educ1, 
                                   smp_data=incomedata, 
                                   lambda=lambda,
                                   transformation="box.cox")$shift
  
  std_transformed_y_bc <-  as.vector(std_data_transformation(income~educ1, 
                                           smp_data=incomedata, 
                                           lambda=lambda,
                                           transformation="box.cox")$income)
  
  expect_equal(transformed_y_bc, as.numeric(as.character(data_bc$y)))
  expect_equal(shift_bc, data_bc$m[1])
  expect_equal(std_transformed_y_bc, as.numeric(as.character(data_bc_std[,1])))
  
  # Log transformation, lambda is NULL for log transformation
  transformed_y_log <- as.vector(data_transformation(income~educ1, 
                                                     smp_data=incomedata, 
                                                     lambda=NULL,
                                                     transformation="log")$transformed_data$income)
  shift_log <- data_transformation(income~educ1, 
                                   smp_data=incomedata, 
                                   lambda=NULL,
                                   transformation="log")$shift
  
  expect_equal(transformed_y_log, as.numeric(as.character(data_log$y)))
  expect_equal(shift_log, data_log$m[1])
  
  # No transformation, lambda is NULL, transformed data needs to be the same as
  # original data
  transformed_y_no <- as.vector(data_transformation(income~educ1, 
                                                    smp_data=incomedata, 
                                                    lambda=NULL,
                                                    transformation="no")$transformed_data$income)
  shift_no <- data_transformation(income~educ1, 
                                  smp_data=incomedata, 
                                  lambda=NULL,
                                  transformation="no")$shift
  
  expect_equal(transformed_y_no, incomedata$income)
  expect_equal(shift_no, NULL)

})






# Test back transformations

test_that("Does back transformation gives sample value?", {
  
  # Load benchmark transformed data or rather y vector and lambda

  lambda <- as.numeric(as.character(ebp_optpar_bc[,"Optpar"]))
  
  
  back_trans_no <- back_transformation(y=incomedata$income, 
                                       transformation="no",
                                       lambda=NULL,
                                       shift=NULL)
  expect_equal(back_trans_no, incomedata$income)
  

  back_trans_log <- back_transformation(y=as.numeric(as.character(data_log$y)), 
                                        transformation="log",
                                        lambda=NULL,
                                        shift=data_log$m[1])
  expect_equal(back_trans_log, incomedata$income)
  
  
  
  
  back_trans_bc <- back_transformation(y=as.numeric(as.character(data_bc$y)), 
                                       transformation="box.cox",
                                       lambda=lambda,
                                       shift=data_bc$m[1])
  expect_equal(back_trans_bc, incomedata$income)

})



