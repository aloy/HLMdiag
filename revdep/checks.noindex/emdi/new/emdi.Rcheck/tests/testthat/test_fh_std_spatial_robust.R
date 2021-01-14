# Test if the same variance, EBLUP and MSE results are obtained as with 
# the sae package


########################## Standard Fay-Herriot model ##########################

# The data that is used for testing is the data from the sae package. 
load("FH/milk.RData")

test_that("Does the fh function in emdi return the same variance, EBLUP and MSE 
          estimates as the function mseFH of package sae?",{
            
            ############################ REML variance estimation ########################
            # Estimation with fh of emdi
            milk$var <- milk$SD^2
            fh_reml <- fh(fixed = yi ~ as.factor(MajorArea), vardir = "var",
                          combined_data = milk, domains = "SmallArea",
                          method = "reml", interval = c(0, 1000), MSE = TRUE, 
                          B = NULL)
            
            # Estimation with mseFH of sae (benchmark)
            fh_reml_sae <- read.csv("FH/fh_reml_sae.csv", sep = ",", 
                                    stringsAsFactors = TRUE) 
            
            # Comparison
            # Variance
            expect_equal(fh_reml$model$variance, fh_reml_sae$variance[1], tolerance = 0.0001)
            # EBLUP
            expect_equal(fh_reml$ind$FH, as.vector(fh_reml_sae$EBLUP), tolerance = 0.0001)
            # MSE
            expect_equal(fh_reml$MSE$FH, fh_reml_sae$MSE, tolerance = 0.00001)
            
            ############################ ML variance estimation ##########################
            # Estimation with fh of emdi
            fh_ml <- fh(fixed = yi ~ as.factor(MajorArea), vardir = "var",
                        combined_data = milk, domains = "SmallArea",
                        method = "ml", interval = c(0, 1000), MSE = TRUE, B = NULL)
            
            # Estimation with mseFH of sae (benchmark)
            fh_ml_sae <- read.csv("FH/fh_ml_sae.csv", sep = ",", 
                                  stringsAsFactors = TRUE) 
            
            # Comparison
            # Variance
            expect_equal(fh_ml$model$variance, fh_ml_sae$variance[1], tolerance = 0.0001)
            # EBLUP
            expect_equal(fh_ml$ind$FH, fh_ml_sae$EBLUP, tolerance = 0.0001)
            # MSE
            expect_equal(fh_ml$MSE$FH, fh_ml_sae$MSE, tolerance = 0.0001)
          })

# Test if the same variance, correlation parameter, EBLUP and MSE results are 
# obtained as with the sae package

########################## Spatial Fay-Herriot model ###########################
# The data that is used for testing is the data from the sae package. 
# Data set
load("FH/grapes.RData")
# Proximity matrix
load("FH/grapesprox.RData")

# Analytical MSE 

test_that("Does the fh function in emdi return the same variance, correlation 
           parameter, EBLUP and MSE estimates as the function mseSFH of package 
           sae?",{
             
             ############################ REML variance estimation ########################
             # Estimation with fh of emdi
             grapes$Domain <- c(1:274)
             fh_spatial_reml_analytical <- fh(fixed = grapehect ~ area + workdays - 1, 
                                              vardir = "var", tol = 0.0001, maxit = 100, 
                                              combined_data = grapes, domains = "Domain", 
                                              method = "reml", correlation = "spatial", 
                                              corMatrix = as.matrix(grapesprox), MSE = TRUE, 
                                              mse_type = "analytical")
             
             
             # Estimation with mseFH of sae (benchmark)
             fh_spatial_reml_analytical_sae <- read.csv("FH/fh_spatial_reml_analytical_sae.csv", 
                                                        sep = ",", 
                                                        stringsAsFactors = TRUE) 
             
             # Comparison
             # Variance
             expect_equal(fh_spatial_reml_analytical$model$variance[1,2], 
                          fh_spatial_reml_analytical_sae$variance[1])
             # Correlation parameter
             expect_equal(fh_spatial_reml_analytical$model$variance[1,1], 
                          fh_spatial_reml_analytical_sae$correlation[1])
             # EBLUP
             expect_equal(fh_spatial_reml_analytical$ind$FH, 
                          as.vector(fh_spatial_reml_analytical_sae$EBLUP))
             # MSE
             expect_equal(fh_spatial_reml_analytical$MSE$FH, 
                          fh_spatial_reml_analytical_sae$MSE)
             
             ############################ ML variance estimation ##########################
             # Estimation with fh of emdi
             fh_spatial_ml_analytical <- fh(fixed = grapehect ~ area + workdays - 1, 
                                            vardir = "var", tol = 0.0001, maxit = 100, 
                                            combined_data = grapes, domains = "Domain", 
                                            method = "ml", correlation = "spatial", 
                                            corMatrix = as.matrix(grapesprox), MSE = TRUE, 
                                            mse_type = "analytical")
             
             # Estimation with mseFH of sae (benchmark)
             fh_spatial_ml_analytical_sae <- read.csv("FH/fh_spatial_ml_analytical_sae.csv", 
                                                      sep = ",", 
                                                      stringsAsFactors = TRUE) 
             
             # Comparison
             # Variance
             expect_equal(fh_spatial_ml_analytical$model$variance[1,2], 
                          fh_spatial_ml_analytical_sae$variance[1])
             # Correlation parameter
             expect_equal(fh_spatial_ml_analytical$model$variance[1,1], 
                          fh_spatial_ml_analytical_sae$correlation[1])
             # EBLUP
             expect_equal(fh_spatial_ml_analytical$ind$FH, 
                          as.vector(fh_spatial_ml_analytical_sae$EBLUP))
             # MSE
             expect_equal(fh_spatial_ml_analytical$MSE$FH, 
                          fh_spatial_ml_analytical_sae$MSE)
           })

# Nonparametric bootstrap MSE 


test_that("Does the fh function in emdi return the same variance, correlation 
          parameter, EBLUP and MSE estimates as the function npbmseSFH of package 
          sae?",{
            suppressWarnings(RNGversion("3.6.3"))  
            ############################ REML variance estimation ########################
            # Estimation with fh of emdi
            grapes$Domain <- c(1:274)
            fh_spatial_reml_npb <- fh(fixed = grapehect ~ area + workdays - 1, 
                                      vardir = "var", tol = 0.0001, maxit = 100, 
                                      combined_data = grapes, domains = "Domain", 
                                      method = "reml", correlation = "spatial", 
                                      corMatrix = as.matrix(grapesprox), MSE = TRUE, 
                                      mse_type = "spatialnonparboot", B = 3, seed = 123)
            
            fh_spatial_reml_npb_bc <- fh(fixed = grapehect ~ area + workdays - 1, 
                                         vardir = "var", tol = 0.0001, maxit = 100, 
                                         combined_data = grapes, domains = "Domain", 
                                         method = "reml", correlation = "spatial", 
                                         corMatrix = as.matrix(grapesprox), MSE = TRUE, 
                                         mse_type = "spatialnonparbootbc", B = 3, seed = 123)
            
            # Estimation with mseFH of sae (benchmark)
            fh_spatial_reml_npb_sae <- read.csv("FH/fh_spatial_reml_npb_sae.csv", 
                                                sep = ",", 
                                                stringsAsFactors = TRUE) 
            
            # Comparison
            # Variance
            expect_equal(fh_spatial_reml_npb$model$variance[1,2], 
                         fh_spatial_reml_npb_sae$variance[1])
            # Correlation parameter
            expect_equal(fh_spatial_reml_npb$model$variance[1,1], 
                         fh_spatial_reml_npb_sae$correlation[1])
            # EBLUP
            expect_equal(fh_spatial_reml_npb$ind$FH, 
                         as.vector(fh_spatial_reml_npb_sae$EBLUP))
            # MSE
            expect_equal(fh_spatial_reml_npb$MSE$FH, 
                         fh_spatial_reml_npb_sae$MSE, tolerance = 0.000001)
            # MSE bias corrected
            expect_equal(fh_spatial_reml_npb_bc$MSE$FH, 
                         fh_spatial_reml_npb_sae$MSE.BC, tolerance = 0.000001)
          })

# Parametric bootstrap MSE 

test_that("Does the fh function in emdi return the same variance, correlation 
          parameter, EBLUP and MSE estimates as the function pbmseSFH of package 
          sae?",{
            
            ############################ REML variance estimation ########################
            # Estimation with fh of emdi
            grapes$Domain <- c(1:274)
            fh_spatial_reml_pb <- fh(fixed = grapehect ~ area + workdays - 1, 
                                     vardir = "var", tol = 0.0001, maxit = 100, 
                                     combined_data = grapes, domains = "Domain", 
                                     method = "reml", correlation = "spatial", 
                                     corMatrix = as.matrix(grapesprox), MSE = TRUE, 
                                     mse_type = "spatialparboot", B = 3, seed = 123)
            
            fh_spatial_reml_pb_bc <- fh(fixed = grapehect ~ area + workdays - 1, 
                                        vardir = "var", tol = 0.0001, maxit = 100, 
                                        combined_data = grapes, domains = "Domain", 
                                        method = "reml", correlation = "spatial", 
                                        corMatrix = as.matrix(grapesprox), MSE = TRUE, 
                                        mse_type = "spatialparbootbc", B = 3, seed = 123)
            
            # Estimation with mseFH of sae (benchmark)
            fh_spatial_reml_pb_sae <- read.csv("FH/fh_spatial_reml_pb_sae.csv", sep = ",", 
                                               stringsAsFactors = TRUE) 
            
            # Comparison
            # Variance
            expect_equal(fh_spatial_reml_pb$model$variance[1,2], 
                         fh_spatial_reml_pb_sae$variance[1])
            # Correlation parameter
            expect_equal(fh_spatial_reml_pb$model$variance[1,1], 
                         fh_spatial_reml_pb_sae$correlation[1])
            # EBLUP
            expect_equal(fh_spatial_reml_pb$ind$FH, 
                         as.vector(fh_spatial_reml_pb_sae$EBLUP))
            # MSE
            expect_equal(fh_spatial_reml_pb$MSE$FH, 
                         fh_spatial_reml_pb_sae$MSE)
            # MSE bias corrected
            expect_equal(fh_spatial_reml_pb_bc$MSE$FH, 
                         fh_spatial_reml_pb_sae$MSE.BC)
            
            ############################ ML variance estimation ##########################
            # Estimation with fh of emdi
            grapes$Domain <- c(1:274)
            fh_spatial_ml_pb <- fh(fixed = grapehect ~ area + workdays - 1, 
                                   vardir = "var", tol = 0.0001, maxit = 100, 
                                   combined_data = grapes, domains = "Domain", 
                                   method = "ml", correlation = "spatial", 
                                   corMatrix = as.matrix(grapesprox), MSE = TRUE, 
                                   mse_type = "spatialparboot", B = 3, seed = 123)
            
            fh_spatial_ml_pb_bc <- fh(fixed = grapehect ~ area + workdays - 1, 
                                      vardir = "var", tol = 0.0001, maxit = 100, 
                                      combined_data = grapes, domains = "Domain", 
                                      method = "ml", correlation = "spatial", 
                                      corMatrix = as.matrix(grapesprox), MSE = TRUE, 
                                      mse_type = "spatialparbootbc", B = 3, seed = 123)
            
            # Estimation with mseFH of sae (benchmark)
            fh_spatial_ml_pb_sae <- read.csv("FH/fh_spatial_ml_pb_sae.csv", sep = ",", 
                                             stringsAsFactors = TRUE) 
            
            # Comparison
            # Variance
            expect_equal(fh_spatial_ml_pb$model$variance[1,2], 
                         fh_spatial_ml_pb_sae$variance[1])
            # Correlation parameter
            expect_equal(fh_spatial_ml_pb$model$variance[1,1], 
                         fh_spatial_ml_pb_sae$correlation[1])
            # EBLUP
            expect_equal(fh_spatial_ml_pb$ind$FH, 
                         fh_spatial_ml_pb_sae$EBLUP)
            # MSE
            expect_equal(fh_spatial_ml_pb$MSE$FH, 
                         fh_spatial_ml_pb_sae$MSE)
            # MSE bias corrected
            expect_equal(fh_spatial_ml_pb_bc$MSE$FH, 
                         fh_spatial_ml_pb_sae$MSE.BC)
          })


# Test if the same variance, EBLUP and MSE results are obtained as with 
# the saeRobust package

########################## Robust area-level model #############################


# Pseudo MSE

test_that("Does the fh function in emdi return the same variance, EBLUP and MSE 
          estimates as the functions rfh and mse (pseudo) of package saeRobust?",{
            
            select <- as.logical(Sys.getenv("_R_TEST_ROBUST_TRUEFALSE_"))
            if (is.na(select)) {
              select <- FALSE
            }
            
            ############################ REBLUP model fitting ############################
            if (select) {
            # Estimation with fh of emdi
              grapes$Domain <- c(1:274)
              fh_robust <- fh(fixed = grapehect ~ area + workdays - 1, vardir = "var",
                              combined_data = grapes, domains = "Domain",
                              method = "reblup", tol = 1e-06, maxit = 100, k = 1.345,
                              MSE = TRUE, mse_type = "pseudo")
              
              # Estimation with mseFH of sae (benchmark)
              fh_robust_saeRobust <- read.csv("FH/fh_robust_saeRobust.csv", sep = ",", 
                                              stringsAsFactors = TRUE) 
              
              # Comparison
              # Variance
              expect_equal(unname(fh_robust$model$variance), 
                           fh_robust_saeRobust$variance[1])
              # EBLUP
              expect_equal(fh_robust$ind$FH, fh_robust_saeRobust$EBLUP)
              # MSE
              expect_equal(fh_robust$MSE$FH, fh_robust_saeRobust$MSE)
              
              ############################ REBLUPBC model fitting ##########################
              
              # Estimation with fh of emdi
              grapes$Domain <- c(1:274)
              fh_robustbc <- fh(fixed = grapehect ~ area + workdays - 1, vardir = "var",
                                combined_data = grapes, domains = "Domain",
                                method = "reblupbc", tol = 1e-06, maxit = 100, k = 1.345, 
                                c = 2, MSE = TRUE, mse_type = "pseudo")
              
              # Estimation with fitRFH of saeRobust (benchmark)
              fh_robustbc_saeRobust <- read.csv("FH/fh_robustbc_saeRobust.csv", sep = ",", 
                                                stringsAsFactors = TRUE) 
              
              # Comparison
              # Variance
              expect_equal(unname(fh_robustbc$model$variance), 
                           fh_robust_saeRobust$variance[1])
              # EBLUP
              expect_equal(fh_robustbc$ind$FH, fh_robustbc_saeRobust$EBLUP)
              # MSE
              expect_equal(fh_robustbc$MSE$FH, fh_robustbc_saeRobust$MSE)
            } else {
              expect_equal(TRUE, TRUE)
            }
          })

# Bootstrap MSE

test_that("Does the fh function in emdi return the same variance, EBLUP and MSE 
          estimates as the functions rfh and mse (boot) of package saeRobust?",{
            
            select <- as.logical(Sys.getenv("_R_TEST_ROBUST_TRUEFALSE_"))
            if (is.na(select)) {
              select <- FALSE
            }
            
            ############################ REBLUP model fitting ############################
            if (select) {
              # Estimation with fh of emdi
              grapes$Domain <- c(1:274)
              fh_robust_boot <- fh(fixed = grapehect ~ area + workdays - 1, vardir = "var",
                                   combined_data = grapes, domains = "Domain",
                                   method = "reblup", tol = 1e-06, maxit = 100, k = 1.345,
                                   MSE = TRUE, mse_type = "boot", B = 3, seed = 123)
              
              # Estimation with fitRFH of saeRobust (benchmark)
              fh_robust_boot_saeRobust <- read.csv("FH/fh_robust_boot_saeRobust.csv", 
                                                   sep = ",", 
                                                   stringsAsFactors = TRUE) 
              
              # Comparison
              # Variance
              expect_equal(unname(fh_robust_boot$model$variance), 
                           fh_robust_boot_saeRobust$variance[1])
              # EBLUP
              expect_equal(fh_robust_boot$ind$FH, fh_robust_boot_saeRobust$EBLUP)
              # MSE
              expect_equal(fh_robust_boot$MSE$FH, fh_robust_boot_saeRobust$MSE)
              
              ############################ REBLUPBC model fitting ##########################
              
              # Estimation with fh of emdi
              grapes$Domain <- c(1:274)
              fh_robustbc_boot <- fh(fixed = grapehect ~ area + workdays - 1, vardir = "var", 
                                     combined_data = grapes, domains = "Domain", 
                                     method = "reblupbc", tol = 1e-06, maxit = 100, k = 1.345, 
                                     c = 2, MSE = TRUE, mse_type = "boot", B = 3, seed = 123)
              
              # Estimation with fitRFH of saeRobust (benchmark)
              fh_robustbc_boot_saeRobust <- read.csv("FH/fh_robustbc_boot_saeRobust.csv", 
                                                     sep = ",", 
                                                     stringsAsFactors = TRUE) 
              
              # Comparison
              # Variance
              expect_equal(unname(fh_robustbc_boot$model$variance), 
                           fh_robust_boot_saeRobust$variance[1])
              # EBLUP
              expect_equal(fh_robustbc_boot$ind$FH, fh_robustbc_boot_saeRobust$EBLUP)
              # MSE
              expect_equal(fh_robustbc_boot$MSE$FH, fh_robustbc_boot_saeRobust$MSE)
            } else {
              expect_equal(TRUE, TRUE)
            }
          })

######################## Robust spatial area-level model #######################

# Pseudo MSE

test_that("Does the fh function in emdi return the same variance, correlation 
           parameter, EBLUP and MSE estimates as the functions rfh and mse 
          (pseudo) of package saeRobust?",{
            
            select <- as.logical(Sys.getenv("_R_TEST_ROBUST_TRUEFALSE_"))
            if (is.na(select)) {
              select <- FALSE
            }
            ############################ REBLUP model fitting ############################
            if (select) {
                
              # Estimation with fh of emdi
              grapes$Domain <- c(1:274)
              fh_robust_spatial <- fh(fixed = grapehect ~ area + workdays - 1, vardir = "var",
                                      combined_data = grapes, domains = "Domain",
                                      method = "reblup", tol = 1e-06, maxit = 100, k = 1.345,
                                      correlation = "spatial", corMatrix = grapesprox,
                                      MSE = TRUE, mse_type = "pseudo")
              
              # Estimation with fitRFH of saeRobust (benchmark)
              fh_robust_spatial_saeRobust <- read.csv("FH/fh_robust_spatial_saeRobust.csv", 
                                                      sep = ",", 
                                                      stringsAsFactors = TRUE)
              # Comparison
              # Variance
              expect_equal(unname(fh_robust_spatial$model$variance[2]), 
                           fh_robust_spatial_saeRobust$variance[1])
              # Correlation parameter
              expect_equal(unname(fh_robust_spatial$model$variance[1]), 
                           fh_robust_spatial_saeRobust$correlation[1])
              # EBLUP
              expect_equal(fh_robust_spatial$ind$FH, fh_robust_spatial_saeRobust$EBLUP)
              # MSE
              expect_equal(fh_robust_spatial$MSE$FH, fh_robust_spatial_saeRobust$MSE)
              
              ############################ REBLUPBC model fitting ##########################
              
              # Estimation with fh of emdi
              grapes$Domain <- c(1:274)
              fh_robust_spatial_bc <- fh(fixed = grapehect ~ area + workdays - 1, vardir = "var",
                                         combined_data = grapes, domains = "Domain",
                                         method = "reblupbc", tol = 1e-06, maxit = 100, 
                                         k = 1.345, c = 2, 
                                         correlation = "spatial", corMatrix = grapesprox, 
                                         MSE = TRUE, mse_type = "pseudo")
              
              # Estimation with fitRFH of saeRobust (benchmark)
              fh_robustbc_spatial_saeRobust <- read.csv("FH/fh_robustbc_spatial_saeRobust.csv", 
                                                        sep = ",", 
                                                        stringsAsFactors = TRUE) 
              
              # Comparison
              # Variance
              expect_equal(unname(fh_robust_spatial_bc$model$variance[2]), 
                           fh_robustbc_spatial_saeRobust$variance[1])
              # Correlation parameter
              expect_equal(unname(fh_robust_spatial_bc$model$variance[1]), 
                           fh_robustbc_spatial_saeRobust$correlation[1])
              # EBLUP
              expect_equal(fh_robust_spatial_bc$ind$FH, 
                           fh_robustbc_spatial_saeRobust$EBLUP)
              # MSE
              expect_equal(fh_robust_spatial_bc$MSE$FH, 
                           fh_robustbc_spatial_saeRobust$MSE)
            } else {
              expect_equal(TRUE, TRUE)
            }
            
          })

# Bootstrap MSE

test_that("Does the fh function in emdi return the same variance, correlation 
           parameter, EBLUP and MSE estimates as the functions rfh and mse 
          (boot) of package saeRobust?",{
            
            select <- as.logical(Sys.getenv("_R_TEST_ROBUST_TRUEFALSE_"))
            if (is.na(select)) {
              select <- FALSE
            }
            ############################ REBLUP model fitting ############################
            
            if (select) {
             
              # Estimation with fh of emdi
              grapes$Domain <- c(1:274)
              fh_robust_spatial_boot <- fh(fixed = grapehect ~ area + workdays - 1, 
                                           vardir = "var", combined_data = grapes, 
                                           domains = "Domain", method = "reblup", 
                                           tol = 1e-06, maxit = 100, k = 1.345,
                                           correlation = "spatial", corMatrix = grapesprox,
                                           MSE = TRUE, mse_type = "boot", B = 3, seed = 123)
              
              # Estimation with fitRFH of saeRobust (benchmark)
              fh_robust_spatial_boot_saeRobust <- read.csv("FH/fh_robust_spatial_boot_saeRobust.csv", 
                                                           sep = ",", 
                                                           stringsAsFactors = TRUE) 
              
              # Comparison
              # Variance
              expect_equal(unname(fh_robust_spatial_boot$model$variance[2]), 
                           fh_robust_spatial_boot_saeRobust$variance[1])
              # Correlation parameter
              expect_equal(unname(fh_robust_spatial_boot$model$variance[1]), 
                           fh_robust_spatial_boot_saeRobust$correlation[1])
              # EBLUP
              expect_equal(fh_robust_spatial_boot$ind$FH, 
                           fh_robust_spatial_boot_saeRobust$EBLUP)
              # MSE
              expect_equal(fh_robust_spatial_boot$MSE$FH, 
                           fh_robust_spatial_boot_saeRobust$MSE)
              
              ############################ REBLUPBC model fitting ##########################
              
              # Estimation with fh of emdi
              grapes$Domain <- c(1:274)
              fh_robust_spatial_bc_boot <- fh(fixed = grapehect ~ area + workdays - 1, 
                                              vardir = "var", combined_data = grapes, 
                                              domains = "Domain", method = "reblupbc", 
                                              tol = 1e-06, maxit = 100, k = 1.345, c = 2, 
                                              correlation = "spatial", corMatrix = grapesprox,
                                              MSE = TRUE, mse_type = "boot", 
                                              B = 3, seed = 123)
              
              # Estimation with fitRFH of saeRobust (benchmark)
              fh_robust_spatial_bc_boot_saeRobust <- read.csv("FH/fh_robust_spatial_bc_boot_saeRobust.csv", 
                                                              sep = ",", 
                                                              stringsAsFactors = TRUE) 
              
              # Comparison
              # Variance
              expect_equal(unname(fh_robust_spatial_bc_boot$model$variance[2]),
                           fh_robust_spatial_bc_boot_saeRobust$variance[1])
              # Correlation parameter
              expect_equal(unname(fh_robust_spatial_bc_boot$model$variance[1]),
                           fh_robust_spatial_bc_boot_saeRobust$correlation[1])
              # EBLUP
              expect_equal(fh_robust_spatial_bc_boot$ind$FH, 
                           fh_robust_spatial_bc_boot_saeRobust$EBLUP)
              # MSE
              expect_equal(fh_robust_spatial_bc_boot$MSE$FH, 
                           fh_robust_spatial_bc_boot_saeRobust$MSE)
            } else {
              expect_equal(TRUE, TRUE)
            }
          })


