# Test if variance estimates, EBLUPs and MSEs remains the same when applying 
# transformations

# Loading data - population and sample data
load("FH/eusilcA_popAgg.RData")
load("FH/eusilcA_smpAgg.RData")

# Combine sample and population data 
combined_data <- combine_data(pop_data = eusilcA_popAgg, pop_domains = "Domain",
                              smp_data = eusilcA_smpAgg, smp_domains = "Domain")

############################## Log transformation ##############################

test_that("Does the fh function with a log transformation return the same 
          results as the ones obtained in simulation studies?", {
            
  ########################## crude backtransformation ##########################
            # Current version 
            fh_log_crude <- fh(fixed = Mean ~ cash + self_empl, 
                               vardir = "Var_Mean", 
                               combined_data = combined_data, domains = "Domain",
                               method = "ml", interval = c(0, 10000000), 
                               transformation = "log", 
                               backtransformation = "bc_crude", MSE = TRUE)
          
            # Status quo (benchmark)
            transf_log_crude <- read.csv("FH/transf_log_crude.csv", sep = ",", 
                                         stringsAsFactors = TRUE)  
            
            # Compare results from current version and benchmark
            # EBLUP
            expect_equal(fh_log_crude$ind[, c("Domain","FH")], 
                         transf_log_crude[, c("Domain","FH")])
            # MSE
            expect_equal(fh_log_crude$MSE$FH, transf_log_crude$MSE)
            # Variance
            expect_equal(fh_log_crude$model$variance, 
                         transf_log_crude$variance[1])
            
  ############################ sm backtransformation ###########################
            # Current version 
            fh_log_sm <- fh(fixed = Mean ~ cash + self_empl, 
                            vardir = "Var_Mean", 
                            combined_data = combined_data, domains = "Domain",
                            method = "ml", interval = c(0, 10000000), 
                            transformation = "log", backtransformation = "bc_sm", 
                            MSE = TRUE)
          
            # Status quo (benchmark)
            transf_log_sm <- read.csv("FH/transf_log_sm.csv", sep = ",", 
                                      stringsAsFactors = TRUE)  
            
            # Compare results from current version and benchmark
            # EBLUP
            expect_equal(fh_log_sm$ind[, c("Domain","FH")], 
                         transf_log_sm[, c("Domain","FH")])
            # MSE
            expect_equal(fh_log_sm$MSE$FH, transf_log_sm$MSE)
            # Variance
            expect_equal(fh_log_sm$model$variance, 
                         transf_log_sm$variance[1])          
            })

############################## arcsin transformation ###########################

test_that("Does the fh function with a arcsin transformation return the same 
          results as the ones obtained in simulation studies?", {
            
  ########################## naive backtransformation ##########################
            # Current version 
            fh_arcsin_naive_jack <- fh(fixed = MTMED ~ cash + age_ben + rent + 
                                         house_allow, 
                                    vardir = "Var_MTMED", 
                                    combined_data = combined_data, 
                                    domains = "Domain",
                                    method = "ml", interval = c(0, 10000000), 
                                    transformation = "arcsin", 
                                    backtransformation = "naive", 
                                    eff_smpsize = "n", MSE = TRUE,
                                    mse_type = "jackknife")
            fh_arcsin_naive_wjack <- fh(fixed = MTMED ~ cash + age_ben + rent + 
                                          house_allow, 
                                       vardir = "Var_MTMED", 
                                       combined_data = combined_data, 
                                       domains = "Domain",
                                       method = "ml", interval = c(0, 10000000), 
                                       transformation = "arcsin", 
                                       backtransformation = "naive", 
                                       eff_smpsize = "n", MSE = TRUE,
                                       mse_type = "weighted_jackknife")
            
            # Status quo (benchmark)
            transf_arcsin_naive <- read.csv("FH/transf_arcsin_naive.csv", 
                                            sep = ",", 
                                            stringsAsFactors = TRUE)  
            
            # Compare results from current version and benchmark
            # EBLUP
            expect_equal(fh_arcsin_naive_jack$ind[, c("Domain","FH")], 
                         transf_arcsin_naive[, c("Domain","FH")])
            # MSE jackknife
            expect_equal(fh_arcsin_naive_jack$MSE$FH, 
                         transf_arcsin_naive$MSE_jack)
            # MSE weighted jackknife
            expect_equal(fh_arcsin_naive_wjack$MSE$FH, 
                         transf_arcsin_naive$MSE_wjack)
            # Variance
            expect_equal(fh_arcsin_naive_jack$model$variance, 
                         transf_arcsin_naive$variance[1])
            
  ############################ sm backtransformation ###########################
            # Current version 
           # fh_arcsin_sm_jack <- fh(fixed = MTMED ~ cash + age_ben + rent + house_allow, 
            #                        vardir = "Var_MTMED", 
            #                        combined_data = combined_data, 
            #                        domains = "Domain",
             #                       method = "ml", interval = c(0, 10000000), 
              #                      transformation = "arcsin", 
              #                      backtransformation = "bc", 
              #                      eff_smpsize = "n", MSE = TRUE,
               #                     mse_type = "jackknife")
           # fh_arcsin_sm_wjack <- fh(fixed = MTMED ~ cash + age_ben + rent + house_allow, 
            #                         vardir = "Var_MTMED", 
             #                        combined_data = combined_data, 
             #                        domains = "Domain",
             #                        method = "ml", interval = c(0, 10000000), 
             #                        transformation = "arcsin", 
             #                        backtransformation = "bc", 
             #                        eff_smpsize = "n", MSE = TRUE,
             #                        mse_type = "weighted_jackknife")
            fh_arcsin_sm_boot <- fh(fixed = MTMED ~ cash + age_ben + rent + house_allow, 
                                    vardir = "Var_MTMED", 
                                    combined_data = combined_data, 
                                    domains = "Domain",
                                    method = "ml", interval = c(0, 10000000), 
                                    transformation = "arcsin", 
                                    backtransformation = "bc", 
                                    eff_smpsize = "n", MSE = TRUE,
                                    mse_type = "boot", B = 3, seed = 123)
           
            # Status quo (benchmark)
            transf_arcsin_sm <- read.csv("FH/transf_arcsin_sm.csv", sep = ",", 
                                         stringsAsFactors = TRUE)  
           
            # Compare results from current version and benchmark
            # EBLUP
            expect_equal(fh_arcsin_sm_boot$ind[, c("Domain","FH")], 
                         transf_arcsin_sm[, c("Domain","FH")])
            # MSE jackknife
           # expect_equal(fh_arcsin_sm_boot$MSE$FH, 
             #            transf_arcsin_sm$MSE_jack)
            # MSE weighted jackknife
          #  expect_equal(fh_arcsin_sm_boot$MSE$FH, 
            #             transf_arcsin_sm$MSE_wjack)
            # MSE bootstrap
            expect_equal(fh_arcsin_sm_boot$MSE$FH, 
                         transf_arcsin_sm$MSE_boot)
            # Variance
            expect_equal(fh_arcsin_sm_boot$model$variance, 
                         transf_arcsin_sm$variance[1])          
            })
