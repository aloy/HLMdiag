context("tests for case_delete") 

#sleepstudy models 
data(sleepstudy, package = 'lme4')
sleep.lmer <- lme4::lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy)
sleep.lme <- nlme::lme(Reaction ~ Days, random =  ~ Days|Subject, data = sleepstudy)

sleep.lmer.case <- case_delete(sleep.lmer)
sleep.lme.case <- case_delete(sleep.lme)

sleep.lmer.caseG <- case_delete(sleep.lmer, level = "Subject")
sleep.lme.caseG <- case_delete(sleep.lme, level = "Subject")

sleep.lmer.caseD <- case_delete(sleep.lmer, delete = c(1, 20, 100))
sleep.lme.caseD <- case_delete(sleep.lmer, delete = c(1, 20, 100))

sleep.lmer.caseGD <- case_delete(sleep.lmer, level = "Subject", delete = "308")
sleep.lme.caseGD <- case_delete(sleep.lme, level = "Subject", delete = "308") 

#chemistry scores models - 3 level 

data(Chem97, package = "mlmRev")  
Chem97 <- Chem97[1:257,]
chem.lmer <- lme4::lmer(score ~ gcsecnt + (1|lea/school), data = Chem97)
chem.lme <- nlme::lme(score ~ gcsecnt, random = ~1|lea/school, data = Chem97)

chem.lmer.case <- suppressMessages(case_delete(chem.lmer)) #this throws isSingular warnings
chem.lme.case <- case_delete(chem.lme)

chem.lmer.caseG <- suppressMessages(case_delete(chem.lmer, level = "lea")) #also throws isSingular warnings 
chem.lme.caseG <- case_delete(chem.lme, level = "lea")

chem.lmer.caseG2 <- suppressMessages(case_delete(chem.lmer, level = "school:lea")) #also throws isSingular warnings
chem.lme.caseG2 <- case_delete(chem.lme, level = "school") #school instead of school:lea 


#passed
test_that("Correct influence diagnostics are returned", {
  expect_equal(length(sleep.lmer.case), 9)
  expect_equal(length(sleep.lme.case), 9)
  expect_equal(length(chem.lmer.case), 9)
  expect_equal(length(chem.lme.case), 9)

})


#passes
test_that("Original fixed effects matches output from fixef",{
  #sleepstudy
  expect_equal(sleep.lmer.case$fixef.original, fixef(sleep.lmer))
  expect_equal(sleep.lme.case$fixef.original, fixef(sleep.lme))
  expect_equal(sleep.lmer.caseG$fixef.original, fixef(sleep.lmer))
  expect_equal(sleep.lme.caseG$fixef.original, fixef(sleep.lme))
  
  #chemistry - still need chemistry group model 
  expect_equal(chem.lmer.case$fixef.original, fixef(chem.lmer))
  expect_equal(chem.lme.case$fixef.original, fixef(chem.lme))
  expect_equal(chem.lmer.caseG$fixef.original, fixef(chem.lmer))
  expect_equal(chem.lme.caseG$fixef.original, fixef(chem.lme))
  expect_equal(chem.lmer.caseG2$fixef.original, fixef(chem.lmer))
  expect_equal(chem.lme.caseG2$fixef.original, fixef(chem.lme))
})



#passes
test_that("Original predicted random effects match output from ranef", {
  
  #sleepstudy
  expect_equal(sleep.lmer.case$ranef.original[[1]], ranef(sleep.lmer)[[1]][,1])
  expect_equal(sleep.lmer.case$ranef.original[[2]], ranef(sleep.lmer)[[1]][,2])
  
  expect_equal(sleep.lme.case$ranef.original[[1]], ranef(sleep.lme)[,1])
  expect_equal(sleep.lme.case$ranef.original[[2]], ranef(sleep.lme)[,2])
  
  expect_equal(sleep.lmer.caseG$ranef.original[[1]], ranef(sleep.lmer)[[1]][,1])
  expect_equal(sleep.lmer.caseG$ranef.original[[2]], ranef(sleep.lmer)[[1]][,2])
  
  expect_equal(sleep.lme.caseG$ranef.original[[1]], ranef(sleep.lme)[,1])
  expect_equal(sleep.lme.caseG$ranef.original[[2]], ranef(sleep.lme)[,2])
  
  #chemistry 
  expect_equal(chem.lmer.case$ranef.original[[1]][,1], ranef(chem.lmer)[[1]][,1]) 
  expect_equal(chem.lmer.case$ranef.original[[2]][,1], ranef(chem.lmer)[[2]][,1]) 
  
  expect_equal(chem.lme.case$ranef.original[[1]][,1], ranef(chem.lme)[[1]][,1]) 
  expect_equal(chem.lme.case$ranef.original[[2]][,1], ranef(chem.lme)[[2]][,1]) 
  
  expect_equal(chem.lmer.caseG$ranef.original[[1]][,1], ranef(chem.lmer)[[1]][,1])
  expect_equal(chem.lmer.caseG$ranef.original[[2]][,1], ranef(chem.lmer)[[2]][,1])
  
  expect_equal(chem.lme.caseG$ranef.original[[1]][,1], ranef(chem.lme)[[1]][,1])
  expect_equal(chem.lme.caseG$ranef.original[[2]][,1], ranef(chem.lme)[[2]][,1])
  
  
  expect_equal(chem.lmer.caseG2$ranef.original[[1]][,1], ranef(chem.lmer)[[1]][,1])
  expect_equal(chem.lmer.caseG2$ranef.original[[2]][,1], ranef(chem.lmer)[[2]][,1])
  
  expect_equal(chem.lme.caseG2$ranef.original[[1]][,1], ranef(chem.lme)[[1]][,1])
  expect_equal(chem.lme.caseG2$ranef.original[[2]][,1], ranef(chem.lme)[[2]][,1])
  
  
})

#passes 
test_that("Original variance-covariance matrix matches output from vcov", {
  #sleepstudy
  expect_equal(sleep.lmer.case$vcov.original, as.matrix(vcov(sleep.lmer)))
  expect_equal(sleep.lme.case$vcov.original, as.matrix(vcov(sleep.lme)))
  
  expect_equal(sleep.lmer.caseG$vcov.original, as.matrix(vcov(sleep.lmer)))
  expect_equal(sleep.lme.caseG$vcov.original, as.matrix(vcov(sleep.lme)))
  
  #chemistry 
  expect_equal(chem.lmer.case$vcov.original, as.matrix(vcov(chem.lmer)))
  expect_equal(chem.lme.case$vcov.original, as.matrix(vcov(chem.lme)))
  
  expect_equal(chem.lmer.caseG$vcov.original, as.matrix(vcov(chem.lmer)))
  expect_equal(chem.lme.caseG$vcov.original, as.matrix(vcov(chem.lme)))
  
  expect_equal(chem.lmer.caseG2$vcov.original, as.matrix(vcov(chem.lmer)))
  expect_equal(chem.lme.caseG2$vcov.original, as.matrix(vcov(chem.lme)))
  
})

#passed
test_that("Variance components matches output from varcomp", {
  #sleepstudy
  expect_equal(sleep.lmer.case$varcomp.original, varcomp.mer(sleep.lmer))
  expect_equal(sleep.lme.case$varcomp.original, varcomp.lme(sleep.lme))
  
  expect_equal(sleep.lmer.caseG$varcomp.original, varcomp.mer(sleep.lmer))
  expect_equal(sleep.lme.caseG$varcomp.original, varcomp.lme(sleep.lme))
  
  #chemistry 
  expect_equal(chem.lmer.case$varcomp.original, varcomp.mer(chem.lmer))
  expect_equal(chem.lme.case$varcomp.original, varcomp.lme(chem.lme))
  
  expect_equal(chem.lmer.caseG$varcomp.original, varcomp.mer(chem.lmer))
  expect_equal(chem.lme.caseG$varcomp.original, varcomp.lme(chem.lme))
  
  expect_equal(chem.lmer.caseG2$varcomp.original, varcomp.mer(chem.lmer))
  expect_equal(chem.lme.caseG2$varcomp.original, varcomp.lme(chem.lme))
  
  
})


#passed
test_that("Dimensions of fixed effects after deletion are correct for single case deletion", {
  #number of rows is number of observations, number of columns is number of fixed effects plus one
  
  #sleepstudy
  expect_equal(nrow(sleep.lmer.case$fixef.delete), nrow(sleep.lmer@frame))
  expect_equal(ncol(sleep.lmer.case$fixef.delete), 1 + length(fixef(sleep.lmer)))
  expect_equal(length(sleep.lmer.caseD$fixef.delete), length(fixef(sleep.lmer)))
  
  expect_equal(nrow(sleep.lme.case$fixef.delete), nrow(sleep.lme$groups))
  expect_equal(ncol(sleep.lme.case$fixef.delete), 1 + length(fixef(sleep.lme)))
  expect_equal(length(sleep.lme.caseD$fixef.delete), length(fixef(sleep.lme)))
  
  #chemistry 
  expect_equal(nrow(chem.lmer.case$fixef.delete), nrow(chem.lmer@frame))
  expect_equal(ncol(chem.lmer.case$fixef.delete), 1 + length(fixef(chem.lmer)))
  expect_equal(nrow(chem.lme.case$fixef.delete), nrow(chem.lme$groups))
  expect_equal(ncol(chem.lme.case$fixef.delete), 1 + length(fixef(chem.lme)))
  
  #add chemistry delete stuff here too? 
  
})


#passed
test_that("Dimensions of fixed effects after deletion are correct for group deletion", {
  #number of rows is number of groups, number of columns is number of fixed effects 
  #sleepstudy
  expect_equal(nrow(sleep.lmer.caseG$fixef.delete), length(unique(sleep.lmer@flist[["Subject"]])))
  expect_equal(ncol(sleep.lmer.caseG$fixef.delete), length(fixef(sleep.lmer)))
  expect_equal(length(sleep.lmer.caseGD$fixef.delete), length(fixef(sleep.lmer)))
  
  expect_equal(nrow(sleep.lme.caseG$fixef.delete), length(unique(sleep.lme$groups$Subject)))
  expect_equal(ncol(sleep.lme.caseG$fixef.delete), length(fixef(sleep.lme)))
  expect_equal(length(sleep.lme.caseGD$fixef.delete), length(fixef(sleep.lme)))
  
  #chemistry
  expect_equal(nrow(chem.lmer.caseG$fixef.delete), length(unique(chem.lmer@flist[["lea"]])))
  expect_equal(ncol(chem.lmer.caseG$fixef.delete), length(fixef(chem.lmer)))
  
  expect_equal(nrow(chem.lme.caseG$fixef.delete), length(unique(chem.lme$groups[["lea"]])))
  expect_equal(ncol(chem.lme.caseG$fixef.delete), length(fixef(chem.lme)))
  
  expect_equal(nrow(chem.lmer.caseG2$fixef.delete), length(unique(chem.lmer@flist[["school:lea"]])))
  expect_equal(ncol(chem.lmer.caseG2$fixef.delete), length(fixef(chem.lmer)))
  
  expect_equal(nrow(chem.lme.caseG2$fixef.delete), length(unique(chem.lme$groups[["school"]])))
  expect_equal(ncol(chem.lme.caseG2$fixef.delete), length(fixef(chem.lme)))
  
})


#passed
test_that("Dimensions of random effects after deletion are correct for single case deletion", {
  #number of rows is number of groups times number of observations, number of columns is two plus number of random effects
  #if delete was set, number of columns is just number of random effects, number of rows is the number of groups 
  
  #sleepstudy 
  expect_equal(nrow(sleep.lmer.case$ranef.delete), nrow(sleep.lmer@frame) * length(unique(sleep.lmer@flist[["Subject"]])))
  expect_equal(ncol(sleep.lmer.case$ranef.delete), 2 + ncol(ranef(sleep.lmer)$Subject))
  expect_equal(nrow(sleep.lmer.caseD$ranef.delete), length(unique(sleep.lmer@flist[["Subject"]])))
  expect_equal(ncol(sleep.lmer.caseD$ranef.delete), ncol(ranef(sleep.lmer)$Subject))
               
  
  expect_equal(nrow(sleep.lme.case$ranef.delete), nrow(sleep.lme$groups) * length(unique(sleep.lme$groups$Subject)))
  expect_equal(ncol(sleep.lme.case$ranef.delete), 2 + ncol(ranef(sleep.lme)))
  expect_equal(nrow(sleep.lme.caseD$ranef.delete), length(unique(sleep.lme$groups$Subject)))
  expect_equal(ncol(sleep.lme.caseD$ranef.delete), ncol(ranef(sleep.lme)))
  
  #chemistry 
  
  expect_equal(nrow(chem.lmer.case$ranef.delete[[1]]), -3 + (nrow(chem.lmer@frame) * length(unique(chem.lmer@flist[["school:lea"]]))))
  expect_equal(ncol(chem.lmer.case$ranef.delete[[1]]), 2 + ncol(ranef(chem.lmer)$'school:lea'))
  
  expect_equal(nrow(chem.lmer.case$ranef.delete[[2]]), nrow(chem.lmer@frame) * length(unique(chem.lmer@flist[["lea"]])))
  expect_equal(ncol(chem.lmer.case$ranef.delete[[2]]), 2 + ncol(ranef(chem.lmer)$lea))
  
  #ranef flips the order for lme 
  expect_equal(nrow(chem.lme.case$ranef.delete[[2]]), -3 + nrow(chem.lme$groups) * length(unique(chem.lme$groups[["school"]])))
  expect_equal(ncol(chem.lme.case$ranef.delete[[2]]), 2 + ncol(ranef(chem.lme)$'school'))
  
  expect_equal(nrow(chem.lme.case$ranef.delete[[1]]), nrow(chem.lme$groups) * length(unique(chem.lme$groups$lea)))
  expect_equal(ncol(chem.lme.case$ranef.delete[[1]]), 2 + ncol(ranef(chem.lme)$lea))
  
  
})

#passed
test_that("Dimensions of random effects after deletion are correct for group deletion ", {
  #number of rows is number of groups times (number of groups minus 1), columns is two plus number of random effects
  
  #sleepstudy
  n <- length(unique(sleep.lmer@flist[["Subject"]]))
  expect_equal(nrow(sleep.lmer.caseG$ranef.delete), n * (n-1))
  expect_equal(ncol(sleep.lmer.caseG$ranef.delete), 2 + ncol(ranef(sleep.lmer)$Subject))
  expect_equal(nrow(sleep.lmer.caseGD$ranef.delete), length(unique(sleep.lmer@flist[["Subject"]])) - 1)
  expect_equal(ncol(sleep.lmer.caseGD$ranef.delete), ncol(ranef(sleep.lmer)$Subject))
  
  n <- length(unique(sleep.lme$groups$Subject))
  expect_equal(nrow(sleep.lme.caseG$ranef.delete), n * (n-1))
  expect_equal(ncol(sleep.lme.caseG$ranef.delete), 2 + ncol(ranef(sleep.lme)))
  expect_equal(nrow(sleep.lme.caseGD$ranef.delete), length(unique(sleep.lme$groups$Subject)) - 1)
  expect_equal(ncol(sleep.lme.caseGD$ranef.delete), ncol(ranef(sleep.lme)))
  
  #chemistry
  n <- length(unique(chem.lmer@flist[["lea"]]))
  expect_equal(nrow(chem.lmer.caseG$ranef.delete[[2]]), n * (n-1))  
  expect_equal(ncol(chem.lmer.caseG$ranef.delete[[2]]), 2 + ncol(ranef(chem.lmer)$lea))
  
  n <- length(unique(chem.lme$groups[["lea"]]))
  expect_equal(nrow(chem.lme.caseG$ranef.delete[[1]]), n * (n-1))
  expect_equal(ncol(chem.lme.caseG$ranef.delete[[1]]), 1 + ncol(ranef(chem.lme)$lea)) #lme doesn't add delete column
  
  n <- length(unique(chem.lmer@flist[["school:lea"]]))
  expect_equal(nrow(chem.lmer.caseG2$ranef.delete[[1]]), n * (n-1))  
  expect_equal(ncol(chem.lmer.caseG2$ranef.delete[[1]]), 2 + ncol(ranef(chem.lmer)$'school:lea'))
  
  n <- length(unique(chem.lme$groups[["school"]]))
  expect_equal(nrow(chem.lme.caseG2$ranef.delete[[2]]), n * (n-1))
  expect_equal(ncol(chem.lme.caseG2$ranef.delete[[2]]), 1 + ncol(ranef(chem.lme)$school))
  
})

#passed
test_that("Dimensions of variance covariance matrices after deletion are correct for single case deletion", {
  #sleepstudy
  expect_equal(length(sleep.lmer.case$vcov.delete), nrow(sleep.lmer@frame))
  #should I check the first one, a random sample, all of them, something else??? For now, just the first
  expect_equal(nrow(sleep.lmer.case$vcov.delete[[1]]), nrow(as.matrix(vcov(sleep.lmer))))
  expect_equal(ncol(sleep.lmer.case$vcov.delete[[1]]), ncol(as.matrix(vcov(sleep.lmer))))
  expect_equal(nrow(sleep.lmer.caseD$vcov.delete), nrow(as.matrix(vcov(sleep.lmer))))
  expect_equal(ncol(sleep.lmer.caseD$vcov.delete), ncol(as.matrix(vcov(sleep.lmer))))
  
  
  expect_equal(length(sleep.lme.case$vcov.delete), nrow(sleep.lme$data))
  expect_equal(nrow(sleep.lme.case$vcov.delete[[1]]), nrow(as.matrix(vcov(sleep.lme))))
  expect_equal(ncol(sleep.lme.case$vcov.delete[[1]]), ncol(as.matrix(vcov(sleep.lme))))
  expect_equal(nrow(sleep.lme.caseD$vcov.delete), nrow(as.matrix(vcov(sleep.lme))))
  expect_equal(ncol(sleep.lme.caseD$vcov.delete), ncol(as.matrix(vcov(sleep.lme))))
  
  #chemistry
  expect_equal(length(chem.lmer.case$vcov.delete), nrow(chem.lmer@frame))
  expect_equal(nrow(chem.lmer.case$vcov.delete[[1]]), nrow(as.matrix(vcov(chem.lmer))))
  expect_equal(ncol(chem.lmer.case$vcov.delete[[1]]), ncol(as.matrix(vcov(chem.lmer))))
  
  expect_equal(length(chem.lme.case$vcov.delete), nrow(chem.lme$data))
  expect_equal(nrow(chem.lme.case$vcov.delete[[1]]), nrow(as.matrix(vcov(chem.lme))))
  expect_equal(ncol(chem.lme.case$vcov.delete[[1]]), ncol(as.matrix(vcov(chem.lme))))
  
})

#passed
test_that("Dimensions of variance covariance matrices after deletion are correct for group deletion", {
  
  #sleepstudy
  expect_equal(length(sleep.lmer.caseG$vcov.delete), length(unique(sleep.lmer@flist[["Subject"]])))
  #should I check the first one, a random sample, all of them, something else??? For now, just the first
  expect_equal(nrow(sleep.lmer.caseG$vcov.delete[[1]]), nrow(as.matrix(vcov(sleep.lmer))))
  expect_equal(ncol(sleep.lmer.caseG$vcov.delete[[1]]), ncol(as.matrix(vcov(sleep.lmer))))
  expect_equal(nrow(sleep.lmer.caseGD$vcov.delete), nrow(as.matrix(vcov(sleep.lmer))))
  expect_equal(ncol(sleep.lmer.caseGD$vcov.delete), ncol(as.matrix(vcov(sleep.lmer))))
  
  
  expect_equal(length(sleep.lme.caseG$vcov.delete), length(unique(sleep.lme$groups$Subject)))
  expect_equal(nrow(sleep.lme.caseG$vcov.delete[[1]]), nrow(as.matrix(vcov(sleep.lme))))
  expect_equal(ncol(sleep.lme.caseG$vcov.delete[[1]]), ncol(as.matrix(vcov(sleep.lme))))
  expect_equal(nrow(sleep.lme.caseGD$vcov.delete), nrow(as.matrix(vcov(sleep.lme))))
  expect_equal(ncol(sleep.lme.caseGD$vcov.delete), ncol(as.matrix(vcov(sleep.lme))))
  
  #chemistry
  expect_equal(length(chem.lmer.caseG$vcov.delete), length(unique(chem.lmer@flist[["lea"]])))
  expect_equal(nrow(chem.lmer.caseG$vcov.delete[[1]]), nrow(as.matrix(vcov(chem.lmer))))
  expect_equal(ncol(chem.lmer.caseG$vcov.delete[[1]]), ncol(as.matrix(vcov(chem.lmer))))
  
  expect_equal(length(chem.lme.caseG$vcov.delete), length(unique(chem.lme$groups[["lea"]])))
  expect_equal(nrow(chem.lme.caseG$vcov.delete[[1]]), nrow(as.matrix(vcov(chem.lme))))
  expect_equal(ncol(chem.lme.caseG$vcov.delete[[1]]), ncol(as.matrix(vcov(chem.lme))))
  
  expect_equal(length(chem.lmer.caseG2$vcov.delete), length(unique(chem.lmer@flist[["school:lea"]])))
  expect_equal(nrow(chem.lmer.caseG2$vcov.delete[[1]]), nrow(as.matrix(vcov(chem.lmer))))
  expect_equal(ncol(chem.lmer.caseG2$vcov.delete[[1]]), ncol(as.matrix(vcov(chem.lmer))))
  
  expect_equal(length(chem.lme.caseG2$vcov.delete), length(unique(chem.lme$groups[["school"]])))
  expect_equal(nrow(chem.lme.caseG2$vcov.delete[[1]]), nrow(as.matrix(vcov(chem.lme))))
  expect_equal(ncol(chem.lme.caseG2$vcov.delete[[1]]), ncol(as.matrix(vcov(chem.lme))))
  
  
  
})

#passed
test_that("Dimensions of fitted values after deletion are correct for single case deletion", {
  #number of rows is number of observations * (number of observations - 1), number of columns is 2 + variables
  
  #sleepstudy
  n <- nrow(sleep.lmer@frame)
  expect_equal(nrow(sleep.lmer.case$fitted.delete), n * (n-1))
  expect_equal(ncol(sleep.lmer.case$fitted.delete), 2 + ncol(sleep.lmer@frame))
  expect_equal(length(sleep.lmer.caseD$fitted.delete), n -3)
  
  n <- nrow(sleep.lme$data)
  expect_equal(nrow(sleep.lme.case$fitted.delete), n * (n-1))
  expect_equal(ncol(sleep.lme.case$fitted.delete), 2 + ncol(sleep.lme$data))
  expect_equal(length(sleep.lme.caseD$fitted.delete), n -3)
  
  #chemistry 
  n <- nrow(chem.lmer@frame)
  expect_equal(nrow(chem.lmer.case$fitted.delete), n * (n-1))
  expect_equal(ncol(chem.lmer.case$fitted.delete), 2 + ncol(chem.lmer@frame))
  
  n <- nrow(chem.lme$data)
  expect_equal(nrow(chem.lme.case$fitted.delete), n * (n-1))
  expect_equal(ncol(chem.lme.case$fitted.delete), 2 + ncol(chem.lmer@frame)) #want to match lme4

  
  
})

#passed 
test_that("Dimensions of fitted values after deletion are correct for group deletion", {
  #number of rows is number of observations * (number of groups - 1), number of columns is 2 + variables
  
  #sleepstudy
  nobs <- nrow(sleep.lmer@frame)
  ngroups <- length(unique(sleep.lmer@flist[["Subject"]]))
  expect_equal(nrow(sleep.lmer.caseG$fitted.delete), nobs * (ngroups - 1))
  expect_equal(ncol(sleep.lmer.caseG$fitted.delete), 2 + ncol(sleep.lmer@frame))
  expect_equal(length(sleep.lmer.caseGD$fitted.delete), nobs - 10)
  
  
  nobs <- nrow(sleep.lme$data)
  ngroups <- length(unique(sleep.lme$groups$Subject))
  expect_equal(nrow(sleep.lme.caseG$fitted.delete), nobs * (ngroups - 1))
  expect_equal(ncol(sleep.lme.caseG$fitted.delete), 2 + ncol(sleep.lme$data)) 
  expect_equal(length(sleep.lme.caseGD$fitted.delete), nobs - 10)
  
  #chemistry 
  nobs <- nrow(chem.lmer@frame)
  ngroups <- length(unique(chem.lmer@flist[["lea"]]))
  expect_equal(nrow(chem.lmer.caseG$fitted.delete), nobs * (ngroups - 1))
  expect_equal(ncol(chem.lmer.caseG$fitted.delete), 2 + ncol(chem.lmer@frame))
  
  nobs <- nrow(chem.lme$data)
  ngroups <- length(unique(chem.lme$groups[["lea"]]))
  expect_equal(nrow(chem.lme.caseG$fitted.delete), nobs * (ngroups - 1))
  expect_equal(ncol(chem.lme.caseG$fitted.delete), 2 + ncol(chem.lmer@frame)) #this should be equal to lmer numbers
  
  ngroups <- length(unique(chem.lme$groups[["school"]]))
  expect_equal(nrow(chem.lme.caseG2$fitted.delete), nobs * (ngroups -1))
  expect_equal(ncol(chem.lme.caseG2$fitted.delete), 2 + ncol(chem.lmer@frame))
})

#passed
test_that("Dimenstions of variance components are correct for single case deletion", {
  
  #sleepstudy 
  expect_equal(length(sleep.lmer.case$varcomp.delete), nrow(sleep.lmer@frame))
  #same issue here, should I check all of them, some of them, only the first one? 
  expect_equal(length(sleep.lmer.case$varcomp.delete[[1]]), length(varcomp.mer(sleep.lmer)))
  expect_equal(length(sleep.lmer.caseD$varcomp.delete), length(varcomp.mer(sleep.lmer)))
  
  expect_equal(length(sleep.lme.case$varcomp.delete), nrow(sleep.lme$data))
  expect_equal(length(sleep.lme.case$varcomp.delete[[1]]), length(varcomp.lme(sleep.lme)))
  expect_equal(length(sleep.lme.caseD$varcomp.delete), length(varcomp.lme(sleep.lme)))
  
  #chemistry 
  expect_equal(length(chem.lmer.case$varcomp.delete), nrow(chem.lmer@frame))
  expect_equal(length(chem.lmer.case$varcomp.delete[[1]]), length(varcomp.mer(chem.lmer)))
  
  expect_equal(length(chem.lme.case$varcomp.delete), nrow(chem.lme$data))
  expect_equal(length(chem.lme.case$varcomp.delete[[1]]), length(varcomp.lme(chem.lme)))
  
})

#passed - here 
test_that("Dimensions of variance components are correct for group deletion", {
  #sleepstudy
  expect_equal(length(sleep.lmer.caseG$varcomp.delete), length(unique(sleep.lmer@flist[["Subject"]])))
  expect_equal(length(sleep.lmer.caseG$varcomp.delete[[1]]), length(varcomp.mer(sleep.lmer)))
  expect_equal(length(sleep.lmer.caseGD$varcomp.delete), length(varcomp.mer(sleep.lmer)))
  
  expect_equal(length(sleep.lme.caseG$varcomp.delete), length(unique(sleep.lme$groups$Subject)))
  expect_equal(length(sleep.lme.caseG$varcomp.delete[[1]]), length(varcomp.lme(sleep.lme)))
  expect_equal(length(sleep.lme.caseGD$varcomp.delete), length(varcomp.lme(sleep.lme)))
  
  #chemistry 
  expect_equal(length(chem.lmer.caseG$varcomp.delete), length(unique(chem.lmer@flist[["lea"]])))
  expect_equal(length(chem.lmer.caseG$varcomp.delete[[1]]), length(varcomp.mer(chem.lmer)))
  
  expect_equal(length(chem.lmer.caseG2$varcomp.delete), length(unique(chem.lmer@flist[["school:lea"]])))
  expect_equal(length(chem.lmer.caseG2$varcomp.delete[[1]]), length(varcomp.mer(chem.lmer)))
  
  expect_equal(length(chem.lme.caseG$varcomp.delete), length(unique(chem.lme$groups[["lea"]])))
  expect_equal(length(chem.lme.caseG$varcomp.delete[[1]]), length(varcomp.lme(chem.lme)))
  
  expect_equal(length(chem.lme.caseG2$varcomp.delete), length(unique(chem.lme$groups[["school"]])))
  expect_equal(length(chem.lme.caseG2$varcomp.delete[[1]]), length(varcomp.lme(chem.lme)))
})



#I should be checking all of these things also when:
#we set level - DONE for sleepstudy 
#we set delete - DONE for sleepstudy
#all tests pass for sleepstudy lmer and lme 

#does delete throw errors at the right spots? (especially for three level models, indices or group names allowed)
#especially check that three level lme models allow the correct names for delete and for level 

#are varcomp.mer and varcomp.lme exported? If I use those functions in the tests, will it work? 



