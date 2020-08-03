context("tests for hlm_influence")

#sleepstudy
sleepstudy <- data(sleepstudy, package = 'lme4')
sleep.lmer <- lme4::lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy)
sleep.lme <- nlme::lme(Reaction ~ Days, random =  ~ Days|Subject, data = sleepstudy)

#chemistry 
data(Chem97, package = "mlmRev")  
Chem97 <- Chem97[1:257,]
chem.lmer <- lme4::lmer(score ~ gcsecnt + (1|lea/school), data = Chem97)
chem.lme <- nlme::lme(score ~ gcsecnt, random = ~1|lea/school, data = Chem97)


#failing right now, probably due to issues with the frame 
test_that("Number of rows and columns is correct for default approximations and full refits", {
  #sleepstudy
  sleep.lmer.infl <- hlm_influence(sleep.lmer)
  expect_equal(ncol(sleep.lmer.infl), 5 + ncol(sleep.lmer@frame))
  expect_equal(nrow(sleep.lmer.infl), nrow(sleep.lmer@frame))
  
  sleep.lme.infl <- hlm_influence(sleep.lme)
  expect_equal(ncol(sleep.lme.infl), 5 + ncol(sleep.lme$data))
  expect_equal(nrow(sleep.lme.infl), nrow(sleep.lme$data))
  
  sleep.lmer.infl2 <- hlm_influence(sleep.lmer, approx = FALSE)
  expect.equal(ncol(sleep.lmer.infl2), 5 + ncol(sleep.lmer@frame) + ncol(rvc(sleep.lmer)))
  
  sleep.lme.infl2 <- hlm_influence(sleep.lme, approx = FALSE)
  expect.equal(ncol(sleep.lme.infl2), 5 + ncol(sleep.lmer@frame) + ncol(rvc(sleep.lme)))
  
  #chemistry
  chem.lmer.infl <- hlm_influence(chem.lmer)
  expect_equal(ncol(chem.lmer.infl), 5 + ncol(chem.lmer@frame)) #should match lme4, this will fail for now 
  expect_equal(nrow(chem.lmer.infl), nrow(chem.lmer@frame))
  
  chem.lme.infl <- hlm_influence(chem.lme)
  expect_equal(ncol(chem.lme.infl), 5 + ncol(chem.lmer@frame)) #same here 
  expect_equal(nrow(chem.lme.infl), nrow(chem.lme$data))
  
  chem.lmer.infl2 <- hlm_influence(chem.lmer, approx = FALSE)
  expect.equal(ncol(chem.lmer.infl2), 5 + ncol(chem.lmer@frame) + ncol(rvc(chem.lmer)))
  
  chem.lme.infl2 <- hlm_influence(chem.lme, approx = FALSE)
  expect.equal(ncol(chem.lme.infl2), 5 + ncol(chem.lmer@frame) + ncol(rvc(chem.lme)))

})

#passes everything that should pass (not last line, until issue is fixed)
test_that("Number of columns is correct when leverage is specified", {
  #sleepstudy
  sleep.lmer.infl <- hlm_influence(sleep.lmer, leverage = c("overall", "fixef", "ranef", "ranef.uc"))
  expect_equal(ncol(sleep.lmer.infl), 8 + ncol(sleep.lmer@frame))
  
  sleep.lme.infl <- hlm_influence(sleep.lme, leverage = c("overall", "fixef", "ranef", "ranef.uc"))
  expect_equal(ncol(sleep.lme.infl), 8 + ncol(sleep.lmer@frame))
  
  #chemistry 
  chem.lmer.infl <- hlm_influence(chem.lmer, leverage = c("overall", "fixef", "ranef", "ranef.uc"))
  expect_equal(ncol(chem.lmer.infl), 8 + ncol(chem.lmer@frame))
  
  chem.lme.infl <- hlm_influence(chem.lme, leverage = c("overall", "fixef", "ranef", "ranef.uc"))
  expect_equal(ncol(chem.lme.infl), 8 + ncol(chem.lmer@frame)) 
})

#passed
test_that("Number of rows is correct when level is specified", {
  #sleepstudy
  sleep.lmer.infl <- hlm_influence(sleep.lmer, level = "Subject")
  expect_equal(nrow(sleep.lmer.infl), length(unique(sleep.lmer@flist[["Subject"]])))
  
  sleep.lme.infl <- hlm_influence(sleep.lme, level = "Subject")
  expect_equal(nrow(sleep.lme.infl), length(unique(sleep.lme$groups[["Subject"]])))
  
  #chemistry
  chem.lmer.infl <- hlm_influence(chem.lmer, level = "lea")
  chem.lme.infl <- hlm_influence(chem.lme, level = "lea")
  
  expect_equal(nrow(chem.lmer.infl), length(unique(chem.lmer@flist[["lea"]])))
  expect_equal(nrow(chem.lme.infl), length(unique(chem.lme$groups[["lea"]])))
  
  chem.lmer.infl2 <- hlm_influence(chem.lmer, level = "school:lea")
  chem.lme.infl2 <- hlm_influence(chem.lme, level = "school")
  
  expect_equal(nrow(chem.lmer.infl2), length(unique(chem.lmer@flist[["school:lea"]])))
  expect_equal(nrow(chem.lme.infl2), length(unique(chem.lme$groups[["school"]])))
})

#passed
test_that("Influence diagnostic columns match output from influence functions", {
  #sleepstudy
  sleep.lmer.infl <- hlm_influence(sleep.lmer)
  expect_equal(sleep.lmer.infl$cooksd, cooks.distance(sleep.lmer))
  expect_equal(sleep.lmer.infl$mdffits, mdffits(sleep.lmer))
  expect_equal(sleep.lmer.infl$covtrace, covtrace(sleep.lmer))
  expect_equal(sleep.lmer.infl$covratio, covratio(sleep.lmer))
  expect_equal(sleep.lmer.infl$leverage.overall, leverage(sleep.lmer)[,1])
  
  sleep.lme.infl <- hlm_influence(sleep.lme)
  expect_equal(sleep.lme.infl$cooksd, cooks.distance(sleep.lme))
  expect_equal(sleep.lme.infl$mdffits, mdffits(sleep.lme))
  expect_equal(sleep.lme.infl$covtrace, covtrace(sleep.lme))
  expect_equal(sleep.lme.infl$covratio, covratio(sleep.lme))
  expect_equal(sleep.lme.infl$leverage.overall, leverage(sleep.lme)[,1])
  
  #chemistry
  chem.lmer.infl <- hlm_influence(chem.lmer)
  expect_equal(chem.lmer.infl$cooksd, cooks.distance(chem.lmer))
  expect_equal(chem.lmer.infl$mdffits, mdffits(chem.lmer))
  expect_equal(chem.lmer.infl$covtrace, covtrace(chem.lmer))
  expect_equal(chem.lmer.infl$covratio, covratio(chem.lmer))
  expect_equal(chem.lmer.infl$leverage.overall, leverage(chem.lmer)[,1])
  
  sleep.lme.infl <- hlm_influence(chem.lme)
  expect_equal(chem.lme.infl$cooksd, cooks.distance(chem.lme))
  expect_equal(chem.lme.infl$mdffits, mdffits(chem.lme))
  expect_equal(chem.lme.infl$covtrace, covtrace(chem.lme))
  expect_equal(chem.lme.infl$covratio, covratio(chem.lme))
  expect_equal(chem.lme.infl$leverage.overall, leverage(chem.lme[,1]))
  
})

#passed
test_that("Number of rows and columns are correct when delete is specified", {
  #sleepstudy
  sleep.lmer.infl <- hlm_influence(sleep.lmer, delete = c(1,10,13))
  expect_equal(nrow(sleep.lmer.infl), 1)
  expect_equal(ncol(sleep.lmer.infl), 4)
  
  sleep.lme.infl <- hlm_influence(sleep.lme, delete = c(1,10,13))
  expect_equal(nrow(sleep.lme.infl), 1)
  expect_equal(ncol(sleep.lme.infl), 4)
  
  expect_warning(hlm_influence(sleep.lmer, delete = c(1,10,13), leverage = "ranef"))
  expect_warning(hlm_influence(sleep.lme, delete = c(1,10,13), leverage = "ranef"))
  
  #chemistry 
  chem.lmer.infl <- hlm_influence(chem.lmer, delete = c(2, 8, 78))
  expect_equal(nrow(chem.lmer.infl), 1)
  expect_equal(ncol(chem.lmer.infl), 4)
  
  chem.lme.infl <- hlm_influence(chem.lme, delete = c(2, 8, 78))
  expect_equal(nrow(chem.lme.infl), 1)
  expect_equal(ncol(chem.lme.infl), 4)
})

