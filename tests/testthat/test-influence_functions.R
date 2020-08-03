context("tests for all influence functions")

#sleepstudy models - 2 level
data(sleepstudy, package = 'lme4')
sleep.lmer <- lme4::lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy)
sleep.lme <- nlme::lme(Reaction ~ Days, random =  ~ Days|Subject, data = sleepstudy)
nsleep <- nrow(sleep.lmer@frame)
sleep.ngroups <- length(unique(sleep.lmer@flist[["Subject"]]))


#chemistry models - 3 level 
data(Chem97, package = "mlmRev")  
Chem97 <- Chem97[1:257,]
chem.lmer <- lme4::lmer(score ~ gcsecnt + (1|lea/school), data = Chem97)
chem.lme <- nlme::lme(score ~ gcsecnt, random = ~1|lea/school, data = Chem97)
nchem <- nrow(chem.lmer@frame)

chem.ngroups2 <- length(unique(chem.lmer@flist[["school:lea"]]))
chem.ngroups3 <- length(unique(chem.lmer@flist[["lea"]]))

test_that("Number of rows from leverage output is equal to number of observations for single case deletion", {
  #sleepstudy
  expect_equal(nrow(leverage(sleep.lmer)), nsleep)
  expect_equal(nrow(leverage(sleep.lme)), nsleep)
  
  #chemistry
  expect_equal(nrow(leverage(chem.lmer)), nchem)
  expect_equal(nrow(leverage(chem.lme)), nchem)
})

test_that("Number of rows from leverage output is equal to number of groups for group deletion", {
  expect_equal(nrow(leverage(sleep.lmer, level = "Subject")), sleep.ngroups)
  expect_equal(nrow(leverage(sleep.lme, level = "Subject")), sleep.ngroups)
  
  expect_equal(nrow(leverage(chem.lmer, level = "lea")), chem.ngroups3)
  expect_equal(nrow(leverage(chem.lmer, level = "school:lea")), chem.ngroups2)
  
  expect_equal(nrow(leverage(chem.lme, level = "lea")), chem.ngroups3)
  expect_equal(nrow(leverage(chem.lme, level = "school")), chem.ngroups2)
})


test_that("Number of rows from cook's distance output is equal to number of observations for single case deletion", {
  expect_equal(length(cooks.distance(sleep.lmer)), nsleep)
  expect_equal(length(cooks.distance(sleep.lme)), nsleep)
  
  expect_equal(length(cooks.distance(chem.lmer)), nchem)
  expect_equal(length(cooks.distance(chem.lme)), nchem)
})

test_that("Number of rows from cook's distance output is equal to number of groups for group deletion", {
  expect_equal(length(cooks.distance(sleep.lmer, level = "Subject")), sleep.ngroups)
  expect_equal(length(cooks.distance(sleep.lme, level = "Subject")), sleep.ngroups)
  
  expect_equal(length(cooks.distance(chem.lmer, level = "lea")), chem.ngroups3)
  expect_equal(length(cooks.distance(chem.lme, level = "lea")), chem.ngroups3)
  
  expect_equal(length(cooks.distance(chem.lmer, level = "school:lea")), chem.ngroups2)
  expect_equal(length(cooks.distance(chem.lme, level = "school")), chem.ngroups2)
})

test_that("Number of rows from MDFFITS output is equal to number of observations for single case deletion", {
  expect_equal(length(mdffits(sleep.lmer)), nsleep)
  expect_equal(length(mdffits(sleep.lme)), nsleep)
  
  expect_equal(length(mdffits(chem.lmer)), nchem)
  expect_equal(length(mdffits(chem.lme)), nchem)
})

test_that("Number of rows from MDFFITS output is equal to number of groups for group deletion", {
  expect_equal(length(mdffits(sleep.lmer, level = "Subject")), sleep.ngroups)
  expect_equal(length(mdffits(sleep.lme, level = "Subject")), sleep.ngroups)
  
  expect_equal(length(mdffits(chem.lmer, level = "lea")), chem.ngroups3)
  expect_equal(length(mdffits(chem.lme, level = "lea")), chem.ngroups3)
  
  expect_equal(length(mdffits(chem.lmer, level = "school:lea")), chem.ngroups2)
  expect_equal(length(mdffits(chem.lme, level = "school")), chem.ngroups2)
})

test_that("Number of rows from covtrace output is equal to number of observations for single case deletion", {
  expect_equal(length(covtrace(sleep.lmer)), nsleep)
  expect_equal(length(covtrace(sleep.lme)), nsleep)
  
  expect_equal(length(covtrace(chem.lmer)), nchem)
  expect_equal(length(covtrace(chem.lme)), nchem)
})

test_that("Number of rows from covtrace output is equal to number of groups for group deletion", {
  expect_equal(length(covtrace(sleep.lmer, level = "Subject")), sleep.ngroups)
  expect_equal(length(covtrace(sleep.lme, level = "Subject")), sleep.ngroups)
  
  expect_equal(length(covtrace(chem.lmer, level = "lea")), chem.ngroups3)
  expect_equal(length(covtrace(chem.lme, level = "lea")), chem.ngroups3)
  
  expect_equal(length(covtrace(chem.lmer, level = "school:lea")), chem.ngroups2)
  expect_equal(length(covtrace(chem.lme, level = "school")), chem.ngroups2)
})

test_that("Number of rows from covratio output is equal to number of observations for single case deletion", {
  expect_equal(length(covratio(sleep.lmer)), nsleep)
  expect_equal(length(covratio(sleep.lme)), nsleep)
  
  expect_equal(length(covratio(chem.lmer)), nchem)
  expect_equal(length(covratio(chem.lme)), nchem)
})

test_that("Number of rows from covratio output is equal to number of groups for group deletion", {
  expect_equal(length(covratio(sleep.lmer, level = "Subject")), sleep.ngroups)
  expect_equal(length(covratio(sleep.lme, level = "Subject")), sleep.ngroups)
  
  expect_equal(length(covratio(chem.lmer, level = "lea")), chem.ngroups3)
  expect_equal(length(covratio(chem.lme, level = "lea")), chem.ngroups3)
  
  expect_equal(length(covratio(chem.lmer, level = "school:lea")), chem.ngroups2)
  expect_equal(length(covratio(chem.lme, level = "school")), chem.ngroups2)
})

