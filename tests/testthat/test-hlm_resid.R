library(dplyr)
library(lme4)
library(nlme)

#2 level, random intercept
bdf <- nlme::bdf
bdf.lmer <- lme4::lmer(IQ.verb ~ ses + aritPOST + langPOST + schoolSES + 
                         (1|schoolNR), data = bdf)
bdf.lme <- nlme::lme(IQ.verb ~ ses + aritPOST + langPOST + schoolSES, 
                     random = ~1|schoolNR, data = bdf)

test_that("correct dimentions, lme4", {
  #check lvl1
  bdf.resids <- hlm_resid(bdf.lmer)
  expect_equal(nrow(bdf.resids), nrow(bdf))
  #check lvl2
  bdf.resids <- hlm_resid(bdf.lmer, level = "schoolNR")
  expect_equal(nrow(bdf.resids), length(unique(bdf$schoolNR)))
})

test_that("correct dimentions, nlme", {
  #check lvl1
  bdf.resids <- hlm_resid(bdf.lme)
  expect_equal(nrow(bdf.resids), nrow(bdf))
  #check lvl2
  bdf.resids <- hlm_resid(bdf.lme, level = "schoolNR")
  expect_equal(nrow(bdf.resids), length(unique(bdf$schoolNR)))
})

test_that("standardize works, lme4", {
  bdf.resids.raw <- hlm_resid(bdf.lmer, standardize = F)
  bdf.resids.std <- hlm_resid(bdf.lmer, standardize = T)
  nc <- ncol(bdf.resids.raw)
  #check column names
  expect_equal(names(bdf.resids.raw[(nc-5):nc]),
               c(".resid", ".fitted", ".ls.resid", ".ls.fitted", ".mar.resid", ".mar.fitted"))
  expect_equal(names(bdf.resids.std[(nc-5):nc]),
               c(".std.resid", ".fitted", ".std.ls.resid", ".ls.fitted", ".chol.mar.resid", ".mar.fitted"))
  
  #check that raw is not the same as std
  expect_false(all(bdf.resids.raw$.resid == bdf.resids.std$.std.resid))
  expect_false(all(bdf.resids.raw$.ls.resid == bdf.resids.std$.std.ls.resid))
  expect_false(all(bdf.resids.raw$.mar.resid == bdf.resids.std$.chol.mar.resid))
  expect_true(all(bdf.resids.raw$.fitted == bdf.resids.std$.fitted))
})

test_that("standardize works, nlme", {
  bdf.resids.raw <- hlm_resid(bdf.lme, standardize = F)
  bdf.resids.std <- hlm_resid(bdf.lme, standardize = T)
  nc <- ncol(bdf.resids.raw)
  #check column names
  expect_equal(names(bdf.resids.raw[(nc-5):nc]),
               c(".resid", ".fitted", ".ls.resid", ".ls.fitted", ".mar.resid", ".mar.fitted"))
  expect_equal(names(bdf.resids.std[(nc-5):nc]),
               c(".std.resid", ".fitted", ".std.ls.resid", ".ls.fitted", ".chol.mar.resid", ".mar.fitted"))
  
  #check that raw is not the same as std
  expect_false(all(bdf.resids.raw$.resid == bdf.resids.std$.std.resid))
  expect_false(all(bdf.resids.raw$.ls.resid == bdf.resids.std$.std.ls.resid))
  expect_false(all(bdf.resids.raw$.mar.resid == bdf.resids.std$.chol.mar.resid))
  expect_true(all(bdf.resids.raw$.fitted == bdf.resids.std$.fitted))
})

test_that("detects level 2 variables, lme4", {
  bdf.resids <- hlm_resid(bdf.lmer, level = "schoolNR")
  #intercept random effect and both fixed effect terms
  expect_equal(ncol(bdf.resids), 4)
  #grabs only two level variables
  expect_equal(names(bdf.resids)[1:2], c("schoolNR", "schoolSES"))
})

test_that("detects level 2 variables, nlme", {
  bdf.resids <- hlm_resid(bdf.lme, level = "schoolNR")
  #intercept random effect and both fixed effect terms
  expect_equal(ncol(bdf.resids), 4)
  #grabs only two level variables
  expect_equal(names(bdf.resids)[1:2], c("schoolNR", "schoolSES"))
})

#3 level, random intercept and slope
class <- read.csv("http://www-personal.umich.edu/~bwest/classroom.csv")
class.lmer <- lmer(mathgain ~ mathkind + minority + ses + housepov + 
                     (mathkind | schoolid/classid), class)
class.lme <- lme(mathgain ~ mathkind + minority + ses + housepov,
                     random = ~mathkind | schoolid/classid, class)

test_that("3 level model tests, lme4", {
  #inner level
  bdf.resids <- hlm_resid(class.lmer, level = "classid:schoolid")
  #intercept and slope random effect and both fixed effect terms
  expect_equal(ncol(bdf.resids), 8)
  #grabs only two level variables
  expect_equal(names(bdf.resids)[1:4], c("group", "classid", "schoolid", "housepov"))
  
  #highest level
  bdf.resids <- hlm_resid(class.lmer, level = "schoolid")
  #intercept and slope random effect and both fixed effect terms
  expect_equal(ncol(bdf.resids), 6)
  #grabs only two level variables
  expect_equal(names(bdf.resids)[1:2], c("schoolid", "housepov"))
})

test_that("3 level model tests, nlme", {
  #inner level
  bdf.resids <- hlm_resid(class.lme, level = "classid")
  #intercept and slope random effect and both fixed effect terms
  expect_equal(ncol(bdf.resids), 8)
  #grabs only two level variables
  expect_equal(names(bdf.resids)[1:4], c("group", "schoolid", "classid", "housepov"))
  
  #highest level
  bdf.resids <- hlm_resid(class.lme, level = "schoolid")
  #intercept and slope random effect and both fixed effect terms
  expect_equal(ncol(bdf.resids), 6)
  #grabs only two level variables
  expect_equal(names(bdf.resids)[1:2], c("schoolid", "housepov"))
})
