#2 level, random intercept
bdf <- nlme::bdf
bdf.lmer <- lme4::lmer(IQ.verb ~ ses + aritPOST + langPOST + schoolSES + 
                         (1|schoolNR), data = bdf)
bdf.lme <- nlme::lme(IQ.verb ~ ses + aritPOST + langPOST + schoolSES, 
                     random = ~1|schoolNR, data = bdf)

bdf.resids.lmer.lvl1 <- hlm_resid(bdf.lmer)
bdf.resids.lmer.lvl2 <- hlm_resid(bdf.lmer, level = "schoolNR")
bdf.resids.lmer.lvl1.std <- hlm_resid(bdf.lmer, standardize = T)

bdf.resids.lme.lvl1 <- hlm_resid(bdf.lme)
bdf.resids.lme.lvl2 <- hlm_resid(bdf.lme, level = "schoolNR")
bdf.resids.lme.lvl1.std <- hlm_resid(bdf.lme, standardize = T)


test_that("correct dimentions, lme4", {
  #check lvl1
  expect_equal(nrow(bdf.resids.lmer.lvl1), nrow(bdf))
  #check lvl2
  expect_equal(nrow(bdf.resids.lmer.lvl2), length(unique(bdf$schoolNR)))
})

test_that("correct dimentions, nlme", {
  #check lvl1
  expect_equal(nrow(bdf.resids.lme.lvl1), nrow(bdf))
  #check lvl2
  expect_equal(nrow(bdf.resids.lme.lvl2), length(unique(bdf$schoolNR)))
})




test_that("standardize works, lme4", {
  nc <- ncol(bdf.resids.lmer.lvl1)
  #check column names
  expect_equal(names(bdf.resids.lmer.lvl1[(nc-5):nc]),
               c(".resid", ".fitted", ".ls.resid", ".ls.fitted", ".mar.resid", ".mar.fitted"))
  expect_equal(names(bdf.resids.lmer.lvl1.std[(nc-5):nc]),
               c(".std.resid", ".fitted", ".std.ls.resid", ".ls.fitted", ".chol.mar.resid", ".mar.fitted"))
  
  #check that raw is not the same as std
  expect_false(all(bdf.resids.lmer.lvl1$.resid == bdf.resids.lmer.lvl1.std$.std.resid))
  expect_false(all(bdf.resids.lmer.lvl1$.ls.resid == bdf.resids.lmer.lvl1.std$.std.ls.resid))
  expect_false(all(bdf.resids.lmer.lvl1$.mar.resid == bdf.resids.lmer.lvl1.std$.chol.mar.resid))
  expect_true(all(bdf.resids.lmer.lvl1$.fitted == bdf.resids.lmer.lvl1.std$.fitted))
})

test_that("standardize works, nlme", {
  nc <- ncol(bdf.resids.lme.lvl1)
  #check column names
  expect_equal(names(bdf.resids.lme.lvl1[(nc-5):nc]),
               c(".resid", ".fitted", ".ls.resid", ".ls.fitted", ".mar.resid", ".mar.fitted"))
  expect_equal(names(bdf.resids.lme.lvl1.std[(nc-5):nc]),
               c(".std.resid", ".fitted", ".std.ls.resid", ".ls.fitted", ".chol.mar.resid", ".mar.fitted"))
  
  #check that raw is not the same as std
  expect_false(all(bdf.resids.lme.lvl1$.resid == bdf.resids.lme.lvl1.std$.std.resid))
  expect_false(all(bdf.resids.lme.lvl1$.ls.resid == bdf.resids.lme.lvl1.std$.std.ls.resid))
  expect_false(all(bdf.resids.lme.lvl1$.mar.resid == bdf.resids.lme.lvl1.std$.chol.mar.resid))
  expect_true(all(bdf.resids.lme.lvl1$.fitted == bdf.resids.lme.lvl1.std$.fitted))
})




test_that("detects level 2 variables, lme4", {
  #intercept random effect and both fixed effect terms
  expect_equal(ncol(bdf.resids.lmer.lvl2), 4)
  #grabs only two level variables
  expect_equal(names(bdf.resids.lmer.lvl2)[1:2], c("schoolNR", "schoolSES"))
})

test_that("detects level 2 variables, nlme", {
  #intercept random effect and both fixed effect terms
  expect_equal(ncol(bdf.resids.lme.lvl2), 4)
  #grabs only two level variables
  expect_equal(names(bdf.resids.lme.lvl2)[1:2], c("schoolNR", "schoolSES"))
})




#3 level, random intercept and slope
class <- read.csv("http://www-personal.umich.edu/~bwest/classroom.csv")
class.lmer <- lme4::lmer(mathgain ~ mathkind + minority + ses + housepov + 
                     (mathkind | schoolid/classid), class)
class.lme <- nlme::lme(mathgain ~ mathkind + minority + ses + housepov,
                     random = ~mathkind | schoolid/classid, class)

test_that("3 level model tests, lme4", {
  #inner level
  class.resids <- hlm_resid(class.lmer, level = "classid:schoolid")
  #intercept and slope random effect and both fixed effect terms
  expect_equal(ncol(class.resids), 8)
  #grabs only two level variables
  expect_equal(names(class.resids)[1:4], c("group", "classid", "schoolid", "housepov"))
  
  #highest level
  class.resids <- hlm_resid(class.lmer, level = "schoolid")
  #intercept and slope random effect and both fixed effect terms
  expect_equal(ncol(class.resids), 6)
  #grabs only two level variables
  expect_equal(names(class.resids)[1:2], c("schoolid", "housepov"))
})

test_that("3 level model tests, nlme", {
  #inner level
  class.resids <- hlm_resid(class.lme, level = "classid")
  #intercept and slope random effect and both fixed effect terms
  expect_equal(ncol(class.resids), 8)
  #grabs only two level variables
  expect_equal(names(class.resids)[1:4], c("group", "schoolid", "classid", "housepov"))
  
  #highest level
  class.resids <- hlm_resid(class.lme, level = "schoolid")
  #intercept and slope random effect and both fixed effect terms
  expect_equal(ncol(class.resids), 6)
  #grabs only two level variables
  expect_equal(names(class.resids)[1:2], c("schoolid", "housepov"))
})
