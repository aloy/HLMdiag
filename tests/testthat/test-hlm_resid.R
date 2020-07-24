library(HLMdiag)
bdf <- nlme::bdf

test_that("detects level 2 variables, lme4", {
  bdf.lmer <- lme4::lmer(IQ.verb ~ ses + aritPOST + langPOST + schoolSES + 
                           (1|schoolNR), data = bdf)
  bdf.resids <- hlm_resid(bdf.lmer, level = "schoolNR")
  expect_equal(ncol(bdf.resids), 4)
  expect_equal(names(bdf.resids)[1:2], c("schoolNR", "schoolSES"))
})

test_that("detects level 2 variables, nlme", {
  bdf.lme <- nlme::lme(IQ.verb ~ sex + ses + denomina + schoolSES, 
                       random = ~1|schoolNR, data = bdf)
  bdf.resids <- hlm_resid(bdf.lme, level = "schoolNR")
  expect_equal(ncol(bdf.resids), 4)
  expect_equal(names(bdf.resids)[1:2], c("schoolNR", "schoolSES"))
})