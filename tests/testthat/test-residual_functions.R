library(lme4)

fm_slopes <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

test_that("mahalanobis_ranef returns one value per group", {
  result <- HLMdiag:::mahalanobis_ranef.lmerMod(fm_slopes)
  expect_length(result, length(unique(sleepstudy$Subject)))
})

test_that("mahalanobis_ranef returns finite non-negative values", {
  result <- HLMdiag:::mahalanobis_ranef.lmerMod(fm_slopes)
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})

test_that("mahalanobis_ranef works for models with more than 2 random effects", {
  fm3re <- suppressWarnings(lmer(Reaction ~ Days + (Days + I(Days^2) | Subject), sleepstudy))
  result <- HLMdiag:::mahalanobis_ranef.lmerMod(fm3re)
  expect_length(result, length(unique(sleepstudy$Subject)))
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
})
