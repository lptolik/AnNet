library(AnNet)
library(testthat)
data(karate,package='igraphdata')
#context("Testing robustness calculations")
test_that("proper algorithm selected",{
  expect_error(getRobustness(karate,'alg',matrix(0,ncol = 1,nrow = 1)),"Membership.*")
  expect_error(getClustering(karate,'alg'),".*lec*")
})
