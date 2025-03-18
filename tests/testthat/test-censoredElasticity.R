test_that("Testing censoredElasticity function.", {

  testthat::skip_on_cran()

  # Reference: https://github.com/NoeNava-USDA/mexSugarTax_replication
  # Nava, Noé J., and Diansheng Dong. 2022. "The impact of taxing sugary-sweetened beverages in México: A censored QUAI demand system approach." Journal of Agricultural and Applied Economics Association, 1(1):1-23:
  # https://onlinelibrary.wiley.com/doi/10.1002/jaa2.6

  testing_data <- readRDS(testthat::test_path('testing_data.rds'))

  s1 <- testing_data$s1
  s2 <- testing_data$s2
  s3 <- testing_data$s3
  s4 <- testing_data$s4
  s5 <- testing_data$s5
  s6 <- testing_data$s6

  lnp1 <- testing_data$lnp1
  lnp2 <- testing_data$lnp2
  lnp3 <- testing_data$lnp3
  lnp4 <- testing_data$lnp4
  lnp5 <- testing_data$lnp5
  lnp6 <- testing_data$lnp6

  age <- testing_data$age
  size <- testing_data$size
  sex <- testing_data$sex
  educ <- testing_data$educ

  # Alpha
  b0 <- rep(0, 5)

  # Beta
  b0 <- c(b0, rep(0.003, 5))

  # Gamma
  b0 <- c(b0,0.01,0,0.01,0,0, 0.01,0,0,0,0.01,0,0,0,0,0.01)

  # Demos
  b0 <- c(b0,rep(0.002, 20))

  # Sigma
  b0 <- c(b0,1,0,1,0,0,1,0,0,0,1,0,0,0,0,1)

  list_etas <- censoredElasticity(
  Params = b ,
  Shares = matrix(c(s1, s2, s3, s4, s5, s6), ncol = 6),
  Prices = matrix(c(lnp1, lnp2, lnp3, lnp4, lnp5, lnp6), ncol = 6),
  Budget = matrix(testing_data$lnw),
  Demographics = matrix(c(age, size, educ, sex), ncol = 4),
  quaids = FALSE,
  func = mean,
  na.rm = TRUE,
  vcov = vcov
  )

  etas <- list_etas$Elasticities
  E_Uobs <- list_etas$E_Uobs

  # Checks for homogeneity
  expect_equal(round(as.numeric(rowSums(etas), 6)), rep(0, nrow(etas)))

  # Engel Aggregation
  expect_equal(round(as.numeric(t(etas[,-7]) %*% E_Uobs),6), round(-E_Uobs, 6))


})
