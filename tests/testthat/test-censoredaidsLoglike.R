test_that("censoredaidsLoglike spits out individual log-likelihood contribution in a matrix of nx1 -- important for gradient descent strategy!", {

  testthat::skip_on_cran()

  # Reference: https://github.com/NoeNava-USDA/mexSugarTax_replication
  # Nava, Noé J., and Diansheng Dong. 2022. "The impact of taxing sugary-sweetened beverages in México: A censored QUAI demand system approach." Journal of Agricultural and Applied Economics Association, 1(1):1-23:
  # https://onlinelibrary.wiley.com/doi/10.1002/jaa2.6

  loglikes <- readRDS(testthat::test_path('loglikes.rds'))
  testing_data <- readRDS(testthat::test_path('testing_data.rds'))

  full_alpha <- readRDS(testthat::test_path('testing_params/full_alpha.rds'))
  full_beta <- readRDS(testthat::test_path('testing_params/full_beta.rds'))
  full_gamma <- readRDS(testthat::test_path('testing_params/full_gamma.rds'))
  full_lamda <- readRDS(testthat::test_path('testing_params/full_lamda.rds'))
  full_theta <- t(readRDS(testthat::test_path('testing_params/full_theta.rds')))
  full_sigma <- readRDS(testthat::test_path('testing_params/full_sigma.rds'))
  full_sigma <- chol(full_sigma)

  li1 <- censoredaidsLoglike(
    Shares = matrix(c(testing_data$s1,
                      testing_data$s2,
                      testing_data$s3,
                      testing_data$s4), ncol = 4),
    Prices = matrix(c(testing_data$lnp1,
                      testing_data$lnp2,
                      testing_data$lnp3,
                      testing_data$lnp4), ncol = 4),
    Budget = matrix(testing_data$lnw, ncol = 1),
    ShareNames = c("SSB", "Juice", "Milk", "Water"),
    Demographics = matrix(c(testing_data$age,
                            testing_data$size,
                            testing_data$educ,
                            testing_data$sex), ncol = 4),
    DemographicNames = c("Age", "Size", "Sex", "Education"),
    Params = c(full_alpha[1:3],
               full_beta[1:3],
               full_gamma[1:3, 1:3][upper.tri(full_gamma[1:3, 1:3], diag = TRUE)],
               full_theta[1:4, 1:3],
               full_lamda[1:3],
               full_sigma[upper.tri(full_sigma, diag = TRUE)]),
    quaids = TRUE
  )
  li1[1:5]

  expect_equal(round(sum(loglikes), 0), round(sum(li1), 0))
})
