test_that("aidsCalculate gives nxm matrix with names, and replicates q shares from JAAEA study", {

  testthat::skip_on_cran()

  # Reference: https://github.com/NoeNava-USDA/mexSugarTax_replication
  # Nava, Noé J., and Diansheng Dong. 2022. "The impact of taxing sugary-sweetened beverages in México: A censored QUAI demand system approach." Journal of Agricultural and Applied Economics Association, 1(1):1-23:
  # https://onlinelibrary.wiley.com/doi/10.1002/jaa2.6

  # qshares <- readRDS(system.file(testthat::test_path("qshares.rds"), package = "censoredAIDS"))
  # testing_data <- readRDS(system.file(testthat::test_path('testing_data.rds'), package = "censoredAIDS"))
  #
  # full_alpha <- readRDS(system.file(testthat::test_path('testing_params/full_alpha.rds'), package = "censoredAIDS"))
  # full_beta <- readRDS(system.file(testthat::test_path('testing_params/full_beta.rds'), package = "censoredAIDS"))
  # full_gamma <- readRDS(system.file(testthat::test_path('testing_params/full_gamma.rds'), package = "censoredAIDS"))
  # full_lamda <- readRDS(system.file(testthat::test_path('testing_params/full_lamda.rds'), package = "censoredAIDS"))
  # full_theta <- readRDS(system.file(testthat::test_path('testing_params/full_theta.rds'), package = "censoredAIDS"))
  # full_sigma <- readRDS(system.file(testthat::test_path('testing_params/full_sigma.rds'), package = "censoredAIDS"))


  qshares <- readRDS(testthat::test_path('qshares.rds'))
  testing_data <- readRDS(testthat::test_path('testing_data.rds'))

  full_alpha <- readRDS(testthat::test_path('testing_params/full_alpha.rds'))
  full_beta <- readRDS(testthat::test_path('testing_params/full_beta.rds'))
  full_gamma <- readRDS(testthat::test_path('testing_params/full_gamma.rds'))
  full_lamda <- readRDS(testthat::test_path('testing_params/full_lamda.rds'))
  full_theta <- readRDS(testthat::test_path('testing_params/full_theta.rds'))
  full_sigma <- readRDS(testthat::test_path('testing_params/full_sigma.rds'))
  full_sigma <- chol(full_sigma)

  qqshares <- aidsCalculate(
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
               c()),
    quaids = TRUE
  )
  colnames(qshares) <- c("SSB", "Juice", "Milk", "Water")
  expect_equal(qshares, qqshares)

})
