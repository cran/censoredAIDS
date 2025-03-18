#' Numerical approximation to censored AI or QUAI demand system elasticities, including demographic variable.
#'
#' @param Prices A matrix of logged prices with (nxm) dimensions where n is the number of observations and m the number of shares.
#' @param Budget A matrix of logged total expenditure/budget with (nx1) dimensions where n is the number of observations.
#' @param ShareNames A vector of strings containing the share names with (mx1) dimensions where m is the number of shares.
#' @param Demographics A matrix of demographic variables with (nxt) dimensions where n is the number of observations and t the number of demographic variables.
#' @param DemographicNames A vector of strings containing the demographic names with (tx1) dimensions where t is the number of demographic variables.
#' @param Params A vector containing the parameters alpha, beta, gamma, and theta and lambda if elected.
#' @param quaids Logical. Should quadratic form be used instead?
#' @param vcov A variance-covariance matrix of the parameters. Must be positive semi-definite and symmetric.
#' @param func A function to be applied to the data at which point estimate elasticities are being evaluated.
#' @param ... Additional arguments to be passed to func.
#'
#' @return A list containing a matrix of price and income elasticities, a matrix with their respective standard errors approximated using the delta method, and a matrix of expected shares that can serve to run further post-estimation analyses.
#' @export
#' @examples
#'
#' \dontrun{
#'
#' testing_data <- censoredAIDS::MexicanHH_foodConsumption
#'
#' # Organizing the data for comfort
#' s1 <- testing_data$s1
#' s2 <- testing_data$s2
#' s3 <- testing_data$s3
#' s4 <- testing_data$s4
#' s5 <- testing_data$s5
#' s6 <- testing_data$s6
#'
#' lnp1 <- testing_data$lnp1
#' lnp2 <- testing_data$lnp2
#' lnp3 <- testing_data$lnp3
#' lnp4 <- testing_data$lnp4
#' lnp5 <- testing_data$lnp5
#' lnp6 <- testing_data$lnp6
#'
#' age <- testing_data$age
#' size <- testing_data$size
#' sex <- testing_data$sex
#' educ <- testing_data$educ
#'
#' # Alpha
#' b0 <- rep(0, 5)
#'
#' # Beta
#' b0 <- c(b0, rep(0.003, 5))
#'
#' # Gamma
#' b0 <- c(b0,0.01,0,0.01,0,0, 0.01,0,0,0,0.01,0,0,0,0,0.01)
#'
#' # Demos
#' b0 <- c(b0,rep(0.002, 20))
#'
#' # Sigma
#' b0 <- c(b0,1,0,1,0,0,1,0,0,0,1,0,0,0,0,1)
#'
#' vcov <- matrix(0, nrow = length(b0), ncol = length(b0))
#' vcov[upper.tri(vcov, diag = TRUE)] <- runif(.5*length(b0)*(1+length(b0)))
#' vcov <- t(vcov) %*% vcov
#'
#' list_etas <- censoredElasticity(
#' Params = b0,
#' Shares = matrix(c(s1, s2, s3, s4, s5, s6), ncol = 6),
#' Prices = matrix(c(lnp1, lnp2, lnp3, lnp4, lnp5, lnp6), ncol = 6),
#' Budget = matrix(testing_data$lnw),
#' Demographics = matrix(c(age, size, educ, sex), ncol = 4),
#' quaids = FALSE,
#' func = mean,
#' na.rm = TRUE,
#' vcov = vcov
#' )
#' }
#'
#'
censoredElasticity <- function(
    Prices = matrix(),
    Budget = matrix(),
    ShareNames = NULL,
    Demographics = matrix(),
    DemographicNames = NULL,
    Params = matrix(),
    quaids = FALSE,
    vcov = matrix(),
    func,
    ...) {

  reps <- 1e+5    # Number of simulation replications
  delta <- 1e-5   # Disturbance factor
  demos <- !all(dim(Demographics) == c(1, 1)) # If demos are present

  # ----: Performing checks :----
  # 1. check psd and symmetry of vcov
  stopifnot({matrixcalc::is.positive.semi.definite(vcov) & matrixcalc::is.symmetric.matrix(vcov)})

  # 2. Check if dim(vcov)[1] == length(b0)
  stopifnot({dim(vcov)[1] == length(Params)})

  # 3. Check that implied number of parameters is equal to length(Params)
  m <- ncol(Prices)       # Number of shares dennoted as m
  n <- nrow(Prices)       # Number of observations dennoted as n
  t <- ncol(Demographics) # Number of demographic variables
  j <- 0.5*(m - 1)*m      # Number of sigma dimensions

  nalpha <- m - 1
  nbeta  <- m - 1
  ngamma <- 0.5*(m - 1)*m

  if (demos) {ntheta <- (m - 1)*t}else{ntheta <- 0}
  if (quaids) {nlamda <- (m - 1)}else{nlamda <- 0}
  stopifnot({length(Params) == (nalpha + nbeta + ngamma + ntheta + nlamda + j)})

  # ----: Get me the point estimates :----
  muPrices <- apply(Prices, 2, func, ...)
  muBudget <- apply(Budget, 2, func, ...)
  muDemogs <- apply(Demographics, 2, func, ...)

  # ----: Creating Sigma & error simulations :----
  Sigmma <- matrix(0, ncol = (m - 1), nrow = (m - 1))
  Sigmma[upper.tri(Sigmma, diag = TRUE)] <- Params[c((length(Params) - j + 1):length(Params))]
  Sigmma <- t(Sigmma) %*% Sigmma

  epsilons <- mvtnorm::rmvnorm(n = reps, sigma = Sigmma, method = "chol")
  epsilons <- cbind(epsilons, -rowSums(epsilons))

  # ----: Amemiya-Tobin mapping of latent shares: f: R->[0, 1] :----
  # Wales and Woodland (1983) and Amemiya and Tobin (1987)
  treatTruncations <- function(x) {

    # Function corrects truncations and guarantees sum of shares equal 1

    x[x < 0] <- 0 # Correcting the left truncations
    rsums <- rowSums(x) # The denominator
    apply(x, 2, function(L) L/rsums) # Correcting the right truncations

  }

  # ----: E(X,b) -- expected share without disturbance :----
  # Constant part no disturbances
  U <- muaidsCalculate(muPrices = muPrices,
                       muBudget = muBudget,
                       muDemographics = muDemogs,
                       Params = Params[-c((length(Params) - j + 1):length(Params))],
                       quaids = quaids,
                       m = m,
                       t = t,
                       dems = demos)

  # Add the disturbances, treat truncations, and obtain the shares
  Ulat <- t(apply(epsilons, 1, function(x) U + x))
  Uobs <- treatTruncations(Ulat)
  E_Uobs <- apply(Uobs, 2, mean)

  # ----: E(X + delta, b) -- expected share with disturbance in variable :----
  # Constant part with disturbance in variable
  U_dx <- sapply(1:(1 + length(muPrices)), function(x) {
    # Note: Creates a matrix of (mxv) where v is the number of prices + 1 and m is the number of shares
    if(x > length(muPrices)) {

      muBudget_delta <- exp(muBudget)
      muBudget_delta <- log(muBudget_delta + delta)
      muaidsCalculate(muPrices = muPrices,
                      muBudget = muBudget_delta,
                      muDemographics = muDemogs,
                      Params = Params[-c((length(Params) - j + 1):length(Params))],
                      quaids = quaids,
                      m = m,
                      t = t,
                      dems = demos)
    }else{

      muPrices_delta <- exp(muPrices)
      muPrices_delta[x] <- muPrices_delta[x] + delta
      muPrices_delta <- log(muPrices_delta)
      muaidsCalculate(muPrices = muPrices_delta,
                      muBudget = muBudget,
                      muDemographics = muDemogs,
                      Params = Params[-c((length(Params) - j + 1):length(Params))],
                      quaids = quaids,
                      m = m,
                      t = t,
                      dems = demos)
    }
  })
  U_dx <- t(U_dx) # Transpose for simplicity

  EUobs_dx <- lapply(split(U_dx,1:nrow(U_dx)), function(u) {

    Ulat_dx <- t(apply(epsilons, 1, function(e) u + e)) # Disturbed latent shares
    Uobs_dx <- treatTruncations(Ulat_dx) # Disturbed observed shares
    apply(Uobs_dx, 2, mean) # Expected disturbed shares

  })

  # ----: E(X, b + delta) -- expected share with disturbance in parameter :----
  # Constant part with disturbance in parameter
  U_db <- sapply(1:length(Params[-c((length(Params) - j + 1):length(Params))]), function(x) {
    # Note: Creates a matrix of (mxp) where v is the number of parameters
    b_delta <- Params[-c((length(Params) - j + 1):length(Params))]
    b_delta[x] <- b_delta[x] + delta

    muaidsCalculate(muPrices = muPrices,
                    muBudget = muBudget,
                    muDemographics = muDemogs,
                    Params = b_delta,
                    quaids = quaids,
                    m = m,
                    t = t,
                    dems = demos)

  })
  U_db <- t(U_db) # Transpose for simplicity

  EUobs_db <- lapply(split(U_db,1:nrow(U_db)), function(u) {

    Ulat_db <- t(apply(epsilons, 1, function(e) u + e)) # Disturbed latent shares
    Uobs_db <- treatTruncations(Ulat_db) # Disturbed observed shares
    apply(Uobs_db, 2, mean) # Expected disturbed shares

  })
  EUobs_db <- do.call(rbind, EUobs_db)

  # ----: E(X + delta, b + delta) -- expected share with disturbance in variable and parameter :----
  # Constant part with disturbance in variable
  U_dxb <- lapply(1:(1 + length(muPrices)), function(x) {
    # Note: Creates a list of matrices of (vxp) where v is the number of prices + 1 and p is the number of parameters
    mat_out <- matrix(0, nrow = m, ncol = length(Params[-c((length(Params) - j + 1):length(Params))]))
    for(p in 1:length(Params[-c((length(Params) - j + 1):length(Params))])) {

      b_delta <- Params[-c((length(Params) - j + 1):length(Params))]
      b_delta[p] <- b_delta[p] + delta

      if(x > length(muPrices)) {

        muBudget_delta <- exp(muBudget)
        muBudget_delta <- log(muBudget_delta + delta)
        out <- muaidsCalculate(muPrices = muPrices,
                               muBudget = muBudget_delta,
                               muDemographics = muDemogs,
                               Params = b_delta,
                               quaids = quaids,
                               m = m,
                               t = t,
                               dems = demos)
        mat_out[,p] <- out
      }else{

        muPrices_delta <- exp(muPrices)
        muPrices_delta[x] <- muPrices_delta[x] + delta
        muPrices_delta <- log(muPrices_delta)
        out <- muaidsCalculate(muPrices = muPrices_delta,
                               muBudget = muBudget,
                               muDemographics = muDemogs,
                               Params = b_delta,
                               quaids = quaids,
                               m = m,
                               t = t,
                               dems = demos)
        mat_out[,p] <- out
      }
    }
    mat_out
  })

  EUobs_dxb <- lapply(1:(1 + length(muPrices)), function(x) {

    focus_dxb <- t(U_dxb[[x]])
    EUobs_dxb <- lapply(split(focus_dxb,1:nrow(focus_dxb)), function(l) {
      Ulat_dxb <- t(apply(epsilons, 1, function(xx) l + xx)) # Disturbed latent shares
      Uobs_dxb <- treatTruncations(Ulat_dxb) # Disturbed observed shares
      apply(Uobs_dxb, 2, mean) # Expected disturbed shares
    })
    do.call(cbind, EUobs_dxb)
  })

  # ----: Elasticity calculations :----
  # Matrix of elasticities are mx(m+1) where rows are quantity and columns are wrt. variables, starting with budget
  m_EUobs_dx <- do.call(cbind, EUobs_dx)
  etas <- sapply(1:m, function(i) {

    # Income Elasticity
    p1 <- (E_Uobs[i] - m_EUobs_dx[i,(m + 1)])/delta
    p2_num <- exp(muBudget) + .5*delta
    p2_den <- E_Uobs[i] + .5*(E_Uobs[i] - m_EUobs_dx[i,(m + 1)])
    p2 <- p2_num/p2_den
    income_eta <- p1*p2 + 1

    # Price Elasticity
    price_eta <- sapply(1:m, function(j) {

      p1 <- (E_Uobs[i] - m_EUobs_dx[i,j])/delta
      p2_num <- exp(muPrices[j]) + .5*delta
      p2_den <- E_Uobs[i] + .5*(E_Uobs[i] - m_EUobs_dx[i,j])
      p2 <- p2_num/p2_den
      p1*p2

    })

    return(c(price_eta, income_eta))

  })

  # ----: Formatting first output: Price and Income elasticities :----
  etas[1:m, 1:m] <- etas[1:m, 1:m] - diag(1, m) # Correcting the diagonal
  etas <- t(etas)
  colnames(etas) <- c(paste0("Price ", as.character(1:m)), "Income")
  rownames(etas) <- paste0("Quantity", as.character(1:m))

  if(!is.null(ShareNames)) {
    colnames(etas)[1:m] <- ShareNames
    rownames(etas) <- ShareNames
  }

  # ----: Standard Errors :----
  # Matrix of elasticities are mx(m+1) where rows are quantity and columns are wrt. variables, starting with budget
  etas_SE <- lapply(1:m, function(i) {

    # Income Elasticity
    m_EUobs_dxb <- t(EUobs_dxb[[(m + 1)]])
    p1 <- (EUobs_db[,i] - m_EUobs_dxb[,i])/delta
    p2_num <- exp(muBudget) + .5*delta
    p2_den <- EUobs_db[,i] + .5*(EUobs_db[,i] - m_EUobs_dxb[,i])
    p2 <- p2_num/p2_den
    SE_income_eta <- p1*p2 + 1

    # Price Elasticity
    SE_price_eta <- sapply(1:m, function(j) {
      m_EUobs_dxb <- t(EUobs_dxb[[j]])
      p1 <- (EUobs_db[,i] - m_EUobs_dxb[,i])/delta
      p2_num <- exp(muPrices[j]) + .5*delta
      p2_den <- EUobs_db[,i] + .5*(EUobs_db[,i] - m_EUobs_dxb[,i])/delta
      p2 <- p2_num/p2_den
      if(i == j) {p1*p2 - 1} else {p1*p2}
    })

    return(cbind(SE_price_eta, SE_income_eta))

  })

  # Helper function to expand rvector into matrix for summations
  expand_rVector <- function(v, M) matrix(v[col(M)], nrow = nrow(M), ncol = ncol(M))

  # Calculation of SEs
  etas_SE <- lapply(1:length(etas_SE), function(x) {

    dy <- etas_SE[[x]]
    f <- expand_rVector(etas[x,], etas_SE[[x]])
    J <- (dy - f) / delta # Jacobian approximation

    # Var(eta(b)) by delta method
    vcov_dims <- length(Params[-c((length(Params) - j + 1):length(Params))])
    etavcov <- t(J) %*% (vcov[1:vcov_dims,1:vcov_dims]/n) %*% J
    sqrt(diag(etavcov)) # Retrieve only the diagonal for S.E.

  })

  # ----: Formatting second output: Elasticity SE:----
  etas_SE <- do.call(rbind, etas_SE)
  colnames(etas_SE) <- c(paste0("Price ", as.character(1:m)), "Income")
  rownames(etas_SE) <- paste0("Quantity", as.character(1:m))

  if(!is.null(ShareNames)) {
    colnames(etas_SE)[1:m] <- ShareNames
    rownames(etas_SE) <- ShareNames
  }

  # ----: Outputs :----
  return(list(Elasticities = etas,
              SE = etas_SE,
              E_Uobs = E_Uobs))

}
