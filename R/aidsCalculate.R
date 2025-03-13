#' Calculate the demand share equations of a AI or QUAI demand system, including demographic variable.
#'
#' @param Prices A matrix of logged prices with (nxm) dimensions where n is the number of observations and m the number of shares.
#' @param Budget A matrix of logged total expenditure/budget with (nx1) dimensions where n is the number of observations.
#' @param ShareNames A vector of strings containing the share names with (mx1) dimensions where m is the number of shares.
#' @param Demographics A matrix of demographic variables with (nxt) dimensions where n is the number of observations and t the number of demographic variables.
#' @param DemographicNames A vector of strings containing the demographic names with (tx1) dimensions where t is the number of demographic variables.
#' @param Params A vector containing the parameters alpha, beta, gamma, and theta and lambda if elected.
#' @param quaids Logical. Should quadratic form be used instead?
#'
#' @return A matrix of estimated shares with (nxm) dimensions where n is the number of observations and m the number of shares.
#' @export
#'
#' @examples
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
#' li1 <- censoredaidsLoglike(
#'   Params = b0,
#'   Shares = matrix(c(s1, s2, s3, s4, s5, s6), ncol = 6),
#'   Prices = matrix(c(lnp1, lnp2, lnp3, lnp4, lnp5, lnp6), ncol = 6),
#'   Budget = matrix(testing_data$lnw),
#'   Demographics = matrix(c(age, size, educ, sex), ncol = 4),
#'   quaids = FALSE
#' )
#'
#'
#'
#'
#'
#'
#'
aidsCalculate <- function(
    Prices = matrix(),
    Budget = matrix(),
    ShareNames = NULL,
    Demographics = matrix(),
    DemographicNames = NULL,
    Params = matrix(),
    quaids = FALSE) {

  # ----: Checks before running the function :----
  # 1. Check correct name dimensions
  if(!is.null(ShareNames)) stopifnot({length(ShareNames) == ncol(Prices)})
  if(!is.null(DemographicNames)) stopifnot({length(DemographicNames) == ncol(Demographics)})

  # 2. Check all matrices have same n
  stopifnot({nrow(Prices) == nrow(Budget)})
  if(all(dim(Demographics) == c(1, 1))) stopifnot({nrow(Budget) == nrow(Demographics)})

  # 3. Check correct column dimensions
  #if(all(dim(Demographics) == c(1, 1))) stopifnot({ncol(Prices) == ncol(Demographics)})

  # 4. Check that implied number of parameters is equal to length(Params)
  m <- ncol(Prices)       # Number of shares dennoted as m
  n <- nrow(Prices)       # Number of observations dennoted as n
  t <- ncol(Demographics) # Number of demographic variables

  nalpha <- m - 1
  nbeta  <- m - 1
  ngamma <- 0.5*(m - 1)*m

  if (!all(dim(Demographics) == c(1, 1))) {ntheta <- (m - 1)*t}else{ntheta <- 0}
  if (quaids) {nlamda <- (m - 1)}else{nlamda <- 0}
  stopifnot({length(Params) == (nalpha + nbeta + ngamma + ntheta + nlamda)})
#browser()
  # ----: Running the function :----
  # Grab the parameters
  Alpha <- Params[1:(m-1)]
  full_alpha <- matrix(c(Alpha,  1 - sum(Alpha)), nrow = m)

  Beta <- Params[m:(2 * (m-1))]
  full_beta <- matrix(c(Beta, - sum(Beta)), nrow = m)

  Gamma <- Params[(2 * (m-1) + 1):((m-1) * (2 + .5 * m))]
  full_gamma <- matrix(0, ncol = (m - 1), nrow = (m - 1))
  full_gamma[upper.tri(full_gamma, diag = TRUE)] <- Gamma
  full_gamma[lower.tri(full_gamma)] <- t(full_gamma)[lower.tri(full_gamma)]
  full_gamma <- cbind(full_gamma, -rowSums(full_gamma))
  full_gamma <- rbind(full_gamma, -colSums(full_gamma))

  nn1 <- (m-1) * (2 + .5 * m) # Helper to make sure we start counting after full_gamma
  if (!all(dim(Demographics) == c(1, 1))) {

      Theta <- Params[(nn1 + 1):(nn1 + t * (m - 1))]
      full_theta <- matrix(0, ncol = m, nrow = t)
      full_theta[1:t, 1:(m - 1)] <- Theta
      full_theta[ , m] <- -rowSums(full_theta)
      full_theta <- t(full_theta)
  }

  if (quaids) {
    #browser()
      if (!all(dim(Demographics) == c(1, 1))) {

          nn2 <- (nn1 + t * (m - 1)) # Helper to make sure we start counting (if) after Demographics
          Lambda <- Params[(nn2 + 1):(nn2 + m - 1)]
          full_lambda <- matrix(c(Lambda, - sum(Lambda)), nrow = m)
      }else {
         nn2 <- nn1
         Lambda <- Params[(nn2 + 1):(nn2 + m - 1)]
         full_lambda <- matrix(c(Lambda, - sum(Lambda)), nrow = m)
      }
  }
#browser()
  # Functions to expand rvector into matrix for summations
  expand_rVector <- function(v, M) matrix(v[col(M)], nrow = nrow(M), ncol = ncol(M))
  expand_cVector <- function(v) matrix(v, nrow = n, ncol = m)

  # Renaming
  Lnp <- Prices
  lnpindex <- Lnp %*% full_alpha + 0.5 * colSums(t(Lnp * t(full_gamma %*% t(Lnp)))) # want this (nx1)
  Lnw <- Budget

  # Creating the shares
  qshare <- expand_rVector(t(full_alpha), Lnp) + Lnp %*% full_gamma

  if (!all(dim(Demographics) == c(1, 1))) {

      Z <- Demographics
      qshare <- qshare + (expand_rVector(full_beta, Lnp) + Z %*% t(full_theta) ) * expand_cVector(Lnw-lnpindex)
  }else {
      qshare <- qshare + (expand_rVector(full_beta, Lnp)) * expand_cVector(Lnw-lnpindex)
  }

  if (quaids) {

      bofp <- exp(Lnp %*% full_beta)  # want this (nx1)
      quaids_term <- (expand_rVector(t(full_lambda), Lnp)/expand_cVector(bofp))*(expand_cVector(Lnw-lnpindex)^2)
      qshare <- qshare + quaids_term #quaids
  }

  if(!is.null(ShareNames)) colnames(qshare) <- ShareNames
  return(qshare)

}
