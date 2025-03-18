#' Calculate the demand share equations of a AI or QUAI demand system, including demographic variable, using point estimates (e.g., mean, median, percentiles).
#'
#' @param muPrices A matrix of logged muPrices with (1xm) dimensions where n is the number of observations and m the number of shares.
#' @param muBudget A matrix of logged total expenditure/muBudget with (1x1) dimensions where n is the number of observations.
#' @param muDemographics A matrix of demographic variables with (1xt) dimensions where n is the number of observations and t the number of demographic variables.
#' @param Params A vector containing the parameters alpha, beta, gamma, and theta and lambda if elected.
#' @param quaids Logical. Should quadratic form be used instead?
#' @param m Number of shares.
#' @param t Number of demographic variables.
#' @param dems Bolean. Should demographic variables be included?
#'
#' @return A matrix of estimated shares with (1xm) dimensions where n is the number of observations and m the number of shares.
#'
#'
muaidsCalculate <- function(
    muPrices = numeric(),
    muBudget = numeric(),
    muDemographics = numeric(),
    Params = matrix(),
    quaids = FALSE,
    m = numeric(),
    t = numeric(),
    dems = FALSE) {

  # -----: Make the demographics correc dimension matrices :-----
  muPrices <- matrix(muPrices, nrow = 1, ncol = m)
  muBudget <- matrix(muBudget, nrow = 1, ncol = 1)
  if (dems) {
    muDemographics <- matrix(muDemographics, nrow = 1, ncol = t)
  }

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
  if (dems) {

    Theta <- Params[(nn1 + 1):(nn1 + t * (m - 1))]
    full_theta <- matrix(0, ncol = m, nrow = t)
    full_theta[1:t, 1:(m - 1)] <- Theta
    full_theta[ , m] <- -rowSums(full_theta)
    full_theta <- t(full_theta)
  }

  if (quaids) {
    #browser()
    if (dems) {

      nn2 <- (nn1 + t * (m - 1)) # Helper to make sure we start counting (if) after muDemographics
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
  expand_cVector <- function(v) matrix(v, nrow = 1, ncol = m)

  # Renaming
  Lnp <- muPrices
  lnpindex <- Lnp %*% full_alpha + 0.5 * colSums(t(Lnp * t(full_gamma %*% t(Lnp)))) # want this (1x1)
  Lnw <- muBudget

  # Creating the shares
  qshare <- expand_rVector(t(full_alpha), Lnp) + Lnp %*% full_gamma

  if (dems) {

    Z <- muDemographics
    qshare <- qshare + (expand_rVector(full_beta, Lnp) + Z %*% t(full_theta) ) * expand_cVector(Lnw-lnpindex)
  }else {
    qshare <- qshare + (expand_rVector(full_beta, Lnp)) * expand_cVector(Lnw-lnpindex)
  }

  if (quaids) {

    bofp <- exp(Lnp %*% full_beta)  # want this (1x1)
    quaids_term <- (expand_rVector(t(full_lambda), Lnp)/expand_cVector(bofp))*(expand_cVector(Lnw-lnpindex)^2)
    qshare <- qshare + quaids_term #quaids
  }

  return(qshare)

}
