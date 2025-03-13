#' Calculates the truncated (censored) log-likehood function of a share system of equations of a AI or QUAI demand system.
#'
#' @param Shares A matrix of shares with (nxm) dimensions where n is the number of observations and m the number of shares.
#' @param Prices A matrix of logged prices with (nxm) dimensions where n is the number of observations and m the number of shares.
#' @param Budget A matrix of logged total expenditure/budget with (nx1) dimensions where n is the number of observations.
#' @param ShareNames A vector of strings containing the share names with (mx1) dimensions where m is the number of shares.
#' @param Demographics A matrix of demographic variables with (nxt) dimensions where n is the number of observations and t the number of demographic variables.
#' @param DemographicNames A vector of strings containing the demographic names with (tx1) dimensions where t is the number of demographic variables.
#' @param Params A vector containing the parameters alpha, beta, gamma, and theta and lambda if elected.
#' @param quaids Logical. Should quadratic form be used instead?
#'
#' @return A numeric vector for individual loglikehood contributions with dimensions (nx1) where n is the number of observations.
#' @export
#'
#' @examples
#'

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
#'   f <- censoredaidsLoglike(
#' Params = b0,
#' Shares = matrix(c(s1, s2, s3, s4, s5, s6), ncol = 6),
#' Prices = matrix(c(lnp1, lnp2, lnp3, lnp4, lnp5, lnp6), ncol = 6),
#' Budget = matrix(testing_data$lnw),
#' Demographics = matrix(c(age, size, educ, sex), ncol = 4),
#' quaids = FALSE
#' )
#'
#'
#'
#'
censoredaidsLoglike <- function(
    Shares = matrix(),
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
  if(!all(dim(Demographics) == c(1, 1))) stopifnot({nrow(Budget) == nrow(Demographics)})

  # 3. Check correct column dimensions
  #if(!all(dim(Demographics) == c(1, 1))) stopifnot({ncol(Prices) == ncol(Demographics)})

  # 4. Check that implied number of parameters is equal to length(Params)
  m <- ncol(Prices)       # Number of shares dennoted as m
  n <- nrow(Prices)       # Number of observations dennoted as n
  t <- ncol(Demographics) # Number of demographic variables
  j <- 0.5*(m - 1)*m      # Number of sigma dimensions

  nalpha <- m - 1
  nbeta  <- m - 1
  ngamma <- 0.5*(m - 1)*m

  if (!all(dim(Demographics) == c(1, 1))) {ntheta <- (m - 1)*t}else{ntheta <- 0}
  if (quaids) {nlamda <- (m - 1)}else{nlamda <- 0}
  stopifnot({length(Params) == (nalpha + nbeta + ngamma + ntheta + nlamda + j)})

  # ----: Calculate shares :----
  U <- aidsCalculate(Prices = Prices,
                     Budget = Budget,
                     Demographics = Demographics,
                     Params = Params[-c((length(Params) - j + 1):length(Params))],
                     quaids = quaids,
                     ShareNames = ShareNames)
  S <- Shares # Renaming shares

  # ----: Truncated log-likelihood function :----

  # Dummy for micro-regimes and index_m re-arragement
  d <- ifelse(S, 1, 0)
  nu_bght <- rowSums(d) # Vector with number of purchased goods per household

  # Log-likelihood values repository
  lf <- matrix(rep(0, n), ncol = 1, nrow = n)

  #### Regime where all goods are purchased ###
  index_m <- nu_bght == m # Index for all households that buy all goods

  # Epsilons = Obs. share - param share (depicted in Eq. 18)
  e <- S[index_m, 1:(m - 1)] - U[index_m, 1:(m - 1)]

  Sigmma <- matrix(0, ncol = (m - 1), nrow = (m - 1))
  Sigmma[upper.tri(Sigmma, diag = TRUE)] <- Params[c((length(Params) - j + 1):length(Params))]
  Sigmma <- t(Sigmma) %*% Sigmma

  lf[index_m] <- mvtnorm::dmvnorm(e, sigma = Sigmma, log = TRUE)

  full_sigma <- cbind(Sigmma, -rowSums(Sigmma))
  full_sigma <- rbind(full_sigma, -colSums(full_sigma))

  for(i in 1:n) {

    bght_index <- which(d[i,] == 1)         # indexing if good is bought
    zero_index <- which(d[i,] == 0)         # indexing if good is not bought
    index_arrn <- c(bght_index, zero_index) # index for re-arragement

    if(nu_bght[i] < m) {
      # First, we re-arrange the our vectors and matrices
      sigma <- full_sigma[index_arrn, index_arrn]   # sorted full sigma
      sigma <- sigma[1:(m - 1), 1:(m - 1)]          # sorted small sigma

      S_a   <- S[i, index_arrn]              # sorted shares
      U_a   <- U[i, index_arrn]              # sorted predicted shares
      U_a   <- matrix(U_a, ncol = 1, nrow = m)

      a <- S_a[1:nu_bght[i]]/S_a[1]    # vector
      a <- matrix(a, ncol = 1, nrow = length(a))

      if(nu_bght[i] == (m - 1)) { ### Regime when all goods except one are bought ###

        AA1   <- diag(array(a))
        omega <- AA1 %*% mnormt::pd.solve(sigma) %*% t(AA1)
        s11   <- omega[1:nu_bght[i], 1:nu_bght[i]]
        s11   <- matrix(s11, nrow = length(1:nu_bght[i]), ncol = length(1:nu_bght[i]))

        U_bar <- U_a[1]
        II    <- matrix(rep(1,nu_bght[i]), ncol = 1, nrow = (m - 1))
        JJ    <- U_a[1:nu_bght[i]]/(a * U_a[1])
        JJ    <- matrix(JJ, ncol = 1, nrow = (m - 1))

        omega11 <- t(II) %*% s11 %*% II
        omega10 <- t(II) %*% s11 %*% JJ
        omega00 <- t(JJ) %*% s11 %*% JJ

        inv_omega11 <- mnormt::pd.solve(omega11)
        U_sta     <- (inv_omega11 %*% omega10)*U_bar
        omg_final <- inv_omega11

        part1 <- t(U_bar) %*% omega00 %*% U_bar - t(U_sta) %*% omega11 %*% U_sta

        part2 <- exp(-0.5*part1)*((2*pi)^(0.5*(1 - nu_bght[i])))*(((det(sigma))^(-0.5))/((det(omg_final))^(-0.5)))

        D_2 <- diag((diag(omg_final)^0.5), m - nu_bght[i])
        RR <- stats::cov2cor(omg_final)

        AA  <- 1/S_a[1]
        CC  <- AA %*% D_2
        PP0 <- 1
        #AA[1] <- -AA[1]
        PP <- PP0 - AA*U_sta
        diagD <- rep(0, m - nu_bght[i])
        for(j in 1:(m - nu_bght[i])) {
          diagD[j] <- (CC[j] %*% RR %*% t(CC[j]))^-0.5
        }
        diagD <- as.data.frame(diagD)
        DD <- diag(diagD, ncol = m - nu_bght[i], nrow = m - nu_bght[i])
        BB <- DD %*% PP
        R_c <- DD %*% CC %*% RR %*% t(CC) %*% DD

        # Depositing log-likelihood values in their positions
        HAJP <- log(
          part2 + 1e-7
        ) + log(
          mvtnorm::pmvnorm(
            upper = as.numeric(-BB),
            sigma = diag(1))[1] + 1e-7
        )
        lf[i] <- HAJP

      }else{ ### All other regimes ###

        ones <- matrix(rep(1, m - nu_bght[i] - 1), nrow = m - nu_bght[i] - 1, ncol = 1)
        AA1 <- rbind(a, ones)
        AA1 <- diag(AA1[1:length(AA1),])

        omega <- AA1 %*% mnormt::pd.solve(sigma) %*% t(AA1)

        s11 <- omega[1:nu_bght[i],1:nu_bght[i]]
        s11 <- matrix(s11,nrow = nu_bght[i], ncol = nu_bght[i])

        s10 <- omega[1:nu_bght[i],(nu_bght[i] +1):(m - 1)]
        s10 <- matrix(s10, nrow = nu_bght[i],
                      ncol = length((nu_bght[i] +1):(m - 1)))

        s00 <- omega[(nu_bght[i]+1):(m - 1),(nu_bght[i]+1):(m - 1)]
        s00 <- matrix(s00, nrow = length((nu_bght[i] +1):(m - 1)),
                      ncol = length((nu_bght[i] +1):(m - 1)))

        Ubar <- U_a[c(1, (nu_bght[i]+1):(m - 1))]
        Ubar <- matrix(Ubar, nrow = length(c(1, (nu_bght[i]+1):(m - 1))), ncol = 1)
        II   <- matrix(rep(1,nu_bght[i]), ncol = 1, nrow = nu_bght[i])
        JJ   <- U_a[1:nu_bght[i]] / (a %*% U_a[1])

        omega11 <- rbind(cbind(t(II) %*% s11 %*% II, t(II) %*% s10), cbind(t(s10) %*% II, s00))
        omega10 <- rbind(cbind(t(II) %*% s11 %*% JJ, t(II) %*% s10), cbind(t(s10) %*% JJ, s00))
        omega00 <- rbind(cbind(t(JJ) %*% s11 %*% JJ, t(JJ) %*% s10), cbind(t(s10) %*% JJ, s00))

        inv_omega11 <- mnormt::pd.solve(omega11)
        U_sta <- inv_omega11 %*% omega10 %*% Ubar
        omg_final <- inv_omega11

        part1 <- t(Ubar) %*% omega00 %*% Ubar - t(U_sta) %*% omega11 %*% U_sta

        part2 <- exp(-0.5*part1)*((2*pi)^(0.5*(1 - nu_bght[i])))*(((det(sigma))^(-0.5))/((det(omg_final))^(-0.5)))

        D_2 <- diag((diag(omg_final)^0.5), nrow = m - nu_bght[i], ncol = m - nu_bght[i])

        RR <- stats::cov2cor(omg_final)

        AA     <- diag(-1*rep(1, m - nu_bght[i]), ncol = m - nu_bght[i], nrow = m - nu_bght[i])
        ones   <- matrix(rep(1, m - nu_bght[i] - 1), ncol = m - nu_bght[i] - 1, nrow = 1)

        AA[1,] <- cbind((1/S_a[1]), ones)
        CC     <- AA %*% D_2
        zeros  <- matrix(rep(0, m - nu_bght[i] - 1), ncol = 1, nrow = m - nu_bght[i] - 1)
        PP0    <- rbind(1, zeros)
        PP     <- PP0 - (AA %*% U_sta)

        diagD <- rep(0, m - nu_bght[i])
        diagD <- matrix(diagD, nrow = 1, ncol = m - nu_bght[i])
        for(j in 1:(m - nu_bght[i])) {
          row <- as.matrix(CC[j,])
          row <- t(row)
          diagD[j] <- (row %*% RR %*% t(row))^-0.5
        }

        DD  <- diag(array(diagD), ncol = length(diagD), nrow = length(diagD))
        BB  <- DD %*% PP
        R_c <- DD %*% CC %*% RR %*% t(CC) %*% DD
        R_c <- as.matrix(Matrix::nearPD(R_c)$mat)
#if(i == 12) browser()
          HAJP <- log(
            part2 + 1e-8
          ) + log(
            mvtnorm::pmvnorm(
              upper = as.numeric(-BB),
              lower = rep(-Inf, nrow(BB)),
              sigma = R_c)[1] + 1e-8)
        }

        lf[i] <- HAJP
      }
    }

  return(lf)

}
