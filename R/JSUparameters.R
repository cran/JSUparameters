JSUparameters <- function(dat) {

  # =====================================================================================
  # MAIN FUNCTION
  # NULL variables to avoid warnings
  sorted <- probs <- theoretical_quantiles <- flag1 <- flag2 <- zeta <- X <- W <- NULL
  # sort the data
  sorted <- sort(dat)
  # theoretical quantiles
  probs <- (1:length(sorted) - 0.5)/length(sorted)
  theoretical_quantiles <- qnorm(probs)

  # set flags manually at the beginning
  flag1 <- 0
  flag2 <- 0   # this may get set in golden but need to just store it for the optimal
  # value of zeta so we will reset it later

  # =====================================================================================
  # CHECK 4 CASES FUNCTION
  Check4Cases = function(sorted, theoretical_quantiles, delta) {
    # NULL variables to avoid warnings
    X <- W <- flag2 <- NULL
    # create the matrix A for the constraints
    A = matrix(c(0, 0, 1, 1, -delta, delta), byrow = F, nrow = 2)

    # JSU distribution -------------------------------------------------
    # calculate the unconstrained minimum
    X <- sorted
    # create the matrix W
    w1 = rep(1, length(sorted))
    w2 = c(delta * (sinh(theoretical_quantiles/delta)))
    w3 = c((delta^2) * (cosh(theoretical_quantiles/delta) - exp(1/(2 * (delta^2)))))
    W <- matrix(c(w1, w2, w3), nrow = length(sorted), ncol = 3, byrow = F)
    # calculate the unconstrained minimum
    inverse = solve(t(W)%*%W, tol = 1e-170)
    betas = inverse %*% t(W) %*% X
    # calculate the SSQ
    ssq = t(X)%*%X - 2*t(betas)%*%t(W)%*%X + t(betas)%*%t(W)%*%W%*%betas

    # set alpha = 0
    alpha = c(0, 0)

    # check conditions with unconstrained minimum
    if (all((alpha >= 0)) & all((A%*%betas >= 0)) & (t(alpha)%*%A%*%betas == 0)) {
      # unconstrained optimum is also a constrained optimum
      CaseID = "JSU distribution"
      # calculate JSU parameters
      gamma = delta * atanh(-((delta * betas[3])/betas[2]))
      lambda = (delta * betas[2])/cosh(gamma/delta)
      if ((betas[2] == 0) | ((abs(lambda) < 10e-13) & all(((1/(2*delta)) * (matrix(c(0, delta, delta, 0, -1, 1), nrow = 2, byrow = T))%*%(t(W)%*%W%*%(c(mean(sorted), 0, 0)) - t(W)%*%X)) <= 10e-13))) {
        CaseID = "Constant distribution"
        constant = mean(sorted)
        SSQ = sum((sorted - constant)^2)
        results = list(CaseID = CaseID, SSQ = SSQ, constant = constant, flags = c(flag1, flag2))
        class(results) = "Constant"
        return(results)
      }
      else {
        xi = -(-betas[1] - ((lambda)*(sinh(gamma/delta))*(exp(1/(2 * (delta^2))))))
        # calculate SSQ
        SSQ = ssq
        # list of JSU parameter estimates
        results = list(CaseID = CaseID, SSQ = SSQ, delta = delta, xi = xi, gamma = gamma, lambda = lambda, flags = c(flag1, flag2))
        class(results) = "JSU"
        return(results)
      }
    }
    else {
      # Lognormal distribution -----------------------------------------
      # define unit vector
      e1 = c(1, 0)
      # solve the first order condition, with alpha being proportional to e1
      # and e1^TAbeta = 0
      # we calculate the analytical solution (which may not be feasible)
      alpha = -((t(e1)%*%A%*%inverse%*%t(W)%*%sorted%*%e1)/(t(e1)%*%A%*%inverse%*%t(A)%*%e1)[1, 1])
      alpha = as.vector(alpha)
      # calculate the corresponding beta_star
      betas = inverse%*%t(W)%*%sorted + inverse%*%t(A)%*%alpha

      # check constraints
      if (all((alpha >= 0)) & ((as.vector(A%*%betas)[2] >= 10e-13)) & (t(alpha)%*%c(0, as.vector(A%*%betas)[2]) == 0)) {
        # constrained optimum is the shifted lognormal distribution
        CaseID = "Shifted lognormal distribution"
        # calculate SSQ
        SSQ = t(X)%*%X - 2*t(betas)%*%t(W)%*%X + t(betas)%*%t(W)%*%W%*%betas
        # transform to lognormal parameters
        log_mod = lm(sorted ~ (exp(theoretical_quantiles/delta)))
        xi = log_mod$coefficients[1][[1]]
        lambda = log_mod$coefficients[2][[1]]
        # check if lambda = 0 (then should be a constant case)
        if ((abs(lambda) < 10e-13) & all(((1/(2*delta)) * (matrix(c(0, delta, delta, 0, -1, 1), nrow = 2, byrow = T))%*%(t(W)%*%W%*%(c(mean(sorted), 0, 0)) - t(W)%*%X)) >= 0)) {
          CaseID = "Constant distribution"
          constant = mean(sorted)
          SSQ = sum((sorted - constant)^2)
          results = list(CaseID = CaseID, SSQ = SSQ, constant = constant, flags = c(flag1, flag2))
          class(results) = "Constant"
          return(results)
        }
        else {
          # check if any data points are beyond this
          if (any(sorted < xi)) {
            flag2 <- 1
          } else {
            flag2 <- 0
          }
          results = list(CaseID = CaseID, SSQ = SSQ, delta = delta, xi = xi, lambda = lambda, flags = c(flag1, flag2))
          class(results) = "Lognormal"
          return(results)
        }
      }
      else {
        # Negative lognormal distribution -----------------------------------------
        # define unit vector
        e2 = c(0, 1)
        # solve the first order condition, with alpha being proportional to e2
        # and e2^TAbeta = 0
        # we calculate the analytic solution (which may not be feasible)
        alpha = -((t(e2)%*%A%*%inverse%*%t(W)%*%sorted%*%e2)/(t(e2)%*%A%*%inverse%*%t(A)%*%e2)[1, 1])
        alpha = as.vector(alpha)
        # calculate the corresponding beta_star
        betas = inverse%*%t(W)%*%sorted + inverse%*%t(A)%*%alpha

        # check constraints
        if (all((alpha >= 0)) & ((as.vector(A%*%betas)[1] >= 0)) & (t(alpha)%*%c(as.vector(A%*%betas)[2], 0) == 0)) {
          # constrained optimum is the shifted negative lognormal distribution
          CaseID = "Shifted negative lognormal distribution"
          # calculate SSQ
          SSQ = t(X)%*%X - 2*t(betas)%*%t(W)%*%X + t(betas)%*%t(W)%*%W%*%betas
          # transform to negative lognormal parameters
          mod_df = data.frame(sorted = sorted, term = exp(-theoretical_quantiles/delta))
          log_mod = lm(sorted ~ term, data = mod_df)
          xi = log_mod$coefficients[1][[1]]
          lambda = -(log_mod$coefficients[2][[1]])
          if ((abs(lambda) < 10e-13) & all(((1/(2*delta)) * (matrix(c(0, delta, delta, 0, -1, 1), nrow = 2, byrow = T))%*%(t(W)%*%W%*%(c(mean(sorted), 0, 0)) - t(W)%*%X)) <= 10e-13)) {
            CaseID = "Constant distribution"
            constant = mean(sorted)
            SSQ = sum((sorted - constant)^2)
            results = list(CaseID = CaseID, SSQ = SSQ, constant = constant, flags = c(flag1, flag2))
            class(results) = "Constant"
            return(results)
          }
          else {
            # check if any data points are beyond this
            if (any(sorted > xi)) {
              flag2 <- 1
            } else {
              flag2 <- 0
            }
            results = list(CaseID = CaseID, SSQ = SSQ, delta = delta, xi = xi, lambda = lambda, flags = c(flag1, flag2))
            class(results) = "NegLognormal"
            return(results)
          }
        }
        else {
          # Constant distribution --------------------------------------------------

          # need to define alpha in some way
          betas = c(mean(sorted), 0, 0)
          mat = matrix(c(0, delta, delta, 0, -1, 1), nrow = 2, byrow = T)
          alpha = (1/(2*delta)) * mat%*%(t(W)%*%W%*%betas - t(W)%*%X)
          # so only need to check if alpha is nonnegative for this case to hold
          if (all((alpha <= 10e-13))) {
            CaseID = "Constant distribution"
            constant = mean(sorted)
            SSQ = sum((sorted - constant)^2)
            results = list(CaseID = CaseID, SSQ = SSQ, constant = constant,flags = c(flag1, flag2))
            class(results) = "Constant"
            return(results)
          }
          else {
            CaseID = "None"
            results = list(CaseID = CaseID, flags = c(flag1, flag2))
            class(results) = "None"
            return(results)
          }
        }
      }
    }
  }
  # =====================================================================================
  # CALCULATE SSQ FUNCTION
  calculate_ssq = function(zeta) {
    if (zeta < 10e-13) {
      # delta = 0
      smallest = sorted[1]
      largest = sorted[length(sorted)]
      average = mean(sorted[2:(length(sorted) - 1)])
      # calculate the SSQ
      ssq = sum(((sorted[2:(length(sorted) - 1)]) - average)^2)
    }
    else if (zeta > 0.9999) {
      # delta = infinity
      df_mod = data.frame(sorted = sorted, theoretical_quantiles = theoretical_quantiles)
      mod = lm(sorted ~ theoretical_quantiles, data = df_mod)
      # calculate SSQ, intercept and slope
      ssq = sum(mod$residuals^2)
    }
    else { # in between
      delta <- res <- NULL
      delta <- zeta/(1 - zeta)
      res <- Check4Cases(sorted, theoretical_quantiles, delta)
      ssq = res$SSQ
    }
    return(ssq)
  }
  # =====================================================================================
  # GOLDEN SECTION SEARCH FUNCTION
  golden = function(func, lower, upper) {
    # NULL variables to avoid warnings
    flag1 <- NULL
    # set up constants
    ratio = 2/(sqrt(5) + 1)
    tol = 10e-15

    # set up initial test points
    x1 = upper - (ratio*(upper - lower))
    x2 = lower + (ratio*(upper - lower))

    # evaluate the function at the test points
    f1 = func(x1)
    f2 = func(x2)

    # check our assumption of one minimum holds based on first four points
    fupper = func(upper)
    flower = func(lower)

    # flag1 <- NULL
    # 4 scenarios are okay and 4 are not - we will return a flag for those not okay
    if (((flower > f1) & (f1 < f2) & (f2 > fupper)) | ((flower < f1) & (f1 < f2) & (f2 > fupper)) | ((flower < f1) & (f1 > f2) & (f2 > fupper)) | ((flower < f1) & (f1 > f2) & (f2 < fupper))) {
      flag1 <- 1
    } else {
      flag1 <- 0
    }

    # continue with GSS
    while (abs(upper - lower) > tol) {

      if (f2 > f1) {
        # minimum is to the left of x2
        # x2 is the new upper bound
        # x1 is the new upper test point
        upper = x2
        x2 = x1
        f2 = f1

        # create new lower test point
        x1 = upper - ratio*(upper - lower)
        f1 = func(x1)
      }
      else {
        # minimum is to the right of x1
        # x1 is the new lower bound
        # x2 is the new lower test point
        lower = x1
        x1 = x2
        f1 = f2

        # create new upper test point
        x2 = lower + ratio*(upper - lower)
        f2 = func(x2)
      }
    }

    # check if the newest upper and lower bounds are 0 or 1
    if (upper == 1) {
      estimated.minimiser = upper
    }
    else if (lower == 0) {
      estimated.minimiser = lower
    }
    else {
      # mid-point of the final interval is the estimate of the minimiser
      estimated.minimiser = (lower + upper)/2
    }
    return(estimated.minimiser)
  }
  # =====================================================================================
  # BACK TO MAIN FUNCTION
  # golden section search to find optimal zeta
  zeta <- golden(calculate_ssq, 0, 1)

  # set flag manually at the beginning
  flag2 <- 0

  # =====================================================================================
  # OPTIMISE GIVEN ZETA FUNCTION
  OptimiseGivenZeta = function(sorted, theoretical_quantiles, zeta) {
    # run checks on the zeta
    if (zeta < 10e-13) {
      # degenerate case
      CaseID = "Degenerate"
      smallest = sorted[1]
      largest = sorted[length(sorted)]
      average = mean(sorted[2:(length(sorted) - 1)])
      # calculate the SSQ
      SSQ = sum(((sorted[2:(length(sorted) - 1)]) - average)^2)
      results = list(CaseID = CaseID, SSQ = SSQ, smallest = smallest, average = average, largest = largest, flags = c(flag1, flag2))
      class(results) = "Degenerate"
    }
    else if (zeta > 0.999) {
      # normal distribution
      CaseID = "Normal distribution"
      # regress xs against zs
      df_mod = data.frame(sorted = sorted, theoretical_quantiles = theoretical_quantiles)
      mod = lm(sorted ~ theoretical_quantiles, data = df_mod)
      # calculate SSQ, intercept and slope
      SSQ = sum(mod$residuals^2)
      intercept = mod$coefficients[1][[1]]
      slope = mod$coefficients[2][[1]]

      # check standard deviation
      if ((slope < 10e-13) & all(((1/(2*delta)) * (matrix(c(0, delta, delta, 0, -1, 1), nrow = 2, byrow = T))%*%(t(W)%*%W%*%(c(mean(sorted), 0, 0)) - t(W)%*%X)) <= 10e-13)) {
        CaseID = "Constant distribution"
        constant = mean(sorted)
        SSQ = sum((sorted - constant)^2)
        results = list(CaseID = CaseID, SSQ = SSQ, constant = constant, flags = c(flag1, flag2))
        class(results) = "Constant"
        return(results)
      }
      else {
        results = list(CaseID = CaseID, SSQ = SSQ, intercept = intercept, slope = slope, flags = c(flag1, flag2))
        class(results) = "Normal"
        return(results)
      }
    }
    else { # in between
      # change back to delta
      delta <- NULL
      delta <- zeta/(1 - zeta)
      # use next function to get the parameters
      results = Check4Cases(sorted, theoretical_quantiles, delta)
    }
    return(results)
  }
  # =====================================================================================
  # BACK TO MAIN FUNCTION
  # supply optimum zeta to next function that will optimise given zeta
  results = OptimiseGivenZeta(sorted, theoretical_quantiles, zeta)

  # =====================================================================================
  # SPECIFIC PRINTING FUNCTIONS
  # JSU
  printJSU = function(results) {
    # print the distribution
    cat("Johnson SU distribution")
    # print the SSQ
    cat("\nSSQ = ", results$SSQ)
    # print the results
    cat("\n\nParameter Estimates:\n")
    cat("delta =", results$delta)
    cat("\nxi =", results$xi)
    cat("\ngamma =", results$gamma)
    cat("\nlambda =", results$lambda)
    if (results$flags[1] == 1) {
      cat("\n\nFlag1: There is a possibility of multiple local minima within the Golden Section Search.")
    }
  }
  # Neg-lognormal
  printNegLognormal = function(results) {
    # print the distribution
    cat("Shifted negative lognormal distribution")
    # print the SSQ
    cat("\nSSQ = ", results$SSQ)
    # print the results
    cat("\n\nParameter Estimates:\n")
    cat("delta =", results$delta)
    cat("\nxi = ", results$xi)
    cat("\nlambda = ", results$lambda)
    if (results$flags[1] == 1) {
      cat("\n\nFlag1: There is a possibility of multiple local minima within the Golden Section Search.")
    }
    if (results$flags[2] == 1) {
      cat("\n\nFlag2: You have some observation(s) in your dataset that exceed the maximum value.")
    }
  }
  # Lognormal
  printLognormal = function(results) {
    # print the distribution
    cat("Shifted lognormal distribution")
    # print the SSQ
    cat("\nSSQ = ", results$SSQ)
    # print the results
    cat("\n\nParameter Estimates:\n")
    cat("delta =", results$delta)
    cat("\nxi = ", results$xi)
    cat("\nlambda = ", results$lambda)
    if (results$flags[1] == 1) {
      cat("\n\nFlag1: There is a possibility of multiple local minima within the Golden Section Search.")
    }
    if (results$flags[2] == 1) {
      cat("\n\nFlag2: You have some observation(s) in your dataset that exceed the minimum value.")
    }
  }
  # Normal
  printNormal = function(results) {
    # print the distribution
    cat("Normal distribution")
    # print the SSQ
    cat("\nSSQ = ", results$SSQ)
    # print the results
    cat("\n\nParameter Estimates:\n")
    cat("mean =", results$mean)
    cat("\nstandard deviation = ", results$std)
    if (results$flags[1] == 1) {
      cat("\n\nFlag1: There is a possibility of multiple local minima within the Golden Section Search.")
    }
  }
  # Constant
  printConstant = function(results) {
    # print the distribution
    cat("Constant distribution")
    # print the SSQ
    cat("\nSSQ = ", results$SSQ)
    # print the results
    cat("\n\nConstant = ", results$constant)
    if (results$flags[1] == 1) {
      cat("\n\nFlag1: There is a possibility of multiple local minima within the Golden Section Search.")
    }
  }
  # Degenerate
  printDegenerate = function(results) {
    # print the distribution
    cat("Degenerate case")
    # print the SSQ
    cat("\nSSQ = ", results$SSQ)
    # print the results
    cat("\n\nSmallest value = ", results$smallest)
    cat("\nLargest value = ", results$largest)
    cat("\nAverage value = ", results$average)
    if (results$flags[1] == 1) {
      cat("\n\nFlag1: There is a possibility of multiple local minima within the Golden Section Search.")
    }
  }
  # None
  printNone = function(results) {
    cat("None of the cases applied to the given dataset.")
  }
  # =====================================================================================
  # BACK TO MAIN FUNCTION
  if (class(results) == "JSU") {
    printJSU(results)
  }
  else if (class(results) == "NegLognormal") {
    printNegLognormal(results)
  }
  else if (class(results) == "Lognormal") {
    printLognormal(results)
  }
  else if (class(results) == "Normal") {
    printNormal(results)
  }
  else if (class(results) == "Constant") {
    printConstant(results)
  }
  else if (class(results) == "Degenerate") {
    printDegenerate(results)
  }
  else if (class(results) == "None") {
    printNone(results)
  }
  invisible(results)
}
