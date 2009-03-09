clarke <- 
function (model1, model2, alpha = 0.05, correction = T)
{
    logliknb <- function(n, th, mu, y) {
        lgamma(th + y) - lgamma(th) - lgamma(y + 1) +
            th * log(th) + y * log(mu + (y == 0)) - (th + y) *
            log(th + mu)
    }
    loglikpoi <- function(n, mu, y, model) {
        -mu + y * model$family$linkfun(mu) - lgamma(y + 1)
    }

    if (is.na(match("ZI.Parameter",names(model1)))) {
      if (is.na(charmatch("Negative Binomial",model1$family$family))==FALSE) {
        mu <- model1$fitted.values
        th <- model1$theta
        Y <- model1$y
        n <- length(Y)
        ll1 <- logliknb(n, th, mu, Y)
        p <- length(model1$coefficients)+1
      }
      if (is.na(charmatch("poisson",model1$family$family))==FALSE) {
        mu <- model1$fitted.values
        Y <- model1$y
        n <- length(Y)
        ll1 <- loglikpoi(n, mu, Y, model1)
        p <- length(model1$coefficients)
      }
    }
    else {
      ll1 <- model1$Loglike.Vuong
      p <- length(model1$Coefficients.Mu) + length(model1$Coefficients.Phi) +
           length(model1$Coefficients.Omega)
    }

    if (is.na(match("ZI.Parameter",names(model2)))) {
      if (is.na(charmatch("Negative Binomial",model2$family$family))==FALSE) {
        mu <- model2$fitted.values
        th <- model2$theta
        Y <- model2$y
        n <- length(Y)
        ll2 <- logliknb(n, th, mu, Y)
        q <- length(model2$coefficients)+1
      }
      if (is.na(charmatch("poisson",model2$family$family))==FALSE) {
        mu <- model2$fitted.values
        Y <- model2$y
        n <- length(Y)
        ll2 <- loglikpoi(n, mu, Y, model2)
        q <- length(model2$coefficients)
      }
    }
    else {
      ll2 <- model2$Loglike.Vuong
      q <- length(model2$Coefficients.Mu) + length(model2$Coefficients.Phi) +
           length(model2$Coefficients.Omega)
    }

    n <- length(ll1)
    m.i <- ll1 - ll2 - ifelse(correction, p/2 * log(n) - q/2 * log(n), 0)/n
    B <- sum(m.i > 0)
    if (B < n/2) {
        cat("Perform lower tail test:", "\n")
        p.value <- pbinom(B, n, 0.5)
        p.value2 <- ifelse(p.value < 10^(-16),
          "<2e-16", as.character(formatC(p.value,
          ifelse(p.value < 10^(-4), 3, 4),
          format=ifelse(p.value < 10^(-4), "g", "f"))))
        cat("H0: 'Models equally good'   vs.   H1: 'Model 2 better than model 1'", "\n")
        cat(ifelse(p.value > alpha, "H0 cannot be rejected on level alpha.",
            "Reject H0 on level alpha."), "\n")
        cat("P value: ", p.value2, "\n")
    }
    if (B >= n/2) {
        cat("Perform upper tail test:", "\n")
        p.value <- 1 - pbinom(B - 1, n, 0.5)
        p.value2 <- ifelse(p.value < 10^(-16),
          "<2e-16", as.character(formatC(p.value,
          ifelse(p.value < 10^(-4), 3, 4),
          format=ifelse(p.value < 10^(-4), "g", "f"))))
        cat("H0: 'Models equally good'   vs.   H1: 'Model 1 better than model 2'", "\n")
        cat(ifelse(p.value > alpha, "H0 cannot be rejected on level alpha.",
            "Reject H0 on level alpha."), "\n")
        cat("P value: ", p.value2, "\n")
    }
}

