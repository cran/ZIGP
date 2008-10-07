loglikelihood.zigp.full <-
function (delta)
{
    n <- get("n", pos=globalenv())
    k <- get("k", pos=globalenv())
    Xsave <- get("Xsave", pos=globalenv())
    Wsave <- get("Wsave", pos=globalenv())
    Zsave <- get("Zsave", pos=globalenv())
    Y <- get("Y", pos=globalenv())
    t.i <- get("t.i", pos=globalenv())

    eta <- double(n)
    s1 <- double(1)
    s2 <- double(1)
    if (k == 1) {
        eta <- delta[3] * Xsave
    }
    else {
        beta <- delta[3:(k + 2)]
        eta <- Xsave %*% beta
    }
    if (is.null(Wsave) == FALSE) { phi <- 1 + exp(delta[1]) }
    else { phi <- 1 }
    if (is.null(Zsave) == FALSE) {
      if (exp(delta[2]) == Inf) { omega <- 1 }
      else { omega <- exp(delta[2])/(1 + exp(delta[2])) }
    }
    else { omega <- 0 }

    if (omega>0) {
      s1 <- sum(ifelse(Y == 0, 1, 0) * log(omega + (1 - omega) *
            exp(-t.i * exp(eta)/phi))) }
    else { s1 <- sum(ifelse(Y == 0, 1, 0) * (-t.i * exp(eta)/phi)) }
    s2 <- sum(ifelse(Y > 0, 1, 0) * (log(1 - omega) + log(t.i) +
        eta + (Y - 1) * log(t.i * exp(eta) + (phi - 1) * Y) -
        Y * log(phi) - (t.i * exp(eta) + (phi - 1) * Y)/phi- lgamma(Y + 1)))
    l <- s1 + s2
    return(-l)
}

