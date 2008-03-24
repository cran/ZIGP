loglikelihood.zigp.full <-
function (delta)
{
    eta <- double(n)
    s1 <- double(1)
    s2 <- double(1)
    if (k == 1) {
        eta <- delta[3] * X
    }
    else {
        beta <- delta[3:(k + 2)]
        eta <- X %*% beta
    }
    if (is.null(W) == FALSE) { phi <- 1 + exp(delta[1]) }
    else { phi <- 1 }
    if (is.null(Z) == FALSE) { omega <- exp(delta[2])/(1 + exp(delta[2])) }
    else { omega <- 0 }
    if (omega>0) {
    s1 <- sum(ifelse(Y == 0, 1, 0) * log(omega + (1 - omega) *
        exp(-t.i * exp(eta)/phi))) }
    else { s1 <- sum(ifelse(Y == 0, 1, 0) * (-t.i * exp(eta)/phi)) }
    s2 <- sum(ifelse(Y > 0, 1, 0) * (log(1 - omega) + log(t.i) +
        eta + (Y - 1) * log(t.i * exp(eta) + (phi - 1) * Y) -
        Y * log(phi) - (t.i * exp(eta) + (phi - 1) * Y)/phi))
    temp2 <- ifelse(Y < 1, 1, Y)
    temp4 <- 0
    for (i in 1:n) {
        temp3 <- c(1:temp2[i])
        temp4 <- temp4 + sum(log(temp3))
    }
    l <- s1 + s2 - temp4
    return(-l)
}