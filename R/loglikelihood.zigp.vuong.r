loglikelihood.zigp.vuong <-
function (delta)
{

    n <- get("n", pos=globalenv())
    k.beta <- get("k.beta", pos=globalenv())
    k.alpha <- get("k.alpha", pos=globalenv())
    k.gamma <- get("k.gamma", pos=globalenv())
    Xsave <- get("Xsave", pos=globalenv())
    Wsave <- get("Wsave", pos=globalenv())
    Zsave <- get("Zsave", pos=globalenv())
    Y <- get("Y", pos=globalenv())
    t.i <- get("t.i", pos=globalenv())

    eta.mu <- double(n)
    eta.phi <- double(n)
    eta.omega <- double(n)
    s1 <- double(n)
    s2 <- double(n)
    if (k.beta == 1) {
        eta.mu <- Xsave * delta[1]
    }
    else {
        beta <- delta[1:k.beta]
        eta.mu <- Xsave %*% beta
    }
    if (is.null(Wsave) == FALSE) {
        if (k.alpha == 1) {
            eta.phi <- Wsave * delta[k.beta + 1]
        }
        else {
            alpha <- delta[(k.beta + 1):(k.beta + k.alpha)]
            eta.phi <- Wsave %*% alpha
        }
    }
    if (is.null(Zsave) == FALSE) {
        if (k.gamma == 1) {
            eta.omega <- Zsave * delta[k.beta + k.alpha + 1]
        }
        else {
            gamma <- delta[(k.beta + k.alpha + 1):(k.beta + k.alpha +
                k.gamma)]
            eta.omega <- Zsave %*% gamma
        }
    }
    mu.i <- t.i * exp(eta.mu)
    if (is.null(Wsave) == FALSE) {
        b.i <- exp(eta.phi)
        phi.i <- 1 + b.i
    }
    else {
        b.i <- rep(0, n)
        phi.i <- rep(1, n)
    }
    if (is.null(Zsave) == FALSE) {
        k.i <- exp(eta.omega)
    }
    else {
        k.i <- rep(0, n)
    }
    temp1 <- (mu.i + b.i * Y)

    if (is.null(Zsave) == FALSE) {
        s1 <- ifelse(Y == 0, 1, 0) * (log(k.i + exp(-1/phi.i *
              mu.i)) - log(1 + k.i))
      }
      else {
        s1 <- ifelse(Y == 0, 1, 0) * (-1/phi.i *
            mu.i)
      }
      s2 <- ifelse(Y > 0, 1, 0) * (-log(1 + k.i) + log(t.i) +
          eta.mu + (Y - 1) * log(temp1) - Y * log(phi.i) - 1/phi.i *
          temp1 - lgamma(Y + 1))
    l <- s1 + s2
    return(l)
}

