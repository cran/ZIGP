loglikelihood.zigp <-
function (delta)
{

    n <- get("n", pos=globalenv())
    k.beta <- get("k.beta", pos=globalenv())
    k.alpha <- get("k.alpha", pos=globalenv())
    k.gamma <- get("k.gamma", pos=globalenv())
    X <- get("X", pos=globalenv())
    W <- get("W", pos=globalenv())
    Z <- get("Z", pos=globalenv())
    Y <- get("Y", pos=globalenv())
    t.i <- get("t.i", pos=globalenv())

    eta.mu <- double(n)
    eta.phi <- double(n)
    eta.omega <- double(n)
    s1 <- double(1)
    s2 <- double(1)
    if (k.beta == 1) {
        eta.mu <- X * delta[1]
    }
    else {
        beta <- delta[1:k.beta]
        eta.mu <- X %*% beta
    }
    if (is.null(W) == FALSE) {
        if (k.alpha == 1) {
            eta.phi <- W * delta[k.beta + 1]
        }
        else {
            alpha <- delta[(k.beta + 1):(k.beta + k.alpha)]
            eta.phi <- W %*% alpha
        }
    }
    if (is.null(Z) == FALSE) {
        if (k.gamma == 1) {
            eta.omega <- Z * delta[k.beta + k.alpha + 1]
        }
        else {
            gamma <- delta[(k.beta + k.alpha + 1):(k.beta + k.alpha +
                k.gamma)]
            eta.omega <- Z %*% gamma
        }
    }
    mu.i <- t.i * exp(eta.mu)
    if (is.null(W) == FALSE) {
        b.i <- exp(eta.phi)
        phi.i <- 1 + b.i
    }
    else {
        b.i <- rep(0, n)
        phi.i <- rep(1, n)
    }
    if (is.null(Z) == FALSE) {
        k.i <- exp(eta.omega)
    }
    else {
        k.i <- rep(0, n)
    }
    temp1 <- (mu.i + b.i * Y)

    if (is.null(Z) == FALSE) {
      s1 <- sum(ifelse(Y == 0, 1, 0) * (log(k.i + exp(-1/phi.i *
          mu.i)) - log(1 + k.i)))
    }
    else {
      s1 <- sum(ifelse(Y == 0, 1, 0) * (-1/phi.i *
            mu.i))
    }
    s2 <- sum(ifelse(Y > 0, 1, 0) * (-log(1 + k.i) + log(t.i) +
        eta.mu + (Y - 1) * log(temp1) - Y * log(phi.i) - 1/phi.i *
        temp1 - lgamma(Y + 1)))
    l <- s1 + s2
    return(-l)
}

