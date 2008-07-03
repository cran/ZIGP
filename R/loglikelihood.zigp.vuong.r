loglikelihood.zigp.vuong <-
function (delta)
{

    n <- get("n")
    k.beta <- get("k.beta")
    k.alpha <- get("k.alpha")
    k.gamma <- get("k.gamma")
    X <- get("X")
    W <- get("W")
    Z <- get("Z")
    Y <- get("Y")
    t.i <- get("t.i")

    eta.mu <- double(n)
    eta.phi <- double(n)
    eta.omega <- double(n)
    s1 <- double(n)
    s2 <- double(n)
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
