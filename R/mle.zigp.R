mle.zigp <-
function (Yin, Xin, Win = NULL, Zin = NULL, Offset = rep(1, length(Yin)),
    init = TRUE)
{

    assign("Y",Yin,.GlobalEnv)
    Y <- get("Y", pos=globalenv())
    assign("X",Xin,.GlobalEnv)
    X <- get("X", pos=globalenv())
    assign("W",Win,.GlobalEnv)
    W <- get("W", pos=globalenv())
    assign("Z",Zin,.GlobalEnv)
    Z <- get("Z", pos=globalenv())
    if (is.matrix(X)) {
        assign("n",dim(X)[1],.GlobalEnv)
        n <- get("n", pos=globalenv())
        assign("k.beta",dim(X)[2],.GlobalEnv)
        k.beta <- get("k.beta", pos=globalenv())
    }
    else {
        assign("n",length(X),.GlobalEnv)
        n <- get("n", pos=globalenv())
        assign("k.beta",1,.GlobalEnv)
        k.beta <- get("k.beta", pos=globalenv())
    }
    if (is.null(W) == FALSE) {
        if (is.matrix(W)) {
            assign("k.alpha",dim(W)[2],.GlobalEnv)
            k.alpha <- get("k.alpha", pos=globalenv())
        }
        else {
            assign("k.alpha",1,.GlobalEnv)
            k.alpha <- get("k.alpha", pos=globalenv())
        }
    }
    else { 
        assign("k.alpha",0,.GlobalEnv)
        k.alpha <- get("k.alpha", pos=globalenv())
    }
    if (is.null(Z) == FALSE) {
        if (is.matrix(Z)) {
            assign("k.gamma",dim(Z)[2],.GlobalEnv)
            k.gamma <- get("k.gamma", pos=globalenv())
        }
        else {
            assign("k.gamma",1,.GlobalEnv)
            k.gamma <- get("k.gamma", pos=globalenv())
        }
    }
    else { 
        assign("k.gamma",0,.GlobalEnv)
        k.gamma <- get("k.gamma", pos=globalenv())
    }
    assign("t.i",Offset,.GlobalEnv)
    t.i <- get("t.i", pos=globalenv())
    if (init) {

        start.delta <- optimized.run(Y, X, W, Z)

        opt <- optim(par = start.delta, fn = loglikelihood.zigp,
            gr = gradient, method = "BFGS")
        # model white noise onto start values if to close to zero (for OD and ZI designs)
        if(sum(start.delta==opt$par)==k.beta+k.alpha+k.gamma) {
          beta <- start.delta[1:k.beta]
          if (k.alpha>0) { alpha <- start.delta[(k.beta + 1):(k.beta + k.alpha)]
                           alpha[alpha<0.01] <- runif(sum(alpha<0.01),-1,1)/sum(alpha<0.01) }
          else { alpha <- NULL }
          if (k.gamma>0) { gamma <- start.delta[(k.beta + k.alpha + 1):(k.beta + k.alpha + k.gamma)]
                           gamma[gamma<0.01] <- runif(sum(gamma<0.01),0,1)/sum(gamma<0.01) }
          else { gamma <- NULL }
          start.delta2 <- c(beta,alpha,gamma)
          test <- loglikelihood.zigp(start.delta2)

          counter <- 0
          while ((is.na(test) | is.nan(test) | test==-Inf | test==Inf) & counter < 30) {
            counter <- counter + 1
            if (k.alpha>0) { alpha <- start.delta[(k.beta + 1):(k.beta + k.alpha)]
                             alpha[alpha<0.01] <- runif(sum(alpha<0.01),0,1)/sum(alpha<0.01) }
            else { alpha <- NULL }
            if (k.gamma>0) { gamma <- start.delta[(k.beta + k.alpha + 1):(k.beta + k.alpha + k.gamma)]
                             gamma[gamma<0.01] <- runif(sum(gamma<0.01),0,1)/sum(gamma<0.01) }
            else { gamma <- NULL }
            start.delta2 <- c(beta,alpha,gamma)
            test <- loglikelihood.zigp(start.delta2)
          }
          if (counter<30) { start.delta <- start.delta2 }
          opt <- optim(par = start.delta, fn = loglikelihood.zigp,
              gr = gradient, method = "BFGS")
        }
    }
    else {
        if (is.null(W) == FALSE) {
            alpha <- rep(-2.5, k.alpha)
        }
        else {
            alpha <- NULL
        }
        if (is.null(Z) == FALSE) {
            gamma <- rep(-2.5, k.gamma)
        }
        else {
            gamma <- NULL
        }
        beta <- summary(glm(Y ~ offset(log(t.i)) + 1 + X, family = poisson(link = log)))$coefficients[,
            1]
        start.delta <- c(beta, alpha, gamma)
        opt <- optim(par = start.delta, fn = loglikelihood.zigp,
            gr = gradient, method = "BFGS")
    }
    delta <- opt$par
    it <- opt$counts
    loglikelihood <- -opt$value
    llvuong <- loglikelihood.zigp.vuong(delta)

    if (is.null(opt$message)) {
        message <- "NULL"
    }
    else {
        message <- opt$message
    }
    AIC <- -2 * loglikelihood + 2 * (k.beta + k.alpha + k.gamma)
    coef <- delta
    beta <- coef[1:k.beta]
    if (is.null(W) == FALSE) {
        alpha <- coef[(k.beta + 1):(k.beta + k.alpha)]
    }
    else {
        alpha <- NULL
    }
    if (is.null(Z) == FALSE) {
        gamma <- coef[(k.beta + k.alpha + 1):(k.beta + k.alpha +
            k.gamma)]
    }
    else {
        gamma <- NULL
    }
    fit <- fit.zigp(delta)
    res <- Y - fit$fit
    chsq <- sum((Y - fit$fit)^2/((1 - fit$omega) * fit$mu * (fit$phi^2 +
        fit$mu * fit$omega)))
    range.mu <- c(min(fit$mu), max(fit$mu))
    if (is.null(Z) == FALSE) {
        if (is.matrix(Z)) {
            eta.omega <- Z %*% gamma
        }
        else {
            eta.omega <- Z * gamma
        }
        omega <- exp(eta.omega)/(1 + exp(eta.omega))
    }
    else {
        omega <- rep(0, n)
    }
    range.omega <- c(min(omega), max(omega))
    if (is.null(W) == FALSE) {
        if (is.matrix(W)) {
            eta.phi <- W %*% alpha
        }
        else {
            eta.phi <- W * alpha
        }
        phi <- 1 + exp(eta.phi)
    }
    else {
        phi <- rep(1, n)
    }
    range.phi <- c(min(phi), max(phi))
    mle.data <- list(ZI.Parameter = omega, Coefficients.Mu = beta,
        Coefficients.Phi = alpha, Coefficients.Omega = gamma,
        Dispersion.Parameter = phi, Range.Mu = range.mu, Range.Phi = range.phi,
        Range.Omega = range.omega, Log.Likelihood = loglikelihood,
        Residuals = res, Pearson = chsq, AIC = AIC, Iterations = it,
        Message = message, Response = Y,
        Fitted.Values = fit$fit, Design.Mu = X, Design.Phi = W,
        Design.Omega = Z,
        Loglike.Vuong = llvuong)
    return(mle.data)
}

