mle.zigp.full.like <-
function (Yein, Xein, Offset = rep(1, length(Y)), summary = TRUE,
    plot = FALSE, method = "CG")
{
    X <<- Xein
    Y <<- Yein
    if (is.matrix(X)) {
        n <<- dim(X)[1]
        k <<- dim(X)[2]
    }
    else {
        n <<- length(X)
        k <<- 1
    }
    t.i <<- Offset
    alpha <- log(sqrt(var(Y)/mean(Y)))
    gamma <- -log(4)
    beta <- summary(glm(Y ~ offset(log(t.i)) + X - 1, family = poisson))$coefficients[,
        1]
    start.delta <- c(alpha, gamma, beta)
    
    loglikelihood.zigp.full(start.delta)
    
    opt <- optim(par = start.delta, fn = loglikelihood.zigp.full,
        gr = gradient1, method = method)
    delta <- opt$par
    it <- opt$counts
    loglikelihood <- -opt$value
    message <- opt$message
    coef <- double(k + 2)
    coef <- delta
    if (is.null(W) == FALSE) {
        phi <- 1 + exp(coef[1])
    }
    else {
        phi <- 1
    }
    if (is.null(Z) == FALSE) {
        omega <- exp(coef[2])/(1 + exp(coef[2]))
    }
    else {
        omega <- 0
    }
    beta <- coef[3:(k + 2)]
    fit <- fit.zigp1(delta)
    res <- double(n)
    RSS <- 0
    res <- Y - fit$fit
    RSS <- sum(res^2)
    AIC <- -2 * loglikelihood + 2 * (k + 2)
    range.mu <- c(min(fit$mu), max(fit$mu))
    ausgabe <<- list(ZI.Parameter = omega, Coefficients = beta,
        AIC = AIC, Dispersion.Parameter = phi, Range.mu = range.mu,
        Log.Likelihood = loglikelihood, Residuals = res, RSS = RSS,
        Iterations = it,
        Message = message, Response = Y,
        Fitted.Values = fit$fit, Design = X)
    if (plot) {
        plot.zigp1(ausgabe)
        if (summary)
            summaryzigp1(ausgabe)
    }
    else {
        if (summary) {
            summaryzigp1(ausgabe)
        }
        else {
            return(ausgabe)
        }
    }
}