mle.zigp.full.like <-
function (Yein, Xein, Offset = rep(1, length(Yein)), summary = TRUE,
    plot = FALSE, method = "CG")
{
    assign("Ysave", Yein, .GlobalEnv)
    Ysave <- get("Ysave", pos=globalenv())
    assign("Xsave", Xein, .GlobalEnv)
    Xsave <- get("Xsave", pos=globalenv())
    Wsave <- get("Wsave", pos=globalenv())
    Zsave <- get("Zsave", pos=globalenv())
    if (is.matrix(Xsave)) {
        assign("n",dim(Xsave)[1],.GlobalEnv)
        n <- get("n", pos=globalenv())
        assign("ksave",dim(Xsave)[2],.GlobalEnv)
        k <- get("ksave", pos=globalenv())
    }
    else {
        assign("n",length(Xsave),.GlobalEnv)
        n <- get("n", pos=globalenv())
        assign("ksave",1,.GlobalEnv)
        k <- get("ksave", pos=globalenv())
    }
    assign("t.i",Offset,.GlobalEnv)
    t.i <- get("t.i", pos=globalenv())
    alpha <- log(sqrt(var(Yein)/mean(Yein))-1)
    gamma <- -log(4)
    beta <- summary(glm(Yein ~ offset(log(t.i)) + Xsave - 1, family = poisson))$coefficients[,
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
    if (is.null(Wsave) == FALSE) {
        phi <- 1 + exp(coef[1])
    }
    else {
        phi <- 1
    }
    if (is.null(Zsave) == FALSE) {
        omega <- exp(coef[2])/(1 + exp(coef[2]))
    }
    else {
        omega <- 0
    }
    beta <- coef[3:(k + 2)]
    fit <- fit.zigp1(delta)
    res <- double(n)
    RSS <- 0
    res <- Yein - fit$fit
    RSS <- sum(res^2)
    AIC <- -2 * loglikelihood + 2 * (k + 2)
    range.mu <- c(min(fit$mu), max(fit$mu))
    assign("ausgabe",list(ZI.Parameter = omega, Coefficients = beta,
        AIC = AIC, Dispersion.Parameter = phi, Range.mu = range.mu,
        Log.Likelihood = loglikelihood, Residuals = res, RSS = RSS,
        Iterations = it,
        Message = message, Response = Yein,
        Fitted.Values = fit$fit, Design = Xsave),.GlobalEnv)
    ausgabe <- get("ausgabe", pos=globalenv())
    #if (plot) {
    #    plot.zigp1(ausgabe)
    #    if (summary)
    #        summaryzigp1(ausgabe)
    #}
    #else {
        if (summary) {
            summaryzigp1(ausgabe)
        }
        else {
            return(ausgabe)
        }
    #}
}
