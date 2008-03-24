optimized.run <-
function (Y, X, W, Z)
{
    out <- mle.zigp.full.like(Y, X, Offset = t.i, summary = FALSE)
    beta.start <- out$Coefficients
    if (is.null(W) == FALSE) {
        phi.first <- out$Dispersion.Parameter
    }
    if (is.null(Z) == FALSE) {
        omega.first <- out$ZI.Parameter
    }
    if (is.null(W) == FALSE) {
        rechte.seite <- rep(log(phi.first - 1), n)
        out <- lm(rechte.seite ~ W - 1)
        alpha.start <- out$coefficients
    }
    else {
        alpha.start <- NULL
    }
    if (is.null(Z) == FALSE) {
        rechte.seite <- rep(log(omega.first) - log(1 - omega.first),
            n)
        out <- lm(rechte.seite ~ Z - 1)
        gamma.start <- out$coefficients
    }
    else {
        gamma.start <- NULL
    }
    start.delta <- c(beta.start, alpha.start, gamma.start)
    rm(beta.start, alpha.start, gamma.start)
    return(start.delta)
}
