optimized.run <-
function (Yein, Xein, Wein, Zein)
{

    n <- get("n", pos=globalenv())
    t.i <- get("t.i", pos=globalenv())

    out <- mle.zigp.full.like(Yein, Xein, Offset = t.i, summary = FALSE)
    beta.start <- out$Coefficients
    if (is.null(Wein) == FALSE) {
        phi.first <- out$Dispersion.Parameter
    }
    if (is.null(Zein) == FALSE) {
        omega.first <- out$ZI.Parameter
    }
    if (is.null(Wein) == FALSE) {
        rechte.seite <- rep(log(phi.first - 1), n)
        out <- lm(rechte.seite ~ Wein - 1)
        alpha.start <- out$coefficients
    }
    else {
        alpha.start <- NULL
    }
    if (is.null(Zein) == FALSE) {
        rechte.seite <- rep(log(omega.first) - log(1 - omega.first),
            n)
        out <- lm(rechte.seite ~ Zein - 1)
        gamma.start <- out$coefficients
    }
    else {
        gamma.start <- NULL
    }
    start.delta <- c(beta.start, alpha.start, gamma.start)
    rm(beta.start, alpha.start, gamma.start)
    return(start.delta)
}
