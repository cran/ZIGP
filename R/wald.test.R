wald.test <-
function (Yin, Xin, Win = NULL, Zin = NULL, Offset = rep(1, length(Yin)),
    init = T)
{
    assign("Y",Yin,.GlobalEnv)
    Y <- get("Y", pos=globalenv())
    assign("Xsave",Xin,.GlobalEnv)
    Xsave <- get("Xsave", pos=globalenv())
    assign("Wsave",Win,.GlobalEnv)
    Wsave <- get("Wsave", pos=globalenv())
    assign("Zsave",Zin,.GlobalEnv)
    Zsave <- get("Zsave", pos=globalenv())
    assign("k.beta",dim(Xsave)[2],.GlobalEnv)
    k.beta <- get("k.beta", pos=globalenv())
    if (is.null(Wsave) == FALSE) {
        if (is.matrix(Wsave)) {
            assign("k.alpha",dim(Wsave)[2],.GlobalEnv)
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
    if (is.null(Zsave) == FALSE) {
        if (is.matrix(Zsave)) {
            assign("k.gamma",dim(Zsave)[2],.GlobalEnv)
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
    assign("n",length(Y),.GlobalEnv)
    n <- get("n", pos=globalenv())
    ausgabe <- mle.zigp(Y, Xsave, Wsave, Zsave, Offset = Offset, init = init)
    hat.beta <- ausgabe$Coefficients.Mu
    if (is.null(Wsave) == FALSE) {
        hat.alpha <- ausgabe$Coefficients.Phi
    }
    else {
        hat.alpha <- NULL
    }
    if (is.null(Zsave) == FALSE) {
        hat.gamma <- ausgabe$Coefficients.Omega
    }
    else {
        hat.gamma <- NULL
    }
    B <- FM(hat.beta, hat.alpha, hat.gamma, Xsave, Wsave, Zsave, Offset = Offset)
    sd.vector <- sqrt(diag(solve(B, tol = 1e-50)))
    hat.sd.beta <- sd.vector[1:k.beta]
    if (is.null(Wsave) == FALSE) {
        hat.sd.alpha <- sd.vector[(k.beta + 1):(k.beta + k.alpha)]
    }
    else {
        hat.sd.alpha <- NULL
    }
    if (is.null(Zsave) == FALSE) {
        hat.sd.gamma <- sd.vector[(k.beta + k.alpha + 1):(k.beta +
            k.alpha + k.gamma)]
    }
    else {
        hat.sd.gamma <- NULL
    }
    z.stat.beta <- hat.beta/hat.sd.beta
    if (is.null(Wsave) == FALSE) {
        z.stat.alpha <- hat.alpha/hat.sd.alpha
    }
    if (is.null(Zsave) == FALSE) {
        z.stat.gamma <- hat.gamma/hat.sd.gamma
    }
    p.value.beta <- 2 * pnorm(-abs(z.stat.beta))
    if (is.null(Wsave) == FALSE) {
        p.value.alpha <- 2 * pnorm(-abs(z.stat.alpha))
    }
    if (is.null(Zsave) == FALSE) {
        p.value.gamma <- 2 * pnorm(-abs(z.stat.gamma))
    }
    glimpse.beta <- rep("", length(z.stat.beta))
    if (is.null(Wsave) == FALSE) {
        glimpse.alpha <- rep("", length(z.stat.alpha))
    }
    else {
        glimpse.alpha <- NULL
    }
    if (is.null(Zsave) == FALSE) {
        glimpse.gamma <- rep("", length(z.stat.gamma))
    }
    else {
        glimpse.gamma <- NULL
    }
    for (i in 1:length(z.stat.beta)) {
        if (p.value.beta[i] < 0.001) {
            glimpse.beta[i] <- "***"
        }
        if (p.value.beta[i] >= 0.001 & p.value.beta[i] < 0.01) {
            glimpse.beta[i] <- "**"
        }
        if (p.value.beta[i] >= 0.01 & p.value.beta[i] < 0.05) {
            glimpse.beta[i] <- "*"
        }
        if (p.value.beta[i] >= 0.05 & p.value.beta[i] < 0.1) {
            glimpse.beta[i] <- "."
        }
    }
    if (is.null(Wsave) == FALSE) {
        for (i in 1:length(z.stat.alpha)) {
            if (p.value.alpha[i] < 0.001) {
                glimpse.alpha[i] <- "***"
            }
            if (p.value.alpha[i] >= 0.001 & p.value.alpha[i] <
                0.01) {
                glimpse.alpha[i] <- "**"
            }
            if (p.value.alpha[i] >= 0.01 & p.value.alpha[i] <
                0.05) {
                glimpse.alpha[i] <- "*"
            }
            if (p.value.alpha[i] >= 0.05 & p.value.alpha[i] <
                0.1) {
                glimpse.alpha[i] <- "."
            }
        }
    }
    if (is.null(Zsave) == FALSE) {
        for (i in 1:length(z.stat.gamma)) {
            if (p.value.gamma[i] < 0.001) {
                glimpse.gamma[i] <- "***"
            }
            if (p.value.gamma[i] >= 0.001 & p.value.gamma[i] <
                0.01) {
                glimpse.gamma[i] <- "**"
            }
            if (p.value.gamma[i] >= 0.01 & p.value.gamma[i] <
                0.05) {
                glimpse.gamma[i] <- "*"
            }
            if (p.value.gamma[i] >= 0.05 & p.value.gamma[i] <
                0.1) {
                glimpse.gamma[i] <- "."
            }
        }
    }
    coef.names.beta <- c("", paste("b", c(0:(k.beta - 1)), sep = ""))
    if (is.null(Wsave) == FALSE) {
        coef.names.alpha <- c("", paste("a", c(0:(k.alpha - 1)),
            sep = ""))
    }
    else {
        coef.names.alpha <- NULL
    }
    if (is.null(Zsave) == FALSE) {
        coef.names.gamma <- c("", paste("g", c(0:(k.gamma - 1)),
            sep = ""))
    }
    else {
        coef.names.gamma <- NULL
    }
    coef.desc.beta <- double(k.beta)
    if (is.null(Wsave) == FALSE) {
        coef.desc.alpha <- double(k.alpha)
    }
    else {
        coef.desc.alpha <- NULL
    }
    if (is.null(Zsave) == FALSE) {
        coef.desc.gamma <- double(k.gamma)
    }
    else {
        coef.desc.gamma <- NULL
    }
    for (i in 1:k.beta) {
        coef.desc.beta[i] <- colnames(Xsave, do.NULL = FALSE)[i]
        if (is.matrix(Xsave)) {
            if (max(Xsave[, i]) == 1 & min(Xsave[, i]) == 1) {
                coef.desc.beta[i] <- "Intercept"
            }
        }
        else {
            if (max(Xsave[i]) == 1 & min(Xsave[i]) == 1) {
                coef.desc.beta[i] <- "Intercept"
            }
        }
    }
    coef.desc.beta <- c("MU REGRESSION", coef.desc.beta)
    if (is.null(Wsave) == FALSE) {
        for (i in 1:k.alpha) {
            coef.desc.alpha[i] <- colnames(Wsave, do.NULL = FALSE)[i]
            if (is.matrix(Wsave)) {
                if (max(Wsave[, i]) == 1 & min(Wsave[, i]) == 1) {
                  coef.desc.alpha[i] <- "Intercept"
                }
            }
            else {
                if (max(Wsave[i]) == 1 & min(Wsave[i]) == 1) {
                  coef.desc.alpha[i] <- "Intercept"
                }
            }
        }
        coef.desc.alpha <- c("PHI REGRESSION", coef.desc.alpha)
    }
    if (is.null(Zsave) == FALSE) {
        for (i in 1:k.gamma) {
            coef.desc.gamma[i] <- colnames(Zsave, do.NULL = FALSE)[i]
            if (is.matrix(Zsave)) {
                if (max(Zsave[, i]) == 1 & min(Zsave[, i]) == 1) {
                  coef.desc.gamma[i] <- "Intercept"
                }
            }
            else {
                if (max(Zsave[i]) == 1 & min(Zsave[i]) == 1) {
                  coef.desc.gamma[i] <- "Intercept"
                }
            }
        }
        coef.desc.gamma <- c("OMEGA REGRESSION", coef.desc.gamma)
    }
    p.value.beta2 <- p.value.beta
    if (is.null(Wsave) == FALSE) { p.value.alpha2 <- p.value.alpha }
    if (is.null(Zsave) == FALSE) { p.value.gamma2 <- p.value.gamma }
    for (i in 1:k.beta) {
      p.value.beta2[i] <- ifelse(p.value.beta[i]<10^(-16),"<2e-16",
              as.character(formatC(p.value.beta[i],
              ifelse(p.value.beta[i]<10^(-4),3,4),
              format=ifelse(p.value.beta[i]<10^(-4),"g","f")))) }
      p.value.beta <- p.value.beta2
    if (is.null(Wsave) == FALSE) {
        hat.alpha <- c("", as.character(formatC(hat.alpha, 5,
            format = "f")))
        hat.sd.alpha <- c("", as.character(formatC(hat.sd.alpha,
            5, format = "f")))
        z.stat.alpha <- c("", as.character(formatC(z.stat.alpha,
            5, format = "f")))
        for (i in 1:k.alpha) {
          p.value.alpha2[i] <- ifelse(p.value.alpha[i]<10^(-16),"<2e-16",
              as.character(formatC(p.value.alpha[i],
              ifelse(p.value.alpha[i]<10^(-4),3,4),
              format=ifelse(p.value.alpha[i]<10^(-4),"g","f")))) }
          p.value.alpha <- c("", p.value.alpha2)
        glimpse.alpha <- c("", glimpse.alpha)
    }
    else {
        hat.alpha <- NULL
        hat.sd.alpha <- NULL
        z.stat.alpha <- NULL
        p.value.alpha <- NULL
    }
    if (is.null(Zsave) == FALSE) {
        hat.gamma <- c("", as.character(formatC(hat.gamma, 5,
            format = "f")))
        hat.sd.gamma <- c("", as.character(formatC(hat.sd.gamma,
            5, format = "f")))
        z.stat.gamma <- c("", as.character(formatC(z.stat.gamma,
            5, format = "f")))
        for (i in 1:k.gamma) {
          p.value.gamma2[i] <- ifelse(p.value.gamma[i]<10^(-16),"<2e-16",
              as.character(formatC(p.value.gamma[i],
              ifelse(p.value.gamma[i]<10^(-4),3,4),
              format=ifelse(p.value.gamma[i]<10^(-4),"g","f")))) }
          p.value.gamma <- c("", p.value.gamma2)
        glimpse.gamma <- c("", glimpse.gamma)
    }
    else {
        hat.gamma <- NULL
        hat.sd.gamma <- NULL
        z.stat.gamma <- NULL
        p.value.gamma <- NULL
    }
    k.alpha <- length(coef.names.alpha)
    k.gamma <- length(coef.names.gamma)
    output <- matrix("", 1 + k.beta + k.alpha + k.gamma + 12,
        7)
    output[2:(1 + k.beta + k.alpha + k.gamma + 1), 1] <- c(coef.names.beta,
        coef.names.alpha, coef.names.gamma)
    output[2:(1 + k.beta + k.alpha + k.gamma + 1), 2] <- c(coef.desc.beta,
        coef.desc.alpha, coef.desc.gamma)
    output[1:(1 + k.beta + k.alpha + k.gamma + 1), 3] <- c("Estimate",
        "", formatC(hat.beta, 5, format = "f"), hat.alpha, hat.gamma)
    output[1:(1 + k.beta + k.alpha + k.gamma + 1), 4] <- c("Std. Error",
        "", formatC(hat.sd.beta, 5, format = "f"), hat.sd.alpha,
        hat.sd.gamma)
    output[1:(1 + k.beta + k.alpha + k.gamma + 1), 5] <- c("z value",
        "", formatC(z.stat.beta, 4, format = "f"), z.stat.alpha,
        z.stat.gamma)
    output[1:(1 + k.beta + k.alpha + k.gamma + 1), 6] <- c("Pr(>|z|)","",
        p.value.beta, p.value.alpha, p.value.gamma)
    output[1:(1 + k.beta + k.alpha + k.gamma + 1), 7] <- c("",
        "", glimpse.beta, glimpse.alpha, glimpse.gamma)
    output[(1 + k.beta + k.alpha + k.gamma + 3), 2] <- "Signif. codes: 0"
    output[(1 + k.beta + k.alpha + k.gamma + 3), 3] <- "`***' 0.001"
    output[(1 + k.beta + k.alpha + k.gamma + 3), 4] <- "`**'  0.01"
    output[(1 + k.beta + k.alpha + k.gamma + 3), 5] <- "`*'  0.05"
    output[(1 + k.beta + k.alpha + k.gamma + 3), 6] <- "`.'  0.1"
    output[(1 + k.beta + k.alpha + k.gamma + 3), 7] <- "` ' 1"
    output[(1 + k.beta + k.alpha + k.gamma + 4), 2] <- "Iterations"
    output[(1 + k.beta + k.alpha + k.gamma + 4), 4] <- ausgabe$Iterations[1]
    output[(1 + k.beta + k.alpha + k.gamma + 5), 2] <- "Log Likelihood"
    output[(1 + k.beta + k.alpha + k.gamma + 5), 4] <- formatC(ausgabe$Log.Likelihood,
        digits = 1, format = "f")
    output[(1 + k.beta + k.alpha + k.gamma + 6), 2] <- "Pearson Chi Squared"
    output[(1 + k.beta + k.alpha + k.gamma + 6), 4] <-
        formatC(ausgabe$Pearson,ifelse(ausgabe$Pearson>10^(6),3,0),format=ifelse(ausgabe$Pearson>10^(6),"g","f"))
    output[(1 + k.beta + k.alpha + k.gamma + 7), 2] <- "AIC"
    output[(1 + k.beta + k.alpha + k.gamma + 7), 4] <- round(ausgabe$AIC)
    output[(1 + k.beta + k.alpha + k.gamma + 8), 2] <- "Range Mu"
    output[(1 + k.beta + k.alpha + k.gamma + 8), 4] <- formatC(ausgabe$Range.Mu[1],
        digits = 2, format = "f")
    output[(1 + k.beta + k.alpha + k.gamma + 8), 5] <- formatC(ausgabe$Range.Mu[2],
        digits = 2, format = "f")
    output[(1 + k.beta + k.alpha + k.gamma + 9), 2] <- "Range Phi"
    output[(1 + k.beta + k.alpha + k.gamma + 9), 4] <- formatC(ausgabe$Range.Phi[1],
        digits = 2, format = "f")
    output[(1 + k.beta + k.alpha + k.gamma + 9), 5] <- formatC(ausgabe$Range.Phi[2],
        digits = 2, format = "f")
    output[(1 + k.beta + k.alpha + k.gamma + 10), 2] <- "Range Omega"
    output[(1 + k.beta + k.alpha + k.gamma + 10), 4] <- formatC(ausgabe$Range.Omega[1],
        digits = 2, format = "f")
    output[(1 + k.beta + k.alpha + k.gamma + 10), 5] <- formatC(ausgabe$Range.Omega[2],
        digits = 2, format = "f")
    output2 <- data.frame(output)
    colnames(output2)[1] <- "#1"
    for (i in 1:dim(output2)[1]) {
        rownames(output2)[i] <- paste("#", i, sep = "")
    }
    return(output2)
}

