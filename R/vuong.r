vuong <- 
function (model1, model2, alpha=0.05, correction=T)
{
    ll1 <- model1$Loglike.Vuong
    ll2 <- model2$Loglike.Vuong
    p <- length(model1$Coefficients.Mu)+length(model1$Coefficients.Phi)+length(model1$Coefficients.Omega)
    q <- length(model2$Coefficients.Mu)+length(model2$Coefficients.Phi)+length(model2$Coefficients.Omega)
    n <- length(ll1)
    m.i <- ll1 - ll2 - ifelse(correction,p/2*log(n)-q/2*log(n),0)/n
    nu <- (sqrt(n) * mean(m.i))/(sqrt((n - 1)/n * var(m.i)))
    print(paste("nu =", formatC(nu, 3, format = "g")))
    if (abs(nu) < qnorm(1-alpha/2)) {
        print("None of the models is favoured.")
    }
    if (nu >= qnorm(1-alpha/2)) {
        print("Favour model 1.")
    }
    if (nu <= -qnorm(1-alpha/2)) {
        print("Favour model 2.")
    }
    print(paste("P-value:", formatC(2 * pnorm(-abs(nu)), 2, format = "f")))
    rm(ll1, ll2)
}
