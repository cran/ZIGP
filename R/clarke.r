clarke <- 
function (model1, model2, alpha=0.05, correction=T)
{
    ll1 <- model1$Loglike.Vuong
    ll2 <- model2$Loglike.Vuong
    p <- length(model1$Coefficients.Mu)+length(model1$Coefficients.Phi)+length(model1$Coefficients.Omega)
    q <- length(model2$Coefficients.Mu)+length(model2$Coefficients.Phi)+length(model2$Coefficients.Omega)
    n <- length(ll1)
    m.i <- ll1 - ll2 - ifelse(correction,p/2*log(n)-q/2*log(n),0)/n

    B <- sum(m.i>0)
    if (B<n/2) {
      print("Perform lower tail test:")
      pval <- pbinom(B,n,.5)
      print("H0: 'Models equally good'   vs.   H1: 'Model 2 better than model 1'")
      print(ifelse(pval>alpha,"H0 cannot be rejected on level alpha.", "Reject H0 on level alpha."))
      print(paste("P value: ",round(pval,5),sep=""))
    }
    if (B>=n/2) {
      print("Perform upper tail test:")
      pval <- 1-pbinom(B-1,n,.5)
      print("H0: 'Models equally good'   vs.   H1: 'Model 1 better than model 2'")
      print(ifelse(pval>1-alpha,"H0 cannot be rejected on level alpha.", "Reject H0 on level alpha."))
      print(paste("P value: ",round(pval,5),sep=""))
    }
}
