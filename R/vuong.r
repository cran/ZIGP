vuong <- function(model1, model2) {
  ll1 <- model1$Loglike.Vuong
  ll2 <- model2$Loglike.Vuong

  m.i <- ll1 - ll2
  #test statistic
  n <- length(ll1)
  nu <- ( sqrt(n) * mean(m.i) ) / ( sqrt((n-1)/n*var(m.i)) )

  print(paste("nu =",formatC(nu,3,format="g")))
  if (abs(nu)<2) { print("None of the models is favoured.") }
  if (nu>=2)     { print("Favour model 1.") }
  if (nu<=-2)    { print("Favour model 2.") }
  print(paste("P-value:",formatC(2*pnorm(-abs(nu)),2,format="f")))
  rm(ll1,ll2)
}

