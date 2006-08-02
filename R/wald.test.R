"wald.test" <-
function(Yin, Xin, Win, Zin, Offset=rep(1,length(Yin)), init=T)

{

Y <<- Yin

X <<- Xin

W <<- Win

Z <<- Zin

k.beta <<- dim(X)[2]

k.alpha <<- dim(W)[2]

k.gamma <<- dim(Z)[2]

n <<- length(Y)



# Estimate coefficients

ausgabe <- mle.zigp(Y, X, W, Z, Offset = Offset, summary = FALSE, init=init)



hat.beta   <- ausgabe$Coefficients.Mu

hat.alpha  <- ausgabe$Coefficients.Phi

hat.gamma  <- ausgabe$Coefficients.Omega



# Compute square root of diagonal elements of FM^-1

B <- FM(hat.beta, hat.alpha, hat.gamma, X, W, Z, Offset = Offset)

sd.vector <- sqrt( diag(solve(B, tol = 1e-20)) )

hat.sd.beta   <- sd.vector[1:k.beta]

hat.sd.alpha  <- sd.vector[(k.beta+1):(k.beta+k.alpha)]

hat.sd.gamma  <- sd.vector[(k.beta+k.alpha+1):(k.beta+k.alpha+k.gamma)]



# Compute T-Statistics

z.stat.beta  <- hat.beta/hat.sd.beta

z.stat.alpha <- hat.alpha/hat.sd.alpha

z.stat.gamma <- hat.gamma/hat.sd.gamma



# Compute P-Values

p.value.beta  <- 2*pnorm(-abs(z.stat.beta ))

p.value.alpha <- 2*pnorm(-abs(z.stat.alpha))

p.value.gamma <- 2*pnorm(-abs(z.stat.gamma))



# Create ***

glimpse.beta  <- rep("",length(z.stat.beta))

glimpse.alpha <- rep("",length(z.stat.alpha))

glimpse.gamma <- rep("",length(z.stat.gamma))

for (i in 1:length(z.stat.beta)) {

if (p.value.beta[i] < 0.001) {glimpse.beta[i] <- "***"}

if (p.value.beta[i] >= 0.001 & p.value.beta[i] < 0.01) {glimpse.beta[i] <- "**"}

if (p.value.beta[i] >= 0.01 & p.value.beta[i] < 0.05) {glimpse.beta[i] <- "*"}

if (p.value.beta[i] >= 0.05 & p.value.beta[i] < 0.1) {glimpse.beta[i] <- "."} }

for (i in 1:length(z.stat.alpha)) {

if (p.value.alpha[i] < 0.001) {glimpse.alpha[i] <- "***"}

if (p.value.alpha[i] >= 0.001 & p.value.alpha[i] < 0.01) {glimpse.alpha[i] <- "**"}

if (p.value.alpha[i] >= 0.01 & p.value.alpha[i] < 0.05) {glimpse.alpha[i] <- "*"}

if (p.value.alpha[i] >= 0.05 & p.value.alpha[i] < 0.1) {glimpse.alpha[i] <- "."} }

for (i in 1:length(z.stat.gamma)) {

if (p.value.gamma[i] < 0.001) {glimpse.gamma[i] <- "***"}

if (p.value.gamma[i] >= 0.001 & p.value.gamma[i] < 0.01) {glimpse.gamma[i] <- "**"}

if (p.value.gamma[i] >= 0.01 & p.value.gamma[i] < 0.05) {glimpse.gamma[i] <- "*"}

if (p.value.gamma[i] >= 0.05 & p.value.gamma[i] < 0.1) {glimpse.gamma[i] <- "."} }



# Create output

coef.names.beta  <- paste("b",c(0:(k.beta-1)),sep="")

coef.names.alpha <- paste("a",c(0:(k.alpha-1)),sep="")

coef.names.gamma <- paste("g",c(0:(k.gamma-1)),sep="")

coef.desc.beta   <- double(k.beta)

coef.desc.alpha  <- double(k.alpha)

coef.desc.gamma  <- double(k.gamma)

for (i in 1:k.beta)  {

coef.desc.beta[i]  <- colnames(X, do.NULL=FALSE)[i]

if (is.matrix(X)) { if (max(X[,i])==1&min(X[,i])==1) {coef.desc.beta[i] <- "Intercept"} }

else{ if (max(X[i])==1&min(X[i])==1) {coef.desc.beta[i] <- "Intercept"} }

}

for (i in 1:k.alpha) {

coef.desc.alpha[i] <- colnames(W, do.NULL=FALSE)[i]

if (is.matrix(W)) { if (max(W[,i])==1&min(W[,i])==1) {coef.desc.alpha[i] <- "Intercept"} }

else{ if (max(W[i])==1&min(W[i])==1) {coef.desc.alpha[i] <- "Intercept"} }

}

for (i in 1:k.gamma) {

coef.desc.gamma[i] <- colnames(Z, do.NULL=FALSE)[i]

if (is.matrix(Z)) { if (max(Z[,i])==1&min(Z[,i])==1) {coef.desc.gamma[i] <- "Intercept"} }

else{ if (max(Z[i])==1&min(Z[i])==1) {coef.desc.gamma[i] <- "Intercept"} }

}



output <- matrix("",1+k.beta+k.alpha+k.gamma+12,7)

output[2:(1+k.beta+k.alpha+k.gamma+3),1] <- c("",coef.names.beta,"", coef.names.alpha,"", coef.names.gamma)

output[2:(1+k.beta+k.alpha+k.gamma+3),2] <- c("MU REGRESSION",coef.desc.beta,"PHI REGRESSION", coef.desc.alpha,"OMEGA REGRESSION", coef.desc.gamma)

output[1:(1+k.beta+k.alpha+k.gamma+3),3] <- c("Estimate","",formatC(hat.beta,5,format="f"),"",formatC(hat.alpha,5,format="f"),"",formatC(hat.gamma,5,format="f"))

output[1:(1+k.beta+k.alpha+k.gamma+3),4] <- c("Std. Error","",formatC(hat.sd.beta,5,format="f"),"",formatC(hat.sd.alpha,5,format="f"),"",formatC(hat.sd.gamma,5,format="f"))

output[1:(1+k.beta+k.alpha+k.gamma+3),5] <- c("z value","",formatC(z.stat.beta,5,format="f"),"",formatC(z.stat.alpha,5,format="f"),"",formatC(z.stat.gamma,5,format="f"))

output[1:(1+k.beta+k.alpha+k.gamma+3),6] <- c("Pr(>|z|)","",formatC(p.value.beta,5,format="f"),"",formatC(p.value.alpha,5,format="f"),"",formatC(p.value.gamma,5,format="f"))

output[1:(1+k.beta+k.alpha+k.gamma+3),7] <- c("","",glimpse.beta,"",glimpse.alpha,"",glimpse.gamma)

output[(1+k.beta+k.alpha+k.gamma+5),2] <- "Signif. codes: 0"

output[(1+k.beta+k.alpha+k.gamma+5),3] <- "`***' 0.001"

output[(1+k.beta+k.alpha+k.gamma+5),4] <- "`**'  0.01"

output[(1+k.beta+k.alpha+k.gamma+5),5] <- "`*'  0.05"

output[(1+k.beta+k.alpha+k.gamma+5),6] <- "`.'  0.1"

output[(1+k.beta+k.alpha+k.gamma+5),7] <- "` ' 1"

output[(1+k.beta+k.alpha+k.gamma+6),2] <- "Iterations"

output[(1+k.beta+k.alpha+k.gamma+6),4] <- ausgabe$Iterations[1]

output[(1+k.beta+k.alpha+k.gamma+7),2] <- "Log Likelihood"

output[(1+k.beta+k.alpha+k.gamma+7),4] <- formatC(ausgabe$Log.Likelihood,digits=1,format="f")

output[(1+k.beta+k.alpha+k.gamma+8),2] <- "Pearson Chi Squared"

output[(1+k.beta+k.alpha+k.gamma+8),4] <- formatC(ausgabe$Pearson,digits=1,format="f")

output[(1+k.beta+k.alpha+k.gamma+9),2] <- "AIC"

output[(1+k.beta+k.alpha+k.gamma+9),4] <- round(ausgabe$AIC)

output[(1+k.beta+k.alpha+k.gamma+10),2] <- "Range Mu"

output[(1+k.beta+k.alpha+k.gamma+10),4] <- formatC(ausgabe$Range.Mu[1],digits=2,format="f")

output[(1+k.beta+k.alpha+k.gamma+10),5] <- formatC(ausgabe$Range.Mu[2],digits=2,format="f")

output[(1+k.beta+k.alpha+k.gamma+11),2] <- "Range Phi"

output[(1+k.beta+k.alpha+k.gamma+11),4] <- formatC(ausgabe$Range.Phi[1],digits=2,format="f")

output[(1+k.beta+k.alpha+k.gamma+11),5] <- formatC(ausgabe$Range.Phi[2],digits=2,format="f")

output[(1+k.beta+k.alpha+k.gamma+12),2] <- "Range Omega"

output[(1+k.beta+k.alpha+k.gamma+12),4] <- formatC(ausgabe$Range.Omega[1],digits=2,format="f")

output[(1+k.beta+k.alpha+k.gamma+12),5] <- formatC(ausgabe$Range.Omega[2],digits=2,format="f")

output2 <- data.frame(output)

return(output2)

}

