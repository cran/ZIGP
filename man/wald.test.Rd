\name{wald.test}
\alias{wald.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting ZIGP(mu(i), phi(i), omega(i)) - Regression Models }
\description{
'wald.test' is used to fit ZIGP(mu(i), phi(i), omega(i)) - Regression Models.
}
\usage{
wald.test(Yin, Xin, Win, Zin, Offset = rep(1, length(Yin)), init = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yin}{ response vector of length n. }
  \item{Xin}{ design matrix of dim (n x p) for mean modelling. }
  \item{Win}{ design matrix of dim (n x r) for overdispersion modelling. }
  \item{Zin}{ design matrix of dim (n x q) for zero inflation modelling. }
  \item{Offset}{ exposure for individual observation lengths. Defaults to a vector of 1.
          The offset MUST NOT be in 'log' scale. }
  \item{init}{ a logical value indicating whether initial optimization values for 
          dispersion are set to -2.5 and values for zero inflation regression 
          parameters are set to -1 (init = F) or are estimated by a 
          ZIGP(mu(i), phi, omega)-model (init = T).  Defaults to 'T'. }
}
\details{
     In order to include an intercept in a design matrix, one has to
     add a vector of ones to the design matrix:
     'Intercept <- rep(1,n)'.
     
     If the output should have variable names additionally to parameter 
     tokens (such as 'b0', 'a0' or 'g0'), create the design matrix by
     'W <- cbind(Intercept, gender, height)'.
}
\examples{
## Number of damages in car insurance.
     
damage <- c(0,1,0,0,0,4,2,0,1,0,1,1,0,2,0,0,1,0,0,1,0,0,0)
Intercept <- rep(1,length(damage))
insurance.year <- c(1,1.2,0.8,1,2,1,1.1,1,1,1.1,1.2,1.3,0.9,1.4,1,1,1,1.2,
1,1,1,1,1)
drivers.age <- c(25,19,30,48,30,18,19,29,24,54,56,20,38,18,23,58,
47,36,25,28,38,39,42)
# for overdispersion: car brand dummy in {1,2,3}, brand = 1 is reference
brand <- c(1,2,1,3,3,2,2,1,1,3,2,2,1,3,1,3,2,2,1,1,3,3,2)
brand2 <- ifelse(brand==2,1,0)
brand3 <- ifelse(brand==3,1,0)
W <- cbind(brand2,brand3)
# abroad: driver has been abroad for longer time (=1)
abroad <- c(0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,1)
Y <- damage
X <- cbind(Intercept, drivers.age)
Z <- cbind(abroad)

wald.test(Yin=Y, Xin=X, Win=W, Zin=Z, Offset = insurance.year, init = FALSE)     

#1                            Estimate Std. Error   z value Pr(>|z|)      
#2           MU REGRESSION                                                
#3  b0           Intercept     1.47148    1.07377   1.37038  0.17057      
#4  b1         drivers.age    -0.05075    0.03907  -1.29897  0.19395      
#5          PHI REGRESSION                                                
#6  a0              brand2    -8.64637 2132.15915  -0.00406  0.99676      
#7  a1              brand3     0.17339    1.50296   0.11536  0.90816      
#8        OMEGA REGRESSION                                                
#9  g0              abroad    -1.10339    2.46771  -0.44713  0.65478      
#10                                                                       
#
#11       Signif. codes: 0 `***' 0.001 `**'  0.01 `*'  0.05 `.'  0.1 ` ' 1
#12             Iterations                     43                         
#13         Log Likelihood                  -23.4                         
#14    Pearson Chi Squared                    2.5                         
#15                    AIC                     57                         
#16               Range Mu                   0.23      2.45               
#17              Range Phi                   1.00      2.19               
#18            Range Omega                   0.25      0.50               

## The function is currently defined as
function(wald.test)

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
}
\keyword{ models }% at least one, from doc/KEYWORDS

