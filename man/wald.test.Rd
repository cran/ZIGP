\name{wald.test}
\alias{wald.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting ZIGP(mu(i), phi(i), omega(i)) - Regression Models }
\description{
'wald.test' is used to fit ZIGP(mu(i), phi(i), omega(i)) - Regression Models.
}
\usage{
wald.test(Yin, Xin, Win=NULL, Zin=NULL, Offset = rep(1, length(Yin)), init = T)
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
     Overall overdispersion and/or zero-inflation
     can be modelled using an Intercept design.
     Setting W to NULL corresponds to modelling a ZIP model.
     Setting Z to NULL corresponds to modelling a GP model.
     Setting W and Z to NULL corresponds to modelling a Poisson GLM.

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
#14    Pearson Chi Squared                   15.1                         
#15                    AIC                     57                         
#16               Range Mu                   0.23      2.45               
#17              Range Phi                   1.00      2.19               
#18            Range Omega                   0.25      0.50               


# approximate equivalence of Poisson-glm and ZIGP-package results
# glm uses IWLS, ZIGP uses numerical maximization of the log-likelihood
# (time series character of the data is neglected)

data(Seatbelts)
DriversKilled <- as.vector(Seatbelts[,1]) # will be response
kms <- as.vector(Seatbelts[,5])           # will be exposure
PetrolPrice <- as.vector(Seatbelts[,6])   # will be covariate 1
law <- as.vector(Seatbelts[,8])           # will be covariate 2

fmla <- DriversKilled ~ PetrolPrice + law
out.glm <- glm(fmla, family=poisson, offset=log(kms))
summary(out.glm)

X <- cbind(rep(1,length(DriversKilled)),PetrolPrice,law)
wald.test(DriversKilled, X, NULL, NULL, Offset = kms)
}
\keyword{ models }% at least one, from doc/KEYWORDS

