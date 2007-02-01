\name{mle.zigp}
\alias{mle.zigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Likelihood Estimates }
\description{
'mle.zigp' is used to calculate the MLEs of the regression parameters for
mean, overdispersion and zero inflation.
}
\usage{
mle.zigp(Yin, Xin, Win, Zin, Offset = rep(1, length(Yin)), summary = TRUE, 
init = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yin}{ response vector of length n. }
  \item{Xin}{ design matrix of dim (n x p) for mean modelling. }
  \item{Win}{ design matrix of dim (n x r) for overdispersion modelling. }
  \item{Zin}{ design matrix of dim (n x q) for zero inflation modelling. }
  \item{Offset}{ exposure for individual observation lengths. Defaults to a vector of 1.
        The offset MUST NOT be in 'log' scale. }
  \item{summary}{ a logical value indicating whether a list summary is shown (=T) or returned
          to the user (=F), i.e. if values of the ML estimation should be extracted, 
          summary has to be set to 'F'. Defaults to 'T'. }
  \item{init}{ a logical value indicating whether initial optimization values for 
          dispersion are set to -2.5 and values for zero inflation regression 
          parameters are set to -1 (init = F) or are estimated by a 
          ZIGP(mu(i), phi, omega)-model (init = T).  Defaults to 'T'. }
}
\value{
If summary is set to 'F', the following values are returned:
     
  \item{ZI.Parameter        }{Zero Inflation estimate}
  \item{Coefficients.Mu     }{Regression Parameters estimates for mean}
  \item{Coefficients.Phi    }{Regression Par. estimates for overdispersion}
  \item{Coefficients.Omega  }{Regression Parameters estimates for ZI}
  \item{Dispersion.Parameter}{Dispersion Parameter Estimate}
  \item{Range.Mu            }{mean range}
  \item{Range.Phi           }{dispersion parameter range}
  \item{Range.Omega         }{ZI parameter range}
  \item{Log.Likelihood      }{log likelihood value at MLE}
  \item{Residuals           }{residuals}
  \item{Pearson             }{Pearson Chi Squared Statistics}
  \item{AIC                 }{Akaike Information Criterion}
  \item{Iterations          }{number of iterations needed}
  \item{Time                }{time needed}
  \item{Message             }{convergence message of underlying 'optim' routine}
  \item{Response            }{response vector}
  \item{Fitted.Values       }{fitted values}
  \item{Design.Mu           }{design matrix for mean}
  \item{Design.Phi          }{design matrix for overdispersion}
  \item{Design.Omega        }{design matrix for ZI}
}
\examples{
## Number of damages in car insurance.

damage <- c(0,1,0,0,0,4,2,0,1,0,1,1,0,2,0,0,1,0,0,1,0,0,0)
Intercept <- rep(1,length(damage))
insurance.year <- c(1,1.2,0.8,1,2,1,1.1,1,1,1.1,1.2,1.3,0.9,1.4,1,1,1,1.2,
1,1,1,1,1)
drivers.age <- c(25,19,30,48,30,18,19,29,24,54,56,20,38,18,23,58,
47,36,25,28,38,39,42)
n <- length(damage)
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

mle.zigp(Yin=Y, Xin=X, Win=W, Zin=Z, Offset = insurance.year, 
summary = TRUE, init = FALSE)     

#[1]   Range for ZI-Parameters:   0.2491062 0.5
#[2]   Range of Dispersion Pars:  1.000176 2.189325
#[3]   Coefficients for mu:       1.471478 -0.05075418
#[4]   Coefficients for phi:      -8.646371 0.1733860
#[5]   Coefficients for omega:    -1.103385
#[6]   Pearson Chi Squared:       15.09779
#[7]   Range of mu:               0.2294054 2.445806
#[8]   Message:                   "NULL"
#[9]   AIC:                       56.88305

## The function is currently defined as
function(Yin, Xin, Win, Zin, Offset = rep(1,length(Yin)), summary = TRUE, init = FALSE)

{

  Y <<- Yin

  X <<- Xin

  W <<- Win

  Z <<- Zin

  if(is.matrix(X)) {

    n <- dim(X)[1]

    k.beta <<- dim(X)[2]

  }

  else {

    n <- length(X)

    k.beta <<- 1

  }

  if(is.matrix(W)) {

    k.alpha <<- dim(W)[2]

    nw <<- dim(W)[1]

  }

  else {

    k.alpha <<- 1

    nw <<- length(W)

  }

  if(is.matrix(Z)) {

    k.gamma <<- dim(Z)[2]

    nz <<- dim(Z)[1]

  }

  else {

    k.gamma <<- 1

    nz <<- length(Z)

  }



  if(nz==n & nw==n){

  t.i <<- Offset



  if(init){   #optimized initial values

    start.delta <- optimized.run(Y,X,W,Z)



    # maximization of log likelihood

    start <- proc.time()[2]

    opt <- optim(par = start.delta, fn = loglikelihood.zigp, gr=gradient,

    method = "BFGS")

    #  method = "CG")

  }

  else

  {

    # simple initial values

    alpha <- rep(-2.5,k.alpha)

    gamma <- rep(-2.5,k.gamma)

    beta <- summary(glm(Y ~ offset(log(t.i)) + 1 + X, family = poisson(link=log)

                    ))$coefficients[, 1]

    start.delta <- c(beta,alpha,gamma)

    # maximization of log likelihood w.r.t. delta s.t. above limits hold

    start <- proc.time()[2]

    opt <- optim(par = start.delta, fn = loglikelihood.zigp, gr=gradient,

    method = "BFGS")

    #  method = "CG")

  }



  delta <- opt$par

  it <- opt$counts

  loglikelihood <-  - opt$value

  if(is.null(opt$message)) {message <- "NULL"}

  else { message <- opt$message }

  AIC <- -2 * loglikelihood + 2 * (k.beta + k.alpha + k.gamma)

  time.needed <- proc.time()[2] - start



# results

  coef <- double(k.beta + k.alpha + k.gamma)

  coef <- delta

  beta <- coef[1 : k.beta]

  alpha <- coef[(k.beta + 1) : (k.beta + k.alpha)]

  gamma <- coef[(k.beta + k.alpha + 1) : (k.beta + k.alpha + k.gamma)]



  fit <- fit.zigp(delta)

  res <- Y - fit$fit

  sum((Y - fit$fit)^2 / ((1-fit$omega)*fit$mu*(fit$phi^2+fit$mu*fit$omega)))



  range.mu <- c(min(fit$mu), max(fit$mu))

# compute zero inflation parameters omega

  if(is.matrix(Z)){ eta.omega <- Z \%*\% gamma }

  else { eta.omega <- Z * gamma }

  omega <- exp(eta.omega)/(1+exp(eta.omega))

  range.omega <- c(min(omega), max(omega))

# compute dispersion parameters phi

  if(is.matrix(W)){ eta.phi <- W \%*\% alpha }

  else { eta.phi <- W * alpha }

  phi <- 1+ exp(eta.phi)

  range.phi <- c(min(phi), max(phi))



# summarize results

  mle.data <- list(ZI.Parameter = omega,

    Coefficients.Mu = beta,

    Coefficients.Phi = alpha,

    Coefficients.Omega = gamma,

    Dispersion.Parameter = phi,

    Range.Mu = range.mu,

    Range.Phi = range.phi,

    Range.Omega = range.omega,

    Log.Likelihood = loglikelihood,

    Residuals = res,

    Pearson = chsq,

    AIC = AIC,

    Iterations = it,

    Time = time.needed,

    Message = message,

    Response = Y,

    Fitted.Values = fit$fit,

    Design.Mu = X,

    Design.Phi = W,

    Design.Omega = Z)

    if(summary) {

      summary.zigp(mle.data)

    }

    else {

      return(mle.data)

    }

  }

  else{

  return("X, W and Z have different numbers of rows!")

  }

  }
}
\keyword{ models }% at least one, from doc/KEYWORDS
