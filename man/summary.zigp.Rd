\name{summary.zigp}
\alias{summary.zigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summary }
\description{
'summary.zigp' summarizes the results of the maximum likelihood estimation.
}
\usage{
summary.zigp(mle.data.in)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mle.data.in}{ a list of result parameters given by function 'mle.zigp' }
}
\examples{
## The function is currently defined as
function(mle.data.in)

{

        X <- mle.data.in$Design.Mu

        if(is.matrix(X)) {

                nx <- dim(X)[1]

                k.beta <- dim(X)[2]

        }

        else {

                nx <- length(X)

                k.beta <- 1

        }

        W <- mle.data.in$Design.Phi

        if(is.matrix(W)) {

                nw <- dim(W)[1]

                k.alpha <- dim(W)[2]

        }

        else {

                nw <- length(W)

                k.alpha <- 1

        }

        Z <- mle.data.in$Design.Omega

        if(is.matrix(Z)) {

                nz <- dim(Z)[1]

                k.gamma <- dim(Z)[2]

        }

        else {

                nz <- length(Z)

                k.gamma <- 1

        }



        # range of ZI parameters

        out0 <- matrix(double(2), 1, 2)

        out0[1, 1] <- min(mle.data.in$ZI.Parameter)

        out0[1, 2] <- max(mle.data.in$ZI.Parameter)

        colnames(out0) <- c("", "")

        rownames(out0) <- c("Range for ZI-Parameters:  ")

        print(out0)

        # range of Dispersion Parameters

        out1 <- matrix(double(2), 1, 2)

        out1[1, 1] <- min(mle.data.in$Dispersion.Parameter)

        out1[1, 2] <- max(mle.data.in$Dispersion.Parameter)

        colnames(out1) <- c("", "")

        rownames(out1) <- c("Range of Dispersion Pars: ")

        print(out1)

        # beta

        out2 <- matrix(double(k.beta), 1, k.beta)

        kopf <- rep("", k.beta)

        for(i in 1:k.beta) {

                out2[1, i] <- mle.data.in$Coefficients.Mu[i]

        }

        colnames(out2) <- kopf

        rownames(out2) <- c("Coefficients for mu:      ")

        print(out2)

        # alpha

        out3 <- matrix(double(k.alpha), 1, k.alpha)

        kopf <- rep("", k.alpha)

        for(i in 1:k.alpha) {

                out3[1, i] <- mle.data.in$Coefficients.Phi[i]

        }

        colnames(out3) <- kopf

        rownames(out3) <- c("Coefficients for phi:     ")

        print(out3)

        # gamma

        out4 <- matrix(double(k.gamma), 1, k.gamma)

        kopf <- rep("", k.gamma)

        for(i in 1:k.gamma) {

                out4[1, i] <- mle.data.in$Coefficients.Omega[i]

        }

        colnames(out4) <- kopf

        rownames(out4) <- c("Coefficients for omega:   ")

        print(out4)

        # Pearson Chi Squared

        out5 <- matrix(double(1), 1, 1)

        out5[1, 1] <- mle.data.in$Pearson

        colnames(out5) <- c("")

        rownames(out5) <- c("Pearson Chi Squared:      ")

        print(out5)

        # range of mu

        out6 <- matrix(double(2), 1, 2)

        out6[1, 1] <- mle.data.in$Range.Mu[1]

        out6[1, 2] <- mle.data.in$Range.Mu[2]

        colnames(out6) <- c("", "")

        rownames(out6) <- c("Range of mu:              ")

        print(out6)

        # optimization message

        out7 <- matrix(double(1), 1, 1)

        out7[1, 1] <- mle.data.in$Message

        colnames(out7) <- c("")

        rownames(out7) <- c("Message:                  ")

        print(out7)

        # AIC

        out8 <- matrix(double(1), 1, 1)

        out8[1, 1] <- mle.data.in$AIC

        colnames(out8) <- c("")

        rownames(out8) <- c("AIC:                     ")

        print(out8)

        return()

  }
}
\keyword{ list }% at least one, from doc/KEYWORDS
