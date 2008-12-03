"fit.zigp1" <-
function(delta)

{

        n <- get("n", pos=globalenv())
        k <- get("ksave", pos=globalenv())
        Xsave <- get("Xsave", pos=globalenv())
        t.i <- get("t.i", pos=globalenv())
        
        eta <- double(n)

        mu <- double(n)

        fit <- double(n)

        omega <- exp(delta[2])/(1+exp(delta[2]))



        if(k == 1) {

                        eta <- delta[3] * Xsave

                }

                else {

                        beta<-delta[3:(k+2)]

                        eta <- Xsave%*%beta

                        }



                mu <- t.i* exp(eta)

                fit <- (1 - omega) * mu



        return(list(fit = fit, mu = mu))

}

