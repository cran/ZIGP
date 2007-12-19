"summaryzigp1" <-
function(mledaten)

{

        X <- mledaten$Design

        if(is.matrix(X)) {

                n <- dim(X)[1]

                k <- dim(X)[2]

        }

        else {

                n <- length(X)

                k <- 1

        }

        out0 <- matrix(double(1), 1, 1)

        out0[1, 1] <- mledaten$ZI.Parameter

        colnames(out0) <- c("")

        rownames(out0) <- c("ZI-Parameter:        ")

        print(out0)

        out1 <- matrix(double(k), 1, k)

        kopf <- factor(rep("", k))

        for(i in 1:k) {

                 out1[1, i] <- mledaten$Coefficients[i]

        }

        colnames(out1) <- kopf

        rownames(out1) <- c("Coefficients:        ")

        print(out1)

        out2 <- matrix(double(1), 1, 1)

        out2[1, 1] <- mledaten$Dispersion.Parameter

        colnames(out2) <- c("")

        rownames(out2) <- c("Dispersion Parameter:")

        print(out2)

        out3 <- matrix(double(1), 1, 1)

        out3[1, 1] <- mledaten$RSS

        colnames(out3) <- c("")

        rownames(out3) <- c("RSS:                 ")

        print(out3)

        out4 <- matrix(double(1), 1, 1)

        out4[1, 1] <- mledaten$AIC

        colnames(out4) <- c("")

        rownames(out4) <- c("AIC:                 ")

        print(out4)

        out5 <- matrix(double(2), 1, 2)

        out5[1, 1] <- mledaten$Range.mu[1]

        out5[1, 2] <- mledaten$Range.mu[2]

        colnames(out5) <- c("", "")

        rownames(out5) <- c("Range of mu:         ")

        print(out5)

        out6 <- matrix(double(1), 1, 1)

        if (is.null(mledaten$Message)) { out6[1, 1] <- "NULL"}

        else {out6[1, 1] <- mledaten$Message}

        colnames(out6) <- c("")

        rownames(out6) <- c("Message:             ")

        print(out6)

        return()

}

