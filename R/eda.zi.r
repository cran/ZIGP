eda.zi <- 
function(x, y, numberclasses=5) {
  Y <- y
  if (length(levels(factor(x)))>2) {
    probs.in <- (0:numberclasses)/numberclasses
    br.x <- quantile(x, probs = probs.in,type=1)

    for (prcounter in 2:numberclasses) {
      if (br.x[prcounter]==br.x[prcounter-1]) {
        qstart <- probs.in[prcounter]
        counter <- 0
        while (qstart+counter*0.01<.99 & quantile(x,qstart+counter*0.01,type=1)==quantile(x,probs.in[prcounter-1],type=1)) counter <- counter + 1
        probs.in[prcounter] <- qstart+counter*0.01
        probs.in[(prcounter+1):numberclasses] <- (1:(numberclasses-prcounter))*(1-probs.in[prcounter])/(numberclasses-prcounter+1)+probs.in[prcounter]
        br.x <- quantile(x, probs = probs.in,type=1)
      }
    }

    probs.in <- probs.in[rank(br.x)-floor(rank(br.x))==0]
    if (probs.in[length(probs.in)]<1) probs.in <- c(probs.in,1)

    xtemp <- x
    xtemp[x==min(x)] <- min(x)-0.00001
    br.x <- quantile(xtemp, probs = probs.in,type=1)

    x.cut <- cut(x, breaks = br.x,include.lowest = FALSE,right=TRUE)

    xlabel <- rep(0,length(br.x)-1)
    for (i in 1:length(br.x)-1) {
    xlabel[i] <- (br.x[i]+br.x[i+1])/2 }
  }
  if (length(levels(factor(x)))==2) {
    br.x <- c(-0.5,0.5,1.5)
    x.cut <- factor(ifelse(x==levels(factor(x))[1],0,1))
    xlabel <- as.double(levels(factor(x)))
  }

  br <- br.x
  no.levels <- length(levels(x.cut))
  label <- "x"

 
  temp <- table(Y, x.cut)
  n <- apply(temp, 2, sum) 
  Anzahl.Nullen <- rep(0,length(n))
  omega.hat <- Anzahl.Nullen
  if (rownames(temp)[1]=="0") {
    Anzahl.Nullen <- sum(temp[1, ])
    omega.hat <- (temp[1, ])/ Anzahl.Nullen * sum(n) / (n*length(levels(factor(x))))
  }
  linie <- log((1/dim(temp)[2]+0.5)/(1-1/dim(temp)[2]+0.5))
  logit <- log((omega.hat + 0.5) / (1-omega.hat + 0.5))

  plot(xlabel,logit,type="b",lty=2,ylab="(shifted) empirical logit", axes=F, xlab=label, lwd=3,
       xlim=c(min(xlabel),max(xlabel)),
       ylim=c(min(logit,linie),max(logit,linie)),
       col="blue")
  axis(2)
  axis(1, at = xlabel, labels = formatC(xlabel,3,format="g"))
  abline(h=linie,lty=3)
}
