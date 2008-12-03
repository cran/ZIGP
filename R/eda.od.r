eda.od <- 
function(x, y, Offset=rep(1,length(y)), numberclasses=5) {
  Y <- y/Offset
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

    # eliminate classes with only 1 observation
    probs.in <- probs.in[c(TRUE,table(x.cut)>1)]
    probs.in[length(probs.in)] <- 1
    br.x <- quantile(xtemp, probs = probs.in,type=1)
    x.cut <- cut(x, breaks = br.x,include.lowest = FALSE,right=TRUE)

    # eliminate classes with only equal observations (variance would be 0)
    probs.in <- probs.in[c(TRUE,(apply(table(Y,x.cut)>0,2,sum)==1)==FALSE)]
    probs.in[length(probs.in)] <- 1
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

  output.means <- rep(0,no.levels)
  output.vars <- rep(0,no.levels)
  for (level in 1:no.levels) {
    output.means[level] <- mean(Y[x.cut==levels(x.cut)[level]])
    output.vars[level] <- var(Y[x.cut==levels(x.cut)[level]])
  }
  output.f <- log(ifelse(sqrt(output.vars/output.means)-1>0,sqrt(output.vars/output.means)-1,0.0001))

  plot(xlabel,output.f,type="b",lty=2,ylab="OD influence", axes=F, xlab=label, lwd=3,
       xlim=c(min(xlabel),max(xlabel)),
       ylim=c(min(output.f),max(output.f)),
       col="blue")
  axis(2)
  axis(1, at = xlabel, labels = formatC(xlabel,3,format="g"))

}

