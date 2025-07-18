require("fitdistrplus")

# ?denscomp
visualize <- FALSE # TRUE for manual test with visualization of plots
nsample <- 1000
nsample <- 10

# (1) Plot various distributions fitted to serving size data
#
data(groundbeef)
serving <- groundbeef$serving
fitW <- fitdist(serving,"weibull")
fitln <- fitdist(serving,"lnorm")
fitg <- fitdist(serving,"gamma")

#sanity checks
try(denscomp("list(fitW, fitln, fitg)",horizontals = FALSE), silent=TRUE)
try(denscomp(list(fitW, fitln, fitg, a=1),horizontals = FALSE), silent=TRUE)

#real call
res <- denscomp(list(fitW, fitln, fitg), probability = TRUE)
str(res)

denscomp(list(fitW, fitln, fitg), probability = FALSE)
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  
  res <- denscomp(list(fitW, fitln, fitg), probability=TRUE, plotstyle = "ggplot")
  str(res)
  
  
  denscomp(list(fitW, fitln, fitg), probability=FALSE, plotstyle = "ggplot")
}

#test ylim argument
denscomp(list(fitW, fitln, fitg), probability=TRUE, ylim=c(0, .05))
denscomp(list(fitW, fitln, fitg), probability=FALSE, ylim=c(0, 100))
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize) {
  denscomp(list(fitW, fitln, fitg), probability=TRUE, ylim=c(0, .05), plotstyle = "ggplot")
  denscomp(list(fitW, fitln, fitg), probability=FALSE, ylim=c(0, 100), plotstyle = "ggplot")
}

#test xlim, legend, main, demp
denscomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
         main="ground beef fits",xlab="serving sizes (g)",
         ylab="F",xlim = c(0,250), xlegend = "topright", demp=TRUE)
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize) {
  denscomp(list(fitW, fitln, fitg), legendtext=c("Weibull","lognormal","gamma"),
           main="ground beef fits",xlab="serving sizes (g)",
           ylab="F",xlim = c(0,250), xlegend = "topright", demp=TRUE, plotstyle = "ggplot")
}


# (2) Plot lognormal distributions fitted by
# maximum goodness-of-fit estimation
# using various distances (data plotted in log scale)
#
data(endosulfan)
ATV <-subset(endosulfan, group == "NonArthroInvert")$ATV
flnMGEKS <- fitdist(ATV,"lnorm",method="mge",gof="KS")
flnMGEAD <- fitdist(ATV,"lnorm",method="mge",gof="AD")
flnMGEADL <- fitdist(ATV,"lnorm",method="mge",gof="ADL")
flnMGEAD2L <- fitdist(ATV,"lnorm",method="mge",gof="AD2L")
llfit <- list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L)

denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
         main="fits of a lognormal dist. using various GOF dist.")
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize) {
  denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
           main="fits of a lognormal dist. using various GOF dist.", plotstyle = "ggplot")
}

denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
         main="fits of a lognormal dist. using various GOF dist.",
         legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"))
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize) {
  denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
           main="fits of a lognormal dist. using various GOF dist.",
           legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"), plotstyle = "ggplot")
}

denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
         main="fits of a lognormal dist. using various GOF dist.",
         legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"),
         fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"))
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize) {
  denscomp(list(flnMGEKS, flnMGEAD, flnMGEADL, flnMGEAD2L),
           main="fits of a lognormal dist. using various GOF dist.",
           legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"),
           fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"), plotstyle = "ggplot")
}

denscomp(llfit,
         main="fits of a lognormal dist. using various GOF dist.",
         legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"),
         fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"),
         datacol="grey")
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize) {
  denscomp(llfit,
           main="fits of a lognormal dist. using various GOF dist.",
           legendtext=c("MGE KS","MGE AD","MGE ADL","MGE AD2L"),
           fitcol=c("black", "darkgreen", "yellowgreen", "yellow2"),
           datacol="grey", plotstyle = "ggplot")
}

denscomp(flnMGEKS, xlim=c(10,100000))
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize)
  denscomp(flnMGEKS, xlim=c(10,100000), plotstyle = "ggplot")


# (3)
#
#

if (visualize)
{
  x1 <- c(6.4,13.3,4.1,1.3,14.1,10.6,9.9,9.6,15.3,22.1,13.4,13.2,8.4,6.3,8.9,5.2,10.9,14.4)
  x <- seq(0, 1.1*max(x1), length=100)
  
  dgumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
  pgumbel <- function(x,a,b) exp(-exp((a-x)/b))
  
  f1 <- mledist(x1,"norm")
  f2 <- mledist(x1,"gumbel", start = list(a = 10, b = 5))
  f3 <- mledist(x1, "exp")
  
  # graph 1
  hist(x1, 10, prob=TRUE)
  lines(x, dnorm(x, f1$estimate[1], f1$estimate[2]), col="red")
  lines(x, dgumbel(x, f2$estimate[1], f2$estimate[2]), col="green")
  lines(x, dexp(x, f3$estimate[1]), col="blue")
  legend("topright", lty=1, leg = c("Normal", "Gumbel", "Exp"), col = c("red", "green", "blue"))
  
  # graph 2
  f1 <- fitdist(x1,"norm")
  f2 <- fitdist(x1,"gumbel", start = list(a = 10, b = 5)) # warning: The pgumbel function should have its first argument named: q as in base R
  f3 <- fitdist(x1, "exp")
  denscomp(list(f1, f2, f3), xlim = c(0, 30), fitlty = 1, legendtext = c("Normal","Gumbel","Exp"))
  
  # graph 3
  if (requireNamespace ("ggplot2", quietly = TRUE))
    denscomp(list(f1, f2, f3), xlim = c(0, 30), fitlty = 1, legendtext = c("Normal","Gumbel","Exp"), breaks = 12, plotstyle = "ggplot")
  
  
}

# (4) normal mixture
#
if (visualize)
{
  
  #mixture of two normal distributions
  #density
  dnorm2 <- function(x, poid, m1, s1, m2, s2)
    poid*dnorm(x, m1, s1) + (1-poid)*dnorm(x, m2, s2)
  #numerical approximate quantile function
  qnorm2 <- function(p, poid, m1, s1, m2, s2)
  {
    L2 <- function(x, prob)
      (prob - pnorm2(x, poid, m1, s1, m2, s2))^2
    sapply(p, function(pr) optimize(L2, c(-1000, 1000), prob=pr)$minimum)
  }
  #distribution function
  pnorm2 <- function(q, poid, m1, s1, m2, s2)
    poid*pnorm(q, m1, s1) + (1-poid)*pnorm(q, m2, s2)
  
  
  #basic normal distribution
  set.seed(1234)
  x2 <- c(rnorm(nsample, 5),  rnorm(nsample, 10))
  #MLE fit
  fit1 <- fitdist(x2, "norm2", "mle", start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2),
                  lower=c(0, 0, 0, 0, 0))
  fit2 <- fitdist(x2, "norm2", "qme", probs=c(1/6, 1/4, 1/3, 1/2, 2/3),
                  start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2),
                  lower=c(0, 0, 0, 0, 0), upper=c(1/2, Inf, Inf, Inf, Inf))
  fit3 <- fitdist(x2, "norm2", "mge", gof="AD",
                  start=list(poid=1/3, m1=4, s1=2, m2=8, s2=2),
                  lower=c(0, 0, 0, 0, 0), upper=c(1/2, Inf, Inf, Inf, Inf))
  
  denscomp(list(fit1, fit2, fit3))
  if (requireNamespace ("ggplot2", quietly = TRUE) & visualize)
    denscomp(list(fit1, fit2, fit3), plotstyle = "ggplot")
}

# (5) large data
#
if (visualize)
{
  n <- 1e4
  x <- rnorm(n)
  f <- fitdist(x, "norm")
  
  denscomp(f)
  denscomp(f, demp=TRUE)
  if (requireNamespace ("ggplot2", quietly = TRUE)) {
    denscomp(f, plotstyle = "ggplot")
    denscomp(f, demp=TRUE, plotstyle = "ggplot")
  }
}


# (6) graphical parameters
#
if (visualize)
{
  # 'graphics' plot style
  denscomp(list(fit1, fit2, fit3), plotstyle = "gr")
  denscomp(list(fit1, fit2, fit3), main = "Fitted distribution")
  denscomp(list(fit1, fit2, fit3), main = "Fitted distribution", addlegend = F, demp = T, dempcol = "purple")
  
  # 'ggplot' plot style
  if (requireNamespace ("ggplot2", quietly = TRUE)) {
    denscomp(list(fit1, fit2, fit3), plotstyle = "gg")
    denscomp(list(fit1, fit2, fit3), plotstyle = "ggplot", breaks = 20, pro = F)
    dcomp <- denscomp(list(fit1, fit2, fit3), plotstyle = "gg", demp = T)
    dcomp + ggplot2::theme_minimal() + ggplot2::ggtitle("Histogram and\ntheoretical densities")
    dcomp + ggplot2::guides(colour = ggplot2::guide_legend("Fitted distribution"), linetype = ggplot2::guide_legend("Fitted distribution")) 
  }
}



# (7) test legend labels
#
if (visualize)
{
  serving <- groundbeef$serving
  fitW <- fitdist(serving,"weibull")
  fitW2 <- fitdist(serving,"weibull", method="qme", probs=c(1/3,2/3))
  fitW3 <- fitdist(serving,"weibull", method="qme", probs=c(1/2,2/3))
  fitln <- fitdist(serving,"lnorm")
  fitg <- fitdist(serving,"gamma")
  
  denscomp(list(fitW, fitln, fitg)) #distrib
  denscomp(list(fitW, fitW2, fitln, fitg)) #distrib+method
  denscomp(list(fitW, fitW2, fitW3, fitln, fitg)) #distrib+method+num
  if (requireNamespace ("ggplot2", quietly = TRUE))
    denscomp(list(fitW, fitW2, fitW3, fitln, fitg), plotstyle = "ggplot") #distrib+method+num
}

# (8) discrete distrib
#

x <- c(rpois(nsample, 5), rbinom(nsample, 12, 2/3))
fpois <- fitdist(x, "pois")
fgeo <- fitdist(x, "geom")
fnbinom <- fitdist(x, "nbinom") 
# Messages d'avis :Messages d'avis :
# 1: Dans cov2cor(varcovar) :
#   diag(V) had non-positive or NA entries; the non-finite result may be dubious
# 2: Dans sqrt(diag(varcovar)) : Production de NaN

par(mar=c(4,4,2,1))
denscomp(list(fpois, fnbinom, fgeo), probability = TRUE)
denscomp(list(fpois, fnbinom, fgeo), probability = FALSE)
denscomp(list(fpois, fnbinom, fgeo), fittype="o")
denscomp(list(fpois, fnbinom, fgeo), fittype="p")
# 'ggplot' plot style
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize) {
  denscomp(list(fpois, fnbinom, fgeo), plotstyle="ggplot", probability = TRUE)
  denscomp(list(fpois, fnbinom, fgeo), plotstyle="ggplot", probability = FALSE)
  denscomp(list(fpois, fnbinom, fgeo), fittype="o", plotstyle="ggplot")
  denscomp(list(fpois, fnbinom, fgeo), fittype="p", plotstyle="ggplot")
}

# test the call to any()
fpois$discrete <- fnbinom$discrete <- FALSE
denscomp(list(fpois, fnbinom, fgeo))
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize)
  denscomp(list(fpois, fnbinom, fgeo), plotstyle="ggplot")
#test the forced usage
fgeo$discrete <- FALSE
denscomp(list(fpois, fnbinom, fgeo), discrete=TRUE)
if (requireNamespace ("ggplot2", quietly = TRUE) & visualize)
  denscomp(list(fpois, fnbinom, fgeo), discrete=TRUE, plotstyle="ggplot")

if (visualize)
{
  x <- c(rpois(nsample, 30), rbinom(nsample, 12, 2/3))
  fpois <- fitdist(x, "pois")
  fgeo <- fitdist(x, "geom")
  fnbinom <- fitdist(x, "nbinom")
  
  #3 types of plot of probability mass function
  par(mar=c(4,4,2,1))
  denscomp(list(fpois, fnbinom, fgeo)) #fittype = "l"
  denscomp(list(fpois, fnbinom, fgeo), fittype = "p")
  denscomp(list(fpois, fnbinom, fgeo), fittype = "o")
  if (requireNamespace ("ggplot2", quietly = TRUE)) {
    denscomp(list(fpois, fnbinom, fgeo), plotstyle="ggplot") #fittype = "l"
    denscomp(list(fpois, fnbinom, fgeo), fittype = "p", plotstyle="ggplot")
    denscomp(list(fpois, fnbinom, fgeo), fittype = "o", plotstyle="ggplot")
  }
}

# (9) examples with user specified regular of irregular breaks in the histogram
#     in probability or not
#
if (visualize)
{
  # two plots with user specified regular breaks in probability or not
  # hist(serving, breaks = seq(0,200,50))
  denscomp(list(fitW, fitln, fitg), addlegend = FALSE,
           main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), 
           breaks = seq(0,200,50))
  denscomp(list(fitW, fitln, fitg), addlegend = FALSE,
           main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), 
           probability = FALSE, breaks = seq(0,200,50))
  # with ggplot2
  denscomp(list(fitW, fitln, fitg), addlegend = FALSE,
           main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), 
           plotstyle = "ggplot", breaks = seq(0,200,50))
  denscomp(list(fitW, fitln, fitg), addlegend = FALSE,
           main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), 
           probability = FALSE, plotstyle = "ggplot", breaks = seq(0,200,50))
  
  
  # two plots with irregular breaks in probability or not
  # hist(serving, breaks = c(0, 20, 50, 100, 200, 300))
  denscomp(list(fitW, fitln, fitg), addlegend = FALSE,
           main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), 
           breaks = c(0, 20, 50, 100, 200, 300))
  # hist(serving, breaks = c(0, 20, 50, 100, 200, 300), probability = FALSE)
  try(denscomp(list(fitW, fitln, fitg), addlegend = FALSE,
               main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), 
               breaks = c(0, 20, 50, 100, 200, 300), probability = FALSE))
  # with ggplot2
  denscomp(list(fitW, fitln, fitg), addlegend = FALSE,
           main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), 
           breaks = c(0, 20, 50, 100, 200, 300), plotstyle = "ggplot")
  ##### ggplot2 does not take into account non-equidistant breaks !!!!!!!!!!!!!!!!
  try(denscomp(list(fitW, fitln, fitg), addlegend = FALSE,
               main = "ground beef fits", xlab = "serving sizes (g)", xlim = c(0, 250), 
               breaks = c(0, 20, 50, 100, 200, 300), 
               probability = FALSE, plotstyle = "ggplot"))
}


# (10) fitlty, fitlwd for discrete
x <- c(rpois(nsample, 30), rbinom(nsample, 12, 2/3))
fpois <- fitdist(x, "pois")
fgeo <- fitdist(x, "geom")
fnbinom <- fitdist(x, "nbinom")
denscomp(list(fpois, fnbinom, fgeo), fitlty = 2, fitlwd = 3:1)
denscomp(list(fpois, fnbinom, fgeo), fittype = "o", fitlwd = 3:1)
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  denscomp(list(fpois, fnbinom, fgeo), plotstyle = "ggplot", fitlty = 2, fitlwd = 3:1)
  denscomp(list(fpois, fnbinom, fgeo), plotstyle = "ggplot", fittype = "o", fitlty = 1, fitlwd = 3:1)
}

# (11) fitlty, fitlwd for non discrete
denscomp(list(fitW, fitln, fitg), fitlty = 1, fitlwd = 3:1)
denscomp(list(fitW, fitln, fitg), fitlty = 1, fitlwd = 1:3, fitcol = c(1:2, 7))
if (requireNamespace ("ggplot2", quietly = TRUE)) {
  denscomp(list(fitW, fitln, fitg), plotstyle = "ggplot", fitlty = 1, fitlwd = 3:1)
  denscomp(list(fitW, fitln, fitg), plotstyle = "ggplot", fitlty = 1, fitlwd = 1:3, fitcol = c(1:2, 7))
}

