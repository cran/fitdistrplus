
library(fitdistrplus)
nbboot <- 100


# (1) basic fit of a gamma distribution by maximum likelihood estimation
#
data(groundbeef)
serving <- groundbeef$serving
fitg <- fitdist(serving, "gamma")
summary(fitg)
plot(fitg)
cdfcomp(fitg, addlegend=FALSE)

#check names
names(fitdist(serving, "gamma", optim.method="Brent", lower=0, upper=10, fix.arg=list(shape=2))$estimate)
names(fitdist(serving, "gamma", optim.method="Nelder-Mead")$estimate)
names(fitdist(serving, "gamma", optim.method="BFGS")$estimate)
# names(fitdist(serving, "gamma", optim.method="CG", start=list(shape=4, rate=1/20))$estimate)
names(fitdist(serving, "gamma", optim.method="L-BFGS-B", lower=0)$estimate)


# (2) use the moment matching estimation (using a closed formula)
#

fitgmme <- fitdist(serving, "gamma", method="mme")
summary(fitgmme)
plot(fitgmme)
quantile(fitgmme)

# (3) fit and comparison of various fits
#
fitW <- fitdist(serving, "weibull")
fitg <- fitdist(serving, "gamma")
fitln <- fitdist(serving, "lnorm")
summary(fitW)
summary(fitg)
summary(fitln)
cdfcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))

# (4) defining your own distribution functions, here for the Gumbel distribution
# for other distributions, see the CRAN task view 
# dedicated to probability distributions
#
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

fitgumbel <- fitdist(serving, "gumbel", start=list(a=10, b=10))
summary(fitgumbel)
plot(fitgumbel)

# (5) fit discrete distributions (Poisson and negative binomial)
#

data(toxocara)
number <- toxocara$number
fitp <- fitdist(number,"pois")
summary(fitp)
plot(fitp)
fitnb <- fitdist(number,"nbinom")
summary(fitnb)
plot(fitnb)

# (6) how to change the optimisation method?
#

oldw <- getOption("warn")
options(warn=-1) #to avoid a high number of warnings when nwarnings > 50
fitdist(serving, "gamma", optim.method="Nelder-Mead")
fitdist(serving, "gamma", optim.method="BFGS") 
fitdist(serving, "gamma", optim.method="SANN")
options(warn=oldw)

# (7) custom optimization function
#

#create the sample
set.seed(1234)
mysample <- rexp(100, 5)
mystart <- list(rate=8)

res1 <- fitdist(mysample, dexp, start= mystart, optim.method="Nelder-Mead")

#show the result
summary(res1)

#the warning tell us to use optimise, because the Nelder-Mead is not adequate.

#to meet the standard 'fn' argument and specific name arguments, we wrap optimize, 
myoptimize <- function(fn, par, ...) 
{
    res <- optimize(f=fn, ..., maximum=FALSE)  
#assume the optimization function minimize
    
    standardres <- c(res, convergence=0, value=res$objective, 
                     par=res$minimum, hessian=NA)
    
    return(standardres)
}

#call fitdist with a 'custom' optimization function
res2 <- fitdist(mysample, dexp, start=mystart, custom.optim=myoptimize, 
interval=c(0, 100))

#show the result
summary(res2)


# (8) custom optimization function - another example with the genetic algorithm
#
if(any(installed.packages()[,"Package"] == "rgenoud"))
{

#set a sample
    fit1 <- fitdist(serving, "gamma")
    summary(fit1)
    
#wrap genoud function rgenoud package
    mygenoud <- function(fn, par, ...) 
    {
        require(rgenoud)
        res <- genoud(fn, starting.values=par, ...)        
        standardres <- c(res, convergence=0, counts=NULL)
        
        return(standardres)
    }
    
#call fitdist with a 'custom' optimization function
    fit2 <- fitdist(serving, "gamma", custom.optim=mygenoud, nvars=2, start=as.list(fit1$estimate),
                    Domains=cbind(c(0, 0), c(10, 10)), boundary.enforcement=1, 
                    print.level=0, hessian=TRUE)
    
    summary(fit2)
}

# (9) estimation of the standard deviation of a gamma distribution 
# by maximum likelihood with the shape fixed at 4 using the argument fix.arg
#
f1c  <- fitdist(serving,"gamma",start=list(rate=0.1),fix.arg=list(shape=4))
summary(f1c)
plot(f1c)


# (10) fit of a Weibull distribution to serving size data 
# by maximum likelihood estimation
# or by quantile matching estimation (in this example 
# matching first and third quartiles)
#
fWmle <- fitdist(serving, "weibull")
summary(fWmle)
plot(fWmle)
gofstat(fWmle)

fWqme <- fitdist(serving, "weibull", method="qme", probs=c(0.25, 0.75))
summary(fWqme)
plot(fWqme)
gofstat(fWqme)


# (11) Fit of a Pareto distribution by numerical moment matching estimation
#
if(any(installed.packages()[,"Package"] == "actuar"))
{

    require(actuar)
    #simulate a sample
    set.seed(1234)
    x4 <- rpareto(1000, 6, 2)
    
    #empirical raw moment
    memp <- function(x, order)
        ifelse(order == 1, mean(x), sum(x^order)/length(x))
    
    #fit
    fP <- fitdist(x4, "pareto", method="mme", order=c(1, 2), memp="memp", 
                  start=list(shape=10, scale=10), lower=1, upper=Inf)
    summary(fP)
    plot(fP)
    
}

# (12) Fit of a Weibull distribution to serving size data by maximum 
# goodness-of-fit estimation using all the distances available
# 

data(groundbeef)
serving <- groundbeef$serving
(f1 <- fitdist(serving, "weibull", method="mge", gof="CvM"))
(f2 <- fitdist(serving, "weibull", method="mge", gof="KS"))
(f3 <- fitdist(serving, "weibull", method="mge", gof="AD"))
(f4 <- fitdist(serving, "weibull", method="mge", gof="ADR"))
(f5 <- fitdist(serving, "weibull", method="mge", gof="ADL"))
(f6 <- fitdist(serving, "weibull", method="mge", gof="AD2R"))
(f7 <- fitdist(serving, "weibull", method="mge", gof="AD2L"))
(f8 <- fitdist(serving, "weibull", method="mge", gof="AD2"))
cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8))
cdfcomp(list(f1, f2, f3, f4, f5, f6, f7, f8), xlogscale=TRUE, xlim=c(8, 250), verticals=TRUE)

# (13) Fit of a uniform distribution using maximum likelihood 
# (a closed formula is used in this special case where the loglikelihood is not defined),
# or maximum goodness-of-fit with Cramer-von Mises or Kolmogorov-Smirnov distance
# 

set.seed(1234)
u <- runif(50, min=5, max=10)

fumle <- fitdist(u, "unif", method="mle")
summary(fumle)
plot(fumle)
gofstat(fumle)

fuCvM <- fitdist(u, "unif", method="mge", gof="CvM")
summary(fuCvM)
plot(fuCvM)
gofstat(fuCvM)

fuKS <- fitdist(u, "unif", method="mge", gof="KS")
summary(fuKS)
plot(fuKS)
gofstat(fuKS)

# (14) scaling problem - too small values
#

x2 <- c(-0.00707717, -0.000947418, -0.00189753,
-0.000474947, -0.00190205, -0.000476077, 0.00237812, 0.000949668,
0.000474496, 0.00284226, -0.000473149, -0.000473373, 0, 0, 0.00283688,
-0.0037843, -0.0047506, -0.00238379, -0.00286807, 0.000478583,
0.000478354, -0.00143575, 0.00143575, 0.00238835, 0.0042847,
0.00237248, -0.00142281, -0.00142484, 0, 0.00142484, 0.000948767,
0.00378609, -0.000472478, 0.000472478, -0.0014181, 0, -0.000946522,
-0.00284495, 0, 0.00331832, 0.00283554, 0.00141476, -0.00141476,
-0.00188947, 0.00141743, -0.00236351, 0.00236351, 0.00235794,
0.00235239, -0.000940292, -0.0014121, -0.00283019, 0.000472255,
0.000472032, 0.000471809, -0.0014161, 0.0014161, -0.000943842,
0.000472032, -0.000944287, -0.00094518, -0.00189304, -0.000473821,
-0.000474046, 0.00331361, -0.000472701, -0.000946074, 0.00141878,
-0.000945627, -0.00189394, -0.00189753, -0.0057143, -0.00143369,
-0.00383326, 0.00143919, 0.000479272, -0.00191847, -0.000480192,
0.000960154, 0.000479731, 0, 0.000479501, 0.000958313, -0.00383878,
-0.00240674, 0.000963391, 0.000962464, -0.00192586, 0.000481812,
-0.00241138, -0.00144963)

for(i in 6:0)
{
	cat("\nscaling", 10^i, "\n")
	res <- try(mledist(x2*10^i, "cauchy"), silent=TRUE)
	if(class(res) == "try-error")
    print(res)
	else
	{
        cat("estimate\n")
        print(res$estimate)
        cat("Hessian\n")
        print(res$hessian)
    }
}



# (15) scaling problem - too big values
#

x1 <- c(
1401928684, 1413455609, 1432458425, 1436910475, 1494883250, 1565770323, 1577486458,
1568908053, 1606424896, 1632264979, 1780495643, 1865525923, 2035689865, 2141429306,
2335443964, 2465661689, 2563368221, 2845012431, 2949890881, 3180645942, 3309009836,
3618581152, 4109197451, 4064662257, 4028375795, 4176781983, 4303024833, 4493470109
)


for(i in 0:5)
{
	cat("\nscaling", 10^(-2*i), "\n")
	res <- mledist(x1*10^(-2*i), "norm")
	Hm1 <- try(solve(res$hessian), silent=TRUE)
	if(class(Hm1) == "try-error")
		print(Hm1)
	else
	{
	    cat("estimate\n")
    	print(res$estimate)
	    cat("Hessian\n")
	    print(res$hessian)
    	cat("inverse Hessian\n")
	    print(Hm1)
    }
}

 fitdist(x1, "norm")
 fitdist(x1*1e-6, "norm")




# (16) Fit of a lognormal distribution on acute toxicity values of endosulfan for
# nonarthropod invertebrates, using maximum likelihood estimation
# to estimate what is called a species sensitivity distribution 
# (SSD) in ecotoxicology, followed by estimation of the 5 percent quantile value of 
# the fitted distribution, what is called the 5 percent hazardous concentration (HC5)
# in ecotoxicology, with its two-sided 95 percent confidence interval calculated by 
# parametric bootstrap
#
data(endosulfan)
ATV <- subset(endosulfan, group == "NonArthroInvert")$ATV
log10ATV <- log10(subset(endosulfan, group == "NonArthroInvert")$ATV)
fln <- fitdist(log10ATV, "norm")

quantile(fln, probs = 0.05)


# (17) Fit of a triangular distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 
if(any(installed.packages()[,"Package"] == "mc2d"))
{
set.seed(1234)
require(mc2d)
t <- rtriang(100,min=5,mode=6,max=10)
fCvM <- fitdist(t,"triang",method="mge",start = list(min=4, mode=6,max=9),gof="CvM")
fKS <- fitdist(t,"triang",method="mge",start = list(min=4, mode=6,max=9),gof="KS")
cdfcomp(list(fCvM,fKS))
}

# (18) gumbel distribution
#
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

data(danishuni)
fitdist(danishuni$Loss, "gumbel", start=list(a=5, b=10))


# (19) check the 'start' argument
#

#create the sample
mysample <- rexp(100, 5)
mystart2 <- list(rate2=8)
mystart3 <- list(8)

try( fitdist(mysample, dexp, start= mystart2, method="mle") ) 
try( fitdist(mysample, dexp, start= mystart3, method="mle") ) 

try( fitdist(mysample, dexp, start= mystart2, method="mme") ) 
try( fitdist(mysample, dexp, start= mystart3, method="mme") ) 

try( fitdist(mysample, dexp, start= mystart2, method="qme", probs=1/2) ) 
try( fitdist(mysample, dexp, start= mystart3, method="qme", probs=1/2) ) 

try( fitdist(mysample, dexp, start= mystart2, method="mge", gof="AD") ) 
try( fitdist(mysample, dexp, start= mystart3, method="mge", gof="AD") ) 

# (20) example with dexgauss
# would require to suggest the package gamlss.dist in the Description file
#
#if(any(installed.packages()[,"Package"] == "gamlss.dist"))
#{
#    require(gamlss.dist)
#    set.seed(1234)
#    a=rexGAUS(100,mu=500,sigma=50,nu=75)
#    fitdist(a,dexGAUS,start=list(mu=median(a),sigma=sqrt(var(a)/2),nu=sqrt(var(a)/2)))
#}


# (21) check the 'keepdata' argument
#

#create the sample
x <- rexp(1e6, 5)
summary(x)
f1 <- fitdist(x, "exp", keepdata=FALSE)
f2 <- fitdist(x, "exp", keepdata=TRUE)

par(mfrow=c(1,2))
cdfcomp(f1)
cdfcomp(f2)

# (22) relevant example for zero modified geometric distribution
#
dzmgeom <- function(x, p1, p2)
{
  p1 * (x == 0) + (1-p1)*dgeom(x-1, p2)
}
pzmgeom <- function(q, p1, p2)
{
  p1 * (q >= 0) + (1-p1)*pgeom(q-1, p2)
}
rzmgeom <- function(n, p1, p2)
{
  u <- rbinom(n, 1, 1-p1) #prob to get zero is p1
  u[u != 0] <- rgeom(sum(u != 0), p2)+1
  u
}

x2 <- rzmgeom(1000, 1/2, 1/10)
table(x2)
#this is the MLE which converges almost surely and in distribution to the true value.
initp1 <- function(x) list(p1=mean(x == 0))

fitdist(x2, "zmgeom", start=list(p1=1/2, p2=1/2))

f2 <- fitdist(x2, "zmgeom", fix.arg=initp1, start=list(p2=1/2))
print(f2)
summary(f2)

f2 <- fitdist(x2, "zmgeom", fix.arg=list(p1=1/2), start=list(p2=1/2))
print(f2)
summary(f2)

# (23) check the use of weights with MLE
#

n <- 1e2
x <- rpois(n, 10)
xtab <- table(x)
xval <- sort(unique(x))
f1 <- fitdist(x, "pois")
f2 <- fitdist(xval, "pois", weights = xtab)

f1$estimate
f2$estimate #should be identical


# (24) check the use of weights with other methods
#

n <- 1e2
x <- rpois(n, 20)
xtab <- table(x)
xval <- sort(unique(x))
(f1 <- fitdist(x, "norm", method = "mle"))
(f2 <- fitdist(xval, "norm", weights = xtab, method = "mle"))

(f1 <- fitdist(x, "norm", method = "mme"))
(f2 <- fitdist(xval, "norm", weights = xtab, method = "mme"))

(f1 <- fitdist(x, "norm", method = "qme", probs=c(1/4, 3/4)))
(f2 <- fitdist(xval, "norm", method = "qme", weights = xtab, probs=c(1/4, 3/4) ))

fitdist(x, "norm", method="mge", gof = "CvM")
try(fitdist(xval, "norm", method="mge", gof = "CvM", weights = xtab)) # not yet developped


# (24b) check the use of weights with qme with a discrete distribution
#
n <- 1e2
x <- rpois(n, 10)
xtab <- table(x)
xval <- sort(unique(x))
(f1 <- fitdist(x, "pois", method = "qme", probs=c(1/2)))
(f2 <- fitdist(xval, "pois", method = "qme", weights = xtab, probs=c(1/2) )) # similar to f1

fitdist(xval, "pois", method = "qme", weights = xtab, probs=c(1/2), optim.method="SANN", control=list(maxit=1000)) #
fitdist(x, "pois", method = "qme", probs=c(1/2), optim.method="SANN", control=list(maxit=1000)) # should be similar


# (25) check the warning messages when using weights in the fit followed by functions
# that do not yet take weights into account
# with an example to be used later to see if weights are well taken into account
#

x3 <- rnorm(100)
x3 <- sort(x3)
(f <- fitdist(x3, "norm", method="mle", weights= c(rep(1, 50), rep(2, 50))))
try(plot(f))
try(cdfcomp(f))
(f2 <- fitdist(x3, "logis", method="mle", weights= c(rep(1, 50), rep(2, 50))))
try(cdfcomp(list(f,f2)))
try(denscomp(f))
try(denscomp(list(f,f2)))
try(ppcomp(f))
try(ppcomp(list(f,f2)))
try(qqcomp(f))
try(qqcomp(list(f,f2)))
try(gofstat(f))
try(gofstat(list(f,f2)))
try(bootdist(f))



