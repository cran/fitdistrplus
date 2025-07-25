# fitdistrplus 1.2-4

BUG FIX

- fix a convergence issue of `fitdist()` in the FAQ 1.11 for U-shaped beta density. Now use the constrained-version Nelder-Mead algorithm.

# fitdistrplus 1.2-3

NEW FEATURES

- update the FAQ 3.2 with a full example.
- correct man page regarding default starting values.
- add a link to the FAQ website in most man pages.
- add defensive programing in `llplot()`.
- add a link to the Software Heritage archive
- update graphics created in ggplot style, in line with the next version of ggplot2

BUG FIX

- bug fix in `mledist()` when Nelder-Mead algorithm was used and infinite log-likelihood was computed. two test files added.
- starting values for inverse distributions as scale or rate parameters was not inverted.
- for positive unbounded distributions, initial shape paramters are now bounded below and above, see `util-startarg.R`.

# fitdistrplus 1.2-2

NEW FEATURES

- a website bringing together all the resources related to the fitdistrplus package now exists on github.io
at the following URL: https://lbbe-software.github.io/fitdistrplus/
- add a logo for the package thanks to Appsilon.

BUG FIXES

- `mgedist()` may suffer a numerical issue for Anderson-Darling GoF metrics. All GoF metrics now take care of numerical issue, such as `log(0)` or `1/0`, and are properly scaled by the sample sized to avoid large sample size issues. Thanks for Ethan Chapman for reporting the bug.
- `mledist()` now takes care of numerical issue, such as `log(0)`.
- the default starting value for the gamma distribution was wrongly computed for the rate parameter. Thanks for Wendy Martin for reporting the bug.
- `mledist()`, `mmedist()`, `qmedist()` may suffer a scaling issue and the objective function is properly scaled by the sample sized to avoid large sample size issues.

# fitdistrplus 1.2-1

NEW FEATURES

- the `fitdistrplus` git repo now belongs to the `lbbe-software` organization
- modify or add a initial value for all univariate distributions provided in `actuar`.
- create a new vignette regarding default initial values.
- add of generic functions `AIC()` and `BIC()` for `fitdist` and `fitdistcens` objects.
- make `gofstat()` work with `fitdistcens()` objects (giving only AIC and BIC values).
- add calculation of the hessian using `optimHess()` within `fitdist()` when it is not given by `optim()`.
- compute the asymptotic covariance matrix with MME : Now the theoretical moments `m<dist>` should be defined up to an order which equals to twice the maximal order given `order`.
- add a new argument `calcvcov()` in order to (dis)able the computation of covariance matrix for any method.
- graphics function `*comp()` now return a list of drawn points and/or lines when `plotstyle == "graphics"`.
- add a `density()` function for `bootdist(cens)()` objects.
- add DOIs in man pages.

BUG FIXES

- when a `scale` parameter was fixed, the `startarg` function also set a `rate` parameter. this leads to an error when calling the density.
- add further sanity check in `plotdistcens()`: the following code `plotdistcens(data.frame(right=smokedfish$right, left=smokedfish$left))` raised an error via `npsurv()`, thanks to R. Pouillot.
- bug fixed in using `breaks` in `plotdist()`.
- solve the extremely long time taking by `lines()` in `descdist()`.
- add defensive programming for input data (check `NA`, `NaN`, `Inf` values).
- correct links in man pages and URL for DOI.
- remove the use of `plot.np` in vignettes.

# fitdistrplus 1.1-11

NEW FEATURES

- add a `print` argument in the `descdist()` function to allow to only plot the skewness-kurtosis graph, without printing the computed parameters.

BUG FIX

- the use of some deprecated `ggplot2` functions is updated.
- the use of some deprecated BibTeX entries is updated.
- bug fixed in drawing CI lines in `CIcdfcplot()` when `ggplot2` is called.
- bug fixed in drawing horizontal lines in `cdfcompcens()`.

# fitdistrplus 1.1-8

WARNING FIX

- update an URL in fitdistrplus.Rd from \href{https://doi.org/10.18637/jss.v064.i04}{} to
\doi{doi.org/10.18637/jss.v064.i04}
- replace `if(class(x) == XX)` by `if(inherits(x, XX))`
- replace all `dontrun` tags by `donttest` in examples in rd files

BUG FIX

- fix an error in `t-detectbound.R` producing "failure: length > 1 in coercion to logical"
reported by Brian Ripley


# fitdistrplus 1.1-6

NEW FEATURES

- new function `Surv2fitdistcens()` to format data for use in `fitdistcens()` from a format
used in the survival package
- new dataset `fremale` in order to illustrate `Surv2fitdistcens()`
- support the use of ggplot2 for `CIcdfplot()`
- add the taxon names to the `endosulfan` dataset
- new argument name.points in `cdfcomp()` and `CIcdfplot()` to add labels next to points


# fitdistrplus 1.1-5

WARNING FIX

- reduce testing times in test files.


# fitdistrplus 1.1-3

NEW FEATURE

- take into account `fix.arg` for uniform distribution.

BUG FIXES

- add the loglikelihood value for uniform distribution (in `mledist()`),
- correct usage of triple dots argument in `llsurface()`,
- fix an error in `ppcomp()` and `qqcomp()` raised for large dataset.


# fitdistrplus 1.1-1

NEW FEATURES

- add of internal functions to cope with problems of lack of maintenance of the package `npsurv` and remove the dependence to this package.
- remove of the deprecated argument Turnbull of `plotdistcens()`.


# fitdistrplus 1.0-14

NEW FEATURES

- add a new estimation method called maximum spacing estimation via `msedist()`.


# fitdistrplus 1.0-13

BUG FIXES

- fix issues coming from the noLD (--disable-long-double) configuration of R.


# fitdistrplus 1.0-12

BUG FIXES

- bug fixed in `qmedist()` and `fitdistcens()` which raised an error in `checkparamlist()`.
- bug fixed in `testdpqfun()` which assumes the first argument of d,p,q,r functions are exactly the same as in base R.


# fitdistrplus 1.0-11

NEW FEATURES

- update the FAQ with beta(a,a).
- improve graphics for discrete distributions in `denscomp()`.
- improve automatic naming of legends in `xxxcomp()`.
- harmonize outputs in `mledist()`, `qmedist()`, `mmedist()`, `mgedist()`, `fitdist()` and 
`fitdistcens()`.
- automatic test of d, p, q functions in `fitdist()` and raise warnings.
- improve test for starting and fixed values.
- add new default starting values for distributions in actuar.
- change of the default CDF plot for censored data, using the Wang NPMLE algorithm provided in the package npsurv (in `plotdistcens()` and `cdfcompcens()`)
- add of two new goodness-of-fit plots (QQ-plot and PP-plot) for censored data (cf. `plotdistcens`, `qqcompcens` and `ppcompcens`). 
- add of a part dedicated to censored datain the FAQ vignette.
- homogeneization of `xlim` and `ylim` default definition in `plotdistcens`.
- Removing of the name of the first argument in calls to dpq functions in order
to make the package compatible with distributions defined with a non classical
name for their first argument (resp. x, q, p for d, p, q functions).
- add the possibility to change the title of the CDF plot in `plotdistcens()` using the argument main.
- support the use of `ggplot2` for `cdfcompcens()`, `qqcompcens()`, `ppcompcens()`.

BUG FIXES

- bug fixed concerning the use of gofstat with a chi squared `df <=0` (error message blocking the other functions)
- bug fix in `mledist()` when bounds were set (so not NULL) for censored MLE.
- enable a correct use of non-equidistant breaks in `denscomp()` for the histogram.
when `plotstyle = "ggplot"`, and prohibit the use of non-equidistant breaks with `probability = FALSE`
(adding a stop in this case).


# fitdistrplus 1.0-9

- update the FAQ with linear inequality constraints.


# fitdistrplus 1.0-8

NEW FEATURES

- support the use of ggplot2 for `cdfcomp()`, `denscomp()`, `qqcomp()`, `ppcomp()`.

BUG FIXES

- correct legend for `qqcomp()` and `ppomp()` on large data.
- correct weights in `mmedist()`.
- correct the name Akaike in `gofstat()`.
- correct the use of trueval in `plot.bootdist()`.
- correct the vignette on truncate (inflated) distributions.


# fitdistrplus 1.0-7

NEW FEATURES

- keep the JSS vignette as a pdf.
- start the FAQ vignette and add datasets (?dataFAQ) for it.
- provide likelihood plot/surface/curve: `llplot()`, `llcurve()`, `llsurface()`.
- provide parallelization of bootstrap in bootdist and bootdistcens.
- provide graphic of (e)cdf with bootstraped confidence interval/area: `CIcdfplot()`.
- allow the use of `constrOptim()` in `mledist()`, `mmedist()`, `mgedist()`, `qmedist()` functions.
- add a possible pre-fitting procedure: prefit.

BUG FIXES

- add `invisible()` for all graphical functions.
- bug fixed concerning the use of weights on censored data.


# fitdistrplus 1.0-6

BUG FIXES

- automatic definition of starting values for distributions `llogis` and `invweibull`
is now working.


# fitdistrplus 1.0-5

NEW FEATURES

- update starting/fixing values in `mledist()`, `mmedist()`, `mgedist()`, `qmedist()` functions.
- update graphics for bootstrap procedure.
- add argument do.points in `cdfcomp()`.
- add argument weights in `mledist()`, `qmedist()`, `mmedist()`, `fitdist()`, `fitdistcens()`.
- add argument keepdata in `fitdist()`, `fitdistcens()`.
- suppress warnings/errors in `fitdist(cens)()`, `bootdist(cens)()`.

BUG FIXES

- defensive programming in `plotdist()`, `cdfcomp()`,...
- simplify plotting curves in `cdfcomp()` where seq(xmin, xmax, by=1) was used.


# fitdistrplus 1.0-4

- release for the JSS publication.


# fitdistrplus 1.0-3

NEW FEATURES

- new generic functions for fitdist(cens): `loglik()`, `vcov()` and `coef()`.
- vignette updated to the version of a paper accepted by the 
Journal of Statistical Software.
- add of an argument discrete in fitdist in order to be able to take into
account non classical discrete distributions while plotting the fit with plot.fitdist or cdfcomp 
and while calculating goodness-of-fit statistics with gofstat (add of an example : fit of a zero
inflate Poisson distribution).
- add of an S3 class for `descdist()` and a print method.

BUG FIXES

- fitdist can handle non invertible Hessian matrices.


# fitdistrplus 1.0-2

NEW FEATURES

- plotdist can plot empirical density as an histogram, a density plot 
or both superimposed.
- a strong warning was added to the documentation of function descdist
about the problematic high variance of skewness and kurtosis.

BUG FIXES

- bug fixed in `bootdistcens()` : argument fix.arg is now correctly passed to mle.


# fitdistrplus 1.0-1

NEW FEATURES

- `gofstat` can handle multiple `fitdist()` objects.
- `plotdist` for discrete data is slightly enhanced.


# fitdistrplus 1.0-0

NEW FEATURES

- update cdfcomp and add denscomp, ppcomp and qqcomp functions.
- add of an argument Turnbull.confint to functions plotdistcens and cdfcompcens in order to 
    draw confidence intervals on the empirical distribution only if requested.
- ppoints now used in `fitdist` for QQ plot, PP plot and cdf plot for continuous data
    (was used only for QQ plot in previous versions) to enable Blom type plotting position
    (using by default Hazen plotting position than can be chanfge using arguments use.ppoints and a.ppoints)
- many changes in the examples given in the reference manual.
- the vignette was removed, to be transformed in a paper that we will soon submit to a journal.
- add of four data sets : `fluazinam`, `salinity`, `danishuni` and `danishmulti`.
- add of functions to calculate quantiles of the fitted distribution, with 95 percent CI calculated by bootstrap :
  quantile generic function is available both for `fitdist()` and `bootdist()` objects and 
  quantile generic function is available both for `fitdistcens()` and `bootdistcens()` objects.

BUG FIXES

- correction the formula for the CvM test for Weibull distribution.
- elimination of CvM and AD tests for normal, lognormal and logistic distributions : 
formulas previously used (given by Stephens 1986) do not use exactly MLE estimates and thus
results were only approximates.
- make arguments xlim and ylim functional in `cdfcompcens()`.
- bug fix in the closed formula in `mmedist()` for lognormal distributions.


# fitdistrplus 0.3-4

NEW FEATURES

- posibility to fix xlegend to a keyword (e.g. `bottomright`) in `cdfcomp()` and `cdfcompdens()`.
- improvement of the new vignette.

BUG FIXES

- correction of the NAMESPACE file in order to enable the correct print of a summary of a `fitdistcens` object
(with the correlation matrix, the loglikelihood and AIC and BIC statistics).


# fitdistrplus 0.3-3

NEW FEATURES

- a new function (`cdfcompcens()`) to plot cumulative distributions corresponding to various fits using a same 
censored data set.
- add an example with scaling problem in man pages.


# fitdistrplus 0.3-2

NEW FEATURES

- new plot of the empirical cdf curve in `plotdistcens()`, using the Turnbull algorithm by a call
to function survfit{survival}.
- new arguments to function `cdfcomp()` : verticals, horizontals and xlim.


# fitdistrplus 0.3-1

NEW FEATURES

- add of a draft of a new version of the vignette.


# fitdistrplus 0.3-0

NEW FEATURES

- a new function (`cdfcomp()`) to plot cumulative distributions corresponding to various fits using a same 
non censored data set.
- add of two data sets : `endosulfan` and `toxocara`.


# fitdistrplus 0.2-2

BUG FIXES

- elimination of NON-ASCII characters in the vignette.


# fitdistrplus 0.2-1

NEW FEATURES

- a new fitting method was implemented for continuous distributions : the maximum goodness-of-fit
estimation (function `mgedist()`) (for the moment only available for non censored data).


# fitdistrplus 0.1-5

NEW FEATURES

- a new goodness-of-fit statistic was added in `gofstat()`, with corresponding test :
the Cramer-von Mises distance.
- a new fitting method has been implemented : the quantile matching estimation (function `qmedist()`).
For the moment, only available for non censored data.
- the moment matching estimation has been extended (in function `mmedist()`) to enable numerical
matching when closed formula are not available.

BUG FIXES

- correction of a bug inserted while adding the argument `fix.arg` which prevent the print of the
results of goodness-of-fit tests.


# fitdistrplus 0.1-4

NEW FEATURES

- a component named dots added to the list returned by `fitdist()` and `fitdistcens()` in order to pass optional arguments
for the control of optimization in `mledist()` to `bootdist()` and `bootdistcens()`. `bootdist()` and `bootdistcens()`
changed to take into account these optional arguments if they are defined in the call to `fitdist()` or `fitdistcens()`.
- an argument added to `fitdist()`, `fitdistcens()` and `mledist()`, named `fix.arg`, and giving the possibility to fix some
of the distribution parameters while maximizing the likelihood. Functions bootdist, bootdistcens and gofstat were also
changed in order to take this new argument into account.
- a new data file of bacterial contamination censored data extracted from Busschaert et al. 2000 and examples
corresponding to analysis of this dataset.

BUG FIXES

- correction of a bug in the print and the plot of bootstraped samples using `bootdist()` or `bootdistcens()`
when there was only one parameter estimated by maximum likelihood.


# fitdistrplus 0.1-3

NEW FEATURES

- new data file `groundbeef` (groundbeef.rda and groundbeef.Rd) and new use of this dataset in some examples.
- new function `gofstat()`.
    Goodness-of-fit statistics are no more computed by fitdist but may computed and printed
    by the use of the function gofstat. In this new function, the whole results computed are not printed :
    results of tests are printed only if the argument print.test==TRUE and for continuous distributions
    only Anderson-Darling and Kolomogorov-Smirnov statistics are printed by default (but complete results
    are returned by gofstat).
- modifications in `descdist()` : three arguments were added in descdist
    1/ method, to choose between unbiased estimations
    of standard deviation, skewness and kurtosis (default choice) and
    sample values.
    2/ obs.col to choose the color used to plot the observed point on the graph.
    3/ boot.col to choose the color used to plot the bootstrap sample of points on the graph.
- modifications in `plotfit()` : minor changes were performed in order to facilitate the use of the argument ...
    to personnalize the plots (see examples in plotdist.Rd)
- modication of the vignette

BUG FIXES

- correction of a bug in `plotdist()` due to the redefinition in the previous version
    of the parameter `"ylim"` for the plot of a histogram with theoretical density function
    (there was a problem with infinite values of theoretical density function).


# fitdistrplus 0.1-2

NEW FEATURES

- deletion of mledistcens and modification of mledist in order
to maximize likelihood for both censored and non censored data.
- possibility to choose the optimization method used for maximum
likelihood estimation (MLE) of distribution parameters using the new
argument `"optim.method"` of `mledist()`.
- possibility to specify contraints on distribution parameters using
the new arguments `"lower"` and `"upper"` of `mledist()`.
- possibility to use a custom optimization function for MLE using the
new argument `"custom.optim"`.
- moment matching estimation is no longer done with argument method
set to `"mom"` but set to `"mme"` in `fitdist()`.
- renaming of `momdist` in `mmedist()`.
- calculation of AIC and BIC criterion after maximum likelihood
estimation of distribution parameters
- change of the default number of iterations from 999 to 1001
for bootstrap in order to avoid interpolation using the quantile function
- use of the argument `"log"` and (resp. `"log.p"`) of density (resp.
distribution) when available to compute the loglikelihood.

BUG FIXES

- omitting the name of the first argument in calls to the density function
during maximization of the likelihood in order to enable the use of a density
function defined with a first parameter (the vector of quantiles) with a name
differing from "x" (classical name for density distributions defined in R),
such as the density function `dexGAUS` from the package gamlss.dist.


# fitdistrplus 0.1-1

- Initial release.
