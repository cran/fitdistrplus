#############################################################################
#   Copyright (c) 2018 Marie Laure Delignette-Muller, Christophe Dutang, 
#                      Aurelie Siberchicot
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the
#   Free Software Foundation, Inc.,
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA
#
#############################################################################
### QQ plot for various fits
### of continuous distribution(s) (fitdistcens results)
### on a same dataset
###
###         R functions
###


qqcompcens <- function(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, fillrect,
                       fitcol, fitlwd, addlegend = TRUE, legendtext, xlegend = "bottomright", ylegend = NULL, 
                       line01 = TRUE, line01col = "black", line01lty = 1, ynoise = TRUE, 
                       NPMLE.method = "Wang", plotstyle = "graphics", ...)
{
  if(inherits(ft, "fitdistcens"))
  {
    ft <- list(ft)
  }else if(!is.list(ft))
  {
    stop("argument ft must be a list of 'fitdistcens' objects")
  }else
  {
    if(any(sapply(ft, function(x) !inherits(x, "fitdistcens"))))        
      stop("argument ft must be a list of 'fitdistcens' objects")
  }
  
  NPMLE.method <- match.arg(NPMLE.method, c("Wang", "Turnbull.intervals", "Turnbull.middlepoints"))
  if (NPMLE.method == "Turnbull.middlepoints")
  {
    warning("The QQcomp plot for censored data is not available with NPMLE.method at Turnbull.middlepoints. 
            Turnbull.intervals will be used instead of Turnbull.middlepoints.")
    NPMLE.method <- "Turnbull.intervals"
  }
  
  # check the 'plotstyle' argument
  plotstyle <- match.arg(plotstyle[1], choices = c("graphics", "ggplot"), several.ok = FALSE)
  
  # In the future developments, it will be necessary to check that all the fits share the same weights
  if(!is.null(ft[[1]]$weights))
    stop("qqcompcens is not yet available when using weights")
  
  censdata <- ft[[1]]$censdata
  
  # check data
  verif.ftidata <- function(fti)
  {
    if (any(fti$censdata$left != censdata$left, na.rm=TRUE) | 
        any(fti$censdata$right != censdata$right, na.rm=TRUE))
      stop("All compared fits must have been obtained with the same dataset")
  }
  l <- lapply( ft, verif.ftidata)
  rm(l)
  
  if (xlogscale != ylogscale)
  {
    xlogscale <- ylogscale <- TRUE
    warning("As a Q-Q plot should use the same scale on x and y axes, 
            both axes were put in a logarithmic scale.")
  }
  logxy <- paste(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""), sep="")
  
  # manage default parameters
  nft <- length(ft)
  if (missing(fitcol)) fitcol <- 2:(nft+1)
  if (missing(fitlwd)) fitlwd <- 1
  fitcol <- rep(fitcol, length.out=nft)
  fitlwd <- rep(fitlwd, length.out=nft)
  if (missing(fillrect)) 
    if ((nft == 1) | plotstyle == "ggplot") fillrect <- "lightgrey" else fillrect <- NA
  
  
  # check legend parameters if added
  if(missing(legendtext)) 
  {
    legendtext <- sapply(ft, function(x) x$distname)
    if(length(legendtext) != length(unique(legendtext)))
      legendtext <- paste(legendtext, sapply(ft, function(x) toupper(x$method)), sep="-")
    if(length(legendtext) != length(unique(legendtext)))
      legendtext <- paste(legendtext, 1:nft, sep="-")
  }
  
  if (missing(xlab))
    xlab <- "Theoretical quantiles"
  if (missing(ylab)) 
    ylab <- "Empirical quantiles"
  if (missing(main)) 
    main <- "Q-Q plot"
  
  # computation from censdata
  f <- npmle(censdata, method = NPMLE.method)
  
  bounds <- c(f$right, f$left)
  finitebounds <- bounds[is.finite(bounds)]
  
  if(missing(xlim) & missing(ylim))
  {
    user.defined.lim <- FALSE
    upper <- max(finitebounds)
    lower <- min(finitebounds)
    width <- upper - lower
    if (xlogscale == TRUE)
    {
      xmin <- lower * (upper / lower)^(-0.1)
      xmax <- upper * (upper / lower)^0.1
      xmininf <- lower * (upper / lower)^(-100) # 100 to be very large
      xmaxinf <- upper * (upper / lower)^100
    } else
    {
      xmin <- lower - width * 0.1
      xmax <- upper + width * 0.1
      xmininf <- lower - width * 100
      xmaxinf <- upper + width * 100
    }
    xlim <- c(xmin, xmax)
    ylim <- c(xmin, xmax)
  } else # at least xlim or ylim are specified
  {
    user.defined.lim <- TRUE
    if (missing(xlim) | missing(ylim))
    {
      warning("By default the same limits are applied to x and y axes.
            You should specify both if you want different x and y limits")
      if (missing(xlim)) xlim <- ylim else ylim <- xlim
    }
    lower <- min(c(xlim, ylim))
    upper <- max(c(xlim, ylim))
    width <- upper - lower
    if (xlogscale == TRUE)
    {
      xmininf <- lower * (upper / lower)^(-100) # 100 to be very large
      xmaxinf <- upper * (upper / lower)^100
    } else
    {
      xmininf <- lower - width * 100
      xmaxinf <- upper + width * 100
    }
  }    
  
  k <- length(f$left)
  Fnpsurv <- cumsum(f$p) 
  Fbefore <- c(0, Fnpsurv[-k])
  df <- data.frame(left = f$left, right = f$right)
  
  # Definition of vertices of each rectangle
  Qi.left <- df$left # dim k
  Qi.left4plot <- Qi.left
  
  # when R is configured with noLD (--disable-long-double), qnorm and other 'q' functions
  # produce NaN values instead of Inf values for 0 and first argument.
  if (is.infinite(Qi.left4plot[1]) | is.nan(Qi.left4plot[1])) Qi.left4plot[1] <- xmininf
  Qi.right <- df$right
  Qi.right4plot <- Qi.right
  if (is.infinite(Qi.right4plot[k]) | is.nan(Qi.right4plot[k])) Qi.right4plot[k] <- xmaxinf
  # keep only 16 significants digits for R configured with noLD (--disable-long-double)
  Pi.low <- signif(Fbefore, 16)
  Pi.up <- signif(Fnpsurv, 16)
  nPi <- length(Pi.low)
  
  lrect <- vector(mode = "list", length = nft)
  theo.xmin <- xlim[1]
  theo.xmax <- xlim[2]
  for(i in 1:nft)
  {
    fti <- ft[[i]]
    para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    qdistname <- paste("q", distname, sep="")
    
    if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois")))
      warning(" Be careful, variables are considered continuous in this function!")
    
    Qitheo.left <- do.call(qdistname, c(list(Pi.low), as.list(para)))
    Qitheo.right <- do.call(qdistname, c(list(Pi.up), as.list(para)))
    theo.xmin <- min(theo.xmin, Qitheo.right[-k])
    Qitheo.left4plot <- Qitheo.left
    theo.xmax <- max(theo.xmax, Qitheo.left[-1])
    if (is.infinite(Qitheo.left4plot[1]) | is.nan(Qitheo.left4plot[1])) Qitheo.left4plot[1] <- xmininf
    Qitheo.right4plot <- Qitheo.right
    if (is.infinite(Qitheo.right4plot[k]) | is.nan(Qitheo.right4plot[k])) Qitheo.right4plot[k] <- xmaxinf
    lrect[[i]] <- data.frame(Qitheo.left4plot = Qitheo.left4plot, 
                             Qi.left4plot = Qi.left4plot, 
                             Qitheo.right4plot = Qitheo.right4plot, 
                             Qi.right4plot = Qi.right4plot, ind = legendtext[i])
  }
  
  # insert here a check of limits in order to enlarge xlim and ylim if needed
  # in order to be sure to visualize each interval, for all the fitted distributions
  if (!user.defined.lim)
  {
    xlim <- c(theo.xmin, theo.xmax)
    ylim <- c(theo.xmin, theo.xmax)
  }
  
  if(plotstyle == "graphics") 
  {
    ######## plot if plotstyle=='graphics' ########
    # main plot
    plot(1, 1, type = "n", main = main, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, log = logxy)
    
    # plot of rectangles
    plot.fti <- function(i, ...)
    {
      Qitheo.left4plot <- lrect[[i]]$Qitheo.left4plot 
      Qi.left4plot <- lrect[[i]]$Qi.left4plot 
      Qitheo.right4plot <- lrect[[i]]$Qitheo.right4plot 
      Qi.right4plot <- lrect[[i]]$Qi.right4plot
      
      if (ynoise & nft > 1)
      {
        if (xlogscale == TRUE)
        {
          noise2mult <- runif(nPi, 0.99, 1.01)
          rect(xleft = Qitheo.left4plot, ybottom = Qi.left4plot * noise2mult, 
               xright = Qitheo.right4plot, 
               ytop = Qi.right4plot * noise2mult, 
               border = fitcol[i], col = fillrect[i], lwd = fitlwd[i])
        }
        else
        {
          noise2add <- runif(nPi, -width*0.01, width*0.01)
          rect(xleft = Qitheo.left4plot, ybottom = Qi.left4plot + noise2add, 
               xright = Qitheo.right4plot, 
               ytop = Qi.right4plot + noise2add, 
               border = fitcol[i], col = fillrect[i], lwd = fitlwd[i])
        }
      } else # ! ynoise
      {
        rect(xleft = Qitheo.left4plot, ybottom = Qi.left4plot, xright = Qitheo.right4plot, 
             ytop = Qi.right4plot, 
             border = fitcol[i], col = fillrect[i], lwd = fitlwd[i])
      }
    }
    s <- sapply(1:nft, plot.fti, ...)
    rm(s)
    
    if(line01)
      abline(0, 1, lty = line01lty, col = line01col)
    
    if (addlegend)
    {
      legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, col=fitcol, lty=1, lwd=fitlwd, ...)
    }
    invisible()
    
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) 
  {
    stop("ggplot2 needed for this function to work with plotstyle = 'ggplot'. Please install it", call. = FALSE)
  } else 
  {
    ######## plot if plotstyle=='ggplot' ########
    drect <-  do.call("rbind", lrect)
    ind <- as.factor(drect$ind)
    fitcol <- rep(fitcol, table(ind))
    fitlwd <- rep(fitlwd, table(ind))
    fillrect <- if(length(fillrect) > 1) {rep(fillrect, table(ind))} else {fillrect}
    
    ggqqcompcens <- ggplot2::ggplot(drect) + 
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)  +
      ggplot2::ggtitle(main) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
      ggplot2::geom_rect(data=drect, mapping=ggplot2::aes(xmin=.data$Qitheo.left4plot, xmax=.data$Qitheo.right4plot, ymin=.data$Qi.left4plot, ymax=.data$Qi.right4plot), colour = fitcol, fill = fillrect, linewidth = fitlwd, alpha=0.5) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      {if(line01) ggplot2::geom_abline(ggplot2::aes(slope = 1, intercept = 0), color = line01col, linetype = line01lty)} +
      {if(xlogscale) ggplot2::scale_x_continuous(trans='log10')} +
      {if(ylogscale) ggplot2::scale_y_continuous(trans='log10')} + 
      ggplot2::facet_wrap(~ind)
    
    return(ggqqcompcens)
  }
}
