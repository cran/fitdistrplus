#############################################################################
#   Copyright (c) 2012 Christophe Dutang, Aurelie Siberchicot
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
### plot density functions for various fits
### of continuous distribution(s) (fitdist results)
### on a same dataset
###
###         R functions
###


qqcomp <- function(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
                   fitpch, fitcol, fitlwd, addlegend = TRUE, legendtext, xlegend = "bottomright", ylegend = NULL, 
                   use.ppoints = TRUE, a.ppoints = 0.5, line01 = TRUE, line01col = "black", line01lty = 1,
                   ynoise = TRUE, plotstyle = "graphics", ...)
{
  if(inherits(ft, "fitdist"))
  {
    ft <- list(ft)
  }else if(!is.list(ft))
  {
    stop("argument ft must be a list of 'fitdist' objects")
  }else
  {
    if(any(sapply(ft, function(x) !inherits(x, "fitdist"))))        
      stop("argument ft must be a list of 'fitdist' objects")
  }
  
  # In the future developments, it will be necessary to check that all the fits share the same weights
  if(!is.null(ft[[1]]$weights))
    stop("qqcomp is not yet available when using weights")
  
  # check the 'plotstyle' argument
  plotstyle <- match.arg(plotstyle[1], choices = c("graphics", "ggplot"), several.ok = FALSE)
  
  # check data
  mydata <- ft[[1]]$data
  verif.ftidata <- function(fti)
  {
    if (any(fti$data != mydata))
      stop("All compared fits must have been obtained with the same dataset")
    invisible()
  }
  lapply(ft, verif.ftidata)
  
  n <- length(mydata)
  sdata <- sort(mydata)
  largedata <- (n > 1e4)
  if (xlogscale != ylogscale)
  {
    warning("As a Q-Q plot should use the same scale on x and y axes, 
            both or none of the axes should be put in a logarithmic scale.")
  }
  logxy <- paste(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""), sep="")
  
  # manage default parameters
  nft <- length(ft)
  if (missing(fitcol)) fitcol <- 2:(nft+1)
  if (missing(fitpch)) fitpch <- ifelse(largedata, 1, 21)
  if (missing(fitlwd)) fitlwd <- 1
  fitcol <- rep(fitcol, length.out=nft)
  fitpch <- rep(fitpch, length.out=nft)
  fitlwd <- rep(fitlwd, length.out=nft)
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
  
  if (use.ppoints)
    obsp <- ppoints(n, a = a.ppoints)
  else
    obsp <- (1:n) / n
  
  # computation of each fitted distribution
  comput.fti <- function(i)
  {
    fti <- ft[[i]]
    para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    qdistname <- paste("q", distname, sep="")
    do.call(qdistname, c(list(obsp), as.list(para)))
  }
  fittedquant <- sapply(1:nft, comput.fti)
  if(NCOL(fittedquant) != nft || NROW(fittedquant) != length(obsp))
    stop("problem when computing fitted CDFs.")
  
  # check limits
  if(missing(xlim))
    xlim <- range(fittedquant)
  if(missing(ylim))
    ylim <- range(mydata)       
  
  if(plotstyle == "graphics") 
  {
    ######## plot if plotstyle=='graphics' ########
    # main plot
    if(!largedata)
      resquant <- plot(fittedquant[,1], sdata, main=main, xlab=xlab, ylab=ylab, log=logxy,
                       pch = fitpch[1], xlim=xlim, ylim=ylim, col=fitcol[1], type="p", ...)
    else
      resquant <- plot(fittedquant[,1], sdata, main=main, xlab=xlab, ylab=ylab, log=logxy,
                       lty = fitpch[1], xlim=xlim, ylim=ylim, col=fitcol[1], type="l", 
                       lwd = fitlwd[1], ...)
    
    #plot of other fitted quantiles
    if(nft > 1 && !ynoise && !largedata)
      for(i in 2:nft)
        points(fittedquant[,i], sdata, pch=fitpch[i], col=fitcol[i], ...)
    if(nft > 1 && ynoise && !largedata)
      for(i in 2:nft)
        if (ylogscale)
        {
          noise2mult <- runif(n, 0.95, 1.05)
          points(fittedquant[,i], sdata*noise2mult, pch=fitpch[i], col=fitcol[i], ...)
        }else
        {
          noise2add <- runif(n, -0.02, 0.02)
          points(fittedquant[,i], sdata+noise2add, pch=fitpch[i], col=fitcol[i], ...)
        }
    
    if(nft > 1 && largedata)
      for(i in 2:nft)
        lines(fittedquant[,i], sdata, col=fitcol[i], lty = fitpch[i], lwd = fitlwd[i], ...)
    
    if(line01)
      abline(0, 1, lty=line01lty, col=line01col)
    
    if (addlegend)
    {
      if(!largedata)
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, col=fitcol, pch = fitpch, ...)
      else
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, col=fitcol, lty = fitpch, lwd = fitlwd, ...)  
    }
    return(invisible(list(obs = sdata, quantiles = fittedquant)))
    
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) 
  {
    stop("ggplot2 needed for this function to work with plotstyle = 'ggplot'. Please install it", call. = FALSE)
  } else 
  {
    ######## plot if plotstyle=='ggplot' ########
    
    # recode the legend position according to available positions in ggplot2
    if(xlegend %in% c("topleft", "bottomleft"))
      xlegend <- "left"
    if(xlegend %in% c("topright", "bottomright"))
      xlegend <- "right"
    
    # structure the fittedquant in a relevant data.frame
    fittedquant <- as.data.frame(fittedquant)
    colnames(fittedquant) <- unlist(lapply(ft, function(X) X["distname"]))
    fittedquant <- stack(fittedquant)
    nfq <- nrow(fittedquant)
    fittedquant$sdata <- sdata   # sdata is recycled in the standard fashion
    fittedquant$ind <- factor(fittedquant$ind, levels = unique(fittedquant$ind))   # reorder levels in the appearance order of the input
    if(nft > 1 && ynoise && !largedata) {
      if (ylogscale)
      {
        noise2mult <- runif(nfq, 0.95, 1.05)
        fittedquant$sdata <- fittedquant$sdata*noise2mult
      }else
      {
        noise2add <- runif(nfq, -0.02, 0.02)
        fittedquant$sdata <- fittedquant$sdata+noise2add
      }
    }
    
    ggqqcomp <-
      ggplot2::ggplot(data = fittedquant, ggplot2::aes(.data$values, .data$sdata, group = .data$ind, colour = .data$ind, shape = .data$ind, linewidth = .data$ind)) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::ggtitle(main) +
      ggplot2::coord_cartesian(xlim = c(xlim[1], xlim[2]), ylim = c(ylim[1], ylim[2])) +
      {if(!largedata) ggplot2::geom_point() else ggplot2::geom_line(ggplot2::aes(linetype = .data$ind, linewidth = .data$ind))} +
      
      {if(addlegend) ggplot2::theme(legend.position = c(xlegend, ylegend), plot.title = ggplot2::element_text(hjust = 0.5)) else ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5))} +
      ggplot2::scale_color_manual(values = fitcol, labels = legendtext) +
      ggplot2::scale_shape_manual(values = fitpch, labels = legendtext) +
      ggplot2::scale_linetype_manual(values = fitpch, labels = legendtext) +
      ggplot2::scale_linewidth_manual(values = fitlwd, labels = legendtext) +
      
      ggplot2::guides(colour = ggplot2::guide_legend(title = NULL)) +
      ggplot2::guides(shape = ggplot2::guide_legend(title = NULL)) +
      ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL)) +
      ggplot2::guides(linewidth = ggplot2::guide_legend(title = NULL)) +
      
      {if(line01) ggplot2::geom_abline(intercept = 0, slope = 1)} +
      {if(xlogscale) ggplot2::scale_x_continuous(trans='log10')} +
      {if(ylogscale) ggplot2::scale_y_continuous(trans='log10')}
    
    return(ggqqcomp)
  }
}
