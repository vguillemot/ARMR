##### Internal functions #####

#' Transform angle in rad to deg
#' @param rad numeric
#' @noRd
rad2deg <- function(rad) {(rad * 180) / (pi)}

#' Transform angle in deg to rad
#' @param deg numeric
#' @noRd
deg2rad <- function(deg) {(deg * pi) / (180)}

#' Computes the arc-cosine in degrees
#' @param lecos numeric
#' @noRd
arcosdeg <- function(lecos){rad2deg(acos(lecos))}
#
#' Plot contributions
#'
#' @param signed.ctr signed contribtions
#' @param col4 colors
#' @param nfac number of dimensions
#' @param print plot the result
#' @param stem title
#' @param font.size font size
#' @param horizontal whether the barplot is horizontal
#'
#' @noRd
plotCtr <- function(signed.ctr, col4, nfac = 1,
                    print = FALSE, stem = "Jset",
                    font.size = 5,
                    horizontal = FALSE){
    leTitre <- paste0("Component ",nfac)
    if (sum(signed.ctr[, nfac] > 0) == nrow(signed.ctr)){
      ylim <- c(0,
                1.2*max(signed.ctr[,nfac]))
    }else if(sum(signed.ctr[, nfac] < 0) == nrow(signed.ctr)){
      ylim <- c(1.2*min(signed.ctr[,nfac]),
                0)
    }else{
      ylim <- c(1.2*min(signed.ctr[,nfac]),
                1.2*max(signed.ctr[,nfac]))
    }

    ctr <- PrettyBarPlot2(signed.ctr[, nfac],
                          threshold = 1 / NROW(signed.ctr),
                          font.size = font.size,
                          signifOnly = FALSE,
                          color4bar = col4,
                          ylab = 'Signed Contributions',
                          ylim = ylim,
                          horizontal = horizontal) +
        ggtitle("Contribution Barplots",
                subtitle = leTitre) +
        theme(plot.title = element_text(
            color = "#3E2E8F", size = 20,
            face = "bold"),
            plot.subtitle = element_text(
                color = "#3E2E8F", size = 16,
                face = "italic"),
            plot.caption =  element_text(
                color = "#3E2E8F", size = 14,
                face = "italic"),
            axis.text =  element_text(
                color = "#3E2E8F", size = 12),
            axis.title.x =  element_text(
                color = "#3E2E8F", size = 16,
                face = "italic"))
    if (print){
        png(paste0(stem,nfac,'_Contributions.png'))
        print(ctr)
        dev.off()
        jpeg(paste0(stem,nfac,'_Contributions.jpeg'))
        print(ctr)
        dev.off()
    }
    return(ctr)

}

#' Boostrap for eigenvalues.
#'
#' @param data the data-set;
#' @param center whether to center the data (default to FALSE), see expo.scale;
#' @param scale whether to scale the data ( default to FALSE), see expo.scale;
#' @param design whether there is a design;
#' @param niter number of bootstrap iterations;
#' @param CI.perc vector of the two levels of confidence for the left and right limits of the bootstrap confidence interval (default to \code{c(0.025, 0.0975)});
#' @param suppressProgressBar whether to show a progress bar or not (default to FALSE).
#'
#' @return list containing four objects
#' @export
#'
#' @examples
#' Boot4Eigs(matrix(rnorm(100), 10, 10))
Boot4Eigs <- function(data,
                      center = FALSE,
                      scale = FALSE,
                      design = NULL,
                      niter = 100,
                      CI.perc = c(0.025, 0.975),
                      suppressProgressBar = TRUE){
  # Boostrap the eigenvalues
  # Resampled within groups
  # Private functions
  GetEigs <- function(input4geteigs, center = TRUE, scale = TRUE, data.centered = TRUE)
  { data4eigs <- ExPosition::expo.scale(DATA = input4geteigs, center = center, scale = scale)
  dim4data <- min((nrow(input4geteigs)-data.centered), ncol(input4geteigs))
  eigs <- eigen(t(data4eigs) %*% data4eigs)$values[1:dim4data]
  return(eigs = eigs)}
  # ********* A function inspired by Derek *****************
  # get the bootstrap index values within groups
  boot.design  <- function(design){# Compute Bootstrap
    # indices according to a factor matrix
    boot.index <- vector()
    ZeGroup = as.factor(as.matrix(design))
    ValG = names(table(ZeGroup)) # values of the group
    nGroups = length(ValG)
    # how many groups of observations do we want to look at
    for(i in 1:nGroups){
      boot.index <- c(boot.index,sample(which(ZeGroup==ValG[i ]),
                                        replace=TRUE))
    }
    return(boot.index)
  }
  # **************************************************************
  # first get the fixed effect eigenvalues
  FixedEigs = GetEigs(data, center = center, scale = scale)

  if (suppressProgressBar != TRUE){
    print('Starting Bootstrap.')
    pb <-txtProgressBar(min = 0, max = niter,
                        initial = 0, char = "=",
                        title = 'Bootstrap Iterations', style = 1)
  }

  nComp = length(FixedEigs) # number of components
  nObs = nrow(data) # number of observations
  # create a niter x nComp empty matrix for the bootstrapped eigenvalues
  ZeMatOfEigs = matrix(NA, nrow = niter, ncol = nComp)

  for (m in 1:niter){ # Bootstrap loop
    if (!is.null(design)){ # if there is a design,
      BootInd = boot.design(design) # the bootstrap indices will be generated within groups
      ZeMatOfEigs[m,] =  GetEigs(data[BootInd,], center = center, scale = scale)[1:nComp] # estimate the eigenvalues
    } else { # if there is no design
      BootInd = sample(c(1:nObs), replace = TRUE) # the bootstrap indices are generated by sampling with replacement
      ZeMatOfEigs[m,] = GetEigs(data[BootInd,], center = center, scale = scale)[1:nComp] # estimate the eigenvalues
    }
    if (suppressProgressBar != TRUE){setTxtProgressBar(pb, m)}
  } # End of loop

  colnames(ZeMatOfEigs) <- paste0('Dimension ',seq(1,nComp))
  # Bootstrapped estimates for the eigenvalues
  BootstrappedEstEigs <- apply(ZeMatOfEigs,2,mean)
  # Bootstrapped confidence intervales for the eigenvalues
  CIind <- c(floor(CI.perc[1]*niter), ceiling(CI.perc[2]*niter))
  # Create a nComp x 2 (i.e., low, high) empty matrix that stores the confidence intervals
  ZeMatOfEigsCI = matrix(NA, nrow = nComp, ncol = 2, dimnames = list(colnames(ZeMatOfEigs),c("LowCI", "HighCI")))

  # The inertia of each iteration
  BootstrappedEstInertia <- rowSums(ZeMatOfEigs)
  # Bootstrapped taus
  ZeMatOfTaus <- ZeMatOfEigs/BootstrappedEstInertia
  # Create a nComp x 2 (i.e., low, high) empty matrix that stores the confidence intervals
  ZeMatOfTausCI = matrix(NA, nrow = nComp, ncol = 2, dimnames = list(colnames(ZeMatOfEigs),c("LowCI", "HighCI")))

  for (cp.walk in 1:nComp){
    # order the bootstrapped samples of eigenvalues of each component
    Eigs.ordered <- ZeMatOfEigs[order(ZeMatOfEigs[,cp.walk]),cp.walk]
    # get the lower and upper CIs
    ZeMatOfEigsCI[cp.walk,] <- Eigs.ordered[CIind]

    # order the bootstrapped samples of taus of each component
    Taus.ordered <- ZeMatOfTaus[order(ZeMatOfTaus[,cp.walk]),cp.walk]
    # get the lower and upper CIs
    ZeMatOfTausCI[cp.walk,] <- Taus.ordered[CIind]

  }

  return.list <- list(BootMatEigs = ZeMatOfEigs,
                      BootMatEigsCI = ZeMatOfEigsCI,
                      FixedEigs = FixedEigs,
                      BootstrappedEstEigs = BootstrappedEstEigs,
                      BootMatTaus = ZeMatOfTaus,
                      BootMatTausCI = ZeMatOfTausCI,
                      BootstrappedEstInertia = BootstrappedEstInertia)
  return(return.list)
}  # End of function Boot4Eigs


#' PlotScreeWithCI ----
#' plot the scree for the eigenvalues
#' of an SVD based multivariate analysis,
#' and add bootstrap confidence intervals.
#'
#' \code{PlotScreeWithCI}: Plot the scree for the eigenvalues
#' of an SVD-based multivariate analysis.
#' Note that the function can recompute the
#' eigen-values when a percentage is given.
#' For example  \code{ExPosition} does not return all ev
#'        but only the requested one. but return all percentage
#'        so if max.ev is specified, it is used to recompute
#'        all eigenvalues.
#'  By default \code{PlotScree}
#'  will not plot the line corresponding to
#'  the average inertia (i.e., Kaiser criterion).
#'  If provided with probabilities,
#'  \code{PlotScree} will
#'  color differently the "significant"
#'  eigenvalues.
#' @author Hervé Abdi with help
#' from Derek Beaton and Ju-Chi Yu.
#' @param ev the eigenvalues to plot.
#' No default.
#' @param p.ev the probabilities
#' associated to the
#' eigen-values, (default = \code{NULL}).
#' @param max.ev
#' the max eigenvalue
#'        needed because \code{ExPosition}
#'        does not always return all
#'        eigenvalues
#'        but sometimes only the requested ones;
#'        however \code{ExPosition} always returns
#'        all percentages i.e., \code{tau}),
#'        so if \code{max.ev} is specified,
#'        it is used to recompute
#'        all eigenvalues.
#' @param ci.ev, confidence intervals for the eigenvalues,
#'        computed with Boot4Eigs. Default to \code{NULL}.
#' @param ci.tau, confidence intervals for the taus,
#'        computed with Boot4Eigs. Default to \code{NULL}.
#' @param polygon.ci, the confidence intervals to plot.
#' \code{'tau'} will plot the bootstrap confidence intervals of taus;
#' \code{'ev'} will plot the bootstrap confidence intervals of eigenvalues;
#' \code{NULL} will not plot the polygon. Default to \code{'tau'}.
#' @param alpha
#' threshold for significance
#'   \code{Default = .05}).
#' @param col.ns color for the non significant
#' eigenvalues. Default is \code{'Green'}.
#' @param col.sig  color for significant
#' eigen-values.
#' Default is \code{'Violet'}.
#' @param title a title for the graph
#' default is
#' \code{"Explained Variance per Dimension"}.
#' @param plotKaiser  when \code{TRUE}
#' plot a line corresponding to the average inertia
#' (Kaiser criterion); do not plot when
#' \code{FALSE} (default).
#' @param color4Kaiser
#' color for Kaiser's
#' line
#' (default is \code{'darkorchid4'})
#' @param lwd4Kaiser lwd value (i.e., width)
#' for Kaiser's criterion line.
#' (default is \code{'2.5'})
#' @examples  # PlotScree(ev)
#' @export
#'
PlotScreeWithCI <- function(ev,
                            p.ev = NULL,
                            max.ev = NULL,
                            ci.ev = NULL, # ADDED - JY ---- # results from Boot4Eigs
                            ci.tau = NULL,
                            polygon.ci = 'tau', # ADDED - JY ---- # TRUE/FALSE statement
                            alpha = .05,
                            col.ns = '#006D2C', col.sig = '#54278F',
                            title = "Explained Variance per Dimension",
                            plotKaiser = FALSE,
                            color4Kaiser = 'darkorchid4',
                            lwd4Kaiser = 2.5
){
  # percentage of inertia
  val.tau = (100*ev / sum(ev))
  # ADDED - JY ---------------------------
  if (polygon.ci == 'ev'){
    val.tau.ci = (100*ci.ev / sum(ev))
  }else if (polygon.ci == 'tau'){
    val.tau.ci = 100*ci.tau
  }

  # percentage of inertia for the CIs
  #---------------------------------------
  Top.y = ceiling(max(val.tau.ci) * .1) * 10
  # if ev is already a percentage convert it back
  if (!is.null(max.ev)){ev = ev * (max.ev / ev[1])}
  #
  par(mar = c(5,6,4,4))
  # plot.window(xlim = c(0, length(val.tau)+5),
  #         ylim = c(0,Top.y),asp = .6)
  plot(x = seq(1, length(val.tau)), y = val.tau, xlab = 'Dimensions',
       ylab = 'Percentage of Explained Variance',
       main = title,
       type = 'l', col = col.ns, lwd = 1,
       xlim = c(1, length(val.tau)),
       ylim = c(0,Top.y) # This was the original line
  )
  # ADDED - JY ------------------------------------------------------
  if(!is.null(ci.ev)){ # plot the confidence intervals if ci.ev exist
    arrows(seq(1, length(val.tau)), val.tau.ci[,1], seq(1, length(val.tau)), val.tau.ci[,2], length=0.05, angle=90, code=3)
  }
  if(!is.null(polygon.ci)){
    polygon(c(seq(1, length(val.tau)),seq(length(val.tau),1)),
            c(val.tau.ci[,1],rev(val.tau.ci[,2])),
            col = adjustcolor(col.ns,alpha.f = 0.2),
            border = FALSE
    )
  }
  #------------------------------------------------------------------
  points(x = seq(1,length(val.tau)),y = val.tau,
         pch = 16,  cex = 1, col = col.ns, lwd = 2.5
  )
  if (!is.null(p.ev) & sum(p.ev < alpha) > 0){# plot the significant vp if exist
    # Plot the significant factors
    signi.vp = which(p.ev < alpha)
    # These are the lines Ju-Chi changed ####
    lines(x = seq(1, max(signi.vp)), y = val.tau[1:max(signi.vp)],
          type = "l", col = col.sig, lwd = 1.5)
    points(x = signi.vp, y = val.tau[signi.vp],
           pch = 16, cex = 1.5, col = col.sig, lwd = 3.5)
    #______________________________________________
  } # end of plot significant vp
  par(new = TRUE)
  par(mar = c(5,6,4,4) + .5)
  le.max.vp = Top.y*(ev[1]/val.tau[1])
  plot(ev, ann = FALSE,axes = FALSE,type = "n",#line=3,
       ylim = c(0,le.max.vp))
  if (plotKaiser){
    abline(h = sum(ev)/length(ev),  col = color4Kaiser, lwd = lwd4Kaiser)
  }
  mtext("Inertia Extracted by the Components", side = 4, line = 3)
  axis(4)
} # end of function




