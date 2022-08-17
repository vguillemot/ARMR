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
    ctr <- PrettyBarPlot2(signed.ctr[, nfac],
                          threshold = 1 / NROW(signed.ctr),
                          font.size = font.size,
                          signifOnly = FALSE,
                          color4bar = col4,
                          ylab = 'Signed Contributions',
                          ylim = c(1.2*min(signed.ctr[,nfac]),
                                   1.2*max(signed.ctr[,nfac])),
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




