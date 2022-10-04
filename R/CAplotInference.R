# Work file for CA output ----
# Created  09/28/2022
# for package ARMR
# Current Version 09/28/2022
# Goal: generates inference graphs and tables for CA
# Author Luke
#
# Documentation -----
#' @title CAplotInference run a CA
#' (with the \code{ExPosition} package)
#' and generates the inference graphs and tables.
#' Note: *Still Under Development*.
#' @description Generates the inference graphs and tables for CA.
#' @param resCA Output from epCA
#' @param data A data frame or a matrix with
#' data suitable for a CA.
#' @param DESIGN  Default: NULL.
#' A design vector (could be factor or character)
#' or (Boolean) matrix used to assign rows
#' to groups.
#' @param make_design_nominal
#' if TRUE (Default) transform
#' the vector from \code{DESIGN} into
#' a Boolean matrix.
#' Passed to
#' \code{InPosition::epCA.inference.battery}.
#' @param k number
#' of factor to keep; when equql to
#' 0  (Default), all factors are kept.
#' Passed to
#' \code{InPosition::epCA.inference.battery}.
#' @param graphs  do we want graphs? At the moment, this parameter is ignored.
#' Current Default is \code{12} which indicates that
#' the graphs are generated for the first
#' 2 components. Note that current version
#' is creating output only for the first two
#' components.
#' @param printGraphs  (Default: FALSE)
#' do we want to print the graphics as \code{.png}?
#' @param col4I  a color vector for
#' plotting the rows (if \code{NULL}
#' Default) use colors from \code{ExPosition::epCA}.
#' @param col4J
#' a color vector for
#' plotting the columns (if \code{NULL}
#' Default) use colors from \code{ExPosition::epCA}.
#' @param niter.boot Number of bootstrap samples to run. Default: 100
#' @param niter.perm Number of permutations to run. Default: 100
#' @param fast Use the fast version of the perm and bootstrap from data4PCCAR
#' @param save2pptx  Default: FALSE
#' @param title4pptx Title of the PPTX, Default:
#' 'CA Inference Results'.
#' @return A list made of two lists
#'
#' @details Work in Progress
#' @author Luke Moraglia
#' @export
#' @import prettyGraphs
#' @import ggplot2 PTCA4CATA data4PCCAR corrplot
#' @importFrom InPosition epCA.inference.battery
#' @importFrom PTCA4CATA PlotScree createFactorMap createxyLabels.gen
#' @importFrom grDevices colorRampPalette  dev.off  jpeg png recordPlot
#' @importFrom stats cor  cov varimax
CAplotInference <- function(
    resCA,
    data,
    DESIGN = NULL,
    make_design_nominal = TRUE,
    k = 0,
    graphs = 12,
    printGraphs = FALSE,
    col4I = NULL,
    col4J = NULL,
    niter.boot = 100,
    niter.perm = 100,
    fast = FALSE,
    save2pptx = FALSE,
    title4pptx = "CA Inference Results"){

  if (is.null(col4I)) {
    if (is.null(DESIGN)){
      col4I <- resCA$Plotting.Data$fi.col
    }else{
      col4I <- createColorVectorsByDesign(makeNominalData(as.matrix(DESIGN)))$oc
    }
  }
  if (is.null(col4J)) {
    col4J <- resCA$Plotting.Data$fj.col
  }

  if ((is.data.frame(DESIGN)|is.matrix(DESIGN)) && ncol(DESIGN) > 1){
    DESIGN <- DESIGN %*% diag(1:ncol(DESIGN)) |> rowSums()
  }else if((is.data.frame(DESIGN)|is.matrix(DESIGN)) && ncol(DESIGN) == 1){
    DESIGN <- as.vector(as.matrix(DESIGN))
  }


  if(fast){
    res_fast_perm <- data4PCCAR::fastPerm4CA(data, niter.perm)
    res_fast_boot <- data4PCCAR::fastBoot4CA(data, nf2keep = 2, nIter = niter.boot)
    p.ev <- res_fast_perm$pEigenvalues
    BR.I <- res_fast_boot$bootRatios.i
    BR.J <- res_fast_boot$bootRatios.j
  }
  else{
    infres.J <- epCA.inference.battery(
      DATA = data,
      DESIGN = DESIGN,
      make_design_nominal = make_design_nominal,
      k = k,
      test.iters = niter.perm,
      graphs = FALSE)
    infres.I <- epCA.inference.battery(
      DATA = t(data),
      DESIGN = NULL,
      make_design_nominal = make_design_nominal,
      k = k,
      test.iters = niter.perm,
      graphs = FALSE)

    p.ev <- infres.J$Inference.Data$components$p.vals
    BR.I <- infres.I$Inference.Data$fj.boots$tests$boot.ratios
    BR.J <- infres.J$Inference.Data$fj.boots$tests$boot.ratios
  }

  ## Scree
  PlotScree(
    ev = resCA$ExPosition.Data$eigs,
    p.ev = p.ev)
  a01.leScree <- recordPlot()

  ## Bootstrap ratios
  rownames(BR.I) <- rownames(data)
  rownames(BR.J) <- colnames(data)

  laDim = 1


  # Plot the bootstrap ratios for Dimension 1
  ba001.BR1.I <- PrettyBarPlot2(BR.I[,laDim],
                                threshold = 2,
                                font.size = 3,
                                color4bar = col4I, # we need hex code
                                main = paste0('CA: Rows Bootstrap ratio ',
                                              laDim),
                                ylab = 'Bootstrap ratios',
                                horizontal = FALSE
  )

  ba002.BR1.J <- PrettyBarPlot2(BR.J[,laDim],
                                threshold = 2,
                                font.size = 3,
                                color4bar = col4J, # we need hex code
                                main = paste0('CA: Columns Bootstrap ratio ',
                                              laDim),
                                ylab = 'Bootstrap ratios',
                                horizontal = FALSE
  )

  laDim = 2


  # Plot the bootstrap ratios for Dimension 1
  ba003.BR2.I <- PrettyBarPlot2(BR.I[,laDim],
                                threshold = 2,
                                font.size = 3,
                                color4bar = col4I, # we need hex code
                                main = paste0('CA: Rows Bootstrap ratio ',
                                              laDim),
                                ylab = 'Bootstrap ratios',
                                horizontal = FALSE
  )

  ba004.BR2.J <- PrettyBarPlot2(BR.J[,laDim],
                                threshold = 2,
                                font.size = 3,
                                color4bar = col4J, # we need hex code
                                main = paste0('CA: Columns Bootstrap ratio ',
                                              laDim),
                                ylab = 'Bootstrap ratios',
                                horizontal = FALSE
  )


  ## graphs ----
  results.graphs <- list(
    scree = a01.leScree,
    BRI1 = ba001.BR1.I,
    BRJ1 = ba002.BR1.J,
    BRI2 = ba003.BR2.I,
    BRJ2 = ba004.BR2.J
  )
  description.graphs <- list(
    scree = "The Eigenvalues Scree Plot",
    BRI1 = "Row Bootstrap Ratios for Dimension 1",
    BRJ1 = "Column Bootstrap Ratios for Dimension 1",
    BRI2 = "Row Bootstrap Ratios for Dimension 2",
    BRJ2 = "Column Bootstrap Ratios for Dimension 2"
  )

  ## list stat & graphs ----
  results <- list(
    results.stats = NULL,
    results.graphs = results.graphs,
    description.graphs = description.graphs
  )

  if (save2pptx) {
    saveAllGraphsInList2pptx(
      list2Save = results.graphs,
      titles4list2Save = description.graphs,
      file2Save.pptx = paste0(title4pptx, ".pptx"),
      title = title4pptx
    )
  }

  return(results)
}
