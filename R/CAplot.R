# Work file for CA output ----
# Created  09/19/2022
# for package ARMR
# Current Version 09/21/2022
# Goal: generates standard graphs and tables for CA
# Author Luke and Hervé
# Entête -----
# install.packages('sinew')
# sinew::makeOxygen(graph4epPCA)
# Preamble ----
# Pass Results from Exposition
# ExPosition
#
# Documentation -----
#' @title CAplot run a CA
#' (with the \code{ExPosition} package)
#' and generates the standard graphs and tables.
#' Note: *Still Under Development*.
#' @description Generates the standard graphs and tables for CA.
#' @param resCA Output from epCA
#' @param data A data frame or a matrix with
#' data suitable for a CA.
#' @param DESIGN  Default: NULL.
#' A design vector (could be factor or character)
#' or (Boolean) matrix used to assign rows
#' to groups.
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
#' @param show.TI whether to plot the tolerance intervals or not. Default: FALSE
#' @param show.CI whether to plot the confidence intervals or not. Default: TRUE
#' @param mean.cex the size of the dots of the means. Default: 3
#' @param mean.textcex the size of the texts of the means. Default: 3
#' @param display.labels.ind If TRUE, the labels of rows will be printed. Default: FALSE.
#' @param display.labels.var If TRUE, the labels of columns will be printed. Default: TRUE.
#' @param display.points.mean If TRUE, the mean row factor scores will be plotted. Default: TRUE.
#' @param mean.constraints A list of the constraints (that include \code{minx}, \code{miny}, \code{maxx}, and \code{maxy})
#' The constraints of the figure that only includes the means. This constraints
#' will be used if \code{only.mean = TRUE}. Default: NULL
#' @param scale.mean.constraints A value used to scale the constraints (by multiplication).
#' This function is used to adjust the constraints when the confidence or the tolerance intervals are outside of the figure.
#' Default: 1.5
#' @param max.n4bar When the number of bars exceed this value, the labels will be hidden. Default: 40.
#' @param save2pptx  Default: FALSE
#' @param title4pptx Title of the PPTX, Default:
#' 'CA Results'.
#' @return A list made of two lists
#'
#' @details Work in Progress
#' @author Luke Moraglia and Hervé Abdi
#' @seealso
#'  \code{\link[ExPosition]{epCA}}
#'  \code{\link[PTCA4CATA]{PlotScree}}, \code{\link[PTCA4CATA]{createFactorMap}}
#' @export
#' @import prettyGraphs
#' @import ggplot2 PTCA4CATA data4PCCAR corrplot
#' @importFrom ExPosition epCA makeNominalData
##  @importFrom
##  @importFrom PTCA4CATA PlotScree createFactorMap createxyLabels.gen
#' @importFrom grDevices colorRampPalette  dev.off  jpeg png recordPlot
#' @importFrom stats cor  cov varimax chisq.test
CAplot <- function(resCA,
                   data,
                   DESIGN = NULL,
                   #make_design_nominal = TRUE,
                   #k = 0,
                   graphs = 12,
                   printGraphs = FALSE,
                   col4I = NULL,
                   col4J = NULL,
                   show.TI = FALSE,
                   show.CI = TRUE,
                   mean.cex = 3,
                   mean.textcex = 3,
                   display.labels.ind = FALSE,
                   display.labels.var = TRUE,
                   display.points.mean = TRUE,
                   mean.constraints = NULL,
                   scale.mean.constraints = 1.5,
                   max.n4bar = 40,
                   save2pptx = FALSE,
                   title4pptx = "CA Results"){

  #Debug graphs
  printTest <- FALSE

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

  # Heatmap of Pearson residuals

  chi2 <- chisq.test(data)
  Inertia.cells <- chi2$residuals / sqrt(sum(data))
  if(nrow(Inertia.cells) > ncol(Inertia.cells)){
    Inertia.cells <- t(Inertia.cells)
  }
  col <-
    colorRampPalette(c("#BB4444", "#EE9988",
                       "#FFFFFF", "#77AADD", "#4477AA"))


  corrplot::corrplot(Inertia.cells, is.cor = FALSE,
                     col = col(200),
                     #method = "color",
                     tl.col = "black",
                     tl.srt = 45
                     )
  a0.residuals <- recordPlot()


  # Scree
  scree <- PTCA4CATA::PlotScree(ev = resCA$ExPosition.Data$eigs,
                                plotKaiser = TRUE)
  a01.leScree <- recordPlot()

  #_________________________________________________
  ## The I-set (plot the rows)----
  #_________________________________________________

  # look at the help for PTCA4CATA::createFactorMap

  Fi   <- resCA$ExPosition.Data$fi
  Fj   <- resCA$ExPosition.Data$fj
  constraintsIJ <- minmaxHelper(mat1 = Fi,
                                mat2  = Fj)


  jolie.ggplot1 <- PTCA4CATA::createFactorMap(
    Fi,
    constraints = constraintsIJ,
    col.points = col4I,
    col.labels = col4I,
    display.labels = display.labels.ind,
    font.face = 'italic'
  )

  jolie.ggplot1.transparent <- PTCA4CATA::createFactorMap(
    Fi,
    constraints = constraintsIJ,
    col.points = col4I,
    col.labels = col4I,
    display.labels = display.labels.ind,
    alpha.points = 0.2,
    font.face = 'italic'
  )

  # Create the labels for Inertia per dimension
  label4Map <- PTCA4CATA::createxyLabels.gen(1,
                                             2,
                                             lambda = resCA$ExPosition.Data$eigs,
                                             tau = resCA$ExPosition.Data$t)
  a3.JolieggMap <- jolie.ggplot1$zeMap +
    label4Map +
    labs(title = 'The Row Map')

  if(printTest){
    print(a3.JolieggMap)
  }

  # To make some changes at the labels:
  #  use and modify this code
  label4Map2 <- list(label4Map,
                     # The standard label from createxyLabels.gen
                     theme(
                       axis.title = element_text(
                         # the new theme
                         color = "darkorchid4",
                         size = rel(1.1),
                         # relative to default
                         # family = 'Times',
                         # "Times", "sans", "Courier"
                         face   = "italic" ,
                         # 'plain','italic', 'bold',
                         # NB: face does not work with current ggplot2
                       ),
                       # end of element_text
                       plot.title = element_text(color = '#5826A3')
                     ))

  a4.JolieggMap <-
    jolie.ggplot1$zeMap + label4Map2
  if (printTest) {
    print(a4.JolieggMap)
  }

  # Final row factor scores graph
  a4.JolieggMap.2 <-
    a3.JolieggMap + label4Map2
  if(printTest){
    print(a4.JolieggMap.2)
  }

  if (printGraphs) {
    png('RowImap.png')
    print(a4.JolieggMap.2)
    dev.off()
  }

  # to look at the map with or without the means

  if (!is.null(DESIGN)) {
    fi.mean <- getMeans(resCA$ExPosition.Data$fi, DESIGN)
    colnames(fi.mean) <- paste0("Dimension ", 1:ncol(fi.mean))
    if (is.null(mean.constraints)) {
      mean.constraints <- lapply(minmaxHelper(fi.mean), '*', scale.mean.constraints)
    }

    grpidx.tmp <- tapply(seq_along(col4I), col4I, identity)[unique(col4I)]
    grpidx <- sapply(grpidx.tmp, "[[", 1)
    col4CI <- names(grpidx)
    names(col4CI) <- DESIGN[grpidx]
    # reorder according to fi.mean
    col4CI <- col4CI[rownames(fi.mean)]

    plot.fi.mean <- createFactorMap(fi.mean,
                                    col.background = NULL,
                                    col.axes = "orchid4",
                                    alpha.axes = 0.5,
                                    #title = "The Observation Map.",
                                    col.points = col4CI[rownames(fi.mean)],
                                    col.labels =  col4CI[rownames(fi.mean)],
                                    constraints = mean.constraints,
                                    display.points = display.points.mean,
                                    cex = mean.cex,
                                    text.cex = mean.textcex,
                                    pch = 17,
                                    alpha.points = 0.8)

    bootCI.res <- Boot4Mean(resCA$ExPosition.Data$fi, DESIGN)
    colnames(bootCI.res$BootCube) <- paste0("Dimension ", 1:ncol(resCA$ExPosition.Data$fi))
    fi.CI <- MakeCIEllipses(bootCI.res$BootCube,
                            col = col4CI,
                            alpha.ellipse = 0.1,
                            line.size = 0.5, alpha.line = 0.2)

    data4TI <- resCA$ExPosition.Data$fi
    colnames(data4TI) <- paste0("Dimension ", 1:ncol(data4TI))
    fi.TI <- MakeToleranceIntervals(data4TI,
                                    design = DESIGN,
                                    axis1 = 1, axis2 = 2,
                                    col = col4CI,
                                    line.size = 1,
                                    alpha.ellipse = 0.05, alpha.line = 0.3,
                                    p.level = .80)

    a5.JolieggMap <- jolie.ggplot1.transparent$zeMap + label4Map2 + plot.fi.mean$zeMap_dots + plot.fi.mean$zeMap_text

    if (show.TI) {
      a5.JolieggMap.Fi.TI <- jolie.ggplot1.transparent$zeMap + label4Map2 + fi.TI + plot.fi.mean$zeMap_dots + plot.fi.mean$zeMap_text
      a5.JolieggMap.Fimean.TI <- jolie.ggplot1.transparent$zeMap_background + label4Map2 + fi.TI + plot.fi.mean$zeMap_dots + plot.fi.mean$zeMap_text
      if (printGraphs) {
        png('RowImapWithTI.png')
        print(a5.JolieggMap.Fi.TI)
        dev.off()

        png('RowImapWithTI(MeansOnly).png')
        print(a5.JolieggMap.Fimean.TI)
        dev.off()
      }
    }
    if (show.CI) {
      a5.JolieggMap.Fi.CI <- jolie.ggplot1.transparent$zeMap + label4Map2 + fi.CI + plot.fi.mean$zeMap_dots + plot.fi.mean$zeMap_text
      a5.JolieggMap.Fimean.CI <- jolie.ggplot1.transparent$zeMap_background + label4Map2 + fi.CI + plot.fi.mean$zeMap_dots + plot.fi.mean$zeMap_text
      if (printGraphs) {
        png('RowImapWithCI.png')
        print(a5.JolieggMap.Fi.CI)
        dev.off()

        png('RowImapWithCI(MeansOnly).png')
        print(a5.JolieggMap.Fimean.CI)
        dev.off()
      }
    }
    if (show.TI == TRUE & show.CI == TRUE) {
      a5.JolieggMap.Fi.CITI <- jolie.ggplot1.transparent$zeMap + label4Map2 + fi.TI + fi.CI + plot.fi.mean$zeMap_dots + plot.fi.mean$zeMap_text
      a5.JolieggMap.Fimean.CITI <- jolie.ggplot1.transparent$zeMap_background + label4Map2 + fi.TI + fi.CI + plot.fi.mean$zeMap_dots + plot.fi.mean$zeMap_text
      if (printGraphs) {
        png('RowImapWithCITI.png')
        print(a5.JolieggMap.Fi.CITI)
        dev.off()

        png('RowImapWithCITI(MeansOnly).png')
        print(a5.JolieggMap.Fimean.CITI)
        dev.off()
      }
    }




  }

  ### I-Cosines -----
  cos4I <- sqrt(resCA$ExPosition.Data$ri) *
    sign(resCA$ExPosition.Data$fi)
  ### 1.Iset - Circle of corr ----
  jolie.ggplot.I <-
    PTCA4CATA::createFactorMap(
      cos4I,
      col.points = col4I,
      col.labels = col4I,
      display.labels = display.labels.ind,
      constraints = list(
        minx = -1,
        miny = -1,
        maxx = 1 ,
        maxy = 1
      ),
      font.face = "italic"
    )
  # draw the circle
  a01.jolieggMap.I <- jolie.ggplot.I$zeMap +
    addCircleOfCor(color = "darkorchid4")
  if (printGraphs) {
    print(a01.jolieggMap.I)
  }
  #  Dot for the I-set, no arrows
  a02.jolieggMap.I <-
    jolie.ggplot.I$zeMap_background +
    jolie.ggplot.I$zeMap_text +
    addCircleOfCor(color = "darkorchid4") +
    jolie.ggplot.I$zeMap_dots  + label4Map2
  if (printTest) {
    print(a02.jolieggMap.I)
  }
  if (printGraphs) {
    png('CircleOfCorrISet.png')
    print(a02.jolieggMap.I)
    dev.off()
  }

  ### J-Cosines -----
  cos4J <- sqrt(resCA$ExPosition.Data$rj) *
    sign(resCA$ExPosition.Data$fj)

  jolie.ggplot.J <-
    PTCA4CATA::createFactorMap(
      cos4J,
      col.points = col4J,
      col.labels = col4J,
      display.labels = display.labels.var,
      constraints = list(
        minx = -1,
        miny = -1,
        maxx = 1 ,
        maxy = 1
      )
    )
  # draw the circle
  b1.jolieggMap.J <- jolie.ggplot.J$zeMap +
    addCircleOfCor(color = "darkorchid4") + label4Map2
  if (printGraphs) {
    png('J-CircleOfCorr_noArrow.png')
    print(b1.jolieggMap.J)
    dev.off()
  }
  #  Add some arrows
  # arrows <-
  #   addArrows(loadings.2, color = col4J)
  # b2.jolieggMap.J <-
  #   jolie.ggplot.J$zeMap_background +
  #   jolie.ggplot.J$zeMap_text +
  #   addCircleOfCor(color = "darkorchid4") +
  #   arrows + label4Map2
  # # print(b2.jolieggMap.J)
  # if (printGraphs) {
  #   png('J-CircleOfCorr.png')
  #   print(b2.jolieggMap.J)
  #   dev.off()
  # }
  ## I & J Circle -----
  b2.jolieggMap.IJ <- b1.jolieggMap.J +
    jolie.ggplot.I$zeMap_text +
    jolie.ggplot.I$zeMap_dots + label4Map2
  #print(b2.jolieggMap.IJ)
  if (printGraphs) {
    png('IJ-CircleOfCorr.png')
    print(b2.jolieggMap.IJ)
    dev.off()
  }

  # Fj factor map

  jolie.ggplot.J.fj <- PTCA4CATA::createFactorMap(
    Fj,
    constraints = constraintsIJ,
    col.points = col4J,
    col.labels = col4J,
    display.labels = display.labels.var
  )

  b3.jolieggMap.J.fj <-
    jolie.ggplot.J.fj$zeMap +
    label4Map2 +
    labs(title = 'The Column Map')
  if(printTest){
    print(b3.jolieggMap.J.fj)
  }

  joint.map <-
    jolie.ggplot1$zeMap + jolie.ggplot.J.fj$zeMap_dots + jolie.ggplot.J.fj$zeMap_text +
    label4Map2

  if(printTest){
    print(joint.map)
  }

  ### Contributions ----
  signed.ctrI <- resCA$ExPosition.Data$ci *
    sign(resCA$ExPosition.Data$fi)
  # if (!is.null(DESIGN)){
  #   signed.ctrI <- signed.ctrI[order(DESIGN),]
  # }
  ### ctr4I -----
  leStem <-  'I-Set'
  nfac  <- 1
  nind <- nrow(signed.ctrI)
  if (nind < 10){
    font.size.ctr <- 5
  }else {
    font.size.ctr = 5-(round((nind-10)/10)-1)
  }
  ctrI1 <- plotCtr(signed.ctrI, col4I,
                   nfac, stem = leStem,
                   font.size = font.size.ctr)

  nfac  <- 2
  ctrI2 <- plotCtr(signed.ctrI, col4I,
                   nfac, stem = leStem,
                   font.size = font.size.ctr)

  if (nind > max.n4bar){
    ctrI1$layers[[3]] <- NULL
    ctrI2$layers[[3]] <- NULL
  }

  # Contributions ----
  signed.ctrJ <- resCA$ExPosition.Data$cj *
    sign(resCA$ExPosition.Data$fj)
  ## plot ctr ----
  # plot contributions for component 1
  col4J.bar <- col4J
  #col4J.bar[3] <- '#2F7D1F'
  nfac  <- 1
  nvar <- nrow(signed.ctrJ)
  if (nvar < 10){
    font.size.ctr <- 5
  }else {
    font.size.ctr = 5-(round((nvar-10)/10)-1)
  }

  ctrJ1 <-
    plotCtr(signed.ctrJ, col4J.bar, nfac,
            font.size = font.size.ctr)
  nfac  <- 2
  ctrJ2 <-
    plotCtr(signed.ctrJ, col4J.bar, nfac,
            font.size = font.size.ctr)
  if (nvar > max.n4bar){
    ctrJ1$layers[[3]] <- NULL
    ctrJ2$layers[[3]] <- NULL
  }


  # return lists ----
  ## stat ----
  results.stats <- list(
    ExPosition.Data = resCA$ExPosition.Data,
    Plotting.Data = resCA$Plotting.Data
  )

  ## graphs ----
  results.graphs <- list(
    # cor and cov mat here
    #covariance = a5.02.covMap,
    #correlation = a5.02.correlationMap,
    chi2.residuals = a0.residuals,
    scree = a01.leScree,
    ctrI.1 = ctrI1,
    ctrI.2 = ctrI2,
    factorScoresI12 = a4.JolieggMap.2,
    factorScoresJ12 = b3.jolieggMap.J.fj,
    #factorScoresJ12.arrow = b3.jolieggMap.J.fj.arrow,
    cosineCircle4I12 = a02.jolieggMap.I,
    cosineCircleJ12  =  b1.jolieggMap.J,
    ctrJ.1 = ctrJ1,
    ctrJ.2 = ctrJ2,
    #cosineCircleArrowJ12  =  b2.jolieggMap.J,
    cosineCircleIJ12 = b2.jolieggMap.IJ,
    #loadings12 = b3.jolieggMap.J.Q,
    #loadings12.arrow = b3.jolieggMap.J.Q.arrow #,
    # biplot12 = e.JolieBiplot
    biplot12 = joint.map
  )

  description.graphs <- list(
    #covariance = "The Covariance Matrix Heat Map",
    #correlation = "The Correlation Matrix Heat Map",
    chi2.residuals = "Chi Square Residuals",
    scree = "The Eigenvalues Scree Plot",
    ctrI.1 = "Rows: Contributions Dimension 1",
    ctrI.2 = "Rows: Contributions Dimension 2",
    factorScoresI12 =  "Rows: Factor Scores 1*2",
    factorScoresJ12 = "Columns: Factor Scores 1*2",
    #factorScoresJ12.arrow = "Variables: Loadings as Inertia 1*2 (with arrows)",
    cosineCircle4I12 = "Rows: Cosine Circle 1*2",
    cosineCircleJ12  = "Columns: Cosine Circle 1*2",
    ctrJ.1 = "Columns: Contributions Dimension 1",
    ctrJ.2 = "Columns: Contributions Dimension 2",
    #cosineCircleArrowJ12  =  "Variables: Correlation Circle 1*2 (with arrows)",
    cosineCircleIJ12 = "Rows and Columns: Cosine Circle 1*2",
    #loadings12 = "Variables: Loadings as Weights  1*2",
    #loadings12.arrow = "Variables: Loadings as Weights  1*2 (with arrows)" #,
    biplot12 = "Row and Column Biplot: Factor Scores 1*2"
  )

  if (!is.null(DESIGN)){
    results.graphs$factorScoresI12design = a5.JolieggMap
    description.graphs$factorScoresI12design = "Rows with design: Factor Scores 1*2"
    if (show.TI) {
      results.graphs$factorScoresI12design.TI = a5.JolieggMap.Fi.TI
      results.graphs$meanfactorScoresI12design.TI = a5.JolieggMap.Fimean.TI
      description.graphs$factorScoresI12design.TI = "Rows with tolerance intervals: Factor Scores 1*2"
      description.graphs$meanfactorScoresI12design.TI = "Rows with tolerance intervals: Mean Factor Scores 1*2"
    }
    if (show.CI) {
      results.graphs$factorScoresI12design.CI = a5.JolieggMap.Fi.CI
      results.graphs$meanfactorScoresI12design.CI = a5.JolieggMap.Fimean.CI
      description.graphs$factorScoresI12design.CI = "Rows with confidence intervals: Factor Scores 1*2"
      description.graphs$meanfactorScoresI12design.CI = "Rows with confidence intervals: Mean Factor Scores 1*2"
    }
    if (show.TI == TRUE & show.CI == TRUE) {
      results.graphs$factorScoresI12design.CITI = a5.JolieggMap.Fi.CITI
      results.graphs$meanfactorScoresI12design.CITI = a5.JolieggMap.Fimean.CITI
      description.graphs$factorScoresI12design.CITI = "Rows with CI and TI: Factor Scores 1*2"
      description.graphs$meanfactorScoresI12design.CITI = "Rows with CI and TI: Mean Factor Scores 1*2"
    }
  }



  results <- list(
    results.stats = results.stats,
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
