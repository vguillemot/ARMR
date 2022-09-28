library(ExPosition)
library(InPosition)
library(ggplot2)
#devtools::install_github('HerveAbdi/PTCA4CATA')
library(PTCA4CATA)
#devtools::install_github('HerveAbdi/data4PCCAR')
library(data4PCCAR)
# **** Get the data ----
data("sixAuthorsPunctuated",
     package = 'data4PCCAR')
## Active set ----
# The active data set
X <- as.matrix(sixAuthorsPunctuated$df.active)

# **** CA Analysis
# Plain symmetric
resCA.sym <- epCA(X, symmetric = TRUE, graphs = FALSE)
# column asymmetric
resCA.asym <- epCA(X, symmetric = FALSE, graphs = FALSE)

res_caplot_sym <- CAplot(resCA.sym, X,
                         display.labels.ind = TRUE, save2pptx = TRUE,
                         show.CI = TRUE, show.TI = TRUE,
                         title4pptx = "CA Results Symmetric")


res_caplot_asym <- CAplot(resCA.asym, X,
                         display.labels.ind = TRUE, save2pptx = TRUE,
                         show.CI = TRUE, show.TI = TRUE,
                         title4pptx = "CA Results Asymmetric")


design = c(rep("G1", 3), rep("G2", 3))

res_caplot_sym_des <- CAplot(resCA.sym, X, DESIGN = design,
                         display.labels.ind = TRUE, save2pptx = TRUE,
                         show.CI = TRUE, show.TI = TRUE,
                         title4pptx = "CA Results Symmetric Design")


res_caplot_asym_des <- CAplot(resCA.asym, X, DESIGN = design,
                          display.labels.ind = TRUE, save2pptx = TRUE,
                          show.CI = TRUE, show.TI = TRUE,
                          title4pptx = "CA Results Asymmetric Design")


excel_output <- saveCAResults2xlsx(res_caplot_sym$results.stats, X, NULL, "CA Results Symmetric.xlsx")
excel_output2 <- saveCAResults2xlsx(res_caplot_sym_des$results.stats, X, design, "CA Results Symmetric Design.xlsx")


res_caplot_sym_inf <- CAplotInference(resCA.sym, X, fast = TRUE,
                                      save2pptx = TRUE, title4pptx = "CA Inference Symmetric.xlsx")
