#' savePCAResults2xlsx
#'
#' Save an Excel Spreadsheet based on the \code{results.stats}
#' from [PCAplot()]
#'
#' @param results.stats The same \code{results.stats} from \code{PCAplot()}
#' @param data_pca The data passed to \code{epPCA()}
#' @param design The design variable passed to \code{epPCA()}
#' @param file2Save.xlsx File name / path to save the Excel sheet
#'
#' @author Luke Moraglia
#'
#' @return a list of data written to the Excel sheet, with class \code{"save2xlsx"}.
#' @export
savePCAResults2xlsx <- function(results.stats,
                                data_pca,
                                design = NULL,
                                file2Save.xlsx = "PCA Results.xlsx"){

  ED <- results.stats$ExPosition.Data

  eigs_tbl <- cbind("Eigenvalues" = ED$eigs, "Percent Inertia" = ED$t)

  list_of_data <- list("PCA Input Data" = as.data.frame(data_pca),
                       "Obs Factor Scores" = as.data.frame(ED$fi),
                       "Obs Signed Contributions" = as.data.frame(ED$ci * sign(ED$fi)),
                       "Var Loadings as Inertia" = as.data.frame(results.stats$loadings.as.inertia),
                       "Var Loadings as Correlation" = as.data.frame(results.stats$loadings.as.correlation),
                       "Var Loadings as Weights" = as.data.frame(results.stats$loadings.as.weights),
                       "Var Signed Contributions" = as.data.frame(ED$cj * sign(ED$fj)),
                       "Eigenvalues" = as.data.frame(eigs_tbl)
  )
  if(!is.null(design)){
    list_of_data <- append(list_of_data, list("Design Variable" = as.data.frame(design)), 1)
  }


  list_to_write <- lapply(list_of_data, function(x) cbind(" " = rownames(x), x))

  res_saveResultsList2xlsx <- saveResultsList2xlsx(list_to_write,
                                                   file2Save.xlsx)

  return(res_saveResultsList2xlsx)


}

#' saveCAResults2xlsx
#'
#' Save an Excel Spreadsheet based on the \code{results.stats}
#' from [CAplot()]
#'
#' @param results.stats The same \code{results.stats} from \code{CAplot()}
#' @param data_ca The data passed to \code{epCA()}
#' @param design The design variable passed to \code{epCA()}
#' @param file2Save.xlsx File name / path to save the Excel sheet
#'
#' @author Luke Moraglia
#'
#' @return a list of data written to the Excel sheet, with class \code{"save2xlsx"}.
#' @export
saveCAResults2xlsx <- function(results.stats,
                                data_ca,
                                design = NULL,
                                file2Save.xlsx = "CA Results.xlsx"){

  ED <- results.stats$ExPosition.Data

  eigs_tbl <- cbind("Eigenvalues" = ED$eigs, "Percent Inertia" = ED$t)

  list_of_data <- list("CA Input Data" = as.data.frame(data_ca),
                       "Row Factor Scores" = as.data.frame(ED$fi),
                       "Row Signed Contributions" = as.data.frame(ED$ci * sign(ED$fi)),
                       "Column Factor Scores" = as.data.frame(ED$fj),
                       "Column Signed Contributions" = as.data.frame(ED$cj * sign(ED$fj)),
                       "Eigenvalues" = as.data.frame(eigs_tbl)
  )
  if(!is.null(design)){
    list_of_data <- append(list_of_data, list("Design Variable" = as.data.frame(design)), 1)
  }


  list_to_write <- lapply(list_of_data, function(x) cbind(" " = rownames(x), x))

  res_saveResultsList2xlsx <- saveResultsList2xlsx(list_to_write,
                                                   file2Save.xlsx)

  return(res_saveResultsList2xlsx)


}

#' saveResultsList2xlsx
#'
#' save a named list of data frames to sheets in an .xlsx file
#'
#' @param list2Save a named list of data frames to save
#' @param file2Save.xlsx file name to save the .xlsx file
#'
#' @author Luke Moraglia
#'
#' @return a list of data written to the Excel sheet, with class \code{"save2xlsx"}.
#' @export
saveResultsList2xlsx <- function(list2Save,
                                 file2Save.xlsx = "Results.xlsx"){
  xlsx.type = 'xlsx'
  if(tools::file_ext(file2Save.xlsx) != xlsx.type){
    file2Save.xlsx <- paste0(file2Save.xlsx, '.', xlsx.type)
  }
  if (file.exists(file2Save.xlsx)){# if file already exists: rename it
    LaDate = substr(as.POSIXlt(Sys.time()),1,10)
    OldFilename = sub(paste0('[.]', xlsx.type),
                      paste0('-',LaDate,'.',xlsx.type),file2Save.xlsx)
    file.rename(from = file2Save.xlsx, to = OldFilename)
    warning(paste0("File: ",file2Save.xlsx,' already exists.\n',
                   ' Oldfile has been renamed: ', OldFilename),
            call. = FALSE)
  }

  writexl::write_xlsx(list2Save,
                      path = file2Save.xlsx,
                      col_names = TRUE,
                      format_headers = TRUE)

  return.list <- structure(
    list(
      listOfSavedData = list2Save,
      nameOfFile = file2Save.xlsx
    ),  class = "save2xlsx")


  return(return.list)
}
