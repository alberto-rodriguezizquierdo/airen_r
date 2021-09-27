#' App structure
#' @name outputGeneration
#' @description Designing output
#' @param results Input of results
#' @param root Working directory
#' @param configFile configFile
#'
#'
#' @author Alberto Rodriguez-Izquierdo, 2021
#'


outputGeneration <- function(results,root){

  dirOutput         <- paste0(root, 'output/results/')

  if (!dir.exists(dirOutput)){

    dir.create(dirOutput)

  }else{

    unlink(dirOutput, recursive = TRUE)

    dir.create(dirOutput)

  }

  write.table(results, file=paste0(dirOutput, 'results_RE.csv'), sep=';')


}
