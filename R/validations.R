#' App structure
#' @name ValidateCharacter
#' @param configFile configFile
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodríguez Izquierdo, 2021


validateCharacter <- function(configFile){

  if (!is.null(configFile)){

    if(configFile == "TRUE"){

      configFile <- TRUE

    }else if(configFile == "FALSE"){

      configFile <- FALSE
    }
  }

  return(configFile)

}


#' App structure
#' @name ValidateNumber
#' @param configFile configFile
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodríguez Izquierdo, 2021



validateNumber <- function(configFile){



  #  Error <- FALSE

  if (is.null(configFile)){

    #    log_error <- paste0('Value ', configFile, ' is null. Please check configFile')
    #    logerror(log_error)
    print('Warning: Node is null')
    #    Error <- TRUE
  }else if (!as.numeric(configFile)){

    #    log_error <- paste0('Value ', configFile, ' is not a number. Please check configFile')
    #    logerror(log_error)
    print('Warning: Node is not numeric char')
    #    Error <- TRUE
  }else{
    configFile <- as.numeric(configFile)
    #  if (any(Error == TRUE)){
    #    error <- Error
    #    stop('Not possible to validate, please check configFile format')
    #  }
  }
  return(configFile)
}

