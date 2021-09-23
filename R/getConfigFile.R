#' App structure
#' @name getConfigFile
#' @param root
#' @description App to start geoapp
#'
#' @import logging
#' @return configFile
#'
#' @author Alberto Rodriguez-Izquierdo, 2021

getConfigFile <- function(root){

  #loginfo('Reading configFile...')
  configFile <- readConfigFile(root)


  #####-------------------Validations----------------------######

  #loginfo('Validating configFile...')


  #------------------Validation nodes----------------------#

  configFile <- nodesValidation(configFile)
  

  ######--------------Validation nodes content --------------######

#  error <- FALSE

  #Validate generalParameters

###-----------------------Validation data-----------------------###
  
  validateDataPath          <- validateCharacter(configFile$data$dataPath)

  configFile$data$dataPath  <- validateDataPath

  validateGenome            <- validateCharacter(configFile$data$genome)

  configFile$data$genome    <- validateGenome

  validateEnzymeDB          <- validateCharacter(configFile$data$enzyme_db)

  configFile$data$enzyme_db <- validateEnzymeDB
  
  
  ###----------------------Validation parameters---------------###

  validateTypeAnalysis                                              <- validateCharacter(configFile$parameters$typeAnalysis)
  
  configFile$parameters$typeAnalysis                                <- validateTypeAnalysis
  
  validateIfRestriction                                             <- validateCharacter(configFile$parameters$ifRestriction)
  
  configFile$parameters$ifRestriction                               <- validateIfRestriction
  
  validateEnzymeSelection                                           <- validateCharacter(configFile$parameters$combination$enzyme_selection)
        
  validateEnzymeSelection                                           <- unlist(strsplit(validateEnzymeSelection, ','))
        
  configFile$parameters$combination$enzyme_selection                <- validateEnzymeSelection
        
  validateUseRepEnzyme                                              <- validateCharacter(configFile$parameters$replicate_enzyme$use_replicate)
  
  configFile$parameters$replicate_enzyme$use_replicate              <- validateUseRepEnzyme
  
  validateRepEnzymeSel                                              <- validateCharacter(configFile$parameters$replicate_enzyme$enzyme_selection)
  
  configFile$parameters$replicate_enzyme$enzyme_selection           <- validateRepEnzymeSel
  
  validateRepEnzymeNb                                               <- validateNumber(configFile$parameters$replicate_enzyme$nb_repeat)
  
  configFile$parameters$replicate_enzyme$nb_repeat                  <- validateRepEnzymeNb
  
  validateRandomNb                                                  <- validateNumber(configFile$parameters$random_genome_fragmentation$nb_fragments)
  
  configFile$parameters$random_genome_fragmentation$nb_fragments    <- validateRandomNb
  
  validateRandomRepe                                                <- validateNumber(configFile$parameters$random_genome_fragmentation$nb_repeat)
  
  configFile$parameters$random_genome_fragmentation$nb_repeat       <- validateRandomRepe
  
  validateUseCalculate                                              <- validateCharacter(configFile$parameters$calculatePosition$use_calculate)
  
  configFile$parameters$calculatePosition$use_calculate             <- validateUseCalculate
  
  validateAlignmentPath                                             <- validateCharacter(configFile$parameters$calculatePosition$alignment_path)
  
  configFile$parameters$calculatePosition$alignment_path            <- validateAlignmentPath
  
  validateGffFile                                                   <- validateCharacter(configFile$parameters$calculatePosition$gffFile)
  
  configFile$parameters$calculatePosition$gffFile                   <- validateGffFile
  
  validateOutputPath                                                <- validateCharacter(configFile$parameters$calculatePosition$outputPath)
  
  configFile$parameters$calculatePosition$outputPath                <- validateOutputPath

  validateCategory                                                  <- validateCharacter(configFile$parameters$calculatePosition$category)
  
  configFile$parameters$calculatePosition$category                  <- validateCategory

  validateMinSize                                                   <- validateNumber(configFile$parameters$min.size)
                
  configFile$parameters$min.size                                    <- validateMinSize
                
  validateMaxSize                                                   <- validateNumber(configFile$parameters$max.size)
                  
  configFile$parameters$max.size                                    <- validateMaxSize
  
  ###-----------------------Validate Output--------------------###
  
  validateUseOutput                      <- validateCharacter(configFile$output$use_output)
  
  configFile$output$use_output         <- validateUseOutput
  
  validateOutput                      <- validateCharacter(configFile$output$outputDir)
  
  configFile$output$outputDir         <- validateOutput
  

  return (configFile)
}


#' App structure
#' @name readConfigFile
#' @param root
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodriguez-Izquierdo, 2021

readConfigFile <- function(root){

  require(XML)

  myDirectoryConfigFile <- paste0(root, 'config/')

  if (!file.exists(paste0(myDirectoryConfigFile,"configFile.xml"))){

    #log_error <- 'File path does not exist! Please review the path'
    #logerror(log_error)

    #stop(log_error)
  }

  configFile <- xmlParse(file = paste0(myDirectoryConfigFile,"configFile.xml"))
  configFile <- xmlToList(configFile)
#  configFile <- data.frame(configFile)

  return (configFile)
}



#' @name nodesValidation
#' @param configFile
#' @return configFile
#'
#' @author Alberto Rodriguez-Izquierdo, 2021


nodesValidation <- function(configFile){

  #Building list with principal and secondary nodes for validation
  
  principalNodes              <- c('data','parameters','output')

  dataNodes                   <- c('dataPath','genome', 'enzyme_db')

  parametersNodes             <- c('typeAnalysis',
                                   'ifRestriction',
                                   'combination',
                                   'replicate_enzyme',
                                   'random_genome_fragmentation',
                                   'calculatePosition',
                                   'max.size', 
                                   'min.size')
  
  combinationNodes            <- c('enzyme_selection')
  
  replicateNodes              <- c('enzyme_selection', 'nb_repeat')
  
  randomGenNodes              <- c('nb_fragments', 'nb_repeat')
  
  calculatePositionNodes      <- c('alignment_path','gffFile','outputPath', 'category')

  outputNodes                 <- c('use_output','outputDir')

  #Validation principal nodes

  generalParametersNodes        <- validateConfigNodes(principalNodes, configFile)

  ValDataNodes                  <- validateConfigNodes(dataNodes, configFile$data)

  ValParametersNodes            <- validateConfigNodes(parametersNodes, configFile$parameters)

  ValCombinationNodes           <- validateConfigNodes(combinationNodes, configFile$parameters$combination)
  
  ValReplicateNodes             <- validateConfigNodes(replicateNodes, configFile$parameters$replicate_enzyme)
  
  valRandomGenNodes             <- validateConfigNodes(randomGenNodes, configFile$parameters$random_genome_fragmentation)
  
  valCalculatePositionNodes     <- validateConfigNodes(calculatePositionNodes, configFile$parameters$calculatePosition)
  
  ValOutputNodes                <- validateConfigNodes(outputNodes, configFile$output)


  return(configFile)
  #loginfo('Nodes validation success!')
}




#' App structure
#' @name validateNodes
#' @param configFile
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodríguez Izquierdo


validateConfigNodes <- function (nodes, configFile){
  
  if (!is.list(nodes)){
    for (x in nodes){
      if (!(x %in% names(configFile))){
        if (!(x %in% configFile)){

#         log_error <- paste0('Node ',x,' does not exist! Please review configFile.xml')
#         logerror(log_error)
          print(paste0('Warning: Node ', nodes,' does not exist'))
#         stop(log_error)
        }
      }
    }
  }else{
    for (x in names(nodes)){
      if (!(x %in% names(configFile))){

        #       log_error <- paste0('Node ',x,' does not exist! Please review configFile.xml')
        #       logerror(log_error)
        print(paste0('Warning: Node ', nodes,' does not exist'))
        #       stop(log_error)
      }
    }
  }
}


#' App structure
#' @name ValidateCharacter
#' @param configFile
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodríguez Izquierdo


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
#' @param configFile
#' @description App to start geoapp
#'
#' @import logging, XML
#' @return configFile
#'
#' @author Alberto Rodríguez Izquierdo



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

