#' @name SimRADoptionsApp
#' @param root
#' @import dplyr, SimRAD, seqinr, XML
#' @author Alberto Rodriguez-Izquierdo, 2021

SimRADoptionsApp <- function(root){

  
  configFile <- getConfigFile(root)

  #########---------Restriction Enzyme simulation----------#########
  if (configFile$parameters$typeAnalysis == 'restriction'){
    
    resultsRestriction <- restrictionSimulation(dnaseq = paste0(configFile$data$dataPath,'/',configFile$data$genome), 
                                                enzyme_db = paste0(configFile$data$dataPath,'/',configFile$data$enzyme_db), 
                                                min.size = configFile$parameters$min.size, 
                                                max.size = configFile$parameters$max.size,
                                                type_analysis = configFile$parameters$ifRestriction,
                                                enzyme_selection = configFile$parameters$combination$enzyme_selection,
                                                nb_repeat=configFile$parameters$replicate_enzyme$nb_repeat,
                                                use_output = configFile$output$use_output,
                                                outputDir = configFile$output$outputDir,
                                                wd=root)
    
    
    
  }
  
  ##########---------Random genome fragmentation------------##########  
  
  if (configFile$parameters$typeAnalysis == 'random'){
    
    for(x in 1:configFile$parameters$random_genome_fragmentation$nb_repeat){
      
      eval(parse(text=paste0('results', x,' <- randomGenomeFragmentation(configFile$data$genome, configFile$parameters$min_size, configFile$parameters$max_size, configFile$parameters$random_genome_fragmentation$nb_fragments)')))
      
      eval(parse(text=paste0("write.fasta(as.list(results",x,"$fragmentRes), results",x,"$V1, file.out=paste0(root,'/output/random_genome_fragmentation/results',x,'.fasta'))")))
      
      eval(parse(text=paste0('rm(results',x,')')))
      
    }
    
  }
  
  if (configFile$parameters$typeAnalysis == 'position'){
    
    calculatePositionFragment(count_dir=configFile$parameters$calculatePosition$outputPath,
                              alignment_path=configFile$parameters$calculatePosition$alignment_path,
                              gffPath=configFile$parameters$calculatePosition$gffFile,
                              category=configFile$parameters$calculatePosition$category)
    
  }
  


  print('App finished successfully!')
}



