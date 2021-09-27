#' App structure
#' @name restrictionSimulation
#' @description Performance of in-silico restriction simulation from reference genome
#' @param dnaseq File genomic data
#' @param enzyme_db Enzyme DB
#' @param min.size min size of length fragment
#' @param max.size max size of length fragment
#' @param type_analysis definition of analysis
#' @param enzyme_selection if selection enzymes, specify
#' @param nb_repeat number of repetitions
#' @param use_output Create output optional
#' @param outputDir output path
#' @import dplyr, XML, tidyverse, SimRAD, seqinr, data.table, phylotools, biomartr
#' @return results
#' @author Alberto Rodriguez-Izquierdo, 2021




restrictionSimulation <- function(dnaseq,
                                  enzyme_db,
                                  min.size,
                                  max.size,
                                  type_analysis,
                                  enzyme_selection = NULL,
                                  nb_repeat=NULL,
                                  use_output = NULL,
                                  outputDir=NULL,
                                  wd=NULL){

  #############---------Check QC data entry--------------#############

  ValDnaseq          <- validateCharacter(dnaseq)

  if (!isTRUE(str_detect(ValDnaseq, '.fasta', negate = FALSE))){

    if (!isTRUE(str_detect(ValDnaseq, '.fa', negate = FALSE))){

      print('File is not a FASTA file! Please check the file (.fa, .fasta)')

      stop(restrictionSimulation)

    }
  }

  dnaseq <- ValDnaseq

  ##############--------------------------------------------#########

  Valenzyme_db <- validateCharacter(Valenzyme_db)

  enzyme_db <- Valenzyme_db


  ##############--------------------------------------------#########

  Valminsize <- validateNumber(min.size)

  if (!is.integer(Valminsize)){

    print('Please input a number')

    stop(restrictionSimulation)

  }

  min.size <- Valminsize

  ##############--------------------------------------------#########

  Valmaxsize <- validateNumber(max.size)

  if (!is.integer(Valmaxsize)){

    print('Please input a number')

    stop(restrictionSimulation)

  }

  max.size <- Valmaxsize


  ##############--------------------------------------------#########

  ValtypeAnalysis <- validateCharacter(type_analysis)

  if (ValtypeAnalysis =='finding'){

    type_analysis <- ValtypeAnalysis

    }else if(ValtypeAnalysis == 'combination'){

        if (!is.null(Valenzyme_selection)){

          if(!length(enzyme_selection)>1){

            Valenzyme_selection <- validateCharacter(enzyme_selection)

          }else{

              valenzyme_selection <- enzyme_selection

          }

          for (enz in Valenzyme_selection){

            if (!str_detect(enzyme_db, enz)){

              print(paste0('Enzyme ', enz, ' not detected in database. Please try other enzyme'))

              stop(restrictionSimulation)

            }

          }

          enzyme_selection <- Valenzyme_selection

        }else{

          print ('Please input the name of the restriction enzyme/s')

          stop(restrictionSimulation)

        }

    }else if(ValtypeAnalysis == 'replicate'){


      if (!is.null(Valenzyme_selection)){

        if(!length(enzyme_selection)>1){

          Valenzyme_selection <- validateCharacter(enzyme_selection)

        }else{

          valenzyme_selection <- enzyme_selection

        }

        for (enz in Valenzyme_selection){

          if (!str_detect(enzyme_db, enz)){

            print(paste0('Enzyme ', enz, ' not detected in database. Please try other enzyme'))

            stop(restrictionSimulation)

          }

        }

        enzyme_selection <- Valenzyme_selection

      }else{

        print ('Please input the name of the restriction enzyme/s')

        stop(restrictionSimulation)

      }

      ValNbRepeat <- validateNumber(nb_repeat)

      if (!is.integer(ValNbRepeat)){

        print('Number of repeats must be integer.')

        stop(restrictionSimulation)

      }

    }else{

      print('Please choose the type of analysis ("finding","combination" or "replicate")')

      stop(restrictionSimulation)


    }




  #############--------Starting function--------------#############

  load(file=enzyme_db)

  enzyme.db     <- data.frame(lapply(enzyme.db, as.character), stringsAsFactors = FALSE)

  root          <- wd

  configFile    <- outputDir

  #Starting restriction simulation
  if (type_analysis == 'finding'){

    print('#########------------Finding the best enzyme----------########')

    for (enzymes in 1:nrow(enzyme.db)){

      row       <- enzyme.db[enzymes,]

      sequences <- as.character(row$restriction.site.sequences)

      eval(parse(text=paste0('simseq', row$enzyme,'.dig <- insilico.digest(ref.DNAseq(dnaseq), cut_site_5prime1="', sequences,'", cut_site_3prime1="", verbose=TRUE)')))

      eval(parse(text=paste0('size.select',row$enzyme,' <- size.select (simseq', row$enzyme,'.dig,min.size = min.size, max.size = max.size, graph = FALSE, verbose= TRUE)')))

      if (isTRUE(use_output)){

        if (eval(parse(text=paste0('!length(size.select', row$enzyme,') == 0')))){

          eval(parse(text=paste0('graph_generator(simseq',row$enzyme,'.dig, row$enzyme, min.size, max.size, configFile, root)')))

        }

      }

      #building exiting file

      if (!exists('results')){

        #rm(list=ls())
        results <- data.frame(row)

        eval(parse(text=paste0('All_digestions_count <- length(simseq', row$enzyme,'.dig)')))

        eval(parse(text=paste0('Sel_digestions_count <- length(size.select', row$enzyme,')')))

        results <- cbind (results, All_digestions_count)

        results <- cbind (results, Sel_digestions_count)

        #eval(parse(text=paste0('results$Number_valid_Rest_fragments <- cbind(results, length(size.select',row$enzyme,'))')))


      }else{

        results1 <- data.frame(row)

        eval(parse(text=paste0('All_digestions_count <- length(simseq', row$enzyme,'.dig)')))

        eval(parse(text=paste0('Sel_digestions_count <- length(size.select', row$enzyme,')')))

        results1 <- cbind (results1, All_digestions_count)

        results1 <- cbind (results1, Sel_digestions_count)

        results  <- rbind (results, results1)

      }

      #cleaning memory

      eval(parse(text=paste0('rm(simseq',row$enzyme,'.dig)')))

      eval(parse(text=paste0('rm(size.select', row$enzyme,')')))
    }

    results[order(-results$Sel_digestions_count),]


########--------------------Use Combination--------------------##########



  }else if(type_analysis == 'combination'){

    print('#########------------Combination of restriction enzymes----------########')

    for (x in enzyme_selection){

      if (exists('enzyme_selectiondb')==FALSE){

        enzyme_selectiondb    <- filter(enzyme.db, enzyme == x)

      }else{

        enzyme_selectiondb_x  <- filter(enzyme.db, enzyme == x)

        enzyme_selectiondb    <- rbind (enzyme_selectiondb, enzyme_selectiondb_x)

      }

    }

    for (enzymes in enzyme_selection){

      for (enzymes_2 in enzyme_selection){

        if (!enzymes_2 == enzymes){

          row1                <- filter(enzyme_selectiondb, enzyme == enzymes)

          row2                <- filter(enzyme_selectiondb, enzyme == enzymes_2)

          sequences1          <- as.character(row1$restriction.site.sequences)

          sequences2          <- as.character(row2$restriction.site.sequences)

          eval(parse(text=paste0('simseq', row1$enzyme,'vs', row2$enzyme,'.dig <- insilico.digest(ref.DNAseq(dnaseq), cut_site_5prime1="', sequences1,'", cut_site_3prime1="", cut_site_5prime2="', sequences2,'", cut_site_3prime2="", verbose=TRUE)')))

          eval(parse(text=paste0('size.select',row1$enzyme,'vs',row2$enzyme,' <- size.select (simseq', row1$enzyme,'vs',row2$enzyme,'.dig,min.size = min.size, max.size = max.size, graph = FALSE, verbose= TRUE)')))

          if (isTRUE(use_output)){

            if (eval(parse(text=paste0('!length(size.select', row1$enzyme,'vs',row2$enzyme,') == 0')))){

              eval(parse(text=paste0('graph_generator(simseq',row1$enzyme,'vs',row2$enzyme,'.dig, "', row1$enzyme,'vs', row2$enzyme,'", min.size, max.size, configFile, root)')))

            }
          }


          #building exiting file

          if (!exists('results')){

            #rm(list=ls())
            results           <- data.frame(row1)

            results           <- cbind(results, data.frame(row2))

            eval(parse(text=paste0('All_digestions_count <- length(simseq', row1$enzyme,'vs', row2$enzyme,'.dig)')))

            eval(parse(text=paste0('Sel_digestions_count <- length(size.select', row1$enzyme,'vs',row2$enzyme,')')))

            results           <- cbind (results, All_digestions_count)

            results           <- cbind (results, Sel_digestions_count)

            #eval(parse(text=paste0('results$Number_valid_Rest_fragments <- cbind(results, length(size.select',row$enzyme,'))')))


          }else{

            results1          <- data.frame(row1)

            results1          <- cbind(results1, data.frame(row2))

            eval(parse(text=paste0('All_digestions_count <- length(simseq', row1$enzyme,'vs', row2$enzyme,'.dig)')))

            eval(parse(text=paste0('Sel_digestions_count <- length(size.select', row1$enzyme,'vs', row2$enzyme,')')))

            results1          <- cbind (results1, All_digestions_count)

            results1          <- cbind (results1, Sel_digestions_count)

            results           <- rbind (results, results1)

          }

          #cleaning memory

          eval(parse(text=paste0('rm(simseq',row1$enzyme,'vs', row2$enzyme,'.dig)')))

          eval(parse(text=paste0('rm(size.select', row1$enzyme,'vs', row2$enzyme,')')))

          }
        }
      }

    results[order(-results$Sel_digestions_count),]


#########-----------Use replicate------------##########


  }else if (type_analysis == 'replicate'){

    print('#########------------Replication restriction selecting enzymes----------########')

    for(nrepeat in 1:nb_repeat){

      enzymes         <- configFile$parameters$replicate_enzyme$enzyme_selection

      row             <- filter(enzyme.db, enzyme == enzymes)

      sequences       <- as.character(row$restriction.site.sequences)

      eval(parse(text=paste0('simseq', row$enzyme,'.dig <- insilico.digest(ref.DNAseq(dnaseq), cut_site_5prime1="', sequences,'", cut_site_3prime1="", verbose=TRUE)')))

      eval(parse(text=paste0('size.select',row$enzyme,' <- size.select (simseq', row$enzyme,'.dig,min.size = min.size, max.size = max.size, graph = FALSE, verbose= TRUE)')))

      if (isTRUE(use_output)){

        if (eval(parse(text=paste0('!length(size.select', row$enzyme,') == 0')))){

          eval(parse(text=paste0('graph_generator(simseq',row$enzyme,'.dig, row$enzyme, min.size, max.size, configFile, root, nrepeat, use_replicate)')))

        }

      }

      #building exiting file

      if (!exists('results')){

        #rm(list=ls())
        results       <- data.frame(row)

        eval(parse(text=paste0('All_digestions_count <- length(simseq', row$enzyme,'.dig)')))

        eval(parse(text=paste0('Sel_digestions_count <- length(size.select', row$enzyme,')')))

        results       <- cbind (results, All_digestions_count)

        results       <- cbind (results, Sel_digestions_count)

        #eval(parse(text=paste0('results$Number_valid_Rest_fragments <- cbind(results, length(size.select',row$enzyme,'))')))


      }else{

        results1        <- data.frame(row)

        eval(parse(text=paste0('All_digestions_count <- length(simseq', row$enzyme,'.dig)')))

        eval(parse(text=paste0('Sel_digestions_count <- length(size.select', row$enzyme,')')))

        results1        <- cbind (results1, All_digestions_count)

        results1        <- cbind (results1, Sel_digestions_count)

        results         <- rbind (results, results1)

      }

      #cleaning memory

      eval(parse(text=paste0('rm(simseq',row$enzyme,'.dig)')))

      eval(parse(text=paste0('rm(size.select', row$enzyme,')')))
    }

    results[order(-results$Sel_digestions_count),]

    }

  if (isTRUE(use_output)){

    outputGeneration(results,root)

  }

  return(results)
}


#' App structure
#' @name graph_generator
#' @description Graph generation from in-silico digestion.
#' @param sequences Sequences to be graphed
#' @param name_enzyme name of the enzyme
#' @param min.size min size to plot
#' @param max.size max size to plot
#' @param configFile configFile
#' @param nrepeat number of repetitions
#' @param use_replicate if use replicate
#' @import dplyr, XML, tidyverse, SimRAD, seqinr, data.table, phylotools, biomartr
#' @author Rodriguez-Izquierdo, Alberto, 2021
#'

graph_generator <- function(sequences, name_enzyme, min.size, max.size, configFile, root, nrepeat=NULL, use_replicate=NULL){

  if (!exists('nrepeat')){

    nrepeat <- NULL

  }

  ssel <- sequences[width(sequences) < max.size & width(sequences) > min.size]

  dirOutput <- paste0(root,'output/',configFile,'/', name_enzyme)

  if (!dir.exists(paste0(root,'output/',configFile))){

    dir.create(paste0(root,'output/',configFile))

    }

  if (!dir.exists(dirOutput)){

    dir.create(dirOutput)

  }

  if(!isTRUE(use_replicate)){

    unlink(dirOutput, recursive = TRUE)

    dir.create(dirOutput)

  }
  if (!is.null(nrepeat)){

    tiff(filename=paste0(dirOutput, '/', name_enzyme,nrepeat,'_histogram.tiff'), units= 'in', width=5, height= 5, res=300)

    bk <- hist(width(sequences), breaks = length(sequences)/20,
               plot = FALSE)$breaks

    hist(width(sequences),
         border = "grey75",
         col = "grey75",
         breaks = bk,
         main = name_enzyme,
         xlab = "Locus size (bp)",
         ylab = "Number of loci")

    hist(width(ssel),
         border = "red",
         col = "red",
         add = TRUE,
         breaks = bk)

    text(mean(c(min.size, max.size)),
         max(hist(width(ssel),breaks = bk, plot = FALSE)$counts),
         pos = 4,
         labels = paste(length(ssel)," loci between ", min.size, " and ", max.size, " bp", sep = ""),
         col = "red",
         cex = 0.9,
         font = 2)

    dev.off()


    sequences_df          <- data.frame(ssel)
    sequences_df$ID       <- seq.int(nrow(sequences_df))

    write.fasta(as.list(sequences_df$ssel), sequences_df$ID, file.out=paste0(dirOutput, '/', name_enzyme,nrepeat, '.fasta'))

  }else{

    tiff(filename=paste0(dirOutput, '/', name_enzyme,'_histogram.tiff'), units= 'in', width=5, height= 5, res=300)

    bk <- hist(width(sequences), breaks = length(sequences)/20,
               plot = FALSE)$breaks

    hist(width(sequences),
         border = "grey75",
         col = "grey75",
         breaks = bk,
         main = name_enzyme,
         xlab = "Locus size (bp)",
         ylab = "Number of loci")

    hist(width(ssel),
         border = "red",
         col = "red",
         add = TRUE,
         breaks = bk)

    text(mean(c(min.size, max.size)),
         max(hist(width(ssel),breaks = bk, plot = FALSE)$counts),
         pos = 4,
         labels = paste(length(ssel)," loci between ", min.size, " and ", max.size, " bp", sep = ""),
         col = "red",
         cex = 0.9,
         font = 2)

    dev.off()


    sequences_df          <- data.frame(ssel)
    sequences_df$ID       <- seq.int(nrow(sequences_df))

    write.fasta(as.list(sequences_df$ssel), sequences_df$ID, file.out=paste0(dirOutput, '/', name_enzyme,'.fasta'))

  }
}
