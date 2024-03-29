#' @title CalculatePos
#' @name calculatePositionFragment
#'
#' @description Function that calculates the position of the founded fragments in sequencing.
#' @param count_dir Output path directory
#' @param alignment_path SAM/BAM files with all alignment positiions
#' @param gffPath GTF/GFF path
#' @param category Category to calculate position.
#' @import stringr
#' @import dplyr
#' @import data.table
#' @import biomartr
#' @import tidyverse
#'
#' @export
#'
#' @author Alberto Rodriguez-Izquierdo



calculatePositionFragment <- function(count_dir,
                                      alignment_path,
                                      gffPath,
                                      category){

#######----------QC entry---------#########

  ValGFF <- validateCharacter(gffPath)

  if(!isTRUE(str_detect(ValGFF, '.gff'))){

    print('Please enter a GFF file')

    stop(calculatePositionFragment)

  }

  gffPath <- ValGFF

  gfffile <- read_gff(gffPath)

  ######------------------------######

  ValAlignmentPath <- validateCharacter(alignment_path)

  if(!dir.exists(ValAlignmentPath)){

    print('Directory path does not exist! Please check the path.')

    stop(calculatePositionFragment)

  }

  alignment_path <- ValAlignmentPath


  ######-------------------------######

  ValCountDir <- validateCharacter(count_dir)

  if(!dir.exists(ValCountDir)){

    print('Directory path does not exist! Please check the path.')

    stop(calculatePositionFragment)

  }

  count_dir <- ValCountDir

  ######------------------------######

  ValCategory <- validateCharacter(category)

  if(!isTRUE(any(str_detect(gfffile$type, ValCategory)))){

    print(paste0('Category ', category, ' is not present in GFF file. Please check availability in the GFF'))

    stop(calculatePositionFragment)

  }

  #########-----------------------------------STARTING APP---------------------------------------###########

  ##selection gene rows in gff

  dirs_alignment          <- list.dirs(path = alignment_path, full.names = FALSE)

  genes                   <- filter(gfffile, type == category)

  print(paste0('Selected category: ', category))

  for(SAMPLE in dirs_alignment){

    if (!SAMPLE==""){

      print(SAMPLE)

      mapping_file        <- paste0(alignment_path, SAMPLE, '/alignment.sam')

      alignment           <- read.table(mapping_file, skip=22, header=FALSE, col.names =paste0('V',seq_len(21)), fill=TRUE)

      if (!file.exists(paste0(count_dir, SAMPLE))){

        dir.create(path=paste0(count_dir,'/',SAMPLE))

        count_res         <- paste0(count_dir,'/',SAMPLE)

        }else{

          unlink(paste0(count_dir, SAMPLE))

          dir.create(path=paste0(count_dir,'/',SAMPLE))

          count_res       <- paste0(count_dir,'/',SAMPLE)

        }

    ##Finding positions inside and outside the genes

      print('###---------------Calculating Fragments--------------###')
      for (positions in 1:nrow(alignment)){

       if (!exists('lengthFragment')){

         lengthFragment         <- nchar(alignment$V10[positions])

         nameFragment           <- alignment$V1[positions]

         chrpos                 <- alignment$V3[positions]

         positionFragmentInit   <- alignment$V4[positions]

         positionFragmentFinal  <- alignment$V4[positions] + lengthFragment

         lengthFragment         <- cbind (nameFragment,
                                          chrpos,
                                          lengthFragment,
                                          positionFragmentInit,
                                          positionFragmentFinal)

       }else{

         lengthFragment1          <- nchar(alignment$V10[positions])

         nameFragment1            <- alignment$V1[positions]

         chrpos1                  <- alignment$V3[positions]

         positionFragmentInit1    <- alignment$V4[positions]

         positionFragmentFinal1   <- alignment$V4[positions] + lengthFragment1

         lengthFragment1          <- cbind (nameFragment1, chrpos1, lengthFragment1, positionFragmentInit1, positionFragmentFinal1)

         lengthFragment           <- rbind(lengthFragment,lengthFragment1)
       }

      }

      lengthFragment_df                     <- data.frame(lengthFragment)

      identifiedGenes                       <- lengthFragment_df %>% filter(grepl('chr', chrpos))

      identifiedGenes                       <- data.frame(identifiedGenes)

      identifiedGenes$positionFragmentInit  <- as.numeric(identifiedGenes$positionFragmentInit)

      identifiedGenes$positionFragmentFinal <- as.numeric(identifiedGenes$positionFragmentFinal)

      print('###---------------Starting Calculate Position from calculated length fragments-------------###')

      for (genesselection in 1:nrow(identifiedGenes)){

        chrselected                         <- filter(genes, seqid == identifiedGenes$chrpos[genesselection])

        value_finding                       <- filter(genes, start == chrselected$start[which.min(abs(chrselected$start - identifiedGenes$positionFragmentInit[genesselection]))])

        if(value_finding$strand == '+'){

           if(identifiedGenes$positionFragmentInit[genesselection] < value_finding$start){

              promotor                      <- value_finding$start - 1000

              if(identifiedGenes$positionFragmentFinal[genesselection] < promotor){

                value <- 'Outside gene'

                }else{

                  value <- 'Promotor'

                  }

              }else if(identifiedGenes$positionFragmentFinal[genesselection] > value_finding$end){

                if (identifiedGenes$positionFragmentFinal[genesselection] > value_finding$end){

                  value <- 'Outside gene'

                  }else{

                    value <- 'Inside gene region'

                    }

                }else{

                  value <- 'Inside gene region'

                  }

          }else{

            if(identifiedGenes$positionFragmentInit[genesselection] > value_finding$end){

              promotor <- value_finding$start + 1000

              if(identifiedGenes$positionFragmentInit[genesselection] > promotor){

                value <- 'Outside gene'

                }else{

                  value <- 'Promotor'

                  }

              }else if(identifiedGenes$positionFragmentInit[genesselection] < value_finding$start){

                if (identifiedGenes$positionFragmentFinal[genesselection] < value_finding$start){

                  value <- 'Outside gene'

                  }else{

                    value <- 'Inside gene region'

                    }

                }else{

                  value <- 'Inside gene region'

                  }

            }

        if (!exists("value_finding1")){

           value_finding1                         <- value_finding

           value_finding1$Results                 <- value

           value_finding1$InitPosFragment         <- identifiedGenes$positionFragmentInit[genesselection]

           value_finding1$FinalPosFragment        <- identifiedGenes$positionFragmentFinal[genesselection]

           }else{

               value_finding2                     <- value_finding

               value_finding2$Results             <- value

               value_finding2$InitPosFragment     <- identifiedGenes$positionFragmentInit[genesselection]

               value_finding2$FinalPosFragment    <- identifiedGenes$positionFragmentFinal[genesselection]

               value_finding1                     <- rbind(value_finding1, value_finding2)

               }

        }

      frequencies_mapping                  <- table(value_finding1$Results)

      frequencies_mapping                  <- data.frame(frequencies_mapping)

      summatory_rows                       <- sum(frequencies_mapping$Freq)

      for (x in 1:nrow(frequencies_mapping)){

        frequencies_mapping$Percentages[x] <- (frequencies_mapping$Freq[x]*100)/summatory_rows

      }

      print('###------------------Writing results-----------------###')

      write.table(value_finding1, paste0 (count_res, "/resMapping.csv"), sep=';')

      write.table(frequencies_mapping, paste0(count_res,'/metrics.csv'), sep=';')

      rm(value_finding1,value_finding2, value_finding, lengthFragment, lengthFragment1)
      }

  }

}
