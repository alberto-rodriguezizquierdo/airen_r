#' @name randomGenomeFragmentation
#' @param a
#' @import seqinr, phylotools, datatable, tydiverse
#' @return fragmentRes
#' @author Alberto Rodriguez-Izquierdo

randomGenomeFragmentation <- function(genome_path, max_size, min_size, nb_repeats){
  
  
  a <- phylotools::read.fasta(file=genome_path)
  
  #, seqtype = 'DNA', as.string=TRUE
  
  # 70000 sequences
  # interval (200,260)
  #Calculate length contig
  length_genome <- nrow(a)
  
  a <- data.frame(a)
  print('------------------Counting bp-----------------')
  for (contigs in 1:nrow(a)){
    
    if (!exists('contig_count')){
      
      contig_count <- nchar(a$seq.text[1])
      
      contig_count <- data.frame(cbind(a$seq.name[1], contig_count))
      
    }else{
      
      contig_count1 <- nchar(a$seq.text[contigs])
      
      contig_count1 <- data.frame(cbind(a$seq.name[contigs], contig_count1))
      
      names(contig_count1) <- names(contig_count)
      
      contig_count <- rbind(contig_count, contig_count1)
      
    }
    
  }
  
  contig_count$contig_count <- as.numeric(contig_count$contig_count)
  
  total_weight_count <- sum(contig_count$contig_count)
  print('------------Calculating weights------------')
  for (weights in 1:nrow(contig_count)){
    
    if (!exists('contigWeights')){
      
      contigWeights <- contig_count$contig_count[weights] / total_weight_count
      
      contigWeights <- data.frame(cbind(contig_count$V1[weights], contigWeights))
      
    }else{
      
      contigWeights_1 <- contig_count$contig_count[weights] / total_weight_count
      
      contigWeights_1 <- data.frame(cbind(contig_count$V1[weights], contigWeights_1))
      
      names(contigWeights_1) <- names(contigWeights)
      
      contigWeights <- rbind(contigWeights, contigWeights_1)
      
    }
    
  }
  rm(contigWeights_1)
  
  contigWeights$contigWeights <- as.numeric(contigWeights$contigWeights)
  
  sum(contigWeights$contigWeights)
  ## Calculate number fragments depending on their weights
  
  fragmentNumber <- nb_repeats
  
  print('------------Weighting fragments------------')
  #rm(weightFragmentNumber1)
  for (eachContig in 1:nrow(contigWeights)){
    
    if(!exists('weightFragmentNumber')){
      
      weightFragmentNumber <- round(contigWeights$contigWeights[eachContig]*fragmentNumber)
      
      weightFragmentNumber <- data.frame(cbind(contigWeights$V1[eachContig], weightFragmentNumber))
      
    }else{
      
      weightFragmentNumber1 <- round(contigWeights$contigWeights[eachContig]*fragmentNumber)
      
      weightFragmentNumber1 <- data.frame(cbind(contigWeights$V1[eachContig], weightFragmentNumber1))
      
      names(weightFragmentNumber1) <- names(weightFragmentNumber)
      
      weightFragmentNumber <- rbind(weightFragmentNumber, weightFragmentNumber1)
      
    }
  
  }
  rm(weightFragmentNumber1)
  
  weightFragmentNumber$weightFragmentNumber <- as.numeric(weightFragmentNumber$weightFragmentNumber)
  
  sum(weightFragmentNumber$weightFragmentNumber)
  
  #loop for calculate contig restriction and position fragment
  
  print('-------Generating position fragments------')
  
  #rm (length_fragment1)
  for (contigsFinal in 1:nrow(weightFragmentNumber)){
    
    numberRestrictionFragments <- weightFragmentNumber$weightFragmentNumber[contigsFinal]
    
    for(numberRest in 1:numberRestrictionFragments){
      
      if(!exists('length_fragment')){
        
        length_fragment <- round(runif(1,min_size,max_size))
        
        position_fragment_min <- round(runif(1,min=0, max=contig_count$contig_count[contigsFinal]))
        
        max_allow_position <- contig_count$contig_count[contigsFinal]-length_fragment
        
        while (position_fragment_min > max_allow_position){
          
          position_fragment_min <- round(runif(1,min=0, max=contig_count$contig_count[contigsFinal]))
          
        }
        
        position_fragment_max <- position_fragment_min + length_fragment
        
        length_fragment <- data.frame(cbind(weightFragmentNumber$V1[contigsFinal],length_fragment, position_fragment_min, position_fragment_max))
        
        
      }else{
        
        length_fragment1 <- round(runif(1,min_size,max_size))
        
        position_fragment_min <- round(runif(1,min=0, max=contig_count$contig_count[contigsFinal]))
        
        max_allow_position <- contig_count$contig_count[contigsFinal]-length_fragment1
        
        while (position_fragment_min > max_allow_position){
          
          position_fragment_min <- round(runif(1,min=0, max=(contig_count$contig_count[contigsFinal])))
          
        }
        
        position_fragment_max <- position_fragment_min + length_fragment1
        
        length_fragment1 <- data.frame(cbind(weightFragmentNumber$V1[contigsFinal], length_fragment1, position_fragment_min, position_fragment_max))
        
        names(length_fragment1) <- names(length_fragment)
        
        length_fragment <- rbind(length_fragment, length_fragment1)
        
        
      }
    }
  }
  
  rm (length_fragment1)
  
  length_fragment$length_fragment <- as.numeric(length_fragment$length_fragment)
  
  length_fragment$position_fragment_min <- as.numeric(length_fragment$position_fragment_min)
  
  length_fragment$position_fragment_max <- as.numeric(length_fragment$position_fragment_max)
  
  
  print('-------Generating fragments...------')
  
  for (fragment in 1:nrow(length_fragment)){
    
    if(!exists('fragmentRes')){
      
      print(fragment)
      
      selectFragment <- a[a$seq.name %like% length_fragment$V1[fragment],]
      
      fragmentRes <- substr(selectFragment$seq.text, length_fragment$position_fragment_min[fragment], length_fragment$position_fragment_max[fragment])
      
      fragmentRes <- data.frame(cbind(length_fragment$V1[fragment],length_fragment$position_fragment_min[fragment],length_fragment$position_fragment_max[fragment], fragmentRes))
      
    }else{
      
      print(fragment)
      
      selectFragment <- a[a$seq.name %like% length_fragment$V1[fragment],]
      
      fragmentRes1 <- substr(selectFragment$seq.text, length_fragment$position_fragment_min[fragment], length_fragment$position_fragment_max[fragment])
      
      fragmentRes1 <- data.frame(cbind(length_fragment$V1[fragment],length_fragment$position_fragment_min[fragment],length_fragment$position_fragment_max[fragment], fragmentRes1))    
      
      names(fragmentRes1) <- names(fragmentRes)
      
      fragmentRes <- rbind(fragmentRes, fragmentRes1)
      
    }
  }
  rm(fragmentRes1)
  
  return(fragmentRes)
}




