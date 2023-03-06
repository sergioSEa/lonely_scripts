Simulate_sequences = function(Sequence_size=10, Mutation_Probability=0.1, Iteration_mutagenesis=10, Number_samples=4, N_strain_prob=c(0.75, 0.2,0.05), Major_strain_frequency_range=c(0.98, 0.8), Error_measurement_frequency=c(0.1, 0 ) ){
  #Simulate sequences

  #Generate strain collection
  Ancestral_sequence = sample(c("A", "T", "G", "C"), size = Sequence_size, replace=T) 
  Strains = list( Ancestral_sequence )
  postion = 2 #index to save the new strain in the list if Strains
  for (i in seq(Iteration_mutagenesis)){
  for (Strain in Strains){
    print(Strain)
    New_seq = Strain
    R = runif(Sequence_size)
    Mutations = R  < Mutation_Probability
    if (sum(Mutations) < 1){ next }
    for (i in seq(Mutations)){
      if( Mutations[i] == F ){ next }
      New = sample(c("A", "T", "G", "C"), size = 1, replace=T)
      New_seq[i] = New
    }
    Strains[[postion]] = New_seq
    postion = postion + 1
  }
}


  #assign strains per person
  list_strains = list(Number_samples)

  for (i in seq(Number_samples)){
    #How many strains will this person have?
    N_strain  = sample(c(1,2,3), prob = N_strain_prob, size = 1)
    #Select the genome of that strain
    sample(Strains, N_strain) -> Strains_person
    #Save strain
    list_strains[[i]] = Strains_person
    #Add frequency of strain
    if (N_strain == 1){
      Major = rep(1, Sequence_size)
      Minor = rep(0, Sequence_size)
    } else {
      runif(1, max=Major_strain_frequency_range[1], min=Major_strain_frequency_range[2]) -> Major
      variation = runif(Sequence_size, max=Error_measurement_frequency[1], min=Error_measurement_frequency[2])
      variation + Major -> Major
      Major[Major > 1] = 1 
      Minor = 1 - Major 
      if (N_strain == 3){ Minor = Minor/2 }
    }
    for (S in seq(list_strains[[i]]) ) {
      if (S == 1){
        names(list_strains[[i]][[S]] ) = Major
      } else {
        names(list_strains[[i]][[S]] ) = Minor
      }
    }
}

  #Calculate strain frequencies
  Strain_df = tibble() #Stores strain sequences
  Frequency_df = tibble() #stores nucleotides frequencies per strain
  Sample_n = 1
  for ( Sample in list_strains ){
    strain_item = 1
    for (Strain in Sample){
      tibble( ID =  Sample_n, Sequence= paste0(Strain, collapse=""), Strain_id= paste0(ID, "_", strain_item) ) %>% rbind(Strain_df, .) -> Strain_df
      strain_item = strain_item + 1
    }
    for (Position in  seq(Strain) ){
      for (Strain in Sample){
        Nucleotide = Strain[Position]
        Nucleotide_f = as.numeric(names(Strain[Position]))
        rbind(Frequency_df, tibble(Position = Position, Nucleotide = Nucleotide, Freq=Nucleotide_f ,ID= Sample_n, Strain=paste0(Strain, collapse="") ) ) ->  Frequency_df
      }
      
    }
    Sample_n = Sample_n + 1
    
  }
  #Group nucleotides from different strains that are identical to obtain the observed frequencies
  Frequency_df %>% group_by(Position, Nucleotide, ID) %>% summarise(Freq = sum(Freq)) %>% ungroup() -> Observed_frequencies
  #Make in wide forma
  Observed_frequencies %>% mutate(SNP = paste0(Position, "_", Nucleotide) ) %>% select(! c(Nucleotide, Position)) %>% spread( SNP, Freq) -> Wide_observed
  Wide_observed[is.na(Wide_observed)] = 0
  #We have several columns per position. Get only one, and use the majority variant as the reference (1) and the minor variant(s) as alternative
  Majority = c()
  for (N in seq(Sequence_size)){
   colnames(Wide_observed)[grepl( paste0(N,"_") , colnames(Wide_observed))] -> Check
    Wide_observed %>% select(Check) -> Check
    if (dim(Check)[2] == 1){ Majority = c(Majority, colnames(Check)) 
    } else{
      Check %>% apply(2, sum) -> C
      Majority = c(Majority,names(which.max(C)) )
    }
  }
  Wide_observed %>% select(c("ID", Majority)) -> Wide_observed_ref
  
  #Calculate distances and do clustering
  
  #Differences using frequencies
  Wide_observed_ref %>% as.data.frame() %>% column_to_rownames("ID") %>% vegan::vegdist("euclidean") -> Dis_1
  hclust(Dis_1)  %>% plot(main="Distance between frequencies")
  #differences using genotypes (0/1/2)
  Wide_observed_ref %>%  as.data.frame() %>% column_to_rownames("ID") %>% apply(2, function(x){ sapply(x, function(y){ ifelse(y==1, 2, ifelse(y==0, 0, 1 ) ) } )  } ) %>% as_tibble() %>% mutate(ID =Wide_observed_ref$ID ) %>%
    as.data.frame() %>% column_to_rownames("ID") %>% vegan::vegdist("euclidean") -> Dis_2
  hclust(Dis_2)  %>% plot(main="Distance between genotypes")
  #differences based in distances between strains
  Dist_matrix_strains = matrix(0, length(Strain_df$Sequence), length(Strain_df$Sequence))
  #Compute hamming distances between strains
  for (S1 in seq(length(Strain_df$Sequence))){
    for (S2 in seq(length(Strain_df$Sequence))){
      sum(  str_split(Strain_df$Sequence[S1], "")[[1]] !=  str_split(Strain_df$Sequence[S2], "")[[1]] ) -> D
      Dist_matrix_strains[S1,S2] = D
    }
  }
  rownames(Dist_matrix_strains) = Strain_df$Strain_id
  hclust(as.dist(Dist_matrix_strains) ) %>% plot(main="Distance between strains")
  
  
  return(list(list_strains, Frequency_df, Observed_frequencies, Wide_observed_ref ))
}
