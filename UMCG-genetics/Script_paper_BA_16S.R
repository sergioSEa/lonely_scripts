#Analysis of 16S bacterial abundance in Mice gut.



setwd("~/Resilio Sync/Transfer/Collaborations/Folkert/")

library(tidyverse)
library(vegan)
library(phyloseq)
library(ape)
library(ggforce)
library(dunn.test)
library(patchwork)
library(nlme)
select = dplyr::select

#Open files
ASVs = read_csv(file = "Matrix_abundance.csv")
str_replace(str_replace(ASVs$ID,"DMP06766_L1_", ""), "_","") -> new_id
str_replace(str_replace(new_id,"DMP06767L1_", ""), "_","") -> new_id
ASVs %>% mutate(ID=new_id) -> ASVs

Change_naming(ASVs) -> ASVs
ASV_tree = read.tree("16S.tree")
Genus_tree = read.tree("16S_Genus.tree")
Family_tree = read.tree("16S_Family.tree")
Order_tree = read.tree("16S_Order.tree")
Phylum_tree = read.tree("16S_Phylum.tree")

Abundance_genus = read_csv(file="Merged_Genus_abundance.csv")
Abundance_family = read_csv(file="Merged_Family_abundance.csv")
Abundance_order = read_csv(file="Merged_Order_abundance.csv")
Abundance_phylum = read_csv(file="Merged_Phylum_abundance.csv")
#Genus_columns = colnames(Abundance_table )[grepl(pattern = "Genus",colnames(Abundance_table ))]
#Genus_table = select(Abundance_table, c(ID,Genus_columns))



metadata = read_tsv("Metadata.txt")

sapply(metadata$Sample, FUN = function(x){ str_split(x, "\\.")[[1]][1] } ) %>% as.vector() ->Mice_names
metadata %>% mutate(Diet = ifelse(grepl("AL", Group), "AL", "CR"), Time = ifelse(grepl(".1", Group), "0","4"), Mouse = Mice_names) -> metadata
batch_mouse = read_tsv("batch_mouse.txt", col_names =F)
metadata %>% mutate(litter = as.factor(batch_mouse$X2)) -> metadata
metadata %>% mutate(Diet2 = ifelse(Time ==  0, "Control", Diet), Group2 = ifelse(grepl("1", Group), "Control", Group)) -> metadata

Preliminary_checks = function(Abundance, metadata){
  #Number of reads per sample
  apply(X = select(Abundance,-ID), MARGIN = 1, FUN = sum) -> Total_number
  Abundance %>% mutate(Total_number = Total_number) -> Abundance_w_number
  ggplot(Abundance_w_number) + geom_bar(aes(x=reorder(ID, -Total_number), y=Total_number), fill="#0080FF",col="black", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1)) -> Total_reads_plot
  #Number of 0s per sample
  apply(X = select(Abundance,-ID), MARGIN = 1, FUN = function(x){ 100*(sum(x==0)/sum(x)) } ) -> Percentage_0
  Abundance_w_number %>% mutate(Percentage_0 = Percentage_0) -> Abundance_w_number
  ggplot(Abundance_w_number) + geom_bar(aes(x=reorder(ID, -Total_number), y=Percentage_0), fill="#00FFFF", col="black", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(size = 7,angle = 90, hjust = 1)) -> Amount_0_plot
  #Both visualizations
  print(Total_reads_plot + Amount_0_plot) -> PLOT
  print(PLOT)
  #Number of species per sample
  apply(select(arrange(Abundance, ID), -ID), MARGIN = 1, FUN= function(x){ length(x[x>0]) }) -> N_ASVs
  metadata %>% filter(Sample %in% Abundance$ID) %>% arrange(Sample) -> metadata
  metadata %>% mutate(N_ASVs) %>% ggplot(aes(x=Group, y=N_ASVs, col=Diet)) + geom_boxplot()  + geom_point() + geom_line(aes(group=Mouse),linetype = 2) +theme_bw()  -> Fig_Abundance                 
  dunn.test::dunn.test(N_ASVs, metadata$Group) 
  print(Fig_Abundance)
  return(Total_number)
}
Filter_abundance = function(M, threshold=0.2, min_abundance = 0){
  #Filter species not seen in at least treshold*100% of samples
  apply(select(M, -ID), MARGIN=2, FUN = function(x){ (sum(x != 0)/length(x)) >= 0.2 }) -> Filter_vector
  F_names = colnames(select(M, -ID))[!Filter_vector]
  M %>% select(!F_names) -> M
  #Filter species with mean relative abundance lower than min_abundance
  apply(select(M, -ID), MARGIN=1, FUN = function(x){ sum(x) }) -> Total_reads_sample
  as_tibble(select(M, -ID)/Total_reads_sample) %>% summarise_all(mean) -> All_means
  colnames(select(M, -ID))[as_vector(All_means) > min_abundance] -> Taxa_not_to_filter
  M %>% select(c(ID, Taxa_not_to_filter)) -> M
  return(M)
}
Alpha_diversity = function(M, metadata){
  #Shanon diversity calculation per sample, plots and non-parametric test
  diversity = vegan::diversity(select(M,-ID), index = "invsimpson")
  metadata %>% mutate(Diversity = diversity) -> metadata
  
  metadata %>% ggplot(aes(x=Group, y=Diversity,col=Time)) + geom_boxplot(outlier.shape = NA) +geom_point() +
    theme_bw() + geom_line(aes(group=Mouse),col="black",linetype=2) -> Fig_alpha
  print(Fig_alpha)
  
  
  summary(lme(Diversity ~ litter + Diet*Time, random=~ 1|Mouse, metadata )) -> Model 
  print(Model)
  
  #dunn.test(x=metadata$Diversity, g= metadata$Group)
  
  #metadata %>% filter(Group == "AL1") -> AL_1
  #metadata %>% filter(Group == "AL2") -> AL_2
  #metadata %>% filter(Group == "CR1") -> CR_1
  #metadata %>% filter(Group == "CR2") -> CR_2
  #delta_AL = AL_1$Diversity - AL_2$Diversity
  #delta_CR = CR_1$Diversity - CR_2$Diversity
  #wilcox.test(delta_AL, delta_CR) -> Delta
  #print(Delta)
}
Beta_diversity = function(Phylobject, Distance_method="bray", SUB=NA){
  meta <- as(sample_data(Phylobject), "data.frame")
  
  if (Distance_method == "Weighted"){
    Ordination = ordinate(Phylobject, method = "PCoA", distance= "unifrac", weighted=T)
    adonis2(distance(Phylobject, method="unifrac", weighted=T) ~ litter + Time * Diet, data = meta) -> adonis_results
    
  } else if (Distance_method == "Unifrac"){
    Ordination = ordinate(Phylobject, method = "PCoA", distance= "unifrac", weighted=F)
    adonis2(distance(Phylobject, method="unifrac", weighted=F) ~ litter + Time * Diet, data = meta) -> adonis_results
    
    
  } else if (Distance_method == "bray") {
    Ordination = ordinate(Phylobject, method = "PCoA", distance= "bray")
    #adonis2(distance(Phylobject, method="bray") ~ litter + Time *  Diet, data = meta) -> adonis_results
    #With repeated
    adonis2(distance(Phylobject, method="bray", weighted=F) ~ litter + Time * Diet + Mouse, data = meta, strata=Mouse) -> adonis_results
  }
  
  plot_ordination(Phylobject, Ordination, color="Time", shape="Diet") + theme_bw() + geom_line(aes(group=Mouse), col="grey", linetype=2, size=0.5) -> PCoA_plot
  plot_ordination(Phylobject, Ordination, color="litter", shape="Diet") + theme_bw() -> PCoA_plot2
  
  print(PCoA_plot)
  print(PCoA_plot2)
  print(adonis_results)
  
}
Prepare_phylobject = function(M, metadata, Tree=NA){
  
  #Prepare abundance object
  M %>% select(-ID) %>% t() %>% as_tibble() %>% mutate(label=colnames(select(M,-ID))) %>% `colnames<-`(c(M$ID, "label")) %>% select(c("label",M$ID)) -> trans_M
  trans_M %>% as.data.frame() %>% column_to_rownames("label") -> ASV_table
  otu_table(ASV_table, taxa_are_rows = TRUE) -> ASV_table
  #Prepare sample name object
  taxmat = matrix(rownames(ASV_table))
  rownames(taxmat) = rownames(ASV_table)
  SampleData = metadata %>% as.data.frame() %>% column_to_rownames("Sample")
  SampleData = sample_data(SampleData)
  
  tax_table(taxmat) -> taxmat
  phyloseq::phyloseq(ASV_table, taxmat, SampleData) -> phylobjetc
  
  if (! is.na(Tree)){ phyloseq::merge_phyloseq(phylobjetc,Tree) -> phylobjetc }
  
  
  return(phylobjetc)
  
}
Rarefy_reads = function(M, metadata){
  Number_reads = apply(select(M, -ID), MARGIN= 1, FUN = function(x){ sum(x) })
  
  ##Check difference in number of species before and after rarefy
  S <- specnumber(select(M, -ID)) # observed number of species
  print(paste("Rarefying to ", min(Number_reads), collapse = " "))
  Srare <- rarefy(x=select(M,-ID),MARGIN=1, sample=min(Number_reads))
  plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
  #rarecurve(select(M,-ID), sample=50317)
  
  ##Perform rarefy
  vegan::rrarefy(x=select(M,-ID), sample=min(Number_reads)) %>% as_tibble() %>% mutate(ID = M$ID) -> Abundance_rarefy
    
  return(Abundance_rarefy)
  
}
Heatmap_composition = function(M, metadata, Top_lineages=NA, ASV = "no"){
  if (ASV == "yes"){
    sapply(str_split(colnames(select(M, -ID)), "\\|"), FUN = function(x){ paste(x[(length(x)-1):length(x)], collapse="|")}) -> N
    colnames(M) = c(N,"ID")
  }
  #Distance among samples
  dist <- vegdist(select(M,-ID),  method = "bray")
  
  #Distance among bacteria
  if (! is.na(Top_lineages)){
    select(M, -ID) %>% summarise_all(sum) %>% as_vector() -> Taxa_total_numbers
    names(sort(desc(Taxa_total_numbers))[1:Top_lineages]) -> Top_bacteria
  } else { Top_bacteria = colnames(select(M, -ID)) }
  
  dist_Bacteria <- vegdist(t(select(M,Top_bacteria)),  method = "bray")
  #Matrix values
  matrix_composition = as.matrix(select(M,-ID))
  rownames(matrix_composition) <- M$ID ; colnames(matrix_composition) <- colnames(select(M,-ID))
  
  Annotation_df = as.data.frame(select(metadata, c(Time, Diet, litter)))
  rownames(Annotation_df) = M$ID
  
  pheatmap(matrix_composition,
           clustering_distance_rows= dist,
           clustering_distance_cols=dist_Bacteria,
           annotation_row = Annotation_df)
  
  
}
Change_naming = function(M){
  sapply(M$ID, FUN = function(x){paste(strsplit(x,"_")[[1]][3],collapse="_")})  -> New_names
  M %>% mutate(ID=as.vector(New_names)) -> M
  return(M)
}
Transformation_composition = function(Count_table){
  #Clr transformation
  library(microbiome)
  SampleID = Count_table$ID
  Count_table %>% select(-ID) -> Counts
  #Transform counts
  colnames(Counts) -> N
  Counts %>% t() -> Check
  Counts_transformed = as_tibble(t(abundances(x=as.data.frame(Check), transform="clr")))
  ###
  #Counts_transformed %>% mutate(ID = SampleID) %>% select(ID,colnames(Counts_transformed)[1:(length(colnames(Counts_transformed))-1)]) -> Counts_transformed
  Counts_transformed %>% mutate(ID = SampleID) %>% select(ID,colnames(Counts_transformed)[1:(length(colnames(Counts_transformed)))]) -> Counts_transformed
  return(Counts_transformed)
}
Test_groups = function(DF, TEST="All"){
  #Statistical test
  if (TEST == "All"){
    Model = dunn.test(x=DF$Dependent, g= DF$Group)
    return(Model$P)
  }else if (TEST == "Paired"){
    DF %>% filter(Group %in% c("AL1", "AL2")) -> DF_1
    wilcox.test(arrange(filter(DF_1, Group=="AL1"), Mouse)$Dependent, arrange(filter(DF_1, Group=="AL2"), Mouse)$Dependent, paired = T) -> Model
    DF %>% filter(Group %in% c("CR1", "CR2")) -> DF_2
    wilcox.test(arrange(filter(DF_2, Group=="CR1"),Mouse)$Dependent, arrange(filter(DF_2, Group=="CR2"), Mouse)$Dependent, paired = T) -> Model2
    
    return(c(Model$p.value, Model2$p.value))
  }
  else if (TEST == "delta"){
    summary(lm(Dependent ~ Diet + Litter, DF)) -> Model
    return(Model$coefficients[14])
  } else if (TEST == "Mixed"){
    library(nlme)
    DF %>% mutate(Group2 = factor(Group2,levels=c("Control", "AL2", "CR2"))) -> DF
    summary(lme(Dependent ~ Group2 + litter, random=~ 1|Mouse, DF )) -> Model 
    DF %>% filter(! Group2 == "Control") -> DF2
    summary(lme(Dependent ~ Group2 + litter, random=~ 1|Mouse, DF2 )) -> Model2 
    Result =  c(Model$tTable[2], Model$tTable[3], Model$tTable[30], Model$tTable[31], Model2$tTable[2], Model2$tTable[26])
    return(Result)
  } 
  
}
Delta_analysis = function(R_abundance, metadata){
  Filter_abundance(R_abundance) -> R_Abundance_f
  
  norm_abundance = Transformation_composition(R_abundance)

  norm_abundance %>% mutate(Group = metadata$Group, Mouse = metadata$Mouse, Litter=metadata$litter) -> Meta_Abundance
  Meta_Abundance %>% filter(Group == "AL1") %>% arrange(Mouse) %>% select(-c(Group,Mouse,Litter,ID)) -> Group_1 
  Meta_Abundance %>% filter(Group == "AL2") %>% arrange(Mouse) %>% select(-c(Group,Mouse,Litter,ID)) -> Group_2 
  Meta_Abundance %>% filter(Group == "CR1") %>% arrange(Mouse) %>% select(-c(Group,Mouse,Litter,ID)) -> Group_1_2 
  Meta_Abundance %>% filter(Group == "CR2") %>% arrange(Mouse) %>% select(-c(Group,Mouse,Litter,ID)) -> Group_2_2 
  Delta_AL = as_tibble(Group_1 - Group_2) %>% mutate(Diet = "AL")
  Delta_CR = as_tibble(Group_1_2 - Group_2_2) %>% mutate(Diet = "CR")
  Delta_abundance = rbind(Delta_AL, Delta_CR)
  Meta_Abundance2 = Meta_Abundance[1:(dim(Meta_Abundance)[1]/2),]
  
  apply(select(Delta_abundance,-Diet), MARGIN = 2, FUN = function(x){tibble(Dependent=x,Diet=Delta_abundance$Diet) %>% group_by(Diet) %>% summarise(N = mean(Dependent)) -> Means ; Means$N[1] - Means$N[2] }) -> Mean_differences
  apply(select(Delta_abundance,-Diet), MARGIN = 2, FUN = function(x){ print(x);Test_groups(tibble(Dependent=x,Diet=Delta_abundance$Diet, Litter=Meta_Abundance2$Litter), TEST="delta")}) -> Model_outcome_rarefy
  Model_outcome_rarefy %>%as_tibble() %>% mutate(Bug=names(Model_outcome_rarefy)) %>%
    mutate(Difference_means = as.vector(Mean_differences)) %>% filter(Bug %in% colnames(R_Abundance_f)) %>%
    mutate(FDR = p.adjust(value, "fdr")) %>% filter(!is.na(FDR)) %>% arrange(FDR) -> Model_outcome_rarefy
  Model_outcome_rarefy %>% ggplot(aes(x=Difference_means, y=-log10(value), col=FDR<0.05)) + geom_point() + theme_bw() -> PLOT
  
  print(PLOT)
  return(Model_outcome_rarefy)
  
}
Check_lineages = function(Stats, Abundance, metadata ,Taxa){
  Stats %>% filter(Bug  %in% Taxa) -> Stats
  for (Bacteria in Taxa){
    Stats %>% filter(Bug == Bacteria) -> Bacteria_stats
    Abundance %>% select(Bacteria) %>% as_vector %>% as.vector() -> Ab
    metadata %>% mutate(Bacterial_abundance = Ab) %>% ggplot(aes(x=Group2, y=Bacterial_abundance)) +
      geom_boxplot() +geom_point() + theme_bw() + geom_line(aes(group=Mouse), col="grey", linetype=2, size=0.5) +
      ggtitle(Bacteria) -> Figure

    print(Bacteria)
    print(Bacteria_stats)
    print(Figure)
  }
  
}
Calculate_fdr = function(norm_abundance, metadata, Real_pvalues){
  FDR_vector = vector()
  #1st, Number of permutations
  Random_distribution  = tibble()
  N_perm = 100
  for (i in seq(N_perm)){
    #Per each permutation rearrange data and fit model
    #rearrange
    sample_n(norm_abundance, size= nrow(norm_abundance)) -> Rearrange
    #Fitmodel
    apply(select(Rearrange, -ID), MARGIN = 2, FUN = function(x){ Test_groups(mutate(metadata, Dependent=x),TEST="Mixed")}) -> H0
    #Save Pvalues in the random distribution of pvalues
    H0 %>% t() %>%as_tibble() %>% mutate(Bug=colnames(Model_outcome_rarefy)) %>%   
      `colnames<-`(c("AL_est", "CR_est", "AL_P", "CR_P", "Diet_est", "Diet_P", "Bug")) %>% gather(Test, Pvalue, c("AL_P", "CR_P","Diet_P")) -> H0
    random_tibble = tibble(AL_P = filter(H0, Test == "AL_P")$Pvalue, CR_P = filter(H0, Test == "CR_P")$Pvalue, Diet_P = filter(H0, Test == "Diet_P")$Pvalue)    
  
    Random_distribution = rbind(Random_distribution, random_tibble)
  }
  #2nd, calculate FDR per each threshold
  H1  = function(Pheno, Distribution, Real_distribution){
    Distribution %>% select(Pheno) %>% as_vector() %>% as.vector() %>% as.numeric() -> Random_distribution
    Real_distribution %>% filter(Test== Pheno) %>% select(Pvalue) %>% arrange(desc(Pvalue)) %>% as_vector() %>% as.vector() %>%as.numeric() -> Real_distribution2
    FDR_vector = c()
    for (Threshold in Real_distribution2){
      Threshold = as.numeric(Threshold)
      P = sum((Real_pvalues <= Threshold) *1)
      FP = sum((Random_distribution <= Threshold)*1)
      FDR = (FP/N_perm)/P
      if (FDR >1){ FDR = 1}
      #If the FDR value of a higher Pvalue is lower, make the FDR of the lower Pvalue at least as low
      if (length(FDR_vector) > 0){ if (FDR > FDR_vector[length(FDR_vector)]){ FDR = FDR_vector[length(FDR_vector)]} }
      FDR_vector = c(FDR_vector,FDR)
    }
    Real_distribution %>% filter(Test== Pheno)  %>% arrange(desc(Pvalue)) %>% mutate(FDR = FDR_vector) -> FDR_vector
    return(FDR_vector)
  }
  lapply(c("AL_P", "CR_P", "Diet_P"), FUN = function(x){ H1(x, Random_distribution,Real_pvalues)  }) -> tibbles
  rbind(tibbles[[1]], tibbles[[2]], tibbles[[3]]) -> New_Real 
  return(New_Real)
}



To_analyze = list(ASVs, Abundance_genus, Abundance_family, Abundance_order, Abundance_phylum)
Trees_list = list(ASV_tree, Genus_tree, Family_tree, Order_tree, Phylum_tree)
Names_list = list("ASV","Genus", "Family", "Order", "Phylum")

set.seed(1194)
library(pheatmap)

#Initial analysis in ASV and Phylum levels. Outcomes: no clear differences in diversity. Beta diversity shows an effect of litter, and an effect of Time:Diet if abundance is used
for (Number_entry in c(1,5)){
  Abundance = To_analyze[[Number_entry]] ; Tree = Trees_list[[Number_entry]] ; Abundance_level = Names_list[[Number_entry]]
  Abundance %>% arrange(ID) -> Abundance ; metadata %>% arrange(Sample) -> metadata
  Preliminary_checks(Abundance, metadata)
  colnames(Abundance)[grepl("Unknown", colnames(Abundance))] -> Remove_counts ;Abundance %>% select(-Remove_counts) -> Abundance
  R_Abundance = Rarefy_reads(Abundance) 
  N_taxa = apply(select(R_Abundance, -ID), MARGIN=1 , FUN = function(x){ length(x[! x==0]) }) ;metadata %>% mutate(Number_taxa = N_taxa) -> metadata_n
  summary(lme(Number_taxa ~ litter + Diet*Time, random=~ 1|Mouse, metadata_n ))
  ggplot(metadata_n,aes(x=Group, y=N_taxa)) + geom_point()+geom_line(aes(group=Mouse)) + theme_bw()
  Heatmap_composition(R_Abundance, metadata, Top_lineages=20, ASV="yes")
  Alpha_diversity(R_Abundance,metadata)
  Prepare_phylobject(R_Abundance,metadata,Tree=Tree) -> phylo_data
  Beta_diversity(phylo_data)
}


for (Number_entry in c(2,5)){
  Abundance = To_analyze[[Number_entry]]
  Tree = Trees_list[[Number_entry]]
  Abundance_level = Names_list[[Number_entry]]
  
  Abundance %>% arrange(ID) -> Abundance
  metadata %>% arrange(Sample) -> metadata
  
  
  
  #Preliminary_checks(Abundance, metadata)
  colnames(Abundance)[grepl("Unknown", colnames(Abundance))] -> Remove_counts
  Abundance %>% select(-Remove_counts) -> Abundance
  
  R_Abundance = Rarefy_reads(Abundance)

    #Target lineages
  colnames(R_Abundance)[grepl("Lactobacillus|Clostridium|Bifidobacterium|Bacteroides|Ruminococcus gnavus",colnames(R_Abundance))] -> BA_bacteria
  #Targeted input
  select(R_Abundance, c("ID",BA_bacteria)) -> BA_bacteria
  #CLR-transformed tests
  Comparisons = c("AL1-AL2","AL1-CR1","AL2-CR1","AL1-CR2","AL2-CR2","CR1-CR2")
  Filter_abundance(R_Abundance) -> R_Abundance_f
  norm_abundance = Transformation_composition(R_Abundance)
  #norm_abundance = Transformation_composition(BA_bacteria)
  
  
  norm_abundance %>% select(colnames(R_Abundance_f)[colnames(R_Abundance_f) %in% colnames(norm_abundance)]) -> norm_abundance
  
  #2. Mixed Model
  apply(select(norm_abundance, -ID), MARGIN = 2, FUN = function(x){ print(x);Test_groups(mutate(metadata, Dependent=x),TEST="Mixed")}) -> Model_outcome_rarefy
  
  ##Calculate FDR by permutatations on each Test independently?
  Model_outcome_rarefy %>% t() %>%as_tibble() %>% mutate(Bug=colnames(Model_outcome_rarefy)) %>%   
    `colnames<-`(c("AL_est", "CR_est", "AL_P", "CR_P", "Diet_est", "Diet_P", "Bug")) %>% gather(Test, Pvalue, c("AL_P", "CR_P","Diet_P")) %>%
    mutate(FDR = p.adjust(Pvalue, "fdr")) %>% filter(!is.na(FDR)) %>% arrange(FDR) -> Model_outcome_rarefy_mixed
  
  Calculate_fdr(norm_abundance, metadata, Model_outcome_rarefy_mixed) -> Model_outcome_rarefy_mixed2
   
   Model_outcome_rarefy_mixed %>% ggplot(aes(x=Test, y=-log10(Pvalue))) + geom_boxplot(outlier.shape = NA
   )+ geom_jitter(aes(col=FDR<0.05)) +theme_bw()+ coord_flip()

   Check_lineages(Stats=Model_outcome_rarefy_mixed2,Abundance=norm_abundance,metadata=metadata,Taxa =  BA_bacteria)
}




