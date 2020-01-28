setwd("~/PhD/WORK/Nijmegen/")

library(tidyverse)
library(vegan)
library(ape)
library(ggforce)
library(ggtree)
library(cluster)

##Load matrix, remove uninfromative row and get records that corrrespond to a given taxonomic level
matrix_taxonomy = read_tsv("Merged_batch.tsv")
matrix_taxonomy[-1,] -> matrix_taxonomy
matrix_taxonomy %>% filter(grepl("s__",ID)) %>% filter(!grepl("t__",ID))  -> species_matrix
matrix_taxonomy %>% filter(grepl("g__",ID) & !grepl("s__",ID)) -> genus_matrix
matrix_taxonomy %>% filter(grepl("f__",ID) & !grepl("g__",ID)) -> family_matrix

##Load metadata of the data
sample_metadata = read_delim("Metadata.csv", delim =";")


sample_metadata %>% mutate(ID_patient = StudieID ,StudieID = paste0(sample_metadata$StudieID,"P"), Source="pocket")-> metadata_P  
sample_metadata %>% mutate(ID_patient = StudieID, StudieID = paste0(sample_metadata$StudieID,"S"), Source="saliva")-> metadata_S  
sample_metadata = rbind(metadata_P,metadata_S)



###Functions

###Metrics Species (raw data with everything before repving anything)
do_Metrics = function(MData,N){
  ID = MData[dim(MData)[2]] 
  MData %>% dplyr::select(-one_of(c("ID"))) -> MData
  number_0 = function(V){
    N_0 = 0
    for (i in V){
      if(i == 0){ N_0 = N_0 + 1}
    }
    return(N_0)
  }
  
  Total_summary = tibble()
  for (C in seq(1:dim(MData)[2])){
  MData[,C] %>% mutate_if(is.character, as.numeric) %>% summarise_all(list(min, max,median,mean,sd,number_0)) %>% mutate(Name = colnames(MData[,C])) -> SUMMARY
  colnames(SUMMARY) = c("Minimum_abundance", "Maximum_abundance", "Median_abundance", "Mean_abundance","Standard_deviation","Number_of_0_counts","Organism")
  rbind(Total_summary, SUMMARY) -> Total_summary
  }
  write_csv(x = Total_summary, path = N)
}

contamination_stats = function(File = "Reads_contaminated.txt.txt", metadata=sample_metadata){
  Contaminated_metrics = read_delim(File, delim=" ")
  Contaminated_metrics %>% mutate(Total= Human_reads+Paired_Reads, Cont_percentage=Human_reads/Total) -> Contaminated_metrics
  Contaminated_metrics %>% ggplot() + geom_density(aes(x=Cont_percentage)) + theme_bw() -> contamination_distribution
  ggplot(Contaminated_metrics,aes(x=Total, y=Cont_percentage)) + geom_point() + theme_bw() + geom_smooth(method='lm', formula= y~x) ->contamination_vs_reads

  metadata %>% filter(StudieID %in% Contaminated_metrics$Sample) %>% arrange(StudieID) -> metadata
  Contaminated_metrics %>% filter(Sample %in% metadata$StudieID) %>% arrange(Sample) -> Contaminated_metrics
  Contaminated_metrics %>% mutate(Source = metadata$Source) -> Contaminated_metrics
  Contaminated_metrics %>%
    ggplot(aes(x=Source, y=Cont_percentage)) + geom_boxplot()+ geom_sina()  + theme_bw()
  
  wilcox.test(filter(Contaminated_metrics, Source=="pocket")$Cont_percentage, filter(Contaminated_metrics, Source=="saliva")$Cont_percentage)
  print(contamination_distribution)
  print(contamination_vs_reads)
  
}

distribution_taxonomy = function(taxa_matrix = matrix_taxonomy){
  t_class = taxa_matrix %>% filter(grepl("t__",ID)) %>% mutate(Level="strain")
  s_class = taxa_matrix %>% filter(grepl("s__",ID)) %>% filter(!grepl("t__",ID)) %>% mutate(Level="species")
  g_class = taxa_matrix %>% filter(grepl("g__",ID) & !grepl("s__",ID)) %>% mutate(Level="genus")
  f_class = taxa_matrix %>% filter(grepl("f__",ID) & !grepl("g__",ID)) %>% mutate(Level="family")
  o_class = taxa_matrix %>% filter(grepl("o__",ID) & !grepl("f__",ID)) %>% mutate(Level="order")
  c_class = taxa_matrix %>% filter(grepl("c__",ID) & !grepl("o__",ID)) %>% mutate(Level="clade")
  p_class = taxa_matrix %>% filter(grepl("p__",ID) & !grepl("c__",ID)) %>% mutate(Level="phylum")
  k_class = taxa_matrix %>% filter(grepl("k__",ID) & !grepl("p__",ID)) %>% mutate(Level="domain")
  
  print("Counting number the species identified per taxonomic level")
  
  raw_numbers = tibble(strain = dim(t_class)[1], species = dim(s_class)[1], 
              genus = dim(g_class)[1], family = dim(g_class)[1], order = dim(o_class)[1],
              clade = dim(c_class)[1], phylum = dim(p_class)[1], domain = dim(k_class)[1])
  raw_numbers %>% gather(Taxonomic_level, Number, strain:domain, factor_key=TRUE) -> tax_numbers
  ggplot(tax_numbers) + geom_bar(aes(x=Taxonomic_level,y=Number), stat="identity") + theme_bw() -> taxonomy_numbers
  print(taxonomy_numbers)
  
  
  print("Counting how many items are there per taxonomic level in each individual")
  classified = tibble()
  for (t_s in list(t_class,s_class,g_class,f_class,o_class,c_class,p_class,k_class)){
    taxa = t_s$ID
    transposed = as_tibble(t(select(t_s,-c(ID,Level))))
    colnames(transposed) = taxa
    transposed  %>% mutate_if(is.character, as.numeric) -> transposed
    (! as.matrix(transposed) == 0)*1 -> transposed
    
    classified_i = tibble()
    for (Indv in seq(1:dim(transposed)[1])){
      presences = sum(as_vector(transposed[Indv,]))
      presences = tibble(N=presences, Level=t_s$Level[1], ID =colnames(t_s)[Indv+1])
      classified_i = rbind(classified_i,presences)
    }
    
    classified = rbind(classified,classified_i)
  }
  
  
  classified %>% mutate(ID= str_replace(ID,"_feature_counts","")) %>% filter(!ID=="Pipelines") ->classified
  
  classified %>% arrange(desc(N)) -> ordered
  classified %>% mutate(Level = factor(Level,levels= c("strain","species","genus","family","order","clade","phylum","domain"))) -> classified
  classified %>% mutate(ID = factor(ID,levels= unique(ordered$ID))) -> classified
  
  classified %>% ggplot(aes(x=ID, y=N)) + geom_bar(stat="identity") + facet_wrap(~Level) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) -> number_classified_by_ID
  print(number_classified_by_ID)
  
  
  print("Loading reads per sample, assesing whether the number of reads affect the taxons identified in each sample")
  Reads = read_delim("Reads_contaminated.txt.txt"," ")
  Reads %>% filter(!Sample== "Pipelines") -> Reads
  Read_number = vector()
  for (entry in arrange(Reads,Sample)$Paired_Reads){
    for (i in c("strain","species","genus","family","order","clade","phylum","domain")){ 
      Read_number = c(Read_number,entry)
    }
  }
  arrange(classified,ID) %>% mutate(Reads= Read_number) -> Classified_reads
  ggplot(Classified_reads, aes(x=N,y=Reads,col=ID)) + geom_point() + facet_wrap(~Level, scales="free") + theme_bw()  + theme(legend.position = "none") -> reads_vs_class
  print(reads_vs_class)
}

Preprocess_f = function(Mtx,Summary_name){
#Transpose matrix
Sample = gsub("_feature_counts","",colnames(Mtx))
Sample = Sample[2:length(Sample)]

Taxon = vector()
list_split = strsplit(Mtx$ID, "\\|")
for (i in seq(1:length(list_split))){
  element =  list_split[[i]]
  Taxon = c(Taxon, element[length(element)])
}
#Taxon = Mtx$ID

Mtx[-1,] -> Mtx2
Mtx2 = as_tibble(t(Mtx)[-1,])
Mtx2%>% mutate_if(is.character, as.numeric) -> Mtx2

colnames(Mtx2) = Taxon

do_Metrics(mutate(Mtx2, ID = Sample), N= Summary_name)
##Remove entries that are not seen in 90% of samples
taxa2 = Mtx2[2:length(Mtx2)]
taxa_filt= taxa2 #taxa2[,((colSums(taxa2 !=0) / nrow(taxa2)) *100 )>20]

mutate(taxa_filt, ID = Sample) -> taxa_filt
taxa_filt[c(dim(taxa_filt)[2], seq(1:(dim(taxa_filt)[2]-1)))] -> taxa_filt

return(taxa_filt)
}

Diversity_analysis = function(taxa_matrix,metadata){
  print("Matching data and metadata IDs")
  metadata %>% filter(StudieID %in% taxa_matrix$ID ) %>% arrange(StudieID) -> metadata  
  taxa_matrix %>% filter(ID %in% metadata$StudieID ) %>% arrange(ID) -> taxa_matrix 

  ###Clustering of the data
  ####Heatmap
  print("Building Heatmap")
  rownames(taxa_matrix) = taxa_matrix$ID
  taxa_matrix %>% select(-ID) %>% as.matrix() %>% heatmap()
  #### Hierarchical clustering
  print("Clustering based on disimilarity matrix")
  dist <- vegdist(taxa_matrix[,2:dim(taxa_matrix)[2]],  method = "bray")
  fit = hclust(dist)
  tree = as.phylo(fit)
  tree$tip.label = metadata$StudieID
  ggtree(tree)  %<+%  metadata + geom_tiplab(aes(col=Source),hjust=-0.3) + geom_tippoint(aes(shape=factor(Group)),size=3) + #aes(col=factor(Group)) + 
    xlim(NA,0.7) +  theme(legend.position=c(.9, .6)) -> Cluster_plot
  print(Cluster_plot)
  
  ### Alpha diversity
  print("Computing in sample alpha diversity and plotting")
  shannon_diversity = diversity(taxa_matrix[,2:dim(taxa_matrix)[2]], index = "shannon")
  simpson_diversity = diversity(taxa_matrix[,2:dim(taxa_matrix)[2]], index = "simpson")
  isimpson_diversity = diversity(taxa_matrix[,2:dim(taxa_matrix)[2]], index = "invsimpson")
  alpha_diversities =  mutate(metadata, Shannon=shannon_diversity, Simpson=simpson_diversity, Inverse_Simpson=isimpson_diversity)
  gather(alpha_diversities, Diversity_index, value, c(Shannon,Simpson,Inverse_Simpson), factor_key=TRUE) -> alpha_diversities
  
  alpha_diversities %>% ggplot(aes(x=factor(Group), y=value,col=Source)) + geom_boxplot() +geom_sina() +
  facet_wrap(~Diversity_index, scales= "free") + theme_bw() -> Fig_alpha
  print(Fig_alpha)
  
  ### Beta Diversity
  print("Calculating Beta diversity and computing PCoA")
  dist <- vegdist(taxa_matrix[,2:dim(taxa_matrix)[2]],  method = "bray")
  PCoA = pcoa(dist)

  ggplot(data=as_tibble(PCoA$values)) + geom_bar(aes(x=seq(1:dim(PCoA$values)[1]), y=Relative_eig), col="black", stat = "identity") + 
    theme_bw() -> Fig_Scree_plot
  print(Fig_Scree_plot)
  Eig_vector = PCoA$vectors[,1:3]
  metadata %>% mutate(PCoA1 = Eig_vector[,1], PCoA2 = Eig_vector[,2], PCoA3=Eig_vector[,3]) -> ordination_data
  ggplot(data= ordination_data) + geom_point(aes(x = PCoA1, y=PCoA2, col=Source, shape=factor(Group))) +
    theme_bw() + labs(x=paste(c("PCoA1(",as.character(round(PCoA$values[1,2],3)*100),"%)"), collapse=""), y=paste(c("PCoA2(",as.character(round(PCoA$values[2,2],3)*100),"%)"), collapse="")) -> Fig_PCoA1
  ggplot(data= ordination_data) + geom_point(aes(x = PCoA2, y=PCoA3, col=Source, shape=factor(Group))) +
    theme_bw() + labs(x=paste(c("PCoA2(",as.character(round(PCoA$values[2,2],3)*100),"%)"), collapse=""), y=paste(c("PCoA3(",as.character(round(PCoA$values[3,2],3)*100),"%)"), collapse="")) -> Fig_PCoA2
  print(Fig_PCoA1)
  print(Fig_PCoA2)
  
  ###Correspondence and canonical correspondence analysis
  #D = decorana(taxa_matrix[,2:dim(taxa_matrix)[2]])
  #plot(D)  
  #cbind(metadata,as_tibble(summary(D)$site.scores)) -> DCA_data
  #ggplot(data=DCA_data) + geom_point(aes(x=DCA1, y=DCA2, col=Source, shape=factor(Group))) + theme_bw()+
  #labs(x=paste(c("CA1(",as.character(round(D$evals[1],3)*100),"%)"), collapse=""), y=paste(c("CA2(",as.character(round(D$evals[2],3)*100),"%)"), collapse=""))
    
  print("Making Canonical Correspondence Analysis")
  ord <- cca(taxa_matrix[,2:dim(taxa_matrix)[2]] ~ Source + factor(Group), data=metadata)
  Info_ord = summary(ord)
  metadata %>% mutate(CCA1= as.vector(Info_ord$sites[,1]), CCA2= as.vector(Info_ord$sites[,2]), CC1_saliva=as.vector(Info_ord$biplot[1,1]), CC2_saliva= as.vector(Info_ord$biplot[1,2]),CC1_patient= as.vector(Info_ord$biplot[2,1]),CC2_patient=as.vector(Info_ord$biplot[2,2])) -> metadata_coa
  
  ggplot(metadata_coa) + geom_point(aes(x=CCA1, y=CCA2, col=Source, shape=factor(Group))) +
    geom_label(label="Saliva",aes(x=CC1_saliva,y=CC2_saliva)) + geom_label(label="Gengivitis",aes(x=CC1_patient,y=CC2_patient)) + theme_bw() -> Fig_CC
  print(Fig_CC)
  
  Reads = read_delim("Reads_contaminated.txt.txt",delim=" ")
  metadata %>% mutate(n_reads = arrange(Reads,Sample)$Paired_Reads) -> metadata
  ###Alpha diversity statistical testing
  print("Using wilcoxon rank test on saliva difference on alpha diversity between patients and controls")
  alpha_diversities %>% filter(Source=="saliva") -> Saliva_diversity
  ggplot(data=Saliva_diversity) + geom_density(aes(x=value, col=factor(Group))) + facet_wrap(~Diversity_index, scales="free") + theme_bw()
  for (Diversity in unique(Saliva_diversity$Diversity_index)){
    Saliva_diversity %>% filter(Diversity_index == Diversity) -> MS
    Test = wilcox.test(filter(MS, Group=="0")$value, filter(MS, Group=="1")$value)
    print(c(Diversity,Test$p.value))
  }
  print("Using wilcoxon rank test on pocket difference on alpha diversity between patients and controls")
  alpha_diversities %>% filter(Source=="pocket") -> Pocket_diversity
  ggplot(data=Pocket_diversity) + geom_density(aes(x=value, col=factor(Group))) + facet_wrap(~Diversity_index, scales="free") + theme_bw()
  for (Diversity in unique(Saliva_diversity$Diversity_index)){
    Pocket_diversity %>% filter(Diversity_index == Diversity) -> MS
    Test = wilcox.test(filter(MS, Group=="0")$value, filter(MS, Group=="1")$value)
    print(c(Diversity,Test$p.value))
  }
  #######
  
  ###Beta diversity Statistical testing
  print("Doing Permanova test on Beta diversity")
  permanova_1 = adonis2(dist~ Source + factor(Group) +n_reads,data=metadata, permutations = 999,strata=Source, method="bray", strata="PLOT")
  
  #Divide by tissue
  metadata %>% filter(Source=="pocket") -> meta_pocket
  metadata %>% filter(Source=="saliva") -> meta_saliva
  taxa_matrix %>% filter(ID %in% meta_saliva$StudieID) -> saliva_matrix
  taxa_matrix %>% filter(ID %in% meta_pocket$StudieID) -> pocket_matrix
  saliva_dist <- vegdist(saliva_matrix[,2:dim(saliva_matrix)[2]],  method = "bray")
  pocket_dist <- vegdist(pocket_matrix[,2:dim(pocket_matrix)[2]],  method = "bray")
  
  permanova_2_saliva = adonis2(saliva_dist~ factor(Group)+n_reads ,data=meta_saliva, permutations = 999,strata=Source, method="bray", strata="PLOT")
  permanova_2_pocket = adonis2(pocket_dist~ factor(Group)+n_reads ,data=meta_pocket, permutations = 999,strata=Source, method="bray", strata="PLOT")
  print("Total permanova")
  print(permanova_1)
  print("In saliva permanova")
  print(permanova_2_saliva)
  print("In pocket permanova")
  print(permanova_2_pocket)
  ###########
}

fit_linear_model = function(abundance_matrix,metadata, covariates){
  result_bugg = tibble()
  for (N in seq(1:dim(abundance_matrix)[2])){
    dependent = as_vector(abundance_matrix[,N])
    if (sum(dependent)==0){next
    }else if (sum((!dependent==0)*1)/sum(dependent) < 0.2){ next }
    
    metadata %>% mutate(Dependent = dependent) -> regression_input
    distribution_bug = ggplot(regression_input) + geom_density(aes(Dependent)) + theme_bw() + facet_wrap(~Source)
    
    if ("source" %in% covariates & length(covariates) == 1){ 
      lm(Dependent ~ Source + factor(Group) + n_reads, data= regression_input ) -> regression_output
      S = summary(regression_output)
      fold_change = mean(filter(regression_input, Group == 1)$Dependent) /  mean(filter(regression_input, Group == 0)$Dependent)
      pvalue = as_tibble(S$coefficients)$`Pr(>|t|)`[3]
   }else if(length(covariates)==0){
      lm(Dependent ~ factor(Group) + n_reads, data= regression_input ) -> regression_output
      S = summary(regression_output)
      fold_change = mean(filter(regression_input, Group == 1)$Dependent) /  mean(filter(regression_input, Group == 0)$Dependent)
      pvalue = as_tibble(S$coefficients)$`Pr(>|t|)`[2]
   }
    
    
    
    #qqnorm(S$residuals)
    #qqline(S$residuals, col = "steelblue", lwd = 2)
    
    shapiro.test(S$residuals) -> normality
    Normality_pval = normality$p.value
    
    f_t = as_tibble(t(c(colnames(abundance_matrix)[N], pvalue, fold_change, Normality_pval)))
    colnames(f_t) = c("Taxa","P-value","Fold_change","Normality_pvalue")
    result_bugg = rbind(result_bugg, f_t)
  }
  result_bugg %>% mutate(FDR=p.adjust(`P-value`,"fdr"), FDR_normality=p.adjust(Normality_pvalue,"fdr")) -> result_bugg
  return(result_bugg)
}

non_parametric_testing = function(abundance_matrix,metadata){
  result_bugg = tibble()
  for (N in seq(1:dim(abundance_matrix)[2])){
    dependent = as_vector(abundance_matrix[,N])
    
    if (sum(dependent)==0){next
    }else if (sum((!dependent==0)*1)/sum(dependent) < 0.2){ next }
    
    metadata %>% mutate(Dependent = dependent) -> regression_input
    lm(Dependent ~  n_reads, data= regression_input ) -> regression_output
    Read_free_abundance = regression_output$residuals
    
    regression_input %>% mutate(Abundance2 = Read_free_abundance) -> regression_input
    wilcox.test(filter(regression_input, Group == 1)$Abundance2, filter(regression_input, Group == 0)$Abundance2) -> non_param
    
    fold_change = mean(filter(regression_input, Group == 1)$Dependent)/mean(filter(regression_input, Group == 0)$Dependent)
    bug = as_tibble(t(c(colnames(abundance_matrix)[N],non_param$p.value, fold_change)))
    colnames(bug) = c("Taxa","P-value","Fold_change")
    
    result_bugg = rbind(result_bugg,bug)
    
  }
  result_bugg %>% mutate(FDR= p.adjust(`P-value`,"fdr")) -> result_bugg
  return(result_bugg)
}


gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }
inverse_rank = function(x){ return(rank(-x)) }
ratio_transform = function(Vector){
  transformed_vector = log(Vector/gm_mean(Vector))
  return(transformed_vector)
}

Differential_abundance = function(taxa_matrix,metadata, transformation = "Square"){
  metadata %>% filter(StudieID %in% taxa_matrix$ID) %>% arrange(StudieID) -> metadata
  Reads = read_delim("Reads_contaminated.txt.txt",delim=" ")
  metadata %>% mutate(n_reads = arrange(Reads,Sample)$Paired_Reads) -> metadata
  taxa_matrix %>% filter(ID %in% metadata$StudieID) %>% arrange(ID) -> taxa_matrix
  abundance_matrix = taxa_matrix[,2:dim(taxa_matrix)[2]]
  
  print("Filtering Columns absent in over 20% of samples")
  abundance_matrix = abundance_matrix[,((colSums(abundance_matrix !=0) / nrow(abundance_matrix)) *100 )>20]
  
  print("Transforming data")
  if (transformation == "Square"){
    abundance_matrix = sqrt(abundance_matrix)
  } else if(transformation == "argsinsquare"){
    abundance_matrix = asin(sqrt(abundance_matrix/100))
  } else if(transformation == "None"){ print("No transformation applied")
  } else if(transformation == "ratio"){ abundance_matrix = as_tibble(sapply(X = abundance_matrix, FUN = ratio_transform)) 
  } else if(transformation == "inverse_rank"){ abundance_matrix = as_tibble(sapply(X = abundance_matrix, FUN = inverse_rank))
  } else if(transformation == "log"){ abundance_matrix = log(abundance_matrix+1)}
  
  print("Iterating by taxa and fitting a linear model")
  Model_1 = fit_linear_model(abundance_matrix, metadata, "source")
  
  #print(Model_1)
  
  ##Tissue-specific models
  as_tibble(abundance_matrix) %>% mutate(ID=metadata$StudieID) -> temporal_matrix
  metadata %>% filter(Source=="pocket") -> meta_pocket
  metadata %>% filter(Source=="saliva") -> meta_saliva
  temporal_matrix %>% filter(ID %in% meta_pocket$StudieID) %>% select(-ID) -> pocket_matrix
  temporal_matrix %>% filter(ID %in% meta_saliva$StudieID) %>% select(-ID) -> saliva_matrix
  
  print("Iterating by taxa and fitting a linear model on pocket samples")
  Model_2 = fit_linear_model(pocket_matrix, meta_pocket, vector())
  Model_21 = non_parametric_testing(pocket_matrix, meta_pocket)
  #print(Model_2)
  print("Iterating by taxa and fitting a linear model on saliva samples")
  Model_3 = fit_linear_model(saliva_matrix, meta_saliva, vector())
  Model_31 = non_parametric_testing(saliva_matrix, meta_saliva)
  #print(Model_3)
  
  print(c(dim(filter(Model_1,FDR_normality > 0.1))[1],dim(filter(Model_2,FDR_normality > 0.1))[1],dim(filter(Model_3,FDR_normality > 0.1))[1], dim(abundance_matrix)[2]))
  
  print(c(dim(filter(Model_21,FDR < 0.1)), dim(filter(Model_31,FDR<0.1))))
  return(Model_3)
}

do_bars = function(taxa_matrix, metadata){
  metadata %>% filter(StudieID %in% taxa_matrix$ID) %>% arrange(StudieID) -> metadata
  taxa_matrix %>% filter(ID %in% metadata$StudieID) %>% arrange(ID) -> taxa_matrix
  taxa_matrix = select(taxa_matrix, -ID)
  
  metadata %>% mutate(N_ID=as.numeric(gsub("[^0-9.-]", "", StudieID))) %>% arrange(N_ID) -> test
  metadata$StudieID = factor(metadata$StudieID,levels(factor(metadata$StudieID))[match(test$StudieID,metadata$StudieID)])
  
  cbind(metadata,taxa_matrix) -> combined_matrix
  
  #tissue
  print("Study taxa composition by tissue")
  combined_matrix %>% gather(Taxa, Abundance, (dim(metadata)[2]+1):dim(combined_matrix)[2], factor_key=TRUE) %>%
  ggplot() + geom_bar(aes(x=StudieID,y=Abundance, fill=Taxa), stat = "identity",position="stack") + facet_wrap(~Source, dir="v")+
   theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",legend.title = element_text(size = 3), legend.text = element_text(size = 6))+
    theme(legend.key.size = unit(0.1, "cm"))
 
  print("Study taxa composition in Patients vs Controls")
  combined_matrix %>% gather(Taxa, Abundance, (dim(metadata)[2]+1):dim(combined_matrix)[2], factor_key=TRUE) %>%
    ggplot() + geom_bar(aes(x=StudieID,y=Abundance, fill=Taxa), stat = "identity",position="stack") + facet_wrap(~factor(Group), dir="v")+
    theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",legend.title = element_text(size = 3), legend.text = element_text(size = 6))+
    theme(legend.key.size = unit(0.1, "cm"))
  
  
  print("Top 5 taxons by tissue")
  combined_matrix %>% filter(Source=="saliva") %>% select(colnames(taxa_matrix)) %>% summarise_all(median) %>%
    gather(Taxa, Median_relative_abundance, 1:length(colnames(taxa_matrix)), factor_key=TRUE) %>%
  arrange(desc(Median_relative_abundance)) %>% head(n=5)
   
  combined_matrix %>% filter(Source=="pocket") %>% select(colnames(taxa_matrix)) %>% summarise_all(median) %>%
    gather(Taxa, Median_relative_abundance, 1:length(colnames(taxa_matrix)), factor_key=TRUE) %>%
    arrange(desc(Median_relative_abundance)) %>% head(n=5)
  
  print("Top 5 taxons by condition")
  
  combined_matrix %>% filter(as.character(Group)=="0") %>% select(colnames(taxa_matrix)) %>% summarise_all(median) %>%
    gather(Taxa, Median_relative_abundance, 1:length(colnames(taxa_matrix)), factor_key=TRUE) %>%
    arrange(desc(Median_relative_abundance)) %>% head(n=5)
  
  combined_matrix %>% filter(as.character(Group)=="1") %>% select(colnames(taxa_matrix)) %>% summarise_all(median) %>%
    gather(Taxa, Median_relative_abundance, 1:length(colnames(taxa_matrix)), factor_key=TRUE) %>%
    arrange(desc(Median_relative_abundance)) %>% head(n=5)
  
  
  
}



### Exploratory analysis
contamination_stats()
distribution_taxonomy()

Preprocess_f(family_matrix,"Summary_family.csv") -> family_matrix2
do_bars(family_matrix2, sample_metadata)
###


### Statistical analysis: Clean, Diversity, indv taxa effect

Preprocess_f(genus_matrix,"Summary_genus.csv") -> genus_matrix2
Diversity_analysis(genus_matrix2,sample_metadata)
Differential_abundance(genus_matrix2,sample_metadata,"Square")

Preprocess_f(species_matrix,"Summary_species.csv") -> species_matrix2
Diversity_analysis(species_matrix2,sample_metadata)
Differential_abundance(species_matrix2,sample_metadata,"argsinsquare")

Diversity_analysis(family_matrix2,sample_metadata)
Differential_abundance(family_matrix2,sample_metadata,"argsinsquare")

#At the family level f__Propionibacteriaceae appears significant both with the linear model (residuals are normal) and in the wilcoxon test
#Fold change = 2.4, meaning that that family appears two times more relative abundant


species_matrix2 %>% arrange(ID) %>% filter(ID %in% sample_metadata$StudieID) -> species_matrix2
sample_metadata %>% arrange(StudieID) %>% filter(StudieID %in% species_matrix2$ID) -> sample_metadata

species_matrix2 %>% mutate(Group = sample_metadata$Group, Source=sample_metadata$Source) %>% filter(Source=="pocket") -> pocket_samples
pocket_samples %>% group_by(Group) %>% summarise_if(is.numeric,median) %>%
  gather(Microorganism, Median_abundance, 2:235, factor_key=TRUE) -> Median_species
Median_species %>% ggplot(aes(x=factor(Group), y=Median_abundance)) + geom_point() + theme_bw()

bigger_changes = tibble()
for (i in unique(Median_species$Microorganism)){
  Median_species %>% filter(Microorganism==i) -> M
  filter(M, Group==1)$Median_abundance/filter(M, Group==0)$Median_abundance -> Change
  ti = tibble(Micro=M$Microorganism[1], Log_fold_Change=log(Change))
  bigger_changes = rbind(bigger_changes,ti)
}
bigger_changes %>% arrange(desc(Log_fold_Change))


Median_species %>% filter(Group==0) %>% arrange(desc(Median_abundance)) %>% print(n=30)
Median_species %>% filter(Group==1) %>% arrange(desc(Median_abundance))

pocket_samples %>% select(s__Porphyromonas_gingivalis, ID, Group) -> por_ging_pocket
wilcox.test(filter(por_ging_pocket, Group==1)$s__Porphyromonas_gingivalis,filter(por_ging_pocket, Group==0)$s__Porphyromonas_gingivalis)




#############################
######PATHWAY ANALYSIS#######
#############################

differential_pathway_analysis = function(Pathway_abundance, metadata){
  metadata %>% filter(StudieID %in% Pathway_abundance$ID) %>% arrange(StudieID) -> metadata
  Pathway_abundance %>% filter(ID %in% metadata$StudieID) %>% arrange(ID) -> Pathway_abundance
  
  Reads = read_delim("Reads_contaminated.txt.txt",delim=" ")
  metadata %>% mutate(n_reads = arrange(Reads,Sample)$Paired_Reads) -> metadata
  
  
  Abundance_matrix = select(Pathway_abundance, -c(ID))
  #Transformation
  Abundance_matrix = asin(sqrt(Abundance_matrix/100))
  
  
  Model_1 = fit_linear_model(abundance_matrix = Abundance_matrix, metadata = metadata, covariates = "source")
  
  ##Tissue-specific models
  as_tibble(Abundance_matrix) %>% mutate(ID=metadata$StudieID) -> temporal_matrix
  metadata %>% filter(Source=="pocket") -> meta_pocket
  metadata %>% filter(Source=="saliva") -> meta_saliva
  temporal_matrix %>% filter(ID %in% meta_pocket$StudieID) %>% select(-ID) -> pocket_matrix
  temporal_matrix %>% filter(ID %in% meta_saliva$StudieID) %>% select(-ID) -> saliva_matrix
  
  print("Iterating by taxa and fitting a linear model on pocket samples")
  Model_2 = fit_linear_model(pocket_matrix, meta_pocket, vector())
  Model_21 = non_parametric_testing(pocket_matrix, meta_pocket)
  #print(Model_2)
  print("Iterating by taxa and fitting a linear model on saliva samples")
  Model_3 = fit_linear_model(saliva_matrix, meta_saliva, vector())
  Model_31 = non_parametric_testing(saliva_matrix, meta_saliva)
  as_tibble(Model_2) %>% group_by(FDR<0.01, FDR_normality>0.01) %>% summarise(n()) -> C
  as_tibble(Model_3)%>% group_by(FDR<0.01, FDR_normality>0.01) %>% summarise(n()) -> C2
  
  print(c(C,C2))
  
}
make_relative = function(Column){
  return(Column/sum(Column)*100)
}


pathway_abundance = read_tsv("Merged_batch_pathabundance.tsv")

pathway_abundance %>% filter(! grepl("\\|",`# Pathway`) )-> pathway_abundance

new_i = vector()
for (i in colnames(pathway_abundance)){
  if (i=="# Pathway"){ 
    i = "Pathway"
  }else{
    i = str_replace(string = i,pattern = "_merged_Abundance",replacement = "")
  }
  new_i = c(new_i,i)
}

colnames(pathway_abundance) = new_i
pathway_abundance %>% select(-Pipelines) -> pathway_abundance
pathway_abundance %>% mutate_if(is_numeric,make_relative) -> pathway_abundance

new_i = colnames(pathway_abundance)
Pathway = pathway_abundance$Pathway
pathway_abundance %>% select(-Pathway) %>% t() -> pathway_abundance2

colnames(pathway_abundance2) = Pathway
pathway_abundance2 %>% as_tibble() %>% mutate(ID = new_i[2:length(new_i)]) -> pathway_abundance2
pathway_abundance2 %>% select(-c(UNMAPPED,UNINTEGRATED)) -> pathway_abundance2
pathway_abundance2 %>% select(c("ID",colnames(select(pathway_abundance2,-ID)))) -> pathway_abundance2



Diversity_analysis(pathway_abundance2,sample_metadata)
differential_pathway_analysis(pathway_abundance2,sample_metadata)


################Compositional Modelling#################

###NOT WORKING
fit_aldex = function(abundance_matrix, metadata, covariates){
  result_bugg = tibble()
  for (N in seq(1:dim(abundance_matrix)[2])){
    dependent = as_vector(abundance_matrix[,N])
    if (sum(dependent)==0){next}
    
    metadata %>% mutate(Dependent = dependent) -> regression_input
    distribution_bug = ggplot(regression_input) + geom_density(aes(Dependent)) + theme_bw() + facet_wrap(~Source)
    
    if ("source" %in% covariates & length(covariates) == 1){ 
      lm(Dependent ~ Source + factor(Group) + n_reads, data= regression_input ) -> regression_output
      S = summary(regression_output)
      fold_change = mean(filter(regression_input, Group == 1)$Dependent) /  mean(filter(regression_input, Group == 0)$Dependent)
      pvalue = as_tibble(S$coefficients)$`Pr(>|t|)`[3]
    }else if(length(covariates)==0){
      abundance_matrix %>% as.data.frame() -> check
      t(check) -> check
      colnames(check) = metadata$StudieID
      mm <- model.matrix(~ factor(Group) + n_reads , data= metadata)
      x <- aldex.clr(check, mm, mc.samples=8, denom="all")
      glm.test <- aldex.glm(x,mm)
      
      
      S = summary(regression_output)
      fold_change = mean(filter(regression_input, Group == 1)$Dependent) /  mean(filter(regression_input, Group == 0)$Dependent)
      pvalue = as_tibble(S$coefficients)$`Pr(>|t|)`[2]
    }
    
    shapiro.test(S$residuals) -> normality
    Normality_pval = normality$p.value
    
    f_t = as_tibble(t(c(colnames(abundance_matrix)[N], pvalue, fold_change, Normality_pval)))
    colnames(f_t) = c("Taxa","P-value","Fold_change","Normality_pvalue")
    result_bugg = rbind(result_bugg, f_t)
  }
  result_bugg %>% mutate(FDR=p.adjust(`P-value`,"fdr"), FDR_normality=p.adjust(Normality_pvalue,"fdr")) -> result_bugg
  return(result_bugg)
}

  

compotisional_diversity = function(taxa_matrix,metadata){
  print("Matching data and metadata IDs")
  metadata %>% filter(StudieID %in% taxa_matrix$ID ) %>% arrange(StudieID) -> metadata  
  taxa_matrix %>% filter(ID %in% metadata$StudieID ) %>% arrange(ID) -> taxa_matrix 
  
  
  ### Beta Diversity
  print("PCA from compositional data")
  dist <- vegdist(taxa_matrix[,2:dim(taxa_matrix)[2]],  method = "euclidean")
  PCA = pcoa(dist)
  
  ggplot(data=as_tibble(PCA$values)) + geom_bar(aes(x=seq(1:dim(PCA$values)[1]), y=Relative_eig), col="black", stat = "identity") + 
    theme_bw() -> Fig_Scree_plot
  print(Fig_Scree_plot)
  Eig_vector = PCA$vectors[,1:3]
  metadata %>% mutate(PCA1 = Eig_vector[,1], PCA2 = Eig_vector[,2], PCA3=Eig_vector[,3]) -> ordination_data
  ggplot(data= ordination_data) + geom_point(aes(x = PCA1, y=PCA2, col=Source, shape=factor(Group))) +
    theme_bw() + labs(x=paste(c("PCA1(",as.character(round(PCA$values[1,2],3)*100),"%)"), collapse=""), y=paste(c("PCA2(",as.character(round(PCoA$values[2,2],3)*100),"%)"), collapse="")) -> Fig_PCoA1
  ggplot(data= ordination_data) + geom_point(aes(x = PCA2, y=PCA3, col=Source, shape=factor(Group))) +
    theme_bw() + labs(x=paste(c("PCA2(",as.character(round(PCA$values[2,2],3)*100),"%)"), collapse=""), y=paste(c("PCA3(",as.character(round(PCoA$values[3,2],3)*100),"%)"), collapse="")) -> Fig_PCoA2
  print(Fig_PCoA1)
  print(Fig_PCoA2)
  
  Reads = read_delim("Reads_contaminated.txt.txt",delim=" ")
  metadata %>% mutate(n_reads = arrange(Reads,Sample)$Paired_Reads) -> metadata
  
  print("Doing Permanova test on Beta diversity")
  permanova_1 = adonis2(dist~ Source + factor(Group) +n_reads,data=metadata, permutations = 999,strata=Source, method="euclidean", strata="PLOT")
  
  #Divide by tissue
  metadata %>% filter(Source=="pocket") -> meta_pocket
  metadata %>% filter(Source=="saliva") -> meta_saliva
  taxa_matrix %>% filter(ID %in% meta_saliva$StudieID) -> saliva_matrix
  taxa_matrix %>% filter(ID %in% meta_pocket$StudieID) -> pocket_matrix
  saliva_dist <- vegdist(saliva_matrix[,2:dim(saliva_matrix)[2]],  method = "euclidean")
  pocket_dist <- vegdist(pocket_matrix[,2:dim(pocket_matrix)[2]],  method = "euclidean")
  
  permanova_2_saliva = adonis2(saliva_dist~ factor(Group)+n_reads ,data=meta_saliva, permutations = 999,strata=Source, method="bray", strata="PLOT")
  permanova_2_pocket = adonis2(pocket_dist~ factor(Group)+n_reads ,data=meta_pocket, permutations = 999,strata=Source, method="bray", strata="PLOT")
  print("Total permanova")
  print(permanova_1)
  print("In saliva permanova")
  print(permanova_2_saliva)
  print("In pocket permanova")
  print(permanova_2_pocket)
  
  total = fit_aldex(select(taxa_matrix,-ID), c("source"))
  saliva = fit_aldex(select(saliva_matrix,-ID), vector())
  pocket = fit_aldex(select(pocket_matrix,-ID), vector())
  
}  




library(compositions)
library(ALDEx2)
species_matrix2 %>% mutate_if(is.numeric, clr) -> species_compositional
compositional_diversity(species_compositional, sample_metadata)
