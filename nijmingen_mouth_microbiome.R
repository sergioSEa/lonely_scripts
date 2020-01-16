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


###Functions
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
taxa_filt=taxa2[,((colSums(taxa2 !=0) / nrow(taxa2)) *100 )>20]

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
  
  permanova_1 = adonis2(dist~ Source + factor(Group),data=metadata, permutations = 999,strata=Source, method="bray", strata="PLOT")
  
  #Divide by tissue
  metadata %>% filter(Source=="pocket") -> meta_pocket
  metadata %>% filter(Source=="saliva") -> meta_saliva
  taxa_matrix %>% filter(ID %in% meta_saliva$StudieID) -> saliva_matrix
  taxa_matrix %>% filter(ID %in% meta_pocket$StudieID) -> pocket_matrix
  saliva_dist <- vegdist(saliva_matrix[,2:dim(saliva_matrix)[2]],  method = "bray")
  pocket_dist <- vegdist(pocket_matrix[,2:dim(pocket_matrix)[2]],  method = "bray")
  
  permanova_2_saliva = adonis2(saliva_dist~ factor(Group) ,data=meta_saliva, permutations = 999,strata=Source, method="bray", strata="PLOT")
  permanova_2_pocket = adonis2(pocket_dist~ factor(Group) ,data=meta_pocket, permutations = 999,strata=Source, method="bray", strata="PLOT")
  
  print(permanova_2_saliva)
  print(permanova_2_pocket)
  ###########
}


Differential_abundance = function(taxa_matrix,metadata, transformation = "Square"){
  abundance_matrix = taxa_matrix[,2:dim(taxa_matrix)[2]]
  print("Transforming data")
  if (transformation == "Square"){
    abundance_matrix = sqrt(abundance_matrix)
  } else if(transformation == "argsinsquare"){
    abundance_matrix = asin(sqrt(abundance_matrix/100))
  } else if(transformation == "None"){}
  
  print("Iterating by taxa and fitting a linear model")
  for (N in seq(1:dim(abundance_matrix)[2])){
    dependent = as_vector(abundance_matrix[,N])
    metadata %>% mutate(Dependent = dependent) -> regression_input
    distribution_bug = ggplot(regression_input) + geom_density(aes(Dependent)) + theme_bw() + facet_wrap(~Source)
    lm(Dependent ~ Source + factor(Group), data= regression_input ) -> regression_output
    S = summary(regression_output)  
    qqnorm(S$residuals)
    qqline(S$residuals, col = "steelblue", lwd = 2)
    shapiro.test(S$residuals) -> normality
    Normality_pval = normality$p.value
    
    print(c(colnames(abundance_matrix)[N], Normality_pval))
  }
  
}

Study_factor_effects = function(taxa_matrix, metadata){
  metadata %>% filter(StudieID %in% taxa_matrix$ID) %>% arrange(StudieID) -> metadata
  taxa_matrix %>% filter(ID %in% metadata$StudieID) %>% arrange(ID) -> taxa_matrix
  taxa_matrix = select(taxa_matrix, -ID)
  
  metadata %>% mutate(N_ID=as.numeric(gsub("[^0-9.-]", "", StudieID))) %>% arrange(N_ID) -> test
  metadata$StudieID = factor(metadata$StudieID,levels(factor(metadata$StudieID))[match(test$StudieID,metadata$StudieID)])
  
  cbind(metadata,taxa_matrix) -> combined_matrix
  
  #tissue
  combined_matrix %>% gather(Taxa, Abundance, (dim(metadata)[2]+1):dim(combined_matrix)[2], factor_key=TRUE) %>%
  ggplot() + geom_bar(aes(x=StudieID,y=Abundance, fill=Taxa), stat = "identity",position="stack") + facet_wrap(~Source, dir="v")+
   theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",legend.title = element_text(size = 3), legend.text = element_text(size = 6))+
    theme(legend.key.size = unit(0.1, "cm"))#+
#guides(fill = F )
  
}


Preprocess_f(genus_matrix,"Summary_genus.csv") -> genus_matrix2
Diversity_analysis(genus_matrix2,sample_metadata)

Preprocess_f(species_matrix,"Summary_species.csv") -> species_matrix2
Diversity_analysis(species_matrix2,sample_metadata)

Preprocess_f(family_matrix,"Summary_family.csv") -> family_matrix2
Study_factor_effects(family_matrix2, sample_metadata)
