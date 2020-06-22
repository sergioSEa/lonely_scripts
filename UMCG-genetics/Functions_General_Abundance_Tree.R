#General alpha and beta diversity analysis
#Differential abundance
#Tree

setwd("~/PhD/WORK/Collaborations/Debby/PPI")
#General and visualization
library(tidyverse)
library(pheatmap)
library(patchwork)
library(vegan)
library(dunn.test)
library(reshape2)
#Tree
library(ape)
library(ggtree)

source("/Users/Sergio/Documents/PhD/WORK/Collaborations/Debby/PPI/Data/ancom_scripts/ancom_v2.1.R")


#Open Data and separate Group in two different features
#Abundance = read_csv("Data/Merged_abundance.csv")
Abundance = read_csv("Data/Merged_Genus_abundance.csv")

Abundance %>% arrange(ID) -> Abundance
metadata = read_tsv("Data/sample.list")
metadata %>% arrange(Sample) -> metadata
metadata %>% mutate(Origin = ifelse(grepl("AGEO",Group), "Oral", "Gut")) %>% mutate(Treatment=ifelse(grepl("CTRL", Group), 0, 1)) -> metadata
#Data used for Tree
Tree = read.tree("Data/16S.tree")
#Tree = read.tree("Data/16S.tree") 

#All_samples = read_csv("Data/Matrix_abundance.csv")#This is the data without merging it into higher taxonomical levels
All_samples = read_csv("Data/Matrix_abundance.csv")
#####################
#####################
Preliminary_checks = function(Abundance,metadata){
  #Number of reads per sample
  apply(X = select(Abundance,-ID), MARGIN = 1, FUN = sum) -> Total_number
  Abundance %>% mutate(Total_number = Total_number) -> Abundance_w_number
  ggplot(Abundance_w_number) + geom_bar(aes(x=reorder(ID, -Total_number), y=Total_number), fill="#0080FF",col="black", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1)) -> Total_reads_plot
  #Number of 0s per sample
  apply(X = select(Abundance,-ID), MARGIN = 1, FUN = function(x){ 100*(sum(x==0)/sum(x)) } ) -> Percentage_0
  Abundance_w_number %>% mutate(Percentage_0 = Percentage_0) -> Abundance_w_number
  ggplot(Abundance_w_number) + geom_bar(aes(x=reorder(ID, -Total_number), y=Percentage_0), fill="#00FFFF", col="black", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(size = 7,angle = 90, hjust = 1)) -> Amount_0_plot
  #Both visualizations
  print(Total_reads_plot + Amount_0_plot)
  return(Total_number)
}
Alpha_diversity = function(M, metadata){
  #Shanon diversity calculation per sample, plots and non-parametric test
  library(ggforce)

  shannon_diversity = vegan::diversity(select(M,-ID), index = "shannon")
  metadata %>% mutate(Diversity = shannon_diversity) -> metadata
  
  metadata %>% ggplot(aes(x=factor(Treatment), y=Diversity,col=Origin)) + geom_boxplot(outlier.size = NA) +geom_sina(aes(size=Reads)) +
    theme_bw() -> Fig_alpha
  print(Fig_alpha)

  #kruskal.test(Diversity ~ Group, data=metadata)
  dunn.test(x=metadata$Diversity, g= metadata$Group)
}
Beta_diversity = function(M, metadata, Input="standard"){
  print("Calculating Beta diversity and computing PCoA")
  if (Input=="CLR"){
    dist = vegdist(select(M,-ID),  method = "euclidean")
  } else{
    dist <- vegdist(select(M,-ID),  method = "bray")
  }
  PCoA = pcoa(dist)
  
  Compare_distance_distributions = function(){
    Distance = as.matrix(dist)
    Distance[lower.tri(Distance)] = 0 -> Distance
    as_tibble(Distance) -> Distance2
  
    Distance2 %>% mutate(Group = metadata$Group) -> Distance2
    Comparison_intra = tibble()
    for (G in unique(metadata$Group)){
      Distance2 %>% filter(Group == G) %>% select(-Group) -> sub_distance
      column_group = which(metadata$Group == G)
      as.vector(as_vector(sub_distance[,column_group])) -> sub_distance
      sub_distance[sub_distance == 0] = NA
      To_tibble = tibble(Group = G, Values = sub_distance)
      Comparison_intra = rbind(Comparison_intra, To_tibble)
    }
    Comparison_intra %>% drop_na() -> Comparison_intra
    ggplot(Comparison_intra,aes(x= Group, y=Values)) + geom_boxplot() + theme_bw()
    dunn.test(Comparison_intra$Values,Comparison_intra$Group)
  }  
  ggplot(data=as_tibble(PCoA$values)) + geom_bar(aes(x=seq(1:dim(PCoA$values)[1]), y=Relative_eig), col="black", stat = "identity") + 
    theme_bw() -> Fig_Scree_plot
  print(Fig_Scree_plot)
  Eig_vector = PCoA$vectors[,1:3]
  metadata %>% mutate(PCoA1 = Eig_vector[,1], PCoA2 = Eig_vector[,2], PCoA3=Eig_vector[,3]) -> ordination_data
  ggplot(data= ordination_data) + geom_point(aes(x = PCoA1, y=PCoA2, col=factor(Treatment), shape=Origin)) +
    theme_bw() + labs(x=paste(c("PCoA1(",as.character(round(PCoA$values[1,2],3)*100),"%)"), collapse=""), y=paste(c("PCoA2(",as.character(round(PCoA$values[2,2],3)*100),"%)"), collapse="")) -> Fig_PCoA1
  ggplot(data= ordination_data) + geom_point(aes(x = PCoA2, y=PCoA3, col=factor(Treatment), shape=Origin)) +
   theme_bw() + labs(x=paste(c("PCoA2(",as.character(round(PCoA$values[2,2],3)*100),"%)"), collapse=""), y=paste(c("PCoA3(",as.character(round(PCoA$values[3,2],3)*100),"%)"), collapse="")) -> Fig_PCoA2
  print(Fig_PCoA1)
  print(Fig_PCoA2)
  
  print("Permanova")
  if (Input=="CLR"){ 
    permanova_1 = adonis2(dist~  Origin + factor(Treatment) + Reads,data=metadata, permutations = 2000, method="euclidean")
  }else{  
    permanova_1 = adonis2(dist~  Origin*factor(Treatment) + Reads,data=metadata, permutations = 2000, method="bray")
  }
  print(permanova_1)
}
Beta_diversity_unifract = function(M=All_samples, metadata=metadata, Tree=Tree){
  library(phyloseq)
  #Prepare Phyloseq
  sapply(str_replace(M$ID,"FDMP19H000760-1a_L1_", ""), FUN = function(x){paste(strsplit(x,"")[[1]][1:8],collapse="")}) -> New_names
  M %>% mutate(ID=as.vector(New_names)) -> M

  M %>% select(-ID) %>% t() %>% as_tibble() %>% mutate(label=colnames(select(M,-ID))) %>% `colnames<-`(c(M$ID, "label")) %>% select(c("label",M$ID)) -> trans_M
  trans_M %>% as.data.frame() %>% column_to_rownames("label") -> ASV_table
  taxmat = matrix(rownames(ASV_table))
  rownames(taxmat) = rownames(ASV_table)
  sampledata = metadata %>% mutate(Treatment = as.factor(Treatment)) %>% as.data.frame() %>% column_to_rownames("Sample")
  
  otu_table(ASV_table, taxa_are_rows = TRUE) -> ASV_table
  tax_table(taxmat) -> taxmat
  phyloseq::phyloseq(ASV_table, taxmat) -> phylobjetc
  phyloseq::merge_phyloseq(phylobjetc,Tree, sample_data(sampledata)) -> phylobjetc
  
  phylobject_norm = rarefy_even_depth(phylobjetc, rngseed = 711, replace=F)
  
  #ordi = ordinate(phylobjetc, method="PCoA", distance="unifrac", weighted=T)
  #ordi2 = ordinate(phylobjetc, method="PCoA", distance="unifrac", weighted=F)
  #ordi3 = ordinate(phylobjetc, method="PCoA", distance="bray")
  
  ordi_1 = ordinate(phylobject_norm, method = "PCoA", distance="unifrac", weighted=T)
  #ordi2_1 = ordinate(phylobject_norm, method="PCoA", distance="unifrac", weighted=F)
  #ordi3_1 = ordinate(phylobject_norm, method="PCoA", distance="bray")
  
  #plot_ordination(phylobjetc, ordi, color="Treatment", shape="Origin") + theme_bw() -> Weighted_unifract_PCoA
  #plot_ordination(phylobjetc, ordi2, color="Treatment", shape="Origin") + theme_bw() -> Unweighted_unifract_PCoA
  #plot_ordination(phylobjetc, ordi3, color="Treatment", shape="Origin") + theme_bw() -> brai_PCoA
  
  plot_ordination(phylobject_norm, ordi_1, color="Treatment", shape="Origin") + theme_bw() -> Weighted_unifract_PCoA_rare
  #plot_ordination(phylobject_norm, ordi2_1, color="Treatment", shape="Origin") + theme_bw() -> Unweighted_unifract_PCoA_rare
  #plot_ordination(phylobject_norm, ordi3_1, color="Treatment", shape="Origin") + theme_bw() -> brai_PCoA_rare
  
  meta <- as(sample_data(phylobject_norm), "data.frame")
  adonis(distance(phylobject_norm, method="unifrac", weighted=T) ~ Origin * Treatment,
         data = meta)
  #adonis(distance(phylobject_norm, method="unifrac", weighted=F) ~ Origin * Treatment,
  #       data = meta)
  
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
  Counts_transformed %>% mutate(ID = SampleID) %>% select(ID,colnames(Counts_transformed)[1:(length(colnames(Counts_transformed))-1)]) -> Counts_transformed
  return(Counts_transformed)
}
Test_groups = function(DF, TEST="Nonparam"){
  #Statistical test
  if (TEST == "Nonparam"){
      Model = dunn.test(x=DF$Dependent, g= DF$Group)
      return(Model$P)
  } else if (TEST == "Normal"){
      Model = lm(Dependent ~ factor(Treatment) + Origin + Reads, data=DF)
      Model = summary(Model)
      return(Model$coefficients[14:15])
      
  }
  
  
}
Filter_abundance = function(M, threshold=0.2){
  #Filter species not seen in at least treshold*100% of samples
  apply(select(M, -ID), MARGIN=2, FUN = function(x){ (sum(x != 0)/length(x)) >= 0.2 }) -> Filter_vector
  F_names = colnames(select(M, -ID)[!Filter_vector])
  M %>% select(! F_names) -> M
  return(M)
}
Phylogenetic_distances_test = function(Distances, bacteria_to_sample_dic){
  library(qdapTools)
  lookup(terms= Distances$Var1, bacteria_to_sample_dic$Bacteria, bacteria_to_sample_dic$ID) ->nVar_1
  lookup(terms= Distances$Var2, bacteria_to_sample_dic$Bacteria, bacteria_to_sample_dic$ID) ->nVar_2
  lookup(terms= nVar_1, metadata$Sample, metadata$Group) -> nVar_12
  lookup(terms= nVar_2, metadata$Sample, metadata$Group) -> nVar_22
  
  Distances %>% mutate(Var1 = nVar_1, Var2=nVar_2, Group_1= nVar_12, Group_2=nVar_22) -> Distances 
  Distances %>% filter(! value == 0) -> Distances_compare
  
  paste(Distances_compare$Group_1, Distances_compare$Group_2, sep="_") -> comp1
  paste(Distances_compare$Group_2, Distances_compare$Group_1, sep="_") -> comp2
  Distances_compare %>% mutate(Comp1 = comp1, Comp2 = comp2, ID_n=seq(dim(Distances_compare)[1])) -> Distances_compare
  
  filter(Distances_compare, comp1=="AGEOCTRL_AGEOOME" | comp2=="AGEOCTRL_AGEOOME") %>% mutate(Distances_between= "Oral_ControlVSOral_Treat") -> Oral_ControlVSOral_Treat 
  filter(Distances_compare, comp1=="AGEFCTRL_AGEFOME" | comp2=="AGEFCTRL_AGEFOME") %>% mutate(Distances_between= "Fecal_ControlVSFecal_Treat") -> Fecal_ControlVSFecal_Treat 
  filter(Distances_compare, comp1=="AGEFCTRL_AGEOCTRL" | comp2=="AGEFCTRL_AGEOCTRL")%>% mutate(Distances_between= "Fecal_ControlVSOral_Control") -> Fecal_ControlVSOral_Control 
  filter(Distances_compare, comp1=="AGEFOME_AGEOCTRL" | comp2=="AGEFOME_AGEOCTRL") %>% mutate(Distances_between= "Fecal_TreatVSOral_Control") -> Fecal_TreatVSOral_Control
  filter(Distances_compare, comp1=="AGEFCTRL_AGEOOME" | comp2=="AGEFCTRL_AGEOOME") %>% mutate(Distances_between= "Fecal_ControlVSOral_Treat")-> Fecal_ControlVSOral_Treat
  filter(Distances_compare, comp1=="AGEFOME_AGEOOME" | comp2=="AGEFOME_AGEOOME") %>% mutate(Distances_between= "Fecal_TreatVSOral_Treat") -> Fecal_TreatVSOral_Treat
  
  filter(Distances_compare, comp1=="AGEOOME_AGEOOME") %>% mutate(Distances_between= "Oral_TreatVSOral_Treat") -> Oral_TreatVSOral_Treat
  filter(Distances_compare, comp1=="AGEOCTRL_AGEOCTRL") %>% mutate(Distances_between= "Oral_ControlVSOral_Control") -> Oral_ControlVSOral_Control
  filter(Distances_compare, comp1=="AGEFOME_AGEFOME")%>% mutate(Distances_between= "Fecal_TreatVSFecal_Treat") -> Fecal_TreatVSFecal_Treat
  filter(Distances_compare, comp1=="AGEFCTRL_AGEFCTRL") %>% mutate(Distances_between= "Fecal_ControlVSFecal_Control") -> Fecal_ControlVSFecal_Control
  
  Distances_compare2 = rbind(Oral_ControlVSOral_Treat,Fecal_ControlVSFecal_Treat,Fecal_ControlVSOral_Control,Fecal_TreatVSOral_Control,Fecal_ControlVSOral_Treat,Fecal_TreatVSOral_Treat,Oral_TreatVSOral_Treat,Oral_ControlVSOral_Control,Fecal_TreatVSFecal_Treat,Fecal_ControlVSFecal_Control)
  ggplot(Distances_compare2, aes(y=value,x=Distances_between))+ theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_boxplot() + geom_point() -> Plot_distances
  Model = dunn.test(x=Distances_compare2$value, g= Distances_compare2$Distances_between)
  return(list(Model, Plot_distances,Distances_compare2))
}
Make_tree_genus = function(Tree, Genus_name, All_samples, metadata){
  #Make phylogenetic tree for a given Name and statistical test on phylogenetic distances
  
  #Get tree for Genus of interest
  Filter_branches = Tree$tip.label[grepl(Genus_name, Tree$tip.label)]
  Subset_Tree = keep.tip(Tree, Filter_branches)
  
  #Get Samples of Genus of interest
  All_samples %>% select(c("ID",Filter_branches)) -> Subset_samples
  sapply(str_replace(Subset_samples$ID,"FDMP19H000760-1a_L1_", ""), FUN = function(x){paste(strsplit(x,"")[[1]][1:8],collapse="")}) -> New_names
  Subset_samples %>% mutate(ID=as.vector(New_names)) -> Subset_samples
  Subset_samples[2:(dim(Subset_samples)[2])] = as.integer(Subset_samples[2:(dim(Subset_samples)[2])] != 0) 
  Subset_samples %>% select(-ID) %>% t() %>% as_tibble() %>% mutate(label=colnames(select(Subset_samples,-ID))) %>% `colnames<-`(c(Subset_samples$ID, "label")) %>% select(c("label",Subset_samples$ID)) -> Subset_samples2
  #Add mouse Number to metadata
  metadata %>% mutate(Mouse = as.vector(sapply(metadata$Sample, FUN= function(X){paste(strsplit(X,"")[[1]][6:8],collapse="")}))) -> metadata
  #Getting distances among each bacteria
  cophenetic(Subset_Tree) -> Matrix_distances
  Matrix_distances[lower.tri(Matrix_distances)] <- 0
  
  melt(Matrix_distances) %>% as_tibble() -> Distances
  #Bacteria to sample connection
  Subset_samples %>% gather(Bacteria, Presence, 2:dim(Subset_samples)[2]) %>% filter(! Presence == 0) -> bacteria_to_sample_dic
  if (length(bacteria_to_sample_dic$Bacteria) == dim(bacteria_to_sample_dic)[1]){
    #If each bacteria is unique per sample, compute non-parametric test between distances
    OUT_f = Phylogenetic_distances_test(Distances, bacteria_to_sample_dic)
    Model = OUT_f[[1]]
    Plot_distances = OUT_f[[2]]
    Plot_distances + ggtitle(paste(c("Phylogenetic distances for ", Genus_name),collapse="")) -> Plot_distances
    print(Plot_distances)
    Distances = OUT_f[[3]]
    
    Results = tibble(Comparison = Model$comparisons, Pvalue=Model$P)
    #Differences between distances between oral and Fecal before and after treatment.
    Results %>% filter(Comparison == "Fecal_TreatVSOral_Treat - Fecal_ControlVSOral_Control" | Comparison =="Fecal_ControlVSOral_Control - Fecal_TreatVSOral_Treat" | Comparison == "Fecal_TreatVSOral_Control - Fecal_ControlVSOral_Control" | Comparison =="Fecal_ControlVSOral_Control - Fecal_TreatVSOral_Control") -> Comparison_interest
    Interest = c("Fecal_ControlVSOral_Control", "Fecal_TreatVSOral_Control", "", "Fecal_TreatVSOral_Treat")
    Distances %>% filter(Distances_between %in% Interest) -> D_interest
    D_interest %>% group_by(Distances_between) %>% summarise(mean(value)) %>% print()
    
    print(Comparison_interest)
  }
  
  #Prepare Phyloseq object for plotting
  Subset_samples2 %>% as.data.frame() %>% column_to_rownames("label") -> ASV_table
  taxmat = matrix(rownames(ASV_table))
  rownames(taxmat) = rownames(ASV_table)
  sampledata = metadata %>% mutate(Treatment = as.factor(Treatment)) %>% as.data.frame() %>% column_to_rownames("Sample")
  
  otu_table(ASV_table, taxa_are_rows = TRUE) -> ASV_table
  tax_table(taxmat) -> taxmat
  phyloseq::phyloseq(ASV_table, taxmat) -> phylobjetc
  merge_phyloseq(phylobjetc,Subset_Tree, sample_data(sampledata)) -> phylobjetc
  #taxa_names
  plot_tree(phylobjetc, color="Origin", shape="Treatment",label.tips="Mouse", ladderize="left", plot.margin=0.3) -> TREE_PLOT
  ggsave(filename = paste("TREE_", Genus_name,".pdf" , collapse="") ,plot = TREE_PLOT, scale = 3)
  #Subset_samples2 %>% gather(Sample, Presence, 2:dim(Subset_samples2)[2]) -> Subset_samples2
  #p = ggtree(Subset_Tree, layout="fan",branch.length="none") + geom_tiplab(size=0.5)
}
DE_seq_preparation = function(M, metadata){
  #Get data into DEseq2
  M %>% as.data.frame() %>% column_to_rownames("ID") -> ASV_table
  phyloseq::otu_table(ASV_table, taxa_are_rows = F) -> ASV_table
  
  taxmat = matrix(colnames(ASV_table))
  rownames(taxmat) = colnames(ASV_table)
  tax_table(taxmat) -> taxmat
  
  sampledata = metadata %>% mutate(Treatment = as.factor(Treatment)) %>% as.data.frame() %>% column_to_rownames("Sample")
  sample_data(sampledata) -> sampledata
  
  phyloseq::phyloseq(ASV_table, taxmat, sampledata) -> phylobjetc
  
  phyloseq_to_deseq2(phylobjetc, design= ~ Group) -> dds
  return(dds)
  
}
DE_seq_normalization = function(dds){
  #Normalize data
  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized=TRUE)
  
  as.data.frame(normalized_counts)  %>%
    t() %>% as_tibble() %>% mutate(ID=colnames(normalized_counts)) %>%
    `colnames<-`(c(rownames(normalized_counts), "ID")) -> Abundance_DEseq_norm
  return(Abundance_DEseq_norm)
}
test_with_ancom = function(Abundance, metadata){
  select(Abundance, -ID) %>% t() %>% `colnames<-`(Abundance$ID) %>% as.data.frame() -> feature_table
  Final_results = tibble()
  for (N in seq(length(unique(metadata$Group)))){
    GROUP1 = unique(metadata$Group)[N]
    if (N == length(unique(metadata$Group))){ next }
    for (M in seq((N+1),length(unique(metadata$Group)))){
      if (M > length(unique(metadata$Group))){ next }
      print(c(N,M))
      GROUP2 = unique(metadata$Group)[M]
      filter(metadata, Group %in% c(GROUP1, GROUP2)) -> meta_group
      feature_table %>% select(meta_group$Sample) -> group_feature_table
      
      feature_table_pre_process(group_feature_table, meta_group, "Sample", group_var = "Group", out_cut = 0.05, zero_cut = 0.90, 0, T) -> list_outcome
      ANCOM(list_outcome[[1]], list_outcome[[2]], list_outcome[[3]], "Group", "fdr", 0.05) -> list_model
    
      list_model[[1]] %>% as_tibble() %>% filter(!W == Inf) %>% mutate(Groups=paste(c(GROUP1, GROUP2), collapse="-")) -> R1
      Final_results = rbind(Final_results, R1)
    }
  }
  return(Final_results)
}

Make_tree_genus2 = function(Tree, Genus_name, All_samples, metadata, Plot = F){
  #Make phylogenetic tree for a given Name and statistical test on phylogenetic distances
  
  if (Genus_name != "ALL"){
    #Get tree for Genus of interest
    Filter_branches = Tree$tip.label[grepl(Genus_name, Tree$tip.label)]
    Subset_Tree = keep.tip(Tree, Filter_branches)
    
    #Get Samples of Genus of interest
    All_samples %>% select(c("ID",Filter_branches)) -> Subset_samples
  } else{  Subset_Tree = Tree ; Subset_samples = All_samples }
  if (length(Subset_Tree$tip.label) < 3){ return(c(NA, NA, NA)) }
  sapply(str_replace(Subset_samples$ID,"FDMP19H000760-1a_L1_", ""), FUN = function(x){paste(strsplit(x,"")[[1]][1:8],collapse="")}) -> New_names
  Subset_samples %>% mutate(ID=as.vector(New_names)) -> Subset_samples
  Subset_samples[2:(dim(Subset_samples)[2])] = as.integer(Subset_samples[2:(dim(Subset_samples)[2])] != 0) 
  Subset_samples %>% select(-ID) %>% t() %>% as_tibble() %>% mutate(label=colnames(select(Subset_samples,-ID))) %>% `colnames<-`(c(Subset_samples$ID, "label")) %>% select(c("label",Subset_samples$ID)) -> Subset_samples2
  #Add mouse Number to metadata
  metadata %>% mutate(Mouse = as.vector(sapply(metadata$Sample, FUN= function(X){paste(strsplit(X,"")[[1]][6:8],collapse="")}))) -> metadata
  arrange(metadata, Sample) -> metadata
  arrange(Subset_samples, ID) %>% mutate(Group = metadata$Group, Mouse=metadata$Mouse, Origin=metadata$Origin, Treatment=metadata$Treatment) -> Data_merg
  Controls_same = vector() ; Controls_other= vector() ; Treated_same = vector() ; Treated_other = vector()
  for (MOUSE in unique(metadata$Mouse)){
    Data_merg %>% filter(Mouse == MOUSE) -> Data_merg2
    Data_merg2 %>% select(-c("Group", "Mouse", "Origin", "Treatment")) -> ID_x
    equals_two = function(x){x==2}
    equals_one = function(x){x==1}
    ID_x %>% summarise_if(is.numeric,sum) %>% as_vector() -> Tip

    sum(Tip == 1 ) -> Other_tip
    sum(Tip == 2 ) -> Same_tip
    
    if (unique(Data_merg2$Treatment) == 0){ 
      Controls_same = c(Controls_same, Same_tip) ; Controls_other = c(Controls_other, Other_tip)  
    }else{ Treated_same = c(Treated_same, Same_tip) ; Treated_other = c(Treated_other, Other_tip)   }
  }
  data.frame(Group_control = c(sum(Controls_same), sum(Controls_other)), Group_treated= c(sum(Treated_same), sum(Treated_other))) -> Table_transfer
  #data.frame(Group_control = c(median(Controls_same), median(Controls_other)), Group_treated= c(median(Treated_same), median(Treated_other))) -> Table_transfer
  rownames(Table_transfer) = c("Same_tip", "Other_tip")
  Test = chisq.test(Table_transfer)
  Summary_results = c(Test$p.value, sum(Controls_same)/sum(Controls_other), sum(Treated_same)/sum(Treated_other))
  
  
  if (Plot != F){
    #Prepare Phyloseq object for plotting
    Subset_samples2 %>% as.data.frame() %>% column_to_rownames("label") -> ASV_table
    taxmat = matrix(rownames(ASV_table))
    rownames(taxmat) = rownames(ASV_table)
    sampledata = metadata %>% mutate(Treatment = as.factor(Treatment)) %>% as.data.frame() %>% column_to_rownames("Sample")
    
    otu_table(ASV_table, taxa_are_rows = TRUE) -> ASV_table
    tax_table(taxmat) -> taxmat
    phyloseq::phyloseq(ASV_table, taxmat) -> phylobjetc
    merge_phyloseq(phylobjetc,Subset_Tree, sample_data(sampledata)) -> phylobjetc
    #taxa_names
    #plot_tree(phylobjetc,  use.edge.length = F, color="Origin", shape="Treatment", ladderize="left", plot.margin=0.3) -> TREE_PLOT
    New_tree(phylobjetc, color="Group", ladderize=F, plot.margin=0.3, LABEL="Mouse", text.size = 3, nodelab= nodeplotblank, base.spacing=0.1) -> TREE_PLOT
    str_replace(Genus_name, "/", "_OR_") -> Genus_name
    trim(Genus_name) -> Genus_name
    FILE = paste(c("Figures/TREE_", Genus_name,".pdf") , collapse="")
    print(FILE)
    print(TREE_PLOT)
    ggsave(filename =  FILE, plot = TREE_PLOT, scale = 3)
    
  }
  return(Summary_results)
}
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
New_tree = function(physeq, method = "sampledodge", nodelabf = NULL, 
                    color = NULL, shape = NULL, size = NULL, min.abundance = Inf, 
                    label.tips = NULL, text.size = NULL, sizebase = 5, base.spacing = 0.02, 
                    ladderize = FALSE, plot.margin = 0.2, title = NULL, treetheme = NULL, 
                    justify = "jagged", LABEL =NULL) {
  library(data.table)
  fix_reserved_vars = function(aesvar) {
    aesvar <- gsub("^abundance[s]{0,}$", "Abundance", 
                   aesvar, ignore.case = TRUE)
    aesvar <- gsub("^OTU[s]{0,}$", "OTU", aesvar, 
                   ignore.case = TRUE)
    aesvar <- gsub("^taxa_name[s]{0,}$", "OTU", 
                   aesvar, ignore.case = TRUE)
    aesvar <- gsub("^sample[s]{0,}$", "Sample", 
                   aesvar, ignore.case = TRUE)
    return(aesvar)
  }
  if (!is.null(label.tips)) {
    label.tips <- fix_reserved_vars(label.tips)
  }
  if (!is.null(color)) {
    color <- fix_reserved_vars(color)
  }
  if (!is.null(shape)) {
    shape <- fix_reserved_vars(shape)
  }
  if (!is.null(size)) {
    size <- fix_reserved_vars(size)
  }
  if (is.null(phy_tree(physeq, FALSE))) {
    stop("There is no phylogenetic tree in the object you have provided.\n", 
         "Try phy_tree(physeq) to see for yourself.")
  }
  if (!inherits(physeq, "phyloseq")) {
    method <- "treeonly"
  }
  treeSegs <- tree_layout(phy_tree(physeq), ladderize = ladderize)
  edgeMap = aes(x = xleft, xend = xright, y = y, yend = y)
  vertMap = aes(x = x, xend = x, y = vmin, yend = vmax)
  p = ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) + 
    geom_segment(vertMap, data = treeSegs$vertDT)
  if (is.null(text.size)) {
    text.size <- manytextsize(ntaxa(physeq))
  }
  if (!is.null(label.tips) & method != "sampledodge") {
    labelDT = treeSegs$edgeDT[!is.na(OTU), ]
    if (!is.null(tax_table(object = physeq, errorIfNULL = FALSE))) {
      taxDT = data.table(tax_table(physeq), OTU = taxa_names(physeq), 
                         key = "OTU")
      labelDT = merge(x = labelDT, y = taxDT, by = "OTU")
    }
    if (justify == "jagged") {
      labelMap <- aes_string(x = "xright", y = "y", 
                             label = label.tips, color = color)
    }
    else {
      labelMap <- aes_string(x = "max(xright, na.rm=TRUE)", 
                             y = "y", label = label.tips, color = color)
    }
    p <- p + geom_text(labelMap, data = labelDT, size = I(text.size), 
                       hjust = -0.1, na.rm = TRUE)
  }
  if (is.null(nodelabf)) {
    nodelabf = howtolabnodes(physeq)
  }
  p = nodelabf(p, treeSegs$edgeDT[!is.na(label), ])
  p = nodelabf(p, treeSegs$vertDT[!is.na(label), ])
  if (is.null(treetheme)) {
    treetheme <- theme(axis.ticks = element_blank(), axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), axis.title.y = element_blank(), 
                       axis.text.y = element_blank(), panel.background = element_blank(), 
                       panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  }
  if (inherits(treetheme, "theme")) {
    p <- p + treetheme
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (method != "sampledodge") {
    return(p)
  }
  dodgeDT = treeSegs$edgeDT[!is.na(OTU), ]
  dodgeDT = merge(x = dodgeDT, y = data.table(psmelt(physeq), 
                                              key = "OTU"), by = "OTU")
  if (justify == "jagged") {
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  if (!is.null(color) | !is.null(shape) | !is.null(size)) {
    setkeyv(dodgeDT, cols = c("OTU", color, shape, 
                              size, LABEL))
  }
  else {
    setkey(dodgeDT, OTU, Sample)
  }
  dodgeDT[, `:=`(h.adj.index, 1:length(xright)), by = OTU]
  if (justify == "jagged") {
    dodgeDT[, `:=`(xdodge, (xright + h.adj.index * 
                              base.spacing * max(xright, na.rm = TRUE)))]
  }
  else {
    dodgeDT[, `:=`(xdodge, max(xright, na.rm = TRUE) + 
                     h.adj.index * base.spacing * max(xright, na.rm = TRUE))]
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  dodgeMap <- aes_string(x = "xdodge", y = "y", 
                         color = color, fill = color, shape = shape, label = LABEL)
  p <- p + geom_text(dodgeMap, data = dodgeDT, na.rm = TRUE, size= text.size)
  if (!is.null(size)) {
    p <- p + scale_size_continuous(trans = log_trans(sizebase))
  }
  if (any(dodgeDT$Abundance >= min.abundance[1])) {
    pointlabdf = dodgeDT[Abundance >= min.abundance[1], ]
    p <- p + geom_text(mapping = aes(xdodge, y, label = Abundance), 
                       data = pointlabdf, size = text.size, na.rm = TRUE)
  }
  if (!is.null(label.tips)) {
    tiplabDT = dodgeDT
    tiplabDT[, `:=`(xfartiplab, max(xdodge)), by = OTU]
    tiplabDT <- tiplabDT[h.adj.index == 1, .SD, by = OTU]
    if (!is.null(color)) {
      if (color %in% sample_variables(physeq, errorIfNULL = FALSE)) {
        color <- NULL
      }
    }
    labelMap <- NULL
    if (justify == "jagged") {
      labelMap <- aes_string(x = "xfartiplab", y = "y", 
                             label = label.tips, color = color)
    }
    else {
      labelMap <- aes_string(x = "max(xfartiplab, na.rm=TRUE)", 
                             y = "y", label = label.tips, color = color)
    }
    p <- p + geom_text(labelMap, tiplabDT, size = I(text.size), 
                       hjust = -0.1, na.rm = TRUE)
  }
  min.x <- -0.01
  max.x <- dodgeDT[, max(xright, na.rm = TRUE)]
  if ("xdodge" %in% names(dodgeDT)) {
    max.x <- dodgeDT[, max(xright, xdodge, na.rm = TRUE)]
  }
  if (plot.margin > 0) {
    max.x <- max.x * (1 + plot.margin)
  }
  p <- p + scale_x_continuous(limits = c(min.x, max.x))
  return(p)
}


#####################
#####################

Genus = colnames(Abundance)[grepl(".Genus",colnames(Abundance))]
Abundance %>% select(c("ID", Genus)) -> Abundance

Total_number = Preliminary_checks(Abundance, metadata)
metadata %>% mutate(Reads = Total_number) -> metadata

Comparisons = c("Fecal_Control-Fecal_Treat", "Fecal_Control-Oral_Control", "Fecal_Treat-Oral_Control", "Fecal_Control-Oral_Treat",  "Fecal_Treat-Oral_Treat", "Oral_Control-Oral_Treat")




Analysis_No_norm = function(){
  #No normalization
  Alpha_diversity(Abundance, metadata)
  Beta_diversity(Abundance, metadata)
  
  Filter_abundance(Abundance) -> Abundance_f
  apply(select(Abundance_f, -ID), MARGIN = 2, FUN = function(x){ Test_groups(mutate(metadata, Dependent=x))}) -> Model_outcome
  Model_outcome %>% as_tibble() %>% mutate(Test=Comparisons) -> Model_outcome
  Model_outcome %>% t() %>%as_tibble() %>% mutate(Bug=colnames(Model_outcome)) %>%   
    `colnames<-`(c(Comparisons, "Bug")) %>% gather(Test, Pvalue, 1:length(Comparisons)) %>%
    mutate(FDR = p.adjust(Pvalue, "fdr")) %>% filter(!is.na(FDR)) -> Model_outcome
  
  write_csv(Model_outcome, "Results/Model_nonorm.csv")
}
Analysis_clr = function(){
  #clr normalization
  Abundance_clr = Transformation_composition(Abundance)
  #Alpha_diversity(Abundance, metadata)
  Beta_diversity(Abundance_clr, metadata, Input="CLR")
  Filter_abundance(Abundance_clr) -> Abundance_clr
  
  apply(select(Abundance_clr, -ID), MARGIN = 2, FUN = function(x){ Test_groups(mutate(metadata, Dependent=x),TEST="Normal")}) -> Model_outcome_clr
  Model_outcome_clr %>% t() %>% as_tibble() %>% mutate(Bug=colnames(Model_outcome_clr)) %>% 
    `colnames<-`(c("Treatment", "Origin", "Bug")) %>% gather(Test, Pvalue, 1:2) %>%
    mutate(FDR = p.adjust(Pvalue, "fdr")) %>% filter(!is.na(FDR)) -> Model_outcome_clr
  write_csv(Model_outcome_clr, "Results/Model_clr.csv")
}
#Rarefraction

##Check difference in number of species before and after rarefy
S <- specnumber(select(Abundance, -ID)) # observed number of species
Srare <- rarefy(x=select(Abundance,-ID),MARGIN=1, sample=min(Total_number))
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

rarecurve(select(Abundance, -ID), step = 10, min(Total_number), xlab = "Sample Size", ylab = "Species",
          label = TRUE) -> rare_curve
ggsave(rare_curve, "Figures/Rarecurve.pdf")

##Perform rarefy
vegan::rrarefy(x=select(Abundance,-ID), sample=min(Total_number)) %>% as_tibble() %>% mutate(ID = Abundance$ID) -> Abundance_rarefy
#Srare <- specnumber(select(Abundance_rarefy, -ID))
Total_number = Preliminary_checks(Abundance_rarefy, metadata)
metadata %>% mutate(Reads = Total_number) -> metadata_rarefy

Alpha_diversity(Abundance_rarefy, metadata_rarefy)


Filter_abundance(Abundance_rarefy) -> Abundance_rarefy
Abundance_rarefy %>% select(-`Unknown|Genus`) -> Abundance_rarefy

Beta_diversity(Abundance_rarefy, metadata_rarefy)



test_with_ancom(Abundance_rarefy, metadata_rarefy) -> ancom_results
write_csv(ancom_results, "Results/Model_rarefy_ancom.csv")
ancom_results  %>% arrange(desc(W)) %>% ggplot(aes(x=Groups, y=W)) + geom_boxplot(outlier.shape = NA
)+ geom_jitter(aes(col=detected_0.7==T)) +theme_bw()+ coord_flip() #theme(axis.text.x = element_text(angle = 90, hjust = 1))



apply(select(Abundance_rarefy, -ID), MARGIN = 2, FUN = function(x){ print(x);Test_groups(mutate(metadata_rarefy, Dependent=x))}) -> Model_outcome_rarefy
Model_outcome_rarefy %>% t() %>%as_tibble() %>% mutate(Bug=colnames(Model_outcome_rarefy)) %>%   
  `colnames<-`(c(Comparisons, "Bug")) %>% gather(Test, Pvalue, 1:length(Comparisons)) %>%
  mutate(FDR = p.adjust(Pvalue, "fdr")) %>% filter(!is.na(FDR)) %>% arrange(FDR) -> Model_outcome_rarefy
write_csv(Model_outcome_rarefy, "Results/Model_rarefy.csv")


Model_outcome_rarefy  %>% ggplot(aes(x=Test, y=-log10(Pvalue))) + geom_boxplot(outlier.shape = NA
)+ geom_sina(aes(col=FDR<0.05)) +theme_bw()+ coord_flip() #theme(axis.text.x = element_text(angle = 90, hjust = 1))

#AGEFCTRL-AGEOCTRL

ancom_results %>% filter(detected_0.7 == T) -> DE_ancom
for (id in unique(DE_ancom$taxa_id)){
  Model_outcome_rarefy %>% filter(FDR < 0.05) %>% filter(Bug == id) -> Matching
  print(c(id, dim(Matching)[1], dim(filter(DE_ancom, taxa_id==id))[1]))
  
}


Analysis_DEseq = funcion(){
  #####DEseq normalization and differential expression
  
  library(DESeq2)
  
  dds = DE_seq_preparation(Abundance,metadata)
  Abundance_DEseq_norm =  DE_seq_normalization(dds) 
  
  Alpha_diversity(M = Abundance_DEseq_norm,metadata = metadata)
  Beta_diversity(M = Abundance_DEseq_norm,metadata = metadata)
  
  #Modelling
  Model_DE = DESeq(dds, test="Wald", fitType="local") 
  plotDispEsts(Model_DE) #Check fit
  
  ###Get results doing different comparisons
  res = results(Model_DE, cooksCutoff = FALSE, contrast = c("Group", "AGEFCTRL", "AGEFOME"))
  as.data.frame(res) %>% rownames_to_column("Test") %>% as_tibble() %>% arrange(padj) %>% mutate(test="FecalControl_vsFecalTreat") -> res_1
  plotMA(res, ylim=c(-2,2))
  
  
  res = results(Model_DE, cooksCutoff = FALSE, contrast = c("Group", "AGEOCTRL", "AGEOOME"))
  as.data.frame(res) %>% rownames_to_column("Test") %>% as_tibble() %>% arrange(padj) %>% mutate(test="OralControl_vsOralTreat") -> res_2
  plotMA(res, ylim=c(-2,2))
  
  res = results(Model_DE, cooksCutoff = FALSE, contrast = c("Group", "AGEOCTRL", "AGEFCTRL"))
  as.data.frame(res) %>% rownames_to_column("Test") %>% as_tibble() %>% arrange(padj) %>% mutate(test="OralControl_vsFecalControl") -> res_3
  plotMA(res, ylim=c(-2,2))
  
  res = results(Model_DE, cooksCutoff = FALSE, contrast = c("Group", "AGEFOME", "AGEOOME"))
  as.data.frame(res) %>% rownames_to_column("Test") %>% as_tibble() %>% arrange(padj) %>% mutate(test="OralTreat_vsFecalTreat") -> res_4
  plotMA(res, ylim=c(-2,2))
  
  rbind(res_1,res_2, res_3, res_4) -> negative_binomial_results
  write_csv(negative_binomial_results, "Results/Model_deseq.csv")
}

##Tree

Result_tranfer_total = tibble()
Abundance_rarefy %>% select(-"Unknown|Genus") -> Abundance_rarefy
for (Genus in colnames(Abundance_rarefy)){
  print(Genus)
  if (Genus == "ID"){ next }
  S = strsplit(Genus, split = "\\|")
  Genus = S[[1]][1]
  Results_transfer =  Make_tree_genus2(Tree, Genus, All_samples, metadata)
  data.frame(t(Results_transfer)) -> Results_transfer
  colnames(Results_transfer) = c("Pvalue", "Proportion_identical_Control", "Proportion_identical_Treated")
  Results_transfer %>% mutate(Genus_name = Genus, Diff= Proportion_identical_Control - Proportion_identical_Treated) -> Results_transfer
  Result_tranfer_total = rbind(Result_tranfer_total,Results_transfer)
}
#Results_transfer =  Make_tree_genus2(Tree, "ALL", All_samples, metadata)

as_tibble(Result_tranfer_total) %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(FDR) -> Result_tranfer_total
write_tsv(Result_tranfer_total,"Results/Transfers.tsv")

Result_tranfer_total %>% filter(Pvalue < 0.05) -> Transfered
sapply(Transfered$Genus_name, FUN = function(x){ Make_tree_genus2(Tree, x, All_samples, metadata, T)} )


Result_tranfer_total %>% filter(Genus_name == "Corynebacterium")
