#Analysis of CMV in LLD

###Packages
library(tidyverse)
library(ape)
library(vegan)
library(viridis)
library(patchwork)
library(car) #cehck inflattion coefficients with vif
###### Set up
setwd("~/Resilio Sync/Transfer/PhD/Immuno_Senescence/")


### Functions
Find_positive_probe = function(Person_ID, CMV_raw,CMV_info){
  CMV_raw %>% filter(ID == Person_ID) %>% t() %>% as.data.frame() %>% rownames_to_column("Probe") %>% as_tibble() %>% filter(V1 == 1) -> Found
  filter(CMV_info, Probe %in% Found$Probe)
}
Clean_MP3 = function(Microbial){
 unknown_rel = Microbial$UNKNOWN
 Microbial %>% select(-UNKNOWN) -> M2
 lapply(colnames(M2)[2:dim(M2)[2]], FUN= function(x){ str_split(x,"\\|")[[1]] -> y ;  return(y[length(y)]) }) %>% as_vector() -> New_names
 colnames(M2) = c("ID", New_names)
#to go back to 100% without unknown we need to apply Back_to_100 in each taxonomic level 
 lapply(str_split(colnames(select(M2, -ID)), "__"), FUN=function(x){ x[1] }) %>% as_vector() ->Taxonomies ; Taxonomies %>% unique() -> Taxonomy_Levels
 Bugs = colnames(select(M2, -ID)) ; tibble(Bug = Bugs, Level = Taxonomies) -> Bug_array
 M3 = tibble(ID = M2$ID)
 for (T_L in Taxonomy_Levels){
   filter(Bug_array, Level == T_L) -> Bugs
   Back_to_100(select(M2, Bugs$Bug)) -> Normalized
   M3 = left_join(M3, mutate(Normalized, ID = M2$ID ),  by= "ID" )
 }
 
 return(list(unknown_rel, M2, M3, Bug_array))
 
}
Back_to_100 = function(Microbial){
  #We can go back to 100%, but this function should be applied per taxonomic level
  apply(Microbial,1, FUN= function(x){ sum(x) -> y ;  return(y) }) %>% as_vector() -> Total_per_sample
  Microbial/Total_per_sample -> Normalized
  return(as_tibble(Normalized))
}
Transformation_clr = function(Count_table){
  SampleID = Count_table$ID
  Count_table %>% select(-ID) -> Counts
  #Transform counts
  
  colnames(Counts) -> N
  Counts %>% t() -> Check
  Counts_transformed = as_tibble(t(microbiome::abundances(x=as.data.frame(Check), transform="clr")))
  ###
  Counts_transformed %>% mutate(ID = SampleID) -> Counts_transformed
  return(Counts_transformed)
}
R2_calc = function(Real,pred){
  rss = sum((pred - Real)^2)
  tss = sum((Real - mean(Real))^2)
  rsq = 1 - rss/tss
  return(rsq)
}
Impute_telomere_data = function(Telomere_markers){
  set.seed(100)
  Tels = c("Lymphocytes","Granulocytes","`Naive T-cells`","`Memory T-cells`", "`B-cells`","`NK-cells`" )#colnames(Telomere_markers)[2: dim(Telomere_markers)[2]]
  Telomere_markers %>% drop_na() -> Complete_telomere
  N = dim(Complete_telomere)[1]
  sample(x= N, replace=F, size= N*0.2) -> Test_samples
  Complete_telomere[Test_samples, 1:(length(Tels)+1)] -> TEST
  Complete_telomere %>% filter(!ID %in% TEST$ID) -> TRAIN
  Model_list = list() ; n =1 ; Acc_tibble = tibble()
  for (i in Tels){
    Formula = as.formula( paste( c( i, "~", paste(Tels[! Tels == i], collapse="+")), collapse= " " )  )
    lm(Formula, TRAIN) -> Model
    predict.lm(Model, TEST) -> Predicted
    Real = as_vector(TEST[colnames(TEST) == str_replace_all(i,"`","") ])
    R2_calc(Real, Predicted) -> pred_r2
    rbind(Acc_tibble, tibble(Telomere= i, R2_test = pred_r2 )) -> Acc_tibble
    Model_list[[n]] = Model
    n = n+1
    
    Predicted = predict.lm( Model , Telomere_markers[is.na(Telomere_markers[colnames(Telomere_markers) == str_replace_all(i,"`","") ]), 1:dim(Telomere_markers)[2]]   )
    Telomere_markers[is.na(Telomere_markers[colnames(Telomere_markers) == str_replace_all(i,"`","") ]),  colnames(Telomere_markers) == str_replace_all(i,"`","")]   = as.numeric(Predicted)
  }
  
  return(list(Telomere_markers, Acc_tibble))
  
}
#####Read_data

#Covariates
Phenotypes <- read_delim("~/Resilio Sync/LLD phenotypes for Sergio/New Phenotypes clean /Merged_Phenotypes_LLD.csv", delim="|")
Phenotypes %>% dplyr::select(ID,Age, Sex) -> Covariates

#CMV data
CMV_raw = read_csv("CMV/Data/CMV_data.csv") ; CMV_info = read_csv("CMV/Data/CMV_probes_info.csv")
apply(select(CMV_raw, -c(X1,ID) ), 1, sum) %>% as.vector() -> CMV_presence
CMV_raw %>% select(ID) %>% mutate(CMV_presence) -> CMV
CMV %>% ggplot() + geom_histogram(aes(x=CMV_presence)) + theme_bw()

Immuno_clusters = read_tsv("CMV/Data/PCs.csv")
Covariates_ab  = read_tsv("CMV/Data/Antibody_Covariates.tsv")
colnames(Immuno_clusters)[4] = "sample_id" 
Immuno_clusters$sample_id = sapply(Immuno_clusters$sample_id, FUN=function(x){ str_split(x, "_")[[1]][2] } )
left_join(Immuno_clusters, Covariates_ab) %>% select(ID, Cluster) -> Immuno_clusters
#There are a bunch of people with few antigens... Check which are those antigens
CMV %>% filter(! CMV_presence == 0) %>% arrange(CMV_presence) %>% print(n=20)
Find_positive_probe("LLDeep_0487",  CMV_raw,CMV_info)

#Bioaging and Telomeric markers

Age_markers <- read_tsv("~/Resilio Sync/LLD phenotypes for Sergio/bio_aging/Combined_data_uncorrected.tsv")
Telomere_markers <- colnames(Age_markers)[grepl("MTL",colnames(Age_markers))]
#ATTENTION: In telomere work I used drop_na, so that all telomeres are completed, here I am not
Telomere_markers = Age_markers %>% dplyr::select(c("Sample",Telomere_markers)) %>% mutate(MTL_gran = as.numeric(MTL_gran), `MTL_CD20+`=as.numeric(`MTL_CD20+`), `MTL_CD45+_CD20-`=as.numeric(`MTL_CD45+_CD20-`),`MTL_CD57+`=as.numeric(`MTL_CD57+`))
Other_age_markers = dplyr::select(Age_markers, c(Sample, Methyl.Age.Hannum., dCT)) %>% mutate(dCT = as.numeric(dCT)) 
colnames(Telomere_markers) = c("ID", "Lymphocytes", "Granulocytes", "Naive T-cells", "Memory T-cells","B-cells","NK-cells")
Telomere_markers[rowSums(is.na(Telomere_markers)) != (ncol(Telomere_markers)-1), ] -> Telomere_markers
Other_age_markers[rowSums(is.na(Other_age_markers)) != (ncol(Other_age_markers)-1), ] -> Other_age_markers
colnames(Other_age_markers)[1] = "ID"
#Microbiome abundance table
###MISSING FOLLOWUP TABLE (Is there in MP3?)
Microbial = read_tsv("/Users/sergio/Resilio Sync/Transfer/PhD/TMAO_project/MSS/Data/MetaPhlan3_table.tsv") %>% filter(grepl("LLD",ID))
List_microbial = Clean_MP3(Microbial) #This removes the unknown column and shorten the taxa name. Argument [1]  are the unknown fractions, the second is the relative abundances without re-assessment, the third are relative abundances after re-assessment, [4] are bugs an taxonomy, used for going for specific taxonomic levels
Unknown_fraction = List_microbial[[1]] ; Microbial = List_microbial[[2]] ; Microbial2 = List_microbial[[3]] ; Bug_level = List_microbial[[4]]
#CLR normalization of microbial abundance, uses the package Microbiome. We perform the Clr normalization per Taxonomic level, and then concatenate the columns together
M4 = tibble(ID = Microbial2$ID) 
for (L in unique(Bug_level$Level) ){
  Bugs = filter(Bug_level, Level == L)
  Microbial2 %>% select(c("ID", Bugs$Bug)) -> Sub_microbial
  Transformation_clr(Sub_microbial) -> Sub_microbial
  left_join(M4, Sub_microbial, by="ID") -> M4
}



#Immuno phenotypes

all_cell_types = read_tsv("~/Resilio Sync/Transfer/PhD/Telomeres_project/BIOS_cell_types_DeconCell_2019-03-08.LLD_subset.txt") #Cell counts from Blood bulkRNA deconvolution
all_cell_types %>% gather(key = var_name, value = value, 2:ncol(all_cell_types)) %>% spread_(key = names(all_cell_types)[1],value = 'value') %>% mutate(ID = var_name) %>% dplyr::select(-var_name) -> all_cell_types


#Inflammation and CVD



########################
####Analysis############
########################

CMV %>% mutate(CMV_presence = ifelse(CMV_presence > 0, 1, 0 )  ) -> CMV_presence
Immuno_clusters %>% mutate(Cluster = abs(Cluster-1)) -> Immuno_clusters
left_join(CMV_presence, select(Immuno_clusters, c(ID, Cluster))) %>% select(-CMV_presence) -> CMV_presence2

#Developer definition
clean_ID = str_replace(CMV$ID, "_F", "")
CMV_presence %>% mutate(Time = ifelse(grepl("_F", ID), "FA", "BL" ), ID=clean_ID )  -> CMV_dev
CMV_dev %>% group_by(ID) %>% summarise(N = n(), D = sum(CMV_presence) ) %>% filter(N>1) %>% filter(D == 1) -> Putative_dev
CMV_dev %>% filter(ID %in% Putative_dev$ID) %>% spread(Time, CMV_presence) %>% filter(FA == 1) -> Developers

CMV_presence2 %>% mutate(Time = ifelse(grepl("_F", ID), "FA", "BL" ), ID=clean_ID )  -> CMV_dev2
CMV_dev2 %>% group_by(ID) %>% summarise(N = n(), D = sum(Cluster) ) %>% filter(N>1) %>% filter(D == 1) -> Putative_dev2
CMV_dev2 %>% filter(ID %in% Putative_dev2$ID) %>% spread(Time, Cluster) %>% filter(FA == 1) -> Developers2

# (!) Interstingly, there are people who "loose" their antigens. Most of them only had 1 to start with. 
CMV %>% mutate(Time = ifelse(grepl("_F", ID), "FA", "BL" ), ID=clean_ID ) %>% filter(ID %in% Putative_dev$ID) %>% spread(Time, CMV_presence)
CMV_dev %>% mutate(Dev = ifelse(ID %in% Developers$ID, 1, 0)) -> CMV_dev

##Age vs prevalence
left_join(CMV_presence,Covariates) %>% ggplot() + geom_histogram(stat = "identity",aes(x=Age, y= CMV_presence )) +
  theme_bw() + geom_histogram(aes(x=Age), stat="count", fill="#00846b", alpha=0.2) + ylab("Number of participants")

##########TELOMERE LENGTH ANALYSIS #########################

##1. PCA of Telomeric length. Are PCs determined by CMV?
#PCA does not accept NAs, so I imputed them using a linear model. Imputed accuracy is the computed R2 based on test data (20% of the complete dataset)
List_telomere_imputation = Impute_telomere_data(Telomere_markers)
Imputed_Telomere_markers = List_telomere_imputation[[1]] ; Imputed_accuracy = List_telomere_imputation[[2]]
#Compute PCA and make two plots. 1- PCA showing age, 2-PCA showing CMV status
PCA <- rda(select(Imputed_Telomere_markers, -ID), scale = FALSE)
Variance_per_axis = as.vector(PCA$CA$eig)/sum(PCA$CA$eig) * 100
Standarized_weights = summary(PCA)$species
PCs = summary(PCA)$sites[,1:6]
left_join( mutate(as_tibble(PCs), ID = Imputed_Telomere_markers$ID) , Covariates ) -> PCs
left_join(PCs, CMV_presence) -> PCs
ggplot() + geom_point(data=PCs, aes(x=PC1, y=PC2, shape= as.factor(Sex) ,col = Age )) + theme_bw() + xlab(paste(c("PC1(", as.character(round(Variance_per_axis[1], 2)), "%)"), collapse="")) + ylab(paste(c("PC2(", as.character(round(Variance_per_axis[2], 2)), "%)"), collapse=""))   + scale_color_viridis() -> PCA_plot
ggplot() + geom_point(data=PCs, aes(x=PC1, y=PC2, shape= as.factor(Sex) ,col = as.factor(CMV_presence) )) + theme_bw() + xlab(paste(c("PC1(", as.character(round(Variance_per_axis[1], 2)), "%)"), collapse="")) + ylab(paste(c("PC2(", as.character(round(Variance_per_axis[2], 2)), "%)"), collapse=""))  + scale_color_viridis(discrete=TRUE)  -> PCA_plot2
PCA_plot + PCA_plot2
#Model PCs, which PCs are mainly driven by Age? Do we see PCs driven by CMV presence?
Results_PC_models = tibble()
#More PC1 means shorter telomeres! So, if positive value, PC2 mainly represents shorter  NK-cell, PC3 represents shorter memory Tcells
Standarized_weights
for (PC in colnames(PCs)[1:6] ){
  Model = as.formula(paste( c(PC, " ~ Age + Sex + CMV_presence"), collapse="" ))
  Model2 = as.formula(paste( c(PC, " ~ Age + Sex * CMV_presence"), collapse="" ))
  Model3 = as.formula(paste( c(PC, " ~ Age * CMV_presence + Sex"), collapse="" ))
  Model4 = as.formula(paste( c(PC, " ~ Age * CMV_presence * Sex"), collapse="" ))
  
  lm(Model, PCs) -> Model1
  lm(Model2, PCs) -> Model2
  lm(Model3, PCs) -> Model3
  lm(Model4, PCs) -> Model4
  #as.data.frame(anova(Model1, Model2))[2,6] -> P1
  #as.data.frame(anova(Model1, Model3))[2,6] -> P2
  which.min(c(BIC(Model1), BIC(Model2), BIC(Model3),BIC(Model4))) -> N
  Model = list(Model1, Model2, Model3, Model4)[[N]]
  
  
  as.data.frame(summary(Model)$coefficients) %>% rownames_to_column("Param") %>% mutate(PC = PC, Model_used= N) %>% as_tibble() -> Result_PC
  rbind(Results_PC_models, Result_PC) -> Results_PC_models
}
#Observation:
#PC1 is fundamentally affected by age (+, so - tel length), Sex (Female, -PC1, so +tel length) but when Females have CMV, their telomere length is specially reduced (no sign effect on guys?) 
#PC2 (mainly affects NK-cell), only age increase while CMV+ is significant, and decreases the telomere length
#PC3 symbolizes mainly shoerter memory Tcell tel length. We see age to have a negative effect somehow. Biggest impact in sex (females also have shorter) and CMV (positive has shorter)


#2. Isteand of PCA, use a Mixed-effect model
Telomere_model_input = left_join(left_join(Telomere_markers, CMV_presence), Covariates)
gather(Telomere_model_input, Cell_pop, Tel_length, Lymphocytes:`NK-cells`, factor_key=TRUE) %>% drop_na() -> Long_Telomere_model_input
Long_Telomere_model_input %>% mutate(Age = Age-min(Age) ) -> Long_Telomere_model_input_centered
nlme::lme(Tel_length ~ Age + Sex + CMV_presence + Age:CMV_presence + Sex:CMV_presence + Age:Sex , random=~1|ID, Long_Telomere_model_input) -> Mixed_model_telLength
summary(Mixed_model_telLength)
nlme::lme(Tel_length ~ Age + Sex + CMV_presence + Age:CMV_presence + Sex:CMV_presence + Age:Sex , random=~1|ID, Long_Telomere_model_input_centered) -> Mixed_model_telLength2
summary(Mixed_model_telLength2)

vif(Mixed_model_telLength) #high inflation values probably due to be using interaction terms
vif(Mixed_model_telLength2)
cor(select(Long_Telomere_model_input, - c(ID, Cell_pop)),method = "spearman" )

#As seen in the PCA analysis, it seems that only NK-cells present a trend for Age+CMV status
lm(Tel_length ~ Age + Sex + CMV_presence + Age:CMV_presence + Sex:CMV_presence + Age:Sex , filter(Long_Telomere_model_input,Cell_pop=="NK-cells")) -> NK_model_telLength
summary(NK_model_telLength)
lm(Tel_length ~ Age + Sex + CMV_presence + Age:CMV_presence + Sex:CMV_presence + Age:Sex , filter(Long_Telomere_model_input_centered,Cell_pop=="NK-cells")) -> NK_model_telLength2
summary(NK_model_telLength2)


PCs %>% gather(Cell_pop, Tel_length, PC1:PC6, factor_key=TRUE) -> Long_PCs

rbind(Long_Telomere_model_input, Long_PCs) -> Long_Telomere_model_input
paste(Long_Telomere_model_input$Sex, Long_Telomere_model_input$CMV_presence) -> Interaction
Long_Telomere_model_input %>% mutate(Sex_CMV = Interaction) %>% drop_na() -> Long_Telomere_model_input


Long_Telomere_model_input %>% 
  ggplot(aes(x=Sex_CMV, y= Tel_length )) + ggforce::geom_sina(alpha=0.5,aes(col = as.factor(Sex), shape =as.factor(CMV_presence))) + geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  facet_wrap(~ Cell_pop, scales= "free") + scale_color_viridis(discrete=TRUE) + scale_x_discrete( labels=c("0 0" = "Male/CNV(-)", "0 1" = "Male/CNV(+)","1 0" = "Female/CNV(-)", "1 1"= "Female/CNV(+)"))  #facet_grid(rows= vars(Cell_pop), scales= "free" )  

Long_Telomere_model_input %>%  mutate( Sex_CMV = ifelse( Sex_CMV == "0 0", "Male/CNV(-)", ifelse(Sex_CMV== "0 1", "Male/CNV(+)", ifelse(Sex_CMV=="1 0", "Female/CNV(-)", ifelse(Sex_CMV=="1 1", "Female/CNV(+)", Sex_CMV))))) %>%
  ggplot(aes(x=Age, y= Tel_length, col = as.factor(Sex_CMV))) + geom_point(alpha=0.5, aes(shape=as.factor(CMV_presence))) + theme_bw() + geom_smooth(method = "lm") +
  facet_wrap(~ Cell_pop, scales= "free")  + scale_color_viridis(discrete=TRUE) #facet_grid(rows= vars(Cell_pop), scales= "free" )  


##########blood cell analysis################
Model_input = left_join(left_join(all_cell_types, CMV_presence), Covariates)
Cell_absolute_results = tibble()
for (Cell in colnames(select(all_cell_types, -ID))){
  Cell = paste(c("`",Cell, "`"), collapse="")
  Model = as.formula(paste( c(Cell, " ~ Age + Sex + CMV_presence"), collapse="" ))
  as.data.frame(summary(lm(Model, Model_input))$coefficients) %>% rownames_to_column("Feature") %>% mutate(Cell = Cell) %>% as_tibble() ->temp_results
  rbind(Cell_absolute_results,temp_results ) -> Cell_absolute_results
}
Cell_absolute_results %>% filter(Feature == "CMV_presence") %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr")) %>% arrange(FDR) -> Cell_absolute_results
#Include more complex models
Cell_absolute_results_interactions = tibble()
for (Cell in colnames(select(all_cell_types, -ID))){
  Cell = paste(c("`",Cell, "`"), collapse="")
  Model = as.formula(paste( c(Cell,  " ~ Age + Sex + CMV_presence"), collapse="" ))
  Model2 = as.formula(paste( c(Cell, " ~ Age + Sex * CMV_presence"), collapse="" ))
  Model3 = as.formula(paste( c(Cell, " ~ Age * CMV_presence + Sex"), collapse="" ))
  Model4 = as.formula(paste( c(Cell, " ~ Age * CMV_presence * Sex"), collapse="" ))

  lm(Model, Model_input) -> Model1
  lm(Model2, Model_input) -> Model2
  lm(Model3, Model_input) -> Model3
  lm(Model4, Model_input) -> Model4
  which.min(c(BIC(Model1), BIC(Model2), BIC(Model3),BIC(Model4))) -> N
  Model = list(Model1, Model2, Model3, Model4)[[N]]

  as.data.frame(summary(Model)$coefficients) %>% rownames_to_column("Feature") %>% mutate(Cell = Cell) %>% as_tibble() ->temp_results
  rbind(Cell_absolute_results_interactions,temp_results ) -> Cell_absolute_results_interactions
}
Cell_absolute_results_interactions %>% filter(`Pr(>|t|)` < 0.05) %>% filter(! Feature %in% c("(Intercept)", "Age", "Sex" )) 
#There are 3 sex interactions:
#`CD4+ Naive CD45RA+ CD27+` are reduced with CMV, but in females this reductions is less strong, same trend is observed in `CD4+ T cells` and `Prol CD4+ Tconv`
#There is an age-CMV interaction in `Prol CD4+ Tconv`, where the decrease is even greater if CMV is presence
Model_input %>% gather(Cell_pop, Abundance, 1:(dim(all_cell_types)[2]-1), factor_key=TRUE) -> Model_input_long
paste(Model_input_long$Sex, Model_input_long$CMV_presence) -> Interaction
Model_input_long %>% mutate(Sex_CMV = as.factor(Interaction) ) ->  Model_input_long
Cell_absolute_results %>% filter(`Pr(>|t|)` < 0.05) %>% mutate( Cell = str_replace_all(Cell, "`", "") ) -> Significant_cell

#Model_input_long %>% filter(Cell_pop %in% Significant_cell$Cell) %>% 
#  ggplot(aes(x=Age, y= Abundance, col = as.factor(Sex))) + geom_point(alpha=0.5) + facet_grid(rows= vars(Cell_pop), scales= "free" ) + theme_bw() + geom_smooth(method = "lm")
Model_input_long %>% filter(Cell_pop %in% Significant_cell$Cell) %>% drop_na() %>% 
  ggplot(aes(x=Sex_CMV, y= Abundance )) + ggforce::geom_sina(alpha=0.5, aes(col = as.factor(Sex), shape =as.factor(CMV_presence))) + geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() + 
  facet_wrap(~ Cell_pop, scales= "free") + scale_color_viridis(discrete=TRUE) + scale_x_discrete( labels=c("0 0" = "Male/CNV(-)", "0 1" = "Male/CNV(+)","1 0" = "Female/CNV(-)", "1 1"= "Female/CNV(+)"))  #facet_grid(rows= vars(Cell_pop), scales= "free" )   #facet_grid(rows= vars(Cell_pop), scales= "free" )  


####### Other bio-aging markers ######
bioage = left_join(left_join(Other_age_markers, CMV_presence), Covariates)
biage_results_interactions = tibble()
for (bioa in colnames(select(Other_age_markers, -ID))){
  Model = as.formula(paste( c(bioa,  " ~ Age + Sex + CMV_presence"), collapse="" ))
  Model2 = as.formula(paste( c(bioa, " ~ Age + Sex * CMV_presence"), collapse="" ))
  Model3 = as.formula(paste( c(bioa, " ~ Age * CMV_presence + Sex"), collapse="" ))
  Model4 = as.formula(paste( c(bioa, " ~ Age * CMV_presence * Sex"), collapse="" ))
  
  lm(Model, bioage) -> Model1
  lm(Model2, bioage) -> Model2
  lm(Model3, bioage) -> Model3
  lm(Model4, bioage) -> Model4
  which.min(c(BIC(Model1), BIC(Model2), BIC(Model3),BIC(Model4))) -> N
  Model = list(Model1, Model2, Model3, Model4)[[N]]
  
  as.data.frame(summary(Model)$coefficients) %>% rownames_to_column("Feature") %>% mutate(Age_marker = bioa) %>% as_tibble() ->temp_results
  rbind(biage_results_interactions,temp_results ) -> biage_results_interactions
}
#CMV is not realted with methylation or dCT markers (!)


######Microbiome
Microbial3 = left_join(left_join(M4, CMV_presence), Covariates)
Microbial3 %>% drop_na() -> Microbial3

Bug_level %>% filter(Level == "s") -> Species
###Beta diversity
Distance = vegdist( select(filter(M4, ID %in% Microbial3$ID), Species$Bug) , method="euclidean") 
PCOA <- pcoa(Distance)
Variance_per_axis2 = PCOA$values["Relative_eig"]* 100
PCs = PCOA$vectors[,1:2]
left_join(left_join( mutate(as_tibble(PCs), ID = Microbial3$ID), CMV_presence) , Covariates) -> PCs2
ggplot() + geom_point(data=PCs2, aes(x=Axis.1, y=Axis.2, col= as.factor(CMV_presence) )) + theme_bw() + xlab(paste(c("PC1(", as.character(round(Variance_per_axis2[1,1], 2)), "%)"), collapse="")) + ylab(paste(c("PC2(", as.character(round(Variance_per_axis2[2,1], 2)), "%)"), collapse=""))  + scale_color_viridis(discrete=TRUE) -> PCA_plot2
adonis2(Distance ~ Age + Sex*CMV_presence , Microbial3 ) -> permanova_results
###Alpha diversity
shannon_diversity = vegan::diversity(select(Microbial2, Species$Bug), index = "shannon")
invsimpson_diversity = vegan::diversity(select(Microbial2, Species$Bug), index = "invsimpson")

left_join(left_join(CMV_presence, Covariates), tibble(ID=Microbial2$ID, Shannon=shannon_diversity, InvSimpson =invsimpson_diversity)) ->Diversity_tibble
summary(lm(Shannon ~ Age + Sex  + CMV_presence  ,Diversity_tibble)) ; summary(lm(InvSimpson ~ Age + Sex  + CMV_presence  ,Diversity_tibble)) 

Diversity_tibble %>% gather(Diversity, Diversity_value, c("Shannon", "InvSimpson")) %>% ggplot(aes(x=as.factor(CMV_presence), y=Diversity_value)) + ggforce::geom_sina() +
  geom_boxplot() + theme_bw() + facet_wrap(~Diversity, scales = "free")



###Association
Prevalence = apply(select(Microbial2, -ID), 2, FUN= function(x){ y = as.numeric(x>0) ; return(100*sum(y)/length(x)) } )
names(Prevalence[Prevalence > 10]) -> Prevalent

Microbial_associations = tibble()
for (Bug in Prevalent){
  Bug2 = paste(c("`",Bug, "`"), collapse="")
  Model = as.formula(paste( c(Bug2,  " ~ Age + Sex + CMV_presence"), collapse="" ))
  Model2 = as.formula(paste( c(Bug2, " ~ Age + Sex * CMV_presence"), collapse="" ))
  Model3 = as.formula(paste( c(Bug2, " ~ Age * CMV_presence + Sex"), collapse="" ))
  Model4 = as.formula(paste( c(Bug2, " ~ Age * CMV_presence * Sex"), collapse="" ))
  
  lm(Model, Microbial3) -> Model1
  lm(Model2, Microbial3) -> Model2
  lm(Model3, Microbial3) -> Model3
  lm(Model4, Microbial3) -> Model4
  which.min(c(BIC(Model1), BIC(Model2), BIC(Model3),BIC(Model4))) -> N
  Model = list(Model1, Model2, Model3, Model4)[[N]]
  
  as.data.frame(summary(Model)$coefficients) %>% rownames_to_column("Feature") %>% mutate(Age_marker = Bug) %>% mutate(Model=N) %>% as_tibble() ->temp_results
  rbind(Microbial_associations,temp_results ) -> Microbial_associations
}
Microbial_associations %>% filter(! Feature %in% c("(Intercept)", "Age", "Sex" )) %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr" )) %>%
  arrange(FDR) 

##Since association will be weak, in order to avoud multiple test correction lost of power, we performe a prediction task using taxonomy balances
Run_balance_analysis = function(Balance_input, Pseudo, with_cov = T, lambda=1){
  set.seed(50)
  Balance_input2 = as.data.frame(select(Balance_input, -c(ID, CMV_presence, Age, Sex)) + Pseudo/2)
  #Split dataset
  tf$random$set_seed(0)
  trainIndex <- sample(1:nrow( Balance_input2), 0.8 * nrow(Balance_input2))
  #Train set
  xTrain <- Balance_input2[trainIndex,]
  yTrain <- as.factor(Balance_input$CMV_presence[trainIndex])
  #Test set
  xTest <- Balance_input2[-trainIndex,]
  yTest <-as.factor(Balance_input$CMV_presence[-trainIndex])
  
  if (with_cov == T){
    partial <- glm( CMV_presence ~ Age + Sex , data= Balance_input[trainIndex,] , family='binomial')
    model=codacore( xTrain ,   yTrain ,offset = predict(partial), logRatioType = 'balances', lambda = lambda) #offset in logit space
    
    print("Train results")
    partialAUC = pROC::auc(pROC::roc(yTrain, predict(partial), quiet=T))
    codacoreAUC = model$ensemble[[1]]$AUC
    cat("Train set AUC =", codacoreAUC, "\n")
    cat("AUC gain:", round(100 * (codacoreAUC - partialAUC)), "%", "\n")
    
    ##Test
    yHat_partial <- predict(partial, newdata = Balance_input[-trainIndex,])
    partialAUC = pROC::auc(pROC::roc(yTest, yHat_partial, quiet=T))
    yHatLogit <- yHat_partial + predict(model, xTest, logits=T)
    #yHat <- yHatLogit > 0
    testAUC <- pROC::auc(pROC::roc(yTest, yHatLogit, quiet=T))
    cat("Test set AUC =", testAUC, "\n")
    cat("AUC gain test set:", round(100 * (testAUC - partialAUC)), "%", "\n")
  } else{
    model=codacore( xTrain ,   yTrain, logRatioType = 'balances', lambda = lambda)
    print("Train results")
    codacoreAUC = model$ensemble[[1]]$AUC
    cat("AUC train:", round(100 * (codacoreAUC)), "%", "\n")
    #Test
    yHat <- predict(model, xTest, logits=F)
    cat("Test set AUC =", pROC::auc(pROC::roc(yTest, yHat, quiet=T)), "\n")
    
  }
  
  Plot_distr = plot(model)
  Plot_roc = plotROC(model)
  Numerator = colnames(xTrain)[getNumeratorParts(model, 1)]
  Denominator = colnames(xTrain)[getDenominatorParts(model, 1)]
  
  cat("Numerator choice: ", Numerator, "\n")
  cat("Denominator choice: ", Denominator, "\n")
  
  as_tibble(getLogRatios(model, xTest)) %>% mutate(CMV_presence = yTest) -> Test_ratios
  as_tibble(getLogRatios(model, xTrain)) %>% mutate(CMV_presence = yTrain) -> Train_ratios
  
  Test_ratios  %>% ggplot(aes(x=CMV_presence, y= `log-ratio1`, fill=CMV_presence) ) +
    geom_boxplot(outlier.shape = NA) + ggforce::geom_sina() + theme_bw() + coord_flip() + ylab("Log-ratio score") + theme(legend.position = "none") + scale_fill_manual(values = c("#FFA500","#ADD8E6")) +ggtitle("Test set")-> Test_distr
  Train_ratios %>% mutate(CMV_presence = yTrain)  %>% ggplot(aes(x=CMV_presence, y= `log-ratio1`, fill=CMV_presence) ) + 
    geom_boxplot(outlier.shape = NA) + ggforce::geom_sina()  + theme_bw() + coord_flip() + ylab("Log-ratio score") + theme(legend.position = "none") + scale_fill_manual(values = c("#FFA500","#ADD8E6")) +ggtitle("Train set") -> Train_distr
  
  wilcox.test(filter(Test_ratios, CMV_presence==1)$`log-ratio1`, filter(Test_ratios, CMV_presence==0)$`log-ratio1` ) -> Dif_test
  wilcox.test(filter(Train_ratios, CMV_presence==1)$`log-ratio1`, filter(Train_ratios, CMV_presence==0)$`log-ratio1` ) -> Dif_train
  Distributions = Train_distr + Test_distr
  print(Distributions)
  cat("Wilcox test of ratio. Test:", Dif_test$p.value, ", Train:", Dif_train$p.value,"\n" )
}

library (codacore)
library("tensorflow")
#For installation I had to use: install_tensorflow(method = 'conda', envname = 'r-reticulate') ; after in creating the envname in the terminal "conda create -n r-reticulate"
#Microbial2 has relative abundances of all taxonomic level, In the examples included all taxonomic levels are also included togeter.
as_vector(select(Microbial2, -ID) ) %>% as_vector() -> Abundances ; Pseudo = min(Abundances[!Abundances == 0])
left_join(left_join(Microbial2, CMV_presence), Covariates) %>% drop_na()  -> Balance_input

#Prevalence filter
left_join(left_join(select(Microbial2, c("ID",Prevalent)), CMV_presence), Covariates) %>% drop_na()  -> Balance_input

Run_balance_analysis(Balance_input, Pseudo)

#Although we find a ratio that seems interesting in the train dataset, it is not replicable in test. 





###Mortality

dead = as_tibble(read.table("~/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/Survival_/LLD_passed.txt",header = T))
Info = read_tsv("~/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/Survival_/baseline.txt")
End_date = as.Date("2021-06-08")
dead   %>% mutate (DATE = as.Date(paste(deathyear,deathmonth,"01",sep="-") )) -> dead
colnames(dead)[2] = "ID"

Info %>% mutate(First_time = as.Date(paste(bl1year,"01-01",sep="-")) ) -> Info #,Second_time = ifelse(LLDEEPID %in% dead$LLDeep_SampleID), dead[LLDeep_SampleID==LLDEEPID,]$DATE,End_date-First_time)
Second_time = c()
for (ID in Info$LLDEEPID){
  Info[Info$LLDEEPID==ID,]$First_time -> First_time
  if(ID %in% dead$LLDeep_SampleID){
    dead[dead$LLDeep_SampleID==ID,]$DATE -First_time -> Datet
  } else { End_date-First_time -> Datet }
  Second_time = c(Second_time, Datet)
}
Info %>% mutate(Second_time = Second_time, Status = ifelse(LLDEEPID %in% dead$ID, 1, 0)) -> Info

library("survival")
select(mutate(Info, ID= LLDEEPID) , c(ID, Second_time, Status))  -> Info2
left_join(left_join(Info2, CMV_presence),Covariates) %>% drop_na() %>% select(-ID) -> Info_mortality
coxph(Surv(Second_time, Status) ~ Age + CMV_presence + Sex, data=Info_mortality) -> res.cox
test.ph <- cox.zph(res.cox) #Doesnt work?
confint(res.cox)["CMV_presence",] -> CI
as.data.frame(summary(res.cox)$coefficients)["CMV_presence",c(1,5)] -> Res_mortality
Res_mortality %>% mutate(`2.5 %`= CI[1], `97.5 %`= CI[2] ) -> Res_mortality
ggplot(Res_mortality) + geom_point(aes(x=as.factor(1), y=coef )) + theme_bw()

