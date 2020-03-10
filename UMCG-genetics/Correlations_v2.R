#TMAO is a metabolite canonically produced by the gut microbiome and that can be obtain through diet or through biosintetic production from dietary products.
#In this script different associations will be made using different OMIC data and TMAO-related metabolites.

#Author: S. Andreu-Sanchez
#28/02/2020
setwd("~/PhD/WORK/TMAO_metagenomics")
##################
###Package load###
##################
#General
library(tidyverse) #Data handling
library(microbiome) #cdr transformation source("Microbiome_function.R")
library(reshape2)
library(patchwork)
library(pheatmap)
#Statistical learning
library(randomForest)
library(glmnet)
library(caret)
library(qdapTools)


###################
###Load Functions##
###################
source("Functions_TMAO.R")

###################
####Load Data######
###################

#TMAO and other metabolites of interest 
Dependent_metabolites = read_tsv("Pheno_file.tsv")
#Covariates of interest (BMI, age, sex)
Covariates =  read_tsv("Cov_file.tsv")
#Clinical and questionary data
Clinical_Questionaries =  read_tsv("Full_proccesed_data_LLD_2017_Jing.txt")
#MSS: Microbial Abundance
Species_abundance = read_tsv("data_LLD_baseline_1135samples_698species.txt")
Genus_abundance = read_tsv("data_LLD_baseline_1135samples_698genus.txt")
Family_abundance = read_tsv("data_LLD_baseline_1135samples_698family.txt")
Phylum_abundance = read_tsv("data_LLD_baseline_1135samples_698phylum.txt")
#MSS: Abundance of genes, gene clusters, pathways...
Gene_abundance = read_tsv(file = "Gene_abundance")
Cluster_abundance = read_tsv(file = "../TMAO_cutC/Total_table.tsv")
Gene_family_abundance = read_tsv("data_LLD_baseline_1135samples_489pathway.txt")

#Diet 
Diet_regressors = read_tsv(file = "20150722_Diet__1135patients.txt")

#Metabolome
Lipids_aa = read_tsv("LLD_nmr_total.txt")
untargeted_metabolome = read_tsv("data_1442samples_LLD_baseline_1183plasma_metabolites.txt")
key_metabolome = read_tsv("key_lld_1183meta_annotation.txt")
####################################

#1. Assess relations among Metabolites of interest.
  ##1.1 Make correlation matrix and distributions
  ###Samples with high levels of TMAO are saved in TMAO_outlayers
TMAO_outlayers = Metabolite_exploration(Dependent_metabolites)
##1.2 Address statistical significance and direction of correlations
  ###Transformation of metabolites to Normality
Normalized_dependent_metabolites = apply(Dependent_metabolites[2:length(Dependent_metabolites)], 2, Rank_norm) %>% as_tibble()
Normalized_dependent_metabolites %>% mutate(ID = Dependent_metabolites$ID) -> Normalized_dependent_metabolites
Datasets = Match_dataset(Metabolites = Normalized_dependent_metabolites, Regressor = Normalized_dependent_metabolites, Covariates = Covariates )

Association1 = Association_metabolites(select(Normalized_dependent_metabolites,-ID), with_Covariates = F)
Association1 %>% mutate(Beta = as.numeric(Beta), Pvalue = as.numeric(Pvalue)) -> Association1
Fig_association1 = Make_heatmap(Association1, "DependentvsDependent_heatmap")
#Assocaition2 = Association_metabolites(select(Datasets[[1]],-ID), with_Covariates = T, Covariates = select(Datasets[[2]],-ID))
Association3 = Association_metabolites(select(Datasets[[1]],-ID), with_Covariates = "BMI", Covariates = select(Datasets[[2]],-ID))
Association3 %>% mutate(Beta = as.numeric(Beta), Pvalue = as.numeric(Pvalue)) -> Association3
Fig_association3 = Make_heatmap(Association3, "DependentvsDependent_adjutsted_heatmap")

#ggsave(filename = "Model_stats/DependentvsDependent_heatmap.pdf", Fig_association1, scale = 3)
#ggsave(filename = "Model_stats/DependentvsDependent_adjutsted_heatmap.pdf", Fig_association3, scale = 3)

####Metabolite are strongly associated among them, with high betas (except TMAO, more complex relation with rest?)


#2. Asses the association of covariates and data
Dependents = select(Datasets[[1]],-ID)
Regressors = select(Datasets[[2]],-ID)

Cov_relation = tibble()
for(Dependent in colnames(Dependents)){
  Dependents %>% select(Dependent) %>% as_vector() -> Dependent_variable
  as_tibble(summary(lm(Dependent_variable ~ rank(Regressors$BMI)))$coefficients) -> Summary_BMI
  as_tibble(summary(lm(Dependent_variable ~ rank(Regressors$Age)))$coefficients) -> Summary_age
  as_tibble(summary(lm(Dependent_variable ~ as.factor(Regressors$Gender)))$coefficients) -> Summary_gender
  for (Summary in seq(1:3)){
    Name = c("BMI", "age", "gender")[Summary]
    Summary = list(Summary_BMI, Summary_age, Summary_gender)[[Summary]]
    Beta = Summary$Estimate[2]
    pvalue = Summary$`Pr(>|t|)`[2]
    stand = Summary$`Std. Error`[2]
    Results = as_tibble(t(c(Dependent, Name, Beta, stand, pvalue)))
    colnames(Results) = c("Regressor", "Covariate", "Beta", "SE", "Pvalue")
    Cov_relation = rbind(Cov_relation, Results)
  }
  
}

  ##Are these associations independent of Kidney function?

Compute_GFR = function(x){
  x %>% mutate(GFR= ifelse(Gender==1, ifelse(Creatinine<= 0.7,144*(Creatinine/0.7)-0.329 *(0.993)^Age, 144*(Creatinine/0.7)-1.209 * (0.993)^Age), NA)) -> x
  x %>% mutate(GFR= ifelse(Gender==0, ifelse(Creatinine<=0.9,141*(Creatinine/0.9)-0.411 * (0.993)^Age, 141*(Creatinine/0.9)-1.209 * (0.993)^Age), GFR)) -> x
  return(x)
}

Datasets = Match_dataset(Regressor = Clinical_Questionaries, Metabolites = Normalized_dependent_metabolites, Covariates )
Datasets[[2]] %>% mutate(Creatinine = Datasets[[3]]$`Creatinine (umol/L)`) -> Covariates
#Covariates %>% mutate(GFR= ifelse(Gender==0, 186 * (Creatinine/88.4)^-1.154 * (Age)^-0.203, 186 * (Creatinine/88.4)^-1.154 * (Age)^-0.203*0.742)) -> Covariates #https://ukidney.com/nephrology-resources/egfr-calculator
Covariates = Compute_GFR(mutate(Covariates, Creatinine = Creatinine*0.0113))
Dependents = select(Datasets[[1]],-ID)
Regressors = select(Covariates,-ID)

Cov_independent_kidney = tibble()
for(Dependent in colnames(Dependents)){
  Dependents %>% select(Dependent) %>% as_vector() -> Dependent_variable
  as_tibble(summary(lm(Dependent_variable ~ rank(Regressors$BMI) + Regressors$GFR + Regressors$Age + Regressors$Gender))$coefficients) -> Summary_BMI
  as_tibble(summary(lm(Dependent_variable ~ rank(Regressors$Age) + Regressors$GFR + Regressors$BMI + Regressors$Gender))$coefficients) -> Summary_age
  as_tibble(summary(lm(Dependent_variable ~ as.factor(Regressors$Gender) + Regressors$GFR + Regressors$BMI + Regressors$Age ))$coefficients) -> Summary_gender
  as_tibble(summary(lm(Dependent_variable ~ rank(Regressors$GFR + Regressors$GFR + Regressors$BMI + Regressors$Gender + Regressors$Age)))$coefficients) -> Summary_GFR
  for (Summary in seq(1:4)){
    Name = c("BMI", "age", "gender","GFR")[Summary]
    Summary = list(Summary_BMI, Summary_age, Summary_gender,Summary_GFR)[[Summary]]
    Beta = Summary$Estimate[2]
    pvalue = Summary$`Pr(>|t|)`[2]
    stand = Summary$`Std. Error`[2]
    Results = as_tibble(t(c(Dependent, Name, Beta, stand, pvalue)))
    colnames(Results) = c("Regressor", "Covariate", "Beta", "SE", "Pvalue")
    Cov_independent_kidney = rbind(Cov_independent_kidney, Results)
  }
  
}

#without GFR
Cov_relation %>% mutate(Metabolite = Covariate ,Pvalue = as.numeric(Pvalue), Beta = as.numeric(Beta)) -> Cov_relation
Fig_covariates = Make_heatmap(Cov_relation, "Covariates_heatmap")

#With GFR
Cov_independent_kidney %>% arrange(Regressor, Covariate) -> Cov_independent_kidney
Cov_independent_kidney %>% mutate(Metabolite = Covariate ,Pvalue = as.numeric(Pvalue), Beta = as.numeric(Beta)) -> Cov_independent_kidney
Fig_covariates_kideney = Make_heatmap(Cov_independent_kidney, "Covariates_kideney_heatmap")

# We observe that in TMAO age is highly associated to higher levels, BMI is not significantly correlated, and GFR and gender are somewhat associated.
#Interestingly Gender is really associated to all other metabolites in the same direction. The little effect on TMAO might be due to higher expression of FMO3 in women (correlate with RNAseq?)
# BMI appear to not affect TMAO, Choline or y-butyrobetaine after corrections.Still associated with L-carnitine. ASsociated with all before adjustments.


### Are TMAO levels different between indifiduals with high Kidney functiona and low kidney fucntion?
Covariates %>% ggplot() + geom_density(aes(GFR))
Covariates %>% mutate(quantile_GFR = ntile(GFR, 4)) -> Quantiles_GFR
First_quantile = filter(Quantiles_GFR,quantile_GFR==1)$ID
Fourth_quantile = filter(Quantiles_GFR,quantile_GFR==4)$ID
ggplot() + geom_violin(aes(x="first",y=filter(Dependent_metabolites, ID %in% First_quantile)$TMAO), fill="red") + geom_violin(aes(x="Fourth",y=filter(Dependent_metabolites, ID %in% Fourth_quantile)$TMAO), fill="blue",alpha=0.5)+ theme_bw()
wilcox.test(filter(Dependent_metabolites, ID %in% First_quantile)$TMAO, filter(Dependent_metabolites, ID %in% Fourth_quantile)$TMAO) 
#3. Association with other clinical features

Regr = Prepare_clinical(Clinical_Questionaries, Metabolites = Dependent_metabolites ,Covariates  = Covariates )
Results_Phenotypes = Iterate_Metabolites(Regr[[1]],Regr[[2]],Regr[[3]], "Clinical_Questionaries",FDR_iter=100)

Fig_phenos = Make_heatmap(Results_Phenotypes, "ClinicalPhenotypes_heatmap" )

#TMAO: Associated with creatinine in serum (<0.05FDR), and with a FDR of 0.2 with excreted Creatinine in opposite direction (consistent)
#Betaine: Many associations. CRP (inflmation marker) is NEGATIVELY correlated. It is also negatively correlated with triglycerides, blood preassure, trhombocytes, Creatinine and creatinine excretion, inactivity, hemoglobine or lymphocytes
#This points to a role of Betaine in protection, which makes sense: Betaine works by preventing the build-up of an amino acid called homocysteine. This amino acid can harm blood vessels and contribute to heart disease, stroke, or circulation problems.
#Choline: Positively associated with creatinine in plasma, with triglyceride levels, negatively associated with inflmmation
#L-carnitine: only strong association with inflammation marker (negative correlation), other positive corelation with tryglycerides and negative to blood pressure
#y-butyro: positive correlated with creatinen in blood, negative correlated with hsCRP, negative correalted pulse.Others: positive with QT, negative with preassure or leykocytes

##All negatively correlated with hsCRP except TMAO (low association)

#4. Association with diet
Regr = Prepare_diet(Diet_regressors, Metabolites = Dependent_metabolites, Covariates = Covariates)
Results_Diet = Iterate_Metabolites(Regr[[1]],Regr[[2]],Regr[[3]], "Diet",FDR_iter=100)
Results_Diet  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_diet
#TMAO associated with how_often_fish, fish and alcohol
#Betaine no significant associations with diet
#Choline neither
#Lcarnitine associated with coffee, alcohol and meat
#y-butyro: no significant association (lowest is meat)

Fig_Diet = Make_heatmap(Results_Diet, "Diet_heatmap")
#ggsave(filename = "Model_stats/Diet_heatmap.pdf", Fig_Diet, scale = 3)



#5. Association with microbes
##5.1 Species
Input_MGS = Prepare_Metagenome_species(Species_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
Results_MGS = Iterate_Metabolites(Input_MGS[[1]],Input_MGS[[2]],Input_MGS[[3]], "MGS_species")
Results_MGS  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_MGS_genus
#TMAO: best associations in 0.14 FDR, positive correlation with Prevotella_copri, Eubacterium_biforme and Paraprevotella
#Betaine: Strongest associations with 0.24 (Prevotella_timonensis)
#Choline: No significant
#L-carnitine: No significant associations
#butyrobetaine: best association FDR 0.09 with PRevotella_histicola

#Does y-butyrobetaine (caiB) or Betaine (lcdH + something) (both produced by moo) come from enzymes from prevotella?
Fig_species = Make_heatmap(Results_MGS, "Species_heatmap")
#ggsave(filename = "Model_stats/Species_heatmap.pdf", Fig_species, scale = 3)


##5.2 Genus
Input_MGS2 = Prepare_Metagenome_genus(Genus_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
Results_MGS2 = Iterate_Metabolites(Input_MGS2[[1]],Input_MGS2[[2]],Input_MGS2[[3]], "MGS_genus")
Results_MGS2  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_MGS_genus

#TMAO: Two positive associations (FDR=0.04) Prevotella and PAraprevotella
#Betaine: No associations
#Choline: Top FDR 0.13, Anaerostipes and LAchnospiraceae
#L-carnitine: no associations
#y-butyrobetaine: no associations

Fig_genus = Make_heatmap(Results_MGS2, "Genus_heatmap")
#ggsave(filename = "Model_stats/Genus_heatmap.pdf", Fig_genus, scale = 3)


##5.2 Family
Input_MGS3 = Prepare_Metagenome_genus(Family_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
Results_MGS3 = Iterate_Metabolites(Input_MGS3[[1]],Input_MGS3[[2]],Input_MGS3[[3]], "MGS_family")
Results_MGS3  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_MGS_genus

Make_heatmap(Results_MGS3, "Family_heatmap") -> Heatmap_family
#ggsave(filename = "Model_stats/Family_heatmap.pdf", Heatmap_family, scale = 3)
#TMAO: Best hit Prevotellaceaea (FDR:0.09), positive correaltion
#Betaine: Best hit is clostridiales, but 0.13 FDR
#Choline: No good hits (0.2 FDR)   No microbial origin
#L-carnitine: No good hits (0.23 FDR) No microbial origin
#y-butyrobetaine: A bunch of ~0.18 FDR


#6 Microbial functional annnotation
##6.1 Gene family 
Input_pt = Prepare_Pathways(Gene_family_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
Results_Pathway = Iterate_Metabolites(Input_pt[[1]], Input_pt[[2]], Input_pt[[3]], "MGS_pathway")
Results_Pathway  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_MGS_pathways


Fig_pathway = Make_heatmap(Results_Pathway, "Pathway_heatmap")
#ggsave(filename = "Model_stats/Pathway_heatmap.pdf", Fig_pathway, scale = 3)


##6.2 Pathways
##6.3 Gene clusters

#7 Gene abunance (shorbred)
Gene_markers = read_tsv(file = "../TMAO_cutC/Total_table.tsv")
Gene_markers %>% group_by(ID) %>% summarise(Count = sum(Count)) -> Gene_markers
Input_shortbred = Prepare_cutC_shortbred(Gene_markers, Metabolites = Dependent_metabolites, Covariates = Covariates)
Results_shortbred = Iterate_Metabolites(Input_shortbred[[1]], Input_shortbred[[2]], Input_shortbred[[3]], "MGS_shortbred2",FDR_iter = 0)

Results_shortbred %>% mutate(log10Pval = -log10(Pvalue)) %>% select(Metabolite, Regressor, log10Pval) %>% spread(Metabolite, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results
pheatmap(wide_results, "Shortbred_heatmap") -> Fig_shortbred
#ggsave(filename = "Model_stats/Shortbred_heatmap.pdf", Fig_shortbred, scale = 3)



#7. Genetics
#....
#8. Metabolomics
##Lipds and amino-acid related metabolites
Input_metabolome = Prepare_metabolome(Lipids_aa, Metabolites = Dependent_metabolites, Covariates=Covariates)
Results_metabolome = Iterate_Metabolites(Input_metabolome[[1]], Input_metabolome[[2]], Input_metabolome[[3]], "Metabolome")
Results_metabolome  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
#Get names from variables somewhere...
#TMAO: 18 associations:  Top four are omega-3, 1 degree of insaturation, correlated to Intermediate denstiy (between low and very low) tryglycerides and negatively correlated to esters
#low density cholesterol tryglycerides ;  positive correlation with creatinine ;  and to Acetoacetate, among others
#Betaine: Role in cholesterol homeostasis/protection : Negative association with all fats except cholesterol in VLDL. Negative correlation with apoa1 (from HDL), also negative total phosphoglycerides 186 associations
#Choline: Positive correlation with fats (with HDL) 104 associations
#L-Carnitine: Negative correlation with HDL, positive correaltions with VLD, 117 (first ones are aminoacids)
# y-butyrobetaine: 52 associations (fist ones amino-acids, creatinine, HDL)

Fig_targeted_metabolites = Make_heatmap(Results_metabolome, "targeted_heatmap")
#ggsave(filename = "Model_stats/targeted_heatmap.pdf", Fig_targeted_metabolites, scale = 3)



##Need some expert in metabolism to check all this...

##Untargeted metabolite approach
untargeted_metabolome
Input_metabolome_un = Prepare_metabolome(untargeted_metabolome, Metabolites = Dependent_metabolites, Covariates=Covariates)
Results_metabolome_un = Iterate_Metabolites(Input_metabolome_un[[1]], Input_metabolome_un[[2]], Input_metabolome_un[[3]], "Metabolome_untargetd")

select(key_metabolome,c(meta, name)) -> Dic_key
lookup(terms=Results_metabolome_un$Regressor, Dic_key$meta, Dic_key$name) -> metabolite_names
Results_metabolome_un$Regressor = metabolite_names


Fig_untargeted_metabolites = Make_heatmap(Results_metabolome_un, "untargeted_heatmap")
#ggsave(filename = "Model_stats/untargeted_heatmap.pdf", Fig_untargeted_metabolites, scale = 3)

##Need some expert in metabolism to check all this...


#X Random Forest


FeatureSelection_RandomForest = function(Input_x, Dependent_metabolites, Covariates, Fold){
  
  set.seed(111)
  ##10-fold cross-validation, using OoB MSE as error measure
  flds <- createFolds(Input_x, k = Fold, list = TRUE, returnTrain = FALSE)
  Min_loss = vector()
  for (i in seq(Fold)){
    Test = Input_x[ flds[[i]], ]
    Train_i = Input_x[-flds[[i]], ]
    
    Train_i %>% select(-c("Age","Gender","BMI","Creatinine","GFR")) -> Train
    Input_train = Match_dataset(Train, Metabolites = Dependent_metabolites, Covariates = Covariates) 
    Results_train = Iterate_Metabolites(select(Input_train[[1]],c(ID,TMAO)),Input_train[[2]],Input_train[[3]], "MGS_species_cv", FDR_iter = 0)
    Results_train %>% filter(Pvalue < 0.01) -> Feature_selection
    #Add covaraites?
    Input_train[[3]] %>% select(Feature_selection$Regressor) %>% cbind(select(Input_train[[2]],-ID)) -> Train
    RF_train = Fit_RandomForest(Train, Input_train[[1]]$TMAO)
    Model_RF = RF_train[[2]]
    
    Input_test = Match_dataset(Test, Metabolites = Dependent_metabolites, Covariates = Covariates) 
    predTest = predict(Model_RF, Input_test[[3]])
    Loss_test = RMSE(predTest,Input_test[[1]]$TMAO)
    if (length(Min_loss) == 0){
      Min_loss = Loss_test
      Model = Model_RF
      Importance = RF_train[[1]]
    }else if (Loss_test < Min_loss){
      Min_loss = Loss_test
      Model = Model_RF
      Importance = RF_train[[1]]
    }
    
  }
  return(list(Min_loss, Model, Importance))
}
#Species
S1 = Prepare_Metagenome_species(Species_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
S1[[3]] %>% cbind(select(S1[[2]],-ID)) %>% as_tibble() %>% drop_na() -> Input_Species
FeatureSelection_RandomForest(Input_x = Input_Species,Dependent_metabolites = Dependent_metabolites ,  Covariates= Covariates, Fold = 10)
#Diet
S1 = Prepare_diet(Diet_regressors, Metabolites = Dependent_metabolites, Covariates = Covariates)
S1[[3]] %>% cbind(select(S1[[2]],-ID)) %>% as_tibble() %>% drop_na() -> Input_Diet
FeatureSelection_RandomForest(Input_x = Input_Diet,Dependent_metabolites = Dependent_metabolites ,  Covariates= Covariates, Fold = 10)
#Genus
S1 = Prepare_Metagenome_genus(Genus_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
S1[[3]] %>% cbind(select(S1[[2]],-ID)) %>% as_tibble() %>% drop_na() -> Input_Diet
FeatureSelection_RandomForest(Input_x = Input_Diet,Dependent_metabolites = Dependent_metabolites ,  Covariates= Covariates, Fold = 10)
#Species + Diet
S1 = Prepare_Metagenome_species(Species_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
S2 = Prepare_diet(Diet_regressors, S1[[3]], Covariates)
Species_diet = S2[[1]] %>% cbind(select(S2[[3]], -ID))
Species_diet %>% cbind(select(S2[[2]],-ID)) %>% as_tibble() %>% drop_na() -> Input_spdi
FeatureSelection_RandomForest(Input_x = Input_spdi,Dependent_metabolites = Dependent_metabolites ,  Covariates= Covariates, Fold = 10)
#Metabolites
S1 = Prepare_metabolome(Lipids_aa, Metabolites = Dependent_metabolites, Covariates=Covariates)
S1[[3]] %>% cbind(select(S1[[2]],-ID)) %>% as_tibble() %>% drop_na() -> Input_metabolome_rf
FeatureSelection_RandomForest(Input_x = Input_metabolome_rf,Dependent_metabolites = Dependent_metabolites ,  Covariates= Covariates, Fold = 10)

#X Lasso
S1 = Prepare_Metagenome_species(Species_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
S2 = Prepare_diet(Diet_regressors, S1[[3]], Covariates)
Species_diet = S2[[1]] %>% cbind(select(S2[[3]], -ID))
Species_diet %>% cbind(select(S2[[2]],-ID)) %>% as_tibble() %>% drop_na() -> Input_spdi
S1[[1]] %>% filter(ID %in% Input_spdi$ID) -> Dependent

#Iterate for each metabolite

for (Metabolite in colnames(select(Dependent, -ID))){
  Dependent %>% select(Metabolite) %>% as_vector() -> Dep
  Lasso_list = Fit_lasso(select(Input_spdi, -ID), Dep)
  
  Lasso_Model = Lasso_list[[1]]
  Betas = Lasso_list[[2]]
  Variance_explained = Lasso_list[[3]]
  
  Non_zero_betas = Betas[which(Betas != 0)]
  Non_zero_betas_names = Betas@Dimnames[[1]][which(Betas != 0 ) ]
  DF_beta = tibble(Feature = Non_zero_betas_names, Beta =  Non_zero_betas)
  write_tsv(x = DF_beta, path = paste(c("Model_stats/Lasso_betas_", Metabolite, ".tsv"), collapse = ""))
  print(paste(c(Variance_explained, Metabolite), collapse = " "))
  Plot = plot(Lasso_Model)
  print(Plot)
}





#Questions 
###What is the effect of Age, Gender, BMI and Kideney function in levels of TMAO.
#####Age is believed to be positively correlated, why? is this kidney independent?
    ### Strong positive correlation with all TMAO metabolites, independnt of any other covariate, even kideney. What is the physiology?
#####Gender appears to have different effects depending on the source. Correct for other covariates and check the gender. Is the expression of FMO3 in liver also dependent of gender?
    ### Association in TMAO independent of any other variables (females have fewer). But, stronger association in all other metabolites. Can it be that even if females have less, there is other mechanism that compensates it? According to GTEX, female expression of FMO3 is higher in liver, which may paliate with the gender differences (A bit). https://gtexportal.org/home/gene/FMO3
#####Does BMI directly have an effect on TMAO levels?
#####What about Kideney function, does it have an independent effect on TMAO levels?

###What are the major sources of TMAO?
#####Are specific bacterial lineages of importance for its presence in blood?
#####What about host genotype?
#####Does diet really affect their levels? Is it moo-dependnt or independent?
#####Do Kideney function levels affect the effect of any of the factors?

##What role plays TMAO?
#####Which metabolites and proteins are linked to its presence in blood?
#######Check correlation of TMAO and bile-acids
#####Are these linked to known TMAO phenotypes?
