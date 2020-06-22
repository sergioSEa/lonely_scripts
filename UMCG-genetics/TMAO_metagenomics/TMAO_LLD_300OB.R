#TMAO correlation with microbiome

setwd("~/PhD/WORK/TMAO_metagenomics")
library(tidyverse) #Data handling
library(microbiome) #cdr transformation source("Microbiome_function.R")
library(reshape2)
library(patchwork)
library(pheatmap)

source("Functions_TMAO.R")

##LLDeep data
Dependent_metabolites = read_tsv("Pheno_file.tsv")
Covariates =  read_tsv("Cov_file.tsv")

Species_abundance = read_tsv("data_LLD_baseline_1135samples_698species.txt")
Genus_abundance = read_tsv("data_LLD_baseline_1135samples_698genus.txt")
Family_abundance = read_tsv("data_LLD_baseline_1135samples_698family.txt")
Phylum_abundance = read_tsv("data_LLD_baseline_1135samples_698phylum.txt")

Gene_abundance = read_tsv(file = "Gene_abundance")
Cluster_abundance = read_tsv(file = "../TMAO_cutC/Total_table.tsv")
Gene_family_abundance = read_tsv("data_LLD_baseline_1135samples_489pathway.txt")

Clinical_Questionaries = read_tsv("LLD_1135subjects_allphenotypes.txt")
Diet_start = which(colnames(Clinical_Questionaries)=="alcohol_products_log") ; Diet_end = which(colnames(Clinical_Questionaries)=="ACE_inhibitor")-1
Diet_regressors = Clinical_Questionaries[,Diet_start:Diet_end] ; Diet_regressors %>% mutate(ID = Clinical_Questionaries$ID) -> Diet_regressors
Clinical_Questionaries = Clinical_Questionaries[,1:Diet_start -1 ]

CutC_LLD = read_tsv("../TMAO_cutC/Total_table.tsv")
##300 OB data
Dependent_metabolites_300O = read_tsv(file = "300OB/330_O_metabolite_measures.tsv")
colnames(Dependent_metabolites_300O)[1] =  "ID"
Dependent_metabolites_300O %>% drop_na() -> Dependent_metabolites_300O
sapply(select(Dependent_metabolites_300O, -ID), FUN= function(x){transform(as.numeric(x))}) %>% as_tibble() %>% mutate(ID = Dependent_metabolites_300O$ID) -> Dependent_metabolites_300O

Species_abundance_300O = read_tsv("300OB/300OB_381sp_unclassified_removed_298samples.txt")
colnames(Species_abundance_300O) = as.vector(sapply(colnames(Species_abundance_300O), FUN = function(x){ if(x=="ID"){ return("ID")} else{str_split(x,"\\|")[[1]][7]} }))
Gene_family_abundance_300O = read_tsv("300OB/300OB_407_generalpath_RTabu_uc_removed_298samples.txt")

Phenotypes_300O = read_tsv("300OB/300OB_phenotypes.txt")
Covariates_300O = select(Phenotypes_300O, c("ID", "Age", "Gender","BMI"))
Covariates_300O %>% mutate(Gender = Gender-1) -> Covariates_300O

CutC_300OB = read_tsv("300OB/Total_table_cutC_300OB.tsv")
#######################
###Analysis############
#######################

Super_plot = vector("list",length(unique(colnames(Dependent_metabolites_300O)))-1)
i = 1
for (Metabolite in unique(colnames(select(Dependent_metabolites_300O, -ID)))){
  Dependent_metabolites_300O %>% select(Metabolite) -> Meta_300OB
  Dependent_metabolites %>% select(Metabolite) -> Meta_LLD
  ggplot() + geom_density(aes_(x=as_vector(Meta_LLD), fill= "red", alpha=0.5)) + geom_density(aes_(x=as_vector(Meta_300OB),fill="blue",alpha=0.5)) + theme_bw() + ggtitle(Metabolite) -> PLOT
  Super_plot[[i]] = PLOT
  i = i + 1
  #print(PLOT)
  wilcox.test(as_vector(Meta_300OB), as_vector(Meta_LLD)) -> Diff_Test
  Ratio = median(as_vector(Meta_300OB)) / median(as_vector(Meta_LLD))
  paste(c(Metabolite,"=", Diff_Test$p.value, ", Ratio(OB/LLD)=", Ratio ), collapse="") -> Diff_R
  print(Diff_R)
}
wrap_plots(Super_plot,nrow = 2) + plot_layout(guides = "collect")


Input_MGS_300OB = Prepare_Metagenome_species(Species_abundance_300O, Metabolites = Dependent_metabolites_300O, Covariates = Covariates_300O)
Input_MGS = Prepare_Metagenome_species(Species_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)


#Age, Sex, BMI associations
Metabolites = c("TMAO", "Choline", "Betaine", "L-Carnitine", "y-butyrobetaine")
apply(select(Input_MGS[[1]],-ID), 2, FUN= function(x){ Y = Rank_norm(x) ; mutate(Covariate,Dep=Y) ; summary(lm(Dep ~ Age + Gender + BMI )) -> Results ; Results$coefficients -> Results ; return(c(Results$Estimate, Results$`Pr(>|t|)`)) }) -> Result_LLD ; apply(select(Input_MGS_300OB[[1]],-ID), 2, FUN= function(x){ Y = Rank_norm(x) ; summary(lm(Y ~ Reg_OB +Covs1_300+Covs2_300 )) -> Results ; as.data.frame(Results$coefficients)["Reg",] -> Results ; return(c(Results$Estimate, Results$`Pr(>|t|)`)) }) -> Result_300

Covariates_results = tibble()
for (Covariate in c("Age", "Gender", "BMI")){
  setdiff(c("Age", "Gender", "BMI"), Covariate) -> Covs
  Input_MGS[[2]] %>% select(Covs[1]) %>% as_vector() -> Covs1_LL ; Input_MGS[[2]] %>% select(Covs[2]) %>% as_vector() -> Covs2_LL ;Input_MGS_300OB[[2]] %>% select(Covs[1]) %>% as_vector() -> Covs1_300 ; Input_MGS_300OB[[2]] %>% select(Covs[2]) %>% as_vector() -> Covs2_300
  
  Input_MGS[[2]] %>% select(Covariate) %>% as_vector() -> Reg
  Input_MGS_300OB[[2]] %>% select(Covariate) %>% as_vector() -> Reg_OB
  apply(select(Input_MGS[[1]],-ID), 2, FUN= function(x){ Y = Rank_norm(x) ; summary(lm(Y ~ Reg +Covs1_LL+Covs2_LL )) -> Results ; as.data.frame(Results$coefficients)["Reg",] -> Results ; return(c(Results$Estimate, Results$`Pr(>|t|)`)) }) -> Result_LLD ; apply(select(Input_MGS_300OB[[1]],-ID), 2, FUN= function(x){ Y = Rank_norm(x) ; summary(lm(Y ~ Reg_OB +Covs1_300+Covs2_300 )) -> Results ; as.data.frame(Results$coefficients)["Reg",] -> Results ; return(c(Results$Estimate, Results$`Pr(>|t|)`)) }) -> Result_300
  as_tibble(Result_LLD) %>% select(Metabolites) %>% t()  %>% as_tibble() %>% mutate(Metabolite = Metabolites, Cohort="LLD", Covariate= Covariate) -> Result_LLD ; as_tibble(Result_300) %>% select(Metabolites) %>% t()  %>% as_tibble() %>% mutate(Metabolite = Metabolites, Cohort="300OB",Covariate= Covariate) -> Result_300
  rbind(Covariates_results, Result_300, Result_LLD) -> Covariates_results
}

Transformation_composition(Input_MGS[[3]])



Filter_abundance(Input_MGS_300OB[[3]],Input_MGS[[3]]) -> Species_abundance_twodata
Input_MGS[[3]] = Species_abundance_twodata[[2]]
Input_MGS_300OB[[3]] = Species_abundance_twodata[[1]]

#LLD species associations
#Results with BMI as covariate
Results_MGS = Iterate_Metabolites(Input_MGS[[1]],Input_MGS[[2]],Transformation_composition(Input_MGS[[3]]), "MGS_species",BMI = T, Correct = "BH")
#Results_MGS_noBMI = Iterate_Metabolites(Input_MGS[[1]],Input_MGS[[2]],Input_MGS[[3]], "MGS_species", BMI = F)

Results_MGS  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_MGS_genus

#300OB species associations
Results_MGS_300OB = Iterate_Metabolites(Input_MGS_300OB[[1]],Input_MGS_300OB[[2]],Transformation_composition(Input_MGS_300OB[[3]]), "MGS_300OB", Correct="BH")
Results_MGS_300OB  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_MGS_genus


full_join(Results_MGS,  Results_MGS_300OB, by= c("Regressor", "Metabolite"), suffix = c("LLD","300OB")) -> Combined_species
write_csv(x = Combined_species, path = "LLD_300OB/Summary_statistics_species.csv")


##Metanalysis
Metanalyze = function(Combined_species){
  library(meta)

  Meta_stats = tibble()
  for (BUG in unique(Combined_species$Regressor)){
    Bug = filter(Combined_species, Regressor == BUG) %>% drop_na()
    for (Metabolite_i in Bug$Metabolite){
      Bug %>% filter(Metabolite == Metabolite_i) -> Entry
      Input_data = tibble(Beta = as.numeric(c(Entry$BetaLLD, Entry$Beta300OB)), SE=as.numeric(c(Entry$SELLD, Entry$SE300OB)), P=c(Entry$PvalueLLD,Entry$Pvalue300OB), N =c(Entry$NLLD, Entry$N300OB)) 
      
      meta::metagen(TE=Beta, seTE=SE,pval=P, data=Input_data, comb.fixed = T, comb.random= T) -> meta_value
      
      Agreement = ifelse(sign(as.numeric(Entry$BetaLLD)) == sign(as.numeric(Entry$Beta300OB)), "Agree", "Desagree")
      Sub_result = tibble(Bug = BUG, Metabolite = Metabolite_i, Concordance = Agreement, beta1= Entry$BetaLLD, beta2 =Entry$Beta300OB, MetaP=meta_value$pval.fixed, MetaW=meta_value$TE.fixed , Meta_random_P=meta_value$pval.random, Meta_random_W=meta_value$TE.random)
      Meta_stats = rbind(Meta_stats, Sub_result)
    }  
  }
  return(Meta_stats)
}

Metanalyze(Combined_species) ->  Meta_stats 
Meta_stats %>% mutate(FDR = p.adjust(MetaP, "fdr")) %>% ggplot(aes(x=MetaW, y = -log10(MetaP), col=FDR<0.05)) + geom_point() + theme_bw() + facet_wrap(~Metabolite)

##



Input_pathway_LLD = Prepare_Pathways(Gene_family_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
Input_pathway_300OB = Prepare_Pathways(Gene_family_abundance_300O, Metabolites = Dependent_metabolites_300O, Covariates = Covariates_300O)
Filter_abundance(Input_pathway_300OB[[3]],Input_pathway_LLD[[3]]) -> pathway_abundance_twodata
Input_pathway_LLD[[3]] = pathway_abundance_twodata[[2]]
Input_pathway_300OB[[3]] = pathway_abundance_twodata[[1]]


###Pathways LLD
Results_pathway_LLD = Iterate_Metabolites(Input_pathway_LLD[[1]],Input_pathway_LLD[[2]],Transformation_composition(Input_pathway_LLD[[3]]), "Results_pathway_LLD",BMI = T,Correct = "BH")
Results_pathway_LLD  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

###Pathways 300OB
Results_pathway_300OB = Iterate_Metabolites(Input_pathway_300OB[[1]],Input_pathway_300OB[[2]],Transformation_composition(Input_pathway_300OB[[3]]), "Input_pathway_300OB",BMI = T, Correct="BH")
Results_pathway_300OB  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

full_join(Results_pathway_LLD,  Results_pathway_300OB, by= c("Regressor", "Metabolite"), suffix = c("LLD","300OB")) -> Combined_pathway
write_csv(x = Combined_pathway, path = "LLD_300OB/Summary_statistics_pathway.csv")

Metanalyze(Combined_pathway) ->  Meta_pathway 
Meta_pathway %>% mutate(FDR = p.adjust(MetaP, "fdr")) %>% ggplot(aes(x=MetaW, y = -log10(MetaP), col=FDR<0.05)) + geom_point() + theme_bw() + facet_wrap(~Metabolite)



####Phenotypic associations

###LLD

Clinical_Questionaries %>% select(-c(antrop_gender.F1M2, antrop_age, antrop_BMI)) -> Clinical_Questionaries
List_output = Match_dataset(Dependent_metabolites, Covariates, Clinical_Questionaries)
Results_Phenos_LLD = Iterate_Metabolites(List_output[[1]],List_output[[2]],List_output[[3]], "Phenotypes_LLD",Correct = "BH")
Results_Phenos_LLD %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.05)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

write_csv(x = Results_Phenos_LLD, path = "LLD_300OB/Phenotypes_LLD.csv")


###300OB


Phenotypes_300O = select(Phenotypes_300O, -c("Age", "Gender","BMI"))
List_output = Match_dataset(Dependent_metabolites_300O, Covariates_300O, Phenotypes_300O)
Results_Phenos_300OB = Iterate_Metabolites(List_output[[1]],List_output[[2]],List_output[[3]], "Phenotypes_300OB",Correct = "BH")
Results_Phenos_300OB  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.05)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

write_csv(x = Results_Phenos_300OB, path = "LLD_300OB/Phenotypes_300OB.csv")


####LLD diet
Regr = Prepare_diet(Diet_regressors, Metabolites = Dependent_metabolites, Covariates = Covariates)
Results_Diet = Iterate_Metabolites(Regr[[1]],Regr[[2]],Regr[[3]], "Diet",Correct="BH")
Results_Diet  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.05)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
write_csv(x = Results_Diet, path = "LLD_300OB/Diet_LLD.csv")


###LLD: Diet and Species
Transformation_composition(Input_MGS[[3]]) -> LLD_species_abundance
Fit_RandomForest(select(LLD_species_abundance, filter(Meta_stats, MetaP < 0.2)$Bug), Input_MGS[[1]]$TMAO) -> RF_species 

for (Metabolite in c("TMAO", "y-butyrobetaine")){
  Regressor = as.vector(as_vector(select(Input_MGS[[1]], Metabolite)))
  Transformation_composition(Input_MGS[[3]]) -> LLD_species_abundance
  Meta_stats %>% filter(Metabolite == Metabolite) %>% filter(MetaP < 0.2) -> Predictors
  
  
  Fit_RandomForest(select(LLD_species_abundance, Predictors$Bug), Regressor) -> RF_species 
  
  
  ######Diet
  Regressor = as.vector(as_vector(select(Regr[[1]], Metabolite)))
  Results_Diet %>% filter(Metabolite == Metabolite) %>% filter(Pvalue < 0.2) -> Predictors
  Fit_RandomForest(select(Regr[[3]], Predictors$Regressor), Regressor) -> RF_diet
}

####Gene abundance

Prepare_cutC_shortbred = function(Cluster_abundance, Metabolites, Covariates){
  Cluster_abundance2 = Cluster_abundance
  Cluster_abundance %>% select(-c(Hits,TotMarkerLength)) %>% spread(Family, Count) -> Cluster_abundance2
  Cluster_abundance2 %>% summarise_if(is_numeric, sum) %>% t() %>% as.data.frame() %>%
  rownames_to_column() %>% filter(V1 > 0) -> Keep_columns 
  Cluster_abundance2 %>% select(c("ID",Keep_columns$rowname)) -> Cluster_abundance2
  
  apply(select(Cluster_abundance2, -ID), 1, FUN = function(x){sum(as.numeric(x))})%>% as.vector() -> SUM_Values
  Cluster_abundance2 %>% mutate(Total = SUM_Values) -> Cluster_abundance2
  List_output = Match_dataset(Metabolites, Covariates, Cluster_abundance2)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Cluster_abundance2=List_output[[3]] 
  IDs = Cluster_abundance2$ID
  Pathway_abundance = Metabolome_transformation(select(Cluster_abundance2,-ID)) #Transformation_composition(Cluster_abundance2)
  as_tibble(Pathway_abundance) %>% mutate(ID= IDs) -> Pathway_abundance
  return(list(Metabolites_meta, Covariates_meta, Pathway_abundance))
}

Prepare_cutC_shortbred(CutC_LLD, Dependent_metabolites, Covariates) -> Input_shortbred
Results_shortbred_LLD = Iterate_Metabolites(Input_shortbred[[1]], Input_shortbred[[2]], Input_shortbred[[3]], "MGS_shortbred_LLD",FDR_iter = 0)
Results_shortbred_LLD  %>% select(Metabolite, Regressor, Pvalue) %>% spread(Metabolite, Pvalue) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results2
Results_shortbred_LLD %>% mutate(Beta = as.numeric(Beta)) %>% select(Metabolite, Regressor, Beta) %>% spread(Metabolite, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results

wide_results_Annotation = wide_results2
wide_results_Annotation[wide_results2 > 0.1] <- ""
wide_results_Annotation[wide_results2 < 0.1] <- "."
wide_results_Annotation[wide_results2 < 0.05] <- "*"
wide_results_Annotation[wide_results2 < 0.01] <- "**"

pheatmap::pheatmap(wide_results,fontsize_number = 20, display_numbers = wide_results_Annotation)


Prepare_cutC_shortbred(Cluster_abundance = CutC_300OB, Metabolites =Dependent_metabolites_300O , Covariates=Covariates_300O) -> Input_shortbred


Results_shortbred = Iterate_Metabolites(Input_shortbred[[1]], Input_shortbred[[2]], Input_shortbred[[3]], "MGS_shortbred_300OB",FDR_iter = 0)
Results_shortbred  %>% select(Metabolite, Regressor, Pvalue) %>% spread(Metabolite, Pvalue) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results2
Results_shortbred %>% mutate(Beta = as.numeric(Beta)) %>% select(Metabolite, Regressor, Beta) %>% spread(Metabolite, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results

wide_results_Annotation = wide_results2
wide_results_Annotation[wide_results2 > 0.1] <- ""
wide_results_Annotation[wide_results2 < 0.1] <- "."
wide_results_Annotation[wide_results2 < 0.05] <- "*"
wide_results_Annotation[wide_results2 < 0.01] <- "**"



pheatmap::pheatmap(wide_results,fontsize_number = 20, display_numbers = wide_results_Annotation)



full_join(Results_shortbred_LLD,  Results_shortbred, by= c("Regressor", "Metabolite"), suffix = c("LLD","300OB")) -> Combined_shortbred
Combined_shortbred %>% select(-c("FDRLLD","FDR300OB")) %>% drop_na() -> Combined_shortbred
Metanalyze(Combined_shortbred) ->  Meta_shortbred
Meta_shortbred %>% mutate(FDR = p.adjust(MetaP, "fdr")) %>% ggplot(aes(x=MetaW, y = -log10(MetaP), col=FDR<0.05)) + geom_point() + theme_bw() + facet_wrap(~Metabolite)
