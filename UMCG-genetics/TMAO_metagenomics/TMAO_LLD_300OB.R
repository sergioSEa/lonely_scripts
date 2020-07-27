#TMAO correlation with microbiome

setwd("~/Documents/PhD/TMAO_project/MSS/")
library(tidyverse) #Data handling
library(microbiome) #cdr transformation source("Microbiome_function.R")
library(reshape2)
library(patchwork)
library(pheatmap)
library(glmnet)
source("Functions_TMAO.R")


################################
####################DATA########
################################

##Metaphlan3
read_tsv("MP3_LLD_300OB_ID.tsv") %>% select(-"X1172") %>% select(-c("UNKNOWN", "UNKNOWN_1")) -> All_taxa
colnames(All_taxa)[1] = "ID"
All_taxa %>% mutate(Cohort = ifelse(grepl("LLD", ID), "LLD", "300OB" )) -> All_taxa

Species_abundance = filter(All_taxa, Cohort=="LLD") %>% select(-Cohort)
Species_abundance_300O = filter(All_taxa, Cohort=="300OB") %>% select(-Cohort)

#Species_abundance2 = read_tsv("Data/LLD/LLD_MP3_250.txt") %>% filter(grepl("s__", ID))%>% select(-X257)
#Species_abundance2 %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% as_tibble() %>% `colnames<-`(c("ID", Species_abundance2$ID)) %>%
#  filter(! ID == "ID")  -> Species_abundance_MP3
#apply(select(Species_abundance_MP3,-ID), 2, FUN= function(x){ as.numeric(x) } ) %>% as_tibble() %>% mutate(ID = Species_abundance_MP3$ID) %>% `colnames<-`(c(colnames(select(Species_abundance_MP3, -ID)), "ID")) -> Species_abundance_MP3





###LLDeep data###
Dependent_metabolites = read_tsv("Data/LLD/Pheno_file.tsv")
Covariates =  read_tsv("Data/LLD/Cov_file.tsv")

Species_abundance = read_tsv("Data/LLD/data_LLD_baseline_1135samples_698species.txt")
Genus_abundance = read_tsv("Data/LLD/data_LLD_baseline_1135samples_698genus.txt")
Family_abundance = read_tsv("Data/LLD/data_LLD_baseline_1135samples_698family.txt")
Phylum_abundance = read_tsv("Data/LLD/data_LLD_baseline_1135samples_698phylum.txt")

#Gene abundance
Gene_abundance = read_tsv(file = "Data/LLD/Gene_abundance")
#Cluster_abundance = read_tsv(file = "Data/LLD/Total_table_cutC_shortbred.tsv")
Gene_family_abundance = read_tsv("Data/LLD/data_LLD_baseline_1135samples_489pathway.txt")


Clinical_Questionaries = read_tsv("Data/LLD/LLD_1135subjects_allphenotypes.txt")
Diet_start = which(colnames(Clinical_Questionaries)=="alcohol_products_log") ; Diet_end = which(colnames(Clinical_Questionaries)=="ACE_inhibitor")-1
Diet_regressors = Clinical_Questionaries[,Diet_start:Diet_end] ; Diet_regressors %>% mutate(ID = Clinical_Questionaries$ID) -> Diet_regressors
Clinical_Questionaries = Clinical_Questionaries[,1:Diet_start -1 ]

CutC_LLD = read_tsv("Data/LLD/Total_table_cutC_shortbred.tsv")

BGC_LLD = read_tsv("Data/LLD/RPKM_BiG-MAP2.tsv")
BGC_LLD_coverage = read_tsv("Data/LLD/coverage_BiG-MAP.tsv")

####300 OB data####
Dependent_metabolites_300O = read_tsv(file = "Data/300OB/330_O_metabolite_measures.tsv")
colnames(Dependent_metabolites_300O)[1] =  "ID"
Dependent_metabolites_300O %>% drop_na() -> Dependent_metabolites_300O
sapply(select(Dependent_metabolites_300O, -ID), FUN= function(x){transform(as.numeric(x))}) %>% as_tibble() %>% mutate(ID = Dependent_metabolites_300O$ID) -> Dependent_metabolites_300O

Species_abundance_300O = read_tsv("Data/300OB/300OB_381sp_unclassified_removed_298samples.txt")
colnames(Species_abundance_300O) = as.vector(sapply(colnames(Species_abundance_300O), FUN = function(x){ if(x=="ID"){ return("ID")} else{str_split(x,"\\|")[[1]][7]} }))
Gene_family_abundance_300O = read_tsv("Data/300OB/300OB_407_generalpath_RTabu_uc_removed_298samples.txt")

Phenotypes_300O = read_tsv("Data/300OB/300OB_phenotypes.txt")
Covariates_300O = select(Phenotypes_300O, c("ID", "Age", "Gender","BMI"))
Covariates_300O %>% mutate(Gender = Gender-1) -> Covariates_300O

CutC_300OB = read_tsv("Data/300OB/Total_table_cutC_300OB.tsv")


BGC_300OB = read_tsv("Data/300OB/RPKM_BiG-MAP_300OB.tsv")
BGC_300OB_coverage = read_tsv("Data/300OB/coverage_BiG-MAP_300OB.tsv")


#Add creatinine as covariate
Clinical_Questionaries %>% mutate(Creatinine = Biochem_Creatinine) %>% select(ID, Creatinine) -> creatinine_LLD
left_join(Covariates, creatinine_LLD) -> Covariates
Phenotypes_300O %>% mutate(Creatinine = kreatinin) %>% select(ID, Creatinine) -> creatinine_300OB
left_join(Covariates_300O, creatinine_300OB) -> Covariates_300O


#######################################################################################
#######################################################################################

#######################
###Analysis############
#######################

#######################
#Analyze distributions#
#######################

Dependent_metabolites_300O %>% gather(Metabolite, Measurement, 1:5, factor_key=TRUE) %>% mutate(Cohort="300OB") -> Info_distr_OB
Dependent_metabolites %>% gather(Metabolite, Measurement, 2:6, factor_key=TRUE) %>% mutate(Cohort="LLD") -> Info_distr_LLD
rbind(Info_distr_OB,Info_distr_LLD) -> Info_distr

#Plot with outliers
Info_distr %>% ggplot() + geom_density(aes(x=Measurement,fill=Cohort), alpha=0.5)+ theme_bw() + facet_wrap(~Metabolite, scales="free") +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1")) -> PLOT

#TMAO has  outlaiers, bibliography sasys levels are around 2/3 micromol/L https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6214130/#:~:text=The%20plasma%20level%20of%20TMAO,atherothrombotic%20cardiovascular%20disease%20(CVD).
##Checking all outliers, don't seem to be the same samples
Info_distr %>% group_by(Metabolite) %>% summarise(Outlier = mean(Measurement) + 3* sd(Measurement)) -> Out_info
Info_distr %>% filter(Metabolite == "TMAO") %>% mutate(Outlier = Measurement > filter(Out_info,Metabolite=="TMAO")$Outlier) %>% filter(Outlier == T) -> Outliers_TMAO
Info_distr %>% filter(Metabolite == "Choline") %>% mutate(Outlier = Measurement > filter(Out_info,Metabolite=="Choline")$Outlier) %>% filter(Outlier == T) -> Outliers_Choline
Info_distr %>% filter(Metabolite == "Betaine") %>% mutate(Outlier = Measurement > filter(Out_info,Metabolite=="Betaine")$Outlier) %>% filter(Outlier == T) -> Outliers_Betaine
Info_distr %>% filter(Metabolite == "L-Carnitine") %>% mutate(Outlier = Measurement > filter(Out_info,Metabolite=="L-Carnitine")$Outlier) %>% filter(Outlier == T) -> Outliers_Lcarnitine
Info_distr %>% filter(Metabolite == "y-butyrobetaine") %>% mutate(Outlier = Measurement > filter(Out_info,Metabolite=="y-butyrobetaine")$Outlier) %>% filter(Outlier == T) -> Outliers_ybutyrobetaine
dim(Outliers_TMAO) ;  dim(Outliers_Choline) ; dim(Outliers_Betaine) ; dim(Outliers_Lcarnitine) ; dim(Outliers_ybutyrobetaine)
match(Outliers_TMAO$ID, Outliers_Betaine$ID)

#Second distribution plot, removing TMAO's outliers
Info_distr %>% filter(Metabolite == "TMAO") %>% mutate(Outlaier = abs(Measurement) > (mean(Measurement) + 2* sd(Measurement))) %>% filter(Outlaier == T) %>% arrange(Measurement)  -> Outlayers_TMAO
Info_distr %>% filter(! (Metabolite == "TMAO" & ID %in% Outlayers_TMAO$ID)) %>% ggplot() + geom_density(aes(x=Measurement,fill=Cohort), alpha=0.5)+ theme_bw() + facet_wrap(~Metabolite, scales="free") +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1"))
#Summary statistics from each cohort, not removing ourliers
Info_distr %>% group_by(Cohort, Metabolite) %>% 
  summarise(Sd = sd(Measurement), Quantile0 = quantile(Measurement)[1],Quantile25 = quantile(Measurement)[2],Quantile50 = quantile(Measurement)[3],Quantile75=quantile(Measurement)[4],Quantile100=quantile(Measurement)[5])
#Test differences in MEAN and SD between both cohorts (after removing outliers)
Info_distr %>% filter(! (Metabolite == "TMAO" & ID %in% Outlayers_TMAO$ID)) -> Info_distr2
Differences = tibble()
Info_distr2 = tibble()
for (Metabol in unique(Info_distr$Metabolite)){
  #Identify outlier per cohort
  Info_distr %>% filter(Metabolite == Metabol) %>% filter(Cohort == "LLD") %>% filter(! abs(Measurement) > (mean(Measurement) + 2* sd(Measurement))) -> LLD_filter
  Info_distr %>% filter(Metabolite == Metabol) %>% filter(Cohort == "300OB") %>% filter(! abs(Measurement) > (mean(Measurement) + 2* sd(Measurement))) -> OB_filter
  
  rbind(LLD_filter, OB_filter) -> Met_filter
  rbind(Info_distr2, Met_filter) -> Info_distr2
  #Test for mean difference
  lm(Measurement ~ Cohort, Met_filter) -> Result
  #Test for sd difference
  var.test(Measurement ~ Cohort, Met_filter, alternative = "two.sided") -> Result2
  
  as.data.frame(summary(Result)$coefficients)["CohortLLD",] %>% as_tibble() %>% mutate(Metabolite = Metabol) ->Diff
  Diff %>% mutate(Ratio_variance_300OBvsLLD = Result2$estimate, Pvalues_variance = Result2$p.value) -> Diff
  Differences = rbind(Differences, Diff)
}

#Plot statistics
Differences %>% ggplot(aes(x=Estimate, y=-log10(`Pr(>|t|)`), shape=Metabolite)) + geom_point() + theme_bw() + scale_y_continuous(trans='log10') +  geom_hline(yintercept = -log10(0.05)) -> test_mean #+ scale_color_manual(values = c(wesanderson::wes_palette("Royal1"),wesanderson::wes_palette("Royal2")))
Differences %>% ggplot(aes(x=Ratio_variance_300OBvsLLD, y=-log10(Pvalues_variance), shape=Metabolite)) + geom_point() + theme_bw() + scale_y_continuous(trans='log10') +  geom_hline(yintercept = -log10(0.05)) -> test_variance
test_mean + test_variance + plot_layout(guides = "collect")

#Plot without outlaiers
Info_distr2 %>% filter(! (Metabolite == "TMAO" & ID %in% Outlayers_TMAO$ID)) %>% ggplot() + geom_density(aes(x=Measurement,fill=Cohort), alpha=0.5)+ theme_bw() + facet_wrap(~Metabolite, scales="free") +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1"))

################################
#Correlations among metabolites#
################################
cor(select(Input_MGS_300OB[[1]],-ID), method = "spear") -> OB_spear
cor(select(Input_MGS[[1]],-ID), method = "spear") -> LL_spear
#Plotting colors
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
cols <- makeColorRampPalette(wesanderson::wes_palette("Royal1",type = "continuous"),0,5)

#Make heatmaps in same order
pheatmap(LL_spear, display_numbers=round(LL_spear,2), fontsize_number = 7, fontsize_col = 7,color = cols) -> Heatmap1 
row.order = Heatmap1$tree_row$order ; col.order = Heatmap1$tree_col$order

values.to.show = OB_spear[c(1,2,3,5,4), c(1,2,3,5,4)]
pheatmap(values.to.show, display_numbers=round(values.to.show,2), fontsize_number = 7, fontsize_col = 7,color = cols, cluster_rows=F, cluster_cols=F) -> Heatmap2

Heatmap1 
Heatmap2 



####Prepare Data###
Input_MGS_300OB = Prepare_Metagenome_species(Species_abundance_300O, Metabolites = Dependent_metabolites_300O, Covariates = Covariates_300O)
Input_MGS = Prepare_Metagenome_species(Species_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
###################

########################################
#Age, Sex, BMI, Creatinine associations#
########################################
#Function
Metabolites = c("TMAO", "Choline", "Betaine", "L-Carnitine", "y-butyrobetaine")
Covariate_association = function(x){
  if (length(x) == dim(Input_MGS[[2]])[[1]]){ INPUT = Input_MGS[[2]]
  }else{ INPUT = Input_MGS_300OB[[2]] }
  Y = Rank_norm(x)
  mutate(INPUT,Dep=Y) -> Input_model
  summary(lm(Dep ~ Age + Gender + BMI + Creatinine, Input_model)) -> Results
  as.data.frame(Results$coefficients) -> Results
  c(Results$Estimate, Results$`Pr(>|t|)`) -> Results_sum
  names(Results_sum) = c("Estimate", "Pvalue")
  return(Results_sum)
}
#LLD
apply(select(Input_MGS[[1]],-ID), 2, FUN= Covariate_association) -> Result_LLD
as_tibble(Result_LLD) %>% mutate(Regressor = c("Inter_B", "Age_B", "Gender_B", "BMI_B","Creatinine_B", "Inter_P", "Age_P", "Gender_P", "BMI_P","Creatinine_P")) -> Result_LLD
#300B
apply(select(Input_MGS_300OB[[1]],-ID), 2, FUN= Covariate_association) -> Result_300OB
as_tibble(Result_300OB) %>% mutate(Regressor = c("Inter_B", "Age_B", "Gender_B", "BMI_B","Creatinine_B", "Inter_P", "Age_P", "Gender_P", "BMI_P","Creatinine_P")) -> Result_300OB
rbind(mutate(Result_LLD, Cohort = "LLD"), mutate(Result_300OB, Cohort = "300OB")) -> Combined_covariates
#Meta-analysis
Sub_result = tibble()
for (Met in Metabolites){
  Combined_covariates %>% select(c(Met, Regressor, Cohort)) -> Metabolite_meta
  Metabolite_meta %>% spread(Regressor,Met) -> Metabolite_meta
  meta::metagen(TE= Metabolite_meta$Age_B, pval=Age_P, data= Metabolite_meta,comb.fixed = T, comb.random= T ) -> meta_value
  Sub_result = rbind(Sub_result,tibble(Met = Met, Covariate = "Age", MetaP=meta_value$pval.fixed, MetaW=meta_value$TE.fixed , Meta_random_P=meta_value$pval.random, Meta_random_W=meta_value$TE.random))
  meta::metagen(TE= Metabolite_meta$BMI_B, pval=BMI_P, data= Metabolite_meta,comb.fixed = T, comb.random= T ) -> meta_value
  Sub_result =rbind(Sub_result,tibble(Met = Met, Covariate = "BMI", MetaP=meta_value$pval.fixed, MetaW=meta_value$TE.fixed , Meta_random_P=meta_value$pval.random, Meta_random_W=meta_value$TE.random))
  meta::metagen(TE= Metabolite_meta$Gender_B, pval=Gender_P, data= Metabolite_meta,comb.fixed = T, comb.random= T ) -> meta_value
  Sub_result =rbind(Sub_result,tibble(Met = Met, Covariate = "Gender", MetaP=meta_value$pval.fixed, MetaW=meta_value$TE.fixed , Meta_random_P=meta_value$pval.random, Meta_random_W=meta_value$TE.random))
  meta::metagen(TE= Metabolite_meta$Creatinine_B, pval=Creatinine_P, data= Metabolite_meta,comb.fixed = T, comb.random= T ) -> meta_value
  Sub_result =rbind(Sub_result,tibble(Met = Met, Covariate = "Creatinine", MetaP=meta_value$pval.fixed, MetaW=meta_value$TE.fixed , Meta_random_P=meta_value$pval.random, Meta_random_W=meta_value$TE.random))

}
#Plot
Sub_result %>% 
  ggplot() + geom_bar(aes(x=Met , y= -log10(MetaP), fill=MetaW>0), stat="identity")+ facet_wrap(~Covariate) + theme_bw() + scale_fill_manual(values = wesanderson::wes_palette("Royal1")) +
  coord_flip() + geom_hline(yintercept = -log10(0.05))
  



###########################
##Microbiome associations##
###########################

#1. Test for effects of the metabolite levels in the overall composition of Beta diversity
Beta_diversity = function(Species, Cov, Metabolite){
  Cov %>% mutate(Metabolite = Metabolite) -> Beta_input
  vegan::vegdist(select(Species, -ID), "bray") -> Distan
  vegan::adonis2(formula = Distan ~ Age + BMI + Gender + Creatinine +  Metabolite, data= Beta_input ) -> adonis_results
  print(adonis_results)
  return(c(adonis_results$R2[5], adonis_results$`Pr(>F)`[5]))
}
Total_beta = tibble()
for (Metabolite in colnames(select(Input_MGS[[1]],-ID))){
  print(Metabolite)
  Input_MGS[[2]] %>% drop_na() -> Cov
  Input_MGS[[1]] %>% filter(ID %in% Cov$ID) %>% select(Metabolite) %>% as_vector() -> Met
  Beta_diversity(Species = filter(Input_MGS[[3]], ID %in% Cov$ID) , Cov = Cov, Metabolite = Met) -> Beta_result
  Beta_result = tibble(R2 = Beta_result[[1]], Pvalue=Beta_result[[2]], Cohort="LLD", Metabolite = Metabolite)
  Total_beta = rbind(Total_beta, Beta_result)
  Input_MGS_300OB[[2]] %>% drop_na() -> Cov
  Input_MGS_300OB[[1]] %>% filter(ID %in% Cov$ID) %>% select(Metabolite) %>% as_vector() -> Met
  Beta_diversity(Species = filter(Input_MGS_300OB[[3]], ID %in% Cov$ID) , Cov = Cov, Metabolite = Met) -> Beta_result
  Beta_result = tibble(R2 = Beta_result[[1]], Pvalue=Beta_result[[2]], Cohort="300OB", Metabolite = Metabolite)
  Total_beta = rbind(Total_beta, Beta_result)
}

Total_beta %>% ggplot(aes(x=Metabolite, y=-log10(Pvalue), fill=Cohort)) + geom_bar(stat="identity",position = "dodge")  + coord_flip() + 
  theme_bw() +  scale_fill_manual(values = c(wesanderson::wes_palette("Royal1"), wesanderson::wes_palette("Royal2"))) + scale_color_manual(values = c("white","black")) + geom_hline(yintercept = -log10(0.05))

#2. Calculation of Alpha diversity before filtering species. This will be added back after filtering
vegan::diversity(select(Input_MGS[[3]],-ID),index = "shannon") -> diversity_LLD
vegan::diversity(select(Input_MGS_300OB[[3]],-ID),index = "shannon") -> diversity_300B

###Prepare abundance data
Filter_02 = function(df){
  apply(select(df, -ID), 2, FUN = function(x){ length(x[x == 0])/length(x) < 0.8 }) -> Keep
  names(Keep[Keep == T]) -> Keep
  
  #df %>% select(Keep) %>%  apply(2, FUN= sum) -> v_total
  #(select(df, Keep) / v_total) %>% as_tibble() %>% mutate(ID = df$ID) %>%
  #  gather(Taxa, Abundance, 1:(ncol(select(df,Keep))), factor_key=TRUE) %>%
  #  group_by(Taxa) %>% summarise(M=median(Abundance)) %>% filter(M > 0.0005) -> Keep2
  #Keep = Keep2$Taxa
  
  return(Keep)
}
#Filter LLD per prevalence
Keep_LLD = Filter_02(Input_MGS[[3]])

Transformation_composition(Input_MGS[[3]]) %>% select(c(ID, Keep_LLD)) ->Input_MGS[[3]]

Keep_300OB = Filter_02(Input_MGS_300OB[[3]])
Transformation_composition(Input_MGS_300OB[[3]]) %>% select(c(ID, Keep_300OB)) -> Input_MGS_300OB[[3]]

Input_MGS[[3]] %>% select(one_of(colnames(Input_MGS_300OB[[3]]))) -> Input_MGS[[3]]
Input_MGS_300OB[[3]] %>% select(one_of(colnames(Input_MGS[[3]]))) -> Input_MGS_300OB[[3]]

Input_MGS[[3]] %>% mutate(Shannon = diversity_LLD) -> Input_MGS[[3]]
Input_MGS_300OB[[3]] %>% mutate(Shannon = diversity_300B) -> Input_MGS_300OB[[3]]
dim(Input_MGS[[3]]) ; dim(Input_MGS_300OB[[3]])

#LLD species associations
#Results with BMI as covariate
Results_MGS = Iterate_Metabolites(Input_MGS[[1]],Input_MGS[[2]], Input_MGS[[3]], "MGS_species",BMI = T, Correct = "BH")
#Results_MGS_noBMI = Iterate_Metabolites(Input_MGS[[1]],Input_MGS[[2]],Input_MGS[[3]], "MGS_species", BMI = F)

Results_MGS  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_MGS_genus

#300OB species associations
Results_MGS_300OB = Iterate_Metabolites(Input_MGS_300OB[[1]],Input_MGS_300OB[[2]], Input_MGS_300OB[[3]], "MGS_300OB", Correct="BH")
Results_MGS_300OB  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.15)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> Fig_MGS_genus


full_join(Results_MGS,  Results_MGS_300OB, by= c("Regressor", "Metabolite"), suffix = c("LLD","300OB")) -> Combined_species
write_csv(x = Combined_species, path = "Summary_statistics_species.csv")


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
Meta_stats %>% mutate(FDR = p.adjust(MetaP, "fdr")) %>% arrange(MetaP) -> Meta_stats
Meta_stats %>% ggplot(aes(x=MetaW, y = -log10(MetaP), col=FDR<0.05)) + geom_point() + theme_bw() + facet_wrap(~Metabolite)

#####Pathwayws: no include

Input_pathway_LLD = Prepare_Pathways(Gene_family_abundance, Metabolites = Dependent_metabolites, Covariates = Covariates)
Input_pathway_300OB = Prepare_Pathways(Gene_family_abundance_300O, Metabolites = Dependent_metabolites_300O, Covariates = Covariates_300O)

Keep_LLD = Filter_02(Input_pathway_LLD[[3]])
Input_pathway_LLD[[3]] %>% select(c(ID, Keep_LLD)) ->Input_pathway_LLD[[3]]
Keep_300OB = Filter_02(Input_pathway_300OB[[3]])
Input_pathway_300OB[[3]] %>% select(c(ID, Keep_300OB)) -> Input_pathway_300OB[[3]]

Input_pathway_LLD[[3]] %>% select(one_of(colnames(Input_pathway_300OB[[3]]))) -> Input_pathway_LLD[[3]]
Input_pathway_300OB[[3]] %>% select(one_of(colnames(Input_pathway_LLD[[3]]))) -> Input_pathway_300OB[[3]] 

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

##############


################################
####Biosynthetic gene clusters#
###############################

Get_BGC_ready =  function(BGC){
  BGC %>% t() %>% as.data.frame() %>% rownames_to_column("gene_clusters") %>% as_tibble() %>% `colnames<-`(c("ID", BGC$gene_clusters)) %>%
    filter(! ID == "gene_clusters") -> BGC
  apply(select(BGC,-ID), 2, FUN= function(x){ as.numeric(x) } ) %>% as_tibble() %>% mutate(ID = BGC$ID) %>% `colnames<-`(c(colnames(select(BGC, -ID)), "ID")) -> BGC
  return(BGC)
}

#LLD
BGC_LLD = Get_BGC_ready(BGC_LLD)
BGC_LLD_coverage = Get_BGC_ready(BGC_LLD_coverage)
Match_dataset(Dependent_metabolites, Covariates, BGC_LLD) -> BGC_match_LLD
#300OB
BGC_300OB = Get_BGC_ready(BGC_300OB)
BGC_300OB_coverage = Get_BGC_ready(BGC_300OB_coverage)
Match_dataset(Dependent_metabolites_300O, Covariates_300O, BGC_300OB) -> BGC_match_300OB


#filter on abundance
Filtering_Median_n_abundance = function(DF, Norm="clr"){
  Above_thresh = function(C){
  C = as_vector(C)
  C2 = sum(C>0)/length(C)
  C2 = C2>0.3
  if (C2 == T){
    return(names(C)[1])  
  }else{ return(NA)
  }
}  
  apply(select(DF, -ID), 2, FUN=Above_thresh) -> Keep
  Keep[!is.na(Keep)] -> Keep
  DF %>% select(c("ID", names(Keep))) -> Filtered
  #filter on median abundance
  apply(select(Filtered, -ID), 2, FUN= sum) -> v_total
  (select(Filtered, -ID) / v_total) %>% as_tibble() %>% mutate(ID = Filtered$ID) %>%
  gather(BGC, Abundance, 2:ncol(Filtered), factor_key=TRUE) %>%
    group_by(BGC) %>% summarise(M=median(Abundance)) %>% filter(M < 0.01) -> Filter_out
  #Normalization on unfiltered data
  if (Norm =="clr"){
    #Normalize comp
    Transformation_composition(DF) -> DF
  } else if(Norm=="log"){
    Pseudocount = min(select(DF,-ID)[select(DF,-ID)>0])/2
    select(DF,-ID) %>% apply(2, FUN=function(x){log10(x+Pseudocount)}) %>% as.data.frame() %>% as_tibble() %>% mutate(ID = DF$ID) ->DF
  }
  #Perform filters
  DF %>% select(c("ID", names(Keep))) -> DF
  DF %>% select(! Filter_out$BGC) -> DF
  return(DF)
}

Filtering_Median_n_abundance(BGC_match_LLD[[3]]) -> BGC_match_LLD[[3]]
Filtering_Median_n_abundance(BGC_match_300OB[[3]]) -> BGC_match_300OB[[3]]

BGC_match_300OB[[3]] %>% select(one_of(colnames(BGC_match_LLD[[3]]))) -> BGC_match_300OB[[3]]
BGC_match_LLD[[3]] %>% select(one_of(colnames(BGC_match_300OB[[3]]))) -> BGC_match_LLD[[3]]
                                                                                                                                    
#LLD association
Iterate_Metabolites(BGC_match_LLD[[1]], BGC_match_LLD[[2]], BGC_match_LLD[[3]], "BGC_LLD",Correct = "BH") -> results_BGCs
results_BGCs %>% mutate(Beta = as.vector(Beta), Pvalue = as.vector(Pvalue)) ->results_BGCs
results_BGCs %>% filter(FDR < 0.05) -> Significant_BGCs
BGC_LLD_coverage %>% select(c(ID, Significant_BGCs$Regressor)) %>% gather(BGC, Coverage, 2:(1+length(unique(Significant_BGCs$Regressor))), factor_key=TRUE) %>%
  group_by(BGC) %>% summarise(median(Coverage))
for (BGC_n in 1:nrow(filter(results_BGCs, FDR<0.1))){
  results_BGCs %>% arrange(FDR) %>% filter(FDR < 0.1) -> BGC_r
  BGC_r[BGC_n,] -> Row_association
  BGC_match[[3]] %>% select(Row_association$Regressor) %>% as_vector() -> X
  BGC_match[[1]] %>% select(Row_association$Metabolite) %>% as_vector() -> Y
  cor.test(X,Y,method = "spear") -> CORR
  print(CORR)
  ggplot() + geom_point(aes(x=X, y=log(Y))) + theme_bw() + xlab(Row_association$Regressor) + ylab(Row_association$Metabolite) -> PLOT
  print(PLOT)
}

#300OB association
Iterate_Metabolites(BGC_match_300OB[[1]], BGC_match_300OB[[2]], BGC_match_300OB[[3]], "BGC_300OB",Correct = "BH") -> results_BGCs_300OB
results_BGCs_300OB %>% mutate(Beta = as.vector(Beta), Pvalue = as.vector(Pvalue)) -> results_BGCs_300OB
results_BGCs_300OB %>% filter(FDR < 0.05) -> Significant_BGCs_300OB
BGC_300OB_coverage %>% select(c(ID, Significant_BGCs_300OB$Regressor)) %>% gather(BGC, Coverage, 2:(1+length(unique(Significant_BGCs$Regressor))), factor_key=TRUE) %>%
  group_by(BGC) %>% summarise(median(Coverage))


full_join(results_BGCs,  results_BGCs_300OB, by= c("Regressor", "Metabolite"), suffix = c("LLD","300OB")) -> Combined_species
Metanalyze(Combined_species) %>% mutate(FDR = p.adjust(MetaP,"fdr")) %>% arrange(MetaP) -> Meta_BGC



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



#######################################
####Non-microbiome correlations########
#######################################

#################
####Phenotypes###
#################

###LLD

Clinical_Questionaries %>% select(-c(antrop_gender.F1M2, antrop_age, antrop_BMI,Biochem_Creatinine)) -> Clinical_Questionaries
#Overview of phenotypes
colnames(select(Clinical_Questionaries, -ID))
List_output = Match_dataset(Dependent_metabolites, Covariates, Clinical_Questionaries)
Results_Phenos_LLD = Iterate_Metabolites(List_output[[1]],List_output[[2]],List_output[[3]], "Phenotypes_LLD",Correct = "BH")
Results_Phenos_LLD %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.05)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
write_csv(x = Results_Phenos_LLD, path = "LLD_300OB/Phenotypes_LLD.csv")

###300OB
Phenotypes_300O = select(Phenotypes_300O, -c("Age", "Gender","BMI", "kreatinin"))
colnames(select(Phenotypes_300O, -ID))
List_output = Match_dataset(Dependent_metabolites_300O, Covariates_300O, Phenotypes_300O)
Results_Phenos_300OB = Iterate_Metabolites(List_output[[1]],List_output[[2]],List_output[[3]], "Phenotypes_300OB",Correct = "BH")
Results_Phenos_300OB  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.05)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

write_csv(x = Results_Phenos_300OB, path = "LLD_300OB/Phenotypes_300OB.csv")


##########
###Diet###
##########

####LLD diet
Regr = Prepare_diet(Diet_regressors, Metabolites = Dependent_metabolites, Covariates = Covariates)
Results_Diet = Iterate_Metabolites(Regr[[1]],Regr[[2]],Regr[[3]], "Diet",Correct="BH")
Results_Diet  %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(Pvalue), col = FDR<0.05)) +
  facet_wrap(~Metabolite, scales="free") + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))
write_csv(x = Results_Diet, path = "LLD_300OB/Diet_LLD.csv")









###Diet, bacteria and TMAO
R2_calc = function(Real,pred){
  rss = sum((pred - Real)^2)
  tss = sum((Real - mean(Real))^2)
  rsq = 1 - rss/tss
  return(rsq)
}

multiomics = function(MET="TMAO", Remove_col = T, alpha=1, More_models= T, Transform=F){ 
  select_to_keep = function(Correlation_bacteria){
    V1 = as.character(Correlation_bacteria[["variable"]])
    V2 = as.character(Correlation_bacteria[["V1"]])
    To_keep %>% filter(Regressor %in% c(V1,V2)) %>% arrange(Pvalue) %>% tail(1) -> Keep
    return(Keep$Regressor)
  }
  Bacteria_composition = Input_MGS[[3]]
  #1. Transform data if not transformed
  if (Transform == T){
    Bacteria_composition = Transformation_composition(Bacteria_composition)
  }
  #2. Filter to bacteria of interest
  Results_MGS %>% filter(Metabolite == MET)  %>% filter(Pvalue<0.2) -> To_keep
  Bacteria_composition %>% select(c(ID,To_keep$Regressor)) -> Bacteria_composition
  #3. Check for correlations
  if (Remove_col == T){
    select(Bacteria_composition, -ID) %>% cor(method = "spear") %>%as.data.frame()%>% rownames_to_column("V1") %>% as_tibble() %>%
      melt(id.vars=c("V1")) %>% filter(value >0.6) %>% filter(! V1 == variable) -> Highly_correlated
    apply(Highly_correlated, 1, FUN=select_to_keep) -> To_remove
    Bacteria_composition %>% select(!unique(To_remove)) -> Bacteria_composition
  }
  #Repeat steps with Diet
  Regr = Prepare_diet(Diet_regressors, Metabolites = Dependent_metabolites, Covariates = Covariates)
  Diet_info = Regr[[3]]
  #1. Filter on Diet of interest
  Results_Diet %>% filter(Metabolite == MET)  %>% filter(Pvalue<0.2) -> To_keep
  Diet_info %>% select(c(ID,To_keep$Regressor)) -> Diet_info
  #2. Check for correlations
  if (Remove_col == T){
    select(Diet_info, -ID) %>% cor(method = "spear") %>%as.data.frame()%>% rownames_to_column("V1") %>% as_tibble() %>%
      melt(id.vars=c("V1")) %>% filter(value >0.6) %>% filter(! V1 == variable) -> Highly_correlated
    apply(Highly_correlated, 1, FUN=select_to_keep) -> To_remove
    Diet_info %>% select(!unique(To_remove)) -> Diet_info
  }  
  #Put diet and bacteria together
  Diet_info %>% filter(ID %in% Bacteria_composition$ID) %>% arrange(ID) -> Diet_info
  Bacteria_composition %>% filter(ID %in% Diet_info$ID) %>% arrange(ID) -> Bacteria_composition
  inner_join(Bacteria_composition,Diet_info, "ID") -> Regressors
  Regressors %>% drop_na() -> Regressors
  #Get Covariates and dependent
  Cov = filter(Regr[[2]], ID %in% Regressors$ID)
  inner_join(Regressors,Cov, by="ID") %>% arrange(ID) -> Regressors
  Regressors %>% drop_na() -> Regressors
  #Remove outlaiers
  Outlier_info = filter(Regr[[1]], ID %in% Regressors$ID) %>% arrange(ID) %>% select(c(MET,ID)) %>% summarise_if(is.numeric, c(sd,mean))
  Dependent = filter(Regr[[1]], ID %in% Regressors$ID) %>% arrange(ID) %>% select(MET) %>% as_vector()
  which(abs(Dependent) < (2*Outlier_info$fn1 + Outlier_info$fn2)) -> Keep
  Dependent[as.vector(Keep)] -> Dependent
  
  Regressors[Keep,]  %>% select(-ID) ->Regressors
  #######
  
  ##Lasso Model##
  print(MET)
  
  
  set.seed(45)
  cvfit = cv.glmnet(y=as.vector(Dependent), x=as.matrix(Regressors), nfolds = 10,type.measure="mse",standardize=T, alpha=alpha)
  Param_best = coef(cvfit, s = "lambda.min")
  Estimation = predict(cvfit, as.matrix(Regressors), s="lambda.min")
  r2 = R2_calc(Dependent, Estimation)
  
  
  if (More_models == F){
    Bac_regressors = Regressors
    Matches = match(colnames(select(Bacteria_composition,-ID)),colnames(Regressors))
    setdiff(seq(1:dim(Bac_regressors)[2]),  Matches) -> To_0
    Bac_regressors[,To_0] = 0 
    EstimationBac = predict(cvfit, as.matrix(Bac_regressors), s="lambda.min")
    r2Bac = R2_calc(Dependent, EstimationBac)
    
    Diet_regress = Regressors
    Diet_positions = match(colnames(select(Diet_info,-ID)),colnames(Regressors))
    setdiff(seq(1:dim(Regressors)[2]),  Diet_positions) -> To_0
    Diet_regress[,To_0] = 0 
    EstimationDiet = predict(cvfit, as.matrix(Diet_regress), s="lambda.min")
    r2Diet = R2_calc(Dependent, EstimationDiet)
  } else{  
    Bac_regressors = as.matrix(select(Regressors,colnames(select(Bacteria_composition,-ID)),))
    cvfitBac = cv.glmnet(y=as.vector(Dependent), x=Bac_regressors, nfolds = 10,type.measure="mse",standardize=T, alpha=alpha)
    Param_bestBac = coef(cvfitBac, s = "lambda.min")
    EstimationBac = predict(cvfitBac, as.matrix(Bac_regressors), s="lambda.min")
    r2Bac = R2_calc(Dependent, EstimationBac)
   
    Diet_regress = as.matrix(select(Regressors,colnames(select(Diet_info,-ID))))
    cvfitDiet = cv.glmnet(y=as.vector(Dependent), x=Diet_regress, nfolds = 10,type.measure="mse",standardize=T, alpha=alpha)
    Param_bestDiet = coef(cvfitDiet, s = "lambda.min")
    EstimationDiet = predict(cvfitDiet, as.matrix(Diet_regress), s="lambda.min")
    r2Diet = R2_calc(Dependent, EstimationDiet)
    

    names(Param_bestBac[abs(Param_bestBac[1:length(Param_bestBac),]) > 0,]) -> Params_bacteria
    if(is.null(Params_bacteria)){ rownames(Param_bestBac[abs(Param_bestBac[1:length(Param_bestBac),]) > 0,]) -> Params_bacteria }
    names(Param_bestDiet[abs(Param_bestDiet[1:length(Param_bestDiet),]) > 0,]) -> Params_diet
    if(is.null(Params_diet)){ rownames(Param_bestDiet[abs(Param_bestDiet[1:length(Param_bestDiet),]) > 0,]) -> Params_diet }
    names(Param_best[abs(Param_best[1:length(Param_best),]) > 0,]) -> Params
    if(is.null(Params)){ rownames(Param_best[abs(Param_best[1:length(Param_best),]) > 0,]) -> Params }
    
    match(Params_bacteria, Params) -> match_bacteria
    match(Params_diet, Params) -> match_diet
    print(paste(c("Params in Bacteria model: ",  length(Params_bacteria), "; Also selected in complete model: ",sum(is.na(match_bacteria) == F)), collapse=""))
    print(paste(c("Params in Diet model: ",  length(Params_diet), "; Also selected in complete model: ",sum(is.na(match_diet) == F)), collapse=""))
    }

  plot(cvfit)
  print(c(r2, r2Bac, r2Diet))

}


for (Metabolite in c("TMAO", "y-butyrobetaine", "Betaine", "Choline", "L-Carnitine")){
  multiomics(MET= Metabolite,Remove_col = F, alpha=1, More_models = T)  
}


#####Bacterial "heritability"
vegan::vegdist(dplyr::select(Input_MGS[[3]], -ID), method = "jaccard") -> Microbial_kinship
Input_MGS[[2]] %>% mutate(Y = as_vector(dplyr::select(Input_MGS[[1]], "TMAO"))) -> Input_model
ans.ADE <- sommer::mmer(Y~Age + Gender + BMI,
                        random=~vs(ID,Microbial_kinship),
                        data=Input_model)
