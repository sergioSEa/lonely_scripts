#Association of 16S data with TMAO-related metabolits

###Package load
library(tidyverse) #Standard packages for data load (readr), visualization (ggplot) and processing (dplyr)
library(microbiome) #Library for clr
library(patchwork)
######For Cluster version, get file path from command line
#Input = commandArgs(trailingOnly = TRUE)
#Count_table = options[1] #Output Step 1.3. Generate summary statistics from https://github.com/alexa-kur/miQTL_cookbook
#Linking_table = options[2] #Column 1: Cohort ID, Column2: Sequencing ID ; No headers
#Covariates = options[3] #Metadata should contain the following columns: ID, Sex, BMI, Age.
#Phenos = option[4] #Should contain a column called ID, the rest will be used as dependent variables in the model

######Running from desktop
setwd("~/PhD/WORK/TMAO_16S")
Count_table =  "LLD_rarefDefault.taxonomyTable.txt.gz"
Linking_table = "coupling_pheno_16s.txt"
Covariates = "Cov_file.tsv"
Phenos = "Pheno_file.tsv"
##########

Count_table = read_tsv(Count_table) #Rows are individuals, columns are Bacteria
Linking_table = read_tsv(Linking_table, col_names=F) 
Covariates = read_tsv(Covariates)
Phenos = read_tsv(Phenos)

#Add ratios
Phenos %>% mutate(TMAO.Choline = TMAO/Choline , TMAO.Betaine = TMAO/Betaine , 
                  TMAO.Butyrobetaine = TMAO/`y-butyrobetaine`, TMAO.Carnitine = TMAO/`L-Carnitine`,
                  Butyrobetain.Carnitine =  `y-butyrobetaine`/`L-Carnitine` ) -> Phenos

#Functions
Choose_replicate = function(Linking_table){
  #Random choice if Same individual's microbiome has been sequenced more than once
  to_remove = vector()
  for (Entry in unique(Linking_table$X1[duplicated(Linking_table$X1)])){
    Linking_table %>% filter(X1 == Entry) %>% sample_n(1) -> Choice
    Linking_table %>% filter(X1 == Entry) %>% filter(!X2 == Choice$X2) -> Out
    to_remove = c(to_remove, Out$X2)
  }
  Linking_table %>% filter(!X2 %in% to_remove) -> Filtered_table
  return(Filtered_table)
}

Filter_unclassified = function(Count_table){
  #NOTAX is the notation for the classifier unclassified samples
  Remove_columns = colnames(Count_table)[grepl("NOTAX", colnames(Count_table))]
  Count_table %>% select(-Remove_columns) %>% select(-rootrank.Root) -> Count_table
  return(Count_table)
}

Transformation = function(Count_table){
  SampleID = Count_table$SampleID
  Count_table %>% select(-SampleID) -> Counts
  #Transform counts
  
  colnames(Counts) -> N
  Counts %>% t() -> Check
  Counts_transformed = as_tibble(t(abundances(x=as.data.frame(Check), transform="clr")))
  ###
  Counts_transformed %>% mutate(SampleID = SampleID) -> Counts_transformed
  return(Counts_transformed)
}
Filter_by_abundance = function(Count_table, threshold=10){
  SampleID = Count_table$SampleID
  Count_table %>% select(-SampleID) -> abundance_matrix
  abundance_matrix = abundance_matrix[,((colSums(abundance_matrix !=0) / nrow(abundance_matrix)) *100 )>threshold]
  abundance_matrix %>% mutate(SampleID = SampleID) -> abundance_matrix
  return(colnames(abundance_matrix))
  #return(abundance_matrix)
}

Fit_model = function(Count_table, Phenotype, Covariates, BMI_model=T){
    Total_results = tibble()
    
    for (Number in seq(1:dim(Count_table)[2])){
      Bug = colnames(Count_table)[Number]
      Bug_abundance = as.vector(as_vector(Count_table[,Number]))
      for (Number_metabolite in seq(1:(dim(Phenotype)[2]-1))){
          Number_metabolite = Number_metabolite + 1
          Metabolite = colnames(Phenotype)[Number_metabolite]
          Metabolite_measurement = as.vector(as_vector(Phenotype[,Number_metabolite]))
          Metabolite_measurement = qnorm((rank(Metabolite_measurement,na.last="keep")-0.5)/sum(!is.na(Metabolite_measurement)))
          ###Fit model
          Covariates %>% mutate(Metabolite_measurement = Metabolite_measurement, Bug_abundance = Bug_abundance ) -> Model_data
          if (length(unique(Model_data$Gender) > 1)){
              if (BMI_model == T){
                MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance + Age + Gender + BMI, data = Model_data)
              }else{ MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance + Age + Gender, data = Model_data) }
              Summary = as_tibble(summary(MODEL_FIT)$coefficients)
              Beta = Summary$Estimate[2]
              pvalue = Summary$`Pr(>|t|)`[2]
              stand = Summary$`Std. Error`[2]
              N_samples = length(MODEL_FIT$residuals)
          }else{
            if (BMI_model == T){
              MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance + Age + BMI, data = Model_data)
            }else{ MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance + Age, data = Model_data) }
            Summary = as_tibble(summary(MODEL_FIT)$coefficients)
            Beta = Summary$Estimate[2]
            pvalue = Summary$`Pr(>|t|)`[2]
            stand = Summary$`Std. Error`[2]
            N_samples = length(MODEL_FIT$residuals)
          }
          Results = as_tibble(t(c(Bug, Metabolite, Beta, stand, pvalue, N_samples)))
          colnames(Results) = c("Bug", "Metabolite", "Beta", "SE", "Pvalue", "N_samples")
          Total_results = rbind(Total_results,Results)    
      }
   
    }
    return(Total_results)  
}

Calculate_fdr = function(Count_table, Phenotype, Covariates, Real_pvalues, BMI_model=T){
  FDR_vector = vector()
  #1st, Number of permutations
  Random_distribution  = vector()
  for (i in seq(100)){
    #Per each permutation rearrange data and fit model
    #rearrange
    sample_n(Count_table, size= nrow(Count_table)) -> Rearrange
    #Fitmodel
    Fit_model(Rearrange, Phenotype, Covariates,BMI_model) -> H0
    #Save Pvalues in the random distribution of pvalues
    Random_distribution = c(Random_distribution,as.numeric(H0$Pvalue))
  }
  #2nd, calculate FDR per each threshold
  Real_pvalues = Real_pvalues$Pvalue 
  for (Threshold in Real_pvalues){
    Threshold = as.numeric(Threshold)
    P = sum((Real_pvalues <= Threshold) *1)
    FP = sum((Random_distribution <= Threshold)*1)
    FDR = (FP/100)/P
    if (FDR >1){ FDR = 1}
    #If the FDR value of a higher Pvalue is lower, make the FDR of the lower Pvalue at least as low
    if (length(FDR_vector) > 0){ if (FDR > FDR_vector[length(FDR_vector)]){ FDR = FDR_vector[length(FDR_vector)]} }
    FDR_vector = c(FDR_vector,FDR)
  }
  return(FDR_vector)
}


###QTL Analysis
Count_table = Filter_unclassified(Count_table)
Linking_table = Choose_replicate(Linking_table)
taxonomic_levels = c("genus","order", "phylum","class","family")
set.seed(111)

for (Level in taxonomic_levels){
  print(paste0(c("Working in",Level, "taxonomic level"),collapse=" "))
  
  Count_table$SampleID -> Sample_ID
  colnames(Count_table)[grepl(Level, colnames(Count_table))] -> Level_filter
  Count_table %>% select(Level_filter) %>% mutate(SampleID =  Sample_ID) -> Taxa_table
  Taxa_table_names_keep = Filter_by_abundance(Taxa_table)
  
  #Merging SeqID and cohort ID
  Linking_table %>% arrange(X2) %>% filter(X2 %in% Taxa_table$SampleID) -> Linking_table_taxa
  Taxa_table %>% arrange(SampleID) %>% filter(SampleID %in% Linking_table$X2) %>% 
    mutate(SampleID = Linking_table_taxa$X1) %>% arrange(SampleID) -> Taxa_table
  
  #Get same number of entries in Covariates and Phenotyes than in the 16S data
  Covariates %>% arrange(ID) %>% filter(ID %in% Taxa_table$SampleID) -> Covariates
  Phenos %>% arrange(ID) %>% filter(ID %in% Taxa_table$SampleID) -> Phenos
  
  Taxa_table %>% filter(SampleID %in% Phenos$ID) -> Taxa_table
  
  
  #Fit each of the models: Females, Males  or All
  for (S in c(unique(Covariates$Gender), 3) ){
    if (is.na(S)){next}
    print(paste0(c("Working in",S, "Gender group"),collapse=" "))
    
    if (S == 3){ Covariates -> Cov_model
    } else{ Covariates %>% filter(Gender == S) -> Cov_model }
    
    Phenos %>% filter(ID %in% Cov_model$ID) -> Pheno_model
    Taxa_table %>% filter(SampleID %in% Cov_model$ID) -> Taxa_model
    
    Transformation(Taxa_model) -> Taxa_model
    select(Taxa_model,Taxa_table_names_keep) -> Taxa_model
    Taxa_model %>% select(-SampleID) -> Taxa_model
    
    for (BMI_model in c(T,F)){
      if (BMI_model == T){
        Output_name = paste(c("16S",Level,S, "tsv"), collapse= ".")
      }else{Output_name = paste(c("16S",Level,S,"noBMI","tsv"), collapse= ".")}
      
      Result_table = Fit_model(Taxa_model, Pheno_model, Cov_model, BMI_model)
  
      Result_table %>% arrange(desc(as.numeric(Pvalue))) -> Result_table
      Calculate_fdr(Taxa_model, Pheno_model, Cov_model, Result_table, BMI_model) -> FDR_vector
      Result_table %>% mutate(FDR = FDR_vector) -> Result_table
      Result_table %>% arrange(FDR) -> Result_table
      write_tsv(Result_table,path = Output_name)
    }
    #p.adjust(Result_table$Pvalue ,"fdr")
  }

  
}

#Put all Results together
#1. BMI
file.names <- dir(".", pattern ="16S.*.tsv")
All_files = tibble()
for(i in 1:length(file.names)){
  if (file.names[i] %in% c("TMAOmetabolites_16S.R", "Results_16S.tsv")){next}
  if (grepl("noBMI",file.names[i])){ next }
  print(file.names[i])
  file <- read_tsv(file.names[i])
  file %>% mutate(Source=file.names[i]) -> O
  if (dim(O)[1] >= 1){ All_files = rbind(All_files, O) }
}
write_tsv(All_files, "Results_16S.tsv") #ICNLUDES BMI

#2. No BMI
file.names <- dir(".", pattern ="16S.*noBMI.tsv")
All_files = tibble()
for(i in 1:length(file.names)){
  if (file.names[i] %in% c("TMAOmetabolites_16S.R", "Results_16S.tsv")){next}
  print(file.names[i])
  file <- read_tsv(file.names[i])
  file %>% mutate(Source=file.names[i]) -> O
  if (dim(O)[1] >= 1){ All_files = rbind(All_files, O) }
}
write_tsv(All_files, "Results_16S_noBMI.tsv") #NOT BMI




############################################################################
#########################METAANALYSIS#######################################
############################################################################


get_level = function(x){
  x = str_split(x, "\\.")[[1]][1]
  return(x)
}


Meta_analyze = function(Gender, Results, BMI=F){
    if (Gender == "Male"){
      Results %>% filter(grepl("\\.0\\.",Source)) %>% select(-c(FDR, Source)) -> Results_LLD
      if (BMI == F){ Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_Jul21_2021_MALES.csv")
      } else {Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_BMI_Jul21_2021_MALES.csv") }
    
    } else if (Gender == "Female"){
      Results %>% filter(grepl("\\.1\\.",Source)) %>% select(-c(FDR, Source)) -> Results_LLD
      if (BMI == F){ Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_Jul21_2021_FEMALES.csv")
      } else {Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_BMI_Jul21_2021_FEMALES.csv") }
      
    } else{
      Results %>% filter(grepl("\\.3\\.",Source)) %>% select(-c(FDR, Source)) -> Results_LLD
      if (BMI == F){ Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_All_Jul21_2021.csv")
      } else {Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_All_BMI_Jul21_2021.csv") }
    }
  
  #Results_Rotter %>% mutate(Metabolite = ifelse(Metabolite == "Deoxycarnitine",  "y-butyrobetaine", Metabolite))
  Results_LLD %>% mutate(Metabolite = ifelse(Metabolite == "y-butyrobetaine",  "Deoxycarnitine", Metabolite)) %>%
                  mutate(Metabolite = ifelse(Metabolite == "L-Carnitine",  "Carnitine", Metabolite) ) -> Results_LLD
  if ("N_samples" %in% colnames(Results_LLD)){ Results_LLD %>% select(-N_samples) ->  Results_LLD}
  
  Results_LLD$Level =  sapply(Results_LLD$Bug, get_level) ; Results_Rotter$Level =  sapply(Results_Rotter$Bug, get_level)
  Results_Rotter %>% select(-X1) %>% filter(!Level == "SampleID") -> Results_Rotter
  
  print("Number of species test")
  length(unique(Results_LLD$Bug)) ; length(unique(Results_Rotter$Bug))
  
  print("Plot summary stats")
  Results_LLD %>% ggplot() + geom_point(aes(x=Beta, y=-log10(Pvalue), col = Pvalue<0.05 )) +
    facet_grid(Metabolite~Level, scales="free") + theme_bw() + ggtitle('LLD') -> Summary_LLD
  Results_Rotter %>% ggplot() + geom_point(aes(x=Beta, y=-log10(Pvalue),  col = Pvalue<0.05)) +
    facet_grid(Metabolite~Level, scales="free") + theme_bw() + ggtitle('Rotterdam') -> Summary_rot
  Summary_LLD + Summary_rot + plot_layout(guides = "collect") -> Summary_cohorts
  
  
    ###########################################
    ##############Meta analysis################
    ###########################################
    print("Iteration through taxa")
    #size1 = 1032
    #size2 = 1427 
    
    Meta_stats = tibble()
    for (BUG in unique(Results_LLD$Bug)){
      if (! BUG %in% Results_Rotter$Bug){ next }
      Results_LLD %>% filter(Bug == BUG) -> Bug_LL
      Results_Rotter %>% filter(Bug == BUG)  -> Bug_R
      
      for (Metabolite_i in Bug_LL$Metabolite){
        Bug_LL %>% filter(Metabolite == Metabolite_i) -> Entry_L ; Bug_R %>% filter(Metabolite == Metabolite_i) -> Entry_R
        
        #p1.z	=	abs(qnorm(Entry_L$Pvalue/2)) ; p2.z	=	abs(qnorm(Entry_R$Pvalue/2))
        #if (Entry_L$Beta<0) p1.z = -1*p1.z
        #if (Entry_R$Beta<0) p2.z	=	-1*p2.z
        #size		= size1+size2
        #w.meta.z		=	(p1.z*size1+p2.z*size2)/sqrt(size1^2+size2^2)
        #w.meta.p		=	(1-pnorm(abs(w.meta.z)))*2
        
        ##Meta function
        #meta::metacont(n.e = size1, mean.e = Entry_L$Beta, sd.e = Entry_L$SE, n.c=size2, mean.c=Entry_R$Beta, sd.c=Entry_R$Pvalue) -> meta_value
        rbind(Entry_L, Entry_R) -> Input_data
        meta::metagen(TE=Beta, seTE=SE, data=Input_data, comb.fixed = T, comb.random= T ) -> meta_value
        
        Agreement = ifelse(sign(Entry_L$Beta) == sign(Entry_R$Beta), "Agree", "Desagree")
        Sub_result = tibble(Bug = BUG, Metabolite = Metabolite_i, MetaP=meta_value$pval.fixed, Treatment_effect=meta_value$TE.fixed,SE_meta=meta_value$seTE.fixed,Heterogeneity_Q=meta_value$Q,Heterogeneity_Pvalue=meta_value$pval.Q,  Concordance = Agreement , Meta_random_P=meta_value$pval.random, Treatment_effect_random=meta_value$TE.random, SE_meta_random=meta_value$seTE.random)
        
        left_join(x=Entry_L, y=Entry_R, by=c("Metabolite", "Bug"), suffix = c("_LLD", "_Rott")) -> Daa
        Sub_result =  left_join(x=Sub_result, y=Daa, by= c("Metabolite", "Bug") )
        
        Meta_stats = rbind(Meta_stats, Sub_result)
      }  
    }
    
    Meta_stats %>% mutate(meta_FDR = p.adjust(MetaP, "fdr")) %>% mutate(Bonferroni_indp = MetaP*134*5) -> Meta_stats
    
    
    #####PLOTS
    print("Uncover taxa associated with more than one metabolite")
    
    Meta_stats %>% filter(meta_FDR < 0.07) %>% filter(Metabolite %in% c("Betaine","Carnitine","Choline", "Deoxycarnitine", "TMAO") ) -> Significant
    Final_figure = NULL
    for (B in unique(Significant$Bug)){
      Significant %>% filter(Bug == B) -> Taxa
      if (dim(Taxa)[1] > 1){
        print(B)
        ggplot(Taxa) + geom_bar(aes_(y=Taxa$Treatment_effect,x=Taxa$Metabolite), stat="identity") + theme_bw() +ggtitle(B) +theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +ylab("Effect") + coord_flip() -> Fig
                                                                                                                                  
                                                                                                                                  
        if (is.null(Final_figure)){ Final_figure = Fig
        }else{ Final_figure = Final_figure  + Fig}
      }
    }
    
    asso_id = function(x,y){
      paste(x,y, sep = "-")
    }
    Assoc_ID = asso_id(Meta_stats$Bug, Meta_stats$Metabolite)
    Meta_stats %>% mutate(Association_ID =  Assoc_ID) -> Meta_stats
    Fig_conc = Meta_stats %>% ggplot(aes(x=Concordance, y=-log10(MetaP) )) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(col = (Heterogeneity_Q<0.05), shape= (meta_FDR<0.05) )) + theme_bw() + geom_hline(yintercept = -log10(0.05), col="grey", linetype="dashed",size=2)
    
    if (BMI == F){ Output_name = paste(c("Metanalysis/Summary_stats_", Gender,"_BMI-",as.character(F) ,".tsv"), collapse="")
    } else{Output_name = paste(c("Metanalysis/Summary_stats_", Gender,".tsv"), collapse="") }
    write_tsv(x= Meta_stats, path = Output_name)
    
    return(list(arrange(Meta_stats, MetaP), Summary_cohorts, Final_figure, Fig_conc) )
   
}


Gender_specific = function(All, Females, Males, Threshold=0.001 ){
  Significant_all = filter(All, MetaP < Threshold)
  #Search for male associations not in All
  Males %>% filter(! Association_ID %in% Significant_all$Association_ID) %>% filter(MetaP < Threshold) -> Bugs_different_male
  #Search for Female associations not in All
  Females %>% filter(! Association_ID %in% Significant_all$Association_ID) %>% filter(MetaP < Threshold) -> Bugs_different_female
  
  #Put together ; left join with Males/Females to get them all. Left join with Bugs_different_male/Bugs_different_female to get the "significant"
  All %>% filter(Association_ID %in% c(Bugs_different_male$Association_ID, Bugs_different_female$Association_ID)) -> Stats_all_in_diff
  #left_join(x=Stats_all_in_diff, y=Bugs_different_male, by=c("Bug", "Metabolite", "Association_ID"), suffix=c("","_Males")) -> Stats_all_in_diff
  #left_join(x=Stats_all_in_diff, y=Bugs_different_female, by=c("Bug", "Metabolite", "Association_ID"), suffix=c("","_Females")) -> Stats_all_in_diff
  left_join(x=Stats_all_in_diff, y=Males, by=c("Bug", "Metabolite", "Association_ID"), suffix=c("","_Males")) -> Stats_all_in_diff
  left_join(x=Stats_all_in_diff, y=Females, by=c("Bug", "Metabolite", "Association_ID"), suffix=c("","_Females")) -> Stats_all_in_diff
  
  #Format to get effect sizes and SE
  #1. Mean
  Stats_all_in_diff %>% select(Bug, Metabolite, MetaP, Treatment_effect, MetaP_Males, Treatment_effect_Males, MetaP_Females, Treatment_effect_Females ) -> Stats_diff
  Stats_diff %>% select(-c(MetaP,MetaP_Females, MetaP_Males)) %>% gather(Group, Effect, c(Treatment_effect,Treatment_effect_Males, Treatment_effect_Females), factor_key=TRUE) -> Stats_diff_long
  #2.Standard Error
  Stats_all_in_diff %>% select(Bug, Metabolite, SE_meta, SE_meta_Females, SE_meta_Males ) %>% gather(Group, StdE, 3:5, factor_key=TRUE) -> Stats_se_long
  #3. Format "Group" uniformely
  Stats_diff_long %>% mutate(Group = ifelse(grepl("Males",Group), "Males", ifelse(grepl("Females", Group), "Females", "All"))) -> Stats_diff_long
  Stats_se_long  %>% mutate(Group = ifelse(grepl("Males",Group), "Males", ifelse(grepl("Females", Group), "Females", "All"))) -> Stats_se_long
  #4. Add together
  left_join(x=Stats_diff_long, y=Stats_se_long, by=c("Bug", "Metabolite", "Group")) -> Stats_diff_long
  #6. Plot
  Stats_diff_long %>% ggplot(aes(x=Bug, y = Effect, col=Group)) + geom_point(position=position_dodge(width=0.3)) + facet_wrap(~Metabolite, scales = "free") +
    geom_errorbar(position=position_dodge(width=0.3), aes(ymin=Effect-(2*StdE), ymax=Effect+(2*StdE)), width=.2)  + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    coord_flip() + theme_bw()-> FPlot
  
  return(FPlot)
  
}


####################################
###############NO BMI ##############
####################################

Data_LLD = read_tsv("Results_16S_noBMI.tsv")
Outcome = Meta_analyze(Gender="All", Results = Data_LLD, BMI=F)
Stats =Outcome[[1]] ;  Fig_cohort = Outcome[[2]] ;  Fig_direct = Outcome[[3]] ;Fig_conc = Outcome[[4]]
Stats %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(7.462687e-05))
Stats %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()


Outcome_m =  Meta_analyze("Male", Data_LLD, BMI=F)
Stats_Male =Outcome_m[[1]] ;  Fig_cohort_Male = Outcome_m[[2]] ;  Fig_direct_Male = Outcome_m[[3]] ; Fig_conc_Male = Outcome_m[[4]]
Stats_Male %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats_Male %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(7.462687e-05))
Stats_Male %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()


Outcome_f = Meta_analyze("Female", Data_LLD, BMI=F)
Stats_Female =Outcome_f[[1]] ;  Fig_cohort_Female = Outcome_f[[2]] ;  Fig_direct_Female = Outcome_f[[3]] ; Fig_conc_Female = Outcome_f[[4]]
Stats_Female %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats_Female %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(7.462687e-05))
Stats_Female %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()


Gender_specific(Stats, Stats_Male, Stats_Female, Threshold = 0.001)

####################################
############### BMI ##############
####################################
Data_LLD = read_tsv("Results_16S.tsv")
Outcome = Meta_analyze("All", Data_LLD, BMI=T)
Stats =Outcome[[1]] ;  Fig_cohort = Outcome[[2]] ;  Fig_direct = Outcome[[3]] ;Fig_conc = Outcome[[4]]
Stats %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(7.462687e-05))
Stats %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()


Outcome_m =  Meta_analyze("Male", Data_LLD, BMI=T)
Stats_Male =Outcome_m[[1]] ;  Fig_cohort_Male = Outcome_m[[2]] ;  Fig_direct_Male = Outcome_m[[3]] ; Fig_conc_Male = Outcome_m[[4]]
Stats_Male %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats_Male %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(7.462687e-05))
Stats_Male %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()



Outcome_f = Meta_analyze("Female", Data_LLD, BMI=T)
Stats_Female =Outcome_f[[1]] ;  Fig_cohort_Female = Outcome_f[[2]] ;  Fig_direct_Female = Outcome_f[[3]]; Fig_conc_Female = Outcome_m[[4]]
Stats_Female %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats_Female %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(7.462687e-05))
Stats_Female %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()


Gender_specific(Stats, Stats_Male, Stats_Female, Threshold = 0.001)


######Gender hetergoneity test####
Stats_Female_noBMI = read_tsv("Metanalysis/Summary_stats_Female_BMI-FALSE.tsv")
Stats_Male_noBMI = read_tsv("Metanalysis/Summary_stats_Male_BMI-FALSE.tsv")

rbind(Stats_Male_noBMI, Stats_Female_noBMI) -> Input_data
paste(Input_data$Bug, Input_data$Metabolite, sep = "-") ->ID_a
Input_data %>% mutate(ID = ID_a) -> Input_data
Hetero_results = tibble()
for (id in unique(ID_a)){
  Input_data %>% filter(ID == id) -> Input
  meta::metagen(TE=Treatment_effect, seTE=SE_meta, data=Input, comb.fixed = T, comb.random= F ) -> meta_value
  N = tibble(ID=id,Q_sex=meta_value$Q, P_Qsex= meta_value$pval.Q)
  rbind(Hetero_results, N) -> Hetero_results
}
unique(Hetero_results) %>% arrange(P_Qsex) %>% mutate(FDR= p.adjust(P_Qsex, "fdr"))
###################################
####Multiple test correction#######
###################################

#This is run once in order to obtain how many independent taxa we have.
#We filter samples as done in association step. Then filter bacteria used for metaanalysis
#B-C PCoA is done and we check how many components are needed to reach 90%

Count_table =  "LLD_rarefDefault.taxonomyTable.txt.gz"
Linking_table = "coupling_pheno_16s.txt"

Count_table = read_tsv(Count_table) #Rows are individuals, columns are Bacteria
Linking_table = read_tsv(Linking_table, col_names=F) 
Count_table = Filter_unclassified(Count_table)
Linking_table = Choose_replicate(Linking_table)

Count_table %>% filter(SampleID %in% Linking_table$X2) %>% select(Stats$Bug) -> Table_PCA


res.pca <- prcomp(Table_PCA, scale = TRUE)
Variability_explained = res.pca$sdev/sum(res.pca$sdev)
SUMA = cumsum(Variability_explained)
Success = SUMA >= 0.9
tibble(component= seq(length(SUMA)), Variability=Variability_explained, Cumulative_v = SUMA, Above_90 = Success) -> Results
Results %>% filter(Success == F) # 134 indp

as.matrix(Table_PCA) %*% res.pca$rotation -> PCs
as_tibble(PCs) %>% ggplot(aes(x=PC1, y=PC2)) + geom_point()

#PCOA, somehow number of dimensions is odd
dist = vegdist(Table_PCA, method = "euclidean")
PCOA = pcoa(dist)
dim(PCOA$vectors) 
as_tibble(PCOA$values) %>% mutate(C = cumsum(Relative_eig)) %>% filter(C<0.9)
