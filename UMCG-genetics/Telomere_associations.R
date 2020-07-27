setwd("/Users/sergio/Documents/PhD/Telomeres_project/")

library(tidyverse)

#Phenotypes and selection of phenotypes to be removed
Phenotypes <- read_delim("Merged_Phenotypes_LLD.csv", delim="|")
Phenotypes %>% select(- colnames(Phenotypes)[grepl("_1",colnames(Phenotypes))]) -> Phenotypes
To_remove <- read_tsv("to_be_excluded.txt", col_names = F) %>% filter(X1 %in% colnames(Phenotypes))

Phenotypes %>% select(! To_remove$X1) -> Phenotypes 
colnames(Phenotypes)[which(grepl("or have they had diabetes mellitus type 2",colnames(Phenotypes)))] = "Parents_diabetes2"

#colnames(Phenotypes)[grepl("smk", colnames(Phenotypes))] = c("Have_you_smoked_a_year?","How_old_when_started_smoking?","father_smk","mother_smk", "mother_smk_pregnancy",	"Number_people_smoke_at_home", "Do_people_smoke_work", "smk_participant")

all_cell_types = read_tsv("BIOS_cell_types_DeconCell_2019-03-08.LLD_subset.txt")
all_cell_types %>% gather(key = var_name, value = value, 2:ncol(all_cell_types)) %>% 
  spread_(key = names(all_cell_types)[1],value = 'value') %>% mutate(ID = var_name) %>% select(-var_name) -> all_cell_types
cell_types_PC = read_tsv("BIOS_cell_types_DeconCell_2019-03-08.LLD_subset_PCA.txt")

##Add new pheno
Parent_smoke = vector()
for (ID in 1:length(Phenotypes$smk13)){
  if( is.na(Phenotypes$smk13[ID]) | is.na(Phenotypes$smk14[ID]) ){
    Parent_smoke = c(Parent_smoke, NA)
  } else if (Phenotypes$smk13[ID] == 1 | Phenotypes$smk14[ID] == 1){
    Parent_smoke = c(Parent_smoke, 1)
  } else{
    Parent_smoke = c(Parent_smoke, 0)
  }
}
Phenotypes %>% mutate(Parental_smoking = Parent_smoke)  -> Phenotypes

Phenotypes %>% mutate(age_m= ifelse(ID == "LLDeep_0968", 24, age_m)) %>% mutate(age_f = ifelse(ID == "LLDeep_1127", 30, age_f)) -> Phenotypes

#All biological aging markers, filter telomeres
Age_markers <- read_tsv("Combined_data_uncorrected.tsv")
Telomere_markers <- colnames(Age_markers)[grepl("MTL",colnames(Age_markers))]
Age_markers %>% select(c("Sample",Telomere_markers)) %>% mutate(MTL_gran = as.numeric(MTL_gran), `MTL_CD20+`=as.numeric(`MTL_CD20+`), `MTL_CD45+_CD20-`=as.numeric(`MTL_CD45+_CD20-`)) %>%
  drop_na() -> Telomere_markers
column_cells = colnames(Phenotypes)[grepl("ytes", colnames(Phenotypes))]
#Make deltas
select(Age_markers, Sample) -> Outcome_delta
for (i in c( "Methyl.Age.Hannum.", "dCT", "MTL_lymp","MTL_gran","MTL_CD45+_CD20-","MTL_CD45-","MTL_CD20+","MTL_CD57+")){
  #Do something with the residuals, they dont have info about the IDs...
  Age_markers %>% select(i) %>% as_vector() %>% as.vector() -> Marker
  Phenotypes %>% select(c("ID", "AgeIncludingMonth", "Sex", column_cells)) %>% mutate(Sample = ID) %>% select(-ID) -> To_add
  left_join(select(Age_markers, Sample), To_add)%>% mutate(Marker = as.numeric(Marker))  %>% drop_na() -> Age_marker
  as.vector(as_vector(lm(Marker ~ ., select(Age_marker, - Sample))$residuals)) -> Delta
  Age_marker %>% mutate(Delta = Delta) %>% select(Sample, Delta) -> Output_delta
  colnames(Output_delta) = c("Sample", paste(c("Delta_",i),collapse="")) 
  left_join(Outcome_delta, Output_delta) -> Outcome_delta
}
write_tsv("~/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/Delta_BioAge.tsv",x = Outcome_delta)

Outcome_delta %>% drop_na() -> CHECK
sum(CHECK$Delta_Methyl.Age.Hannum.^2)
Age_markers %>% select(c(`Methyl.Age.Hannum.`,Age)) %>% drop_na() -> CHECK2
sum((CHECK2$Methyl.Age.Hannum. -  CHECK2$Age)^2)


#Proteomics 
read_tsv("CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt") -> Proteomics ; Proteomics %>% select(-Protein) %>% t() %>% as_tibble() %>% mutate(ID= colnames(Proteomics)[2:length(colnames(Proteomics))]) %>%
  `colnames<-`(c(Proteomics$Protein, "ID")) %>% arrange(ID) -> Proteomics
read_tsv("~/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/age_gender_smk_contrac_cell_counts.txt")  -> Protein_covariates
as_tibble(cbind(nms = names(Protein_covariates), t(Protein_covariates))) %>%  `colnames<-`(c("ID",Protein_covariates$personid)) %>%
  filter(!ID == "personid") %>% apply(2,FUN=as.numeric) %>% as_tibble() %>% mutate(ID = colnames(Protein_covariates)[2:length(colnames(Protein_covariates))]) -> Protein_covariates

#################
###FUNCTIONS#####
#################
Linear_regression = function(Dependent, Regressor,Covariates, Regressor_name, Dependent_name, Plot = F,cells=F, BMI=F){
  Regressor = as.numeric(Regressor)
  if (BMI == T){  Model_data = tibble(Dependent=Dependent, Regressor=Regressor, Sex= factor(Covariates$Sex), Age= as.numeric(Covariates$Age), ID = Covariates$ID, BMI= Covariates$`Body Mass Index (kg/M^2)`)
  }else{Model_data = tibble(Dependent=Dependent, Regressor=Regressor, Sex= factor(Covariates$Sex), Age= as.numeric(Covariates$Age), ID = Covariates$ID)}
  Model_data %>% drop_na() -> Model_data
  if (BMI == T){
    left_join(Model_data, select(Covariates, c("ID",names(cell_types))),by="ID") -> Model_data
    Model_data %>% drop_na() -> Model_data
    paste(c("Dependent ~ Sex + Age + BMI",cell_types, "Regressor"), collapse = "+") -> Model_f
  }
  else if (cells==T){
    left_join(Model_data, select(Covariates, c("ID",names(cell_types))),by="ID") -> Model_data
    Model_data %>% drop_na() -> Model_data
    paste(c("Dependent ~ Sex + Age",cell_types, "Regressor"), collapse = "+") -> Model_f
  }else{
    c("Dependent ~Sex + Age + Regressor") -> Model_f
  }
  Model = lm(Model_f , Model_data )
  Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Regressor",])
  
  Beta = Summary$Estimate
  pvalue = Summary$`Pr(>|t|)`
  stand = Summary$`Std. Error`
  N = length(summary(Model)$residuals)
  Results = as_tibble(t(c(Regressor_name, Dependent_name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
  colnames(Results) = c("Regressor", "Dependent", "Beta", "SE", "Pvalue", "N")
  
  if (Plot == T){
    RES = lm(Regressor ~ Age ,Model_data)$residuals
    
    Model_data %>% mutate(RES = RES) -> Model_data
    if (length(unique(Dependent)) > 8){
      paste(c(Regressor_name, Dependent_name), collapse= " ") -> NAME
      Figure1 = ggplot(data=Model_data,aes(x=RES, y=Dependent)) + geom_point() + theme_bw() + geom_smooth(method ="lm") + ggtitle(NAME)
    } else{ Figure1 = ggplot(data=Model_data,aes(x=factor(Dependent), y=RES)) + geom_boxplot() + theme_bw() + geom_smooth(method ="lm") }
    print(Figure1)
    #try(suppressMessages(ggsave(plot=Figure1,filename=paste(c("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Plots/",Regressor_name,"vs",Dependent_name,".png"),collapse=""))))
  }
  
  
  return(Results)
}

lm_proteomics = function(Dependent, Regressor,Covariates, Regressor_name, Dependent_name){
  Regressor = as.numeric(Regressor)
  Covariates %>% mutate(Regressor = Regressor, Dependent = Dependent) -> Input_data
  Model = lm(Dependent ~ ., data=select(Input_data, -ID))
  Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Regressor",])
  
  Beta = Summary$Estimate
  pvalue = Summary$`Pr(>|t|)`
  stand = Summary$`Std. Error`
  N = length(summary(Model)$residuals)
  Results = as_tibble(t(c(Regressor_name, Dependent_name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
  colnames(Results) = c("Regressor", "Dependent", "Beta", "SE", "Pvalue", "N")
  return(Results)
}


Logistic_regression = function(Dependent, Regressor,Covariates, Regressor_name, Dependent_name, Plot=T, cells=F){
  Regressor = as.numeric(Regressor)
  #If instead of a binary 0-1 it is 1-2, make it 0-1
  if (length(which(Dependent == 2)) > 0){ Dependent = as.numeric(as.factor(Dependent)) - 1}
  Model_data = tibble(Dependent=Dependent, Regressor=Regressor, Sex= factor(Covariates$Sex), Age= as.numeric(Covariates$Age), ID = Covariates$ID)
  Model_data %>% filter(! is.na(Dependent)) %>% filter(! is.na(Regressor)) -> Model_data
  
  
  if (cells==T){
    left_join(Model_data, select(Covariates, c("ID",names(cell_types))),by="ID") -> Model_data
    Model_data %>% drop_na() -> Model_data
    paste(c("Dependent ~ Sex + Age",cell_types,"Regressor"), collapse = "+") -> Model_f
  }else{
    c("Dependent ~ Sex + Age + Regressor") -> Model_f
  }

  Model = glm(Model_f, Model_data, family=binomial(link='logit'))

  Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Regressor",])
  Beta = Summary$Estimate
  pvalue = Summary$`Pr(>|z|)`
  stand = Summary$`Std. Error`
  N = length(summary(Model)$deviance.resid)
  Results = as_tibble(t(c(Regressor_name, Dependent_name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
  colnames(Results) = c("Regressor", "Dependent", "Beta", "SE", "Pvalue", "N")
  if (Plot == T ){
    RES = lm(Regressor ~ Age ,Model_data)$residuals
    Figure1 = ggplot(data=Model_data, aes(x=as.factor(Dependent), y=RES)) + geom_boxplot() + geom_point() + theme_bw()
    print(Figure1)
    #try(suppressMessages(ggsave(plot=Figure1,filename=paste(c("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Plots/",Regressor_name,"vs",Dependent_name,".png"),collapse=""))))
  }
  
  
  return(Results)
}

clean_name = function(Name){
  Name= str_trim(Name)
  Name = gsub(" ", "_",Name, fixed= T)
  Name = gsub("(", "",Name, fixed= T)
  Name = gsub(")", "",Name, fixed= T)
  Name = gsub("%", "",Name, fixed= T)
  Name = gsub("/", "",Name, fixed= T)
  return(Name)
}

Model_phenotypes = function(Phenotypes, Data, Covariates, cells=NULL, BMI=F, Proteomics = F ){ 
  Overall = tibble()
  for (Col_n in 1:ncol(Phenotypes)){
    Name = colnames(Phenotypes)[Col_n]
    print(Name)
    Name= clean_name(Name)
    Dependent = as.vector(as_vector((Phenotypes[Col_n])))
    #If it is a binary outcome, perform a logistic regression. Otherwise, OLS
    if (length(unique(Dependent[!is.na(Dependent)])) == 2){ 
      Categorical = T
      Dependent = as.numeric(as.factor(Dependent))
    }else if (length(unique(Dependent[!is.na(Dependent)])) >= 2){ 
      Categorical = F
      Dependent = as.numeric(Dependent)
    }else{ next }
    
    for (Col_m in 1:ncol(Data)){
      Explanatory  = as.vector(as_vector(Data[,Col_m]))
      Name_2 = colnames(Data)[Col_m]
      Name_2= clean_name(Name_2)
      To_plot = c()#c("age_f", "age_m", "Body_Mass_Index_kgM^2", "Do_your_parents,_siblings_or_children_have_diabetes_mellitus_type_2,_or_have_they_had_diabetes_mellitus_type_2?_No.", "Have_you_ever_had_poorly_healing_wounds_on_your_feet?","Parental_smoking","smk13","Waist_circumference_in_cm" )
      if (Name %in% To_plot){ PL = T
      } else{ PL = F}
      if (Categorical == T){
        if (is.null(cells)){
          Results = Logistic_regression(Dependent, Explanatory, Covariates, Name_2, Name, Plot= PL)
        } else{ Results = Logistic_regression(Dependent, Explanatory, Covariates, Name_2, Name, Plot= PL, cells =T) }
      }else{
        if (Proteomics == T){ 
          Results = lm_proteomics(Dependent,Explanatory, Covariates, Name_2, Name)
        }else if (is.null(cells)){
          Results = Linear_regression(Dependent, Explanatory, Covariates,Name_2, Name, Plot = PL)
        }else{ Results = Linear_regression(Dependent, Explanatory, Covariates, Name_2, Name, Plot= PL, cells =T, BMI=BMI) }  
      }

      #Count how many entries per level
      if (length(unique(Dependent[!is.na(Dependent)])) < 10){ 
        LEVEL =  paste(unique(Dependent[!is.na(Dependent)]),collapse=",")
        N_LEVEL = c()
        for (L in str_split(LEVEL, ",")[[1]]){
          N_LEVEL = c(N_LEVEL,length(which(as.character(Dependent) == L)))
        }
        N_LEVEL = paste(N_LEVEL, collapse=",")
      }else{ 
        LEVEL =  "Numeric" 
        N_LEVEL = length(Dependent[!is.na(Dependent)])
      }
      
      Results %>% mutate(Levels = LEVEL, Levels_n = N_LEVEL) -> Results
      Overall = rbind(Overall, Results)
      }
    }
  return(Overall)     
}

Plot_phenotypes = function(Phenotypes, Data, Covariates){
  Phenos = c("age_f", "age_m", "Body Mass Index (kg/M^2)", "Do your parents, siblings or children have diabetes mellitus type 2, or have they had diabetes mellitus type 2? No.", "Have you ever had poorly healing wounds on your feet?","Parental_smoking","smk13","Waist circumference in cm" )
  Phenotypes %>%  select(c("ID",Phenos)) -> Phenotypes
  #left_join(Phenotypes, Covariates) -> Model_input
  Data  %>% mutate(ID = Sample) %>% select(-Sample) -> Data
  left_join(Covariates, Data, by="ID") -> Model_input
  for (Phenotype in colnames(Phenotypes)[2:length(colnames(Phenotypes))]){
    Phenotypes %>% select(Phenotype) %>% as_vector() %>% as.vector() -> D
    Model_input %>% mutate(Dependent = D) -> Model
    
    
    Model %>% gather(Cell_line, Telomere_length, colnames(select(Data, -ID)) ) -> Model
    
    if (Phenotype %in% c("Do your parents, siblings or children have diabetes mellitus type 2, or have they had diabetes mellitus type 2? No.", "Have you ever had poorly healing wounds on your feet?")){
      Model %>% drop_na() %>% ggplot(aes(x=as.factor(Dependent), y=as.numeric(Telomere_length))) + ggforce::geom_sina() + theme_bw() + facet_wrap(~Cell_line) + geom_smooth(method = "lm") -> PLOT
    } else{
      Model %>% ggplot(aes(y=as.numeric(Dependent), x=as.numeric(Telomere_length))) + geom_point() + theme_bw() + facet_wrap(~Cell_line) + geom_smooth(method = "lm", formula = y~x) -> PLOT
    }
      print(PLOT)
  }
  
}


Calculate_fdr = function(Phenotypes, Data, Covariates, Real_pvalues, cells=NULL, BMI=F){
  FDR_vector = vector()
  #1st, Number of permutations
  Random_distribution  = vector()
  for (i in seq(100)){
    #Per each permutation rearrange data and fit model
    #rearrange
    sample_n(Phenotypes, size= nrow(Phenotypes)) -> Rearrange
    #Fitmodel
    Model_phenotypes(Rearrange, Data, Covariates,cells, BMI) -> H0
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

Compare_outlaiers_bioage = function(DF){
  AGE = DF$AgeIncludingMonth
  DF_all_shared = tibble()
  for (BioAge in colnames(DF)){
    if (BioAge %in% c("ID","AgeIncludingMonth")){ next }
    DF %>% select(BioAge) %>% as_vector() %>% as.vector() -> bioage
    lm(bioage ~ AGE) -> Model
    
    residual_sd = summary(Model)$sigma
    predict(Model) -> Mean_value
    Upper_threshold = Mean_value + 1.2*residual_sd #1.2sd
    Lower_threshold = Mean_value - 1.2*residual_sd
    
    DF$ID[as.vector(which(bioage > Upper_threshold))] -> outliers_upper
    DF$ID[as.vector(which(bioage < Lower_threshold))] -> outliers_lower
    DF %>% select(ID) %>% mutate(Up_outlier = ifelse(ID %in% outliers_upper, T, F), Down_outlier = ifelse(ID %in% outliers_lower, T, F), Marker=BioAge ) -> DF_shared
    if (BioAge == "Methyl.Age.Hannum."){ DF %>% select(ID) %>% mutate(Up_outlier = ifelse(ID %in% outliers_lower, T, F), Down_outlier = ifelse(ID %in% outliers_upper, T, F), Marker=BioAge ) -> DF_shared  }
    DF_all_shared = rbind(DF_all_shared,DF_shared)  
  }
  DF_all_shared %>% filter(Up_outlier == T | Down_outlier == T) %>%  group_by(Marker) %>% summarise(n())

  DF_all_shared %>% mutate(Outlier = Up_outlier+Down_outlier) %>%  select(-c(Up_outlier, Down_outlier)) %>% 
    spread(Marker, Outlier) -> Upset_input
  upset(as.data.frame(Upset_input), text.scale = 2,  nsets = 20) -> upset_shared_outliers
  
  #Do the same for negative and positive
  DF_all_shared %>% mutate(Outlier = Up_outlier*1) %>%  select(-c(Up_outlier, Down_outlier)) %>% 
    spread(Marker, Outlier) -> Upset_input_up
  DF_all_shared %>% mutate(Outlier = Down_outlier*1) %>%  select(-c(Up_outlier, Down_outlier)) %>% 
    spread(Marker, Outlier) -> Upset_input_down 
  
  upset(as.data.frame(Upset_input_up), text.scale = 2,nintersects=NA, nsets = 20) -> upset_shared_outliers_SlowAging
  upset(as.data.frame(Upset_input_down), text.scale = 2,nintersects=NA,  nsets = 20) -> upset_shared_outliers_FastAging
  print(upset_shared_outliers_SlowAging)
  print(upset_shared_outliers_FastAging)
}


Compare_outlaiers_bioage2 = function(DF = Age_markers, Phenos = Phenotypes ){
  DF_all_shared = tibble()
  apply(select(DF, -Sample),2,FUN=as.numeric) %>% as_tibble() %>% mutate(Sample = DF$Sample) -> DF
  for (BioAge in colnames(DF)[grepl("MTL",colnames(DF))]){
    if (BioAge %in% c("ID","AgeIncludingMonth")){ next }
    DF %>% select(c(Sample,BioAge)) %>% drop_na() -> temp_df
    Phenos %>% filter(ID %in% temp_df$Sample) %>% arrange(ID) %>% drop_na(AgeIncludingMonth) -> AGE  
    temp_df %>% filter(Sample %in% AGE$ID) %>% arrange(Sample) %>% select(BioAge) %>% as_vector() %>% as.vector() -> bioage
    lm(bioage ~ AGE$AgeIncludingMonth) -> Model
    
    residual_sd = summary(Model)$sigma
    predict(Model) -> Mean_value
    Upper_threshold = Mean_value + 1*residual_sd #1.2sd
    Lower_threshold = Mean_value - 1*residual_sd
    
    AGE %>% mutate(bioage) %>% filter(bioage > Upper_threshold) %>% select(ID) %>% as_vector() %>% as.vector() ->outliers_upper
    AGE %>% mutate(bioage) %>% filter(bioage < Lower_threshold) %>% select(ID) %>% as_vector() %>% as.vector() ->outliers_lower
    
    DF %>% select(Sample) %>% mutate(ID=Sample) %>% select(-Sample) %>% mutate(Up_outlier = ifelse(ID %in% outliers_upper, T, F), Down_outlier = ifelse(ID %in% outliers_lower, T, F)) %>%
      mutate(Marker = BioAge) -> DF_all
    rbind(DF_all_shared, DF_all) -> DF_all_shared
  }
  DF_all_shared %>% filter(Up_outlier == T | Down_outlier == T) %>%  group_by(Marker) %>% summarise(n())
  DF_all_shared %>% mutate(Status = Up_outlier+Down_outlier) %>% group_by(ID) %>% summarise(N=sum(Status)) %>% arrange(desc(N)) %>% filter(N>=2) ->
    Outliers
  DF_all_shared  %>% filter(ID %in% Outliers$ID) %>% mutate(Up_outlier = ifelse(Up_outlier == TRUE, 1, 0), Down_outlier = ifelse(Down_outlier == TRUE, -1, 0)) %>% filter(!Down_outlier+Up_outlier == 0) %>% 
    mutate(Outlier = Up_outlier+Down_outlier) %>% group_by(ID) %>% summarise(N= sum(Outlier), C=n() ) %>% filter(! abs(N) != C ) ->
    Groups 
  
  Down_group = filter(Groups, N<0 )$ID
  Up_group =  filter(Groups, N>0 )$ID
  
  select(DF, Sample, Methyl.Age.Hannum.) %>% mutate(ID=Sample) %>% select(-Sample) %>% mutate(Group = ifelse(ID %in% Up_group, "1-Long",ifelse(ID %in% Down_group, "0-Short",NA))) %>% drop_na() -> Methylation
  left_join(Methylation,Phenos) -> Methylation
  summary(glm(as.factor(Group) ~ Methyl.Age.Hannum. + AgeIncludingMonth + Sex, Methylation, family = binomial(link="logit")))
  Methylation %>% select(Group,  Methyl.Age.Hannum.,AgeIncludingMonth,Sex, colnames(Methylation)[grepl("cytes", colnames(Methylation))] ) -> Methylation2
  summary(glm(as.factor(Group) ~ . , Methylation2, family = binomial(link="logit")))
  
  select(DF, Sample, dCT) %>% mutate(ID=Sample) %>% select(-Sample) %>% mutate(Group = ifelse(ID %in% Up_group, "1-Long",ifelse(ID %in% Down_group, "0-Short",NA))) %>% drop_na() -> TRECS
  left_join(TRECS,Phenos) -> TRECS
  summary(glm(as.factor(Group) ~ dCT + AgeIncludingMonth + Sex, TRECS, family = binomial(link="logit")))
  TRECS %>% select(Group,  dCT,AgeIncludingMonth,Sex, colnames(TRECS)[grepl("cytes", colnames(TRECS))] ) -> TRECS2
  summary(glm(as.factor(Group) ~ . , TRECS2, family = binomial(link="logit")))
  
  
  
}



SAVE = function(NAME, PLOT, Width = 6.3, Height = 6.3, Dpi = 300, Scale = 1){
  ggsave(filename = NAME, PLOT,
         width = Width, height = Height, dpi = Dpi, units = "in", scale = Scale)
}
#################################
##Prepare input for Analysis#####
#################################

Phenotypes %>% filter(ID %in% Telomere_markers$Sample) %>% arrange(ID) %>%
filter(! is.na(AgeIncludingMonth)) %>% filter(! is.na(Sex)) -> Phenotypes 

Telomere_markers %>% filter(Sample %in% Phenotypes$ID) %>% arrange(Sample) -> Telomere_markers

Phenotypes %>% select(ID, AgeIncludingMonth, AgeInDays, Sex) %>% mutate(Age =AgeIncludingMonth)-> Cov
Phenotypes %>% select(-c(AgeIncludingMonth, AgeInDays, Sex)) -> Phenotypes

cell_types = colnames(Phenotypes)[grepl("10E9/L",colnames(Phenotypes))]
cell_types[!cell_types == "Leukocytes (10E9/L)" ] -> cell_types
c(cell_types,"Erythrocytes (10E12/L)") -> cell_types
left_join(Cov, select(Phenotypes, c("ID", cell_types))) -> Cov 
sapply(cell_types, FUN=function(x){paste(c("`",x,"`"), collapse="")}) -> cell_types


cell_types_PC %>% arrange(`-`) %>% filter(`-` %in% Cov$ID) -> cell_types_PC
all_cell_types %>% arrange(ID) %>% filter(ID %in% Cov$ID) -> all_cell_types

###########################
##Exploratory analysis#####
###########################
library(RColorBrewer)
library(wesanderson)
library(pheatmap)
left_join(mutate(Telomere_markers, ID = Sample), select(Cov,c("ID", "AgeIncludingMonth" , "Sex"))) %>% gather(Cell_line, Telomere_length, 2:7) %>% select(-Sample) %>% mutate(Telomere_length = as.numeric(Telomere_length)) -> Exploration_telomeres
#distribution of telomere length by gender
Exploration_telomeres %>% ggplot(aes(x=Telomere_length, fill=as.factor(Sex) )) + geom_histogram() + theme_bw() + facet_wrap(~ Cell_line) +
  scale_fill_manual(values = wesanderson::wes_palette("Royal1"),name="Sex",labels = c("Male", "Female")) -> Distribution_plot
SAVE(NAME =  "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Distribution_Telomeres.pdf",PLOT =Distribution_plot)
#Number of Males/Females
Exploration_telomeres %>% group_by(Sex,Cell_line) %>% summarise(n())
#correlation between telomere lengths
Exploration_telomeres %>% spread(Cell_line, Telomere_length) %>% drop_na() -> Exploration_telomeres_wide
cor(Exploration_telomeres_wide[4:9],method = "pearson") -> Correlation_matrix
pheatmap(Correlation_matrix, color=RColorBrewer::brewer.pal(name ="YlOrRd", n=9)) -> Pheatmap_correlation_Telomeres
SAVE(NAME = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Heatmap_Correlation_Telomeres.pdf",PLOT = Pheatmap_correlation_Telomeres)
pheatmap(abs(Correlation_matrix), color=RColorBrewer::brewer.pal(name ="YlOrRd", n=9)) -> Pheatmap_correlation_Telomeres_abs
SAVE(NAME = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Heatmap_Correlation_Telomeres_abs.pdf",PLOT = Pheatmap_correlation_Telomeres_abs)


#Correlation with Age and with Gender
Exploration_telomeres %>% ggplot(aes(x=AgeIncludingMonth, y=Telomere_length, col=as.factor(Sex) )) + geom_point() + theme_bw() + facet_wrap(~ Cell_line) +
  scale_color_manual(values = wesanderson::wes_palette("Royal1"),name="Sex",labels = c("Male", "Female")) + geom_smooth(method="lm", formula=y~x) + scale_x_continuous(name="Age") -> Age_vs_Telomere
SAVE(NAME = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Age_vs_Telomere.pdf",PLOT = Age_vs_Telomere)


model_age = function(x){ 
  Exploration_telomeres_wide[1:3] %>% mutate(Dependent=as_vector(x)) ->Input_m
  summary(lm(Dependent ~ AgeIncludingMonth + as.factor(Sex), Input_m)) -> Model
  as.data.frame(Model$coefficients) -> Param
  return(c(Param[2,1],Param[3,1],Param[2,4],Param[3,4] ))
}
apply(Exploration_telomeres_wide[,4:9], 2, FUN = model_age) %>% as_tibble() -> df
as_tibble(cbind(nms = names(df), t(df))) %>% `colnames<-`(c("Cell_line","Estimate_Age", "Estimate_Sex", "P_Ave", "P_sex")) -> df
apply(df[2:dim(df)[2]],2, FUN=as.numeric) %>% as_tibble() %>% mutate(ID = df$Cell_line) -> df
df[,c(1,3,5)]  %>% `colnames<-`(c("Estimate","Pvalue","Cell_line")) %>% mutate(Phenotype="Age") -> Age_df ; df[,c(2,4,5)]  %>% `colnames<-`(c("Estimate","Pvalue","Cell_line")) %>% mutate(Phenotype= "Sex") -> Sex_df
rbind(Age_df, Sex_df) %>% ggplot(aes(x=Cell_line , y=-log10(Pvalue) , fill = Phenotype, col=Estimate>0 )) + geom_bar(stat = "identity",position = "dodge",size=1) + theme_bw() + 
  scale_fill_manual(values = wesanderson::wes_palette("Royal1")[3:4]) + coord_flip() + geom_hline(yintercept = -log10(0.05)) +
  scale_color_manual(values = wesanderson::wes_palette("Royal1")[1:2]) -> Covariates_VS_telomeres
SAVE(NAME = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Covariates_vs_Telomere.pdf",PLOT = Covariates_VS_telomeres)

#Relation of Telomeres with other Biological aging values
Age_markers %>% select(c(Sample, dCT, Methyl.Age.Hannum.)) %>% mutate(sjTREC = as.numeric(dCT), ID=Sample) %>% select(-c(Sample,dCT)) -> Other_markers
left_join(Exploration_telomeres_wide, Other_markers, by="ID") %>% drop_na() %>% select(-Sex) -> Subset_all
cor(select(Subset_all, -ID),method = "pearson") -> Correlation_matrix
pheatmap(Correlation_matrix, color=rev(RColorBrewer::brewer.pal(name ="RdYlBu", n=9))) -> Corre_bioages
SAVE(NAME = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Correlation_BioAges.pdf",PLOT = Corre_bioages)

#Are samples aging too-fast or too-slow the same in the different bio-aging measurements?
Compare_outlaiers_bioage(Subset_all)

#Questionaires available
Study_phenos = function(Pheno){
  Pheno = as.numeric(Pheno)
  Pheno[is.na(Pheno)] -> MIA
  sum(length(MIA)) -> Number_NA
  Pheno[!is.na(Pheno)] -> Pheno
  max(Pheno) -> Maximum ; min(Pheno) -> Minimum ; sd(Pheno) -> SD ; mean(Pheno) -> Average_value ; median(Pheno) -> Median_value
  OUT = c(MAX=Maximum, MIN=Minimum, SD=SD, AVERAGE=Average_value, MEDIAN=Median_value, NA_number=Number_NA, N= length(Pheno) )
  return(OUT)
}

Phenotypes %>% filter(ID %in% Exploration_telomeres_wide$ID) %>% select(-ID) %>% apply(MARGIN=2, FUN=Study_phenos) %>% as_tibble() -> df
as_tibble(cbind(nms = names(df), t(df))) %>% `colnames<-`(c("Phenotype","Max","Min","SD", "Average","Median","NA", "Number") ) -> df
write_tsv(df, "/Users/sergio/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/Summary_phenotypes.tsv")


#############################
#####Phenotype association###
#############################


##Association without correcting for cell number
Model_phenotypes(select(Phenotypes,-"ID"), select(Telomere_markers, -"Sample"), Cov) %>% mutate(Pvalue=as.numeric(Pvalue)) -> Model_results
Model_results %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Model_results1
write_tsv(path = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Summary_stats/Phenos_vs_Telomere.tsv", x=Model_results)

Model_results %>% filter(FDR < 0.05)

##Association while controlling for cell number
Model_phenotypes(select(Phenotypes,-"ID"), select(Telomere_markers, -"Sample"), Cov, cells=T) %>% mutate(Pvalue=as.numeric(Pvalue)) -> Model_results
Model_results %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Model_results
write_tsv(path = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Summary_stats/Phenos_vs_Telomere_cellcorrected.tsv", x=Model_results)

left_join(select(Model_results, c(Regressor, Dependent, ,Beta, Pvalue, FDR, Levels_n)),select(Model_results1, c(Regressor, Dependent,Beta, Pvalue, FDR, Levels_n)), by=c("Regressor", "Dependent", "Levels_n"),suffix=c(".cellcorrec",".notcorrec") ) -> Table_results
Table_results %>% filter(FDR.cellcorrec < 0.05 | FDR.notcorrec<0.05 ) %>% arrange(Pvalue.cellcorrec) %>% mutate(Beta.cellcorrec=as.numeric(Beta.cellcorrec), Beta.notcorrec=as.numeric(Beta.notcorrec)) -> Table_results_sig
write_tsv(Table_results_sig,path = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Summary_stats/Phenos_vs_Telomere_allSig.tsv")

##Plot associatiations
library(UpSetR)
Model_results %>% filter(FDR<0.05) %>% select(Regressor, Dependent) %>% mutate(Value = 1) %>% spread(Dependent, Value) -> Upset_input
Upset_input[is.na(Upset_input)] = 0
as_tibble(cbind(nms = names(Upset_input), t(Upset_input))) -> Upset_input2 ; Upset_input2 %>% filter(! nms == "Regressor") %>% `colnames<-`(c("Regressor",Upset_input$Regressor)) ->Upset_input2
apply(select(Upset_input2, -Regressor),2, FUN=as.numeric) %>% as_tibble() %>% mutate(Regressor = Upset_input2$Regressor) -> Upset_input2
#upset(as.data.frame(Upset_input))
upset(as.data.frame(Upset_input2),sets = colnames(select(Upset_input2, -Regressor)), text.scale = 2,nsets = 10) -> upset_pheno_associations
Model_results %>% filter(FDR<0.05) %>% group_by(Dependent) %>% summarise(n())

pdf(file = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Summary_associations_phenoVStelomere.pdf",width = 6.3, height = 6.3)
upset_pheno_associations
dev.off()



#Check the signiciant hits.
# Are associations to BMI independent to parental smoking associations?
Phenotypes %>% select(`Waist circumference in cm`,`Body Mass Index (kg/M^2)`) -> Phenotype_BMI
Cov %>% select(-c(Age, AgeInDays)) %>% mutate(Parental_smoke = Phenotypes$Parental_smoking, Age_f = Phenotypes$age_f, Age_m=Phenotypes$age_m) -> Cov_BMI
select(Telomere_markers, -"Sample") -> Data_BMI
BMI_associations = tibble()
for (Col_n in 1:ncol(Phenotype_BMI)){
  Dependent = as.vector(as_vector((Phenotype_BMI[,Col_n])))
  Name = colnames(Phenotype_BMI)[Col_n]
  Name= clean_name(Name)
  for (Col_m in 1:ncol(Data_BMI)){
    Explanatory  = as.vector(as_vector(Data_BMI[,Col_m]))
    Name_2 = colnames(Data_BMI)[Col_m]
    Name_2= clean_name(Name_2)
    Cov_BMI %>% select(-ID) %>% mutate(Dependent = Dependent, Explanatory = as.numeric(Explanatory)) -> Model_input
    lm(Dependent ~ Explanatory + ., Model_input) -> Model

    Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Explanatory",])
    Beta = Summary$Estimate
    pvalue = Summary$`Pr(>|t|)`
    stand = Summary$`Std. Error`
    N = length(summary(Model)$residuals)
    Results = as_tibble(t(c(Name_2, Name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
    colnames(Results) = c("Regressor", "Dependent", "Beta", "SE", "Pvalue", "N")
    BMI_associations = rbind(BMI_associations, Results)
  }
}

inner_join(Model_results, BMI_associations, c("Regressor", "Dependent"), suffix=c("Regular", "Parentsmk_corrected")) -> Check_BMI
Check_BMI %>%
  ggplot(aes(x=as.numeric(BetaRegular), y=as.numeric(BetaParentsmk_corrected), col=Dependent)) + geom_point() + theme_bw() -> BMI_1

Check_BMI %>%
  ggplot(aes(x=-log10(as.numeric(PvalueRegular)), y=-log10(as.numeric(PvalueParentsmk_corrected)),col=Dependent)) + geom_point() + theme_bw() + geom_hline(yintercept = -log10(0.05)) -> BMI_2
BMI_1 + BMI_2 + plot_layout(guides = "collect")


# Are associations to Diabetes independent to parental age?
Phenotypes %>% select(Parents_diabetes2,`Have you ever had poorly healing wounds on your feet?`) -> Phenotype_Diabetes
Cov %>% select(-c(Age, AgeInDays)) %>% mutate(Parental_smoke = Phenotypes$Parental_smoking, Age_f = Phenotypes$age_f, Age_m=Phenotypes$age_m) -> Cov_Diabetes
select(Telomere_markers, -"Sample") -> Data_Diab
Diab_associations = tibble()
for (Col_n in 1:ncol(Phenotype_Diabetes)){
  Dependent = as.vector(as_vector((Phenotype_Diabetes[,Col_n])))
  Name = colnames(Phenotype_Diabetes)[Col_n]
  Name= clean_name(Name)
  for (Col_m in 1:ncol(Data_Diab)){
    Explanatory  = as.vector(as_vector(Data_Diab[,Col_m]))
    Name_2 = colnames(Data_Diab)[Col_m]
    Name_2= clean_name(Name_2)
    Cov_Diabetes %>% select(-ID) %>% mutate(Dependent = Dependent, Explanatory = as.numeric(Explanatory)) -> Model_input
    glm(Dependent ~ Explanatory + ., Model_input,family = binomial(link='logit')) -> Model
    
    Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Explanatory",])
    Beta = Summary$Estimate
    pvalue = Summary$`Pr(>|z|)`
    stand = Summary$`Std. Error`
    N = length(summary(Model)$residuals)
    Results = as_tibble(t(c(Name_2, Name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
    colnames(Results) = c("Regressor", "Dependent", "Beta", "SE", "Pvalue", "N")
    Diab_associations = rbind(Diab_associations, Results)
  }
}


inner_join(Model_results, Diab_associations, c("Regressor", "Dependent"), suffix=c("Regular", "ParentalAge_corrected")) -> Check_Diab
Check_Diab %>%
  ggplot(aes(x=as.numeric(BetaRegular), y=as.numeric(BetaParentalAge_corrected),col=Dependent)) + geom_point() + theme_bw() -> Diabetes_1

Check_Diab %>%
  ggplot(aes(x=-log10(as.numeric(PvalueRegular)), y=-log10(as.numeric(PvalueParentalAge_corrected)), col=Dependent)) + geom_point() + theme_bw() + geom_hline(yintercept = -log10(0.05)) -> Diabetes_2
Diabetes_1 + Diabetes_2 + plot_layout(guides = "collect")

#Show the effect of diabetes
Cov %>% mutate(Telomere=Telomere_markers$`MTL_CD20+`, Diabetes = Phenotypes$Parents_diabetes2) %>% ggplot(aes(x=as.factor(Diabetes), y=Telomere)) + geom_boxplot() + geom_point()
Cov %>% mutate(Telomere=Telomere_markers$MTL_gran, Diabetes=Phenotypes$type2diabetes) %>% ggplot(aes(x=as.factor(Diabetes), y=Telomere)) + geom_boxplot() + geom_point()
#diabetes
Cov %>% mutate(Telomere=Telomere_markers$MTL_gran, smk= Phenotypes$smk13, BMI = Phenotypes$`Body Mass Index (kg/M^2)`, Diabetes=Phenotypes$type2diabetes, Glucose = Phenotypes$`HbA1c (mmol/mol)`) -> Cov
#diabetes_parents
Cov %>% mutate(Telomere=Telomere_markers$`MTL_CD20+`,smk= Phenotypes$smk13, BMI = Phenotypes$`Body Mass Index (kg/M^2)`, Diabetes = Phenotypes$Parents_diabetes2,Age_father=Phenotypes$age_f,Age_mother=Phenotypes$age_m, Glucose = Phenotypes$`HbA1c (mmol/mol)`) -> Cov
Cov %>% mutate(Age_father = (Age_father + Age), Age_mother = (Age_mother + AgeIncludingMonth)) -> Cov
#Cov$Diabetes[is.na(Cov$Diabetes)] = 0
summary(glm(formula = as.factor(Diabetes) ~ Telomere + Age_father + Age_mother+ AgeIncludingMonth + Sex ,data=Cov, family = binomial(link="logit") ))

lm(formula= Telomere ~ as.factor(Diabetes) + AgeIncludingMonth + Sex, data=Cov) -> check_model

#Is The polygenenic risk score of diabetes correlatd with telomeres? Common genetic architecture
PRS_BMI = as_tibble(read.table("Diabetes_risk_score_BMIcorrected.tsv", sep="", header=T))
PRS = as_tibble(read.table("Diabetes_risk_score.tsv", sep="", header=T))

colnames(PRS_BMI)[2] = "ID" ; colnames(PRS)[2] = "ID" 
left_join(Cov, select(PRS_BMI, c("ID","SCORESUM"))) %>% mutate(Diabetes_risk_bmi = SCORESUM ) %>% select(-SCORESUM) -> Cov_diabetes
left_join(Cov_diabetes, select(PRS, c("ID","SCORESUM")), by="ID") %>% mutate(Diabetes_risk = SCORESUM ) -> Cov_diabetes
summary(lm(formula= Telomere ~  Diabetes_risk + AgeIncludingMonth + Sex, data=Cov_diabetes)) -> check_model
summary(lm(formula= Telomere ~  Diabetes_risk_bmi + AgeIncludingMonth + Sex, data=Cov_diabetes)) -> check_model

PRS_BMI %>% arrange(ID) %>% filter(ID %in% Telomere_markers$Sample) %>% mutate(PRS_BMI = SCORESUM) %>% select(PRS_BMI) -> Phenotypes_PRS
cbind(Phenotypes_PRS,PRS %>% arrange(ID) %>% filter(ID %in% Telomere_markers$Sample) %>% mutate(PRS = SCORESUM) %>% select(PRS,ID) ) %>% as_tibble -> Phenotypes_PRS

Model_phenotypes(select(Phenotypes_PRS,-"ID"), select(filter(Telomere_markers, Sample %in% Phenotypes_PRS$ID), -"Sample"), filter(Cov,ID %in%Phenotypes_PRS$ID) , cells=T) %>% arrange(Pvalue)


#glucosa levels between diabetes/non diabetes
lm(Glucose ~ Diabetes + AgeIncludingMonth,Cov) 
###Association between phenotypes?
Phenotypes %>% select(`Body Mass Index (kg/M^2)`,Parental_smoking) %>% drop_na() -> Subhypothesis
cor.test(Subhypothesis$`Body Mass Index (kg/M^2)`, Subhypothesis$Parental_smoking,method = "spearman")

Phenotypes %>% select(`Do your parents, siblings or children have diabetes mellitus type 2, or have they had diabetes mellitus type 2? No.`,age_f, age_m) %>% drop_na() -> Subhypothesis
cor.test(Subhypothesis$`Do your parents, siblings or children have diabetes mellitus type 2, or have they had diabetes mellitus type 2? No.`, Subhypothesis$age_f,method = "spearman")
cor.test(Subhypothesis$`Do your parents, siblings or children have diabetes mellitus type 2, or have they had diabetes mellitus type 2? No.`, Subhypothesis$age_m,method = "spearman")

####

#How do the telomeres associate to all cell  counts
Cov %>% select(-c(AgeIncludingMonth, AgeInDays)) %>% filter(ID %in% all_cell_types$ID) -> Cov_2
Telomere_markers %>% filter(Sample %in% all_cell_types$ID) %>% select(-"Sample") -> Data
cell_associations = tibble()

for (Col_n in 1:(ncol(all_cell_types)-1)){
  Dependent = as.vector(as_vector((all_cell_types[,Col_n])))
  Name = colnames(all_cell_types)[Col_n]
  Name= clean_name(Name)
  for (Col_m in 1:ncol(Data)){
    Explanatory  = as.vector(as_vector(Data[,Col_m]))
    Name_2 = colnames(Data)[Col_m]
    Name_2= clean_name(Name_2)
    print(c(Name,Name_2))
    Cov_2 %>% select(-ID) %>% mutate(Dependent = Dependent, Explanatory = as.numeric(Explanatory)) -> Model_input
    lm(Dependent ~ Explanatory + ., Model_input) -> Model
    Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Explanatory",])
    Beta = Summary$Estimate
    pvalue = Summary$`Pr(>|t|)`
    stand = Summary$`Std. Error`
    N = length(summary(Model)$residuals)
    Results = as_tibble(t(c(Name_2, Name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
    colnames(Results) = c("Regressor", "Dependent", "Beta", "SE", "Pvalue", "N")
    cell_associations = rbind(cell_associations, Results)
  }
}
cell_associations %>% arrange(Pvalue) %>% mutate(FDR = p.adjust(as.numeric(Pvalue), "fdr")) %>% ggplot() + geom_point(aes(x=as.numeric(Beta), y=-log10(as.numeric(Pvalue)), col=FDR<0.05)) + theme_bw()
cell_associations %>% arrange(Pvalue) %>% mutate(FDR = p.adjust(as.numeric(Pvalue), "fdr")) -> Results_cell
write_tsv(x = Results_cell, path = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Phenotype_summaries/Cell_counts.tsv")
##

##Check correlations of smoking phenotypes
Model_results %>% filter(grepl("smk",Dependent) ) %>% ggplot(aes(x= Dependent, y=as.numeric(Beta), fill=-log10(Pvalue), col=Pvalue<0.05 )) + geom_bar(size=1, stat = "identity") +
  theme_bw() + facet_wrap(~Regressor) + coord_flip() + scale_color_manual(values = c("white", "black")) + scale_fill_distiller(palette="YlOrRd" ,direction=1) +
  scale_x_discrete(labels=c("smk1" = "Have_you_smoked_a_year?", "smk2" ="How_old_when_started_smoking?","smk13"= "father_smk", "smk14" = "mother_smk", "smk15"="mother_smk_pregnancy","smk10"="Number_people_smoke_at_home","smk11"= "Do_people_smoke_work","smk_now"="smk_participant")) -> Smoking_Associations
SAVE(NAME = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Smoking_associations.pdf",PLOT=Smoking_Associations)




#compare effect size of boys and girls
Model_results %>% filter(grepl("smk",Dependent) ) %>% filter(Dependent %in% c("smk13", "smk14")) -> Man_vs_woman
summary(lm(Beta ~ Dependent, Man_vs_woman))
left_join(filter(Man_vs_woman,Dependent=="smk13"),filter(Man_vs_woman,Dependent=="smk14"), by="Regressor") %>% ggplot(aes(x=as.numeric(Beta.x), y=as.numeric(Beta.y))) + geom_point()+theme_bw() +geom_abline()



#####Combine all phenotypes and calculate total variance explained
set.seed(2020)
R2_calc = function(Real,pred){
  rss = sum((pred - Real)^2)
  tss = sum((Real - mean(Real))^2)
  rsq = 1 - rss/tss
  return(rsq)
}
#Add telomere Polygenic risk score
PRS_telomeres = read_tsv("PRS_telomere.tsv")
PRS_telomeres %>% filter(ID %in% Phenotypes$ID) -> PRS_totest
Variance_explained = tibble()
clean_name(colnames(Phenotypes)) -> colnames(Phenotypes)
Phenotypes %>% select(c(unique(filter(Model_results,FDR<0.05)$Dependent),ID)) -> Phenotypes_significant
for (Cell_line in colnames(Telomere_markers)){
  if (Cell_line == "Sample"){ next }
  Telomere_markers %>% select(Cell_line) %>% as_vector() %>% as.vector() %>% as.numeric() -> Dependent
  #PRS_totest %>% select(ID ,Cell_line) -> PRS_cell
  select(Cov, -c(AgeInDays, Age))[,1:10] -> Cov_to_add
  Phenotypes_significant %>% mutate(Dependent = Dependent) -> Entry_model
  #Remove colinear
  left_join(Entry_model, Cov_to_add, by="ID") -> Entry_model
  
  #Entry_model %>% select(-c(smk13, age_m,Waist_circumference_in_cm)) -> Entry_model
  #lm(Dependent ~ ., select(Entry_model, -ID)) -> Model_summary
  #Complete_model = summary(Model_summary)$adj.r.squared
  #lm(Dependent ~ ., select(Entry_model, -one_of(c(ID, colnames(Phenotypes_significant) )))) -> Model0_summary
  #Null_model = summary(Model0_summary)$adj.r.squared

  Entry_model %>% drop_na() -> Entry_model
  cvfit = cv.glmnet(y=as.vector(Entry_model$Dependent), x=as.matrix(select(Entry_model,-c(Dependent,ID))), nfolds = 10,type.measure="mse",standardize=T, alpha=1)
  Param_best = coef(cvfit, s = "lambda.min")
  Estimation = predict(cvfit, as.matrix(select(Entry_model,-c(Dependent,ID))), s="lambda.min")
  r2 = R2_calc(as.vector(Entry_model$Dependent), Estimation)
  
  as.matrix(Param_best) %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>%  filter(! `1` == 0) %>% select(rowname) %>% as_vector() %>% as.vector() -> Keep_covariates
  Entry_model %>% select(one_of(c("ID",Keep_covariates))) -> Telomere_param_used
  write_tsv(x = Telomere_param_used, path = paste(c("Param_chosen/covariates_",Cell_line,".tsv"),collapse=""))
  
  cvfit_null = cv.glmnet(y=as.vector(Entry_model$Dependent), x=as.matrix(select(Entry_model,-one_of(c("ID", colnames(Phenotypes_significant),"Dependent" )))), nfolds = 10,type.measure="mse",standardize=T, alpha=1)
  Param_best_null = coef(cvfit_null, s = "lambda.min")
  Estimation_null = predict(cvfit_null, as.matrix(select(Entry_model,-one_of(c("ID", colnames(Phenotypes_significant),"Dependent" )))), s="lambda.min")
  r2_null = R2_calc(as.vector(Entry_model$Dependent), Estimation_null)
  
  lm(Dependent ~ AgeIncludingMonth, Entry_model) ->  Null_0
  R2_null0 = R2_calc(as.vector(Entry_model$Dependent), predict(Null_0))
  
  #colnames(PRS_cell)[2] = "PRS"
  #left_join(Entry_model, PRS_cell) -> Complete_model
  #Complete_model %>% drop_na() -> Complete_model
  #cvfit = cv.glmnet(y=as.vector(Complete_model$Dependent), x=as.matrix(select(Complete_model,-c(Dependent,ID))), nfolds = 10,type.measure="mse",standardize=T, alpha=1)
  #Param_best = coef(cvfit, s = "lambda.min")
  #Estimation = predict(cvfit, as.matrix(select(Complete_model,-c(Dependent,ID))), s="lambda.min")
  #r2_prs = R2_calc(as.vector(Complete_model$Dependent), Estimation)
  
  rbind(Variance_explained,tibble(Telomere=Cell_line, R2_environment= r2, R2_Covariates = r2_null, R2_age=R2_null0, R2_genetics=r2_prs)) -> Variance_explained
  #rbind(Variance_explained,tibble(Telomere=Cell_line, R2_environment= r2, R2_Covariates = r2_null, R2_age=R2_null0, R2_genetics=r2_prs)) -> Variance_explained
    
}


Variance_explained %>% mutate(R2_genet9cs = R2_genetics -R2_environment ,R2_environment = R2_environment-R2_Covariates, R2_Covariates = R2_Covariates-R2_age) %>% gather("Model","r2", c("R2_genetics","R2_environment","R2_Covariates","R2_age")) %>% mutate(Model = fct_relevel(Model,c("R2_genetics","R2_environment", "R2_Covariates", "R2_age"))) %>%
  ggplot(aes(x=Telomere, y=r2, fill=Model)) + geom_bar(position="stack", stat="identity") +
  coord_flip() + theme_bw() + scale_fill_manual(values = wesanderson::wes_palette("Royal1")) 



################################
##Analysis with Proteomics######
################################

Protein_covariates %>% arrange(ID) %>% filter(ID %in% Telomere_markers_protein$Sample) -> Cov_protein
Telomere_markers %>% filter(Sample %in% Cov_protein$ID) -> Telomere_markers_protein
#Cov %>% filter(ID %in% Proteomics$ID) -> Cov_protein
Proteomics %>% filter(ID %in% Cov_protein$ID) %>% arrange(ID) -> Proteomics




##Model2
Model_phenotypes(select(Proteomics, -"ID"), select(Telomere_markers_protein, - "Sample"), Cov_protein,cells = T, Proteomics = T) %>% mutate(Pvalue=as.numeric(Pvalue)) -> Model_results_olink
Model_results_olink %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Model_results_olink1
write_tsv(path = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Summary_stats/Telomere_Olink_CellCorrecter.tsv", x=Model_results_olink1)

Cov_protein %>% mutate(BMI = filter(arrange(Phenotypes,ID), ID %in% Cov_protein$ID)$`Body Mass Index (kg/M^2)`) -> Cov_protein2
Model_phenotypes(select(Proteomics, -"ID"), select(Telomere_markers_protein, - "Sample"), Cov_protein2, Proteomics = T) %>% mutate(Pvalue=as.numeric(Pvalue)) -> Model_results_olink_BMI
Model_results_olink_BMI %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Model_results_olink_BMI


#Model_phenotypes(select(Proteomics, -"ID"), select(Telomere_markers_protein, - "Sample"), Cov_protein) %>% mutate(Pvalue=as.numeric(Pvalue)) -> Model_results_olink
#Model_results_olink %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Model_results_olink2
#write_tsv(path = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Summary_stats/Telomere_Olink.tsv", x=Model_results_olink2)


Model_results_olink1 %>% filter(FDR<0.05)

###are RARRES and FABP4 concounded by BMI?
left_join(Cov_protein, select(Phenotypes, c("ID","Body Mass Index (kg/M^2)"))) -> Cov_with_BMI
Cov_with_BMI %>%  mutate(RARRES2 = Proteomics$RARRES2, FABP4 = Proteomics$FABP4) -> Input_model
corrected_stats = tibble()
for (Cell_line in colnames(select(Telomere_markers_protein, -Sample))){
  Telomere_markers_protein %>% select(Cell_line) %>% as_vector() %>% as.vector() -> Regressor
  Input_model %>% mutate(Cell_line=as.numeric(Regressor)) -> Input_model
  paste(c("`",paste(colnames(Input_model)[c(2,4,6:13)], collapse = "`+`"), "`"), collapse="") -> Covariates_BMI
  as.formula(paste(c("RARRES2 ~ as.numeric(Regressor) + ",Covariates_BMI),collapse="")) -> formula1
  as.formula(paste(c("FABP4  ~ as.numeric(Regressor) + ",Covariates_BMI),collapse="")) -> formula2
  
  summary(lm(formula1, data=Input_model)) -> Results_model
  as.data.frame(Results_model$coefficients)["as.numeric(Regressor)", c("Estimate","Pr(>|t|)")] %>% as_tibble() %>% mutate(Cell_line =Cell_line, Protein="RARRES2") -> Tibble_r
  rbind(corrected_stats,Tibble_r) -> corrected_stats
  summary(lm(formula2, data=Input_model)) -> Results_model
  as.data.frame(Results_model$coefficients)["as.numeric(Regressor)", c("Estimate","Pr(>|t|)")] %>% as_tibble() %>% mutate(Cell_line =Cell_line, Protein="FABP4") -> Tibble_r
  rbind(corrected_stats,Tibble_r) -> corrected_stats
  
}

Model_results_olink1 %>% filter(Dependent %in% c("RARRES2", "FABP4")) %>% arrange(Regressor, Dependent) -> uncorrected_stats

corrected_stats %>% arrange(Cell_line, Protein) ->  corrected_stats
cbind(uncorrected_stats, corrected_stats) %>% as_tibble() -> plot_input

plot_input %>% ggplot(aes(x=as.numeric(Beta), y = as.numeric(Estimate), col= Protein)) + geom_point() + theme_bw()+ geom_abline()
plot_input %>% ggplot(aes(x=-log10(as.numeric(Pvalue)), y = -log10(as.numeric(`Pr(>|t|)`)), col= Protein)) + geom_point() + theme_bw()+ geom_abline()

Model_phenotypes(select(Proteomics, -"ID"), select(Telomere_markers_protein, - "Sample"), Cov_with_BMI,cells = T, BMI=T) %>% mutate(Pvalue=as.numeric(Pvalue)) -> Model_results_olink
Model_results_olink %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Model_results_olink3
write_tsv(path = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Summary_stats/Telomere_Olink_CellBMI_corrected.tsv", x=Model_results_olink3)

##what about TFPI and smoking
Cov %>% mutate(smk_now = Phenotypes$smk_now) %>% filter(ID %in% Proteomics$ID) -> Cov_TFPI
colnames(Cov_TFPI) = clean_name(colnames(Cov_TFPI))
Cov_TFPI %>% mutate(Dependent = filter(Proteomics, ID %in% Cov_TFPI$ID)$TFPI) ->  model_input
MODEL = paste(c("Dependent ~ Sex + AgeIncludingMonth +", paste(colnames(Cov_TFPI)[6:13], collapse="+"), "+Telomere"),collapse="") 
Output_TFPI = tibble()
for (cell in colnames(Telomere_markers_protein)){
    if(cell == "Sample"){next}
  Telomere_markers_protein %>% select(cell) %>% as_vector() %>% as.numeric() -> Telo_length
  model_input %>% mutate(Telomere = Telo_length) -> model_input2
  as.data.frame(summary(lm(MODEL,model_input2))$coefficients)["Telomere",] %>% as_vector() %>% as.vector() -> out_model
  tibble(Regressor=cell, Dependent="TFPI", Beta=out_model[1], Pvalue=out_model[4]) -> out_model
  colnames(out_model) = c("Regressor","Dependent", "Beta", "Pvalue")
  Output_TFPI = rbind(Output_TFPI, out_model)
}
Output_TFPI
Model_results_olink %>% mutate(Beta = as.numeric(Beta), Pvalue=as.numeric(Pvalue)) -> Model_results_olink
inner_join(Model_results_olink, Output_TFPI, by=c("Regressor","Dependent")) %>% ggplot(aes(x=Beta.x, y = Beta.y)) + theme_bw() + geom_point()
inner_join(Model_results_olink, Output_TFPI, by=c("Regressor","Dependent")) %>% ggplot(aes(x=Pvalue.x, y = Pvalue.y)) + theme_bw() + geom_point()

##




##Heritability
for (Cell_line in colnames(Telomere_markers, -ID)){
  Telomere_markers %>%mutate(ID=Sample) %>% select("ID", Cell_line) -> Telomere_entry
  colnames(Telomere_entry)[2] = "Dependent"
  PRS_totest  %>% select("ID", Cell_line) -> PRS_entry
  colnames(PRS_entry)[2] = "Regressor"
  Cov %>% select(ID,AgeIncludingMonth, Sex) -> Cov_entry
  left_join(left_join(Telomere_entry, PRS_entry), Cov_entry) -> H_input
  H_input %>% drop_na() -> H_input
  summary(lm(as.numeric(Dependent) ~ as.numeric(Regressor) + AgeIncludingMonth + Sex ,select(H_input, -ID)))
}
heritbility_summary = read_tsv("Heritability_stuff/Heritability_total.tsv")
heritbility_summary %>% select(-c(Pvalue, h2_LD, Sd_LD)) %>% mutate(method="GREML") -> GREML
heritbility_summary %>% select(-c(Pvalue, h2, Sd)) %>% mutate(method = "LD") -> LD_regression
colnames(LD_regression) = colnames(GREML)
rbind(GREML, LD_regression) -> heritability_summary2

heritability_summary2 %>% mutate(Error = ifelse((h2+Sd) > 1, 1, h2+Sd )) %>%  ggplot(aes(x = Cell_line, y = h2, fill=method)) + geom_bar(stat="identity",position = "dodge",color="black") + coord_flip() +
  geom_errorbar(aes(ymin=h2, ymax=Error), width=.2,position=position_dodge(.9)) + theme_bw() + scale_fill_manual(values = wesanderson::wes_palette("Royal1")[1:2]) -> Heritability_plot
heritability_summary2 %>% group_by(method) %>% summarise(median(h2))                                                                                                                                                        
SAVE(NAME = "~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Figures/Heritability.pdf" ,PLOT = Heritability_plot)



#########################
#Prediction of mortality#
#########################
#Prediction (association) of death 
dead = as_tibble(read.table("~/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/Survival_/LLD_passed.txt",header = T))
Info = read_tsv("~/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/Survival_/baseline.txt")
End_date = as.Date("2020-05-01")
dead   %>% mutate (DATE = as.Date(paste(deathyear,deathmonth,"01",sep="-") )) -> dead


Info %>% mutate(First_time = as.Date(paste(bl1year,"01-01",sep="-")) ) -> Info #,Second_time = ifelse(LLDEEPID %in% dead$LLDeep_SampleID), dead[LLDeep_SampleID==LLDEEPID,]$DATE,End_date-First_time)

Second_time = c()
for (ID in Info$LLDEEPID){
  Info[Info$LLDEEPID==ID,]$First_time -> First_time
  if(ID %in% dead$LLDeep_SampleID){
    dead[dead$LLDeep_SampleID==ID,]$DATE -First_time -> Datet
  } else { End_date-First_time -> Datet }
  Second_time = c(Second_time, Datet)
}
Info %>% mutate(Second_time = Second_time, Status = ifelse(LLDEEPID %in% dead$LLDeep_SampleID, 1, 0)) -> Info

library("survival")
library("ggfortify")

Cov %>%  select(-c(AgeInDays,Age)) -> Cov
Predictors = Telomere_markers

Cov %>% filter(ID %in% Info$LLDEEPID) %>% arrange(ID) -> Cov
Predictors %>% filter(Sample %in% Info$LLDEEPID) %>% arrange(Sample) -> Predictors

Results = tibble()
for( Phenotype in colnames(select(Predictors,-c("Sample")))){
  print(Phenotype)
  Regressor = Predictors %>% select(Phenotype) %>% as_vector() %>% as.numeric()
  Info %>% filter(LLDEEPID %in% Predictors$Sample) -> Info2 ;  Cov %>% filter(ID  %in% Info2$LLDEEPID) -> Cov2
  cbind(Info2 %>% select(Second_time,Status), select(Cov2, - ID)) %>% as_tibble() %>% mutate(Pheno=Regressor) ->  All_info
  All_info %>% group_by(Status) %>% summarise(n()) %>% print()
  coxph(Surv(Second_time, Status) ~ ., data=All_info) -> res.cox
  autoplot(survfit(res.cox), surv.linetype = 'dashed', surv.colour = 'blue',
           conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.5, censor = FALSE)
  as.data.frame(summary(res.cox)$coefficients)["Pheno",] %>% as_tibble() %>% mutate(Regressor = Phenotype) -> res.cox
  rbind(res.cox, Results) -> Results
}
arrange(Results, `Pr(>|z|)`)
ggplot(Results) + geom_bar(aes(x=Regressor, y= -log10(`Pr(>|z|)`) , fill = `exp(coef)` > 1), stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept = -log10(0.05))


##Plot presence
library(UpSetR)
Info_avail = read.table("~/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/Genetic_and_methylation_availability.tsv",header = T)
Info_avail[! Info_avail==0] = 1 ; Info_avail[Info_avail==0] = 0 ; Info_avail[is.na(Info_avail)] = 0
apply(as.data.frame(Info_avail), 2, as.numeric) -> Info_avail
colnames(Info_avail) = c("LLDeep", "Methylation", "Genetics", "Onlink", "Telomere")
upset(as.data.frame(Info_avail),nsets = 10,text.scale = 2)


Info_avail = as_tibble(read.table("~/Resilio Sync/LLD phenotypes for Sergio/Phenotypes/Genetic_and_methylation_availability.tsv",header = T))

as_tibble(cbind(nms = names(Info_avail), t(Info_avail))) -> Info_avail
colnames(Info_avail) = Info_avail[1,]
Info_avail[2:5,] -> Info_avail
