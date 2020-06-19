library(tidyverse)



Read_data = function(){

  #Read phenotypes. Many record have a "Record" written in the name. Remove those and arrange by ID
  Phenotypes = read_csv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Merged_Phenotypes_LLD.csv")
  Rm = colnames(Phenotypes)[grepl("Record", colnames(Phenotypes))]
  Phenotypes %>% select(-Rm) -> Phenotypes 
  Phenotypes %>% arrange(ID) -> Phenotypes
  #Extra phenos
  How_are_you_Q = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Phenotype_Howareyou.tsv")
  How_are_you_Q %>% filter(ID %in% Phenotypes$ID) %>% arrange(ID) -> How_are_you_Q
  left_join(Phenotypes, filter(How_are_you_Q, Time == "Baseline assessment (1A)"), "ID") -> Phenotypes 


  #Read the age in months file, arrange by LLDEEPID
  Age_pheno = read_tsv(file = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Age.txt")
  Age_pheno %>% arrange(LLDEEPID) -> Age_pheno
  #Read the biological age file, arrange, remove records without age. Arrange by sample
  Data = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Combined_data_uncorrected.tsv")
  Data %>% arrange(Sample) %>% filter(Sample %in% Age_pheno$LLDEEPID) -> Data
  #Read proteomics data, transpose and arrange by ID
  Proteomics = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt")
  select(Proteomics, -Protein) %>% t() %>% as_tibble() %>% mutate(ID= colnames(Proteomics)[2:length(colnames(Proteomics))]) %>%
    `colnames<-`(c(Proteomics$Protein, "ID")) %>% arrange(ID) -> Proteomics2
  Calories = read_delim("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/calories__Diet_LLD.csv", delim=";")
  Calories %>% arrange(UMCGIBDResearchIDorLLDeepID) %>% filter(UMCGIBDResearchIDorLLDeepID %in% Phenotypes$ID) -> Calories
  #Change the answer of an ID who doesnt know how to read
  Phenotypes %>% mutate(`Do you know the outcome of these urine tests? Protein in the urine.(SZ:selfreported, 1-yes, 0-no orNA or nr)!!.` = ifelse(ID=="LLDeep_1422", "1", `Do you know the outcome of these urine tests? Protein in the urine.(SZ:selfreported, 1-yes, 0-no orNA or nr)!!.`)) -> Phenotypes
  colnames(Data)[grepl("INT",colnames(Data))] -> Remove_INT
  Data %>% select(-Remove_INT) -> Data
  
  Combine_fm_smoke = function(Phenotypes){
    Parent_smoke = vector()
    for (ID in 1:length(Phenotypes$fathersmk)){
      if( is.na(Phenotypes$fathersmk[ID]) | is.na(Phenotypes$mothersmk[ID]) ){
        Parent_smoke = c(Parent_smoke, NA)
      } else if (Phenotypes$fathersmk[ID] == 1 | Phenotypes$mothersmk[ID] == 1){
        Parent_smoke = c(Parent_smoke, 1)
      } else{
        Parent_smoke = c(Parent_smoke, 0)
      }
    }
    return(Parent_smoke)
  }
  
  Phenotypes %>% mutate(Parent_smoker = Combine_fm_smoke(Phenotypes)) -> Phenotypes
  Phenotypes %>% filter(! ID == "#N/A") -> Phenotypes
  
  #Dasha's Age file
  LLD_age_smk = read_tsv("/Users/Sergio/Resilio Sync/LLD phenotypes for Sergio/Results/aging_vs_expression_methylation/parental_age_smk_final.txt")
  LLD_age_smk %>% filter(LLDEEPID %in% Phenotypes$ID) %>% arrange(LLDEEPID) -> LLD_age_smk
  Phenotypes %>% mutate(AGE_FATHER == LLD_age_smk$age_f, AGE_MOTHER = LLD_age_smk$age_m) -> Phenotypes
  
  return(list(Phenotypes, Data, Age_pheno, Calories, Proteomics))
}  
  


Read_data() -> List_inputs
Phenotypes <- List_inputs[[1]] ; Data <- List_inputs[[2]] ; Age_pheno <- List_inputs[[3]] ; Calories <- List_inputs[[4]] ; Proteomics <- List_inputs[[5]]

D = as.matrix(select(Data, -c(Sample,Age)))
D[is.na(D)] = 0
D[D != 0] = 1
apply(D, 2,as.numeric) -> D

library(UpSetR)


as_tibble(D) %>% mutate(Sample = Data$Sample) -> D
upset(as.data.frame(D),sets= colnames(select(D, -Sample)),  text.scale = c(1.3, 1, 1, 1, 1, 2))
   
as_tibble(D) %>% mutate(Sample = Data$Sample) %>% filter(`Methyl.Age.Hannum.` == 1 & `Methyl.Age.Weidner.` == 0)  %>% select(Sample)                                                    


##Matching data

#Filter Phenotypes so that there are no individuals who lack biological age measures
Phenotypes %>% filter(ID %in% Calories$UMCGIBDResearchIDorLLDeepID) -> Phenotypes #%>% mutate(Calorie = Calories$SUMOFKCAL) -> Phenotypes
Phenotypes = cbind(Phenotypes, select(Calories, -UMCGIBDResearchIDorLLDeepID))
Phenotypes %>% filter(ID %in% Data$Sample) %>% filter(ID %in% Age_pheno$LLDEEPID) -> Phenotypes

#Filter Age pheno so that they are in the filtered biologicala ge
Age_pheno %>% filter(LLDEEPID %in% Phenotypes$ID) -> Age_pheno
#Create covariates: Age, Sex, ID of individuals with phenotypes, data and age
Covariates = tibble(Age = Age_pheno$AgeIncludingMonth, Sex = Phenotypes$Sex, ID = Age_pheno$LLDEEPID)
#Removed phenotupes we do not want to test
Phenotypes %>% select(-c(ID,Sex,Age)) -> Phenotypes 

#Counts of smokers
Phenotypes %>% group_by(Parent_smoker, smoking) %>% summarise(n()) %>% drop_na() %>% mutate(Group = paste(Parent_smoker, smoking, sep="-")) %>% ggplot(aes(x=Group, y = `n()`)) + geom_bar(stat="identity") + theme_bw()




#####Age vs Bilogical age: examples
Data %>% select(-Age) %>% mutate(Age = Covariates$Age) %>% gather(Biological_age, Measure, 2:11) %>% mutate(Measure = as.numeric(Measure)) %>%
ggplot(aes(y=Measure, x=Age))+ geom_point() + theme_bw()  + geom_smooth(se=F, method=lm) + facet_wrap(~Biological_age, scales="free")


#Function that outputs the residuals of removing "Sex" and removing or not Age
Get_residuals = function(){
    #Go through each variable in Data, adjust to models in which Sex (and Age) are the regressor.
    #Get residuals
    Results_c = tibble(ID = Data$Sample)
    Results_u = tibble(ID = Data$Sample)
    for (Col_n in 2:ncol(Data)){
      Name = colnames(Data)[Col_n]
      print(Name)
      #Name= clean_name(Name)
      
      Dependent = as.numeric(as.vector(as_vector((Data[Col_n]))))
      
      Model_data = tibble(Dependent=Dependent, Sex= factor(Covariates$Sex), Age= as.numeric(Covariates$Age), ID = Data$Sample)
      Model_data %>% drop_na() -> Model_data
      Model = lm(as.numeric(Dependent) ~ as_factor(Sex) + as.numeric(Age), Model_data)
      Model_un = lm(as.numeric(Dependent) ~ as_factor(Sex), Model_data)
      
      Residuals_age = summary(Model)$residuals
      Residuals_noage = summary(Model_un)$residuals
      
      Corrected = tibble(ID = Model_data$ID, Residuals_Age = Residuals_age)
      colnames(Corrected) = c("ID", Name)
      
      UnCorrected = tibble(ID = Model_data$ID, Residuals_noage = Residuals_noage)
      colnames(UnCorrected) = c("ID", Name)
      #Residuals are added as columns to the rows with the same IDs
      Results_c = left_join(Results_c, Corrected)
      Results_u = left_join(Results_u, UnCorrected)  
        }
    #Save residuals  
    write_csv(x = Results_c, path = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/BiologicalAge_Sex_Age.tsv")
    write_csv(x = Results_u, path = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/BiologicalAge_Sex.tsv")

}
#Functions for correlations of PHENO  ~ Bio_Age + COV_AGE + COV_SEX
Linear_regression = function(Dependent, Regressor,Covariates, Regressor_name, Dependent_name, Plot = F, wAge ="yes"){
  Regressor = as.numeric(Regressor)
  Model_data = tibble(Dependent=Dependent, Regressor=Regressor, Sex= factor(Covariates$Sex), Age= as.numeric(Covariates$Age))
  Model_data %>% drop_na() -> Model_data
  if (wAge == "No"){
    Model = lm(Dependent ~Regressor + Sex, Model_data )
  }else{
  Model = lm(Dependent ~Regressor + Sex + Age, Model_data )
  }
  Summary = as_tibble(summary(Model)$coefficients)
  
  Beta = Summary$Estimate[2]
  pvalue = Summary$`Pr(>|t|)`[2]
  stand = Summary$`Std. Error`[2]
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
Logistic_regression = function(Dependent, Regressor,Covariates, Regressor_name, Dependent_name, Plot=T){
  Regressor = as.numeric(Regressor)
  #If instead of a binary 0-1 it is 1-2, make it 0-1
  if (length(which(Dependent == 2)) > 0){ Dependent = as.numeric(as.factor(Dependent)) - 1}
  Model_data = tibble(Dependent=Dependent, Regressor=Regressor, Sex= factor(Covariates$Sex), Age= as.numeric(Covariates$Age))
  Model_data %>% filter(! is.na(Dependent)) %>% filter(! is.na(Regressor)) -> Model_data
  Model = glm(Dependent ~Regressor + Sex + Age, Model_data, family=binomial(link='logit'))
  Summary = as_tibble(summary(Model)$coefficients)
  Beta = Summary$Estimate[2]
  pvalue = Summary$`Pr(>|z|)`[2]
  stand = Summary$`Std. Error`[2]
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

#Iteration through Phenotypes, in each phenotype adjust a model per Biological Age after adjusting for Covariates
Model_phenotypes = function(){ 
    
    Overall = tibble()
    Overall_age = tibble() #Biological age comparison among them
    Overall_age_noage = tibble() #Biological age comparison among them with no age correction
    DONE = F
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
      
      for (Col_m in 2:ncol(Data)){
        Explanatory  = as.vector(as_vector(Data[,Col_m]))
        Name_2 = colnames(Data)[Col_m]
        Name_2= clean_name(Name_2)
        To_plot = c("fathersmk", "AGE_FATHER", "AGE_MOTHER")
        #To_plot = c("Anti-CCP (U/mL)", "P","PQ", "QRS", "QT","QTC","T")
        #To_plot = c("Anti-CCP_UmL")
        if (Name %in% To_plot){ PL = T
        } else{ PL = F}
        if (Categorical == T){
          Results = Logistic_regression(Dependent, Explanatory, Covariates, Name_2, Name, Plot= PL)
        }else{
          Results = Linear_regression(Dependent, Explanatory, Covariates,Name_2, Name, Plot = PL)
        }
        
        #Count how many entreis per level
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
        #If it is the first iteration, check also how similar the measurements are among them
        if (DONE == T){ next }
        #if ( (Col_m+1) > ncol(Data)){next}
        #for (Col_m2 in (Col_m+1):ncol(Data)){
        for (Col_m2 in 2:ncol(Data)){
          if(Col_m2 == Col_m){next}
          Dependent_n = as.numeric(as.vector(as_vector(Data[,Col_m2])))
          Name_n2 = colnames(Data)[Col_m2]
          Name_n2= clean_name(Name_n2)
          Results =  Linear_regression(scale(Dependent_n), scale(as.numeric(Explanatory)),Covariates, Name_2, Name_n2, Plot=F)
          Results_noage =  Linear_regression(scale(Dependent_n), scale(as.numeric(Explanatory)),Covariates, Name_2, Name_n2, Plot=F, wAge ="No")
          
          LEVEL =  "Numeric" 
          N_LEVEL = length(Dependent[!is.na(Dependent_n)])
           
          Results %>% mutate(Levels = LEVEL, Levels_n = N_LEVEL) -> Results
          Results_noage %>% mutate(Levels = LEVEL, Levels_n = N_LEVEL) -> Results_noage
          
          Overall_age = rbind(Overall_age, Results)
          Overall_age_noage =  rbind(Overall_age_noage,Results_noage)
          
          
        }
      }
      DONE = T
    }
 
    return(list(Overall, Overall_age,Overall_age_noage))     
}

Model_phenotypes() -> list_results
list_results[[1]] %>% mutate(Beta= as.numeric(Beta),Pvalue = as.numeric(Pvalue),FDR = p.adjust(Pvalue,"fdr")) %>% arrange(Pvalue) -> Overall
list_results[[2]] %>% mutate(Beta= as.numeric(Beta),Pvalue = as.numeric(Pvalue), FDR = p.adjust(Pvalue,"fdr")) %>% arrange(Pvalue) -> Overall_age
list_results[[3]] %>% mutate(Beta= as.numeric(Beta),Pvalue = as.numeric(Pvalue), FDR = p.adjust(Pvalue,"fdr")) %>% arrange(Pvalue) -> Overall_age_notcorrected

Overall %>% filter(Dependent %in% c("Parent_smoker", "fathersmk", "mothersmk"))





Overall %>% mutate(Significance = ifelse(FDR < 0.05, "FDR" ,ifelse(Pvalue < 0.01, "Nominal_0.01", "Not_significant"))) -> Overall
Overall %>% ggplot(aes(x=Beta, y=-log10(Pvalue), col=Significance)) + geom_point() + theme_bw() + facet_wrap(~ Regressor, scales="free")
Overall %>% group_by(Regressor, Significance) %>% summarise(N=n()) %>% drop_na() %>% ggplot() + geom_bar(aes(x=Regressor, y=N, fill=Significance), stat="identity", position ="dodge") + theme_bw() + scale_y_continuous(trans='log10') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

Overall %>% mutate(Significant_FDR = ifelse(FDR< 0.05, 1, 0)) %>% select(c(Regressor, Dependent, Significant_FDR)) %>% spread(Dependent, Significant_FDR) -> matrix_presence
  t(matrix_presence) %>% as_tibble() %>% `colnames<-`(matrix_presence$Regressor)  -> matrix_presence
  matrix_presence %>% filter(!Age == "Age") -> matrix_presence
  apply(matrix_presence,2, FUN = as.numeric)  -> matrix_presence
  #apply(matrix_presence, 1, FUN= function(x){ sum(as.numeric(x))})
upset(as.data.frame(matrix_presence),sets= colnames(matrix_presence), text.scale = c(1, 1, 1, 1, 1, 2))



write_tsv(Overall,"C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Results_associations.tsv")
write_tsv(Overall_age,"C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Results_associations_comparisons.tsv")


Overall_age_notcorrected %>% select(Regressor, Dependent, Beta) %>% spread(Dependent, Beta) %>% filter(! is.na(Regressor)) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results2
Overall_age_notcorrected %>% select(Regressor, Dependent, Pvalue) %>% mutate(Pvalue = round(-log10(Pvalue), 2)) %>% spread(Dependent, Pvalue) %>% filter(! is.na(Regressor)) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results

pheatmap::pheatmap(wide_results2, display_numbers = wide_results,cluster_rows=FALSE, cluster_cols=FALSE)

Overall_age %>% select(Regressor, Dependent, Beta) %>% spread(Dependent, Beta) %>% filter(! is.na(Regressor)) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results2
Overall_age %>% select(Regressor, Dependent, Pvalue) %>% mutate(Pvalue = round(-log10(Pvalue), 2)) %>% spread(Dependent, Pvalue) %>% filter(! is.na(Regressor)) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results

wide_results2[wide_results < -log10(0.05)] = 0
wide_results[wide_results < -log10(0.05)] = ""

pheatmap::pheatmap(wide_results2, display_numbers = wide_results,cluster_rows=FALSE, cluster_cols=FALSE)

##########Proteomics

Model_proteomics = function(Proteomics2 = Proteomics){
  Proteomics2$Protein -> Proteins 
  t(select(Proteomics2, -Protein)) %>% as_tibble() %>% mutate(ID = colnames(select(Proteomics2, -Protein))) -> Proteomics2
  colnames(Proteomics2) = c(Proteins,"ID")
  
  Age_pheno = read_tsv(file = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Age.txt")
  Age_pheno %>% arrange(LLDEEPID) -> Age_pheno
  Data = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Combined_data_uncorrected.tsv")
  Data %>% filter(Sample %in% Proteomics2$ID) %>% arrange(Sample) -> D
  Phenotypes = read_csv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Merged_Phenotypes_LLD.csv")
  
  Proteomics2 %>% filter(ID %in% D$Sample) %>% arrange(ID) %>% select(-ID) -> Proteomics2
  Age_pheno %>% filter(LLDEEPID %in% D$Sample) -> Age_pheno
  Phenotypes %>% filter(ID %in% D$Sample) -> Phenos
  Covariates = tibble(Age = Age_pheno$AgeIncludingMonth, Sex = Phenos$Sex, ID = Age_pheno$LLDEEPID)
  
  Overall = tibble()
  Overall2 = tibble()
  for (Col_n in 1:ncol(Proteomics2)){
    Name = colnames(Proteomics2)[Col_n]
    print(Name)
    Dependent = as.vector(as_vector((Proteomics2[Col_n])))
    for (Col_m in 2:ncol(D)){
      Explanatory  = as.vector(as_vector(D[,Col_m]))
      Name_2 = colnames(D)[Col_m]
      Name_2= clean_name(Name_2)
      Results = Linear_regression(Dependent, Explanatory, Covariates,Name_2, Name, Plot = F)
      Results2 = Linear_regression(Dependent, Explanatory, Covariates,Name_2, Name, Plot = F, wAge="No")
      Overall = rbind(Overall, Results)
      Overall2 = rbind(Overall2, Results2)
    }
  }
  Overall %>% mutate(Pvalue = as.numeric(Pvalue),Beta= as.numeric(Beta), FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Results_proteomics
  Overall2 %>% mutate(Pvalue = as.numeric(Pvalue),Beta= as.numeric(Beta), FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Results_proteomics2
  return(list(Results_proteomics, Result_proteomics2))  
}

Model_proteomics() -> Results_proteomics
Results_proteomics_age = Results_proteomics[[2]]
Results_proteomics = Results_proteomics[[1]]
write_tsv(x = Results_proteomics, path="C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Proteomics_associations.tsv")
Results_proteomics %>% mutate(Pvalue = as.numeric(Pvalue), Beta = as.numeric(Beta)) %>% filter(! Regressor == "Age" &  ! grepl("INT", Regressor)) -> Results_proteomics

Results_proteomics %>% mutate(log10Pval = round(-log10(Results_proteomics$Pvalue),2)) %>% dplyr::select(c(Regressor, Dependent, log10Pval)) %>% spread(Regressor, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results
Results_proteomics %>% select(Dependent, Regressor, Beta) %>% spread(Regressor, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results2

#wide_results_Annotation = wide_results2 ; wide_results_Annotation[wide_results_Annotation < "0"] <- "-" ;wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
wide_results_Annotation = wide_results
wide_results_Annotation[wide_results < -log10(0.05)] <- "" ; wide_results_Annotation[! wide_results_Annotation == ""] = "*"
wide_results2[wide_results < -log10(0.05)] <- 0

pheatmap::pheatmap(wide_results2, display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 9, fontsize_row = 5 ) #-> Fig_0


Results_proteomics %>% mutate(Significance = ifelse(FDR < 0.05, "FDR" ,ifelse(Pvalue < 0.01, "Nominal_0.01", "Not_significant"))) -> Results_proteomics
Results_proteomics %>% group_by(Regressor, Significance) %>% summarise(N=n()) %>% drop_na() %>% ggplot() + geom_bar(aes(x=Regressor, y=N, fill=Significance), stat="identity", position ="dodge") + theme_bw() + scale_y_continuous(trans='log10') + theme(axis.text.x = element_text(angle = 90, hjust = 1))

Results_proteomics %>% mutate(Significant_FDR = ifelse(FDR< 0.05, 1, 0)) %>% select(c(Regressor, Dependent, Significant_FDR)) %>% spread(Dependent, Significant_FDR) -> matrix_presence
t(matrix_presence) %>% as_tibble() %>% `colnames<-`(matrix_presence$Regressor) %>% filter(!dCT == "dCT")  -> matrix_presence
apply(matrix_presence,2, FUN = as.numeric)  -> matrix_presence
#apply(matrix_presence, 1, FUN= function(x){ sum(as.numeric(x))})
upset(as.data.frame(matrix_presence),sets= colnames(matrix_presence), text.scale = c(1, 1, 1, 1, 1, 2))




#Proteomics set that represents aging:
Proteomics = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt")
Aging_proteins = colnames(read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/UMCG_senescence_set.txt"))
change_names = function(x){
  if(grepl("IL", x)){ paste(str_split(x,"L")[[1]],collapse="L-")
  }else if(grepl("MMP", x)){ paste(str_split(x,"MP")[[1]],collapse="MP-")
    
  } else{ return(x)}
}

sapply(Aging_proteins, FUN = change_names) -> Aging_proteins2
names(Aging_proteins2) = NULL

Aging_proteins2 <- read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/gene2protein.txt")
Proteomics %>% filter(Protein %in% Aging_proteins2$protein) -> Proteomics_age

Results_proteomics = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Proteomics_associations.tsv")
Results_proteomics %>% filter(Dependent %in% Proteomics_age$Protein) %>% filter(!grepl("INT",Regressor)) -> Results_proteomics_aging_markers 
Results_proteomics_Age %>% filter(Dependent %in% Proteomics_age$Protein) %>% filter(!grepl("INT",Regressor)) -> Results_proteomics2_aging_markers 


Results_proteomics_aging_markers %>% mutate(log10Pval = -log10(Results_proteomics_aging_markers$Pvalue)) %>% dplyr::select(c(Regressor, Dependent, log10Pval)) %>% spread(Regressor, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results
Results_proteomics_aging_markers %>% select(Dependent, Regressor, Beta) %>% spread(Regressor, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results2
wide_results_Annotation = wide_results2 ; wide_results_Annotation[wide_results_Annotation < "0"] <- "-" ;wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
wide_results_Annotation[wide_results < -log10(0.05)] <- ""
pheatmap::pheatmap(wide_results, display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 9, fontsize_row = 5 ) #-> Fig_0

Results_proteomics2_aging_markers %>% mutate(log10Pval = -log10(Results_proteomics_aging_markers$Pvalue)) %>% dplyr::select(c(Regressor, Dependent, log10Pval)) %>% spread(Regressor, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results
Results_proteomics2_aging_markers %>% select(Dependent, Regressor, Beta) %>% spread(Regressor, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results2
wide_results_Annotation = wide_results2 ; wide_results_Annotation[wide_results_Annotation < "0"] <- "-" ;wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
wide_results_Annotation[wide_results < -log10(0.05)] <- ""
pheatmap::pheatmap(wide_results, display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 9, fontsize_row = 5 ) #-> Fig_0






####Part 2 correlate methylation in CpG islands of telomerase with Telomere-related phenotypes

#1.Correlation of CpG islands in TERT and smoking/parental age
Methylation = function(){

    Phenotypes = read_csv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Merged_Phenotypes_LLD.csv")
    Phenotypes %>% select(c(ID,AGE_FATHER, fathersmk, `Has your father ever smoked regularly during your childhood?`, AGE_MOTHER)) -> Smoke_and_age
    TERT_methylation = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/TERT_CpG.tsv")
    TERT_methylation %>% summarise_if(is_numeric, sum)  -> Total
    TERT_methylation %>% summarise_if(is_numeric, mean)  -> Total2
    
    TERC_methylation = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/TERC_CpG.tsv")
    TERT_methylation %>% summarise_if(is_numeric, mean)  -> Total2_2
    
    rbind(TERT_methylation, select(mutate(Total, `-` = "Total"), c(`-`, colnames(Total)))) -> TERT_methylation
    rbind(TERT_methylation, select(mutate(Total2, `-` = "Total2"), c(`-`, colnames(Total2)))) -> TERT_methylation
    rbind(TERC_methylation, select(mutate(Total2, `-` = "Total"), c(`-`, colnames(Total2_2)))) -> TERC_methylation
    
    
    
    
    Smoke_and_age %>% filter(ID %in% colnames(TERT_methylation)) %>% arrange(ID) -> Smoke_and_age
    TERT_methylation %>% select(c("-", Smoke_and_age$ID)) -> TERT_methylation
    t(TERT_methylation)[2:dim(TERT_methylation)[2],] %>% as_tibble() %>% mutate(ID= colnames(TERT_methylation)[2:length(colnames(TERT_methylation))]) %>%
      `colnames<-`(c(TERT_methylation$`-`, "ID")) %>% arrange(ID) -> Methylation_info
    
    t(TERC_methylation)[2:dim(TERC_methylation)[2],] %>% as_tibble() %>% mutate(ID= colnames(TERC_methylation)[2:length(colnames(TERC_methylation))]) %>%
      `colnames<-`(c(TERC_methylation$`-`, "ID")) %>% arrange(ID) -> Methylation_info2
    
    
    Age_pheno = read_tsv(file = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Age.txt")
    Age_pheno %>% filter(LLDEEPID %in% Methylation_info$ID) %>% arrange(LLDEEPID) %>% select(AgeIncludingMonth) -> Cov_age
    
    Results_met_total = tibble()
    Results_met_TERC = tibble()
    
    Results_met_total_effect = tibble()
    Results_met_TERC_effect = tibble()
    
    for (Dependent in colnames(Methylation_info)){
        if (Dependent == "ID"){ next }
        Dependent = select(Methylation_info, Dependent)
        cbind(Smoke_and_age, dependent = as.numeric(as.vector(as_vector(Dependent)))) %>% as_tibble() -> Model_info
        Model_info %>% filter(ID  %in% Age_pheno$LLDEEPID) %>% mutate(Age = Cov_age$AgeIncludingMonth) -> Model_info
    
        #ggplot(Model_info) + geom_density(aes(x=dependent))
        #ggplot(Model_info) + geom_point(aes(y= dependent, x=fathersmk))
        summary(lm(dependent ~AGE_FATHER + Age, Model_info)) -> Model_ageF
        summary(lm(dependent ~AGE_MOTHER + Age, Model_info)) -> Model_ageM
        as.vector(lm(dependent ~ Age, Model_info)$residuals) -> Residuals
        Model_info %>% mutate(dependent = Residuals) -> Model_info 
        if (colnames(Dependent) %in% c("Total2")){
          Model_info %>% drop_na() %>% ggplot(aes(y= dependent, x=as.factor(fathersmk))) + geom_violin()  + geom_boxplot() + theme_bw() -> PLOT
          print(PLOT)
          Model_info %>% drop_na() %>% ggplot(aes(y= dependent, x=AGE_FATHER)) + geom_point() + theme_bw() + geom_smooth(method="lm" ) -> PLOT
          print(PLOT)
        }
        #wilcox.test(filter(Model_info, fathersmk==0)$dependent, filter(Model_info, fathersmk==1)$dependent) -> Model_smokeF
        summary(lm(dependent ~ fathersmk, data=Model_info)) -> Model_smokeF
        tibble(CG= colnames(Dependent), fathersmk= Model_smokeF$coefficients[8], Age_Father = Model_ageF$coefficients[11], Age_Mother = Model_ageM$coefficients[11]) -> Results_met
        tibble(CG= colnames(Dependent), fathersmk= Model_smokeF$coefficients[2], Age_Father = Model_ageF$coefficients[2], Age_Mother = Model_ageM$coefficients[2]) -> Results_met_effect
        
        Results_met_total = rbind(Results_met_total, Results_met)
        Results_met_total_effect = rbind(Results_met_total_effect, Results_met_effect)
    }
    return(list(Results_met_total,Results_met_total_effect))
}
Methylation() -> list_methylation
list_methylation[[1]] -> Results_met_total ; list_methylation[[2]] -> Results_met_effects
Results_met_total %>% gather(Variable, Pvalue, 2:4, factor_key=TRUE)-> Results_met_total_2 

write_tsv(path = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/TERT_CpG.tsv", x=Results_met_effects)

as.data.frame(Results_met_total) %>% column_to_rownames(var = "CG") -> wide_results
as.data.frame(Results_met_effects) %>% column_to_rownames(var = "CG") -> wide_results_Annotation
wide_results_Annotation[wide_results_Annotation < "0"] <- "-"
wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
pheatmap::pheatmap(-log10(wide_results), display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 10) #-> Fig_0





#2. Correlation of CpG islands in TERC and smoking/parental age
Methylation_TERC = function(TERC_methylation = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/TERC_CpG.tsv")){
  
  Phenotypes = read_csv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Merged_Phenotypes_LLD.csv")
  Phenotypes %>% select(c(ID,AGE_FATHER, fathersmk, `Has your father ever smoked regularly during your childhood?`, AGE_MOTHER)) -> Smoke_and_age
  
  TERC_methylation %>% summarise_if(is_numeric, mean)  -> Total2_2
  
  rbind(TERC_methylation, select(mutate(Total2_2, `-` = "Total"), c(`-`, colnames(Total2_2)))) -> TERC_methylation
  
  Smoke_and_age %>% filter(ID %in% colnames(TERC_methylation)) %>% arrange(ID) -> Smoke_and_age
  
  TERC_methylation %>% select(c("-", Smoke_and_age$ID)) -> TERC_methylation
  t(TERC_methylation)[2:dim(TERC_methylation)[2],] %>% as_tibble() %>% mutate(ID= colnames(TERC_methylation)[2:length(colnames(TERC_methylation))]) %>%
    `colnames<-`(c(TERC_methylation$`-`, "ID")) %>% arrange(ID) -> Methylation_info
  
  Age_pheno = read_tsv(file = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Age.txt")
  Age_pheno %>% filter(LLDEEPID %in% Methylation_info$ID) %>% arrange(LLDEEPID) %>% select(AgeIncludingMonth) -> Cov_age
  
  Results_met_total = tibble()
  Results_met_total_effect = tibble()
  for (Dependent in colnames(Methylation_info)){
    if (Dependent == "ID"){ next }
    Dependent = select(Methylation_info, Dependent)
    cbind(Smoke_and_age, dependent = as.numeric(as.vector(as_vector(Dependent)))) %>% as_tibble() -> Model_info
    Model_info %>% filter(ID  %in% Age_pheno$LLDEEPID) %>% mutate(Age = Cov_age$AgeIncludingMonth) -> Model_info
    
    summary(lm(dependent ~AGE_FATHER + Age, Model_info)) -> Model_ageF
    summary(lm(dependent ~AGE_MOTHER + Age, Model_info)) -> Model_ageM
    as.vector(lm(dependent ~ Age, Model_info)$residuals) -> Residuals
    Model_info %>% mutate(dependent = Residuals) -> Model_info 
    if (colnames(Dependent) %in% c("cg25090302")){
      Model_info %>% drop_na() %>% ggplot(aes(y= dependent, x=as.factor(fathersmk))) + geom_violin()  + geom_boxplot() + theme_bw() -> PLOT
      print(PLOT)
      Model_info %>% drop_na() %>% ggplot(aes(y= dependent, x=AGE_FATHER)) + geom_point() + theme_bw() + geom_smooth(method="lm" ) -> PLOT
      print(PLOT)
    }
    
    summary(lm(dependent ~ fathersmk, data=Model_info)) -> Model_smokeF
    tibble(CG= colnames(Dependent), fathersmk= Model_smokeF$coefficients[8], Age_Father = Model_ageF$coefficients[11], Age_Mother = Model_ageM$coefficients[11]) -> Results_met
    tibble(CG= colnames(Dependent), fathersmk= Model_smokeF$coefficients[2], Age_Father = Model_ageF$coefficients[2], Age_Mother = Model_ageM$coefficients[2]) -> Results_met_effect
    
    Results_met_total = rbind(Results_met_total, Results_met)
    Results_met_total_effect = rbind(Results_met_total_effect, Results_met_effect)
    
    
    
    
  }
  return(list(Results_met_total, Results_met_total_effect))
}
Methylation_TERC() -> list_TERC
Results_met_TERC <- list_TERC[[1]] ; Results_met_TERC_effects <- list_TERC[[2]]
Results_met_TERC %>% gather(Variable, Pvalue, 2:4, factor_key=TRUE)-> Results_met_total_TERC
as.data.frame(Results_met_TERC) %>% column_to_rownames(var = "CG") -> wide_results
as.data.frame(Results_met_TERC_effects) %>% column_to_rownames(var = "CG") -> wide_results_Annotation
wide_results_Annotation[wide_results_Annotation < "0"] <- "-" ; wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
pheatmap::pheatmap(-log10(wide_results), display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 10) #-> Fig_0

write_tsv(path = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/TERC_CpG.tsv", x=Results_met_total_TERC)



#Correlation of the eQTMs
TERT_significant = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/TERT_significant_CpG.tsv")
TERT_significant %>% select(-"X822") -> TERT_significant
Methylation_TERC(TERT_significant)


Data = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Combined_data_uncorrected.tsv")
Data %>% arrange(Sample) %>% filter(Sample %in% colnames(TERT_significant)) -> Data

Age_pheno = read_tsv(file = "C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Age.txt")
Age_pheno %>% arrange(LLDEEPID) %>% filter(LLDEEPID %in% Data$Sample ) -> Age_pheno

TERT_significant %>% filter(`-` == "cg00968538") %>% select(Data$Sample) %>% t() %>% as_vector() -> TERT_significant_d

for (Telomere in c("MTL_lymp","MTL_gran","MTL_CD45+_CD20-","MTL_CD45-","MTL_CD20+","MTL_CD57+")){
  Data %>% select(Telomere) %>% as_vector %>% as.numeric() -> Value
  print(Telomere)
  print(summary(lm(Value ~ TERT_significant_d + Age_pheno$AgeIncludingMonth)))
  }
TERT_significant %>% filter(`-` == "cg09951201") %>% select(Data$Sample) %>% t() %>% as_vector() -> TERT_significant_d
for (Telomere in c("MTL_lymp","MTL_gran","MTL_CD45+_CD20-","MTL_CD45-","MTL_CD20+","MTL_CD57+")){
  Data %>% select(Telomere) %>% as_vector %>% as.numeric() -> Value
  print(Telomere)
  print(summary(lm(Value ~ TERT_significant_d + Age_pheno$AgeIncludingMonth)))
}




##Study Smoking phenotypes

Smoke_pheno = read_csv2(file = "C:/Users/Sergio/Resilio Sync/LLD phenotypes for Sergio/Smoking phentoypes_/Smoking_LLD.csv")
Smoke_pheno2 = read_tsv("C:/Users/Sergio/Resilio Sync/LLD phenotypes for Sergio/Smoking phentoypes_/Current_smoker.txt")
Smoke_pheno %>% filter(LLDEEPID %in% Smoke_pheno2$ID) %>% mutate(Smk_now = as.factor(arrange(filter(Smoke_pheno2, ID %in% Smoke_pheno$LLDEEPID), ID)$smk)) ->Smoke_pheno

Smoke_pheno %>% select(c("smk13", "smk14", "smk15", "smk1", "smk10", "smk11", "Smk_now", "LLDEEPID")) -> Smoke_pheno
Smoke_pheno %>% mutate(Parental_smoking = ifelse(smk13 == 1 | smk14== 1, 1, ifelse(is.na(smk13) | is.na(smk14), NA, 0))) -> Smoke_pheno

Telomeres_data_n = colnames(Data)[grepl("MTL",colnames(Data))]
Data %>% select(c("Sample", Telomeres_data_n)) %>% filter(Sample %in% Smoke_pheno$LLDEEPID) %>% arrange(Sample) -> Telomeres_data
Covariates %>% filter(ID %in% Telomeres_data$Sample) %>% arrange(ID) %>% select(-ID) -> Covariates
Smoke_pheno %>% filter(LLDEEPID %in% Telomeres_data$Sample ) %>% arrange(LLDEEPID) %>% select(-LLDEEPID) -> Smoke_pheno

Parent_smoke = vector()
for (ID in 1:length(Smoke_pheno$smk13)){
  if( is.na(Smoke_pheno$smk13[ID]) | is.na(Smoke_pheno$smk14[ID]) ){
    Parent_smoke = c(Parent_smoke, NA)
  } else if (Smoke_pheno$smk13[ID] == 1 | Smoke_pheno$smk14[ID] == 1){
    Parent_smoke = c(Parent_smoke, 1)
  } else{
    Parent_smoke = c(Parent_smoke, 0)
  }
}
#Smoke_pheno %>% mutate(Parent_smoker = Parent_smoke) -> Smoke_pheno
Overall_smk = tibble()
for (Col_n in 2:ncol(Telomeres_data)){
  Name = colnames(Telomeres_data)[Col_n]
  print(Name)
  Name= clean_name(Name)
  Dependent = as.vector(as_vector((Telomeres_data[Col_n])))
  Smoke_pheno %>% mutate(Dependent = as.numeric(Dependent)) -> Input_cov
  Input_cov %>% drop_na() -> Input_cov
  summary(lm(Dependent ~ ., Input_cov)) -> Complete_model
  print(Complete_model)
  
  
  for (Col_m in 1:ncol(Smoke_pheno)){
    Explanatory  = as.numeric(as_vector(Smoke_pheno[,Col_m]))
    Name_2 = colnames(Smoke_pheno)[Col_m]
    
    if (Name_2 %in% c("smk15") & Name == "MTL_lymp"){ P = T
    } else { P = F}
    #if (length(unique(Explanatory)) < 4){
    #Results = Logistic_regression(Explanatory, Dependent, Covariates, Name_2, Name, Plot= F)
    #} else{
      Results = Linear_regression(Explanatory, Dependent, Covariates, Name_2, Name, Plot= P)  
    #}
    LEVEL =  paste(unique(Explanatory[!is.na(Explanatory)]),collapse=",")
    N_LEVEL = c()
    for (L in str_split(LEVEL, ",")[[1]]){
      N_LEVEL = c(N_LEVEL,length(which(as.character(Explanatory) == L)))
    }
      N_LEVEL = paste(N_LEVEL, collapse=",")
      
      Results %>% mutate(Levels = LEVEL, Levels_n = N_LEVEL) -> Results
      Overall_smk = rbind(Overall_smk, Results)
      
  }    

}

Overall_smk %>% mutate(Pvalue = as.numeric(Pvalue), Beta = as.numeric(Beta)) -> Overal_smk

write_tsv(Overal_smk,"C:/Users/Sergio/Resilio Sync/LLD phenotypes for Sergio/Results/Smoking_vs_Telomeres.csv")
Overall_smk %>% mutate(log10Pval = round(-log10(Overal_smk$Pvalue),2)) %>% dplyr::select(c(Regressor, Dependent, log10Pval)) %>% spread(Regressor, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results
Overall_smk %>% select(Dependent, Regressor, Beta) %>% spread(Regressor, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results2

#wide_results_Annotation = wide_results2
#wide_results_Annotation[wide_results_Annotation < "0"] <- "-"
#wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
wide_results_Annotation = wide_results
wide_results = apply(wide_results2,2, FUN = as.numeric)
rownames(wide_results) = rownames(wide_results2)
wide_results_Annotation[wide_results_Annotation < -log10(0.05)] = ""
wide_results[wide_results_Annotation < -log10(0.05)] = 0


breaksList = seq(7, 30, by = 1)
library(RColorBrewer)
pheatmap::pheatmap(wide_results, display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 10,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(length(breaksList))) #-> Fig_0



Partial_correlation = function(Telomeres, Covariates, Phenos){
  #Partial correlaton AGE father vs Smoking in telomere length
  Telomeres = colnames(Data)[grepl("MTL_", colnames(Data))]
  Data %>% select(Telomeres) -> Data_telomeres
  Phenotypes %>% select(c(AGE_FATHER, fathersmk, `Has your father ever smoked regularly during your childhood?`, AGE_MOTHER)) -> Smoke_and_age
  
  for (Col_m in 1:ncol(Data_telomeres)){
    Explanatory  = as.numeric(as.vector(as_vector(Data_telomeres[,Col_m])))
    Name_2 = colnames(Data_telomeres)[Col_m]
    Name_2= clean_name(Name_2)
    
    Partial_correlation_input = tibble(Telomere= Explanatory, Age_father= Smoke_and_age$AGE_FATHER, Smoking = Smoke_and_age$fathersmk,Sex= Covariates$Sex, Age=Covariates$Age)
    Partial_correlation_input %>% filter(! is.na(Telomere)) %>% filter( ! is.na(Age)) -> Partial_correlation_input
    Partial_correlation_input %>% ggplot(aes(y=Telomere,x=as.factor(Smoking))) + geom_boxplot()
    
    
    Explanatory2 = summary(lm(Telomere ~ Age, Partial_correlation_input))$residuals
    
    Partial_correlation_input %>% mutate(Telomere = Explanatory2) -> Partial_correlation_input
    
    Model1 = lm(Telomere ~ Age_father  + Sex , Partial_correlation_input)
    Model2 = lm(Telomere ~ as_factor(Smoking)  + Sex , Partial_correlation_input)
    Model3 = lm(Telomere ~ Age_father + as_factor(Smoking) + Sex , Partial_correlation_input)
    
    Partial_correlation_input %>% ggplot(aes(y=Telomere,x=Age_father)) + geom_point() + geom_smooth()
    Partial_correlation_input %>% ggplot(aes(y=Telomere,x=as.factor(Smoking))) + geom_boxplot()
  }  
  
  
}



#Study associations to Telomeres including Genotypes
Covariates = tibble(Age = Age_pheno$AgeIncludingMonth, Sex = filter(Phenotypes, ID %in% Age_pheno$LLDEEPID)$Sex, ID = Age_pheno$LLDEEPID)

Overall = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Results_associations.tsv")
Overall %>% filter(grepl("MTL",Regressor)) %>% filter(FDR< 0.05) -> FDR_Telo
Phenotypes = read_csv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Merged_Phenotypes_LLD.csv") ; Rm = colnames(Phenotypes)[grepl("Record", colnames(Phenotypes))]
colnames(Phenotypes) = clean_name(colnames(Phenotypes))
Parent_smoke = vector()
for (ID in 1:length(Phenotypes$fathersmk)){
  if( is.na(Phenotypes$fathersmk[ID]) | is.na(Phenotypes$mothersmk[ID]) ){
    Parent_smoke = c(Parent_smoke, NA)
  } else if (Phenotypes$fathersmk[ID] == 1 | Phenotypes$mothersmk[ID] == 1){
    Parent_smoke = c(Parent_smoke, 1)
  } else{
    Parent_smoke = c(Parent_smoke, 0)
  }
}
Phenotypes %>% mutate(Parent_smoker = Parent_smoke) -> Phenotypes
Phenotypes %>% arrange(ID) %>% select(c("ID",unique(FDR_Telo$Dependent))) -> Phenotypes
Data = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Combined_data_uncorrected.tsv") ; Data %>% arrange(Sample) %>% filter(Sample %in% Age_pheno$LLDEEPID) -> Data
Data %>% select(c("Sample", colnames(Data)[grepl("MTL", colnames(Data))])) -> Data_telomeres
Data_telomeres %>% filter(Sample %in% Phenotypes$ID) %>% arrange(Sample) -> Data_telomeres
Genotypes = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Genotypes_TERT_LLD.tsv")
Genotypes %>% filter(ID %in% Data_telomeres$Sample) %>% arrange(ID) -> Genotypes

Data_telomeres #Dependent
Phenotypes #Regressor
Genotypes #Covariate
Covariates  #Sex and Age

Data_telomeres %>% filter(Sample %in% Genotypes$ID) -> Data_telomeres
Phenotypes %>% filter(ID %in% Genotypes$ID) -> Phenotypes
Covariates %>% filter(ID %in% Genotypes$ID) -> Covariates

Covariates %>% mutate(Genotype = Genotypes$GG*2 + 0*Genotypes$AA + 1*Genotypes$GA) -> Covariates

Results_with_Geno = tibble()
for (Col_n in 2:ncol(Data_telomeres)){
  Name = colnames(Data_telomeres)[Col_n]
  print(Name)
  Name= clean_name(Name)
  Dependent = as.vector(as_vector((Data_telomeres[Col_n])))
  
  for (Col_m in 2:ncol(Phenotypes)){
    Explanatory  = as.numeric(as.vector(as_vector(Phenotypes[,Col_m])))
    Name_2 = colnames(Phenotypes)[Col_m]
    Input_model = tibble(Dependent = as.numeric(Dependent),Explanatory = Explanatory, Age= Covariates$Age, Sex=Covariates$Sex, Genotype=Covariates$Genotype)
    Model = lm(Dependent ~Explanatory + Age + Sex + Genotype, Input_model)
    Model0 = lm(Dependent ~Explanatory + Age + Sex, Input_model)
    Results = tibble(Pvalue = summary(Model)$coefficients[17], Beta = summary(Model)$coefficients[2], Dependent = Name, Explanatory = Name_2, Model ="geno")
    Results0 = tibble(Pvalue = summary(Model0)$coefficients[14], Beta = summary(Model0)$coefficients[2], Dependent = Name, Explanatory = Name_2, Model ="regular")
    Results = left_join(Results, Results0, by=c("Dependent", "Explanatory"), suffix= c("Geno", "NoGeno"))
    Results_with_Geno = rbind(Results_with_Geno, Results)
  
  }
}

Results_with_Geno %>% ggplot(aes(x=-log10(PvalueGeno), y=-log10(PvalueNoGeno))) + geom_point() +theme_bw() + geom_abline()
Results_with_Geno %>% ggplot(aes(x=-log10(PvalueGeno), y=-log10(PvalueNoGeno))) + geom_point() +theme_bw() + geom_abline()

write_tsv(Results_with_Geno,path ="C:/Users/Sergio/Resilio Sync/LLD phenotypes for Sergio/Results/Telomere_associations_GenotypeCov.csv")

Results_with_Geno %>% gather(Variable, Pvalue, 2:4, factor_key=TRUE)-> wide_results


Results_with_Geno %>% mutate(log10Pval = -log10(Pvalue)) %>% dplyr::select(c(Explanatory, Dependent, log10Pval)) %>% spread(Explanatory, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results
Results_with_Geno %>% select(Dependent, Explanatory, Beta) %>% spread(Explanatory, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results2
wide_results_Annotation = wide_results2 ; wide_results_Annotation[wide_results_Annotation < "0"] <- "-" ;wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
pheatmap::pheatmap(wide_results, display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 9, fontsize_row = 5 ) #-> Fig_0


#Predictions


##Functions
Train_and_test = function(Training_set, Test_set, Phenotype, Predictor,Sex){
  Phenotypes_followup %>% arrange(ID) %>% select(c(ID, Time,Phenotype)) -> Phenos_interest
  
  Phenos_interest %>% filter(ID %in% Train_ID) -> Train_ys
  Phenos_interest %>% filter(ID %in% Test_set$Sample) -> Test_ys
  #Train
  
  filter(Train_ys, ! Time ==  "Baseline assessment (1A)") -> Y    #Use first time point as part of the model
  filter(Train_ys,  Time ==  "Baseline assessment (1A)") -> Cov
  
  Cov %>% filter(ID %in% Y$ID) -> Cov ; Y %>% filter(ID %in% Cov$ID) -> Y #Match time points (some samples lack one)
  Sex %>% arrange(ID) %>% filter(ID %in% Cov$ID) -> Cov_sex ; Cov_sex = as.numeric(Cov_sex$Sex)
  Cov %>% select(Phenotype) %>% as_vector() -> Cov  #Select the Phenotype to predict
  Training_set %>% filter(Sample %in% Y$ID) -> X
  Y %>% select(Phenotype) %>% as_vector() -> Y
  X %>% select(Predictor) %>% as_vector() %>% as.numeric() -> X #Select predictor
  
  lm(Y ~ X + Cov + Cov_sex) -> Model_training
  
  #Test
  filter(Test_ys, ! Time ==  "Baseline assessment (1A)") -> Y
  filter(Test_ys,  Time ==  "Baseline assessment (1A)") -> Cov
  Cov %>% filter(ID %in% Y$ID) -> Cov ; Y %>% filter(ID %in% Cov$ID) -> Y ; Sex %>% arrange(ID) %>% filter(ID %in% Cov$ID) -> Cov_sex ; Cov_sex = as.numeric(Cov_sex$Sex)
  Cov %>% select(Phenotype) %>% as_vector() -> Cov ; Test_set %>% filter(Sample %in% Y$ID) -> X ; Y %>% select(Phenotype) %>% as_vector() -> Y ; X %>% select(Predictor) %>% as_vector() %>% as.numeric() -> X
  
  
  New_data = tibble(Cov = Cov, X=X, Cov_sex=Cov_sex)
  predict(Model_training, newdata = New_data ) -> y_hat
  Data_prediction = tibble(Y_hat = y_hat, Y = Y) %>% drop_na()
  Error = RMSE(Data_prediction$Y_hat, Data_prediction$Y)
  R2 = calculate_R2(Data_prediction$Y, Data_prediction$Y_hat)
  ggplot(Data_prediction, aes(x=Y_hat, y=Y)) + geom_point() + theme_bw() -> Plot
  return(list(Error, R2))
}
RMSE = function(Y_hat, Y){ sqrt(sum((Y_hat - Y)^2)/length(Y)) }
error = function(Y_hat, Y){ 1 - (sum(Y_hat == Y)/length(Y))}
calculate_R2 = function(real, pred){
  rss = sum((pred - real)^2)
  tss = sum((real - mean(real))^2)
  rsq = 1 - rss/tss
  rsq
}
Fit_model_ANOVA = function(Set, Phenotype, Predictor, Cov){
  Phenotypes_followup %>% arrange(ID) %>% select(c(ID, Time,Phenotype)) -> Phenos_interest
  Cov %>% arrange(ID) -> Cov
  
  Phenos_interest %>% filter(ID %in% Cov$ID) %>% filter(ID %in% Set$Sample) -> Phenos_interest
  filter(Phenos_interest, ! Time ==  "Baseline assessment (1A)") -> Y    #Use first time point as part of the model
  filter(Phenos_interest,  Time ==  "Baseline assessment (1A)") -> Baseline0
  Y %>% filter(ID %in% Baseline0$ID) %>% arrange(ID) -> Y ; Baseline0 %>% filter(ID %in% Y$ID) %>% arrange(ID) -> Baseline0
  Cov %>% filter(ID %in% Baseline0$ID) -> Cov
  
  Baseline0 %>% select(Phenotype) %>% as_vector() -> Baseline  #Select the Phenotype to predict
  Y %>% select(Phenotype) %>% as_vector() -> Y
  Set  %>% filter(Sample %in% Baseline0$ID) %>% select(Predictor) %>% as_vector() %>% as.numeric() -> X #Select predictor
  
  Input = tibble(Y = Y, Age= Cov$Age, Sex=Cov$Sex, Baseline = Baseline, X = X) %>% drop_na()
  lm(Y ~ Age + Sex + Baseline + X, Input) -> Model_B
  lm(Y ~ Age + Sex + Baseline, Input) -> Model_A
  P = anova(Model_A,Model_B,test="Chisq")$`Pr(>Chi)`
  return(P[2])
}





Cov = read_csv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Merged_Phenotypes_LLD.csv") %>% select(ID, Age, Sex)
Pheno_P = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Phenotype_Others.tsv")
Pheno_H = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Phenotype_Howareyou.tsv")
left_join(Pheno_P, Pheno_H, c("ID", "Time")) -> Phenotypes_followup

Read_data() -> data_entries
Predictors <- data_entries[[2]]

#Divide in Train/Test IDs
Predictors %>% filter(Sample %in% Phenotypes_followup$ID) -> Predictors
N = length(Predictors$Sample)
set.seed(800)
Train_n = 0.8 *N ; Train_ID = sample(Predictors$Sample, size = Train_n, replace = F)
Test_n = 0.2*N

Predictors %>% filter(Sample %in% Train_ID) -> Training_set
Predictors %>% filter(! Sample %in% Train_ID) -> Test_set

Overall_result = tibble()
Best = tibble()
for (Pheno_to_predict in colnames(Phenotypes_followup) ){
  Best_prediction = ""
  Best_error = 100
  for (Predictor_pheno in colnames(select(Predictors, -Sample))){
    if (Pheno_to_predict %in% c("ID", "Time")){ next }
    Errors = Train_and_test(Training_set, Test_set, Pheno_to_predict, Predictor_pheno, Cov)
    Pvalue_comp = Fit_model_ANOVA(Predictors,Pheno_to_predict, Predictor_pheno, Cov)
    Error = Errors[[1]]
    R2 = Errors[[2]]
    Result = tibble(Phenotype = Pheno_to_predict, Predictor = Predictor_pheno, Error = Error, R2=R2, Pval=Pvalue_comp )
    Overall_result =rbind(Overall_result, Result)
    if (Best_error > Error){ Best_error = Error ; Best_prediction = Predictor_pheno  }
  }
  Best = rbind(Best, tibble(Phenotype = Pheno_to_predict, Best_predictor=Best_prediction ))
}
library(patchwork)
Total_plot = NULL
for (Phenotypes in unique(Overall_result$Phenotype)){
  Overall_result  %>% filter(Phenotype == Phenotypes) -> OR
  OR %>% mutate(Error = R2) -> OR
  Error_age = filter(OR, Predictor=="Age")$Error
  OR %>% mutate(Smaller_than_age = ifelse(Error < Error_age, "Smaller", ifelse(Error ==Error_age, "Equal", "Greater") )) -> OR
  ggplot() + geom_bar(aes_(x=OR$Predictor,y=OR$Error,fill=OR$Smaller_than_age), stat = "identity")  + theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(Phenotypes)  -> PLOT
  if (! is.null(Total_plot)){
    Total_plot = Total_plot +PLOT
  } else { Total_plot = PLOT }
}
Total_plot + plot_layout(guides="collect")
Overall_result %>% spread(Predictor, Error)

Best %>% group_by(Best_predictor) %>% summarise(n())




#Association in the second timepoint
Cov 
Phenotypes_followup
Predictors %>% mutate(dCT = as.numeric(dCT)) -> Predictors
Result_f = tibble() ; Result_f2 = tibble()

for (Predictor in colnames(select(Predictors, -Sample))){
  Predictors %>% select(c("Sample", Predictor)) %>% drop_na() %>% arrange(Sample) %>% filter(Sample %in% Phenotypes_followup$ID) %>% filter(Sample %in% Cov$ID) -> Pred_info
  Cov %>% filter(ID %in% Pred_info$Sample) %>% arrange(ID) -> Pred_cov
  Phenotypes_followup %>% filter(ID %in% Pred_info$Sample) %>% arrange(ID) -> Pred_phenos
  for (To_predict in colnames(select(Pred_phenos, -c("ID", "Time")))){
    Pred_phenos %>% select(c("ID", "Time", To_predict)) %>% spread(Time, To_predict) %>% drop_na() -> Pred_phenos2
    print(c(Predictor, To_predict))
    Pred_phenos2$`Baseline assessment (1A)` -> Covariate0
    Pred_phenos2$`Second assessment (2A)` -> Y
    
    Pred_info %>% filter(Sample %in% Pred_phenos2$ID) %>% select(Predictor) %>% as_vector() -> Regressor
    Pred_cov %>% filter(ID %in% Pred_phenos2$ID) %>% select(Sex) %>% as_vector() -> Covariate
    Pred_cov %>% filter(ID %in% Pred_phenos2$ID) %>% select(Age) %>% as_vector() -> Age
    
    
    tibble(Regressor = scale(as.numeric(Regressor)), Covariate0 = scale(as.numeric(Covariate0)), Covariate = as.numeric(Covariate), Y = scale(as.numeric(Y)), Age = as.numeric(Age)) %>% drop_na() -> Input
    summary(lm(Y ~ Covariate0 + Covariate + Regressor, data = Input))$coefficients %>% as.data.frame() -> Result
    
    
    Result["Regressor",] %>% as_tibble() %>% mutate(Regressor = Predictor, Dependent=To_predict ) -> Result
    rbind(Result_f, Result) -> Result_f
    
    if (Predictor == "Age"){ next }
    lm(Y ~ Covariate0 + Covariate + Age + Regressor, data = Input)-> Result2
    lm(Y ~ Covariate0 + Covariate + Age, data = Input) -> Result2_0
    anova(Result2, Result2_0,test="Chisq")$`Pr(>Chi)` -> Pvalue
    
    
    
    as.data.frame(summary(Result2)$coefficients)["Regressor",] %>% as_tibble() %>% mutate(Regressor = Predictor, Dependent=To_predict, P=Pvalue[2] ) -> Result2
    rbind(Result_f2, Result2) -> Result_f2
  }
  
    
}


breaksList = seq(7, 30, by = 1)


Result_f %>% mutate(log10Pval = round(-log10(`Pr(>|t|)`), 2)) %>% dplyr::select(c(Regressor, Dependent, log10Pval)) %>% spread(Regressor, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results
Result_f %>% select(Dependent, Regressor, Estimate) %>% spread(Regressor, Estimate) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results2
wide_results_Annotation = wide_results
wide_results_Annotation[wide_results_Annotation < -log10(0.05)] = ""
wide_results2[wide_results_Annotation < -log10(0.05)] = 0

pheatmap::pheatmap(wide_results2, display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 9, fontsize_row = 5,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))) #-> Fig_0


Result_f2 %>% mutate(log10Pval = round(-log10(`Pr(>|t|)`), 2)) %>% dplyr::select(c(Regressor, Dependent, log10Pval)) %>% spread(Regressor, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results
Result_f2 %>% select(Dependent, Regressor, Estimate) %>% spread(Regressor, Estimate) %>% as.data.frame() %>% column_to_rownames(var = "Dependent") -> wide_results2
wide_results_Annotation = wide_results
wide_results_Annotation[wide_results_Annotation < -log10(0.05)] = ""
wide_results2[wide_results_Annotation < -log10(0.05)] = 0
pheatmap::pheatmap(wide_results2, display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 9, fontsize_row = 5, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))) #-> Fig_0



#Prediction (association) of death 
dead = as_tibble(read.table("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/LLD_passed.txt",header = T))
Info = read_tsv("/Users/Sergio/Resilio Sync/LLD phenotypes for Sergio/Survival_/baseline.txt")

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
library("survminer")


Cov = read_csv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Merged_Phenotypes_LLD.csv") %>% select(ID, Age, Sex)
Predictors = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Combined_data_uncorrected.tsv")

Cov %>% filter(ID %in% Info$LLDEEPID) %>% arrange(ID) -> Cov
Predictors %>% filter(Sample %in% Info$LLDEEPID) %>% arrange(Sample) -> Predictors



Results = tibble()
for( Phenotype in colnames(select(Predictors,-c("Sample")))){
  print(Phenotype)
  Regressor = Predictors %>% select(Phenotype) %>% as_vector() %>% as.numeric()
  Info %>% filter(LLDEEPID %in% Predictors$Sample) -> Info2 ;  Cov %>% filter(ID  %in% Info2$LLDEEPID) -> Cov2
  Info2 %>% mutate(Age = Cov2$Age, Sex = Cov2$Sex, Pheno = Regressor) -> All_info
  coxph(Surv(Second_time, Status) ~ Age + Sex + Pheno, All_info) -> res.cox
  #coxph(Surv(Second_time, Status) ~ Sex + Pheno, All_info,) -> res.cox
  ggsurvplot(survfit(res.cox), color = "#2E9FDF",data=All_info,
             ggtheme = theme_minimal()) + ggtitle(Phenotype)-> FIG
  #print(FIG)
  
  as.data.frame(summary(res.cox)$coefficients)["Pheno",] %>% as_tibble() %>% mutate(Regressor = Phenotype) -> res.cox
  rbind(res.cox, Results) -> Results
  
}
arrange(Results, `Pr(>|z|)`)
ggplot(Results) + geom_bar(aes(x=Regressor, y= -log10(`Pr(>|z|)`) , fill = `exp(coef)` > 1), stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(yintercept = -log10(0.05))


#Correlation with cancer development
Read_data() -> List_inputs
Phenotypes <- List_inputs[[1]] ; Data <- List_inputs[[2]] ; Age_pheno <- List_inputs[[3]] ; Calories <- List_inputs[[4]] ; Proteomics <- List_inputs[[5]]
Predictors = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/Cancer_changes.tsv")
#Predictors = read_tsv("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Phenotype_Correlation/selected_phenos_ss_ts_ageing_prediction_final_with_changes")


Predictors %>% filter(LLDEEPID %in% Data$Sample) %>% arrange(LLDEEPID) -> Predictors
Phenotypes %>% select(ID, Sex, Age) %>% filter(ID %in% Predictors$LLDEEPID)  %>% arrange(ID)-> Cov
Data %>% filter(Sample %in% Predictors$LLDEEPID) %>% arrange(Sample) %>% select(-Sample)  -> Data2

fit_logit = function(x,y,Cov){
  tibble(y = y, x = scale(as.numeric(x)), Cov$Sex, Cov$Age) -> DATA
  #DATA %>% mutate(y = ifelse(is.na(y), "No", y))  -> DATA
  DATA %>% drop_na() %>% mutate(y = as.factor(y)) -> DATA
  glm(y ~ x , DATA, family = binomial) -> model

  as.data.frame(summary(model)$coefficients)["x",] -> model
  
  return(c(model$Estimate, model$`Pr(>|z|)`))
}


apply(Data2, 2, FUN = function(x){ fit_logit(x, Predictors$Cancer_change, Cov) } ) %>% as_tibble()  -> Cancer_results

Cancer_results %>% t() %>% as_tibble() %>% mutate(Biological_age = colnames(Cancer_results))  %>% ggplot() + geom_bar(aes(x=Biological_age, y=V2, fill= V1>0 ), stat="identity") + theme_bw() + coord_flip() + ylab("-log10(Pvalue)")
