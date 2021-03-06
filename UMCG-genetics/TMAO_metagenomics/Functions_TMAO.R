#############################
######General functions######
#############################
Rank_norm = function(Variable){
  qnorm((rank(Variable,na.last="keep")-0.5)/sum(!is.na(Variable)))
}



Association_metabolites = function(Normalized_dependent_metabolites, with_Covariates=F,Covariates="None"){
  Summary_results = tibble()
  for (Metabolite in colnames(Normalized_dependent_metabolites)){
    Dependent = Normalized_dependent_metabolites %>% select(Metabolite)
    Dependent = as_vector(Dependent)
    Regressors = Normalized_dependent_metabolites %>% select(-Metabolite)
    for (Metabolite2 in colnames(Regressors)){
      Regressor = Regressors %>% select(Metabolite2) %>% as_vector()
      Regressor = as_vector(Regressor)
      if (with_Covariates == F){
        Summary = as_tibble(summary(lm(Dependent ~ Regressor))$coefficients)
      } else if (with_Covariates == T){
        Summary = as_tibble(summary(lm(Dependent ~ Regressor + Covariates$Age + Covariates$Gender))$coefficients)
      } else if (with_Covariates == "BMI"){
        Summary = as_tibble(summary(lm(Dependent ~ Regressor + Covariates$Age + Covariates$Gender + Covariates$BMI))$coefficients)
        
      }
      
      Beta = Summary$Estimate[2]
      pvalue = Summary$`Pr(>|t|)`[2]
      stand = Summary$`Std. Error`[2]
      Results = as_tibble(t(c(Metabolite2, Metabolite, Beta, stand, pvalue)))
      colnames(Results) = c("Regressor", "Metabolite", "Beta", "SE", "Pvalue")
      Summary_results = rbind(Summary_results, Results) 
    }
  }
  return(Summary_results)
}
Metabolite_exploration = function(Phenos){
  # Ditributions with colours
  Phenos %>% melt() %>% ggplot() + geom_density(aes(x=value, fill=variable)) + facet_wrap(~variable, scales="free") + theme_bw() -> Distributions
  # Distribution of TMAO and focus on outlayers
  Phenos %>% melt() %>% filter(variable=="TMAO") -> TMAO_measures
  TMAO_measures %>% ggplot() + geom_histogram(aes(x=value)) + theme_bw()
  TMAO_measures %>% filter(value > 40) %>% arrange(value) -> TMAO_outlayers
  
  #Correlation table 
  Phenos %>% select(-ID) -> Check
  round(cor(Check,method = "spear"),2) %>% melt() -> correlation_phenos
  correlation_phenos %>% ggplot() + geom_tile(aes(x=Var1, y=Var2, fill=value)) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> correlation_matrix
  print(correlation_matrix)
  print(Distributions)
  return(list(TMAO_outlayers, correlation_phenos)) 
}
Linear_regression = function(Dependent, Regressor, Covariates, Regressor_name, Dependent_name, BMI=T){
  Dependent = Rank_norm(Dependent)
  Covariates %>% mutate(Regressor = Regressor, Dependent = Dependent) -> Model_variables
  if (BMI == T){
    Model = lm(Dependent ~ Regressor + Age + BMI + Gender + Creatinine, Model_variables)
  } else{ Model = lm(Dependent ~ Regressor + Age + Gender + Creatinine, Model_variables) }
  Summary = as_tibble(summary(Model)$coefficients)
  Beta = Summary$Estimate[2]
  pvalue = Summary$`Pr(>|t|)`[2]
  stand = Summary$`Std. Error`[2]
  N_samples = length(Model$fitted.values)
  Results = as_tibble(t(c(Regressor_name, Dependent_name, Beta, stand, pvalue, N_samples)))
  colnames(Results) = c("Regressor", "Metabolite", "Beta", "SE", "Pvalue", "N")
  return(Results)
}
Match_dataset = function(Metabolites, Covariates, Regressor){
  IDs = intersect(intersect(Metabolites$ID, Covariates$ID), Regressor$ID)
  Metabolites %>% filter(ID %in% IDs) %>% arrange(ID) -> Metab_filter
  Covariates  %>% filter(ID %in% IDs) %>% arrange(ID) -> Cov_filter
  Regressor   %>% filter(ID %in% IDs) %>% arrange(ID) -> Regr_filter
  return(list(Metab_filter, Cov_filter,Regr_filter))
}
Iteration_Regressors_and_model = function(Dependent, Regressors, Covariates, Dependent_name, BMI=T){
  Stats_per_Metabolite = tibble()
  for (Column_number in seq(dim(Regressors)[2])){
    Regressor_name = colnames(Regressors)[Column_number]
    Regressor = as.vector(as_vector(Regressors[,Column_number]))
    Summary_stats = Linear_regression(Dependent = Dependent, Regressor = Regressor, Covariates = Covariates, Regressor_name = Regressor_name, Dependent_name = Dependent_name, BMI=BMI)
    Stats_per_Metabolite = rbind(Stats_per_Metabolite, Summary_stats)
  }
  return(Stats_per_Metabolite)
}  
FDR_correction = function(Regressors, Pvalues, Dependent, Covariates, Dependent_name, Permutations = 100, BMI=T){
  Permuted_Pvalues = vector()
  for (Permutation_number in seq(Permutations)){
    Permuted = sample_n(Regressors, size= nrow(Regressors))
    Summary_stats = Iteration_Regressors_and_model(Dependent, Permuted, Covariates, Dependent_name, BMI)
    Permuted_Pvalues = c(Permuted_Pvalues,as.numeric(Summary_stats$Pvalue))
  }
  FDR_vector = vector()
  #Pvalues = sort(Pvalues) 
  for (Threshold in Pvalues){
    Threshold = as.numeric(Threshold)
    P = sum((Pvalues <= Threshold) *1)
    FP = sum((Permuted_Pvalues <= Threshold)*1)
    FDR = (FP/Permutations)/P
    if (FDR >1){ FDR = 1}
    #If there are already elements in the vector, dont let a higher pvalue have lower FDR than a lower pvalue. Make them equal.
    if (length(FDR_vector) > 0){ if (FDR > FDR_vector[length(FDR_vector)]){ FDR = FDR_vector[length(FDR_vector)]} }
    FDR_vector = c(FDR_vector,FDR)
  }
  return(FDR_vector)
}
Iterate_Metabolites = function(Metabolites, Covariates, Regressors,Analname,FDR_iter=100, BMI=T, Correct = "Permutations"){
  Total_results = tibble()
  Metabolites = select(Metabolites, -ID)
  Regressors = select(Regressors, -ID)
  for (Metabolite_column in seq(length(colnames(Metabolites)))){
    Dependent = as.vector(as_vector(Metabolites[,Metabolite_column]))
    Dependent_name = colnames(Metabolites)[Metabolite_column]
    #Make a model per Regressor in a given Metablite
    Summary_stats = Iteration_Regressors_and_model(Dependent, Regressors,Covariates,Dependent_name, BMI)
    #Calculate False discovery Rates
    Summary_stats %>% mutate(Pvalue = as.numeric(Pvalue)) %>% arrange(desc(Pvalue)) -> Summary_stats
    Summary_stats %>% filter(Pvalue < 0.05) -> DA
    Remove_null_distribution = DA$Regressor
    Summary_stats$Pvalue -> Pvalues
    if (Correct == "Permutations"){
      if (length(Remove_null_distribution) == dim(Regressors)[2] || FDR_iter == 0){ FDRs = sample(c(NA),size = dim(Regressors)[2],replace = T) 
      } else { FDRs = FDR_correction(select(Regressors, -Remove_null_distribution), Pvalues, Dependent, Covariates, Dependent_name, FDR_iter, BMI)
      }
    } else{ FDRs = p.adjust(Pvalues, "fdr") }  
    Summary_stats %>% mutate(FDR = FDRs) -> Summary_stats
    #Prepare Output
    Output = paste(c("Model_stats/",Dependent_name,"_",Analname,".tsv"), collapse="")
    #write_tsv(x= Summary_stats,path=Output)
    Total_results = rbind(Total_results, Summary_stats)
  }
  return(Total_results)
}

Transformation_composition = function(Count_table){
  #Clr transformation
  SampleID = Count_table$ID
  Count_table %>% select(-ID) -> Counts
  #Transform counts
  colnames(Counts) -> N
  Counts %>% t() -> Check
  Counts_transformed = as_tibble(t(abundances(x=as.data.frame(Check), transform="clr")))
  ###
  Counts_transformed %>% mutate(ID = SampleID) %>% select(ID,colnames(Counts_transformed)[1:(length(colnames(Counts_transformed)))]) -> Counts_transformed
  return(Counts_transformed)
}
Filter_name_meta = function(Abundance_names){
  Taxon = vector()
  list_split = strsplit(Abundance_names, "\\.")
  for (i in seq(1:length(list_split))){
    element =  list_split[[i]]
    Taxon = c(Taxon, element[length(element)])
  }
  return(Taxon)
}
Transformation_logarithmic = function(Features){
  for (Feature_n in seq(dim(Features)[2])){
    Feature = Features[,Feature_n]
    if (colnames(Feature) == "ID"){ next }
    Feature = as.vector(as_vector(Feature))
    Feature_na = Feature[! is.na(Feature)]
    pseudocount = min(Feature_na[Feature_na > 0])/3
    Transformed_feature = log10(Feature+pseudocount)
    Features[,Feature_n] = Transformed_feature
  }
  return(Features)
}
Make_genus = function(Abundances, Level=1){
  Genus_Vector = vector()
  
  for (Species_n in seq(dim(Abundances)[2])){
    Species_name = colnames(Abundances)[Species_n]
    if (Species_name == "ID"){next}
    vector_taxon = strsplit(Species_name, "\\.")
    Genus = vector_taxon[[1]][length(vector_taxon[[1]])-Level]
    Genus_Vector = c(Genus_Vector, Genus)
  }
  Total_genus = tibble()
  for (Indv_n in seq(dim(Abundances)[1])){
    Indv = Abundances[Indv_n,]
    tibble_Indv = tibble(ID=Indv$ID)
    Indv %>% gather(Species, Abundance, 2:dim(Abundances)[2], factor_key=TRUE) -> Indv
    for (Genus in unique(Genus_Vector)){
      Indv %>% group_by(grepl(pattern = Genus, Species)) %>% summarise(Value=sum(Abundance)) %>% filter(`grepl(pattern = Genus, Species)` == T) %>%
        select(Value) -> Value
      Entry = tibble(Value$Value)
      colnames(Entry) = Genus
      tibble_Indv = cbind(tibble_Indv,Entry)
    }
    Total_genus = rbind(Total_genus, as_tibble(tibble_Indv))
  }  
  return(Total_genus)
}
Gene_abundance_process = function(Gene_abundance, sum_up=F){
  Gene_abundance %>% filter(!grepl(x=`# Gene Family`, pattern="\\|")) -> Gene_abundance
  
  if (sum_up == T){
    read_delim("Gene_ID_metadata.txt", delim="=", col_names = F) %>% select(X1) -> Meta
    sapply(X = Meta$X1, function(x){strsplit(x = x,split=" ")}) -> Name_and_gene
    Meta_proper = tibble()
    for (Entry in Name_and_gene){
      Name = Entry[1]
      Gene = paste(Entry[2:(length(Entry)-1)],collapse="_")
      NG = tibble(Name, Gene)
      colnames(NG) = c("ID","Gene")
      Meta_proper = rbind(Meta_proper,NG)
    }
    Meta_proper %>% filter(ID %in% Gene_abundance$`# Gene Family`) -> Meta_proper
    Gene_abundance %>% mutate(`# Gene Family`=Meta_proper$Gene) %>% group_by(`# Gene Family`) %>% summarise_all(sum) -> Gene_abundance
  }
  
  Header = colnames(Gene_abundance) 
  Header = as.vector(sapply(Header,function(x){strsplit(x ,"_")[[1]][1]}))
  New_header = c(Gene_abundance$`# Gene Family`, "ID")
  Transposed_genes = as_tibble(t(Gene_abundance[,2:length(Header)]))
  Transposed_genes %>% mutate(ID = Header[2:length(Header)]) -> Transposed_genes
  colnames(Transposed_genes) = New_header
  
  
  Translation = read_tsv(file = "ID_barcode.txt")
  Translation %>% filter(!(sampleID=="LLDeep_1442" & Project=="G83489"))-> Translation
  Translation$feceID = sapply(Translation$feceID,function(x){return(strsplit(x ,"_")[[1]][2])})
  Translation %>% arrange(feceID) %>% filter(feceID %in% Transposed_genes$ID) -> Translation
  Transposed_genes %>% arrange(ID) %>% filter(ID %in% Translation$feceID) %>% mutate(ID = Translation$sampleID) -> Transposed_genes
  
  Transposed_genes %>% filter(ID=="LLDeep_1433") %>% sample_n(1) -> Choice
  Transposed_genes %>% filter(!(ID=="LLDeep_1433")) -> Transposed_genes
  Tranposed_genes = rbind(Transposed_genes,Choice)
  
  
  return(Transposed_genes)
}
Metabolome_transformation = function(x, samples_row=T, method="hfmin", missing_filter=0){
  x[x==0] <- NA
  log10(x) -> x
  x <- as_tibble(scale(x))
  
  #if samples are in columns transpose
  if (!samples_row){ x=as.data.frame(t(x)) }
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>1){ stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 1") }
  col_ori=ncol(x)
  col_keep=colnames(x)[colSums(!is.na(x))/nrow(x)>missing_filter]
  x=x[,col_keep]
  my_num_removed=col_ori-ncol(x)
  warning (paste(my_num_removed, "columns removed due to many missing values"))
  if (method=="zero"){ x[is.na(x)]=0 
  } else if (method=="min"){ x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], function(a) ifelse(is.na(a), min(a, na.rm = TRUE), a)) 
  } else if (method=="hfmin"){
    x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], function(a) ifelse(is.na(a), min(a, na.rm = TRUE)/2, a))
  } else if (method=="mean"){
    x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], function(a) ifelse(is.na(a), mean(a, na.rm = TRUE), a))
  } else if (method=="median"){
    x[sapply(x, is.numeric)] <- lapply(x[sapply(x, is.numeric)], function(a) ifelse(is.na(a), median(a, na.rm = TRUE), a))
  } else if (method=="none"){
    x=x
  } else{
    stop("\n Hey! \n Method parameters for imputing missing values per column are: \n min -> min (col), hfmin -> min(col)/2, mean -> mean (col), median -> median(col), mice -> using mice package [Multivariate Imputation By Chained Equations] defaults \n")
  }
  return(as.data.frame(x))
}
Make_heatmap = function(Results, Name, Scale=3, size_rows=5){
  Results %>% mutate(Beta=as.numeric(Beta)) -> Results
  Results %>% mutate(log10Pval = -log10(Pvalue)) %>% select(Metabolite, Regressor, log10Pval) %>% spread(Metabolite, log10Pval) %>%filter(! is.na(Regressor)) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results
  Results %>% select(Metabolite, Regressor, Beta) %>% spread(Metabolite, Beta) %>% filter(! is.na(Regressor)) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results2
  
  wide_results_Annotation = wide_results2
  
  wide_results_Annotation[wide_results_Annotation < "0"] <- "-"
  wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
  
  
  # set colors
  #breaksList = seq(0, +max(abs(wide_results)),by=2 * max(abs(wide_results))/55)
  #cols <- colorRampPalette(c("Darkblue","White","Darkred"))(length(breaksList))
  #cols[seq(length(cols)/2-1,length(cols)/2+1)] <- "White"

  pheatmap(wide_results, display_numbers = wide_results_Annotation, fontsize_number = 7, fontsize_col = 7,fontsize_row = size_rows) -> Fig_0#,border_color = "#EEEEEE",na_col = "white",
           #color = cols,treeheight_row = 0, treeheight_col = 0,legend_labels = "log(p-value)",
           #breaks = breaksList ) -> Fig_0
  ggsave(filename = paste(c("Model_stats/", Name, ".pdf"), collapse = ""), Fig_0, scale = Scale)
  
  return(Fig_0)
}
Compute_GFR = function(x){
  x %>% mutate(GFR= ifelse(Gender==1, ifelse(Creatinine<= 0.7,144*(Creatinine/0.7)^-0.329 *(0.993)^Age, 144*(Creatinine/0.7)^-1.209 * (0.993)^Age), NA)) -> x
  x %>% mutate(GFR= ifelse(Gender==0, ifelse(Creatinine<=0.9,141*(Creatinine/0.9)^-0.411 * (0.993)^Age, 141*(Creatinine/0.9)^-1.209 * (0.993)^Age), GFR)) -> x
  return(x)
}
Do_enrichment = function(DE_genes, Metab, ensembl){
  DE = filter(DE_genes, Metabolite == Metab)
  Genes = DE$Regressor
  
  foo <- getBM(attributes=c('ensembl_gene_id',
                            'external_gene_name'),
               filters = 'ensembl_gene_id',
               values = Genes,
               mart = ensembl)
  Genes = foo$external_gene_name
  
  
  
  #dbs <- listEnrichrDbs()
  dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
  enriched <- enrichr(Genes, dbs)
  Process = as_tibble(enriched[["GO_Biological_Process_2018"]])
  Component = as_tibble(enriched[["GO_Cellular_Component_2018"]])
  Function = as_tibble(enriched[["GO_Molecular_Function_2018"]])
  
  return(list(Process, Component, Function))
  
}
Presence <- function(x, na.rm = FALSE){ sapply(X = x, FUN= function(x){ if(as.numeric(x)>0){ return(1) }else{return(0)}}) }
Percentage <- function(x, na.rm = FALSE){ x/dim(Species_abundance)[[1]] }
Selection <- function(x, Threshold=0.1, na.rm = FALSE){ x > Threshold }

##############################
#########ML functions#########
##############################
RMSE = function(Y_hat, Y){ sqrt(sum((Y_hat - Y)^2)/length(Y)) }
error = function(Y_hat, Y){ 1 - (sum(Y_hat == Y)/length(Y))}
Split_Train_Test = function(X,Y){
  round(dim(X)[1]*0.7) -> Training_size
  dim(X)[1] -> End
  #Predictors
  X[1:Training_size,] -> Train
  X[(Training_size+1):End,] -> Test
  #Dependent
  Y[1:Training_size] -> Label_Train
  Y[(Training_size+1):End] -> Label_Test
  
  return(list(Train, Test, Label_Train, Label_Test))
}
Fit_RandomForest = function(X,Y){
  X %>% mutate(Metabolite = Y) -> Train
  #Fit Model
  model1 <- randomForest(Metabolite ~ ., data = Train, importance = TRUE, mtry=(max(floor(ncol(X)/3), 1)*2))
  predTrain = predict(model1, Train)
  Loss_train = RMSE(predTrain,Train$Metabolite)
  
  as.data.frame(importance(model1)) %>% rownames_to_column() %>% as_tibble() %>%
    arrange(desc(`%IncMSE`)) -> Variable_Importance
  ggplot(Variable_Importance) + geom_bar(aes(x=rowname,y=`%IncMSE`),stat = "identity") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(list(Variable_Importance,model1))
  
}
Fit_RandomForest_classification = function(X,Y, quartiles){
  X %>% mutate(Metabolite = Y) -> Train
  #Train$quartile <- with(Train, cut(Metabolite, 
  #                               breaks=quantile(Metabolite, probs=seq(0,1, by=0.25), na.rm=TRUE), 
  #                                include.lowest=TRUE))
  Train %>% mutate(quartile = as.factor(quartiles)) %>% select(-Metabolite) -> Train
  model1 <- randomForest(quartile ~ ., data = Train, importance = TRUE)
  predTrain = predict(model1, Train)
  Loss_train = error(predTrain,Train$quartile)
  
  as.data.frame(importance(model1)) %>% rownames_to_column() %>% as_tibble() %>%
  arrange(desc(MeanDecreaseGini)) -> Variable_Importance
  return(list(Variable_Importance,model1))
  
}
Fit_lasso = function(X,Y){
  #findCorrelation(cor(X,method = "spearman"), cutoff = 0.75) -> Variables_correlated
  #X[-Variables_correlated] -> X
  
  #List_split = Split_Train_Test(X,Y)
  #Train = List_split[[1]]
  #Train_y = List_split[[3]]
  #Test = List_split[[2]]
  #Test_y = List_split[[4]]
  Train = X
  X = as.matrix(Train)
  Train_y = Y
  #X_test = as.matrix(Test)
  
  #Cross Validation 
  cvfit = cv.glmnet(X, Train_y, nfolds = 20,type.measure="mse",standardize=T, alpha=1)
  #Param
  Param_best = coef(cvfit, s = "lambda.min")
  #Fit in Test
  #Prediction = predict(cvfit, newx = X_test, s = "lambda.min")
  #Loss = RMSE(Prediction, Test_y)
  
  #Variation explained Given the lambda that minimizes the error function
  Explained_variance = cvfit$glmnet.fit$dev.ratio[which(cvfit$glmnet.fit$lambda == cvfit$lambda.min)]
  return(list(cvfit, Param_best, Explained_variance))
}

FeatureSelection_RandomForest = function(Input_x, Dependent_metabolites, Covariates, Fold, Method="Regression", Min_Pvalue = 0.2){
  RETURN = NULL
  set.seed(111)
  Iteration = colnames(Dependent_metabolites)[2:length(colnames(Dependent_metabolites))]
  Iteration[!Iteration=="quartile"] -> Iteration
  for (Metabolite in Iteration){
    
    print(Metabolite)
    if (Metabolite == "TMAO"){ Dependent_metabolites = filter(Dependent_metabolites, TMAO<20)  }
    
    Quantile_input = select(Dependent_metabolites,Metabolite) 
    colnames(Quantile_input) = c("Meta")
    Quantile_input$quartile <- with(Quantile_input, cut(Meta, breaks=quantile(Meta, probs=seq(0,1, by=0.25), na.rm=TRUE), include.lowest=TRUE))
    Dependent_metabolites %>% mutate(quartile = Quantile_input$quartile) -> Dependent_metabolites
    ##10-fold cross-validation, using OoB MSE as error measure
    flds <- createFolds(Input_x, k = Fold, list = TRUE, returnTrain = FALSE)
    Min_loss = vector()
    
    for (i in seq(Fold)){
      Test = Input_x[ flds[[i]], ]
      Train_i = Input_x[-flds[[i]], ]
      
      Train_i %>% select(-c("Age","Gender","BMI","Creatinine","GFR")) -> Train
      Input_train = Match_dataset(Train, Metabolites = Dependent_metabolites, Covariates = Covariates) 
      Results_train = Iterate_Metabolites(select(Input_train[[1]],c(ID,Metabolite)),Input_train[[2]],Input_train[[3]], "MGS_species_cv", FDR_iter = 0)
      #Results_train = Iterate_Metabolites(select(Input_train[[1]],c(ID,`y-butyrobetaine`)),Input_train[[2]],Input_train[[3]], "MGS_species_cv", FDR_iter = 0)
      Results_train %>% filter(Pvalue < Min_Pvalue) -> Feature_selection
      #Add covaraites?
      Input_train[[3]] %>% select(Feature_selection$Regressor) %>% cbind(select(Input_train[[2]],-ID)) -> Train
      
      if (Method == "Regression"){
        RF_train = Fit_RandomForest(Train, as_vector(select(Input_train[[1]], Metabolite)))
        Model_RF = RF_train[[2]]
        Input_test = Match_dataset(Test, Metabolites = Dependent_metabolites, Covariates = Covariates) 
        predTest = predict(Model_RF, Input_test[[3]])
        Loss_test = RMSE(predTest,as_vector(select(Input_test[[1]],Metabolite)))
      } else {
        RF_train = Fit_RandomForest_classification(Train, as_vector(select(Input_train[[1]], Metabolite)), as_vector(select(Input_train[[1]], quartile)))
        Model_RF = RF_train[[2]]
        Input_test = Match_dataset(Test, Metabolites = Dependent_metabolites, Covariates = Covariates) 
        predTest = predict(Model_RF, Input_test[[3]])
        Loss_test = error(predTest,as_vector(select(Input_test[[1]], Metabolite)))
        
      }
      
      if (length(Min_loss) == 0){
        Min_loss = Loss_test
        Model = Model_RF
        Importance = RF_train[[1]]
        Variables = Feature_selection$Regressor
      }else if (Loss_test < Min_loss){
        Min_loss = Loss_test
        Model = Model_RF
        Importance = RF_train[[1]]
        Variables = Feature_selection$Regressor
      }
      
    }
    if (Method == "Regression"){
      Input_x %>% select(Variables,ID) -> Final_set
      Input_final = Match_dataset(Final_set, Metabolites = Dependent_metabolites, Covariates = Covariates) 
      
      RF = Fit_RandomForest(select(Input_final[[3]],-ID), as_vector(select(Input_final[[1]], Metabolite)))
      Model_RF = RF[[2]]
      Importance = RF[[1]]
      
    }else{
      Input_x %>% select(Variables,ID) -> Final_set
      Input_final = Match_dataset(Final_set, Metabolites = Dependent_metabolites, Covariates = Covariates) 
      RF = Fit_RandomForest_classification(select(Input_final[[3]], -ID), as_vector(select(Input_final[[1]], Metabolite)), as_vector(select(Input_final[[1]], quartile)))
      Model_RF = RF[[2]]
      Importance = RF[[1]]
      
    }
    
    OUT = list(RF, Importance, Metabolite)
    RETURN = list(RETURN, OUT)
  }
  return(RETURN)
}


####################################
#########Preprocess functions#######
####################################
Prepare_Metagenome_species = function(Species_abundance, Metabolites, Covariates){
  List_output = Match_dataset(Metabolites, Covariates, Species_abundance)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Species_abundance=List_output[[3]] 
  colnames(Species_abundance) = Filter_name_meta(colnames(Species_abundance))
  return(list(Metabolites_meta, Covariates_meta, Species_abundance))
}
Prepare_Metagenome_genus = function(Genus_abundance, Metabolites, Covariates){
  List_output = Match_dataset(Metabolites, Covariates, Genus_abundance)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Genus_abundance=List_output[[3]] 
  Genus_abundance = Transformation_composition(Genus_abundance)
  colnames(Genus_abundance) = Filter_name_meta(colnames(Genus_abundance))
  return(list(Metabolites_meta, Covariates_meta, Genus_abundance))
}
Prepare_Metagenome_family = function(Family_abundance, Metabolites, Covariates){
  List_output = Match_dataset(Metabolites, Covariates, Family_abundance)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Family_abundance=List_output[[3]] 
  Family_abundance = Transformation_composition(Family_abundance)
  colnames(Family_abundance) = Filter_name_meta(colnames(Family_abundance))
  return(list(Metabolites_meta, Covariates_meta, Family_abundance))
}
Prepare_Metagenome_phylum = function(Phylum_abundance, Metabolites, Covariates){
  List_output = Match_dataset(Metabolites, Covariates, Phylum_abundance)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Phylum_abundance=List_output[[3]] 
  Phylum_abundance = Transformation_composition(Phylum_abundance)
  colnames(Phylum_abundance) = Filter_name_meta(colnames(Phylum_abundance))
  return(list(Metabolites_meta, Covariates_meta, Phylum_abundance))
}

Prepare_Pathways = function(Pathway_abundance, Metabolites, Covariates){
  Pathway_abundance %>% select(-c(UNMAPPED, UNINTEGRATED)) -> Pathway_abundance
  List_output = Match_dataset(Metabolites, Covariates, Pathway_abundance)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Pathway_abundance=List_output[[3]] 
  return(list(Metabolites_meta, Covariates_meta, Pathway_abundance))
}
Prepare_diet = function(Diet_regressors, Metabolites, Covariates){
  
  #Include Diet
  #Diet_regressors %>% mutate(ID = LL_ID_corrected) %>% arrange(ID) %>%  select(-LL_ID_corrected) -> Diet_regressors
  List_output = Match_dataset(Metabolites, Covariates, Diet_regressors)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Regressors_trs=List_output[[3]] 
  return(list(Metabolites_meta, Covariates_meta, Regressors_trs))
}
Prepare_clinical = function(Clinical_Questionaries, Metabolites, Covariates){
  Remove = colnames(Clinical_Questionaries)[grepl(pattern = "\\?", colnames(Clinical_Questionaries))]
  Remove = c(Remove,colnames(Clinical_Questionaries)[grepl(pattern = "[Aa]dded", colnames(Clinical_Questionaries))])
  Remove = c(Remove, colnames(Clinical_Questionaries)[grepl(pattern = "technical", colnames(Clinical_Questionaries))])
  Remove = c(Remove, colnames(Clinical_Questionaries)[grepl(pattern = "_", colnames(Clinical_Questionaries))])
  Remove = c(Remove, c("age at first visit", "female=1; male=0", "birthyear", "Allergy", "Body Mass Index (kg/M^2)", "type1diabetes", "type2diabetes"))
  Remove = unique(Remove)
  Clinical_Questionaries %>% select(-Remove) %>% select_if(is.numeric) %>% mutate(ID=Clinical_Questionaries$ID)  -> Regressors
  
  
  IDs = Regressors$ID
  Regressors_trs = Metabolome_transformation(select(Regressors,-ID)) #Transformation_composition(Cluster_abundance2)
  as_tibble(Regressors_trs) %>% mutate(ID= IDs) -> Regressors_trs
  
  
  #Regressors_trs = Transformation_logarithmic(Regressors)
  List_output = Match_dataset(Metabolites, Covariates, Regressors_trs)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Regressors_trs=List_output[[3]] 
  return(list(Metabolites_meta, Covariates_meta, Regressors_trs))
  
}

Prepare_metabolome = function(Metabolome, Metabolites, Covariates){
  IDs = Metabolome$ID
  Metabolome = Metabolome_transformation(select(Metabolome,-ID), samples_row=T, method="hfmin", missing_filter=0.3)
  List_output = Match_dataset(Metabolites, Covariates, mutate(Metabolome,ID=IDs))
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Regressors_trs=List_output[[3]] 
  return(list(Metabolites_meta, Covariates_meta, Regressors_trs))
}

Prepare_cutC_shortbred = function(Cluster_abundance, Metabolites, Covariates){
  Cluster_abundance2 = Cluster_abundance
  #Cluster_abundance %>% select(-c(Hits,TotMarkerLength)) %>% spread(Family, Count) -> Cluster_abundance2
  #Cluster_abundance2 %>% summarise_if(is_numeric, sum) %>% t() %>% as.data.frame() %>%
  #  rownames_to_column() %>% filter(V1 > 0) -> Keep_columns 
  # Cluster_abundance2 %>% select(c("ID",Keep_columns$rowname)) -> Cluster_abundance2
  
  List_output = Match_dataset(Metabolites, Covariates, Cluster_abundance2)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Cluster_abundance2=List_output[[3]] 
  IDs = Cluster_abundance2$ID
  Pathway_abundance = Metabolome_transformation(select(Cluster_abundance2,-ID)) #Transformation_composition(Cluster_abundance2)
  as_tibble(Pathway_abundance) %>% mutate(ID= IDs) -> Pathway_abundance
  return(list(Metabolites_meta, Covariates_meta, Pathway_abundance))
}
Prepare_humman2_gene = function(Gene_abundance, Dependent_metabolites, Covariates, SUM=F){
  if (SUM == T){
    Gene_abundance2 = Gene_abundance_process(Gene_abundance, sum_up = T)
    IDs = Gene_abundance2$ID
    Gene_abundance2 = Metabolome_transformation(select(Gene_abundance2,-ID))
    mutate(Gene_abundance2, ID = IDs) -> Gene_abundance2
  }else{
    Gene_abundance2 = Gene_abundance_process(Gene_abundance, sum_up = F)
    Gene_abundance2 %>% summarise_if(is.numeric, sum) %>% t() ->Total_gene_abundance 
    Total_gene_abundance %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% filter(V1 == 0) -> Total_gene_abundance_zero
    Gene_abundance2 %>% select(-Total_gene_abundance_zero$rowname) -> Gene_abundance2
    IDs = Gene_abundance2$ID
    Gene_abundance2 = Metabolome_transformation(select(Gene_abundance2,-ID))
    mutate(Gene_abundance2, ID = IDs) -> Gene_abundance2
  }
  
  List_output = Match_dataset(Dependent_metabolites, Covariates, Gene_abundance2)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Gene_abundance2=List_output[[3]] 
  return(list(Metabolites_meta, Covariates_meta, Gene_abundance2))
}

Prepare_expression = function(Expression_Counts, Dependent_metabolites, Covariates){
  colnames(Expression_counts)[2:length(colnames(Expression_counts))] -> Column_ID 
  Expression_counts$ID -> New_colnames
  Expression_counts %>% select(-ID) %>% t() %>% as_tibble()  -> Expression_counts2
  colnames(Expression_counts2) = New_colnames
  Expression_counts2[,1:51800] -> Expression_counts2
  Metabolome_transformation(Expression_counts2,missing_filter = 0.2)  -> Expression_counts2
  Expression_counts2 %>% as_tibble() %>% mutate(ID = Column_ID) -> Expression_counts2
  
  List_output = Match_dataset(Dependent_metabolites, Covariates, Expression_counts2)
  Metabolites_meta = List_output[[1]]; Covariates_meta=List_output[[2]]; Regressors_trs=List_output[[3]] 
}

Filter_abundance = function(M,M2, threshold=0.2, min_abundance = 0){
  #Filter species not seen in at least treshold*100% of samples
  apply(select(M, -ID), MARGIN=2, FUN = function(x){ (sum(x != 0)/length(x)) >= min_abundance }) -> Filter_vector1
  apply(select(M2, -ID), MARGIN=2, FUN = function(x){ (sum(x != 0)/length(x)) >= min_abundance }) -> Filter_vector2
  Filter_vector = c(unique(c(names(Filter_vector1), names(Filter_vector2))), "ID")
    
  M %>% select(one_of(Filter_vector)) -> M
  M2 %>% select(one_of(Filter_vector)) -> M2
  
  #Filter species with mean relative abundance lower than min_abundance
  apply(select(M, -ID), MARGIN=1, FUN = function(x){ sum(x) }) -> Total_reads_sample ; apply(select(M2, -ID), MARGIN=1, FUN = function(x){ sum(x) }) -> Total_reads_sample2
  as_tibble(select(M, -ID)/Total_reads_sample) %>% summarise_all(mean) -> All_means ; as_tibble(select(M2, -ID)/Total_reads_sample2) %>% summarise_all(mean) -> All_means2
  colnames(select(M, -ID))[as_vector(All_means) > min_abundance] -> Taxa_not_to_filter ; colnames(select(M2, -ID))[as_vector(All_means) > min_abundance] -> Taxa_not_to_filter2
  M %>% select(c(ID, Taxa_not_to_filter)) -> M ; M2 %>% select(c(ID, Taxa_not_to_filter2)) -> M2
  
  
  return(list(M,M2))
}


############################################
####Compare MSS abundances with 16S#########
############################################

Choose_replicate = function(Linking_table){
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
  Remove_columns = colnames(Count_table)[grepl("NOTAX", colnames(Count_table))]
  Count_table %>% select(-Remove_columns) %>% select(-rootrank.Root) -> Count_table
  return(Count_table)
}
Change_format = function(Names){
  strsplit(x = Names,split = "\\.")[[1]] ->y
  if (y[2] != "" ){
    return(y[2])
  }else{
    return(y[3])}
}
Compare_16S_MSS = function(Count_table, Level, MSS_Abundance, Transform = "no", correlation="kendall"){
  if (Transform == "yes"){
    MSS_Abundance = Transformation_composition(MSS_Abundance)
  }
  ###Part1, process of 16S as in TMAOmetabolites_16S
  Count_table$SampleID -> Sample_ID
  colnames(Count_table)[grepl(Level, colnames(Count_table))] -> Level_filter
  Count_table %>% select(Level_filter) %>% mutate(SampleID =  Sample_ID) -> Taxa_table
  #Merging SeqID and cohort ID
  Linking_table %>% arrange(X2) %>% filter(X2 %in% Taxa_table$SampleID) -> Linking_table_taxa
  Taxa_table %>% arrange(SampleID) %>% filter(SampleID %in% Linking_table$X2) %>% 
    mutate(SampleID = Linking_table_taxa$X1) %>% arrange(SampleID) -> Taxa_table
  
  ###Part2, match taxa and individuals in 16S and MSS data
  
  #Change Taxa name format in 16S
  ID = Taxa_table$SampleID
  select(Taxa_table,-SampleID) -> Taxa_table
  sapply(colnames(Taxa_table), Change_format) %>% as.vector() -> Columns
  colnames(Taxa_table) = Columns
  #Taxa_table[-which(Columns==paste(c("unknown",Level),collapse=""))] -> Taxa_table
  Remove = which(grepl(Level, Columns,ignore.case=T))
  if (length(Remove > 0)){ Taxa_table[-Remove] -> Taxa_table }
  Taxa_table %>% mutate(ID = ID) -> Taxa_table
  if (Transform == "yes"){
    Taxa_table = Transformation_composition(Taxa_table)
  }
  #Match ID and format MSS
  MSS_Abundance %>% filter(ID %in% Taxa_table$ID) -> MSS_Abundance
  colnames(MSS_Abundance) = c("ID",as.vector(sapply(colnames(MSS_Abundance), function(x){strsplit(x,"__")[[1]][2]}))[2:length(colnames(MSS_Abundance))])
  
  #Just keep the Taxa name that are identical
  #Keep in 16S
  match(colnames(MSS_Abundance),colnames(Taxa_table)) -> Keep_columns
  Taxa_table[Keep_columns[!is.na(Keep_columns)]] -> Taxa_table
  Taxa_table %>% filter(ID %in% MSS_Abundance$ID) -> Taxa_table
  #Keep in MSS
  match(colnames(Taxa_table),colnames(MSS_Abundance)) -> Keep_columns
  MSS_Abundance[Keep_columns[!is.na(Keep_columns)]] -> MSS_Abundance
  
  #Merging in a single Data structure
  MSS_Abundance %>% gather(Taxa, measurement, 2:dim(MSS_Abundance)[2], factor_key=TRUE) %>% mutate(Source = "MSS") -> Ready_MSS
  Taxa_table %>% gather(Taxa, measurement, 2:dim(MSS_Abundance)[2], factor_key=TRUE) %>% mutate(Source = "16S") -> Ready_16S
  rbind(Ready_MSS,Ready_16S) -> Comparative_data
  Comparative_data %>% spread(Source,measurement) -> Comparative_data
  
  ###Part3 compute correlations
  All_correlations = tibble()
  for (TAXA in unique(Comparative_data$Taxa)){
    Comparative_data %>% filter(Taxa== TAXA) -> Comp_d
    if (Transform == "yes"){ 
      cor(x=Comp_d$`16S`,y=Comp_d$MSS) -> C
    }else{
      cor(x=Comp_d$`16S`,y=Comp_d$MSS, method = correlation) -> C
    }
    Result = tibble(Taxa = TAXA, Correlation=C)
    All_correlations = rbind(All_correlations, Result)
  }
  return(All_correlations)
}
