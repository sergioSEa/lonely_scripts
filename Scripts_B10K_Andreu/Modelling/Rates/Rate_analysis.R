setwd("~/Documents/Paper/Research_Assitant/dnds")
library(tidyverse)
library(ggrepel)
library(ape)
library(phytools)
library(stringr)
library(nlme)
library(ggforce)
library(ggpubr)
library(scales)
library(cowplot)
library(caper)
require(MASS)
require(dplyr)
library(RColorBrewer)
#Input: TSV long format, one line for syn, other line for non-syn.
#Three sets: 13 Mitogenes, 13 Mitogenes in Passeriformes, 100 Background
#From tip and from root


##########All Organisms########
#Tip
Data_tips = read_tsv("hyphy_tip.tsv")
iqtree_tips = read_tsv("iqtree_tips.tsv")
rbind(Data_tips,iqtree_tips) -> Data_tips
#Root-tip
Data_root = read_tsv("hyphy_root.tsv")
iqtree_root = read_tsv("iqtree_root.tsv")
rbind(Data_root,iqtree_root) -> Data_root


##########All Organisms - IQTREE codon########
#Tip
Data_tips2 = read_tsv("hyphy_tip.tsv")
iqtree_tips2 = read_tsv("Complex/iqtree_tips.tsv")
rbind(Data_tips2,iqtree_tips2) -> Data_tips2
#Root-tip
Data_root2 = read_tsv("hyphy_root.tsv")
iqtree_root2 = read_tsv("Complex/iqtree_root.tsv")
rbind(Data_root2,iqtree_root2) -> Data_root2

##########All Organisms - IQTREE no third position########
#Tip
Data_tips3 = read_tsv("hyphy_tip.tsv")
iqtree_tips3 = read_tsv("Complex/remove_third/iqtree_tips.tsv")
rbind(Data_tips3,iqtree_tips3) -> Data_tips3
#Root-tip
Data_root3 = read_tsv("hyphy_root.tsv")
iqtree_root3 = read_tsv("Complex/remove_third/iqtree_root.tsv")
rbind(Data_root3,iqtree_root3) -> Data_root3



##########Passeriformes data##########
#Tip
Data_pass_tips = read_tsv("Passeriformes/hyphy_merged_tip.tsv")
iqtree_pass_tips = read_tsv("Passeriformes/iqtree_tips.tsv")
rbind(Data_pass_tips,iqtree_pass_tips) -> Data_pass_tips
#Root-tip
Data_pass_root = read_tsv("Passeriformes/hyphy_merged_tiproot.tsv")
iqtree_pass_root = read_tsv("Passeriformes/iqtree_root.tsv")
rbind(Data_pass_root,iqtree_pass_root) -> Data_pass_root


##########NonPasseriformes data##########
#Tip
Data_nopass_tips = read_tsv("NoPasseriformes/hyphy_tip.tsv")
iqtree_nopass_tips = read_tsv("NoPasseriformes/iqtree_tips.tsv")
rbind(Data_nopass_tips,iqtree_nopass_tips) -> Data_nopass_tips
#Root-tip
Data_nopass_root = read_tsv("NoPasseriformes/hyphy_root.tsv")
iqtree_nopass_root = read_tsv("NoPasseriformes/iqtree_root.tsv")
rbind(Data_nopass_root,iqtree_nopass_root) -> Data_nopass_root


#######Background##############
#Tip
Data_back_tips = read_tsv("Background/Background_tip_rate.tsv")
iqtree_back_tips = read_tsv("Background/iqtree_tips.tsv")
rbind(Data_back_tips,iqtree_back_tips) -> Data_back_tips
#Root-tip
Data_back_root = read_tsv("Background/Background_root_rate.tsv")
iqtree_back_root = read_tsv("Background/iqtree_root.tsv")
rbind(Data_back_root,iqtree_back_root) -> Data_back_root

##########Background Passeriformes data##########
#Tip
Data_backpass_tips = read_tsv("Passeriformes/Background/hyphy_merged_tip.tsv")
iqtree_backpass_tips = read_tsv("Passeriformes/Background/iqtree_tips.tsv")
rbind(Data_backpass_tips,iqtree_backpass_tips) -> Data_backpass_tips
#Root-tip
Data_backpass_root = read_tsv("Passeriformes/Background/hyphy_merged_tiproot.tsv")
iqtree_backpass_root = read_tsv("Passeriformes/Background/iqtree_root.tsv")
rbind(Data_backpass_root,iqtree_backpass_root) -> Data_backpass_root




##########Background NoPasseriformes data##########
#Tip
Data_nobackpass_tips = read_tsv("NoPasseriformes/Background/rate_tip.tsv")
iqtree_nobackpass_tips = read_tsv("NoPasseriformes/Background/iqtree_tips.tsv")
rbind(Data_nobackpass_tips,iqtree_nobackpass_tips) -> Data_nobackpass_tips
#Root-tip
Data_nobackpass_root = read_tsv("NoPasseriformes/Background/rate_root.tsv")
iqtree_nobackpass_root = read_tsv("NoPasseriformes/Background/iqtree_root.tsv")
rbind(Data_nobackpass_root,iqtree_nobackpass_root) -> Data_nobackpass_root











#####################################################ANALYSIS
Perform_analisys = function(data,outcome,Set){
  #Distribution
  data %>% ggplot() + geom_density(aes(x=nosyn)) + facet_wrap(~Protein) + theme_bw()
  data %>% ggplot() + geom_density(aes(x=syn)) + facet_wrap(~Protein) + theme_bw()
  data %>% ggplot() + geom_density(aes(x=ML_branch)) + facet_wrap(~Protein) + theme_bw()
  #Plot of Mass and latitude effects
  ###Plots Mass
  data  %>%  ggplot(aes(y=nosyn, x=log(unsexed_mass))) + geom_point() + facet_wrap(~Protein) + theme_bw() -> nosyn_mass
  data  %>%  ggplot(aes(y=syn, x=log(unsexed_mass))) + geom_point() + facet_wrap(~Protein) + theme_bw() -> syn_mass
  data  %>%  ggplot(aes(y=ML_branch, x=log(unsexed_mass))) + geom_point() + facet_wrap(~Protein) + theme_bw() -> subs_mass
  ###Plots latitude
  data  %>%  ggplot(aes(y=nosyn, x=log(Latitude))) + geom_point() + facet_wrap(~Protein) + theme_bw()  -> nosyn_lat
  data  %>%  ggplot(aes(y=syn, x=log(Latitude))) + geom_point() + facet_wrap(~Protein) + theme_bw() -> syn_lat
  data  %>%  ggplot(aes(y=ML_branch, x=log(Latitude))) + geom_point() + facet_wrap(~Protein) + theme_bw() -> subs_lat
  
  ggsave(plot = nosyn_mass,filename = paste(outcome,"_nosyn_vs_mass.pdf"))
  ggsave(plot =syn_mass,filename =paste(outcome,"_syn_vs_mass.pdf"))
  ggsave(plot =subs_mass,filename =paste(outcome,"_rates_vs_mass.pdf"))
  ggsave(plot =nosyn_lat,filename =paste(outcome,"_nosyn_vs_lat.pdf"))
  ggsave(plot =syn_lat,filename =paste(outcome,"_syn_vs_lat.pdf"))
  ggsave(plot =subs_lat,filename =paste(outcome,"_rates_vs_lat.pdf"))
  
  ##PGLS
  ALL_ID = read_csv("../Analysis_Mtc_number/info/NAME-ID.txt")
  t2 <- read.tree("../Analysis_Mtc_number/info/B10K.new")

  RESULTS = tibble("Lambda" = numeric(), "Latitude_Estimation" = numeric(), "Mass_Estimation"= numeric(),"Latitude_pvalue" = numeric(),"Mass_pvalue"=numeric(), "Model" =character(), "Protein"  = character())
  RESULTS = as.data.frame(RESULTS)
  
  
  ####Multiple regression PGLS
  for (Protein_n in unique(data$Protein)){
    Data = filter(data, Protein == Protein_n)
    Data %>%  dplyr::select(ID, syn, nosyn,Latitude,unsexed_mass,ML_branch) %>% as.data.frame() -> mm
    comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
    #MAKE LOG SO RESIDUALS ARE NORMAL
    Model_syn = pgls(syn ~ log(abs(Latitude)) + log(unsexed_mass), comp, lambda=1)
    Model_nosyn = pgls(nosyn ~ log(abs(Latitude)) + log(unsexed_mass), comp, lambda=1)
    Model_branch = pgls(ML_branch ~ log(abs(Latitude)) + log(unsexed_mass), comp, lambda=1)
    n =0
    L <- vector("list",3)
    L[[1]] = Model_syn
    L[[2]] = Model_nosyn
    L[[3]] = Model_branch
    for (model_n in 1:length(L)){
      
      model = L[[model_n]]
      lambda = model$param[[2]]
      Latitude_est = summary(model)$coefficients[2]
      mass_est = summary(model)$coefficients[3]
      Latitude_pval =  summary(model)$coefficients[11]
      Mass_pval =  summary(model)$coefficients[12]
      n = n +1
      if (n == 1){
        Model = "Synonymous"
      }else if (n ==2){
        Model = "NonSynonymous"
      }else{
        Model = "ML_branch"
      }
      Line = data.frame(lambda,Latitude_est,mass_est,Latitude_pval,Mass_pval,Model,Protein_n)
      names(Line) = c("Lambda","Latitude_Estimation","Mass_Estimation", "Latitude_pvalue", "Mass_pvalue", "Model", "Protein")
      RESULTS = rbind(RESULTS,Line)  
    }
    
  }
  
  RESULTS_hy = as_tibble(RESULTS)
  #write.csv(RESULTS_hy,paste(ourcome,"_models.csv"))

  
  col_wes = brewer.pal(12, "Paired")
  col_wes = c(col_wes,"#000000","#C0C0C0")
  
  L_corrected = p.adjust(RESULTS_hy$Latitude_pvalue,method="fdr")
  M_corrected = p.adjust(RESULTS_hy$Mass_pvalue,method="fdr")
  RESULTS_hy %>% mutate(FDR_lat= L_corrected ,FDR_mass= M_corrected) -> RESULTS_hy
  
  
  if (Set == "Background"){
    RESULTS_hy %>% ggplot(aes(y=-log(Latitude_pvalue),x=Model)) + geom_boxplot() + geom_point(size=3) + 
      theme_bw() + geom_hline(yintercept=-log(0.05), linetype="dashed",  color = "blue", size=1) +  scale_color_manual(values=col_wes)+
      geom_hline(yintercept=-log(0.05/dim(RESULTS_hy)[1]), linetype="dashed",  color = "red", size=1) -> Latitude_dN_dS
    
    RESULTS_hy %>% ggplot(aes(y=-log(Mass_pvalue),x=Model)) + geom_boxplot() + geom_point(size=3) +
      theme_bw() + geom_hline(yintercept=-log(0.05), linetype="dashed",  color = "blue", size=1) + scale_color_manual(values=col_wes)+
      geom_hline(yintercept=-log(0.05/dim(RESULTS_hy)[1]), linetype="dashed",  color = "red", size=1) -> mass_dN_dS
    
    RESULTS_hy %>% ggplot(aes(y=Latitude_Estimation,x=Model)) + geom_sina(aes(col=FDR_lat < 0.05),size=3) + 
      theme_bw() + geom_hline(yintercept= 0, linetype="dashed",  color = "black", size=1)  -> Latitude_effect
    RESULTS_hy %>% ggplot(aes(y=Mass_Estimation,x=Model)) + geom_sina(aes(col=FDR_mass < 0.05),size=3) + 
      theme_bw() + geom_hline(yintercept= 0, linetype="dashed",  color = "black", size=1)  -> Mass_effect

  }else{
  
  RESULTS_hy %>% ggplot(aes(y=-log(Latitude_pvalue),x=Model)) + geom_boxplot() + geom_point(aes(col=Protein,shape=Protein=="TOTAL"),size=3) + 
    theme_bw() + geom_hline(yintercept=-log(0.05), linetype="dashed",  color = "blue", size=1) +  scale_color_manual(values=col_wes)+
    geom_hline(yintercept=-log(0.05/dim(RESULTS_hy)[1]), linetype="dashed",  color = "red", size=1) -> Latitude_dN_dS
  
  RESULTS_hy %>% ggplot(aes(y=-log(Mass_pvalue),x=Model)) + geom_boxplot() + geom_point(aes(col=Protein,shape=Protein=="TOTAL"),size=3) + 
    theme_bw() + geom_hline(yintercept=-log(0.05), linetype="dashed",  color = "blue", size=1) + scale_color_manual(values=col_wes)+
    geom_hline(yintercept=-log(0.05/dim(RESULTS_hy)[1]), linetype="dashed",  color = "red", size=1) -> mass_dN_dS
  

  RESULTS_hy %>% ggplot(aes(y=Latitude_Estimation,x=Model)) + geom_point(aes(col=FDR_lat < 0.05, shape=Protein=="TOTAL"),size=3) + 
    theme_bw() + geom_hline(yintercept= 0, linetype="dashed",  color = "black", size=1)  -> Latitude_effect
  RESULTS_hy %>% ggplot(aes(y=Mass_Estimation,x=Model)) + geom_point(aes(col=FDR_mass < 0.05,shape=Protein=="TOTAL"),size=3) + 
    theme_bw() + geom_hline(yintercept= 0, linetype="dashed",  color = "black", size=1)  -> Mass_effect
  }
  

  RESULTS_hy %>% ggplot(aes(x=Mass_Estimation,y=-log(Mass_pvalue))) + geom_point(aes(col=FDR_mass < 0.05),size=3) + 
    theme_bw() + geom_hline(yintercept= -log(0.05), linetype="dashed",  color = "black", size=1) + 
    geom_vline(xintercept= 0, linetype="dashed",  color = "black", size=1) + facet_wrap(~Model, scales = "free") -> effect_pval_mass
  RESULTS_hy %>% ggplot(aes(x=Latitude_Estimation,y=-log(Latitude_pvalue))) + geom_point(aes(col=FDR_lat < 0.05),size=3) + 
    theme_bw() + geom_hline(yintercept= -log(0.05), linetype="dashed",  color = "black", size=1) + 
    geom_vline(xintercept= 0, linetype="dashed",  color = "black", size=1) + facet_wrap(~Model, scales = "free")  -> effect_pval_lat
  
  
   
  ggsave(plot = Latitude_dN_dS,filename = paste(outcome,"_Latitude_pval.pdf",collapse=""))
  ggsave(plot = mass_dN_dS,filename = paste(outcome,"_mass_pval.pdf",collapse=""))

  
  ggsave(plot = effect_pval_lat,filename = paste(outcome,"_Latitude_effectpval.pdf",collapse=""))
  ggsave(plot = effect_pval_mass,filename = paste(outcome,"_mass_effectpval.pdf",collapse=""))
  
  ggsave(plot = Latitude_effect,filename = paste(outcome,"_Latitude_effect.pdf",collapse=""))
  ggsave(plot = Mass_effect,filename = paste(outcome,"_mass_effect.pdf",collapse=""))
  
 
  
  print(Latitude_dN_dS)
  print(mass_dN_dS)
  print(Latitude_effect)
  print(Mass_effect)
  print(effect_pval_mass)
  print(effect_pval_lat)
  
}
##################################################




############### COMPLETE
Data_root %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_root
Data_root %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_root

Data_tips %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_tips
Data_tips %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_tips
##############

############## PASSERIFORMES
Data_pass_root %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_pass_root
Data_pass_root %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_pass_root

Data_pass_tips %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_pass_tips
Data_pass_tips %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_pass_tips
#############

############## NOPASSERIFORMES
Data_nopass_root %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_nopass_root
Data_nopass_root %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_nopass_root

Data_nopass_tips %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_nopass_tips
Data_nopass_tips %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_nopass_tips
#############


############# Background
Data_back_root %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_back_root
Data_back_root %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_back_root

Data_back_tips %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_back_tips
Data_back_tips %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_back_tips

############


############# Background Passeriformes
Data_backpass_root %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_backpass_root
Data_backpass_root %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_backpass_root

Data_backpass_tips %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_backpass_tips
Data_backpass_tips %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_backpass_tips

############

############# Background Passeriformes
Data_nobackpass_root %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_nobackpass_root
Data_nobackpass_root %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_nobackpass_root

Data_nobackpass_tips %>%dplyr::select(ID_organism, Protein, Rate, Tree, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(Tree, Rate) -> Data_nobackpass_tips
Data_nobackpass_tips %>% mutate(dNdS_ratio= nosyn/syn, syn= as.numeric(syn), nosyn= as.numeric(nosyn)) -> Data_nobackpass_tips

############




Perform_analisys(Data_root,"../Figures/Rates/Total_root","MTC")
Perform_analisys(Data_tips,"../Figures/Rates/Total_tips","MTC")


Perform_analisys(Data_pass_root,"../Figures/Rates/Pass_root", "MTC")
Perform_analisys(Data_pass_tips,"../Figures/Rates/Pass_tips", "MTC")

Perform_analisys(Data_nopass_root,"../Figures/Rates/nopass_root","MTC")
Perform_analisys(Data_nopass_tips,"../Figures/Rates/nopass_tips", "MTC")

Perform_analisys(Data_back_root,"../Figures/Rates/Back_root","Background")
Perform_analisys(Data_back_tips,"../Figures/Rates/Back_tips", "Background")


Perform_analisys(Data_backpass_root,"../Figures/Rates/Backpass_root","Background")
Perform_analisys(Data_backpass_tips,"../Figures/Rates/Backpass_tips", "Background")


Perform_analisys(Data_nobackpass_root,"../Figures/Rates/noBackpass_root","Background")
Perform_analisys(Data_nobackpass_tips,"../Figures/Rates/noBackpass_tips", "Background")




###############################################################
##Differences in IQTREE substitution models with complex models
iqtree_root %>% mutate(IQTREE_v = "Regular") -> iqtree_root
iqtree_root3 %>% mutate(IQTREE_v = "No_third") -> iqtree_root3
iqtree_root2 %>% mutate(IQTREE_v = "Complex") -> iqtree_root2
rbind(iqtree_root, iqtree_root2, iqtree_root3) -> all_iqtree
all_iqtree %>% dplyr::select(ID_organism, Protein, Rate, IQTREE_v, unsexed_mass, Latitude, Order, Latin_name,ID) %>% spread(IQTREE_v, Rate) -> all_iqtree
all_iqtree %>% ggplot(aes(x=Regular, y=Complex, col=Order)) + geom_point() + theme_bw()+ guides(col=FALSE) -> p1
all_iqtree %>% ggplot(aes(x=Regular, y=No_third, col=Order)) + geom_point() + theme_bw()+ guides(col=FALSE) -> p2
all_iqtree %>% ggplot(aes(x=Complex, y=No_third, col=Order)) + geom_point() + theme_bw() -> p3
ggarrange(p1,p2,p3, labels = "AUTO", ncol = 3, nrow = 1, common.legend = T, legend= "bottom")

cor.test(all_iqtree$Regular, all_iqtree$Complex, method="pearson")
cor.test(all_iqtree$Regular, all_iqtree$No_third, method="pearson")
cor.test(all_iqtree$No_third, all_iqtree$Complex, method="pearson")

##PGLS
ALL_ID = read_csv("../Analysis_Mtc_number/info/NAME-ID.txt")
t2 <- read.tree("../Analysis_Mtc_number/info/B10K.new")

RESULTS = tibble("Lambda" = numeric(), "Latitude_Estimation" = numeric(), "Mass_Estimation"= numeric(),"Latitude_pvalue" = numeric(),"Mass_pvalue"=numeric(), "Model" =character(), "Protein"  = character())
RESULTS = as.data.frame(RESULTS)

####Multiple regression PGLS
for (Protein_n in unique(all_iqtree$Protein)){
  Data = filter(all_iqtree, Protein == Protein_n)
  Data %>%  dplyr::select(ID,Latitude,unsexed_mass, Regular,Complex,No_third) %>% as.data.frame() -> mm
  comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
  #MAKE LOG SO RESIDUALS ARE NORMAL
  Model_regular = pgls(Regular ~ log(abs(Latitude)) + log(unsexed_mass), comp, lambda=1)
  Model_complex = pgls(Complex ~ log(abs(Latitude)) + log(unsexed_mass), comp, lambda=1)
  Model_Third = pgls(No_third ~ log(abs(Latitude)) + log(unsexed_mass), comp, lambda=1)
  n =0
  L <- vector("list",3)
  L[[1]] = Model_regular
  L[[2]] = Model_complex
  L[[3]] = Model_Third
  for (model_n in 1:length(L)){
    
    model = L[[model_n]]
    lambda = model$param[[2]]
    Latitude_est = summary(model)$coefficients[2]
    mass_est = summary(model)$coefficients[3]
    Latitude_pval =  summary(model)$coefficients[11]
    Mass_pval =  summary(model)$coefficients[12]
    n = n +1
    if (n == 1){
      Model = "Regular"
    }else if (n ==2){
      Model = "Complex"
    }else{
      Model = "No_third"
    }
    Line = data.frame(lambda,Latitude_est,mass_est,Latitude_pval,Mass_pval,Model,Protein_n)
    names(Line) = c("Lambda","Latitude_Estimation","Mass_Estimation", "Latitude_pvalue", "Mass_pvalue", "Model", "Protein")
    RESULTS = rbind(RESULTS,Line)  
  }
  
}
RESULTS_hy = as_tibble(RESULTS)
col_wes = brewer.pal(12, "Paired")
col_wes = c(col_wes,"#000000","#C0C0C0")

L_corrected = p.adjust(RESULTS_hy$Latitude_pvalue,method="fdr")
M_corrected = p.adjust(RESULTS_hy$Mass_pvalue,method="fdr")
RESULTS_hy %>% mutate(FDR_lat= L_corrected ,FDR_mass= M_corrected) -> RESULTS_hy

RESULTS_hy %>% ggplot(aes(y=-log(Latitude_pvalue),x=Model)) + geom_boxplot() + geom_point(size=3) + 
    theme_bw() + geom_hline(yintercept=-log(0.05), linetype="dashed",  color = "blue", size=1) +  scale_color_manual(values=col_wes)+
    geom_hline(yintercept=-log(0.05/dim(RESULTS_hy)[1]), linetype="dashed",  color = "red", size=1) -> Latitude_IQTREE
  
RESULTS_hy %>% ggplot(aes(y=-log(Mass_pvalue),x=Model)) + geom_boxplot() + geom_point(size=3) +
    theme_bw() + geom_hline(yintercept=-log(0.05), linetype="dashed",  color = "blue", size=1) + scale_color_manual(values=col_wes)+
    geom_hline(yintercept=-log(0.05/dim(RESULTS_hy)[1]), linetype="dashed",  color = "red", size=1) -> mass_IQTREE
  
RESULTS_hy %>% ggplot(aes(y=Latitude_Estimation,x=Model)) + geom_sina(aes(col=FDR_lat < 0.05),size=3) + 
    theme_bw() + geom_hline(yintercept= 0, linetype="dashed",  color = "black", size=1)  -> Latitude_effect
RESULTS_hy %>% ggplot(aes(y=Mass_Estimation,x=Model)) + geom_sina(aes(col=FDR_mass < 0.05),size=3) + 
    theme_bw() + geom_hline(yintercept= 0, linetype="dashed",  color = "black", size=1)  -> Mass_effect

RESULTS_hy %>% ggplot(aes(x=Mass_Estimation,y=-log(Mass_pvalue))) + geom_point(aes(col=FDR_mass < 0.05),size=3) + 
  theme_bw() + geom_hline(yintercept= -log(0.05), linetype="dashed",  color = "black", size=1) + 
  geom_vline(xintercept= 0, linetype="dashed",  color = "black", size=1) + facet_wrap(~Model, scales = "free") -> effect_pval_mass
RESULTS_hy %>% ggplot(aes(x=Latitude_Estimation,y=-log(Latitude_pvalue))) + geom_point(aes(col=FDR_lat < 0.05),size=3) + 
  theme_bw() + geom_hline(yintercept= -log(0.05), linetype="dashed",  color = "black", size=1) + 
  geom_vline(xintercept= 0, linetype="dashed",  color = "black", size=1) + facet_wrap(~Model, scales = "free")  -> effect_pval_lat


#############################################
###########################################


##Passeriformes show a lower effect on mass because have a distribution of smaller values than the rest
iqtree_tips %>% filter(! Order == "Passeriformes") %>% arrange(unsexed_mass) %>% distinct(ID_organism, unsexed_mass) %>% mutate(Pass=F) -> no_pass
iqtree_tips %>% filter(Order == "Passeriformes") %>% arrange(unsexed_mass) %>% distinct(ID_organism,unsexed_mass) %>% mutate(Pass=T)  -> pass

R = rbind(no_pass,pass)
R %>% ggplot() + geom_histogram(aes(x=log(unsexed_mass),fill=Pass)) + theme_bw()       





#################################################







### The correlation of Rates and Mitochondrial copy number or BMR

#####################
data = read_delim("~/Documents/Paper/Research_Assitant/Analysis_Mtc_number/F_Quantification.tsv",delim = "\t")
dplyr::select(data,-ID) -> data
names(data)[1] = "ID"
###Naming###
MIXED_TISSUE = c(NA,"Heart/liver/muscle", "multiple embryos", "Heart (possibly muscle too)", "Liver (possibly heart too)", "muscle; mixed tissue sample?", "unknown tissue sample", "tissue", "muscle (10); liver (11)", "Unknown","blood/muscle/liver", "heart (possibly muscle too)", "liver (possibly heart too)","heart, liver, muscle","liver, muscle (breast)")
MUSCLES = c("muscle, pectoral", "muscle-crop", "Muscle", "muscle, leg","tissue is presumed to be muscle, but not recorded in catalog")
LIVER = c("liver")
BLOOD = c("hematocrit")
HEART= c("Heart")
data %>% mutate(`Tissue type` = ifelse(`Tissue type` %in% MIXED_TISSUE, "mixed tissue sample", `Tissue type`)) %>%
  mutate(`Tissue type` = ifelse(`Tissue type` %in% MUSCLES , "muscle", `Tissue type`)) %>%
  mutate(`Tissue type` = ifelse(`Tissue type` %in% LIVER , "Liver", `Tissue type`)) %>%
  mutate(`Tissue type` = ifelse(`Tissue type` %in% HEART , "heart", `Tissue type`)) %>%
  mutate(`Tissue type` = ifelse(`Tissue type` %in% BLOOD , "blood", `Tissue type`)) -> data
data %>% filter(Coverage_ncl>=0.9 & max_sixe > 10000) -> data_soft
#######################
data_soft %>% arrange(ID) %>% filter(ID %in% Data_tips$ID_organism) -> data_soft

##PGLS
ALL_ID = read_csv("../Analysis_Mtc_number/info/NAME-ID.txt")
t2 <- read.tree("../Analysis_Mtc_number/info/B10K.new")
#Compare MTC number estimation against all mitochondrial proteins
for (Protein_n in unique(Data_root$Protein)){
  Data = filter(Data_root, Protein == Protein_n)
  Add_data = filter(data_soft, ID %in% Data$ID_organism)
  Data %>% arrange(ID_organism) %>% filter(ID_organism %in% data_soft$ID)%>% mutate(Mitochondrial_number=Add_data$`Ratio of Median Depth`) -> Data
  
  Data %>%  dplyr::select(Mitochondrial_number, ML_branch,ID) %>% as.data.frame() -> mm
  comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
  #MAKE LOG SO RESIDUALS ARE NORMAL
  Model_p = pgls(log(Mitochondrial_number) ~ log(ML_branch), comp, lambda="ML")
  p_p = summary(Model_p)$coefficients[8]
  p_v = summary(Model_p)$coefficients[2]
  p_r = summary(Model_p)$r.squared
  figure = ggplot(data=Data,aes(x=log(ML_branch), y=log(Mitochondrial_number))) + geom_point() + theme_bw() + geom_smooth(method='lm', formula= y~x)
  print(c(Protein_n,p_p,p_v, p_r))
  plot(figure)
  }

for (Protein_n in unique(Data_root$Protein)){
  Data = filter(Data_root, Protein == Protein_n)
  Add_data = filter(data_soft, ID %in% Data$ID_organism)
  Data %>% arrange(ID_organism) %>% filter(ID_organism %in% data_soft$ID)%>% mutate(Mitochondrial_number=Add_data$`Ratio of Median Depth`) -> Data
  
  Data %>%  dplyr::select(Mitochondrial_number, ML_branch,ID) %>% as.data.frame() -> mm
  comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
  #MAKE LOG SO RESIDUALS ARE NORMAL
  Model_p = pgls(log(ML_branch) ~ log(Mitochondrial_number) , comp, lambda=0.01)
  p_p = summary(Model_p)$coefficients[8]
  p_v = summary(Model_p)$coefficients[2]
  p_r = summary(Model_p)$r.squared
  figure = ggplot(data=Data,aes(y=log(ML_branch), x=log(Mitochondrial_number))) + geom_point() + theme_bw() + geom_smooth(method='lm', formula= y~x)
  print(c(Protein_n,p_p,p_v, p_r))
  plot(figure)
}







###Compare against BMR
#1
BMR_measures = read_tsv("~/Documents/Paper/Research_Assitant/Bibliography/Bird_BMR-Fristoe_et_al.tsv")
#2
BMR_measures2 = read_csv("~/Documents/Paper/Research_Assitant/Bibliography/BMR_uyeda_et_al.csv")
#3
BMR_measures3 = read_tsv("~/Documents/Paper/Research_Assitant/Bibliography/Lodono_et_al_BMR.tsv")
#4
BMR_measures4 = read_delim("~/Documents/Paper/Research_Assitant/Bibliography/Bushuev_BRM.csv",delim=";")
change = remove_sd(BMR_measures4$`BMR (ml O2/min)`)
change2 = remove_sd(BMR_measures4$`M (g)`)
BMR_measures4 %>% mutate(`BMR (ml O2/min)`=change, `M (g)` = change2) -> BMR_measures4
####
BMR_measures %>% filter(Species %in% gsub("_"," ",Data_tips$Latin_name)) %>% arrange(Species) ->b_overl
BMR_measures2 %>% filter(X1 %in% Data_tips$Latin_name) %>% arrange(X1) -> b_overl2
BMR_measures3 %>% filter(Species %in%gsub("_"," ",Data_tips$Latin_name)) %>% arrange(Species) -> b_overl3
BMR_measures4 %>% filter(Species %in%gsub("_"," ",Data_tips$Latin_name)) %>% arrange(Species) -> b_overl4

b_overl %>% mutate(BMR_watts = conversion(`BMR(mlO2/hour)`)) -> b_overl
b_overl4  %>% mutate("BMR(watts)" = conversion(`BMR (ml O2/min)`)*60) -> b_overl4

b_overl3 %>% mutate("Mass(g)"=`Mean mass (g)`, "BMR(watts)"=`BMR (watts; = 20.1 J/ml O2)`, "Source"="Londono") %>%
  dplyr::select(`Mass(g)`,`BMR(watts)`, Species,Source) %>% arrange(Species) -> overlap_londono
b_overl %>% mutate("BMR(watts)"=BMR_watts, "Source"="Fristoe") %>%  dplyr::select(`Mass(g)`,`BMR(watts)`, Species,Source)%>% arrange(Species) %>% 
  filter(!Species %in% overlap_londono$Species)-> overlap_Fristoe
b_overl4 %>% mutate("Source"="Bushuev", "Mass(g)"=`M (g)` ) %>%  dplyr::select(`Mass(g)`,`BMR(watts)`, Species,Source)%>% arrange(Species) %>%
  filter(! Species %in% overlap_londono$Species) %>% filter(! Species %in% overlap_Fristoe$Species) -> overlap_bushuev

rbind(overlap_londono,overlap_Fristoe,overlap_bushuev) -> data_BMR
data_BMR[-c(14,27),] -> data_BMR
data_BMR %>% ggplot(aes(x=log(`Mass(g)`),y=log(`BMR(watts)`))) + geom_point()



for (Protein_n in unique(Data_root$Protein)){
  Data = filter(Data_root, Protein == Protein_n)
  Data %>% arrange(Latin_name) %>% filter( gsub("_"," ",Latin_name) %in% data_BMR$Species) -> rc

  data_BMR %>% filter(Species %in%  gsub("_"," ",rc$Latin_name)) -> data_BMR2 
  Data %>% arrange(Latin_name)  %>% filter( gsub("_"," ",Latin_name) %in% data_BMR2$Species) %>% mutate(BMR=data_BMR2$`BMR(watts)`, Mass=data_BMR2$`Mass(g)`) -> root_check
  
  regular_model = lm(log(ML_branch) ~ log(BMR/Mass), root_check)
  p_r = summary(regular_model)$coefficients[8]
  
  root_check %>%  dplyr::select(BMR,ID, ML_branch,Mass) %>% as.data.frame() -> mm
  comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
  #MAKE LOG SO RESIDUALS ARE NORMAL
  Model_p = pgls(log(ML_branch) ~ log(BMR/Mass), comp, lambda="ML")
  p_p = summary(Model_p)$coefficients[8]
  figure = ggplot(data=root_check,aes(x=log(ML_branch), y=log(BMR/Mass))) + geom_point() + theme_bw() + geom_smooth(method='lm', formula= y~x)
  print(c(Protein_n,p_p,p_r))
  plot(figure)
}


for (Protein_n in unique(Data_root$Protein)){
  Data = filter(Data_root, Protein == Protein_n)
  Data %>% arrange(Latin_name) %>% filter( gsub("_"," ",Latin_name) %in% data_BMR$Species) -> rc
  
  data_BMR %>% filter(Species %in%  gsub("_"," ",rc$Latin_name)) -> data_BMR2 
  Data %>% arrange(Latin_name)  %>% filter( gsub("_"," ",Latin_name) %in% data_BMR2$Species) %>% mutate(BMR=data_BMR2$`BMR(watts)`, Mass=data_BMR2$`Mass(g)`) -> root_check
  
  regular_model = lm(log(ML_branch) ~ log(BMR/Mass), root_check)
  p_r = summary(regular_model)$coefficients[8]
  
  root_check %>%  dplyr::select(BMR,ID, ML_branch,Mass) %>% as.data.frame() -> mm
  comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
  #MAKE LOG SO RESIDUALS ARE NORMAL
  Model_p = pgls(log(BMR/Mass) ~ log(ML_branch), comp, lambda="ML")
  p_p = summary(Model_p)$coefficients[8]
  figure = ggplot(data=root_check,aes(y=log(ML_branch), x=log(BMR/Mass))) + geom_point() + theme_bw() + geom_smooth(method='lm', formula= y~x)
  print(c(Protein_n,p_p,p_r))
  plot(figure)
}


