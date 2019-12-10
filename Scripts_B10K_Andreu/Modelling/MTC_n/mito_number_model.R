setwd("~/Documents/Paper/Research_Assitant/Analysis_Mtc_number")
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


#size text axis
TS=15


#DATA, after running the script Fuse_trait_n_Quan.py
data = read_delim("F_Quantification.tsv",delim = "\t")
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

high_avail = c("Aq","Vt","Om","Sc")
low_avail = c("Fr","Ne","Gr","Pl","In") #Invertebrades have lower metabolism
data %>% mutate(main_diet = ifelse(main_diet %in% high_avail, "CARN", ifelse(main_diet %in% low_avail, "LEAF", main_diet))) -> data

#####

#Only muscle data, 166 samples out of the 166 muscle samples available in the consortium's data
data %>% group_by(`Tissue type`) %>% summarise(n())

#see the effect of Coverage and length
data %>% ggplot(aes(y=`Ratio of Median Depth`, x=Coverage_mtc)) + geom_point() + theme_bw()+ 
  labs(y='Estimation mitochondrial number', x= "Mitochondrial genome breadth of coverage") -> mtc_effect
data %>% ggplot(aes(y=`Ratio of Median Depth`, x=Coverage_ncl)) + geom_point() + theme_bw()   + labs(y='Estimation mitochondrial number', x= "Nuclear genome breath of coverage") +
  geom_vline(xintercept = 0.9, linetype="dashed",color = "red", size=1.5) -> ncl_effect
data %>% ggplot(aes(y=`Ratio of Median Depth`, x=max_sixe)) + geom_point() + theme_bw()  + 
  geom_vline(xintercept = 10000, linetype="dashed",color = "red", size=1.5)+ 
  labs(y='Estimation mitochondrial number', x= "Mitochondrial size (bp)") -> size_effect

ggarrange(mtc_effect, ncl_effect, size_effect, labels = c("A", "B", "C"), ncol = 3, nrow = 1) -> Fig0
ggsave(filename = "Biases.pdf", plot = Fig0)
ggsave(filename = "Biases.png", plot = Fig0)

#Filter: 26 smaples.  Other filter: mtc, ncl>0.8, size> 10K --> 162/166
data %>% filter(Coverage_mtc>=0.8 & Coverage_ncl>=0.9 & max_sixe > 14000 & max_sixe < 19000 ) -> data_soft
data %>% filter(Coverage_ncl>=0.9 & max_sixe > 10000) -> data_soft
#183/176



#data %>% ggplot(aes(x=factor(migration),y=`Ratio of Median Depth`)) + geom_boxplot() + theme_bw()

###Non-phylogeny

#Mass trend
data_soft %>% ggplot(aes(y=log(`Ratio of Median Depth`),x=log(unsexed_mass))) + geom_point(size=4) +
theme_bw() + geom_smooth(method='lm',formula=y~x) + theme(axis.title.x = element_text(size=TS), axis.title.y = element_text(size=TS),axis.text=element_text(size=12))+
  labs(y='Estimation mitochondrial number', x= "ln(Body mass (g))")  -> mass_trend

summary(lm(log(`Ratio of Median Depth`) ~ log(unsexed_mass), data))


#Latitude trend
data_soft %>% ggplot(aes(y=log(`Ratio of Median Depth`),x=log(abs(Latitude)))) + geom_point(size=4) +
  theme_bw() + geom_smooth(method='lm',formula=y~x) + theme(axis.title.x = element_text(size=TS), axis.title.y = element_text(size=TS),axis.text=element_text(size=12)) + 
  labs(y='Estimation mitochondrial number', x= "ln(abs(Latitude (degrees)))")-> latitude_trend
summary(lm(log(`Ratio of Median Depth`) ~ log(abs(Latitude)), data))

summary(lm(log(`Ratio of Median Depth`) ~ log(abs(Latitude)):factor(Order), data))

#Temperature trend
data_soft %>% ggplot(aes(y=log(`Ratio of Median Depth`),x=Temperature)) + geom_point(size=4) +
  theme_bw() + geom_smooth(method='lm',formula=y~x) + theme(axis.title.x = element_text(size=TS), axis.title.y = element_text(size=TS),axis.text=element_text(size=12))-> temperature_trend
summary(lm(log(`Ratio of Median Depth`) ~ Temperature, data))


#Temperature and latitude correlate

data_soft %>% ggplot(aes(x=abs(Latitude),y=Temperature)) + geom_point(size=4) + theme_bw() + geom_smooth(method='lm',formula=y~x) + theme(axis.title.x = element_text(size=TS), axis.title.y = element_text(size=TS),axis.text=element_text(size=12))-> sup_T_Latitude

ggsave("T_vs_Latitude.pdf",sup_T_Latitude)
ggsave("T_vs_Latitude.png",sup_T_Latitude)
#Multiple regression
summary(lm(log(`Ratio of Median Depth`) ~ log(abs(Latitude)) + log(unsexed_mass) , data))

summary(lm(log(`Ratio of Median Depth`) ~ log(unsexed_mass) + Temperature , data))

#Colinearity
data_soft %>% ggplot(aes(y=log(unsexed_mass),x=log(abs(Latitude)))) + geom_point(size=4) +
  theme_bw() + geom_smooth(method='lm',formula=y~x) + theme(axis.title.x = element_text(size=TS), axis.title.y = element_text(size=TS),axis.text=element_text(size=12))-> colinearity

summary(lm(log(unsexed_mass) ~ log(abs(Latitude)) + log(unsexed_mass) , data))
#R square 0.001551
###

ggarrange(mass_trend,latitude_trend, labels = c("A", "B"), ncol = 2, nrow = 1) -> Fig1 
#ggsave("../Figures/Factors_miton_number.pdf", Fig1)
#ggsave("../Figures/Factors_miton_number.png", Fig1)


ggsave("../Figures/Sup/Sup_colinearity.pdf", colinearity)

##############################
##PGLS
ALL_ID = read_csv("info/NAME-ID.txt")
t2 <- read.tree("info/B10K.new")

library(caper)
require(MASS)
require(dplyr)

####Mass PGLS
data_soft %>%  dplyr::select(ID, `Ratio of Median Depth`,unsexed_mass) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_0 = pgls(log(`Ratio of Median Depth`) ~ log(unsexed_mass), comp, lambda="ML")
summary(Model_0)

####Latitude PGLS
data_soft %>%  dplyr::select(ID, `Ratio of Median Depth`,Latitude) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_1 = pgls(log(`Ratio of Median Depth`) ~ log(abs(Latitude)), comp, lambda="ML")
summary(Model_1)

####Multiple regression PGLS
data_soft %>%  dplyr::select(ID, `Ratio of Median Depth`,Latitude,unsexed_mass) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_2 = pgls(log(`Ratio of Median Depth`) ~ log(abs(Latitude)) + log(unsexed_mass), comp, lambda="ML")
summary(Model_2)
#Adjusted R-squared: 0.1011

###MR with other Latitude values
data_soft %>%  dplyr::select(ID, `Ratio of Median Depth`,Latitude,unsexed_mass, Max_latitude, Min_latitude) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_2.1 = pgls(log(`Ratio of Median Depth`) ~ log(abs(Max_latitude)) + log(unsexed_mass), comp, lambda="ML")
summary(Model_2.1)

data_soft %>%  dplyr::select(ID, `Ratio of Median Depth`,Latitude,unsexed_mass, Max_latitude, Min_latitude) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_2.2 = pgls(log(`Ratio of Median Depth`) ~ log(abs(Min_latitude)) + log(unsexed_mass), comp, lambda="ML")
summary(Model_2.2)



##Multiple with Temperature

data_soft %>%  dplyr::select(ID, `Ratio of Median Depth`,Latitude,unsexed_mass,Temperature) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_25 = pgls(log(`Ratio of Median Depth`) ~ log(unsexed_mass) + Temperature, comp, lambda="ML")
summary(Model_25)


data_soft %>%  dplyr::select(ID, `Ratio of Median Depth`,Latitude,unsexed_mass) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_3 = pgls(log(`Ratio of Median Depth`) ~ log(abs(Latitude)):log(unsexed_mass), comp, lambda="ML")
summary(Model_3)


####Multiple regression only in Passeriformes####
data_soft %>% filter(Order=="Passeriformes") %>% dplyr::select(ID, `Ratio of Median Depth`,Latitude,unsexed_mass) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_4 = pgls(log(`Ratio of Median Depth`) ~ log(abs(Latitude))+log(unsexed_mass), comp, lambda="ML")
summary(Model_4)


#non pass
data_soft %>% filter(!Order=="Passeriformes") %>% dplyr::select(ID, `Ratio of Median Depth`,Latitude,unsexed_mass) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_5 = pgls(log(`Ratio of Median Depth`) ~ log(abs(Latitude))+log(unsexed_mass), comp, lambda="ML")
summary(Model_5)


###Test for normality in the residuals
shapiro.test(Model_0$phyres)
shapiro.test(Model_1$phyres)
shapiro.test(Model_2$phyres)
shapiro.test(Model_4$phyres)




#Mass trend
Mass_data = tibble(Y=Model_0$y, X=Model_0$x[,2], Fitted=Model_0$fitted )
Mass_data %>% ggplot(aes(y=Y,x=X)) + geom_point(size=4) + theme_bw() + geom_line(aes(x=X, y=Fitted),col="blue") +
  theme(axis.title.x = element_text(size=TS), axis.title.y = element_blank()) + labs(x="ln(Body Mass (grams))",y='Estimation mitochondrial number') -> mass_trend
#Latitude trend
Lat_data = tibble(Y=Model_1$y, X=Model_1$x[,2], Fitted=Model_1$fitted )
Lat_data %>% ggplot(aes(y=Y,x=X)) + geom_point(size=4) + theme_bw() + geom_line(aes(x=X, y=Fitted),col="blue") + 
  theme(axis.title.x = element_text(size=TS), axis.title.y = element_blank()) + labs(x="ln(Latitude (degrees))",y='Estimation mitochondrial number') ->Lat_trend


ggarrange(mass_trend,Lat_trend, labels = c("A", "B"), ncol = 1, nrow = 2) -> Fig1 
annotate_figure(Fig1, left=text_grob("ln(Estimation of mitochondria number)", size=TS,rot = 90)) -> Fig1
ggsave("../Figures/Factors_miton_number.pdf", Fig1)

ggsave("../Figures/Factors_miton_number.png", Fig1)




###########################Distribtuion data###############
Lat_data = read_csv("/Users/Sergio/Documents/Paper/Research_Assitant/Distribution_data2.csv")
Lat_data %>% distinct(Sp,.keep_all = TRUE) %>% arrange(Sp) -> Lat_data
dim(data_soft) 
data_soft %>% filter(Latin_name %in% Lat_data$Sp) %>% arrange(Latin_name) -> data_soft
Lat_data %>% filter(Sp %in% data_soft$Latin_name) -> Lat_data
data_soft %>% mutate(Max_latitude=Lat_data$mx, Min_latitude=Lat_data$mn, Lat_Range= abs(Max_latitude-Min_latitude)) -> data_soft


data_soft %>% filter(Lat_Range < 30) -> s_data_soft

s_data_soft %>%  dplyr::select(ID, `Ratio of Median Depth`,Latitude,unsexed_mass) %>% as.data.frame() -> mm
comp <- comparative.data(force.ultrametric(t2, method=c("nnls","extend")), mm, ID, vcv=TRUE, vcv.dim=3)
Model_2 = pgls(log(`Ratio of Median Depth`) ~ log(abs(Latitude)) + log(unsexed_mass), comp, lambda="ML")
summary(Model_2)


##############################################################################
######################## BMR correlation with estimator ######################
##############################################################################


#Using BMR data from Fristoe et al.
BMR_measures = read_tsv("~/Documents/Paper/Research_Assitant/Bibliography/Bird_BMR-Fristoe_et_al.tsv")

#Mass-BMR and mass-specific trends
summary(lm(log(`BMR(mlO2/hour)`) ~ log(`Mass(g)`) , data=BMR_measures))
summary(lm(log(`BMR(mlO2/hour)`/`Mass(g)`) ~ log(`Mass(g)`) , data=BMR_measures)) 

#Comparing to data project
#Get Data from Modellin.T
data_soft %>% filter(Latin_name_1 %in% BMR_measures$Species) %>% arrange(Latin_name_1) -> d_overl
BMR_measures %>% filter(Species %in% data_soft$Latin_name_1) %>% arrange(Species) -> b_overl

b_overl %>% mutate(Mass_mine=d_overl$unsexed_mass, Estimation=d_overl$`Ratio of Median Depth`,latitude = d_overl$Latitude, tissue= d_overl$`Tissue type`, Order= d_overl$Order, Temperature = d_overl$Temperature) -> b_overl

b_overl %>% ggplot() + geom_point(aes(x=log(Mass_mine),y=log(Estimation))) + facet_grid(~tissue)
b_overl %>% ggplot() + geom_point(aes(x=log(`BMR(mlO2/hour)`/`Mass(g)`),y=log(Estimation))) + facet_grid(~tissue)

b_overl %>% ggplot() + geom_point(aes(x=log(`BMR(mlO2/hour)`),y=log(abs(latitude)))) 

summary(lm(log(`BMR(mlO2/hour)`/`Mass(g)`) ~ log(`Mass(g)`) , data=b_overl)) -> model_Fristoe

summary(lm(model_Fristoe$residuals ~ log(abs(b_overl$latitude))))



#Using BMR data from Uyeda et al. The evolution of energetic scaling across the vertebrate tree of life
BMR_measures2 = read_csv("~/Documents/Paper/Research_Assitant/Bibliography/BMR_uyeda_et_al.csv")
#Mass-BMR and mass-specific trends
summary(lm(lnBMR ~ lnMass, data=BMR_measures2)) 
#Compare to estimates
data_soft %>% filter(Latin_name %in% BMR_measures2$X1) %>% arrange(Latin_name) -> d_overl2
BMR_measures2 %>% filter(X1 %in% data_soft$Latin_name) %>% arrange(X1) -> b_overl2

b_overl2 %>% mutate(Mass_mine=d_overl2$unsexed_mass, Estimation=d_overl2$`Ratio of Median Depth`,latitude = d_overl2$Latitude, tissue= d_overl2$`Tissue type`, Order =d_overl2$Order, Other_naming=d_overl2$Latin_name_1, Temperature = d_overl2$Temperature) -> b_overl2

b_overl2 %>% ggplot() + geom_point(aes(x=lnBMR,y=log(Estimation))) + facet_grid(~tissue)


#Using BMR from Lodono et al.

BMR_measures3 = read_tsv("~/Documents/Paper/Research_Assitant/Bibliography/Lodono_et_al_BMR.tsv")
#Mass-BMR and mass-specific trends
summary(lm(log(`BMR (watts; = 20.1 J/ml O2)`/`Mean mass (g)`) ~ log(`Mean mass (g)`), data=BMR_measures3)) 
#Compare to estimates
data_soft %>% filter(Latin_name_1 %in% BMR_measures3$Species) %>% arrange(Latin_name) -> d_overl3
BMR_measures3 %>% filter(Species %in% data_soft$Latin_name_1) %>% arrange(Species) -> b_overl3

b_overl3 %>% mutate(Mass_mine=d_overl3$unsexed_mass, Estimation=d_overl3$`Ratio of Median Depth`,latitude = d_overl3$Latitude, tissue= d_overl3$`Tissue type`, Order=d_overl3$Order, Temperature = d_overl3$Temperature) -> b_overl3

b_overl3  %>% ggplot(aes(x=log(`BMR (watts; = 20.1 J/ml O2)`/`Mean mass (g)`),y=log(Estimation),col=Order)) + geom_point() +
  theme_bw() +  geom_smooth(method='lm',formula=y~x, se=F) -> Londono_comparison1
b_overl3  %>% filter(Order=="Passeriformes") %>% ggplot(aes(x=log(`BMR (watts; = 20.1 J/ml O2)`/`Mean mass (g)`),y=log(Estimation),col=Order)) + geom_point() +
  theme_bw() +  geom_smooth(method='lm',formula=y~x, se=F) -> Londono_comparison2
ggarrange(Londono_comparison1,Londono_comparison2, labels = c("A", "B"), ncol = 2, nrow = 1) -> BMR_comparison_londono
ggsave("BMR_comparison_londono.pdf",BMR_comparison_londono)



#Using BMR from Bushuev
remove_sd = function(DATA){
  output = vector()
  for (item in DATA){
      v = strsplit(item,"[\\\\]|[^[:print:]]")[[1]][[1]] #REMOVE COMA by POINT
      if (str_count(v,"\\.") == 2){ v = strsplit(v,"\\.")[[1]][[1]] }
      v = as.numeric(gsub(",",".",v))
      output = c(output,v)}
  return(output)
}

BMR_measures4 = read_delim("~/Documents/Paper/Research_Assitant/Bibliography/Bushuev_BRM.csv",delim=";")
change = remove_sd(BMR_measures4$`BMR (ml O2/min)`)
change2 = remove_sd(BMR_measures4$`M (g)`)
BMR_measures4 %>% mutate(`BMR (ml O2/min)`=change, `M (g)` = change2) -> BMR_measures4
summary(lm(log(`BMR (ml O2/min)`/`M (g)` ) ~ log(`M (g)`), data=BMR_measures4)) 


data_soft %>% filter(Latin_name_1 %in% BMR_measures4$Species) %>% arrange(Latin_name) -> d_overl4
BMR_measures4 %>% filter(Species %in% data_soft$Latin_name_1) %>% arrange(Species) -> b_overl4


b_overl4 %>% mutate(Mass_mine=d_overl4$unsexed_mass, Estimation=d_overl4$`Ratio of Median Depth`,latitude = d_overl4$Latitude, tissue= d_overl4$`Tissue type`, Order=d_overl4$Order,Temperature = d_overl4$Temperature) -> b_overl4




#Matching sets

#x = list("Fristoe" = b_overl$Species, "Uyeda" = b_overl2$Other_naming, "Londono" = b_overl3$Species)
#library("VennDiagram")
#venn.diagram(x, "venn_BMR.tiff")


#O2 consumption to watts
conversion = function(x){
  #Jules = x *5/1000* 4186 * 1/36000
  #watts = Jules/20.1
  watts = x *1/3600*5/1000*4186
  return(watts)
}

b_overl %>% mutate(BMR_watts = conversion(`BMR(mlO2/hour)`)) -> b_overl

b_overl %>% arrange(Species) %>% filter(Species %in% b_overl3$Species) -> check
b_overl3 %>% arrange(Species) %>% filter(Species %in% b_overl$Species) -> check2
check %>% mutate(BMR2 = check2$`BMR (watts; = 20.1 J/ml O2)`) -> check
check %>% ggplot(aes(x=BMR2, y=BMR_watts)) + geom_point()

names(b_overl)
names(b_overl3)

b_overl3 %>% mutate("Mass(g)"=`Mean mass (g)`, "BMR(watts)"=`BMR (watts; = 20.1 J/ml O2)`, "Source"="Londono") %>%
  dplyr::select(`Mass(g)`,`BMR(watts)`, Species, Estimation, latitude, Order,Source,Temperature) %>% arrange(Species) -> overlap_londono
b_overl %>%mutate("BMR(watts)"=BMR_watts, "Source"="Fristoe") %>%  dplyr::select(`Mass(g)`,`BMR(watts)`, Species, Estimation, latitude,Order,Source,Temperature)%>% arrange(Species) %>% 
  filter(!Species %in% overlap_londono$Species)-> overlap_Fristoe
b_overl4  %>% mutate("BMR(watts)" = conversion(`BMR (ml O2/min)`)*60) -> b_overl4 
b_overl4 %>% mutate("Source"="Bushuev","Mass(g)"=`M (g)`) %>%  dplyr::select(`Mass(g)`,`BMR(watts)`, Species, Estimation, latitude,Order,Source, Temperature)%>% arrange(Species) %>%
  filter(! Species %in% overlap_londono$Species) %>% filter(! Species %in% overlap_Fristoe$Species) -> overlap_bushuev



rbind(overlap_londono,overlap_Fristoe,overlap_bushuev) -> Overlap

###BMR###
summary(lm(log(`BMR(watts)`/`Mass(g)`) ~ log(`Mass(g)`) + I(Temperature^-1), Overlap))
#####

Overlap %>% ggplot(aes(x=`BMR(watts)`/`Mass(g)`,y=Estimation)) + geom_point() + geom_smooth(method='lm',formula=y~x, se=F)
Overlap %>% filter(Order=="Passeriformes") %>% ggplot(aes(x=`BMR(watts)`/`Mass(g)`,y=Estimation)) + geom_point() + 
  geom_smooth(method='lm',formula=y~x, se=F) + theme_bw()

Overlap %>% mutate(Mass_specific_BMR =`BMR(watts)`/`Mass(g)`) -> Overlap
Overlap %>% filter(Order=="Passeriformes") -> Overlap_pass
summary(lm(Estimation ~ Mass_specific_BMR, Overlap))
summary(lm(Estimation ~ Mass_specific_BMR, Overlap_pass))
cor.test(Overlap$Mass_specific_BMR, Overlap$Estimation, method="pearson")
cor.test(Overlap_pass$Mass_specific_BMR, Overlap_pass$Estimation, method="pearson")


#######################







######Check unsexed's mass correlation with Actual mass##########



replace_usnm = function(INPUT){
  if (grepl("USNM", INPUT)){
    I = gsub("USNM-", "", INPUT)
  }else{
    I = INPUT
  }
  return(I)
}
make_numeric = function(INPUT){
  kilos =F
  if (grepl(";", INPUT)){
    INPUT = strsplit(INPUT, ";")[[1]][1]
  }
  if (grepl("kg", INPUT)){
    kilos = T
    }
  I = as.numeric(gsub("[^0-9.-]", "", INPUT))
  if (kilos == T){ I = I*1000}
  return(I)
}

N = as.vector(sapply(data$`Source specimen ID code`,replace_usnm))
data %>% mutate(`Source specimen ID code` = N) -> data

smith = read_csv("info/nmnhsearch-20191008093347.csv")
smith %>% filter(`Catalog Number - USNM` %in% data$`Source specimen ID code`) %>% arrange(`Catalog Number - USNM`) -> smith

data %>% filter(`Source specimen ID code` %in% smith$`Catalog Number - USNM`) %>% arrange(`Source specimen ID code`) -> data

N = as.vector(sapply(smith$Measurements,make_numeric))
smith %>% mutate(Measurements = N) -> smith

data %>% mutate("Actual_mass" = smith$Measurements) -> data

data %>% ggplot(aes(x=log(unsexed_mass),y=log(Actual_mass))) + geom_point()  + theme_bw() -> sup_mass
ggsave("Supplementary_Mass.pdf",sup_mass)

summary(lm(unsexed_mass ~ Actual_mass, data))




#####Correlation of Estimator with Rates and of BMR with Rates:

data_hyphy = read_tsv(file = "../dNdS/hyphy_merged.tsv")
data_hyphy %>% filter(protein == "TOTAL") %>% filter(`original name` %in% data_soft$ID) %>% arrange(`original name`) ->
  data_hyphy2
data_iqtree = read_tsv("../EvoRates/combined_merged.tsv")
data_iqtree %>% filter(ID_organism %in% data_soft$ID) %>% arrange(ID_organism) -> data_iqtree2


data_soft %>% arrange(ID) %>% filter(ID %in% data_iqtree$ID_organism) %>%
  mutate(IQTREE=data_iqtree2$Rate, dNdS=data_hyphy2$nonsynonymous/data_hyphy2$synonymous,
dN = data_hyphy2$nonsynonymous, dS= data_hyphy2$synonymous) -> data_total

data_total %>% ggplot(aes(x=`Ratio of Median Depth`,y=dN)) + geom_point() + theme_bw()
data_total %>% ggplot(aes(x=`Ratio of Median Depth`,y=dS)) + geom_point() + theme_bw()
data_total %>% ggplot(aes(x=`Ratio of Median Depth`,y=dNdS)) + geom_point() + theme_bw()
data_total %>% ggplot(aes(x=`Ratio of Median Depth`,y=IQTREE)) + geom_point() + theme_bw()


####Correlation with BMR

BMR_measures3 = read_tsv("~/Documents/Paper/Research_Assitant/Bibliography/Lodono_et_al_BMR.tsv")
BMR_measures3$Species

BMR_measures3 %>% filter(Species %in% data_iqtree$Latin_name_1) %>% arrange(Species) %>% filter(!Species == "Gallus gallus") %>% filter(!Species=="Picoides pubescens") -> BMR3_2
data_hyphy %>% filter(protein == "TOTAL") %>% filter(Latin_name_1 %in% BMR3_2$Species) %>% arrange(Latin_name_1) -> data_hyphy
data_iqtree %>% filter(Latin_name_1 %in% BMR3_2$Species) %>% arrange(Latin_name_1) -> data_iqtree

BMR3_2 %>% mutate(IQTREE=data_iqtree$Rate, dNdS=data_hyphy$nonsynonymous/data_hyphy$synonymous,
                  dN = data_hyphy$nonsynonymous, dS= data_hyphy$synonymous) -> final_bmr

final_bmr %>% ggplot(aes(x=`BMR (watts; = 20.1 J/ml O2)`/`Mean mass (g)`,y=dN)) + geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)
final_bmr %>% ggplot(aes(x=`BMR (watts; = 20.1 J/ml O2)`/`Mean mass (g)`,y=dS)) + geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)
final_bmr %>% ggplot(aes(x=`BMR (watts; = 20.1 J/ml O2)`/`Mean mass (g)`,y=dNdS)) + geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)
final_bmr %>% ggplot(aes(x=`BMR (watts; = 20.1 J/ml O2)`/`Mean mass (g)`,y=IQTREE)) + geom_point() + theme_bw() + geom_smooth(method='lm',formula=y~x)


final_bmr %>% mutate(mass_specific_bm= `BMR (watts; = 20.1 J/ml O2)`/`Mean mass (g)`) -> final_bmr
summary(lm(log(mass_specific_bm) ~ log(dN), final_bmr))
summary(lm(log(mass_specific_bm) ~ log(dS), final_bmr))
summary(lm(log(mass_specific_bm) ~ log(dNdS), final_bmr))
summary(lm(log(mass_specific_bm) ~ log(IQTREE), final_bmr))
