#ND3 tree
setwd("~/Documents/Paper/ND3_insertion")
library("treeio")
library(tidyverse)
library(ape)
library(phytools)
library(ggtree)
library(rotl)
library(googlesheets)
library(caper)
library(ggimage)
#ALL_ID = read_csv("../Research_Assitant/Analysis_Mtc_number/info/NAME-ID.txt")
dataset = read_tsv("Position_nucleotide.tsv",col_names=F)

#Need information from which order are they from in order to differentiate Passeriformes from non-passeriformes
#gs_title("Genomes for the family phylogeny") -> sheet
#data_order <- gs_read(ss = sheet, ws = "Genomes_phase2")
#data_order %>% dplyr::select(Latin_name,Order) %>% filter(Order=="Passeriformes") -> data_order
#sub_names = gsub(" ", "_", data_order$Latin_name)
#data_order %>% mutate(Order=sub_names) -> data_order
#####

#dataset %>% mutate(Passeriformes = ifelse( X1 %in% data_order$Order, "yes","no")) -> dataset





#############################################################
#######################Get tree##############################
taxa = gsub(" ","_",dataset$X1)
resolved_names = vector()

resolved_names = tibble(search_string=character(),unique_name=character(),approximate_match=logical(),ott_id=numeric(),is_synonym=logical(),flags=character(),number_matches=numeric())
resolved_names = as.data.frame(resolved_names)
#346 --> Edolisoma coeruslescens ?? Edolisoma mindanense
#555 --> hypsiglena_sp
#609 --> Leptocoma_aspasia ?? Leptocoma sericea
#979 "Rhadina_sibilatrix" ?? Phylloscopus sibillatrix
#1043 Strigops_habroptilus
for (item in  1044:length(taxa)){
  print(item)
  Name = tnrs_match_names(taxa[item])
  resolved_names = rbind(resolved_names, Name)
  }
Name = tnrs_match_names(c("Edolisoma mindanense","Leptocoma sericea","Phylloscopus sibillatrix"))
resolved_names = rbind(resolved_names, Name)

my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)


redo_names = vector()
for (n in my_tree$tip.label){
  new_name = paste(strsplit(n, "_")[[1]][1:2],collapse="_")
  if (new_name == "Edolisoma_mindanense"){ new_name = "Edolisoma_coeruslescens"
  } else if(new_name == "Leptocoma_sericea"){ new_name = "Leptocoma_aspasia"
  } else if(new_name == "Phylloscopus_sibillatrix"){ new_name = "Rhadina_sibilatrix"}
  redo_names =c(redo_names,new_name)
}

my_tree$tip.label = redo_names
#575 tip of Nestor notabilis

bind.tip(tree = my_tree, where = 575,tip.label = "Strigops_habroptilus") -> my_tree
write.tree(my_tree, "Total_tree.tree")

##################################################
Tree = read.tree("Total_tree.tree")

vector1 = vector()
for (i in dataset$X1){
  if(! i %in% Tree$tip.label){
    vector1 = c(vector1,i)
  }
}
vector2 = vector()
for (i in Tree$tip.label){
  if(! i %in% dataset$X1){
    vector2 = c(vector2,i)
  }
}

names = as_tibble(resolved_names)
names$unique_name = gsub(" ","_",names$unique_name)
names %>% filter(unique_name %in% vector2) %>% filter(is_synonym==F) -> toss
drop.tip(Tree,toss$unique_name) -> Tree2

names %>% filter(unique_name %in% vector2) %>% filter(is_synonym==T)  -> switch
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
firstup(switch$search_string) -> switch$search_string
n = match(switch$unique_name,Tree2$tip.label)
Tree2$tip.label[n] = switch$search_string

vector2 = vector()
for (i in Tree2$tip.label){
  if(! i %in% dataset$X1){
    vector2 = c(vector2,i)
  }
}
Tree2$tip.label[Tree2$tip.label=="Edolisoma_coeruslescens"] = "Edolisoma_coerulescens"
Tree3 = drop.tip(Tree2,vector2)

vector1 = vector()
for (i in dataset$X1){
  if(! i %in% Tree3$tip.label){
    vector1 = c(vector1,i)
  }
}
dataset %>% filter(X1 %in% vector1) -> dataset_notused
write_tsv(dataset_notused,"discarted.tsv")
######################################################################
######################################################################






################################################################
################################################################

write.tree(Tree3, "Total_tree2.tree")
write_tsv(names,"naming.tsv")
##################################################################
Tree3 = read.tree("Total_tree2.tree")
dataset_notused = read_tsv("discarted.tsv")

dataset %>% filter(! X1 %in% dataset_notused$X1) -> dataset


#####################Ancestral estate reconstruction
library(castor)

dataset %>% arrange(X1) -> dataset
states = dataset$X2
####Check match
n = match(Tree3$tip.label,dataset$X1)
#####
states = states[n]
states[states=="-"]=1
states[states=="a"]=2
states[states=="t"]=3
states[states=="c"]=4
states[states=="g"]=5
states = as.integer(states)

#reconstruction = hsp_empirical_probabilities(tree = Tree3,tip_states = states) #crude method
reconstruction = hsp_mk_model(tree = Tree3,tip_states = states, rate_model="ER",Nthreads = 3)
#reconstruction = hsp_max_parsimony(tree= Tree3,tip_states = states)

Likelihoods = reconstruction$likelihoods

###Checking weather the leafs matches the given input
Leaf = Likelihoods[0:length(states),]

Check_leaf = function(item){
      if (item[1] == 1){ check_leaf = 1
      }else if (item[2] == 1){check_leaf = 2
      }else if (item[3] == 1){check_leaf = 3
      }else if (item[4] == 1){check_leaf = 4
      }else if (item[5] == 1){check_leaf = 5}
return(check_leaf)
}

MM=apply(Leaf, 1, Check_leaf)
Outcome = MM == states
p =which(Outcome == F)
############
leaf_size = length(Tree3$tip.label)

LL_nodes = Likelihoods[(length(states)+1):dim(Likelihoods)[1],]

Turtle = which(grepl("Testudines",Tree3$node.label))  #Turtles, Testudines
Serpentis = which(grepl("mrcaott1662ott20148", Tree3$node.label)) #Serpentis 
Bifurcata = which(grepl("ott4945781", Tree3$node.label)) #Squamata
Crocodyla = which(grepl("mrcaott31216ott35864", Tree3$node.label)) #crocodylia
Aves = which(grepl("Aves", Tree3$node.label)) #Aves
Passeriformes = which(grepl("Passeriformes", Tree3$node.label)) #passserines

Ancestor_all = which(grepl("mrcaott59ott1662", Tree3$node.label))

Testudines2 = which(grepl("mrcaott59ott3697",Tree3$node.label))
Testudines3 = which(grepl("Pleurodira",Tree3$node.label))
Testudines4 = which(grepl("mrcaott2982ott235762",Tree3$node.label))
Testudines5 = which(grepl("mrcaott3697ott810244",Tree3$node.label))
Testudines6 = which(grepl("mrcaott59ott14158",Tree3$node.label))

data_piechart = LL_nodes[c(Ancestor_all,Turtle,Serpentis,Bifurcata,Crocodyla,Aves,Passeriformes,Testudines2,Testudines3,Testudines4,Testudines4,Testudines6),]
data_piechart = cbind(data_piechart, c(Ancestor_all,Turtle,Serpentis,Bifurcata,Crocodyla,Aves,Passeriformes,Testudines2,Testudines3,Testudines4,Testudines4,Testudines6)+leaf_size)
colnames(data_piechart) = c("-", "A","T","C","G", "node")
data_piechart = as_tibble(data_piechart)



###########



Turtle_ID= "3827b9bf-ca2f-4613-aa9c-28c4101e04d8"
Serp_ID = "b294378c-0b6c-426e-aed4-999b90ad4ee9"
Squamata_ID = "83053aee-0f56-4cf3-bbfa-9207e6f13f46"
Croco_ID= "dffda000-77cb-4251-b837-0cd2ab21ed5b"
Ave_Id= "d5194798-0b50-4a79-8b3d-e6efe468c5ab"
Pass_id ="dfdfb59e-8126-44e1-a7a9-1bf698113e1c"
Phylopic =c(Turtle_ID, Serp_ID,Squamata_ID, Croco_ID,Ave_Id, Pass_id)

test = tibble(node=c(Turtle,Serpentis,Bifurcata,Crocodyla,Aves,Passeriformes) + leaf_size, phylopic = Phylopic)

######################radial with images #####################
PAL =  c("#DD8D29", "#E2D200", "#46ACC8", "#9986A5", "#B40F20")
dataset %>% mutate(Insertion = ifelse(X2 == "-", "Absent", "Present")) -> ds
p = ggtree(Tree3, layout="radial",branch.length="none") 

p %<+% ds  + theme(legend.position="bottom") + geom_tippoint(aes(x=x+1,color=factor(Insertion)),na.rm=TRUE, size= 0.5)+
  labs(color = "ND3-174") +scale_color_manual(values=c(PAL[2],PAL[5])) -> p2
  
p2 + geom_cladelabel(Turtle+leaf_size, "Testudines", offset=7, barsize=2, angle=0, offset.text=1.5, hjust=0.5, fontsize=3) + 
  geom_cladelabel(Serpentis+ leaf_size, "Serpentis", offset=2, barsize=2, angle=-35, offset.text=2, hjust=0.5, fontsize=3) +
  geom_cladelabel(Bifurcata+ leaf_size, "Squamata", offset=7, barsize=2, angle=0, offset.text=4.5, hjust=0.5, fontsize=3) +
  geom_cladelabel(Crocodyla+ leaf_size, "Crocodilia", offset=7, barsize=2, angle=0, offset.text=3.5, hjust=0.5, fontsize=3) +
  geom_cladelabel(Aves+ leaf_size, "Aves", offset=7, barsize=2, angle=0, offset.text=3.5, hjust=0.5, fontsize=3) +
  geom_cladelabel(Passeriformes+ leaf_size, "Passeriformes", offset=2, barsize=2, angle=205, offset.text=2, hjust=0.5, fontsize=3)  ->p3

p3 %<+% test + geom_nodelab(aes(image=phylopic),alpha=1, geom="phylopic", color=PAL[3], size=0.04) -> p4

p4 + labs(colour = "Status Position 174") -> p4

ggsave("Tree.pdf",p4)
ggsave("Tree_paper.pdf",p4)
ggsave("Tree_paper.png",p4)
###############################################################
##########Rectangular with piecharts################
p = ggtree(Tree3, layout="rectangular",branch.length="none") 

p %<+% dataset  + theme(legend.position="bottom") + #xlim(NA, 50) + 
  geom_tippoint(aes(x=x+1,color=factor(X2)),na.rm=TRUE) -> p2

p2 + geom_cladelabel(3+leaf_size, "Testudines", offset=5, barsize=2, angle=0, offset.text=5, hjust=0.5, fontsize=3) + 
  geom_cladelabel(874+ leaf_size, "Serpentis", offset=3, barsize=2, angle=0, offset.text=6, hjust=0.5, fontsize=3) +
  geom_cladelabel(870+ leaf_size, "Squamata", offset=5, barsize=2, angle=0, offset.text=5, hjust=0.5, fontsize=3) +
  geom_cladelabel(861+ leaf_size, "Crocodylia", offset=5, barsize=2, angle=0, offset.text=5, hjust=0.5, fontsize=3) +
  geom_cladelabel(87+ leaf_size, "Aves", offset=5, barsize=2, angle=0, offset.text=5, hjust=0.5, fontsize=3) +
  geom_cladelabel(97+ leaf_size, "Passeriformes", offset=3, barsize=2, angle=0, offset.text=6, hjust=0.5, fontsize=3)  ->p3
pies = nodepie(data_piechart, cols=1:5)
ggtree::inset(p3, pies,width=1, height=1) 
ggsave("Tree_pie.pdf",p4)
##########################################################



dataset %>% group_by(X2) %>% summarise(Counts= n())


getDescendants(Tree3,node=Passeriformes+leaf_size) -> Descendant_passerines
leafs_pass = Descendant_passerines[Descendant_passerines <= leaf_size]
Name_pass = Tree3$tip.label[leafs_pass]

dataset %>% mutate(is.pass = ifelse(X1 %in% Name_pass, T,F)) -> dataset

dataset %>% filter(is.pass == T) %>% group_by(X2) %>% summarise(n())
dataset %>% filter(is.pass == T) %>% filter(!X2 == "-")





extract.clade(Tree3, node=(Turtle + leaf_size)) -> All_testudinest
extract.clade(Tree3,node=Aves+leaf_size) -> All_birdst
extract.clade(Tree3, node=Crocodyla+leaf_size) -> All_crocodylat
extract.clade(Tree3, node=Passeriformes+leaf_size) -> All_Passeriformest
extract.clade(Tree3, node= leaf_size+Bifurcata) -> All_lizardst
extract.clade(Tree3, node= leaf_size+Serpentis) -> All_snakest

dataset %>% filter(X1 %in% All_testudinest$tip.label) -> Data_testudines
dataset %>% filter(X1 %in% All_birdst$tip.label) -> Data_birds
dataset %>% filter(X1 %in% All_Passeriformest$tip.label) -> Data_Passeriformes
dataset %>% filter(X1 %in% All_lizardst$tip.label) -> Data_lizards
dataset %>% filter(X1 %in% All_snakest$tip.label) -> Data_snakes
dataset %>% filter(X1 %in% All_crocodylat$tip.label) -> Data_crocodile



Data_testudines %>% mutate(Order="Testudines", Gap = ifelse(X2=="-", "yes","no")) -> Com1
Data_birds %>% mutate(Order="Birds", Gap = ifelse(X2=="-", "yes","no")) -> Com2
Data_lizards %>% mutate(Order="Lizards", Gap = ifelse(X2=="-", "yes","no")) -> Com3
Data_crocodile %>% mutate(Order="Crocodile", Gap = ifelse(X2=="-", "yes","no")) -> Com4


rbind(Com1,Com2,Com3,Com4) -> complete_meta
write.csv(x = complete_meta, file = "Metadata_sp.csv")


dim(Data_testudines)
dim(Data_birds)
dim(Data_Passeriformes)
dim(Data_lizards)
dim(Data_snakes)
dim(Data_crocodile)
103+807+202+20

dataset %>% group_by(X2) %>% summarise(n())
Data_testudines %>% group_by(X2) %>% summarise(n())
Data_birds %>% group_by(X2) %>% summarise(n())
Data_birds %>% filter(!X1 %in% Data_Passeriformes$X1) %>% group_by(X2) %>% summarise(n())
Data_Passeriformes %>% group_by(X2) %>% summarise(n())
Data_lizards%>% group_by(X2) %>% summarise(n())
Data_snakes%>% group_by(X2) %>% summarise(n())
Data_crocodile%>% group_by(X2) %>% summarise(n())
Data_birds %>% filter(!X1 %in% Data_Passeriformes$X1) -> Data_nopass

list_subsets = list(Data_testudines, Data_birds, Data_nopass, Data_Passeriformes)
vector_nodes = c(Turtle, NA, NA, NA)

library(phangorn)
for (i in seq(1,length(list_subsets))){
  Subset = list_subsets[[i]]
  
  
  keep.tip(Tree3, Subset$X1) -> Tr
  p = ggtree(Tr, layout="radial",branch.length="none") 
  p %<+% Subset  + theme(legend.position="bottom") +
    geom_tippoint(aes(x=x+1,color=factor(X2)),na.rm=TRUE) -> p2
  print(p2)
  
  if (is.na(vector_nodes[i])){
    next
  }
  Desc = Descendants(Tree3, vector_nodes[i]+leaf_size, type = "all")
  Desc = Desc[Desc>leaf_size]
  Desc = c(Desc,vector_nodes[i]+leaf_size)
  Desc = sort(Desc)
  
  
 
  
  if (i==2){
    as.data.frame(Likelihoods) %>% rownames_to_column() %>% as_tibble -> LL_nodes
    LL_nodes = LL_nodes[Desc,]
    LL_nodes %>% mutate(rowname = as.numeric(rowname)) %>% mutate(rowname_s=rowname - min(rowname) + 1 +length(Tr$tip.label)) -> LL_nodes
    Passeriformes_list = Descendants(Tree3, Passeriformes+leaf_size, type = "all")
    LL_nodes %>% filter(! as.numeric(rowname) %in% Passeriformes_list) -> LL_nodes
    LL_nodes %>% filter(V4 < 0.6) -> LL_nodes
    data_piechart = LL_nodes[,c(2,3,4,5,6,7)]
    
  } else{
    LL_nodes = Likelihoods[Desc,]
    LL_nodes = as.data.frame(LL_nodes) 
    data_piechart = cbind(LL_nodes, seq(dim(Subset)[1]+1,dim(Subset)[1]+Tr$Nnode))
  }
  
  
  colnames(data_piechart) = c("-", "A","T","C","G", "node")
  data_piechart = as_tibble(data_piechart)
  
  p = ggtree(Tr, layout="rectangular",branch.length="none") 
  p %<+% Subset  + theme(legend.position="bottom") + labs(color = "ND3-174") +
    geom_tippoint(aes(x=x+1,color=factor(X2)),na.rm=TRUE) + 
    scale_color_manual(values=c(PAL[2], PAL[1], PAL[3], PAL[4], PAL[5])) -> p2
  pies = nodepie(data_piechart, cols=1:5, color =c(PAL[2], PAL[1], PAL[3], PAL[4], PAL[5]))
  ggtree::inset(p2, pies,width=1, height=1) -> p3
  ggsave(plot = p3,filename = "Supplmentary_Turtles.pdf")
  ggsave(plot = p3, filename="Supplmentary_Turtles.png")
}




##Bird orders
orders_info = read_tsv("All_bird_orders.tsv", col_names=F)
orders_info %>% arrange(X1) %>% filter(X1 %in% Data_birds$X1) -> orders_info
Data_birds %>% mutate(Family = orders_info$X2, Order= orders_info$X3) -> Data_birds

Data_birds %>% group_by(Order, X2) %>% summarise(Number = n()) %>% spread(X2, Number) -> Information_birds_orders
Information_birds_orders %>% mutate(c = ifelse(is.na(c), 0, c),g = ifelse(is.na(g), 0, g),t = ifelse(is.na(t), 0, t), `-`=ifelse(is.na(`-`), 0, `-`) ) %>%
  mutate(Total = `-`+c + g + t, Frq_gap = ifelse(is.na(`-`),0, `-`/Total)) -> Order_insertion
Order_insertion %>% ggplot(aes(x=Order,y=Frq_gap)) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  geom_text(aes(y= -0.01, label=Total)) -> Freq_insertion_by_order

ggsave(plot = Freq_insertion_by_order, filename = "Birds_Insertion_Order.pdf")
ggsave(plot = Freq_insertion_by_order, filename = "Birds_Insertion_Order.png")
#Kernel "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine"




Jarvis = read.tree("/Users/Sergio/Documents/Paper/Research_Assitant/Chronogram01.TENT.ExAML.tre")
INFO = read_delim("/Users/Sergio/Documents/Paper/Research_Assitant/Analysis_Mtc_number/info/Trait_data2.csv", delim = ";")
INFO$Latin_name

new_lab = vector()
n = 0
for (label in Jarvis$tip.label){
  l = strsplit(label,"_")[[1]][1:2]
  nl = paste(l, collapse = "_")
  if (nl == "Tinamus_major"){
    nl = "Tinamus_guttatus"
  } else if(nl == "Manacus_vitellinus"){
    nl = "Manacus_manacus"
  }
  INFO %>% filter(Latin_name == nl) -> S
  if (S$order %in% new_lab){
    S = paste(c(S$order,as.character(n)), collapse="")
  }else{ S = S$order}
  n = n+1
  new_lab = c(new_lab,S)
  
}

Jarvis$tip.label = new_lab
Jarvis$tip.label[Jarvis$tip.label=="Sphenisciformes24"] = "Apterygiformes"
drop.tip(Jarvis, c("Caprimulgiformes14","Caprimulgiformes15","Passeriformes35","Passeriformes36","Passeriformes37","Passeriformes34")) -> p

match(Order_insertion$Order, Jarvis$tip.label)
Order_insertion %>% filter(Order %in% Jarvis$tip.label) -> Order_insertion
keep.tip(Jarvis, Order_insertion$Order) -> Jarvis
Order_insertion[match(Jarvis$tip.label,Order_insertion$Order),] -> Order_insertion



library(ggstance)

ggtree(Jarvis, layout = "rectangular")  + theme_tree2() + geom_tiplab(offset = -20,vjust = 1 )  -> p  #+ labs(caption = "Divergence time")
revts(p) + theme(plot.margin=margin(60,60,60,60)) -> p #+ xlim(-105,14) -> p
facet_plot(p, panel="Frequency Insertion", data=Order_insertion, geom= geom_barh, mapping=aes(x=Order_insertion$Frq_gap),stat="identity") -> p2
facet_plot(p2,panel="Frequency Insertion", data=Order_insertion, geom=geom_text, mapping=aes(x= -0.05, label=Order_insertion$Total) ) -> p3

ggsave(plot = p3, filename = "TreeBirds_Insertion_Order.pdf")
ggsave(plot = p3, filename = "TreeBirds_Insertion_Order.png")



##Entrophy

data_entrophy = read_tsv("Shanon.tsv")
data_entrophy %>% mutate(Information = log2(4) - Entropy, Information_codon = log2(4) - codon_bin) -> data_entrophy
data_entrophy %>% filter(Position == 174)

BIN = 1
data_entrophy %>% mutate(Gap_status = ifelse(Gap_status=="yes", "Gap", "Insertion")) -> data_entrophy



######################
#data_entrophy %>% ggplot(aes(x=Position, y=Entropy)) + geom_bar(stat="identity") + facet_grid(Group ~ Gap_status)

#data_entrophy %>% ggplot(aes(x=Position, y=codon_bin)) + geom_bar(stat="identity", width=BIN, aes(fill=Information_codon)) + facet_grid(Group ~ Gap_status) + 
#  geom_vline(xintercept= 174 + (-BIN+1), linetype="dashed",  color = "blue", size=0.1) +  coord_flip() + theme_bw()

#data_entrophy %>% ggplot(aes(x=Position, y=codon_bin)) + geom_bar(stat="identity", width=BIN) + facet_grid(Group ~ Gap_status) + geom_bar(aes(y=-0.2,x=Position,fill= Information^3), stat="identity") +
#  geom_vline(xintercept= 174 + (-BIN+1), linetype="dashed",  color = "blue", size=0.1) +  coord_flip() + theme_bw() +
#viridis::scale_fill_viridis()

#data_entrophy %>% ggplot(aes(x=Position, y=codon_bin)) + geom_bar(stat="identity", width=BIN) + facet_grid(Group ~ Gap_status) + geom_bar(aes(y=-0.2,x=Position,fill= Information_codon), stat="identity") +
#  geom_vline(xintercept= 174 + (-BIN+1), linetype="dashed",  color = "blue", size=0.1) +  coord_flip() + theme_bw() +
#  viridis::scale_fill_viridis()

#data_entrophy %>% ggplot(aes(x=Gap_status, y=1, fill=Information^3)) + geom_bar(stat="identity") + facet_grid(~Group) +
#  theme_bw() +geom_hline(yintercept= 174, linetype="dashed",  color = "red", size=0.1) + 
#  geom_hline(yintercept= 164, linetype="dashed",  color = "black", size=0.5) + 
#  geom_hline(yintercept= 180, linetype="dashed",  color = "black", size=0.5)+
#  viridis::scale_fill_viridis() + coord_flip()
#########################

Legend = c(NULL)
for (n in seq(1,dim(data_entrophy)[1])){
  entry = data_entrophy[n,]
  L = paste(c(entry$Gap_status,entry$Group), collapse="_")
  Legend = c(Legend,L)
}

Legend
data_entrophy %>% mutate(Name = Legend) -> data_entrophy
data_entrophy %>% ggplot(aes(x=Name, y=1, fill=Information^3)) + geom_bar(stat="identity")  +
  theme_bw() +geom_hline(yintercept= 174, linetype="dashed",  color = "red", size=0.1) + 
  geom_hline(yintercept= 164, linetype="dashed",  color = "black", size=0.5) + 
  geom_hline(yintercept= 180, linetype="dashed",  color = "black", size=0.5)+ coord_flip()+
  labs(fill = "Information content (to the cube)") + ylab(label = "ND3 position") -> Information_content
ggsave(plot = Information_content,filename = "Information.pdf")
ggsave(plot = Information_content,filename = "Information.png")
