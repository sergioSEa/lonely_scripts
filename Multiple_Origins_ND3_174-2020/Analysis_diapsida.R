#Script used in order to analyze Diapsida records.

#The output of Prepare_dataset.sh generates a table of Species and nucleotide found in position 174 and a MSA of those records
#In this  script we will only use the records that correspond to the phylogenetic group of Diapsida

#Parts:
#1. Ancestral State reconstruction: Requires dataset of diapsida (subset of Prepare_dataset.sh) and phylogenetic tree (can be obtained from Rotl, check Make_otol_tree.R)
#2. Transition count: Requires state likelihoods produced at #1
#3. Study of codon frequency by bird order. Requires states in birds (can be obtained at #1), a table of Orders per bird and, for visualization, the topology of ED Jarvis 2015 Science
#4. Entrophy. Requires entrophy metrics. Requires to run Calculate_diversity.py, which input is both the output from #1 and a MSA only of diapsida.
#5. Codon analysis: Requieres output from study_codons_insertion2.py (input generated from Divide.py)



library(tidyverse)
library(ape)
library(phytools)
library(ggtree)
library(castor)
library(ggimage)
library(phangorn)
library(ggstance)
library(Biostrings)
library(ggrepel)


Ancestral_state = function(Path_data = "Diapsida_table.tsv", Path_tree="Total_tree.tree"){
  #Dataset needed. Table of states and phylogenetic tree
  dataset = read_tsv(Path_data,col_names=F)
  Tree3 = ape::read.tree(Path_tree)
  
  print("Filtering out species from dataset not found in the tree")
  #If any UNVERIFIED record present, remove
  dataset %>% filter(! grepl("UNVERIFIED",X1)) -> dataset
  #Homogenize records
  print(paste(c("Filtered records:",dim(filter(dataset, ! X1 %in% Tree3$tip.label))[1]), collapse= " "))
  dataset %>% filter(X1 %in% Tree3$tip.label) -> dataset
  keep.tip(Tree3, dataset$X1) -> Tree3
  
  
  ##############################################
  #######Ancestral state reconstruction#########
  ##############################################
  
  print("Filtering ambiguous status at ND3-174")
  print(paste(c("Filtered records:",dim(filter(dataset, (X2 == "n" | is.na(X2))))[1]),collapse=" "))
  
  dataset %>% arrange(X1) %>% filter(! X2 == "n") %>% filter(! is.na(X2)) -> dataset
  keep.tip(Tree3, dataset$X1) -> Tree3
  states = dataset$X2
  
  n = match(Tree3$tip.label,dataset$X1)
  
  #####Change encoding for package
  print("Fitting a 5 state model with a ML equal-rate params")
  states = states[n]
  states[states=="-"]=1
  states[states=="a"]=2
  states[states=="t"]=3
  states[states=="c"]=4
  states[states=="g"]=5
  states = as.integer(states)
  Tree3$edge.length = NULL #Otol tree has 0s in the branch lengths. The model fails if that is the case.
  #Rate model is chosen to equal rates. 
  reconstruction = hsp_mk_model(tree = Tree3,tip_states = states, rate_model="ER",Nthreads = 3)
  #All likelihoods
  Likelihoods = reconstruction$likelihoods
  
  #Looking for Likelihoods only on nodes, positions after the leafs
  leaf_size = length(Tree3$tip.label)
  LL_nodes = Likelihoods[(length(states)+1):dim(Likelihoods)[1],]
  
  
  #Check different nodes
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
  Birds_and_Croco = which(grepl("mrcaott246ott31216", Tree3$node.label))
  
  print("Assessed Likelihood at each clade ancestor")
  colnames(Likelihoods) = c("-","a","t","c","g")
  Likelihoods[c(Ancestor_all, Bifurcata, Serpentis, Turtle, Crocodyla, Aves, Passeriformes, Birds_and_Croco)+leaf_size,] -> Ancestor_info
  rownames(Ancestor_info) = c("Ancestor_all", "Bifurcata", "Serpentis", "Turtle", "Crocodyla", "Aves", "Passeriformes", "Birds-Croco")
  round(Ancestor_info,2)
  
  #Use the LL of interest for making piecharts
  data_piechart = LL_nodes[c(Ancestor_all,Turtle,Serpentis,Bifurcata,Crocodyla,Aves,Passeriformes,Testudines2,Testudines3,Testudines4,Testudines4,Testudines6),]
  data_piechart = cbind(data_piechart, c(Ancestor_all,Turtle,Serpentis,Bifurcata,Crocodyla,Aves,Passeriformes,Testudines2,Testudines3,Testudines4,Testudines4,Testudines6)+leaf_size)
  colnames(data_piechart) = c("-", "A","T","C","G", "node")
  data_piechart = as_tibble(data_piechart)
  #Get images for the tree
  Turtle_ID= "3827b9bf-ca2f-4613-aa9c-28c4101e04d8"
  Serp_ID = "b294378c-0b6c-426e-aed4-999b90ad4ee9"
  Squamata_ID = "83053aee-0f56-4cf3-bbfa-9207e6f13f46"
  Croco_ID= "dffda000-77cb-4251-b837-0cd2ab21ed5b"
  Ave_Id= "d5194798-0b50-4a79-8b3d-e6efe468c5ab"
  Pass_id ="dfdfb59e-8126-44e1-a7a9-1bf698113e1c"
  Phylopic =c(Turtle_ID, Serp_ID,Squamata_ID, Croco_ID,Ave_Id, Pass_id)
  test = tibble(node=c(Turtle,Serpentis,Bifurcata,Crocodyla,Aves,Passeriformes) + leaf_size, phylopic = Phylopic)
  #Make tree plot
  print("Making plot of diapsida tree")
  
  PAL =  c("#DD8D29", "#E2D200", "#46ACC8", "#9986A5", "#B40F20")
  dataset %>% mutate(Insertion = ifelse(X2 == "-", "Absent", "Present")) -> ds
  p = ggtree(Tree3, layout="fan",branch.length="none") 
  p %<+% ds  + theme(legend.position="bottom") + aes(color=factor(Insertion)) + labs(color = "ND3-174") + geom_tiplab(aes(label=label, angle=angle), size= 0.7) -> p2
  p2 + geom_cladelabel(Turtle+leaf_size, "Testudines", offset=13, barsize=1, angle=-50, offset.text=2.5, hjust=0.5, fontsize=2) +  geom_cladelabel(Serpentis+ leaf_size, "Serpentis", offset=7, barsize=0.5, angle=-65, offset.text=2, hjust=0.5, fontsize=2) + geom_cladelabel(Bifurcata+ leaf_size, "Squamata", offset=13, barsize=1, angle=-78, offset.text=4.5, hjust=0.5, fontsize=2) + geom_cladelabel(Crocodyla+ leaf_size, "Crocodilia", offset=13, barsize=1, angle=-45, offset.text=6, hjust=0.5, fontsize=2) + geom_cladelabel(Aves+ leaf_size, "Aves", offset=13, barsize=1, angle=113, offset.text=3.5, hjust=0.5, fontsize=2) + geom_cladelabel(Passeriformes+ leaf_size, "Passeriformes", offset=7, barsize=0.5, angle=160, offset.text=2, hjust=0.5, fontsize=2)  ->p3
  p3 %<+% test + geom_nodelab(aes(image=phylopic),alpha=1, geom="phylopic", color=PAL[3], size=0.04) -> p4
  p4 + labs(colour = "Status Position 174") -> p4
  
  print("Printing labeled tree")
  print(p4)
  ggsave("FIGURES/Labeled_tree.pdf", p4, scale=4)
  ggsave("FIGURES/Labeled_tree.png", p4, scale=4)
  
  p = ggtree(Tree3, layout="fan",branch.length="none") 
  p %<+% ds  +  theme(legend.position="bottom") + aes(color=factor(Insertion)) +
    labs(color = "ND3-174") -> p2
  p2 + geom_cladelabel(Turtle+leaf_size, "Testudines", offset=13, barsize=1, angle=-50, offset.text=2.5, hjust=0.5, fontsize=2) +  geom_cladelabel(Serpentis+ leaf_size, "Serpentis", offset=7, barsize=0.5, angle=-65, offset.text=2, hjust=0.5, fontsize=2) + geom_cladelabel(Bifurcata+ leaf_size, "Squamata", offset=13, barsize=1, angle=-78, offset.text=4.5, hjust=0.5, fontsize=2) + geom_cladelabel(Crocodyla+ leaf_size, "Crocodilia", offset=13, barsize=1, angle=-45, offset.text=6, hjust=0.5, fontsize=2) + geom_cladelabel(Aves+ leaf_size, "Aves", offset=13, barsize=1, angle=113, offset.text=3.5, hjust=0.5, fontsize=2) + geom_cladelabel(Passeriformes+ leaf_size, "Passeriformes", offset=7, barsize=0.5, angle=160, offset.text=2, hjust=0.5, fontsize=2)  ->p3
  p3 %<+% test + geom_nodelab(aes(image=phylopic),alpha=1, geom="phylopic", color=PAL[3], size=0.04) -> p4_2
  p4_2 + labs(colour = "Status Position 174") -> p4_2
  
  print("Printing unlabeled tree")
  print(p4_2)
  ggsave("FIGURES/UnLabeled_tree.pdf", p4_2, scale=4)
  ggsave("FIGURES/UnLabeled_tree.png", p4_2, scale=4)
  print("Splitting entries per taxonomy")
  #Label taxonomy for all Speecies
  extract.clade(Tree3, node=(Turtle + leaf_size)) -> All_testudinest
  extract.clade(Tree3,node=Aves+leaf_size) -> All_birdst
  extract.clade(Tree3, node=Crocodyla+leaf_size) -> All_crocodylat
  extract.clade(Tree3, node=Passeriformes+leaf_size) -> All_Passeriformest
  extract.clade(Tree3, node= leaf_size+Bifurcata) -> All_lizards
  extract.clade(Tree3, node= leaf_size+Serpentis) -> All_snakest
  dataset %>% filter(X1 %in% All_testudinest$tip.label) -> Data_testudines
  dataset %>% filter(X1 %in% All_birdst$tip.label) -> Data_birds
  dataset %>% filter(X1 %in% All_Passeriformest$tip.label) -> Data_Passeriformes
  dataset %>% filter(X1 %in% All_lizards$tip.label) -> Data_lizards
  dataset %>% filter(X1 %in% All_snakest$tip.label) -> Data_snakes
  dataset %>% filter(X1 %in% All_crocodylat$tip.label) -> Data_crocodile
  
  
  Data_testudines %>% mutate(Order="Testudines", Gap = ifelse(X2=="-", "yes","no")) -> Com1
  Data_birds %>% mutate(Order="Birds", Gap = ifelse(X2=="-", "yes","no")) -> Com2
  Data_lizards %>% mutate(Order="Lizards", Gap = ifelse(X2=="-", "yes","no")) -> Com3
  Data_crocodile %>% mutate(Order="Crocodile", Gap = ifelse(X2=="-", "yes","no")) -> Com4
  print("Saving metadata with group information at 'Metadata_sp.csv', this information can be use for running the Shanon script")
  rbind(Com1,Com2,Com3,Com4) -> complete_meta
  write.csv(x = complete_meta, file = "Metadata_sp.csv") 
  
  
  ############################################################################
  #########Find all cases where transitions happen in birds and turtles#######
  ###########################################################################
  print("Finding transitions between states in Birds and turtles")
  find_transitions = function(Parent_node,threshold, STATES=c("-", "A","T","C","G")){
    #Funtion for finding transition cases. From a parental node, e.g the bird ancestor, takes al descendants and go one by one checking if their most recent ancestor is predicted with a likelihood of at least threshold to have the insertion.
    Transition =vector()
    #Get all descendants from a given node
    Desc = sort(c(Descendants(Tree3, Parent_node, type = "all"),Parent_node))
    #Get likelihoods for all those descendents
    LLH = Likelihoods[Desc,]
    as_tibble(LLH) %>% mutate(N=Desc) -> LLH
    #Go descendant by descendant and check if current state did change from its closest ancestor
    for (n in Desc){
      Node = Ancestors(Tree3,n)
      if (length(Node) == 0){ next}
      Node = Node[1]
      Prob = filter(LLH,N==n)
      Prob_ancest = filter(LLH,N==Node)
      if (dim(Prob_ancest)[1] == 0){ next}
      P1 =0
      P2 =0
      for (m in seq(1:5)){
        P1_a = Prob[m]
        P2_a = Prob_ancest[m]
        if (P1_a > P1){
          P1 = P1_a
          state = STATES[m]
        }
        if (P2_a > P2){
          P2 = P2_a
          state_a = STATES[m]
        }
      }
      if (P1 < threshold){ next }
      if (P2 < threshold){ next }
      if(! state_a == state ){
        if(n > length(Tree3$tip.label)){
          Name = Tree3$node.label[n-length(Tree3$tip.label)]
        }else{ Name = Tree3$tip.label[n]}
        Transition = c(Transition, paste(c(state_a,state),collapse=""))
        if(n<length(Tree3$tip.label)+1){
          print(paste(c("Transition",state_a, "to", state, "in node",n, "named", Name),collapse=" "))}
      }
      
    }
    return(Transition)
  }
  #Fint the transitions for turtles and for birds in two different minimal threshold of marginal likelihood.
  print("Tranisitions in turtle above 0.6 likelihood")
  t1 = find_transitions(Parent_node = Turtle + leaf_size, threshold = 0.6)
  print("Tranisitions in turtle above 0.9 likelihood")
  t12 = find_transitions(Parent_node = Turtle + leaf_size, threshold = 0.9)
  print("Tranisitions in Birds above 0.6 likelihood")
  t2 = find_transitions(Parent_node = Aves + leaf_size, threshold = 0.6)
  print("Tranisitions in Birds above 0.9 likelihood")
  t22 = find_transitions(Parent_node = Aves + leaf_size, threshold = 0.9)
  
  Transitions_df = tibble(Transition = t1, species = "Turtle", threshold = 0.6)
  Transitions_df2 = tibble(Transition = t12, species = "Turtle", threshold = 0.9)
  Transitions_df3 = tibble(Transition = t2, species = "Bird", threshold = 0.6)
  Transitions_df4 = tibble(Transition = t22, species = "Bird", threshold = 0.9)
  
  Transitions_df_total = rbind(rbind(rbind(Transitions_df, Transitions_df2),Transitions_df3), Transitions_df4)
  print("Summarizing number of transitions per clade")
  Transitions_df_total %>% group_by(Transition,species, threshold) %>% summarise(n()) -> SAVE_transitions
  Transitions_df_total %>% group_by(Transition,species, threshold) %>% summarise(n()) %>% print()
  
  print("Producing testudines and aves phylogenetic tree")
  
  
 ############################
 #########TESTDUINES TREE####
 ############################
  Subset = Data_testudines
  keep.tip(Tree3, Subset$X1) -> Tr
  
  Desc = Descendants(Tree3, Turtle+leaf_size, type = "all")
  Turtle_leafs = Desc[Desc<leaf_size]
  Desc = Desc[Desc>leaf_size]
  Desc = c(Desc,Turtle+leaf_size)
  Desc = sort(Desc)
  #Desc = c(Desc, Turtle_leafs)
  
  LL_nodes = Likelihoods[Desc,]
  LL_nodes = as.data.frame(LL_nodes) 
  data_piechart = cbind(LL_nodes, seq(dim(Subset)[1]+1,dim(Subset)[1]+Tr$Nnode))
  
  STATES = c("-", "A","T","C","G")
  colnames(data_piechart) =c(STATES, "node")
  data_piechart = as_tibble(data_piechart)
  
  print("Unlabelled turtle tree")
  p = ggtree(Tr, layout="rectangular",branch.length="none") 
  p %<+% Subset  + theme(legend.position="bottom") + labs(color = "ND3-174") +
    geom_tippoint(aes(x=x+1,color=factor(X2)),size=0.3,na.rm=TRUE) + 
    scale_color_manual(values=c(PAL[2], PAL[1], PAL[3], PAL[4], PAL[5])) -> p2
  pies = nodepie(data_piechart, cols=1:5, color =c(PAL[2], PAL[1], PAL[3], PAL[4], PAL[5]))
  ggtree::inset(p2, pies,width=0.06, height=0.06) -> p3_turtle ###THE ERROR LINE
  print(p3_turtle)
  ggsave(plot = p3_turtle, filename  = "FIGURES/Testudines_tree.png", scale = 2)
  ggsave(plot = p3_turtle, filename  = "FIGURES/Testudines_tree.pdf", scale = 2)
  
  
  print("Labelled turtle tree")
  p = ggtree(Tr, layout="rectangular",branch.length="none") 
  p %<+% Subset  + theme(legend.position="bottom") + labs(color = "ND3-174") +
    geom_tiplab(aes(color = factor(X2)), geom = "text", size=0.75) + 
    scale_color_manual(values=c(PAL[2], PAL[1], PAL[3], PAL[4], PAL[5])) -> p2
  pies = nodepie(data_piechart, cols=1:5, color =c(PAL[2], PAL[1], PAL[3], PAL[4], PAL[5]))
  ggtree::inset(p2, pies,width=.06, height=.06) -> p3_turtle
  
  ggsave(plot = p3_turtle, filename  = "FIGURES/Testudines_tree_labeled.png", scale= 2)
  ggsave(plot = p3_turtle, filename  = "FIGURES/Testudines_tree_labeled.pdf", scale = 2)
  
  #Prepare tips
  #select(unique(Subset), X2) %>% as_vector() -> x
  #names(x) = unique(Subset)$X1
  #Prepare nodes
  #y = as.matrix(data_piechart[,1:5])
  #names(y) = data_piechart[6]
  
  #Tr
  #setNames(c("red","yellow","green", "blue","purple"),c("-","a", "t", "c","g")) -> PAL_cutre
  #dotTree(Tr,as.matrix(x)[Tr$tip.label,], legend = F, fsize=0.7,ftype="i",use.edge.length = FALSE, colors=PAL_cutre)
  #nodelabels(pie=y,adj=c(0.5,0.9),cex=0.5)
  
  
    
  print("Finishing function. Returning tibble of birds and tibble of passeriformes")
  return(list(Data_birds, Data_Passeriformes, SAVE_transitions))  
}

Ancestral_state_parsimony = function(Path_data = "Diapsida_table.tsv",Path_tree="Total_tree.tree"){
  dataset = read_tsv(Path_data,col_names=F)
  Tree3 = phyloseq::read_tree("Total_tree.tree")
  
  #If any UNVERIFIED record present, remove
  dataset %>% filter(! grepl("UNVERIFIED",X1)) -> dataset
  #Homogenize records
  dataset %>% filter(X1 %in% Tree3$tip.label) -> dataset
  keep.tip(Tree3, dataset$X1) -> Tree3
  
  #Ancestral reconstruction
  
  dataset %>% arrange(X1) %>% filter(! X2 == "n") %>% filter(! is.na(X2)) -> dataset
  keep.tip(Tree3, dataset$X1) -> Tree3
  states = dataset$X2
  
  n = match(Tree3$tip.label,dataset$X1)
  
  #####Change encoding for package
  states = states[n]
  states[states=="-"]=1
  states[!states=="1"]=2
  
  states = as.integer(states)
  Tree3$edge.length = NULL #Otol tree has 0s in the branch lengths. The model fails if that is the case.
  
  print("Fitting maximum parsimony model of ND3-174 ancestral state")
  #reconstruction = castor::hsp_max_parsimony(tree = Tree3,tip_states = states)
  #Phangorn
  #phyDat(as.matrix(states2), type ="USER", levels=c("2", "1") ) -> CHE
  #phangorn::ancestral.pars(Tree3,CHE,type="MPR", return="prob") -> answer
  #Ape
  Out =  Tree3$tip.label[length(Tree3$tip.label)]
  MPR(phy=unroot(multi2di(Tree3)), x=states, outgroup = Out) -> Predictions
  as.data.frame(Predictions) %>% rownames_to_column("Node")    %>% as_tibble() %>% filter(Node %in% Tree3$node.label) %>%
    mutate(V1 = ifelse(lower == 1,1,0), V2 = ifelse(upper == 1,0,1), V1 = ifelse(upper!=lower, 0.5, V1 ),  V2= ifelse(upper!=lower, 0.5,V2)) -> Predictions
  
  Predictions[match(Tree3$node.label, Predictions$Node),] %>% select(V1,V2) -> Predictions
  tibble(V1 = states == 1, V2 = states== 2) %>% mutate(V1 = ifelse(V1==1, 1, 0), V2 = ifelse(V2==2, 1, 0)) -> tip_status
  Predictions[1,] = c(1,0)
  rbind(tip_status, Predictions) -> Likelihoods

  print("Fitting maximum likelihood model of ND3-174 ancestral state for two states")
  reconstruction2 = castor::hsp_mk_model(tree = Tree3,tip_states = states, rate_model="ER",Nthreads = 3)
  
  #All likelihoods
  #Likelihoods = reconstruction$likelihoods
  Likelihoods_mk = reconstruction2$likelihoods
  #Looking for Likelihoods only on nodes, positions after the leafs
  leaf_size = length(Tree3$tip.label)
  LL_nodes = Likelihoods[(length(states)+1):dim(Likelihoods)[1],]
  LL_nodes_mk = Likelihoods_mk[(length(states)+1):dim(Likelihoods_mk)[1],]
  
  #Cases were the 0.9 threshold does not agree with parsimony 30/3091
  print("Comparison with 2-state ML model")
  as_tibble(LL_nodes_mk) %>% mutate(V1 = ifelse(V1 > 0.9, 1, V1), V1 = ifelse(V1 < 0.1, 0, V1),V2 = ifelse(V2 > 0.9, 1, V2), V2 = ifelse(V2 < 0.1, 0, V2) ) -> hsp_check
  as_tibble(LL_nodes) -> parsimony_check
  sum((hsp_check$V1 == parsimony_check$V1) == F)
  
  
  #Check different nodes
  Turtle = which(grepl("Testudines",Tree3$node.label))  #Turtles, Testudines
  Serpentis = which(grepl("mrcaott1662ott20148", Tree3$node.label)) #Serpentis 
  Bifurcata = which(grepl("ott4945781", Tree3$node.label)) #Squamata
  Crocodyla = which(grepl("mrcaott31216ott35864", Tree3$node.label)) #crocodylia
  Aves = which(grepl("Aves", Tree3$node.label)) #Aves
  Passeriformes = which(grepl("Passeriformes", Tree3$node.label)) #passserines
  Ancestor_all = which(grepl("mrcaott59ott1662", Tree3$node.label))
  
  print("Assessed Likelihood at each clade ancestor")
  colnames(Likelihoods) = c("-","insertion")
  Likelihoods[c(Ancestor_all, Bifurcata, Serpentis, Turtle, Crocodyla, Aves, Passeriformes)+leaf_size,] -> Ancestor_info
  rownames(Ancestor_info) = c("Ancestor_all", "Bifurcata", "Serpentis", "Turtle", "Crocodyla", "Aves", "Passeriformes")
  round(Ancestor_info,2)
  
  
  
  ############################################################################
  #########Find all cases where transitions happen in birds and turtles#######
  ###########################################################################
  find_transitions = function(Parent_node,threshold, STATES=c("-", "ATGC"), Like = Likelihoods){
    Transition =vector()
    #Get all descendants from a given node
    Desc = sort(c(Descendants(Tree3, Parent_node, type = "all"),Parent_node))
    #Get P for all those descendents
    LLH = Like[Desc,]
    as_tibble(LLH) %>% mutate(N=Desc) -> LLH
    Desc = c(Desc[Desc>leaf_size], Desc[Desc<leaf_size])
    #Go descendant by descendant and check if current state did change from its closest ancestor
    for (n in Desc){
      Node = Ancestors(Tree3,n, type="parent")
      if (length(Node) == 0){ next}
      Prob = filter(LLH,N==n)
      Prob_ancest = filter(LLH,N==Node)
      if (dim(Prob_ancest)[1] == 0){ next}
      P1 =0
      P2 =0
      for (m in seq(1:length(STATES))){
        P1_a = Prob[m]
        P2_a = Prob_ancest[m]
        if (P1_a > P1){
          P1 = P1_a
          state = STATES[m]
        }
        if (P2_a > P2){
          P2 = P2_a
          state_a = STATES[m]
        }
      }
      if (P1 < threshold){ next }
      if (P2 < threshold){ next }
      if(! state_a == state ){
        if(n > length(Tree3$tip.label)){
          Name = Tree3$node.label[n-length(Tree3$tip.label)]
        }else{ Name = Tree3$tip.label[n]}
        Transition = c(Transition, paste(c(state_a,state),collapse=""))
        #if(n<length(Tree3$tip.label)+1){
        print(paste(c("Transition",state_a, "to", state, "in node",n, "named", Name),collapse=" "))}
      #}
      
    }
    return(Transition)
  }
  
  t1 = find_transitions(Parent_node = Turtle + leaf_size, threshold = 0.9)
  lapply(t1, FUN = function(x){ substr(x,1,1) == "-" }) %>% as_vector() -> gains ; sum(gains) ; length(t1) - sum(gains)
  
  t12 = find_transitions(Parent_node = Turtle + leaf_size, threshold = 0.9, Like =Likelihoods_mk)
  #Turtles have the same transitions
  t2 = find_transitions(Parent_node = Aves + leaf_size, threshold = 0.9)
  lapply(t2, FUN = function(x){ substr(x,1,1) == "-" }) %>% as_vector() -> gains ; sum(gains)  ; length(t2) - sum(gains)
  
  t22 = find_transitions(Parent_node = Aves + leaf_size, threshold = 0.9, Like=Likelihoods_mk)
  lapply(t22, FUN = function(x){ substr(x,1,1) == "-" }) %>% as_vector() -> gains ; sum(gains)  ; length(t22) - sum(gains)
  #Birds have only 5 transitions which are different
  
  #All transitions
  t3 = find_transitions(Parent_node = Ancestor_all + leaf_size, threshold = 0.9)
  lapply(t3, FUN = function(x){ substr(x,1,1) == "-" }) %>% as_vector() -> gains ; sum(gains) ; length(t3) - sum(gains)
  return(t3)
}


Frequency_insertion_bird_order = function(Data_birds, Data_Passeriformes, orders_info = "Howard_Moore_taxonomy.tsv", Jarvis_tree="Chronogram01.TENT.ExAML.tre"){
  #Data_birds and Data_Passeriformes can be obtained by returning them from the function Ancestral_state()
  
  
  #########################################################
  ####################Study bird orders####################
  #########################################################
  
  orders_info = read_tsv(orders_info, col_names=F) #Table of three columns. Column 1 bird taxonomy, 2. Bird Family, 3. Bird Order
  
  #Add taxonomical information for records that have NA. 
  orders_info %>% mutate(X3 = ifelse(X1 %in% Data_Passeriformes$X1, "Passeriformes",X3)) -> orders_info
  orders_info %>% mutate(X3 = ifelse(grepl("Anthropoides",X1), "Gruiformes", X3)) -> orders_info
  
  big_extinct = c("Anomalopteryx_didiformis","Megalapteryx_benhami","Pachyornis_australis","Pachyornis_elephantopus","Pachyornis_mappini","Megalapteryx_didinus","Dinornis_giganteus","Dinornis_robustus","Dinornis_novaezealandiae","Emeus_crassus")
  Rhea = c("Pterocnemia_pennata")
  
  orders_info %>% mutate(X3 = ifelse(X1 %in% big_extinct, "Dinornithiformes", X3)) -> orders_info
  orders_info %>% mutate(X3 = ifelse(X1 %in% Rhea, "Rheiformes", X3)) -> orders_info
  orders_info %>% mutate(X3 = ifelse(X1 == "Nannopterum_brasilianus", "Pelecaniformes", X3)) -> orders_info 
  orders_info %>% mutate(X3 = ifelse(X1 == "Megalaima_virens", "Piciformes", X3)) -> orders_info 
  orders_info %>% mutate(X3 = ifelse(X1 == "Eurynorhynchus_pygmeus", "Charadriiformes", X3)) -> orders_info 
  orders_info %>% mutate(X3 = ifelse(X1 == "Dupetor_flavicollis", "Pelecaniformes", X3)) -> orders_info 
  orders_info %>% mutate(X3 = ifelse(X1 == "Harpagornis_moorei", "Accipitriformes", X3)) -> orders_info	
  orders_info %>% mutate(X3 = ifelse(X1 == "Sarothrura_rufa", "Gruiformes", X3)) -> orders_info	
  orders_info %>% mutate(X3 = ifelse(X1 == "Mullerornis_agilis", "Aepyornithiformes", X3)) -> orders_info	
  
  #Add order information to the table with ND3-174 status on birds
  unique(orders_info) %>% arrange(X1) %>% filter(X1 %in% Data_birds$X1) -> orders_info
  unique(Data_birds) %>% arrange(X1) %>% filter() %>% mutate(Family = orders_info$X2, Order= orders_info$X3) -> Data_birds
  
  ####Get information about Frequency of insertion per bird order
  Data_birds %>% group_by(Order, X2) %>% summarise(Number = n()) %>% spread(X2, Number) -> Information_birds_orders
  Information_birds_orders %>% mutate(c = ifelse(is.na(c), 0, c),g = ifelse(is.na(g), 0, g),t = ifelse(is.na(t), 0, t), `-`=ifelse(is.na(`-`), 0, `-`) ) %>%
    mutate(Total = `-`+c + g + t, Frq_gap = ifelse(is.na(`-`),0, `-`/Total)) -> Order_insertion
  
  Palaeognathae = c("Casiariiformes","Rheiformes","Tinamiformes","Strurhioniformes","Apterygiformes","Casuariiformes","Aepyornithiformes")
  Order_insertion %>% filter(Order %in% Palaeognathae) -> Pale
  Frq_gap2= sum(Pale$`-`)/sum(Pale$Total)
  Pale2 = tibble(Order="Palaeognathae","-"=NA,"c"=NA,"g"=NA,"t"=NA,Total=sum(Pale$Total),Frq_gap=Frq_gap2)
  bind_rows(Order_insertion,Pale2)-> Order_insertion
  ################################################
  #####Make figure based on Jarvis topology#######
  ################################################
  Jarvis = read.tree(Jarvis_tree)
  INFO = read_delim("Howard_Moore_taxonomy.csv", delim = ";") #Howard&Moore 
  
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
  Jarvis$tip.label[Jarvis$tip.label=="Sphenisciformes24"] = "Sphenisciformes"#"Apterygiformes"
  drop.tip(Jarvis, c("Caprimulgiformes14","Caprimulgiformes15","Passeriformes35","Passeriformes36","Passeriformes37","Passeriformes34")) -> p
  Jarvis$tip.label[Jarvis$tip.label=="Struthioniformes"] = "Palaeognathae"
  
  Order_insertion %>% mutate(ifelse(Order== "Apterygiformes","Sphenisciformes", Order )) -> Order_insertion
  
  match(Order_insertion$Order, Jarvis$tip.label)
  Order_insertion %>% filter(Order %in% Jarvis$tip.label) -> Order_insertion
  keep.tip(Jarvis, Order_insertion$Order) -> Jarvis
  
  d= fortify(Jarvis)
  d = subset(d, isTip)
  MATCH = with(d, label[order(y, decreasing=T)])
  
  Order_insertion[match(MATCH,Order_insertion$Order),] -> Order_insertion
  Order_insertion[seq(dim(Order_insertion)[1],1),] -> Order_insertion
  
  
  ggtree(Jarvis, layout = "rectangular")  + theme_tree2() + geom_tiplab(offset = -30,vjust = 0.2,size=2 )  -> p  #+ labs(caption = "Divergence time")
  revts(p) + theme(plot.margin=margin(60,60,60,60)) -> p #+ xlim(-105,14) -> p
  facet_plot(p, panel="Frequency absence of insertion", data=Order_insertion, geom= geom_barh, mapping=aes(x=Order_insertion$Frq_gap),stat="identity") -> p2
  facet_plot(p2,panel="Frequency absence of insertion", data=Order_insertion, geom=geom_text, mapping=aes(x= -0.05, label=Order_insertion$Total), size=1.5 ) -> p3
  print(p3)
  ggsave(plot = p3, "FIGURES/Bird_orders.pdf", scale=2)
  ggsave(plot = p3, "FIGURES/Bird_orders.png", scale=2)
}

Process_entrophy = function(Shanon_file = "Shanon.tsv"){
  #Shannon.tsv is the output from Calculate_diversity.py
  ################################
  #####Plot entrophy #############
  ###############################
  data_entrophy = read_tsv(Shanon_file)
  data_entrophy %>% mutate(Information = log2(4) - Entropy, Information_codon = log2(4) - codon_bin) -> data_entrophy
  data_entrophy %>% filter(Position == 174)
  
  BIN = 1
  data_entrophy %>% mutate(Gap_status = ifelse(Gap_status=="yes", "Gap", "Insertion")) -> data_entrophy
  
  R2_gap = data_entrophy %>% filter(Gap_status == "Gap") %>% filter(Position %in% c(163:180)) %>% drop_na()
  R2_insertion = data_entrophy %>% filter(Gap_status == "Insertion") %>% filter(Position %in% c(163:180))
  
  ggplot(R2_insertion) + geom_density(aes(x=Information))
  
  wilcox.test(R2_gap$Information, R2_insertion$Information) -> Test
  mean(R2_gap$Information)/mean(R2_insertion$Information) -> Fold_change
  print("Testing difference in conservation between sequences with insertion and without it")
  print(Test)
  print(Fold_change)
  
  Legend = c(NULL)
  for (n in seq(1,dim(data_entrophy)[1])){
    entry = data_entrophy[n,]
    L = paste(c(entry$Gap_status,entry$Group), collapse="_")
    Legend = c(Legend,L)
  }
  
  data_entrophy %>% mutate(Name = Legend) -> data_entrophy
  data_entrophy %>% ggplot(aes(x=Name, y=1, fill=Information^3)) + geom_bar(stat="identity")  + scale_fill_gradient(low="black", high="white")+
    theme_bw() +geom_hline(yintercept= 174, linetype="dashed",  color = "red", size=0.1) + 
    geom_hline(yintercept= 164, linetype="dashed",  color = "black", size=0.5) + 
    geom_hline(yintercept= 180, linetype="dashed",  color = "black", size=0.5)+ coord_flip()+
    labs(fill = "Information content (to the cube)") + ylab(label = "ND3 position") -> Information_content
  
  print(Information_content)
  ggsave("FIGURES/Information_content.pdf", Information_content, scale= 4)
  ggsave("FIGURES/Information_content.png",Information_content,scale = 4)
}

Study_codon_frequencies = function(Path_file = "codons.tsv"){
  #condons.tsv is obtained from Study_conds_interstion2.py
  ###############################
  #######Codon analysis##########
  ###############################
  
  codons = read_tsv(Path_file)
  codons %>% filter(! Codon == "---") -> codons
  
  codons %>% mutate(Codon_position = ifelse(Sequence=="regular",Codon_position-4,Codon_position)) -> codons
  codons %>% mutate(Codon_position = ifelse(Sequence=="frameshift",Codon_position-1,Codon_position)) -> codons
  codons %>% mutate(Sequence = ifelse(Sequence=="frameshift","Corrected frame","0-frame")) -> codons
  
  codons$Codon_position = factor(as.character(codons$Codon_position))
  codons$Codon_position = factor(codons$Codon_position,levels(codons$Codon_position)[c(3,2,1,4,5,6,7,8)])
  
  
  
  N=names(getGeneticCode("SGC1"))
  V = as.vector(getGeneticCode("SGC1"))
  
  tibble(codon = N, AminoAcid=V) -> amino_translation
  v = vector()
  for (value in codons$Codon){
    amino_translation %>% filter(codon == toupper(value)) -> Value
    AA = Value$AminoAcid
    if (length(AA) == 0){v = c(v,NA)}
    else{ v = c(v,AA)}
  }
  
  codons %>% mutate(AminoAcid = v) -> codons
  codons %>% ggplot(aes(group=factor(Codon),x=Codon_position,y=Counts, fill=AminoAcid,label=Codon)) + 
    geom_bar(stat="identity",position='dodge',col="black") + facet_wrap(~Sequence, scales="free",dir="v") +
    theme_bw() + geom_text_repel(aes(x=Codon_position), size= 3,position = position_dodge(width = 1),direction="y") -> fig_codon
  ggsave("FIGURES/Codon_freq.pdf")
  ggsave("FIGURES/Codon_freq.png")
  print(fig_codon)
}


setwd("/Users/Sergio/Resilio Sync/Transfer/Master/GigaScience/GigaScience/Second_round")

Ancestral_state(Path_data ="Diapsida_table.tsv") -> bird_frequency #Works
Ancestral_state_parsimony(Path_data ="Diapsida_table.tsv") -> Total_transitions #Works

Frequency_insertion_bird_order(bird_frequency[[1]], bird_frequency[[2]],orders_info = "Howard_Moore_taxonomy.tsv") #Missing files on birds
Process_entrophy() #Works, need to do the logo
Study_codon_frequencies(Path_file="scripts/codons.tsv") #Works


D = read_tsv("Diapsida_table.tsv", col_names=F)
Tree = read.tree("Total_tree.tree")
D %>% filter(X1 %in% Tree$tip.label) %>% group_by(X3) %>% summarise(n())
D %>% filter(X1 %in% Tree$tip.label) %>% group_by(X2 == "-")  %>% summarise(n())
D %>% filter(X1 %in% Tree$tip.label) %>% group_by(X2)  %>% summarise(n())


D %>% filter(X1 %in% Tree$tip.label) %>% group_by(X2 == "-", X3)  %>% summarise(n())
