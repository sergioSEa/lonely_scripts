library(tidyverse)
library(ape)
library(phytools)
library(rotl)

dataset = read_tsv("Postion_nucleotide_diapsida.tsv",col_names=F)


#############################################################
#######################Get tree##############################
#############################################################



taxa = gsub(" ","_",dataset$X1)
resolved_names = vector()

resolved_names = tibble(search_string=character(),unique_name=character(),approximate_match=logical(),ott_id=numeric(),is_synonym=logical(),flags=character(),number_matches=numeric())
resolved_names = as.data.frame(resolved_names)

###For each species, get otol name
for (item in 1:length(taxa)){
  tryCatch({
   Name = tnrs_match_names(taxa[item])
    resolved_names = rbind(resolved_names, Name)
  }, error=function(e){cat(paste(c("ERROR: ",item, taxa[item], collapse=F)) , "\n")})  
 
}
write_csv(x = resolved_names, path = "R_names_itol.csv")


#Tree with all vertebrades
my_tree  = read.tree("opentree11.4_tree/labelled_supertree/labelled_supertree_ottnames.tre")
#Extract diapsida
my_tree = extract.clade(my_tree, "mrcaott59ott1662")


add_tonumber = function(vector_id){ paste(c("ott",vector_id), collapse="") }

v_d = vector()
#Prepare names in teh same format than seen in the Otol tree. Save them in v_d  vector
for (i in seq(dim(resolved_names)[1])){
  I = resolved_names[i,]
  name = I$search_string
  name = paste(c(toupper(substring(name, 1, 1)), substring(name, 2)), collapse="")
  vector_id =  I$OTT_name
  d = paste(c(name,vector_id), collapse="_")
  v_d = c(v_d,d)
}
naming = sapply(resolved_names$ott_id,FUN=add_tonumber)

#Check for names in the Otol tree
resolved_names %>% mutate(OTT_name = naming, Match_name = v_d) -> resolved_names
resolved_names %>% filter(! Match_name %in% my_tree$tip.label) -> not_found


get_N = function(NAME, OUT){
  N = strsplit(NAME,split = "_")
  ott = N[[1]][3]
  N = paste(N[[1]][1:2],collapse=" ")
  if (OUT == "ott"){return(ott)
    }else{ return(N)}
}

###Try to match the unmatched either by NAME or OTT, then change the tip names so they all look the same.
Name = unlist(lapply(FUN = get_N,X = my_tree$tip.label, OUT="name"))
OTTs = unlist(lapply(FUN = get_N,X = my_tree$tip.label, OUT="ott"))

resolved_names %>% mutate(unique_name = ifelse(search_string=="chinemys_reevesi","Chinemys_reevesi",unique_name)) -> resolved_names

Remove = which((toupper(resolved_names$search_string) == str_replace(toupper(resolved_names$unique_name)," ","_")) == F)
resolved_names = resolved_names[-Remove,]

resolved_names %>% filter(unique_name %in% Name) -> RN
resolved_names %>% filter(! unique_name %in% Name) %>% filter(!search_string == "chinemys_reevesi") -> not_found2
#Species not found in the tree somehow
dim(not_found2)


####RN are the ones found
RN  %>% filter(search_string == "chinemys_reevesi") -> INCLUDE
RN %>% filter(is_synonym==F) -> RN
RN = rbind(RN, INCLUDE)
my_tree$tip.label = Name

###There are unique names which are repeated length(unique(RN$unique_name)) != length(RN$unique_name)

my_tree = keep.tip(my_tree, RN$unique_name)

write.tree(my_tree, "Total_tree.tree")

resolved_names %>% filter(! search_string %in% RN$search_string) -> dataset_notused
write_tsv(dataset_notused,"discarted.tsv")
write_tsv(RN,"nondiscarted.tsv")

