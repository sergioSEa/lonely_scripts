#Rarefication

setwd("~/PhD/WORK/Collaborations/Debby/Diet")
library(vegan)
library(tidyverse)
library(ape)

set.seed(1205)

Genus_abundance = read_csv("Merged_Genus_abundance.csv")
ASV_abundance  = read_csv("Matrix_abundance.csv")
read.tree("16S_ASV.tree") -> ASV_tree
metadata = read_tsv("Metadata.txt")


Rarefy_reads = function(M){
  Number_reads = apply(select(M, -ID), MARGIN= 1, FUN = function(x){ sum(x) })
  ggplot() + geom_bar(aes(x=M$ID, y=Number_reads), stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ##Check difference in number of species before and after rarefy
  S <- specnumber(select(M, -ID)) # observed number of species
  Srare <- rarefy(x=select(M,-ID),MARGIN=1, sample=min(Number_reads))
  plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
  #rarecurve(select(M,-ID), sample=50317)
  
  ##Perform rarefy
  vegan::rrarefy(x=select(M,-ID), sample=min(Number_reads)) %>% as_tibble() %>% mutate(ID = M$ID) -> Abundance_rarefy
  
  return(Abundance_rarefy)
}
Make_tree_genus2 = function(Tree, Genus_name, All_samples, metadata, Plot = F){
  #Make phylogenetic tree for a given Name and statistical test on phylogenetic distances
  
  if (Genus_name != "ALL"){
    #Get tree for Genus of interest
    Filter_branches = Tree$tip.label[grepl(Genus_name, Tree$tip.label)]
    Subset_Tree = keep.tip(Tree, Filter_branches)
    
    #Get Samples of Genus of interest
    All_samples %>% select(c("ID",Filter_branches))  -> Subset_samples
  } else{  Subset_Tree = Tree ; Subset_samples = All_samples }
  if (length(Subset_Tree$tip.label) < 3){ return(c(NA, NA, NA)) }
  #sapply(Subset_samples$ID, FUN = function(x){ str_split(x, "_")[[1]][[3]] }) -> New_names
  #Subset_samples %>% mutate(ID=as.vector(New_names)) %>% filter(ID %in% metadata$`Sample Name`) -> Subset_samples
  Subset_samples %>% filter(ID %in% metadata$`Sample Name`) -> Subset_samples
  Subset_samples[2:(dim(Subset_samples)[2])] = as.integer(Subset_samples[2:(dim(Subset_samples)[2])] != 0) 
  Subset_samples %>% select(-ID) %>% t() %>% as_tibble() %>% mutate(label=colnames(select(Subset_samples,-ID))) %>% `colnames<-`(c(Subset_samples$ID, "label")) %>% select(c("label",Subset_samples$ID)) -> Subset_samples2
  #Add mouse Number to metadata
  metadata %>% mutate(Mouse = as.vector(sapply(metadata$`Sample Name`, FUN= function(X){ str_replace(strsplit(X,"\\.")[[1]][[2]],"F|O", "") }))) -> metadata
  arrange(metadata, `Sample Name`) %>% filter(`Sample Name` %in% Subset_samples$ID) -> metadata
  arrange(Subset_samples, ID) %>% mutate(Group = metadata$`Group Name`, Mouse=metadata$Mouse, Origin=metadata$`*Tissue`, Diet=as.numeric(as.factor(metadata$Diet))) -> Data_merg
  Controls_same = vector() ; Controls_other= vector() ; Treated_same = vector() ; Treated_other = vector()
  
  for (MOUSE in unique(metadata$Mouse)){
    Data_merg %>% filter(Mouse == MOUSE) -> Data_merg2
    Data_merg2 %>% select(-c("Group", "Mouse", "Origin", "Diet")) -> ID_x
    equals_two = function(x){x==2}
    equals_one = function(x){x==1}
    ID_x %>% summarise_if(is.numeric,sum) %>% as_vector() -> Tip
    
    sum(Tip == 1 ) -> Other_tip
    sum(Tip == 2 ) -> Same_tip
    

    if (unique(Data_merg2$Diet) == 1){ 
      Controls_same = c(Controls_same, Same_tip) ; Controls_other = c(Controls_other, Other_tip)  
    }else{ Treated_same = c(Treated_same, Same_tip) ; Treated_other = c(Treated_other, Other_tip)   }
  }
  data.frame(Group_control = c(sum(Controls_same), sum(Controls_other)), Group_treated= c(sum(Treated_same), sum(Treated_other))) -> Table_transfer
  #data.frame(Group_control = c(median(Controls_same), median(Controls_other)), Group_treated= c(median(Treated_same), median(Treated_other))) -> Table_transfer
  #data.frame(Group_control = c(sum(Controls_same), sum(Controls_other))) -> Table_transfer
  
  rownames(Table_transfer) = c("Same_tip", "Other_tip")
  Test = chisq.test(Table_transfer)
  Summary_results = c(Test$p.value, sum(Controls_same)/sum(Controls_other), sum(Treated_same)/sum(Treated_other))
  #Summary_results = c(sum(Controls_same)/sum(Controls_other))
  
  if (Plot != F){
    library(phyloseq)
    #Prepare Phyloseq object for plotting
    Subset_samples2 %>% as.data.frame() %>% column_to_rownames("label") -> ASV_table
    taxmat = matrix(rownames(ASV_table))
    rownames(taxmat) = rownames(ASV_table)
    metadata %>% mutate(Group = `Group Name`) -> metadata
    sampledata = metadata %>% mutate(Treatment = as.factor(Diet)) %>% as.data.frame() %>% column_to_rownames("Sample Name")
    otu_table(ASV_table, taxa_are_rows = TRUE) -> ASV_table
    tax_table(taxmat) -> taxmat
    phyloseq::phyloseq(ASV_table, taxmat) -> phylobjetc
    merge_phyloseq(phylobjetc,Subset_Tree, sample_data(sampledata)) -> phylobjetc
    #taxa_names
    #plot_tree(phylobjetc,  use.edge.length = F, color="Origin", shape="Treatment", ladderize="left", plot.margin=0.3) -> TREE_PLOT
    New_tree(phylobjetc, color="Group", ladderize=F, plot.margin=0.3, LABEL="Mouse", text.size = 3, nodelab= nodeplotblank, base.spacing=0.1) -> TREE_PLOT
    str_replace(Genus_name, "/", "_OR_") -> Genus_name
    trim(Genus_name) -> Genus_name
    FILE = paste(c("Figures/TREE_", Genus_name,".pdf") , collapse="")
    print(FILE)
    print(TREE_PLOT)
    #ggsave(filename =  FILE, plot = TREE_PLOT, scale = 3)
    
  }
  return(Summary_results)
}
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
New_tree = function(physeq, method = "sampledodge", nodelabf = NULL, 
                    color = NULL, shape = NULL, size = NULL, min.abundance = Inf, 
                    label.tips = NULL, text.size = NULL, sizebase = 5, base.spacing = 0.02, 
                    ladderize = FALSE, plot.margin = 0.2, title = NULL, treetheme = NULL, 
                    justify = "jagged", LABEL =NULL) {
  library(data.table)
  fix_reserved_vars = function(aesvar) {
    aesvar <- gsub("^abundance[s]{0,}$", "Abundance", 
                   aesvar, ignore.case = TRUE)
    aesvar <- gsub("^OTU[s]{0,}$", "OTU", aesvar, 
                   ignore.case = TRUE)
    aesvar <- gsub("^taxa_name[s]{0,}$", "OTU", 
                   aesvar, ignore.case = TRUE)
    aesvar <- gsub("^sample[s]{0,}$", "Sample", 
                   aesvar, ignore.case = TRUE)
    return(aesvar)
  }
  if (!is.null(label.tips)) {
    label.tips <- fix_reserved_vars(label.tips)
  }
  if (!is.null(color)) {
    color <- fix_reserved_vars(color)
  }
  if (!is.null(shape)) {
    shape <- fix_reserved_vars(shape)
  }
  if (!is.null(size)) {
    size <- fix_reserved_vars(size)
  }
  if (is.null(phy_tree(physeq, FALSE))) {
    stop("There is no phylogenetic tree in the object you have provided.\n", 
         "Try phy_tree(physeq) to see for yourself.")
  }
  if (!inherits(physeq, "phyloseq")) {
    method <- "treeonly"
  }
  treeSegs <- tree_layout(phy_tree(physeq), ladderize = ladderize)
  edgeMap = aes(x = xleft, xend = xright, y = y, yend = y)
  vertMap = aes(x = x, xend = x, y = vmin, yend = vmax)
  p = ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) + 
    geom_segment(vertMap, data = treeSegs$vertDT)
  if (is.null(text.size)) {
    text.size <- manytextsize(ntaxa(physeq))
  }
  if (!is.null(label.tips) & method != "sampledodge") {
    labelDT = treeSegs$edgeDT[!is.na(OTU), ]
    if (!is.null(tax_table(object = physeq, errorIfNULL = FALSE))) {
      taxDT = data.table(tax_table(physeq), OTU = taxa_names(physeq), 
                         key = "OTU")
      labelDT = merge(x = labelDT, y = taxDT, by = "OTU")
    }
    if (justify == "jagged") {
      labelMap <- aes_string(x = "xright", y = "y", 
                             label = label.tips, color = color)
    }
    else {
      labelMap <- aes_string(x = "max(xright, na.rm=TRUE)", 
                             y = "y", label = label.tips, color = color)
    }
    p <- p + geom_text(labelMap, data = labelDT, size = I(text.size), 
                       hjust = -0.1, na.rm = TRUE)
  }
  if (is.null(nodelabf)) {
    nodelabf = howtolabnodes(physeq)
  }
  p = nodelabf(p, treeSegs$edgeDT[!is.na(label), ])
  p = nodelabf(p, treeSegs$vertDT[!is.na(label), ])
  if (is.null(treetheme)) {
    treetheme <- theme(axis.ticks = element_blank(), axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), axis.title.y = element_blank(), 
                       axis.text.y = element_blank(), panel.background = element_blank(), 
                       panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  }
  if (inherits(treetheme, "theme")) {
    p <- p + treetheme
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (method != "sampledodge") {
    return(p)
  }
  dodgeDT = treeSegs$edgeDT[!is.na(OTU), ]
  dodgeDT = merge(x = dodgeDT, y = data.table(psmelt(physeq), 
                                              key = "OTU"), by = "OTU")
  if (justify == "jagged") {
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  if (!is.null(color) | !is.null(shape) | !is.null(size)) {
    setkeyv(dodgeDT, cols = c("OTU", color, shape, 
                              size, LABEL))
  }
  else {
    setkey(dodgeDT, OTU, Sample)
  }
  dodgeDT[, `:=`(h.adj.index, 1:length(xright)), by = OTU]
  if (justify == "jagged") {
    dodgeDT[, `:=`(xdodge, (xright + h.adj.index * 
                              base.spacing * max(xright, na.rm = TRUE)))]
  }
  else {
    dodgeDT[, `:=`(xdodge, max(xright, na.rm = TRUE) + 
                     h.adj.index * base.spacing * max(xright, na.rm = TRUE))]
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  dodgeMap <- aes_string(x = "xdodge", y = "y", 
                         color = color, fill = color, shape = shape, label = LABEL)
  p <- p + geom_text(dodgeMap, data = dodgeDT, na.rm = TRUE, size= text.size)
  if (!is.null(size)) {
    p <- p + scale_size_continuous(trans = log_trans(sizebase))
  }
  if (any(dodgeDT$Abundance >= min.abundance[1])) {
    pointlabdf = dodgeDT[Abundance >= min.abundance[1], ]
    p <- p + geom_text(mapping = aes(xdodge, y, label = Abundance), 
                       data = pointlabdf, size = text.size, na.rm = TRUE)
  }
  if (!is.null(label.tips)) {
    tiplabDT = dodgeDT
    tiplabDT[, `:=`(xfartiplab, max(xdodge)), by = OTU]
    tiplabDT <- tiplabDT[h.adj.index == 1, .SD, by = OTU]
    if (!is.null(color)) {
      if (color %in% sample_variables(physeq, errorIfNULL = FALSE)) {
        color <- NULL
      }
    }
    labelMap <- NULL
    if (justify == "jagged") {
      labelMap <- aes_string(x = "xfartiplab", y = "y", 
                             label = label.tips, color = color)
    }
    else {
      labelMap <- aes_string(x = "max(xfartiplab, na.rm=TRUE)", 
                             y = "y", label = label.tips, color = color)
    }
    p <- p + geom_text(labelMap, tiplabDT, size = I(text.size), 
                       hjust = -0.1, na.rm = TRUE)
  }
  min.x <- -0.01
  max.x <- dodgeDT[, max(xright, na.rm = TRUE)]
  if ("xdodge" %in% names(dodgeDT)) {
    max.x <- dodgeDT[, max(xright, xdodge, na.rm = TRUE)]
  }
  if (plot.margin > 0) {
    max.x <- max.x * (1 + plot.margin)
  }
  p <- p + scale_x_continuous(limits = c(min.x, max.x))
  return(p)
}
Filter_abundance = function(M, threshold=0.2, min_abundance = 0){
  #Filter species not seen in at least treshold*100% of samples
  apply(select(M, -ID), MARGIN=2, FUN = function(x){ (sum(x != 0)/length(x)) >= 0.2 }) -> Filter_vector
  F_names = colnames(select(M, -ID))[!Filter_vector]
  M %>% select(!F_names) -> M
  #Filter species with mean relative abundance lower than min_abundance
  apply(select(M, -ID), MARGIN=1, FUN = function(x){ sum(x) }) -> Total_reads_sample
  as_tibble(select(M, -ID)/Total_reads_sample) %>% summarise_all(mean) -> All_means
  colnames(select(M, -ID))[as_vector(All_means) > min_abundance] -> Taxa_not_to_filter
  M %>% select(c(ID, Taxa_not_to_filter)) -> M
  return(M)
}




Rarefy_reads(Genus_abundance) -> rarefied_genus
#write_csv(rarefied_genus, "Merged_Genus_abundance_rarefy.csv")
Rarefy_reads(ASV_abundance) -> rarefied_ASV



sapply(rarefied_ASV$ID, FUN = function(x){ str_split(x, "_")[[1]][[3]] }) -> ID
rarefied_ASV$ID = as.vector(ID)

str_replace(metadata$`*Sample Name`, "_", ".") -> NAMES
metadata %>% mutate(`Sample Name`= NAMES) -> metadata
metadata %>% mutate(Time_Diet = paste(`Period Ttreatment`, Diet, sep="-")) -> metadata

Result_tranfer_total = tibble()
rarefied_genus %>% select(-"Unknown|Genus") -> Abundance_rarefy
Comparisons = list(c("W4-LFD", "W4-HFD"), c("W4-HFD", "W4-HFC"),c("W8-LFD", "W8-HFD"), c("W8-HFD", "W8-HFC"),c("W12-LFD", "W12-HFD"), c("W12-HFD", "W12-HFC") )
#Treatment 1 is gonna be always the one in the right

for (T_n_D in Comparisons){
  Meta = arrange(filter(metadata, Time_Diet %in% T_n_D), `Sample Name`)
  Abundance_rarefy %>% arrange(ID) %>% filter(ID %in% Meta$`Sample Name`) -> Abundance_data
  Abundance_data = Filter_abundance(Abundance_data)
  for (Genus in colnames(Abundance_data)){
    print(Genus)
    if (Genus == "ID"){ next }
    S = strsplit(Genus, split = "\\|")
    Genus = S[[1]][1]
  
    Results_transfer =  Make_tree_genus2(Tree = ASV_tree, Genus_name = Genus, All_samples = rarefied_ASV, metadata = Meta,Plot = F)
    data.frame(t(Results_transfer)) -> Results_transfer
    colnames(Results_transfer) = c("Pvalue", "Proportion_identical_Treated", "Proportion_identical_Control")
    Results_transfer %>% mutate(Comparison = paste(T_n_D, collapse="|")) -> Results_transfer
    Results_transfer %>% mutate(Genus_name = Genus, Diff= Proportion_identical_Control - Proportion_identical_Treated) -> Results_transfer
    Result_tranfer_total = rbind(Result_tranfer_total,Results_transfer)
  }
}
write_tsv(Result_tranfer_total,"All_transfer_results.tsv")

as_tibble(Result_tranfer_total) %>% arrange(Pvalue) %>% mutate(Time = ifelse(grepl("4", Comparison), 4, ifelse(grepl("8", Comparison), 8,12))) -> Result_tranfer_total
Result_tranfer_total %>% mutate(Diet = ifelse(grepl("HFC", Comparison), "HFC", "HFD")) -> Result_tranfer_total
Result_tranfer_total %>% filter(Genus_name == "Lachnospira") %>% ggplot() + geom_point(aes(x=Time, y=Diff, col=Diet)) + theme_bw()
                                  