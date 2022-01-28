library(tidyverse)
library('org.Hs.eg.db')


read_delim(file = "~/Desktop/panel_senescence.txt", delim = " ") -> Genes
Summary_stats = readxl::read_excel("~/Resilio Sync/LLD phenotypes for Sergio/Results/Manuscript/Manuscript_2_0620201/Supplementary/Supplementary_Files.xlsx", sheet = "Table_S7.2")
Other_names = read_tsv("~/Resilio Sync/LLD phenotypes for Sergio/scripts/Pathway_enrichment/All_gene_info_broad.tsv")
Other_names %>% mutate(primerid = `Original name`) %>% dplyr::select(entrezgene_id, primerid) -> Other_names
left_join(Summary_stats, Other_names) %>% distinct() -> Summary_stats

#Get entrez ids in mayo clinic senescence panel
mapIds(org.Hs.eg.db, Genes$`Gene(human)`, 'ENTREZID', 'SYMBOL') -> entrez
Genes %>% mutate(Entrez = entrez) %>% distinct() ->  Genes

#22/125
Genes %>% dplyr::filter(Entrez %in% Summary_stats$entrezgene_id)

#GSEA analysis

set.seed(999)
All_cells_results = tibble()
Plots = list()
for (i in seq(length(unique(Summary_stats$Cell_type)))){
  I = unique(Summary_stats$Cell_type)[i]
  Pathways = list()
  Pathways$SENSMAYO = c( unique(Genes$Entrez) )

  Summary_stats %>% filter(Cell_type == I) %>% mutate(z= as.numeric(z)) %>% filter(!is.na(z)) %>% arrange(z) %>% mutate(ID = primerid) -> Ranked
  Ranked %>% dplyr::select(entrezgene_id, z) %>% drop_na() -> Ranked
  as.numeric(Ranked$z) -> Ranked_v
  names(Ranked_v) = Ranked$entrezgene_id
  
  fgsea(pathways = Pathways, stats= Ranked_v )   -> Results
  
  as_tibble(Results)  %>% mutate(Cell = I) -> Results
  rbind(All_cells_results , Results) -> All_cells_results
  
  #Complains about ties because several Z values are identical
  plotEnrichment(Pathways[["SENSMAYO"]], Ranked_v) + labs(title=paste(c("SENSMAYO", I), collapse=":" )) -> PLOT
  Plots[[i]] = PLOT
  #plotGseaTable(pathways = Pathways["SENSMAYO"], stats = Ranked_v, fgseaRes = Results)
  
}
