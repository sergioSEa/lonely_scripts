format_metaphlanoutput_df = function(path, ID_name = 'ID' ){
  DB = readLines(path, n = 1)
  df = readr::read_tsv(path, skip=1)
  df %>% as.data.frame() %>% tibble::column_to_rownames('clade_name') %>% t() %>% as.data.frame() %>%
    tibble::rownames_to_column(ID_name) %>% tibble::as_tibble() -> df
  return( list(abundance_table = df, mp_version = DB ) )
}

extract_taxonomy <- function(taxa_list) {
  taxonomy_levels <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__")
  column_names <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "SGB")
  
  parsed <- lapply(taxa_list, function(entry) {
    values <- sapply(taxonomy_levels, function(prefix) {
      match <- grep(paste0("^", prefix), entry, value = TRUE)
      if (length(match) > 0) sub(prefix, "", match) else NA
    })
    setNames(values, column_names)
  })
  
  as.data.frame(do.call(rbind, parsed), stringsAsFactors = FALSE)
}


retrieve_taxonomic_level = function(df, taxonomy='t', clean_names=T,sep='\\|'){
  if (! grepl('_$', taxonomy) ){
    taxonomy = paste0(taxonomy, '__')
  }
  Taxa = colnames(df)
  clean = Taxa %>% str_split(sep) %>% sapply( function(taxa){ taxa[length(taxa)] } )
  Keep = grepl(taxonomy, clean )
  Keep[1] = T #ID column
  Keep[2] = T #UNCLASSIFIED column
  
  df %>% select( Taxa[Keep] ) -> df
  if(clean_names == T){
    colnames(df) = clean[Keep]
  }
  complete_taxonomy = Taxa[Keep] %>% str_split(sep) %>% extract_taxonomy() %>% as_tibble() %>% drop_na()
  
  return( list(abundance_table = df, taxonomy_table = complete_taxonomy ) )
}

add_spname_to_sgb = function(output_retrieve_taxonomy){
  #Make sure input is a list of two
  if ( length(output_retrieve_taxonomy) != 2 ){ return(NULL) }
  #Make sure taxonomy level is SGB (t__)
  df = output_retrieve_taxonomy[[1]]
  SGBs = grepl('^t__', colnames(df) )
  if (  sum( SGBs )   == 0  ){ return(NULL) }
  #Get species names
  Species_table =  output_retrieve_taxonomy[[2]]
  colnames(df)[SGBs] %>% sapply( function(x){ Species_table$species[Species_table$SGB== str_split(x, '__')[[1]][2] ]    }  ) -> Sp_names
  
  colnames(df)[SGBs] = paste(Sp_names, names(Sp_names), sep = "|")
  return(df)
}
make_barplots = function(df, id_col = 'ID', top_taxa = 10, Keep_taxa= c(), rm_cols = c()) {
  if (length(df) != 0) {
    df <- dplyr::select(df, -all_of(rm_cols))
  }
  
  # Get top taxa by average abundance
  taxa_avg_ab <- df %>%
    dplyr::select(-all_of(id_col)) %>%
    apply(2, mean)
  taxa_do <- sort(taxa_avg_ab, decreasing = TRUE)[1:top_taxa]
  taxa_do = c(names(taxa_do), Keep_taxa) %>% unique() 
  # Subset dataframe to include only top taxa and ID column
  df_top <- df %>%
    dplyr::select(all_of(c(id_col, names(taxa_avg_ab))))
  
  # Convert to long format
  df_long <- df_top %>%
    tidyr::pivot_longer(-all_of(id_col), names_to = "Taxa", values_to = "Abundance")
  
  # Label non-top taxa as "Other"
  df_long <- df_long %>%
    dplyr::mutate(Taxa = ifelse(Taxa %in% taxa_do, Taxa, "Other"))
  
  # Normalize to relative abundance
  df_long <- df_long %>%
    dplyr::group_by(!!sym(id_col)) %>%
    dplyr::mutate(RelAbundance = Abundance / sum(Abundance)) %>%
    dplyr::ungroup()
  
  # Plot
  Plot = ggplot2::ggplot(df_long, ggplot2::aes(x = !!sym(id_col), y = RelAbundance, fill = Taxa)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = "Relative Abundance", x = id_col, fill = "Taxa") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))  + ggplot2::scale_fill_manual(values=c25)
  
  return(Plot)
}
