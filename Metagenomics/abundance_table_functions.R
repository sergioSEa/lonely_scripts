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


retrieve_taxonomic_level = function(df, taxonomy='t', clean_names=T){
  if (! grepl('_$', taxonomy) ){
    taxonomy = paste0(taxonomy, '__')
  }
  Taxa = colnames(df)
  clean = Taxa %>% str_split('\\|') %>% sapply( function(taxa){ taxa[length(taxa)] } )
  Keep = grepl(taxonomy, clean )
  Keep[1] = T #ID column
  Keep[2] = T #UNCLASSIFIED column
  
  df %>% select( Taxa[Keep] ) -> df
  if(clean_names == T){
    colnames(df) = clean[Keep]
  }
  complete_taxonomy = Taxa[Keep] %>% str_split('\\|') %>% extract_taxonomy() %>% as_tibble() %>% drop_na()
  
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
