
#' Convert a phyloseq object to an agglomerated (melted) dataframe.
#' @param ps a phyloseq object
#' @param sample.col.name the sample_data() column corresponding rownames(otu_table()), as a string
#' @param otu_col the desired name of the otu id column
#' @param count_col the desired name of the count column
#' @import reshape2
phyloseq_to_agglomerated <- function(ps, sample_col, otu_col, count_col) {
  otus <- otu_table(ps)
  otus <- otus[, colSums(otus) > 0]
  taxa <- tax_table(ps)
  taxa <- taxa[rownames(taxa) %in% colnames(otus), ]
  taxa.df <- as.data.frame(taxa@.Data)
  taxa.df$otu_id <- rownames(taxa)
  agg <- reshape2::melt(otus, varnames=c(sample_col, otu_col), value.name=count_col) %>%
    left_join(sample_data(ps)) %>%
    left_join(taxa.df)
}

#' Convert an agglomerated dataframe to a phyloseq object.
#' @param agg the agglomerated dataframe
#' @param otu_col the unquoted name of the otu id column in agg
#' @param sample_col the unquoted name of the sample id column in agg
#' @param count_col the unquoted name of the counts column in agg
#' @param ... remaining unquoted names of taxonomic ranks (i.e. Kingdom:Species)
#' @import phyloseq
#' @import dplyr
agglomerated_to_phyloseq <- function(agg, otu_col, sample_col, count_col, ...) {
  otu_col <- enquo(otu_col)
  sample_col <- enquo(sample_col)
  count_col <- enquo(count_col)
  taxa_ranks <- quos(...)

  # Rebuild count matrix
  otu <- select(agg, !! otu_col, !! sample_col, !! count_col) %>%
    tidyr::spread(!! sample_col, !! count_col, fill=0) %>%
    filter(!is.na(!! otu_col)) %>%
    as.data.frame()
  rownames(otu) <- otu[[quo_name(otu_col)]]
  otu[[quo_name(otu_col)]] <- NULL
  otu.mat <- t(as.matrix(otu))

  # Rebuild taxonomy table
  taxa <- select(agg, !! otu_col, !!! taxa_ranks) %>%
    distinct() %>%
    as.data.frame()
  rownames(taxa) <- taxa[[quo_name(otu_col)]]
  taxa[[quo_name(otu_col)]] <- NULL
  taxa <- as.matrix(taxa)
  taxa <- taxa[match(rownames(taxa), colnames(otu.mat)),]

  # Rebuild sample table
  samples <- select(agg, -c(!! otu_col, !! count_col, !!! taxa_ranks)) %>%
    distinct() %>%
    as.data.frame()
  rownames(samples) <- samples[[quo_name(sample_col)]]
  samples <- samples[match(rownames(samples), rownames(otu.mat)), ]

  phyloseq(
    otu_table(otu.mat, taxa_are_rows = FALSE),
    sample_data(samples),
    tax_table(taxa)
  )
}
