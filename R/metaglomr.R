#' Example feature counts.
"features"

#' Example sample metadata.
"samples"

#' Example taxonomy annotations
"taxa"


#' Create an agglomerated dataframe.
#'
#' Creates an agglomerated dataframe from feature counts, sample metadata, and taxonomy annotations.
#'
#' @param features a matrix-like object with features as columns and samples as rows
#' @param samples an optional dataframe containing sample metadata, with each sample a separate row
#' @param taxonomy an optional dataframe containing taxonomy information, with each taxa a row and
#'   rownames providing taxa id
#' @param sample_col a string giving the name of the column in `samples` corresponding to the rows
#'   in the features table
#' @param otu_col a string giving the desired name of the otu column in the resulting dataframe
#' @param count_col a string giving the desired name of the counts column in the resulting dataframe
#' @return A dataframe with each row representing a unique feature+sample+count observation, plus
#'   any metadata about the sample and any taxonomic annotation about the feature. This results in a
#'   lot of repetition, but is very easy to manipulate with dplyr and the tidyverse.
#' @export
agglomerated <-
  function(features,
           samples,
           taxonomy,
           sample_col,
           otu_col = "otu_id",
           count_col = "count") {
    agg <- reshape2::melt(features, varnames=c(sample_col, otu_col), value.name=count_col)
    agg[[sample_col]] <- as.character(agg[[sample_col]])
    if (!missing(samples)) {
      samples[[sample_col]] <- as.character(samples[[sample_col]])
      agg <- dplyr::left_join(agg, samples, by=sample_col)
    }
    if (!missing(taxonomy)) {
      if (!otu_col %in% colnames(taxonomy)) {
        taxonomy[[otu_col]] <- rownames(taxonomy)
      } else {
        taxonomy[[otu_col]] <- as.character(taxonomy[[otu_col]])
      }
      agg[[otu_col]] <- as.character(agg[[otu_col]])
      agg <- dplyr::left_join(agg, taxonomy, by=otu_col)
    }
    as_tibble(agg)
  }

#' Extract sample metadata from an agglomerated dataframe
#' @param agg agglomerated dataframe
#' @param sample_col unquoted name of the sample id column
#' @param otu_col unquoted name of the otu id column
#' @param count_col unquoted name of the count column
#' @param ... unquoted names of the taxonomic ranks (i.e. Kingdom:Species)
#' @return the sample metadata, deduplicated
#' @import dplyr
#' @export
get_samples <- function(agg, sample_col, otu_col, count_col, ...) {
  if (missing(sample_col) | missing(otu_col) | missing(count_col)) {
    stop("Must specify sample, otu, and count columns")
  }
  otu_col <- enquo(otu_col)
  sample_col <- enquo(sample_col)
  count_col <- enquo(count_col)
  taxa_ranks <- quos(...)

  samples <- select(agg, -c(UQ(otu_col), UQ(count_col), UQS(taxa_ranks)))
  samples <- as.data.frame(distinct(samples))
  rownames(samples) <- samples[[quo_name(sample_col)]]
  samples
}

#' Extract feature count matrix from agglomerated dataframe
#' @param agg agglomerated dataframe
#' @param sample_col unquoted name of the sample id column
#' @param otu_col unquoted name of the otu id column
#' @param count_col unquoted name of the count column
#' @return the feature count matrix, invisibly
#' @import dplyr
#' @export
get_features <- function(agg, sample_col, otu_col, count_col) {
  if (missing(sample_col) | missing(otu_col) | missing(count_col)) {
    stop("Must specify sample, otu, and count columns")
  }
  otu_col <- enquo(otu_col)
  sample_col <- enquo(sample_col)
  count_col <- enquo(count_col)

  # Rebuild count matrix
  otu <- select(agg, UQ(otu_col), UQ(sample_col), UQ(count_col))
  otu <- tidyr::spread(otu, UQ(sample_col), UQ(count_col), fill=0)
  otu <- filter(otu, !is.na(UQ(otu_col)))
  otu <- as.data.frame(otu)
  rownames(otu) <- otu[[quo_name(otu_col)]]
  otu[[quo_name(otu_col)]] <- NULL
  otu.mat <- t(as.matrix(otu))
  invisible(otu.mat)
}

#' Extract taxonomic annotations from agglomerated dataframe
#' @param agg agglomerated dataframe
#' @param otu_col unquoted name of the otu id column
#' @param ... unquoted names of taxonomic ranks, i.e. Kingdom:Species
#' @return a character matrix of taxonomic annotations
get_taxa <- function(agg, otu_col, ...) {
  if (missing(otu_col) | missing(...)) {
    stop("Must specify otu and taxonomic rank columns")
  }
  otu_col <- enquo(otu_col)
  tax_ranks <- quos(...)

  # Rebuild taxonomy table
  taxa <- select(agg, UQ(otu_col), UQS(tax_ranks))
  taxa <- as.data.frame(distinct(taxa))
  rownames(taxa) <- taxa[[quo_name(otu_col)]]
  taxa[[quo_name(otu_col)]] <- NULL
  taxa <- as.matrix(taxa)
}

#' Convert a phyloseq object to an agglomerated dataframe.
#' @param ps a phyloseq object
#' @param sample_col the sample_data() column corresponding to rownames(otu_table()), as a string
#' @param otu_col the desired name of the otu id column
#' @param count_col the desired name of the count column
#' @export
phyloseq_to_agglomerated <- function(ps, sample_col, otu_col="otu_id", count_col="count") {
  otus <- phyloseq::otu_table(ps)
  otus <- otus[, colSums(otus) > 0]
  taxa <- phyloseq::tax_table(ps)
  taxa <- taxa[rownames(taxa) %in% colnames(otus), ]
  taxa.df <- as.data.frame(taxa@.Data)
  taxa.df$otu_id <- rownames(taxa)
  agg <- reshape2::melt(otus, varnames=c(sample_col, otu_col), value.name=count_col)
  agg <- dplyr::left_join(agg, phyloseq::sample_data(ps))
  dplyr::left_join(agg, taxa.df)
}

#' Convert an agglomerated dataframe to a phyloseq object.
#' @param agg the agglomerated dataframe
#' @param sample_col the name of the sample id column in agg
#' @param otu_col the name of the otu id column in agg
#' @param count_col the unquoted name of the counts column in agg
#' @param ... remaining unquoted names of taxonomic ranks (i.e. Kingdom:Species)
#' @import dplyr
#' @export
agglomerated_to_phyloseq <- function(agg, sample_col, otu_col, count_col, ...) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop(
      "Package \"phyloseq\" needed for this function to work. Please install it from Bioconductor.",
      call. = FALSE
    )
  }
  otu_col <- enquo(otu_col)
  sample_col <- enquo(sample_col)
  count_col <- enquo(count_col)
  tax_ranks <- quos(...)

  # Rebuild count matrix
  features <- get_features(agg, UQ(sample_col), UQ(otu_col), UQ(count_col))

  # Rebuild taxonomy table
  taxa <- get_taxa(agg, UQ(otu_col), UQS(tax_ranks))
  taxa <- taxa[match(rownames(taxa), colnames(features)),]

  # Rebuild sample table
  samples <- get_samples(agg, UQ(sample_col), UQ(otu_col), UQ(count_col), UQS(tax_ranks))
  samples <- samples[match(rownames(samples), rownames(features)), ]

  phyloseq::phyloseq(
    phyloseq::otu_table(features, taxa_are_rows = FALSE),
    phyloseq::sample_data(samples),
    phyloseq::tax_table(taxa)
  )
}


