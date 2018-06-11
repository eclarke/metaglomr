metaglomr: tools for tidy metagenomics
================
Erik Clarke

This is a package to combine counts, sample metadata, and taxonomic annotations into one long-form tidy dataframe for further analysis. It offers functions to switch between phyloseq objects and these 'agglomerated' dataframes.

[Read the documentation.](https://eclarke.github.io/metaglomr)

Installation
------------

``` r
devtools::install_github("eclarke/metaglomr")
```

Short usage overview
--------------------

``` r
suppressPackageStartupMessages(library(tidyverse))
library(metaglomr)

data(features)  # A counts matrix, with rows as samples and columns as features
data(samples)   # A dataframe with sample metadata
data(taxa)      # A dataframe or matrix with taxonomic annotations for features

# Combine all of these datasets into one long-form dataframe
(agg <- agglomerated(features, samples, taxa, "sample_id"))
## # A tibble: 1,024 x 12
##    sample_id otu_id count study_group sample_type Kingdom  Phylum   Class 
##    <chr>     <chr>  <dbl> <fct>       <fct>       <fct>    <fct>    <fct> 
##  1 S_1       f_1       25 case        A           Bacteria Bactero… Bacte…
##  2 S_2       f_1       78 case        B           Bacteria Bactero… Bacte…
##  3 S_3       f_1       34 case        C           Bacteria Bactero… Bacte…
##  4 S_4       f_1       61 case        D           Bacteria Bactero… Bacte…
##  5 S_5       f_1      258 case        A           Bacteria Bactero… Bacte…
##  6 S_6       f_1       54 case        B           Bacteria Bactero… Bacte…
##  7 S_7       f_1       56 case        C           Bacteria Bactero… Bacte…
##  8 S_8       f_1       58 case        D           Bacteria Bactero… Bacte…
##  9 S_9       f_1       48 control     A           Bacteria Bactero… Bacte…
## 10 S_10      f_1        7 control     B           Bacteria Bactero… Bacte…
## # ... with 1,014 more rows, and 4 more variables: Order <fct>,
## #   Family <fct>, Genus <fct>, Species <fct>
```
