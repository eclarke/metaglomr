metaglomr: tools for tidy metagenomics
================
Erik Clarke

This is a package to combine counts, sample metadata, and taxonomic annotations into one long-form tidy dataframe for further analysis. It offers functions to switch between phyloseq objects and these 'agglomerated' dataframes.

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

Motivation
----------

In metagenomics, we often work with three distinct datasets: a feature count table, a sample table, and a taxonomic annotation table.

### 1. Feature count table

The feature count table is is a matrix of features × samples, and the cells of the matrix are the times that feature was seen in that sample. Features can be OTUs, ASVs, species, genomes, or whatever.

``` r
data(features)
features[1:5, 1:15]
##     f_1 f_2 f_3 f_4 f_5 f_6 f_7 f_8 f_9 f_10 f_11 f_12 f_13 f_14 f_15
## S_1  25 112 130 123 101  61 220  72  43   19  105   97   52   76  258
## S_2  78 402  67 188  69  89 150  44 129  236   56   83  284   15  112
## S_3  34 728 293 342   4   6  19 122 110   19  387    5  116   17  148
## S_4  61  84  41   5 106  39  97  51  97   58  322   84  132    1   35
## S_5 258  40  65 160  38  31  52  88  36   88   41   91  128  140  105
```

### 2. Sample metadata

The sample metadata can be anything about the samples, but frequently contains at least a description of the sample type and study group of each sample. This is frequently referred to as a mapping file in QIIME.

``` r
data(samples)
as_tibble(samples)
## # A tibble: 16 x 3
##    sample_id study_group sample_type
##    <fct>     <fct>       <fct>      
##  1 S_1       case        A          
##  2 S_2       case        B          
##  3 S_3       case        C          
##  4 S_4       case        D          
##  5 S_5       case        A          
##  6 S_6       case        B          
##  7 S_7       case        C          
##  8 S_8       case        D          
##  9 S_9       control     A          
## 10 S_10      control     B          
## 11 S_11      control     C          
## 12 S_12      control     D          
## 13 S_13      control     A          
## 14 S_14      control     B          
## 15 S_15      control     C          
## 16 S_16      control     D
```

### 3. Taxonomic annotations

This contains taxonomic annotations for all the features (OTUs or otherwise), split by rank (i.e. Kingdom, Phylum, etc).

``` r
data(taxa)
taxa[1:5, ]
##      Kingdom         Phylum               Class            Order
## f_1 Bacteria  Bacteroidetes         Bacteroidia    Bacteroidales
## f_2 Bacteria Proteobacteria Alphaproteobacteria      Rhizobiales
## f_3 Bacteria     Firmicutes       Negativicutes  Selenomonadales
## f_4 Bacteria     Firmicutes          Clostridia    Clostridiales
## f_5 Bacteria  Bacteroidetes         Bacteroidia Flavobacteriales
##                 Family                     Genus Species
## f_1 Porphyromonadaceae             Porphyromonas    <NA>
## f_2  Xanthobacteraceae              Pseudolabrys    <NA>
## f_3    Veillonellaceae               Veillonella    <NA>
## f_4  Defluviitaleaceae Defluviitaleaceae_UCG-011    <NA>
## f_5      Weeksellaceae          Chryseobacterium    <NA>
```

### Agglomerating the three datasets

For users of the [tidyverse](https://tidyverse.org), it's frequently easiest to work with long-form melted datasets, where each row is a unique observation or data point. The unique reference in these three datasets is the count of a feature in a particular sample. Therefore, we can create a dataframe where each row is this unique combination of feature + sample + count, with additional columns describing the sample and feature further.

``` r
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

While this may seem overly repetitive (as the metadata is duplicated in lots of rows), R and dplyr actually handle this pretty well. Things only start breaking down with feature tables that have more than &gt; 100,000,000 cells. What this buys you is the ability to use standard tidyverse verbs and operations easily.

Examples
--------

### Subsetting

Subsetting is easy through the use of the `filter` verb:

``` r
filter(agg, study_group == "case")
## # A tibble: 512 x 12
##    sample_id otu_id count study_group sample_type Kingdom  Phylum  Class  
##    <chr>     <chr>  <dbl> <fct>       <fct>       <fct>    <fct>   <fct>  
##  1 S_1       f_1       25 case        A           Bacteria Bacter… Bacter…
##  2 S_2       f_1       78 case        B           Bacteria Bacter… Bacter…
##  3 S_3       f_1       34 case        C           Bacteria Bacter… Bacter…
##  4 S_4       f_1       61 case        D           Bacteria Bacter… Bacter…
##  5 S_5       f_1      258 case        A           Bacteria Bacter… Bacter…
##  6 S_6       f_1       54 case        B           Bacteria Bacter… Bacter…
##  7 S_7       f_1       56 case        C           Bacteria Bacter… Bacter…
##  8 S_8       f_1       58 case        D           Bacteria Bacter… Bacter…
##  9 S_1       f_2      112 case        A           Bacteria Proteo… Alphap…
## 10 S_2       f_2      402 case        B           Bacteria Proteo… Alphap…
## # ... with 502 more rows, and 4 more variables: Order <fct>, Family <fct>,
## #   Genus <fct>, Species <fct>
filter(agg, Phylum == "Bacteroidetes")
## # A tibble: 176 x 12
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
## # ... with 166 more rows, and 4 more variables: Order <fct>, Family <fct>,
## #   Genus <fct>, Species <fct>
```

Here's how to convert your counts to proportions:

``` r
agg <- agg %>%
  group_by(sample_id) %>%
  mutate(proportion = count/sum(count))

# Showing just a subset of the data:
select(agg, sample_id, otu_id, count, proportion)
## # A tibble: 1,024 x 4
## # Groups:   sample_id [16]
##    sample_id otu_id count proportion
##    <chr>     <chr>  <dbl>      <dbl>
##  1 S_1       f_1       25    0.00429
##  2 S_2       f_1       78    0.0104 
##  3 S_3       f_1       34    0.00510
##  4 S_4       f_1       61    0.0101 
##  5 S_5       f_1      258    0.0434 
##  6 S_6       f_1       54    0.00811
##  7 S_7       f_1       56    0.00874
##  8 S_8       f_1       58    0.0120 
##  9 S_9       f_1       48    0.00635
## 10 S_10      f_1        7    0.00116
## # ... with 1,014 more rows
```

Aggregate based on taxonomic rank:

``` r
agg %>% 
  group_by(sample_id, Phylum) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  # re-add sample data that got lost in the summarizing
  left_join(get_samples(agg, sample_id, study_group, sample_type))
## Joining, by = "sample_id"
## # A tibble: 80 x 5
##    sample_id Phylum         count study_group sample_type
##    <chr>     <fct>          <dbl> <fct>       <fct>      
##  1 S_1       Actinobacteria   363 case        A          
##  2 S_1       Bacteroidetes   1007 case        A          
##  3 S_1       Firmicutes      3616 case        A          
##  4 S_1       Proteobacteria   699 case        A          
##  5 S_1       Spirochaetes     136 case        A          
##  6 S_10      Actinobacteria   627 control     B          
##  7 S_10      Bacteroidetes   1169 control     B          
##  8 S_10      Firmicutes      2869 control     B          
##  9 S_10      Proteobacteria  1043 control     B          
## 10 S_10      Spirochaetes     323 control     B          
## # ... with 70 more rows
```

Or find the most prevalent phyla in your study groups:

``` r
agg %>% 
  ungroup() %>%
  group_by(study_group, Phylum) %>%
  summarize(mean_proportion = mean(proportion)) %>%
  top_n(1, mean_proportion)
## # A tibble: 2 x 3
## # Groups:   study_group [2]
##   study_group Phylum        mean_proportion
##   <fct>       <fct>                   <dbl>
## 1 case        Firmicutes             0.0172
## 2 control     Bacteroidetes          0.0173
```
