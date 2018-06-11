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
agg <- agglomerated(features, samples, taxa, "sample_id")
str(agg)
```

Classes 'tbl\_df', 'tbl' and 'data.frame': 1024 obs. of 12 variables: $ sample\_id : chr "S\_1" "S\_2" "S\_3" "S\_4" ... $ otu\_id : chr "f\_1" "f\_1" "f\_1" "f\_1" ... $ count : num 25 78 34 61 258 54 56 58 48 7 ... $ study\_group: Factor w/ 2 levels "case","control": 1 1 1 1 1 1 1 1 2 2 ... $ sample\_type: Factor w/ 4 levels "A","B","C","D": 1 2 3 4 1 2 3 4 1 2 ... $ Kingdom : Factor w/ 1 level "Bacteria": 1 1 1 1 1 1 1 1 1 1 ... $ Phylum : Factor w/ 5 levels "Actinobacteria",..: 2 2 2 2 2 2 2 2 2 2 ... $ Class : Factor w/ 10 levels "Actinobacteria",..: 4 4 4 4 4 4 4 4 4 4 ... $ Order : Factor w/ 20 levels "Actinomycetales",..: 3 3 3 3 3 3 3 3 3 3 ... $ Family : Factor w/ 35 levels "0319-6G20","Actinomycetaceae",..: 24 24 24 24 24 24 24 24 24 24 ... $ Genus : Factor w/ 50 levels "Acinetobacter",..: 31 31 31 31 31 31 31 31 31 31 ... $ Species : Factor w/ 3 levels "intermedia","jejuensis",..: NA NA NA NA NA NA NA NA NA NA ...

Motivation
----------

In metagenomics, we often work with three distinct datasets: a feature count table, a sample table, and a taxonomic annotation table.

### 1. Feature count table

The feature count table is is a matrix of features × samples, and the cells of the matrix are the times that feature was seen in that sample. Features can be OTUs, ASVs, species, genomes, or whatever.

``` r
data(features)
knitr::kable(features[1:5, 1:10], format = "markdown")
```

|      |  f\_1|  f\_2|  f\_3|  f\_4|  f\_5|  f\_6|  f\_7|  f\_8|  f\_9|  f\_10|
|:-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|------:|
| S\_1 |    25|   112|   130|   123|   101|    61|   220|    72|    43|     19|
| S\_2 |    78|   402|    67|   188|    69|    89|   150|    44|   129|    236|
| S\_3 |    34|   728|   293|   342|     4|     6|    19|   122|   110|     19|
| S\_4 |    61|    84|    41|     5|   106|    39|    97|    51|    97|     58|
| S\_5 |   258|    40|    65|   160|    38|    31|    52|    88|    36|     88|

### 2. Sample metadata

The sample metadata can be anything about the samples, but frequently contains at least a description of the sample type and study group of each sample. This is frequently referred to as a mapping file in QIIME.

``` r
data(samples)
knitr::kable(samples, format="markdown")
```

| sample\_id | study\_group | sample\_type |
|:-----------|:-------------|:-------------|
| S\_1       | case         | A            |
| S\_2       | case         | B            |
| S\_3       | case         | C            |
| S\_4       | case         | D            |
| S\_5       | case         | A            |
| S\_6       | case         | B            |
| S\_7       | case         | C            |
| S\_8       | case         | D            |
| S\_9       | control      | A            |
| S\_10      | control      | B            |
| S\_11      | control      | C            |
| S\_12      | control      | D            |
| S\_13      | control      | A            |
| S\_14      | control      | B            |
| S\_15      | control      | C            |
| S\_16      | control      | D            |

### 3. Taxonomic annotations

This contains taxonomic annotations for all the features (OTUs or otherwise), split by rank (i.e. Kingdom, Phylum, etc).

``` r
data(taxa)
knitr::kable(taxa[1:5, 1:5], format="markdown")
```

|      | Kingdom  | Phylum         | Class               | Order            | Family             |
|:-----|:---------|:---------------|:--------------------|:-----------------|:-------------------|
| f\_1 | Bacteria | Bacteroidetes  | Bacteroidia         | Bacteroidales    | Porphyromonadaceae |
| f\_2 | Bacteria | Proteobacteria | Alphaproteobacteria | Rhizobiales      | Xanthobacteraceae  |
| f\_3 | Bacteria | Firmicutes     | Negativicutes       | Selenomonadales  | Veillonellaceae    |
| f\_4 | Bacteria | Firmicutes     | Clostridia          | Clostridiales    | Defluviitaleaceae  |
| f\_5 | Bacteria | Bacteroidetes  | Bacteroidia         | Flavobacteriales | Weeksellaceae      |

### Agglomerating the three datasets

For users of the [tidyverse](https://tidyverse.org), it's frequently easiest to work with long-form melted datasets, where each row is a unique observation or data point. The unique reference in these three datasets is the count of a feature in a particular sample. Therefore, we can create a dataframe where each row is this unique combination of feature + sample + count, with additional columns describing the sample and feature further.

``` r
agg <- agglomerated(features, samples, taxa, "sample_id")
print(agg)
```

A tibble: 1,024 x 12
====================

sample\_id otu\_id count study\_group sample\_type Kingdom Phylum Class <chr> <chr> <dbl> <fct> <fct> <fct> <fct> <fct> 1 S\_1 f\_1 25 case A Bacteria Bactero… Bacte… 2 S\_2 f\_1 78 case B Bacteria Bactero… Bacte… 3 S\_3 f\_1 34 case C Bacteria Bactero… Bacte… 4 S\_4 f\_1 61 case D Bacteria Bactero… Bacte… 5 S\_5 f\_1 258 case A Bacteria Bactero… Bacte… 6 S\_6 f\_1 54 case B Bacteria Bactero… Bacte… 7 S\_7 f\_1 56 case C Bacteria Bactero… Bacte… 8 S\_8 f\_1 58 case D Bacteria Bactero… Bacte… 9 S\_9 f\_1 48 control A Bacteria Bactero… Bacte… 10 S\_10 f\_1 7 control B Bacteria Bactero… Bacte… \# ... with 1,014 more rows, and 4 more variables: Order <fct>, \# Family <fct>, Genus <fct>, Species <fct>

While this may seem overly repetitive (as the metadata is duplicated in lots of rows), R and dplyr actually handle this pretty well. Things only start breaking down with feature tables that have more than &gt; 100,000,000 cells. What this buys you is the ability to use standard tidyverse verbs and operations easily. Here's how to convert your counts to proportions:

``` r
agg <- agg %>%
  group_by(sample_id) %>%
  mutate(proportion = count/sum(count))
select(agg, sample_id, otu_id, count, proportion)
```

A tibble: 1,024 x 4
===================

Groups: sample\_id \[16\]
=========================

sample\_id otu\_id count proportion <chr> <chr> <dbl> <dbl> 1 S\_1 f\_1 25 0.00429 2 S\_2 f\_1 78 0.0104 3 S\_3 f\_1 34 0.00510 4 S\_4 f\_1 61 0.0101 5 S\_5 f\_1 258 0.0434 6 S\_6 f\_1 54 0.00811 7 S\_7 f\_1 56 0.00874 8 S\_8 f\_1 58 0.0120 9 S\_9 f\_1 48 0.00635 10 S\_10 f\_1 7 0.00116 \# ... with 1,014 more rows

Or easily find the most prevalent phyla in your study groups:

``` r
agg %>% 
  ungroup() %>%
  group_by(study_group, Phylum) %>%
  summarize(mean_proportion = mean(proportion)) %>%
  top_n(1, mean_proportion)
```

A tibble: 2 x 3
===============

Groups: study\_group \[2\]
==========================

study\_group Phylum mean\_proportion <fct> <fct> <dbl> 1 case Firmicutes 0.0172 2 control Bacteroidetes 0.0173
