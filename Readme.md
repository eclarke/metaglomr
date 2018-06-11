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
knitr::kable(head(agg))
```

| sample\_id | otu\_id |  count| study\_group | sample\_type | Kingdom  | Phylum        | Class       | Order         | Family             | Genus         | Species |
|:-----------|:--------|------:|:-------------|:-------------|:---------|:--------------|:------------|:--------------|:-------------------|:--------------|:--------|
| S\_1       | f\_1    |     25| case         | A            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_2       | f\_1    |     78| case         | B            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_3       | f\_1    |     34| case         | C            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_4       | f\_1    |     61| case         | D            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_5       | f\_1    |    258| case         | A            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_6       | f\_1    |     54| case         | B            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |

Motivation
----------

In metagenomics, we often work with three distinct datasets: a feature count table, a sample table, and a taxonomic annotation table.

### 1. Feature count table

The feature count table is is a matrix of features × samples, and the cells of the matrix are the times that feature was seen in that sample. Features can be OTUs, ASVs, species, genomes, or whatever.

``` r
data(features)
knitr::kable(features[1:5, 1:15], format = "markdown")
```

|      |  f\_1|  f\_2|  f\_3|  f\_4|  f\_5|  f\_6|  f\_7|  f\_8|  f\_9|  f\_10|  f\_11|  f\_12|  f\_13|  f\_14|  f\_15|
|:-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|------:|------:|------:|------:|------:|------:|
| S\_1 |    25|   112|   130|   123|   101|    61|   220|    72|    43|     19|    105|     97|     52|     76|    258|
| S\_2 |    78|   402|    67|   188|    69|    89|   150|    44|   129|    236|     56|     83|    284|     15|    112|
| S\_3 |    34|   728|   293|   342|     4|     6|    19|   122|   110|     19|    387|      5|    116|     17|    148|
| S\_4 |    61|    84|    41|     5|   106|    39|    97|    51|    97|     58|    322|     84|    132|      1|     35|
| S\_5 |   258|    40|    65|   160|    38|    31|    52|    88|    36|     88|     41|     91|    128|    140|    105|

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
knitr::kable(head(agg))
```

| sample\_id | otu\_id |  count| study\_group | sample\_type | Kingdom  | Phylum        | Class       | Order         | Family             | Genus         | Species |
|:-----------|:--------|------:|:-------------|:-------------|:---------|:--------------|:------------|:--------------|:-------------------|:--------------|:--------|
| S\_1       | f\_1    |     25| case         | A            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_2       | f\_1    |     78| case         | B            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_3       | f\_1    |     34| case         | C            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_4       | f\_1    |     61| case         | D            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_5       | f\_1    |    258| case         | A            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |
| S\_6       | f\_1    |     54| case         | B            | Bacteria | Bacteroidetes | Bacteroidia | Bacteroidales | Porphyromonadaceae | Porphyromonas | NA      |

While this may seem overly repetitive (as the metadata is duplicated in lots of rows), R and dplyr actually handle this pretty well. Things only start breaking down with feature tables that have more than &gt; 100,000,000 cells. What this buys you is the ability to use standard tidyverse verbs and operations easily. Here's how to convert your counts to proportions:

``` r
agg <- agg %>%
  group_by(sample_id) %>%
  mutate(proportion = count/sum(count))

# Showing just a subset of the data:
knitr::kable(head(select(agg, sample_id, otu_id, count, proportion)))
```

| sample\_id | otu\_id |  count|  proportion|
|:-----------|:--------|------:|-----------:|
| S\_1       | f\_1    |     25|   0.0042948|
| S\_2       | f\_1    |     78|   0.0104278|
| S\_3       | f\_1    |     34|   0.0050959|
| S\_4       | f\_1    |     61|   0.0100660|
| S\_5       | f\_1    |    258|   0.0433832|
| S\_6       | f\_1    |     54|   0.0081093|

Aggregate based on taxonomic rank:

``` r
agg %>% 
  group_by(sample_id, Phylum) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  # re-add sample data that got lost in the summarizing
  left_join(get_samples(agg, sample_id, study_group, sample_type))
```

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

Or find the most prevalent phyla in your study groups:

``` r
agg %>% 
  ungroup() %>%
  group_by(study_group, Phylum) %>%
  summarize(mean_proportion = mean(proportion)) %>%
  top_n(1, mean_proportion) %>%
  knitr::kable()
```

| study\_group | Phylum        |  mean\_proportion|
|:-------------|:--------------|-----------------:|
| case         | Firmicutes    |         0.0171578|
| control      | Bacteroidetes |         0.0172942|
