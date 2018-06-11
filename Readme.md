---
title: 'metaglomr: tools for tidy metagenomics'
author: "Erik Clarke"
output: 
  html_document: 
    keep_md: yes
---



This is a package to combine counts, sample metadata, and taxonomic annotations into one long-form tidy dataframe for further analysis. It offers functions to switch between phyloseq objects and these 'agglomerated' dataframes.

## Installation


```r
devtools::install_github("eclarke/metaglomr")
```

## Short usage overview


```r
library(metaglomr)
library(tidyverse)
```

```
## ── Attaching packages ──────────────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 2.2.1     ✔ purrr   0.2.5
## ✔ tibble  1.4.2     ✔ dplyr   0.7.5
## ✔ tidyr   0.8.1     ✔ stringr 1.3.1
## ✔ readr   1.1.1     ✔ forcats 0.3.0
```

```
## ── Conflicts ─────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
data(features)  # A counts matrix, with rows as samples and columns as features
data(samples)   # A dataframe with sample metadata
data(taxa)      # A dataframe or matrix with taxonomic annotations for features

# Combine all of these datasets into one long-form dataframe
agg <- agglomerated(features, samples, taxa, "sample_id")
str(agg)
```

```
## Classes 'tbl_df', 'tbl' and 'data.frame':	1024 obs. of  12 variables:
##  $ sample_id  : chr  "S_1" "S_2" "S_3" "S_4" ...
##  $ otu_id     : chr  "f_1" "f_1" "f_1" "f_1" ...
##  $ count      : num  25 78 34 61 258 54 56 58 48 7 ...
##  $ study_group: Factor w/ 2 levels "case","control": 1 1 1 1 1 1 1 1 2 2 ...
##  $ sample_type: Factor w/ 4 levels "A","B","C","D": 1 2 3 4 1 2 3 4 1 2 ...
##  $ Kingdom    : Factor w/ 1 level "Bacteria": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Phylum     : Factor w/ 5 levels "Actinobacteria",..: 2 2 2 2 2 2 2 2 2 2 ...
##  $ Class      : Factor w/ 10 levels "Actinobacteria",..: 4 4 4 4 4 4 4 4 4 4 ...
##  $ Order      : Factor w/ 20 levels "Actinomycetales",..: 3 3 3 3 3 3 3 3 3 3 ...
##  $ Family     : Factor w/ 35 levels "0319-6G20","Actinomycetaceae",..: 24 24 24 24 24 24 24 24 24 24 ...
##  $ Genus      : Factor w/ 50 levels "Acinetobacter",..: 31 31 31 31 31 31 31 31 31 31 ...
##  $ Species    : Factor w/ 3 levels "intermedia","jejuensis",..: NA NA NA NA NA NA NA NA NA NA ...
```


## Motivation

In metagenomics, we often work with three distinct datasets: a feature count table, a sample table, and a taxonomic annotation table.

### 1. Feature count table
The feature count table is is a matrix of features $\times$ samples, and the cells of the matrix are the times that feature was seen in that sample. 
Features can be OTUs, ASVs, species, genomes, or whatever.


```r
data(features)
knitr::kable(features[1:5, 1:10], format = "markdown")
```



|    | f_1| f_2| f_3| f_4| f_5| f_6| f_7| f_8| f_9| f_10|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|----:|
|S_1 |  25| 112| 130| 123| 101|  61| 220|  72|  43|   19|
|S_2 |  78| 402|  67| 188|  69|  89| 150|  44| 129|  236|
|S_3 |  34| 728| 293| 342|   4|   6|  19| 122| 110|   19|
|S_4 |  61|  84|  41|   5| 106|  39|  97|  51|  97|   58|
|S_5 | 258|  40|  65| 160|  38|  31|  52|  88|  36|   88|

### 2. Sample metadata

The sample metadata can be anything about the samples, but frequently contains at least a description of the sample type and study group of each sample. 
This is frequently referred to as a mapping file in QIIME.


```r
data(samples)
knitr::kable(samples, format="markdown")
```



|sample_id |study_group |sample_type |
|:---------|:-----------|:-----------|
|S_1       |case        |A           |
|S_2       |case        |B           |
|S_3       |case        |C           |
|S_4       |case        |D           |
|S_5       |case        |A           |
|S_6       |case        |B           |
|S_7       |case        |C           |
|S_8       |case        |D           |
|S_9       |control     |A           |
|S_10      |control     |B           |
|S_11      |control     |C           |
|S_12      |control     |D           |
|S_13      |control     |A           |
|S_14      |control     |B           |
|S_15      |control     |C           |
|S_16      |control     |D           |

### 3. Taxonomic annotations
This contains taxonomic annotations for all the features (OTUs or otherwise), split by rank (i.e. Kingdom, Phylum, etc).


```r
data(taxa)
knitr::kable(taxa[1:5, 1:5], format="markdown")
```



|    |Kingdom  |Phylum         |Class               |Order            |Family             |
|:---|:--------|:--------------|:-------------------|:----------------|:------------------|
|f_1 |Bacteria |Bacteroidetes  |Bacteroidia         |Bacteroidales    |Porphyromonadaceae |
|f_2 |Bacteria |Proteobacteria |Alphaproteobacteria |Rhizobiales      |Xanthobacteraceae  |
|f_3 |Bacteria |Firmicutes     |Negativicutes       |Selenomonadales  |Veillonellaceae    |
|f_4 |Bacteria |Firmicutes     |Clostridia          |Clostridiales    |Defluviitaleaceae  |
|f_5 |Bacteria |Bacteroidetes  |Bacteroidia         |Flavobacteriales |Weeksellaceae      |

### Agglomerating the three datasets

For users of the [tidyverse](https://tidyverse.org), it's frequently easiest to work with long-form
melted datasets, where each row is a unique observation or data point. 
The unique reference in these three datasets is the count of a feature in a particular sample. 
Therefore, we can create a dataframe where each row is this unique combination of feature + sample + count, with additional columns describing the sample and feature further.


```r
agg <- agglomerated(features, samples, taxa, "sample_id")
print(agg)
```

```
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

While this may seem overly repetitive (as the metadata is duplicated in lots of rows), R and dplyr actually handle this pretty well. 
Things only start breaking down with feature tables that have more than > 100,000,000 cells. 
What this buys you is the ability to use standard tidyverse verbs and operations easily. Here's how to convert your counts to proportions:


```r
agg <- agg %>%
  group_by(sample_id) %>%
  mutate(proportion = count/sum(count))
select(agg, sample_id, otu_id, count, proportion)
```

```
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

Or easily find the most prevalent phyla in your study groups:


```r
agg %>% 
  ungroup() %>%
  group_by(study_group, Phylum) %>%
  summarize(mean_proportion = mean(proportion)) %>%
  top_n(1, mean_proportion)
```

```
## # A tibble: 2 x 3
## # Groups:   study_group [2]
##   study_group Phylum        mean_proportion
##   <fct>       <fct>                   <dbl>
## 1 case        Firmicutes             0.0172
## 2 control     Bacteroidetes          0.0173
```

