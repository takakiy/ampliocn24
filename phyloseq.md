# ANALYSIS OF PHYLOSEQ


## 準備

```
$ R

library(tidyverse)

library(phyloseq)
library(MicrobeR)
library(ggplot2)
library(ggsci)

```

## Qiime2からphyloseqクラスの作成

```
library(qiime2R)

 otus <- read_qza("table-dada2-nochim.qza")
 tree <- read_qza("rooted-tree.qza")
 taxonomy <- read_qza("taxonomy.qza")
 tax_table <- do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), ";"))

 colnames(tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
 rownames(tax_table) <- taxonomy$data$Feature.ID
 metadata <- read.table("sample-metadata.tsv", sep='\t', header=T, row.names=1, comment="")
 physeq <- phyloseq(
    otu_table(otus$data, taxa_are_rows = T),
    phy_tree(tree$data),
    tax_table(tax_table),
    sample_data(metadata)
    )

```

## 1. 細菌占有率グラフの作成(Familyレベル)


## Qiime2からphyloseqクラスの作成


## Qiime2からphyloseqクラスの作成


## Qiime2からphyloseqクラスの作成


## Qiime2からphyloseqクラスの作成


## Qiime2からphyloseqクラスの作成


