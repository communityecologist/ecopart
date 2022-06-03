# Partitioning the temporal changes in beta diversity
The R package `ecopart` ("Ecological COmmunity PARTitioning" or "Extinction and COlonization PARTitioning") allows one to partition the temporal changes in beta diversity into multiple components. The components represent biotic homogenization and differentiation (i.e., decreases and increases in beta diversity) that result from local extinctions and colonizations of species or the losses and gains in species abundances. The pacakge consists of two functions: `ecopart.pair()` and `ecopart.multi()`.

## Installation
```{r}
#library(remotes)
remotes::install_github("communityecologist/ecopart")
library(ecopart)
```

## Usage
```{r}
ecopart.pair(d1, d2, index = "sorensen", components="four")
```
- `d1` : A matrix or dataframe at time 1. Rows are a pair of sites (sites 1 and 2), columns are species, and elements are presence-absence (01) or abundances of species.
- `d2` : A matrix or dataframe at time 2. Note that d1 and d2 must have exactly the same sites and species in the same order.
- `index` : Type of dissimilarity measure. Options are "jaccard", "sorensen", "ruzicka", and "bray-curtis".
- `components` : Types of components into which the total change in beta diversity is partitioned. Options are "two", "four", "six", and "sp".

```{r}
ecopart.multi(d1, d2, index = "whittaker", components="four")
```
- `d1` : A matrix or dataframe at time 1. Rows are sites, columns are species, and elements are presence-absence (01) or abundances of species.
- `d2` : A matrix or dataframe at time 2. Note that d1 and d2 must have exactly the same sites and species in the same order.
- `index` : Type of dissimilarity measure. Options are "whittaker" and "baselga".
- `components` : Types of components into which the total change in beta diversity is partitioned. Options are "two", "four", and "sp".

Run `?ecopart.pair()` and `?ecopart.multi()` for detail.

## Citations
* [Tatsumi S, Iritani R, Cadotte MW (2021) Temporal changes in spatial variation: partitioning the extinction and colonisation components of beta diversity. *Ecology Letters* 24(5): 1063–1072.](https://onlinelibrary.wiley.com/doi/10.1111/ele.13720)
* Tatsumi S, Iritani R, Cadotte MW (2022) Partitioning the temporal changes in abundance-based beta diversity into loss and gain components. *Methods in Ecology and Evolution*, in press.
