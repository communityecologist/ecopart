# Partitioning the temporal changes in beta diversity
The R package `ecopart` ("Ecological COmmunity PARTitioning" or "Extinction and COlonization PARTitioning") allows partitioning of the temporal changes in beta diversity into multiple components. The components represent biotic homogenization and differentiation (i.e., decrease and increase in beta diversity) that result from local extinctions and colonizations of species or the losses and gains in species abundances. The pacakge consists of two functions: `ecopart.pair()` and `ecopart.multi()`.

## Installation
```{r}
library(remotes)
remotes::install_github(“communityecologist/ecopart”)
```

## Usage
###### ecopart.pair(d1, d2, index="ruzicka", components="two")
- `d1` : Community matrix at time 1, where rows are sites, columns are species, and elements are species presence-abacence (01) or abundances.
- `d2` : Community matrix at time 2, where rows are sites, columns are species, and elements are species presence-abacence (01) or abundances.
- `index` : Type of beta diversity measure. Partial match of "ruzicka, "bray-curtis", or "harrison-baselga".
- `components` : Types of components into which Δβ is partitioned. Partial match of "two", "four", "six", or "sp".

###### ecopart.multi(d1, d2, components="two")
- `d1` : Community matrix at time 1, where rows are sites, columns are species, and elements are species presence-abacence (01) or abundances.
- `d2` : Community matrix at time 2, where rows are sites, columns are species, and elements are species presence-abacence (01) or abundances.
- `components` : Types of components into which Δβ is partitioned. Partial match of "two", "four", "six", or "sp".

See Appendix 3 in Tatsumi et al. (2021) *Ecology Letters* for detail.

## Citations
* [Tatsumi S, Iritani R, Cadotte MW (2021) Temporal changes in spatial variation: partitioning the extinction and colonisation components of beta diversity. *Ecology Letters* 24(5): 1063–1072.](https://onlinelibrary.wiley.com/doi/10.1111/ele.13720)
* [Tatsumi S, Iritani R, Cadotte MW (2021) Partitioning the temporal changes in abundance-based beta diversity into loss and gain components. *bioRxiv*, 2021.05.28.446095.](https://www.biorxiv.org/content/10.1101/2021.05.28.446095v1)
