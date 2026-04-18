# RCMap
Group Concept Mapping in R

RCMap is an open-source implementation of group concept mapping, implemented in R. It provides a menu-driven user interface to guide users through the concept mapping analytical process.

## Requirements

- R >= 4.0.0
- Dependencies: `smacof`, `crayon`, `ape`, `tcltk`, `stringdist`, `cluster`

## Installation

```r
# install.packages("remotes")
remotes::install_github("haimbar/RCMap")
```

## Quick Start

```r
library(RCMap)
RCMapMenu()
```

This launches the interactive menu-driven interface for the full concept mapping workflow.

## Documentation

For detailed documentation and examples, see the [Documentation](https://haimbar.github.io/RCMap/) page (web-only resource, not a package vignette).

**Data**: The repository contains two sample datasets: *sampleproject* contains a synthetic dataset with 10 sorters/raters. *DataNCATS* contains a real dataset with 74 (anonymized) sorters/raters.