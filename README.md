# partitions
Performs entity resolution for data bases with categorical fields using partition-based Bayesian clustering models. Includes two new microclustering prior models for random partitions, and the traditional Dirichlet and Pitman-Yor process priors.

## Installation

```r
devtools::install_github("bbetanc1/partitions", build_vignettes = TRUE)
```

## Citation

This package implements the methods introduced in the following paper:

> Betancourt, Brenda, Giacomo Zanella and Rebecca C. Steorts  (2020). "Random Partitions Models for Microclustering tasks".


## Background

Entity resolution (record linkage or de-duplication) is used to join multiple databases to remove duplicate entities. Recent methods tackle the entity resolution problem as a clustering task. While traditional Bayesian random partition models assume that the size of each cluster grows linearly with the number of data points, this assumption is not appropriate for applications such as entity resolution. This problem requires
models that yield clusters whose sizes grow sublinearly with the total number of data
points -- `the microclustering property`. The `partitions` package includes two random partition models that satisfy the microclustering property and implements entity resolution with categorical data.

## Main functions

The package contains the implemetation of four random partition models that are used to perform entitiy resolution as a clustering task: 

* Two traditional random partition models: Dirichlet process (DP) mixtures and Pitman–Yor process (PY) mixtures. 
* Two random partition models that exhibit the microclustering property, referred to as: The ESCNB model and the ESCD model.

The two main functions in `partitions` are `SampleCluster` and `SampleClusterESCD`, which perform entity resolution for the four models. Additionally, we have added the function `calcError` to evaluate the record linkage performance when ground truth is available.

```r
help(partitions)
```

For more extensive documentation of the use of this package, please see the vignette.

```r
vignette("partitions")
```

## Acknowledgements

This work was partially supported by the National Science Foundation through
