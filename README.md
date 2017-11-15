
<!-- README.md is generated from README.Rmd. Please edit that file -->
boclust
=======

The goal of boclust is to provide a new normalization method for sparse data by a feature boosting strategy with the latent representation, especially for scRNA-seq data consisted of many zeros. Based on the normalization, a new measure of similarity is defined for the following clustering algorithm. Unlike other unsupervised cluster methods, boclust provides the suggestion K to determine the number of clusters. In this way, it may be unsuitable for low-dimentional data.

There are three major functions:

-   BossaSimi: to calculate the similarity matrix and normalized data.
-   BossaClust: the main function which provide an object including clustering result for shiny.
-   bossa\_interactive: a shiny framework to show the clustering result.

Installation
------------

You can install boclust from github with:

``` r
# install.packages("devtools")
devtools::install_github("TinyOpen/boclust")
```

Example
-------

``` r
# generate sparse data from the toy model of CIDR
sparse.data <- data.frame(g.1 = c(0, 5, 0, 6, 8, 6, 7, 7), 
                          g.2 = c(5, 0, 0, 0, 5, 7, 5, 7)) 
bossa.change <- BossaSimi(sparse.data, is.pca = FALSE) # with low-dimensional data, pca is uncessary
data.after <- bossa.change$U.score.non.pca # data after normalization
```

You can check after normalization, the first 4 cells which are actually from the same cluster 
are more closer. The seperation between the first 4 cells and the last 4 cells is large enough to 
get the correct clustering result.

``` r
d3heatmap(sparse.data) ## show heatmap of original data
d3heatmap(data.after) ## show heatmap of bossa-normalized data 
```

Now, when it comes to your **high-dimentional** data, which is the target which `boclust` is designed for. You can either use `BossaClust` to get the final result:

``` r
object <- BossaClust(high.dim.data) # do normalization and clustering at the same time
bossa_interactive(object) # use shiny frame to show the result
```

Or, you can store the normalized data first, which is obtained from function `BossaSimi`, and then do the rest work.

``` r
pre.object <- BossaSimi(high.dim.data)
object <- BossaClust(data = high.dim.data, data.pre = pre.object) # do normalization and clustering at the same time
bossa_interactive(object) # use shiny frame to show the result
```
