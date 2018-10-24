Variance components for pairs involving Hnf4a
================
Frederick Boehm
10/12/2018

We want to infer variance components for all 139 pairs that involve
*Hnf4a* in the Chromosome 2 hotspot in the Keller et al 2018 data.

``` r
load("../data-to-ignore/Attie_DO378_eQTL_viewer_v1.Rdata")
readRDS("../data/hotspot_expr_tib2_keller_chr2.rds") -> hotspot_expr
readRDS("../data/keller2018-chr2-local-expr.rds") -> local_expr
```

``` r
# what is the ENSEMBL id for Hnf4a?
hnf4a_id <- "ENSMUSG00000017950"
```

``` r
library(qtl2)
library(qtl2pleio)
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.0.0     ✔ purrr   0.2.5
    ## ✔ tibble  1.4.2     ✔ dplyr   0.7.6
    ## ✔ tidyr   0.8.1     ✔ stringr 1.3.1
    ## ✔ readr   1.1.1     ✔ forcats 0.3.0

    ## ── Conflicts ──────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter()   masks stats::filter()
    ## ✖ dplyr::lag()      masks stats::lag()
    ## ✖ readr::read_csv() masks qtl2::read_csv()

``` r
# convert to needed phenotypes input format for scan1
as.matrix(local_expr[, 1:21]) -> local_expr_mat
rownames(local_expr_mat) <- local_expr$mouse_id %>% unlist() %>% as.character()
#
as.matrix(hotspot_expr[, 1:139]) -> hotspot_expr_mat
rownames(hotspot_expr_mat) <- hotspot_expr$mouse_id %>% unlist() %>% as.character()
#
# isolate additive covariates
dataset.islet.rnaseq$covar -> addcovar
addcovar2 <- addcovar[rownames(addcovar) %in% hotspot_expr$mouse_id, ]
```

``` r
# prepare inputs for calc_covs()
k2 <- K$`2`[rownames(K$`2`) %in% rownames(hotspot_expr_mat), colnames(K$`2`) %in% rownames(hotspot_expr_mat)]
```

``` r
get_vc <- function(e2, # n by 1 expression matrix for gene 1
                   e1, # n by 1 expression matrix for gene 2
                   kinship, #kinship matrix
                   covariates = NULL # matrix
                   ){
  pheno <- cbind(e1, e2)
  shared <- intersect(rownames(e1),rownames(e2))
  id2keep <- intersect(shared, rownames(kinship))
  #phea <-subset_input(pheno, id2keep = id2keep)
  #print(phea)
  #subset_kinship(kinship, id2keep) -> ka
  if (!is.null(covariates)) {
    subset_input(covariates, id2keep) -> cova
  }
  else {
    cova <- NULL
  }
  cc_out <- calc_covs(pheno = pheno, kinship = kinship, covariates = cova)
  Vg_vec <- as.vector(cc_out$Vg)
  Ve_vec <- as.vector(cc_out$Ve)
  out <- tibble::as_tibble(t(c(Vg_vec, Ve_vec)))
  out2 <- dplyr::rename(out, Vg11 = V1, Vg12 = V2, Vg21 = V3, Vg22 = V4, Ve11 = V5, Ve12 = V6, Ve21 = V7, Ve22 = V8)
  return(out2)
}
```

``` r
get_vc(e1 = local_expr_mat[ , colnames(local_expr_mat) == hnf4a_id, drop = FALSE],
       e2 = hotspot_expr_mat[ , 1, drop = FALSE],
       kinship = k2,
       covariates = addcovar2[, 1:4]
       )
```

    ## # A tibble: 1 x 8
    ##    Vg11  Vg12  Vg21  Vg22  Ve11  Ve12  Ve21  Ve22
    ##   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 0.962 0.460 0.460  1.78 0.520 0.143 0.143 0.169

``` r
hotspot_expr_mat[, 1:3] %>%
  as_tibble() %>%
  purrr::map_dfr(get_vc,
                 e1 = local_expr_mat[ , colnames(local_expr_mat) == hnf4a_id,
                                      drop = FALSE],
                 kinship = k2,
                 covariates = NULL
                 )
```

    ## # A tibble: 3 x 8
    ##    Vg11  Vg12  Vg21  Vg22  Ve11  Ve12  Ve21  Ve22
    ##   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 0.857 0.466 0.466  1.72 0.613 0.159 0.159 0.198
    ## 2 0.862 0.333 0.333  1.14 0.611 0.126 0.126 0.485
    ## 3 0.851 0.542 0.542  1.32 0.615 0.300 0.300 0.436
