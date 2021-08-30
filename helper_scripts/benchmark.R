# benchmark

library(microbenchmark)
library(tidyverse)
library(Matrix)
library(Rcpp)


############### SET UP TEST EXAMPLE ###################
nrow_ <- 96
ncol_ <- 546

x <- matrix(runif(nrow_ * ncol_), nrow = nrow_, ncol = ncol_)
x2 <- apply(x, 1, function(r) {
    r[sample(c(1:ncol_), floor(ncol_ / 10))] <- NA
    return(r)
}) %>% t()


c_full <- cor(x2, use = "pairwise.complete.obs", method = "spearman")
c <- c_full[seq_len(12), seq_len(12)]

names_ <- sapply(
    c("A", "B", "C", "D"),
    function(chr) str_c(chr, c("x","y","z"), sep = "--")
) %>% as.character()
colnames(c) <- names_
rownames(c) <- names_
############### SET UP TEST EXAmPLE ###################

get_grp_names_from_matrix <- function(M, unique_cs, sep_ = "--") {
    grp_names <- str_split(unique_cs, pattern = sep_, simplify = T)[, 1]
    return(grp_names)
}

