library(readr)
library(dplyr)
library(purrr)
library(doParallel)

library(data.table)

# read all files in imputed_expected_matches into a dataframe d
d <- list.files("../computed_matches", full.names = TRUE) %>%
    map_df(read_csv, col_types = cols()) %>% setDT()

registerDoParallel(cores = 3)

tc_allele_tests <- foreach(i = 1:100) %dopar% {
    wilcox.test(d[d$rep == i, ]$called_allele_matches, d[d$rep == i, ]$true_allele_matches, paired = TRUE, conf.int = TRUE)
}

tc_full_tests <- foreach(i = 1:100) %dopar% {
    wilcox.test(d[d$rep == i, ]$called_full_matches, d[d$rep == i, ]$true_full_matches, paired = TRUE, conf.int = TRUE)
}

tc_partial_tests <- foreach(i = 1:100) %dopar% {
    wilcox.test(d[d$rep == i, ]$called_partial_matches, d[d$rep == i, ]$true_partial_matches, paired = TRUE, conf.int = TRUE)
}

tm_allele_tests <- foreach(i = 1:100) %dopar% {
    wilcox.test(d[d$rep == i, ]$mean_allele_matches, d[d$rep == i, ]$true_allele_matches, paired = TRUE, conf.int = TRUE)
}

tm_full_tests <- foreach(i = 1:100) %dopar% {
    wilcox.test(d[d$rep == i, ]$mean_full_matches, d[d$rep == i, ]$true_full_matches, paired = TRUE, conf.int = TRUE)
}

tm_partial_tests <- foreach(i = 1:100) %dopar% {
    wilcox.test(d[d$rep == i, ]$mean_partial_matches, d[d$rep == i, ]$true_partial_matches, paired = TRUE, conf.int = TRUE)
}

tc_allele_medians <- sapply(tc_allele_tests, function(x) x$estimate)
tc_full_medians <- sapply(tc_full_tests, function(x) x$estimate)
tc_partial_medians <- sapply(tc_partial_tests, function(x) x$estimate)
tm_allele_medians <- sapply(tm_allele_tests, function(x) x$estimate)
tm_full_medians <- sapply(tm_full_tests, function(x) x$estimate)
tm_partial_medians <- sapply(tm_partial_tests, function(x) x$estimate)

tc_allele_pvals <- sapply(tc_allele_tests, function(x) x$p.value)
tc_full_pvals <- sapply(tc_full_tests, function(x) x$p.value)
tc_partial_pvals <- sapply(tc_partial_tests, function(x) x$p.value)
tm_allele_pvals <- sapply(tm_allele_tests, function(x) x$p.value)
tm_full_pvals <- sapply(tm_full_tests, function(x) x$p.value)
tm_partial_pvals <- sapply(tm_partial_tests, function(x) x$p.value)

tc_allele_confint_l <- sapply(tc_allele_tests, function(x) x$conf.int[1])
tc_allele_confint_u <- sapply(tc_allele_tests, function(x) x$conf.int[2])
tc_full_confint_l <- sapply(tc_full_tests, function(x) x$conf.int[1])
tc_full_confint_u <- sapply(tc_full_tests, function(x) x$conf.int[2])
tc_partial_confint_l <- sapply(tc_partial_tests, function(x) x$conf.int[1])
tc_partial_confint_u <- sapply(tc_partial_tests, function(x) x$conf.int[2])
tm_allele_confint_l <- sapply(tm_allele_tests, function(x) x$conf.int[1])
tm_allele_confint_u <- sapply(tm_allele_tests, function(x) x$conf.int[2])
tm_full_confint_l <- sapply(tm_full_tests, function(x) x$conf.int[1])
tm_full_confint_u <- sapply(tm_full_tests, function(x) x$conf.int[2])
tm_partial_confint_l <- sapply(tm_partial_tests, function(x) x$conf.int[1])
tm_partial_confint_u <- sapply(tm_partial_tests, function(x) x$conf.int[2])

wilcoxon_test_results <- data.frame(
    tc_allele_medians, tc_allele_pvals, tc_allele_confint_l, tc_allele_confint_u,
    tc_full_medians, tc_full_pvals, tc_full_confint_l, tc_full_confint_u,
    tc_partial_medians, tc_partial_pvals, tc_partial_confint_l, tc_partial_confint_u,
    tm_allele_medians, tm_allele_pvals, tm_allele_confint_l, tm_allele_confint_u,
    tm_full_medians, tm_full_pvals, tm_full_confint_l, tm_full_confint_u,
    tm_partial_medians, tm_partial_pvals, tm_partial_confint_l, tm_partial_confint_u
)

write_csv(wilcoxon_test_results, "wilcoxon_tests.csv")
