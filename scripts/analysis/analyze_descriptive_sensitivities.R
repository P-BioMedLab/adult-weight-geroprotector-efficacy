suppressPackageStartupMessages(library(survival))
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
outdir <- file.path(package, "data", "outputs")

km_median <- function(z) {
  fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z, conf.type = "none")
  as.numeric(summary(fit)$table["median"])
}

empirical_bins <- function(v, k) {
  cut(v, breaks = quantile(v, probs = 0:k / k, type = 7),
      include.lowest = TRUE, labels = FALSE)
}

rank_bins <- function(v, k) {
  cut(rank(v, ties.method = "average") / length(v),
      breaks = seq(0, 1, length.out = k + 1),
      include.lowest = TRUE, labels = FALSE)
}

summarize_bins <- function(d, analysis, k, rank_based = FALSE) {
  d$bin <- ave(d$body_weight, d$site, FUN = function(v) {
    if (rank_based) rank_bins(v, k) else empirical_bins(v, k)
  })
  rows <- lapply(seq_len(k), function(q) {
    z <- d[d$bin == q, ]; c <- z[z$Tx == 0, ]; t <- z[z$Tx == 1, ]
    mc <- km_median(c); mt <- km_median(t)
    data.frame(analysis = analysis, bins = k, bin = q,
               mean_control_weight = mean(c$body_weight),
               control_median = mc, treated_median = mt,
               gain_pct = 100 * (mt - mc) / mc)
  })
  do.call(rbind, rows)
}

lower <- read.csv(file.path(outdir, "predosing_weight_lowering_27arm_comparisons.csv"),
                  stringsAsFactors = FALSE)
lower <- lower[!duplicated(lower$mouse), ]
nonlower <- read.csv(file.path(outdir, "predosing_weight_nonlowering_29arm_unique.csv"),
                     stringsAsFactors = FALSE)

rows <- list(
  summarize_bins(lower, "weight-lowering conventional quartiles", 4, FALSE),
  summarize_bins(nonlower, "non-lowering conventional quartiles", 4, FALSE),
  summarize_bins(lower, "weight-lowering rank tertiles", 3, TRUE),
  summarize_bins(lower, "weight-lowering rank quartiles", 4, TRUE),
  summarize_bins(lower, "weight-lowering rank fifths", 5, TRUE)
)
out <- do.call(rbind, rows)
write.csv(out, file.path(outdir, "descriptive_sensitivities_reproduced.csv"), row.names = FALSE)
print(out, row.names = FALSE, digits = 6)
