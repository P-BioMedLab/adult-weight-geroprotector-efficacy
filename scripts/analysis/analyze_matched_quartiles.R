suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
input <- if (length(args) >= 1) args[1] else file.path(package, "data", "outputs", "predosing_weight_lowering_27arm_comparisons.csv")
outdir <- if (length(args) >= 2) args[2] else file.path(package, "data", "outputs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input, stringsAsFactors = FALSE)
d$Tx <- as.integer(d$Tx)
d$dead <- as.integer(d$dead)
d$matched_quartile <- as.integer(d$matched_quartile)
d$q_factor <- factor(d$matched_quartile, levels = 1:4, labels = paste0("Q", 1:4))
d$arm <- factor(d$arm)
d$stratum <- factor(d$stratum)
d$mouse <- factor(d$mouse)

weighted_median <- function(z) {
  fit <- survfit(
    Surv(entry, lifespan_days, dead) ~ 1,
    data = z,
    weights = standardization_weight,
    conf.type = "none"
  )
  as.numeric(summary(fit)$table["median"])
}

quartile_rows <- lapply(1:4, function(q) {
  z <- d[d$matched_quartile == q, ]
  control <- z[z$Tx == 0, ]
  treated <- z[z$Tx == 1, ]
  mc <- weighted_median(control)
  mt <- weighted_median(treated)
  fit <- coxph(
    Surv(entry, lifespan_days, dead) ~ Tx + strata(stratum) + cluster(mouse),
    data = z,
    ties = "efron"
  )
  s <- summary(fit)$coefficients["Tx", ]
  data.frame(
    quartile = q,
    comparison_rows = nrow(z),
    unique_mice = length(unique(z$mouse)),
    mean_control_weight = weighted.mean(control$body_weight, control$standardization_weight),
    standardized_control_median = mc,
    standardized_treated_median = mt,
    gain_days = mt - mc,
    gain_pct = 100 * (mt - mc) / mc,
    treatment_HR = exp(s["coef"]),
    treatment_lo = exp(s["coef"] - 1.96 * s["robust se"]),
    treatment_hi = exp(s["coef"] + 1.96 * s["robust se"]),
    treatment_p = s["Pr(>|z|)"]
  )
})
quartiles <- do.call(rbind, quartile_rows)

# Shared ordinal treatment-by-quartile trend while allowing each randomized
# arm to retain its own treatment main effect. Quartiles were assigned within
# arm x site x treatment, so every comparison contributes across the gradient.
d$q_centered <- d$matched_quartile - 1
d$tx_q <- d$Tx * d$q_centered
trend_fit <- coxph(
  Surv(entry, lifespan_days, dead) ~ 0 + arm:Tx + q_factor + tx_q +
    strata(stratum) + cluster(mouse),
  data = d,
  ties = "efron",
  singular.ok = TRUE
)
trend_s <- summary(trend_fit)$coefficients["tx_q", ]
trend <- data.frame(
  analysis = "within-arm-site matched quartile trend",
  arms = length(unique(d$arm)),
  comparison_rows = nrow(d),
  unique_mice = length(unique(d$mouse)),
  HR_per_quartile = exp(trend_s["coef"]),
  lo = exp(trend_s["coef"] - 1.96 * trend_s["robust se"]),
  hi = exp(trend_s["coef"] + 1.96 * trend_s["robust se"]),
  p = trend_s["Pr(>|z|)"]
)

write.csv(quartiles, file.path(outdir, "matched_quartile_standardized_survival.csv"), row.names = FALSE)
write.csv(trend, file.path(outdir, "matched_quartile_trend.csv"), row.names = FALSE)
print(quartiles, row.names = FALSE, digits = 5)
print(trend, row.names = FALSE, digits = 5)
