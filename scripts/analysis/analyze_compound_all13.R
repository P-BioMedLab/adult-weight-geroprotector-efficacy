# Descriptive quartile analysis for every single compound in the frozen
# pre-treatment weight-lowering set. Paths are package-relative.
suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
input <- file.path(root, "data", "outputs", "predosing_weight_lowering_27arm_comparisons.csv")
output <- file.path(root, "data", "outputs", "compound_all13_descriptive.csv")

d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
d$Tx <- as.integer(d$Tx); d$dead <- as.integer(d$dead)
d$matched_quartile <- as.integer(d$matched_quartile)
raw <- sub("^[^|]+\\|", "", d$arm)
d$compound <- ifelse(grepl("^17aE2", raw), "17alpha-estradiol",
  ifelse(grepl("^ACA_", raw), "Acarbose",
  ifelse(grepl("^Cana", raw), "Canagliflozin",
  ifelse(grepl("^[Rr]apa(_|$)", raw), "Rapamycin",
  ifelse(raw == "CC", "Candesartan",
  ifelse(raw == "DMF_9", "Dimethyl fumarate",
  ifelse(raw == "Fis_Cyc", "Fisetin",
  ifelse(raw == "GGA", "Geranylgeranylacetone",
  ifelse(raw == "Gly", "Glycine",
  ifelse(raw == "MIF098", "MIF098",
  ifelse(raw == "Mec", "Meclizine",
  ifelse(raw == "NDGA", "NDGA",
  ifelse(raw == "OH-EST", "16alpha-hydroxyestradiol", NA)))))))))))))
d <- d[!is.na(d$compound), ]

weighted_median <- function(z) {
  fit <- survfit(Surv(entry, lifespan_days, dead) ~ 1, data = z,
                 weights = standardization_weight, conf.type = "none")
  as.numeric(summary(fit)$table["median"])
}
estimate <- function(z) {
  cmed <- weighted_median(z[z$Tx == 0, ])
  tmed <- weighted_median(z[z$Tx == 1, ])
  c(control = cmed, treated = tmed, gain_pct = 100 * (tmed - cmed) / cmed)
}

rows <- lapply(sort(unique(d$compound)), function(nm) {
  z <- d[d$compound == nm, ]; lo <- estimate(z[z$matched_quartile == 1, ])
  ov <- estimate(z); hi <- estimate(z[z$matched_quartile == 4, ])
  data.frame(compound = nm, arms = length(unique(z$arm)), unique_mice = length(unique(z$mouse)),
    lowest_control_median = lo["control"], lowest_treated_median = lo["treated"], lowest_gain_pct = lo["gain_pct"],
    overall_control_median = ov["control"], overall_treated_median = ov["treated"], overall_gain_pct = ov["gain_pct"],
    highest_control_median = hi["control"], highest_treated_median = hi["treated"], highest_gain_pct = hi["gain_pct"])
})
out <- do.call(rbind, rows)
out <- out[order(-out$overall_gain_pct), ]
lo <- estimate(d[d$matched_quartile == 1, ]); ov <- estimate(d); hi <- estimate(d[d$matched_quartile == 4, ])
out <- rbind(out, data.frame(compound = "Pooled single-agent set", arms = length(unique(d$arm)),
  unique_mice = length(unique(d$mouse)), lowest_control_median = lo["control"], lowest_treated_median = lo["treated"],
  lowest_gain_pct = lo["gain_pct"], overall_control_median = ov["control"], overall_treated_median = ov["treated"],
  overall_gain_pct = ov["gain_pct"], highest_control_median = hi["control"], highest_treated_median = hi["treated"],
  highest_gain_pct = hi["gain_pct"]))
stopifnot(nrow(out) == 14, tail(out$unique_mice, 1) == 5523,
          sum(out$highest_gain_pct[-14] > out$lowest_gain_pct[-14]) == 12)
write.csv(out, output, row.names = FALSE)
print(out, row.names = FALSE, digits = 4)
