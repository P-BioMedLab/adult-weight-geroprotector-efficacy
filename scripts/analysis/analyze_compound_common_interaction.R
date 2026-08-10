# Common treatment-by-weight interaction in the frozen 24-arm, 13-compound
# pre-treatment weight-lowering population.
# Paths are package-relative when installed under scripts/analysis/.
suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
input <- file.path(root, "data", "outputs", "predosing_weight_lowering_27arm_comparisons.csv")
output <- file.path(root, "data", "outputs", "compound_common_interaction.csv")

d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
d$Tx <- as.integer(d$Tx); d$dead <- as.integer(d$dead)
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
d$arm <- factor(d$arm)
d$Txz <- d$Tx * d$weight_z_matched

base_formula <- Surv(entry, lifespan_days, dead) ~ 0 +
  Tx:arm + weight_z_matched:arm + Txz + strata(stratum) + cluster(mouse)
fit <- coxph(base_formula, data = d, ties = "efron", singular.ok = TRUE)
s <- summary(fit)$coefficients["Txz", ]
ci <- summary(fit)$conf.int["Txz", ]

out <- data.frame(
  analysis = "pooled_common_interaction",
  estimate = unname(ci["exp(coef)"]),
  lower_95 = unname(ci["lower .95"]),
  upper_95 = unname(ci["upper .95"]),
  p = unname(s["Pr(>|z|)"]),
  compounds = length(unique(d$compound)), arms = length(unique(d$arm)),
  unique_mice = length(unique(d$mouse))
)
print(out, row.names = FALSE, digits = 7)
stopifnot(nrow(out) == 1, out$compounds[1] == 13, out$arms[1] == 24,
          out$unique_mice[1] == 5523,
          abs(out$estimate[1] - 0.9059) < 0.001)
write.csv(out, output, row.names = FALSE)
