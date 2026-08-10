suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript prepare_developmental_figure.R input.csv output.csv")
x <- read.csv(args[1], stringsAsFactors = FALSE)

fit_row <- function(data, formula, term, panel, label, sex) {
  f <- coxph(formula, data = data, ties = "efron")
  s <- summary(f)$coefficients[term, ]
  data.frame(
    panel = panel, label = label, sex = sex, n = nrow(data),
    HR = exp(s["coef"]),
    lo = exp(s["coef"] - 1.96 * s["se(coef)"]),
    hi = exp(s["coef"] + 1.96 * s["se(coef)"]),
    p = s["Pr(>|z|)"]
  )
}

rows <- list()
for (sx in c("m", "f")) {
  z <- subset(x, sex == sx & is.finite(weight_42d) & age > 183)
  z$w14z <- ave(z$weight_42d, z$cohort, z$site,
                FUN = function(v) (v - mean(v)) / sd(v))
  z <- subset(z, is.finite(w14z))
  rows[[length(rows) + 1]] <- fit_row(
    z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + strata(cohort, site),
    "w14z", "age", "1.4", sx
  )
}

z <- subset(x, sex == "m" & is.finite(weight_42d) & is.finite(weight_6m))
z$w14z <- ave(z$weight_42d, z$cohort, z$site,
              FUN = function(v) (v - mean(v)) / sd(v))
z$w6z <- ave(z$weight_6m, z$cohort, z$site,
             FUN = function(v) (v - mean(v)) / sd(v))
z <- subset(z, age > 183 & is.finite(w14z) & is.finite(w6z))

rows[[length(rows) + 1]] <- fit_row(
  z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + strata(cohort, site),
  "w14z", "decomposition", "1.4-month alone", "m"
)
rows[[length(rows) + 1]] <- fit_row(
  z, Surv(rep(183, nrow(z)), age, dead) ~ w6z + strata(cohort, site),
  "w6z", "decomposition", "6-month alone", "m"
)
rows[[length(rows) + 1]] <- fit_row(
  z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + w6z + strata(cohort, site),
  "w14z", "decomposition", "1.4-month adjusted", "m"
)
rows[[length(rows) + 1]] <- fit_row(
  z, Surv(rep(183, nrow(z)), age, dead) ~ w14z + w6z + strata(cohort, site),
  "w6z", "decomposition", "6-month adjusted", "m"
)

out <- do.call(rbind, rows)
write.csv(out, args[2], row.names = FALSE)
print(out, row.names = FALSE, digits = 5)
