suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
input <- file.path(root, "data", "inputs", "control_weight_gradient_input.csv")
output <- file.path(root, "data", "outputs", "control_weight_gradient_by_age.csv")

d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
d$dead <- as.integer(d$dead)
d$lifespan_days <- as.numeric(d$lifespan_days)
d <- d[d$group == "Control" & is.finite(d$dead) & is.finite(d$lifespan_days), ]

extract <- function(fit, term) {
  estimate <- summary(fit)$coefficients[term, ]
  c(
    HR = unname(exp(estimate["coef"])),
    lo = unname(exp(estimate["coef"] - 1.96 * estimate["se(coef)"])),
    hi = unname(exp(estimate["coef"] + 1.96 * estimate["se(coef)"])),
    p = unname(estimate["Pr(>|z|)"])
  )
}

rows <- list()
for (month in c(6, 12, 18, 24)) {
  column <- paste0("bw_", month)
  entry <- month * 30.4
  d[[column]] <- as.numeric(d[[column]])
  d[[column]][d[[column]] < 10 | d[[column]] > 80] <- NA
  for (sex in c("m", "f")) {
    frame <- d[d$sex == sex & is.finite(d[[column]]) & d$lifespan_days > entry, ]
    frame$weight_z <- ave(frame[[column]], frame$cohort, frame$site, FUN = function(values) {
      (values - mean(values)) / sd(values)
    })
    frame <- frame[is.finite(frame$weight_z), ]
    # Baseline hazards vary across both cohort and site. Retaining every
    # eligible control group also prevents the result from depending on an
    # undocumented site exclusion.
    frame$cohort_site <- interaction(frame$cohort, frame$site, drop = TRUE)
    fit <- coxph(
      Surv(rep(entry, nrow(frame)), lifespan_days, dead) ~
        weight_z + strata(cohort_site),
      data = frame, ties = "efron"
    )
    estimate <- extract(fit, "weight_z")
    rows[[length(rows) + 1]] <- data.frame(
      weight_age = month, sex = sex, n = nrow(frame),
      HR_per_SD = estimate["HR"], lo = estimate["lo"],
      hi = estimate["hi"], p = estimate["p"]
    )
  }
}
result <- do.call(rbind, rows)
write.csv(result, output, row.names = FALSE)
print(result, row.names = FALSE, digits = 5)
