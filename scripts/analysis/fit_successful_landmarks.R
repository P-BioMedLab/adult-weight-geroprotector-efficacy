suppressPackageStartupMessages(library(survival))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
input <- file.path(package, "data", "outputs", "successful_landmarks_long.csv")
outdir <- file.path(package, "data", "outputs")
d <- read.csv(input, stringsAsFactors = FALSE)
d$dead <- as.integer(d$dead); d$Tx <- as.numeric(d$Tx)
d$positive_m <- as.logical(d$positive_m); d$positive_f <- as.logical(d$positive_f)
d$interaction <- d$Tx * d$weight_z

fit_pool <- function(x, label) {
  x$armstage <- factor(x$arm_stage)
  x$stratum <- factor(paste(x$arm_stage, x$site, x$sex, sep = "|"))
  x$mouse <- factor(x$mouse)
  fit <- coxph(
    Surv(entry, lifespan_days, dead) ~ 0 + armstage:Tx + armstage:weight_z + interaction +
      strata(stratum) + cluster(mouse),
    data = x, ties = "efron", singular.ok = TRUE
  )
  s <- summary(fit)$coefficients["interaction", ]
  data.frame(
    analysis = label, landmark = unique(x$landmark), sex = unique(x$sex),
    arms = length(unique(x$arm)), rows = nrow(x), unique_mice = length(unique(x$mouse)),
    HR = exp(s["coef"]), lo = exp(s["coef"] - 1.96 * s["robust se"]),
    hi = exp(s["coef"] + 1.96 * s["robust se"]), p = s["Pr(>|z|)"]
  )
}

fit_grid <- function(data, lag_days = 0) {
  if (lag_days > 0) {
    data <- data[data$lifespan_days > data$entry + lag_days, ]
    data$entry <- data$entry + lag_days
  }
  rows <- list()
  for (sx in c("m", "f")) for (lm in sort(unique(data$landmark))) {
    positive <- if (sx == "m") data$positive_m else data$positive_f
    x <- data[data$sex == sx & data$landmark == lm & positive, ]
    good <- names(which(tapply(x$Tx, x$arm_stage, function(z) length(unique(z)) == 2 && min(table(z)) >= 15)))
    x <- x[x$arm_stage %in% good, ]
    if (length(unique(x$arm)) >= 1) rows[[length(rows) + 1]] <- fit_pool(x, ifelse(lag_days, "90-day lag", "published extenders"))
  }
  do.call(rbind, rows)
}

main <- fit_grid(d, 0)
lag90 <- fit_grid(d, 90)
write.csv(main, file.path(outdir, "successful_landmarks_reproduced.csv"), row.names = FALSE)
write.csv(lag90, file.path(outdir, "successful_landmarks_lag90_reproduced.csv"), row.names = FALSE)

x0 <- d[d$sex == "m" & d$landmark == 12 & d$positive_m, ]
loo <- list()
for (drop in c("none", sort(unique(x0$cohort)))) {
  x <- if (drop == "none") x0 else x0[x0$cohort != drop, ]
  good <- names(which(tapply(x$Tx, x$arm_stage, function(z) length(unique(z)) == 2 && min(table(z)) >= 15)))
  x <- x[x$arm_stage %in% good, ]
  z <- fit_pool(x, paste("drop", drop)); z$drop_cohort <- drop
  loo[[length(loo) + 1]] <- z
}
loo <- do.call(rbind, loo)
write.csv(loo, file.path(outdir, "successful_landmarks_loo12_reproduced.csv"), row.names = FALSE)
print(main, row.names = FALSE, digits = 5)
