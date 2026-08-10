suppressPackageStartupMessages(library(survival))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
input <- file.path(package, "data", "inputs", "early_weight_change_input.csv")
output <- file.path(package, "data", "outputs", "early_weight_change_reproduced.csv")
d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
d$is_control <- tolower(as.character(d$is_control)) == "true"

pooled_within_sd <- function(x, cell) {
  centered <- x - ave(x, cell, FUN = mean)
  sqrt(sum(centered^2) / (length(x) - length(unique(cell))))
}

prepare <- function(sex_value, timing) {
  w <- d[d$sex == sex_value, ]
  if (timing == "primary") {
    w <- w[w$is_control | (!is.na(w$init_months) & w$init_months < 12), ]
  } else if (timing == "treated_primary") {
    w <- w[!w$is_control & !is.na(w$init_months) & w$init_months < 12, ]
  } else if (timing == "after_6") {
    w <- w[w$is_control | (!is.na(w$init_months) & w$init_months > 6), ]
  } else if (timing == "at_or_after_12") {
    w <- w[w$is_control | (!is.na(w$init_months) & w$init_months >= 12), ]
  } else {
    stop(timing)
  }
  event_count <- tapply(w$dead, w$cell, sum)
  valid <- names(event_count)[event_count >= 12]
  w <- w[w$cell %in% valid, ]
  w$change_z <- (w$pct_6_12 - ave(w$pct_6_12, w$cell, FUN = mean)) /
    pooled_within_sd(w$pct_6_12, w$cell)
  w$w6_z <- (w$bw_6 - ave(w$bw_6, w$cell, FUN = mean)) /
    pooled_within_sd(w$bw_6, w$cell)
  w
}

fit_one <- function(sex_value, timing, adjusted) {
  w <- prepare(sex_value, timing)
  form <- if (adjusted) {
    Surv(lifespan_days - 365, dead) ~ change_z + w6_z + strata(cell) + cluster(cell)
  } else {
    Surv(lifespan_days - 365, dead) ~ change_z + strata(cell) + cluster(cell)
  }
  fit <- coxph(form, data = w, ties = "efron", model = FALSE, x = FALSE, y = FALSE)
  tab <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  row <- data.frame(
    sex = sex_value,
    timing = timing,
    adjusted_for_w6 = adjusted,
    n_mice = nrow(w),
    events = sum(w$dead),
    cells = length(unique(w$cell)),
    cohorts = length(unique(w$cohort)),
    hr = unname(ci["change_z", "exp(coef)"]),
    ci_low = unname(ci["change_z", "lower .95"]),
    ci_high = unname(ci["change_z", "upper .95"]),
    p = unname(tab["change_z", "Pr(>|z|)"]),
    w6_hr = NA_real_,
    w6_ci_low = NA_real_,
    w6_ci_high = NA_real_,
    w6_p = NA_real_
  )
  if (adjusted) {
    row$w6_hr <- unname(ci["w6_z", "exp(coef)"])
    row$w6_ci_low <- unname(ci["w6_z", "lower .95"])
    row$w6_ci_high <- unname(ci["w6_z", "upper .95"])
    row$w6_p <- unname(tab["w6_z", "Pr(>|z|)"])
  }
  print(row)
  row
}

specs <- list(
  c("m", "primary", FALSE),
  c("m", "primary", TRUE),
  c("f", "primary", FALSE),
  c("f", "primary", TRUE),
  c("m", "treated_primary", TRUE),
  c("m", "after_6", TRUE),
  c("m", "at_or_after_12", TRUE)
)

rows <- lapply(specs, function(s) fit_one(s[[1]], s[[2]], as.logical(s[[3]])))
result <- do.call(rbind, rows)
write.csv(result, output, row.names = FALSE)
cat("Saved", output, "\n")
