suppressPackageStartupMessages(library(survival))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
input <- file.path(package, "data", "outputs", "itp_allarm_withinarm_long_reproduced.csv")
output <- file.path(package, "data", "outputs", "itp_allarm_primary_reproduced.csv")

d <- read.csv(input, stringsAsFactors = FALSE)
d$dead <- as.integer(d$dead); d$Tx <- as.numeric(d$Tx); d$wz <- as.numeric(d$wz)
d$interaction <- d$Tx * d$wz
rows <- list()
for (sx in c("m", "f")) {
  x <- d[d$sex == sx, ]
  x$armsex <- factor(x$armsex); x$stratum <- factor(x$stratum); x$mouse <- factor(x$mouse)
  fit <- coxph(
    Surv(entry, lifespan_days, dead) ~ 0 + armsex:Tx + armsex:wz + interaction +
      strata(stratum) + cluster(mouse),
    data = x, ties = "efron", singular.ok = TRUE
  )
  s <- summary(fit)$coefficients["interaction", ]
  rows[[length(rows) + 1]] <- data.frame(
    sex = sx, arms = length(unique(x$armsex)), rows = nrow(x),
    unique_mice = length(unique(x$mouse)), deaths = sum(x$dead),
    HR = exp(s["coef"]), lo = exp(s["coef"] - 1.96 * s["robust se"]),
    hi = exp(s["coef"] + 1.96 * s["robust se"]), p = s["Pr(>|z|)"]
  )
}
out <- do.call(rbind, rows)
write.csv(out, output, row.names = FALSE)
print(out, row.names = FALSE, digits = 6)
