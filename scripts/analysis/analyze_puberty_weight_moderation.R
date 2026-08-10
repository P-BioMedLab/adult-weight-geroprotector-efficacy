suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
input <- file.path(root, "data", "inputs", "itp_gn_puberty_weight_data.csv")
output <- file.path(root, "data", "outputs", "itp_gn_puberty_weight_moderation.csv")

d <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
d$Tx <- as.numeric(d$Tx)
d$dead <- as.integer(d$dead)
d$positive_m <- tolower(as.character(d$positive_m)) %in% c("true", "1")
d$positive_f <- tolower(as.character(d$positive_f)) %in% c("true", "1")
d$interaction <- d$Tx * d$weight_42d_z

fit_sex <- function(x, sex) {
  x$armsex <- factor(paste(x$arm, x$sex, sep = "|"))
  x$stratum <- factor(paste(x$arm, x$site, x$sex, sep = "|"))
  x$mouse <- factor(x$mouse)
  fit <- coxph(
    Surv(entry, age, dead) ~ 0 + armsex:Tx + armsex:weight_42d_z +
      interaction + strata(stratum) + cluster(mouse),
    data = x, ties = "efron", singular.ok = TRUE
  )
  estimate <- summary(fit)$coefficients["interaction", ]
  data.frame(
    label = paste("pubertal weight", sex), sex = sex,
    arms = length(unique(x$arm)), unique_mice = length(unique(x$mouse)),
    HR = exp(estimate["coef"]),
    lo = exp(estimate["coef"] - 1.96 * estimate["robust se"]),
    hi = exp(estimate["coef"] + 1.96 * estimate["robust se"]),
    p = estimate["Pr(>|z|)"]
  )
}

rows <- list()
for (sex in c("m", "f")) {
  eligible <- if (sex == "m") d$positive_m else d$positive_f
  frame <- d[d$sex == sex & eligible, ]
  valid_arms <- names(which(tapply(frame$Tx, frame$arm, function(values) {
    length(unique(values)) == 2 && min(table(values)) >= 15
  })))
  rows[[length(rows) + 1]] <- fit_sex(frame[frame$arm %in% valid_arms, ], sex)
}
result <- do.call(rbind, rows)
write.csv(result, output, row.names = FALSE)
print(result, row.names = FALSE, digits = 4)
