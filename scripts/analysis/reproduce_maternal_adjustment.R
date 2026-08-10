# Standalone reproduction of the maternal-weight adjustment analysis.
# Paths are resolved relative to this script, so it is independent of cwd.
suppressPackageStartupMessages(library(survival))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
root <- normalizePath(file.path(dirname(script), "..", ".."), winslash = "/")
input <- file.path(root, "data", "inputs", "itp_gn_earlylife_controls.csv")
output <- file.path(root, "data", "outputs", "itp_gn_maternal_adjustment.csv")

d <- read.csv(input, check.names = FALSE)
names(d) <- tolower(names(d))

# Frozen input columns are used explicitly; fail clearly if the package changes.
required <- c("sex", "age", "dead", "site", "cohort", "weight_6m", "maternal_age_days", "parity")
missing <- setdiff(required, names(d))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

extract <- function(fit, term) {
  b <- unname(coef(fit)[term])
  se <- unname(sqrt(vcov(fit)[term, term]))
  c(HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))
}

rows <- list()
for (sx in c("m", "f")) {
  z <- d[d$sex == sx & complete.cases(d[, required]), ]
  z <- z[z$age > 183, ]
  z$w6z <- as.numeric(scale(z$weight_6m))
  z$maz <- as.numeric(scale(z$maternal_age_days))
  z$pz <- as.numeric(scale(z$parity))
  f0 <- coxph(Surv(rep(183, nrow(z)), age, dead) ~ w6z + strata(cohort, site), data = z, ties = "efron")
  f1 <- coxph(Surv(rep(183, nrow(z)), age, dead) ~ w6z + maz + pz + strata(cohort, site), data = z, ties = "efron")
  a <- extract(f0, "w6z")
  b <- extract(f1, "w6z")
  rows[[length(rows) + 1]] <- data.frame(
    sex = sx, n = nrow(z),
    HR_weight = a["HR"], lo = a["lo"], hi = a["hi"], p = a["p"],
    HR_adjusted_maternal = b["HR"], adj_lo = b["lo"], adj_hi = b["hi"], adj_p = b["p"]
  )
}
out <- do.call(rbind, rows)
write.csv(out, output, row.names = FALSE)
print(out)
