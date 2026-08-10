suppressPackageStartupMessages(library(survival))
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
x <- read.csv(file.path(package, "data", "inputs", "itp_gn_earlylife_controls.csv"))
x <- subset(x, sex == "m" & is.finite(weight_42d) & is.finite(weight_6m))
x$w42z <- ave(x$weight_42d, x$cohort, x$site,
              FUN=function(z)(z-mean(z))/sd(z))
x$w6z <- ave(x$weight_6m, x$cohort, x$site,
             FUN=function(z)(z-mean(z))/sd(z))
x <- subset(x, age > 183 & is.finite(w42z) & is.finite(w6z))
f <- coxph(Surv(rep(183,nrow(x)),age,dead) ~ w42z + w6z +
             strata(cohort,site), data=x, ties="efron")
print(summary(f)$coefficients)
print(exp(cbind(coef(f),confint(f))))
s <- summary(f)$coefficients
ci <- exp(confint(f))
out <- data.frame(
  term = rownames(s), n = nrow(x), HR = exp(s[, "coef"]),
  lo = ci[, 1], hi = ci[, 2], p = s[, "Pr(>|z|)"]
)
write.csv(out, file.path(package, "data", "outputs", "puberty_joint_reproduced.csv"), row.names = FALSE)
