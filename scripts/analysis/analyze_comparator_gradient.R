script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
input <- file.path(package, "data", "outputs", "comparator_effect_cells.csv")
output <- file.path(package, "data", "outputs", "comparator_gradient_sensitivity.csv")

d <- read.csv(input, stringsAsFactors = FALSE)
stopifnot(nrow(d) == 426, sum(d$sex == "m") == 219, sum(d$sex == "f") == 207)
d$cluster <- paste(d$cohort, d$site, sep = "|")

fit_gradient <- function(z, cohort_fe = FALSE, clustered = FALSE) {
  x <- z$ctl_bw6
  y <- z$log_hr
  w <- 1 / z$se^2
  X <- cbind(intercept = 1, weight = x)
  if (cohort_fe) {
    X <- cbind(X, model.matrix(~ factor(cohort), data = z)[, -1, drop = FALSE])
  }
  bread <- solve(crossprod(X, w * X))
  beta <- bread %*% crossprod(X, w * y)
  residual <- as.vector(y - X %*% beta)
  if (clustered) {
    groups <- unique(z$cluster)
    meat <- matrix(0, ncol(X), ncol(X))
    for (group in groups) {
      selected <- z$cluster == group
      score <- crossprod(X[selected, , drop = FALSE], w[selected] * residual[selected])
      meat <- meat + score %*% t(score)
    }
    covariance <- bread %*% meat %*% bread * length(groups) / (length(groups) - 1)
    degrees <- length(groups) - 1
    critical <- qt(0.975, degrees)
    p <- 2 * pt(-abs(beta[2] / sqrt(covariance[2, 2])), degrees)
  } else {
    covariance <- bread
    critical <- qnorm(0.975)
    p <- 2 * pnorm(-abs(beta[2] / sqrt(covariance[2, 2])))
  }
  slope_se <- sqrt(covariance[2, 2])
  weight_sd <- sd(x)
  data.frame(
    cells = nrow(z), control_clusters = length(unique(z$cluster)),
    hr_per_sd = exp(beta[2] * weight_sd),
    ci_low = exp((beta[2] - critical * slope_se) * weight_sd),
    ci_high = exp((beta[2] + critical * slope_se) * weight_sd),
    p = as.numeric(p)
  )
}

rows <- list()
for (sx in c("m", "f")) {
  z <- d[d$sex == sx, ]
  rows[[length(rows) + 1]] <- cbind(sex = sx, analysis = "within_cohort",
                                    fit_gradient(z, cohort_fe = TRUE))
  rows[[length(rows) + 1]] <- cbind(sex = sx, analysis = "shared_control_clustered",
                                    fit_gradient(z, clustered = TRUE))
}
result <- do.call(rbind, rows)
write.csv(result, output, row.names = FALSE)

expected <- c(m_within_cohort = 0.008, f_within_cohort = 0.63,
              m_shared_control_clustered = 0.078, f_shared_control_clustered = 0.37)
for (i in seq_len(nrow(result))) {
  key <- paste(result$sex[i], result$analysis[i], sep = "_")
  tolerance <- if (expected[key] < 0.1) 0.0015 else 0.015
  stopifnot(abs(result$p[i] - expected[key]) <= tolerance)
}
print(result, row.names = FALSE, digits = 6)
