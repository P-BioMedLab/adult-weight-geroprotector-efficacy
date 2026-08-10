script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
input <- file.path(package, "data", "outputs", "comparator_effect_cells.csv")
output <- file.path(package, "data", "outputs", "comparator_split_reproduced.csv")
d <- read.csv(input, stringsAsFactors = FALSE)

# Tertiles and other displayed splits are empirical quantile bins of the
# observed control-group mean weights, defined separately by sex. This is the
# specification used for Table 4 (not rank/n bins, which assign weight ties to
# different cut points and answer a slightly different question).
rows <- list()
for (sx in c("m", "f")) {
  x <- d[d$sex == sx, ]
  for (k in c(3, 4, 5, 10)) {
    breaks <- unique(quantile(x$ctl_bw6, probs = 0:k / k, type = 7, na.rm = TRUE))
    if (length(breaks) != k + 1) stop("Non-unique empirical quantile cut points")
    x$bin <- cut(x$ctl_bw6, breaks = breaks, include.lowest = TRUE,
                 labels = paste0("Q", seq_len(k)))
    for (b in levels(x$bin)) {
      z <- x[x$bin == b, ]; w <- 1 / z$se^2
      beta <- sum(z$log_hr * w) / sum(w); se <- sqrt(1 / sum(w))
      rows[[length(rows) + 1]] <- data.frame(
        sex = sx, split = k, bin = b, cells = nrow(z),
        mean_control_weight = mean(z$ctl_bw6), HR = exp(beta),
        reduction_pct = 100 * (1 - exp(beta)),
        lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se)
      )
    }
  }
}
out <- do.call(rbind, rows)
write.csv(out, output, row.names = FALSE)
print(out[out$split == 3 & out$bin %in% c("Q1", "Q3"), ], row.names = FALSE, digits = 6)
