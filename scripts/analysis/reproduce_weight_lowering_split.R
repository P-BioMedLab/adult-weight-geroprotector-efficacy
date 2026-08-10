suppressPackageStartupMessages(library(survival))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
package <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
d <- read.csv(file.path(package, "data", "outputs", "itp_allarm_withinarm_long_reproduced.csv"), stringsAsFactors = FALSE)
w <- read.csv(file.path(package, "data", "outputs", "weight_effect_by_arm.csv"), stringsAsFactors = FALSE)
meta <- read.csv(file.path(package, "data", "outputs", "withinarm_interactions.csv"), stringsAsFactors = FALSE)

w$post_dosing_flag <- tolower(as.character(w$post_dosing)) == "true"
w <- w[w$post_dosing_flag, ]
w$weighted_diff <- w$diff_g * w$n_tx
agg <- aggregate(cbind(weighted_diff, n_tx) ~ cohort + group, data = w, sum)
agg$mean_diff_g <- agg$weighted_diff / agg$n_tx
agg$lowering <- agg$mean_diff_g < -1

meta$armsex <- paste(meta$cohort, meta$group, meta$sex, sep = "|")
meta <- merge(meta, unique(w[c("cohort", "group", "init")]), by = c("cohort", "group"), all.x = TRUE)
meta$predosing <- is.finite(meta$init) & meta$weight_age < meta$init

d$cohort <- sub("\\|.*$", "", d$armsex)
d$remainder <- sub("^[^|]+\\|", "", d$armsex)
d$group <- sub("\\|[^|]+$", "", d$remainder)
d <- merge(d, agg[c("cohort", "group", "lowering")], by = c("cohort", "group"), all.x = TRUE)
d$lowering[is.na(d$lowering)] <- FALSE
d <- merge(d, meta[c("armsex", "predosing")], by = "armsex", all.x = TRUE)
d$dead <- as.integer(d$dead)
d$Tx <- as.numeric(d$Tx)
d$wz <- as.numeric(d$wz)
d$interaction <- d$Tx * d$wz

fit_subset <- function(x, label) {
  x$armsex <- factor(x$armsex)
  x$stratum <- factor(x$stratum)
  x$mouse <- factor(x$mouse)
  f <- coxph(
    Surv(entry, lifespan_days, dead) ~ 0 + armsex:Tx + armsex:wz +
      interaction + strata(stratum) + cluster(mouse),
    data = x, ties = "efron", singular.ok = TRUE
  )
  s <- summary(f)$coefficients["interaction", ]
  data.frame(
    analysis = label,
    arms = length(unique(x$armsex)),
    rows = nrow(x),
    unique_mice = length(unique(x$mouse)),
    HR = exp(s["coef"]),
    lo = exp(s["coef"] - 1.96 * s["robust se"]),
    hi = exp(s["coef"] + 1.96 * s["robust se"]),
    p = s["Pr(>|z|)"]
  )
}

fit_contrast <- function(x, label) {
  x$armsex <- factor(x$armsex)
  x$stratum <- factor(x$stratum)
  x$mouse <- factor(x$mouse)
  x$lowering_num <- as.numeric(x$lowering)
  f <- coxph(
    Surv(entry, lifespan_days, dead) ~ 0 + armsex:Tx + armsex:wz +
      interaction + interaction:lowering_num + strata(stratum) + cluster(mouse),
    data = x, ties = "efron", singular.ok = TRUE
  )
  s <- summary(f)$coefficients["interaction:lowering_num", ]
  data.frame(
    analysis = label,
    arms = length(unique(x$armsex)),
    rows = nrow(x),
    unique_mice = length(unique(x$mouse)),
    HR = exp(s["coef"]),
    lo = exp(s["coef"] - 1.96 * s["robust se"]),
    hi = exp(s["coef"] + 1.96 * s["robust se"]),
    p = s["Pr(>|z|)"]
  )
}

male <- d[d$sex == "m", ]
female <- d[d$sex == "f", ]
male_pre <- male[male$predosing %in% TRUE, ]

out <- do.call(rbind, list(
  fit_subset(male[male$lowering, ], "male_weight_lowering"),
  fit_subset(male[!male$lowering, ], "male_nonlowering"),
  fit_contrast(male, "male_class_contrast"),
  fit_subset(female[female$lowering, ], "female_weight_lowering"),
  fit_subset(female[!female$lowering, ], "female_nonlowering"),
  fit_subset(male_pre, "male_predosing_all"),
  fit_subset(male_pre[male_pre$lowering, ], "male_predosing_weight_lowering"),
  fit_subset(male_pre[!male_pre$lowering, ], "male_predosing_nonlowering"),
  fit_contrast(male_pre, "male_predosing_class_contrast")
))

write.csv(out, file.path(package, "data", "outputs", "weight_lowering_split_reproduced.csv"), row.names = FALSE)
print(out, row.names = FALSE, digits = 6)
