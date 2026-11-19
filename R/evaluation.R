set.seed(20261118)
kf <- 10
folds <- sample(rep(1:kf, length.out = nrow(svi_pct)))

delta_r2 <- numeric(kf)
for (k in 1:kf) {
  tr <- svi_pct[folds != k, ]
  te <- svi_pct[folds == k, ]
  
  # rebuild ranks on tr, map to te
  ecdfs <- lapply(pct_vars, function(v) ecdf(tr[[v]]))
  names(ecdfs) <- pct_vars
  rankify <- function(df) {
    as.data.frame(lapply(pct_vars, function(v) ecdfs[[v]](df[[v]]))) |>
      setNames(paste0("r_", pct_vars))
  }
  tr_r <- cbind(tr, rankify(tr))
  te_r <- cbind(te, rankify(te))
  
  # fit LASSO on tr (constrained), build composite, rank, calibrate
  x_tr <- as.matrix(tr_r[, paste0("r_", pct_vars)])
  fit  <- cv.glmnet(x_tr, tr_r$diabetes, alpha = 1, lower.limits = 0, nfolds = 5)
  w    <- as.numeric(coef(fit, s = "lambda.min")[-1, 1])
  tr_r$ci <- as.numeric(x_tr %*% w)
  te_r$ci <- as.numeric(as.matrix(te_r[, paste0("r_", pct_vars)]) %*% w)
  
  e <- ecdf(tr_r$ci)
  tr_r$ci_rank <- e(tr_r$ci); te_r$ci_rank <- e(te_r$ci)
  
  # calibrate + eval
  m_c <- lm(diabetes ~ ci_rank, data = tr_r)
  m_s <- lm(diabetes ~ RPL_THEMES, data = tr_r)
  te_r$pc <- predict(m_c, te_r); te_r$ps <- predict(m_s, te_r)
  
  r2c <- 1 - sum((te_r$diabetes - te_r$pc)^2) / sum((te_r$diabetes - mean(te_r$diabetes))^2)
  r2s <- 1 - sum((te_r$diabetes - te_r$ps)^2) / sum((te_r$diabetes - mean(te_r$diabetes))^2)
  delta_r2[k] <- r2c - r2s
}
cat("ΔR² mean:", round(mean(delta_r2), 3),
    "| 95% CI:", round(quantile(delta_r2, c(0.025, 0.975)), 3), "\n")

county_results <- test_r %>%
  group_by(county) %>%
  summarise(
    n = n(),
    r2_svi = 1 - sum((diabetes - svi_fit)^2) / sum((diabetes - mean(diabetes))^2),
    r2_con = 1 - sum((diabetes - con_fit)^2) / sum((diabetes - mean(diabetes))^2),
    delta  = r2_con - r2_svi,
    .groups = "drop"
  ) %>%
  arrange(delta)

print(county_results, n = 60)
summary(county_results$delta)
cat("Counties where custom wins:", sum(county_results$delta > 0), "of", nrow(county_results), "\n")

county_results_clean <- county_results %>%
  filter(is.finite(delta), n >= 10) %>%
  arrange(delta)

cat("Counties with n>=10:", nrow(county_results_clean), "\n")
cat("Custom wins:", sum(county_results_clean$delta > 0),
    "of", nrow(county_results_clean), "\n")
cat("Median delta:", round(median(county_results_clean$delta), 3), "\n")
cat("Min delta:", round(min(county_results_clean$delta), 3),
    "in county", county_results_clean$county[1], "\n")
cat("Max delta:", round(max(county_results_clean$delta), 3),
    "in county", county_results_clean$county[nrow(county_results_clean)], "\n")

county_results_clean %>%
  mutate(wins = delta > 0) %>%
  summarise(
    n_counties = n(),
    n_wins = sum(wins),
    pct_wins = round(100 * mean(wins), 1),
    median_delta = round(median(delta), 3),
    weighted_mean_delta = round(weighted.mean(delta, w = n), 3)
  )
