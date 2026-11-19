library(glmnet)

x_train <- as.matrix(select(train_r, starts_with("r_")))
y_train <- train_r$diabetes

x_test  <- as.matrix(select(test_r,  starts_with("r_")))
y_test  <- test_r$diabetes

set.seed(20261118)
cv_fit <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10)

cat("lambda.min:", cv_fit$lambda.min, "\n")
cat("lambda.1se:", cv_fit$lambda.1se, "\n")
plot(cv_fit)

coef_mat <- coef(cv_fit, s = "lambda.min")
weights  <- as.numeric(coef_mat[-1, 1])   # drop intercept
names(weights) <- rownames(coef_mat)[-1]

cat("Non-zero weights:\n"); print(round(weights[weights != 0], 3))
cat("Num non-zero:", sum(weights != 0), "of", length(weights), "\n")

build_composite <- function(df, w) {
  r_cols <- paste0("r_", pct_vars)
  # weighted sum of percentile ranks; negative LASSO weights are allowed
  # (they encode "protective" variables, same as a regression)
  as.numeric(as.matrix(df[, r_cols]) %*% w)
}

train_r$custom_index <- build_composite(train_r, weights)
test_r$custom_index  <- build_composite(test_r,  weights)

# SVI baseline: RPL_THEMES is already a 0-1 percentile rank, comparable scale
summary(train_r$custom_index)
summary(train_r$RPL_THEMES)
cor(train_r$custom_index, train_r$RPL_THEMES)   # how similar are they?

cv_fit_pos <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10, lower.limits = 0)

coef_pos  <- coef(cv_fit_pos, s = "lambda.min")
w_pos     <- as.numeric(coef_pos[-1, 1])
names(w_pos) <- rownames(coef_pos)[-1]
cat("Constrained non-zero:\n"); print(round(w_pos[w_pos != 0], 3))
cat("Constrained zeros:", sum(w_pos == 0), "of", length(w_pos), "\n")

# build + rank both composites
train_r$custom_index_pos <- as.numeric(as.matrix(train_r[, paste0("r_", pct_vars)]) %*% w_pos)
test_r$custom_index_pos  <- as.numeric(as.matrix(test_r[,  paste0("r_", pct_vars)]) %*% w_pos)

# rank-transform (matches SVI construction: 0-1 percentile)
rank_on_train <- function(tr, te, col) {
  e <- ecdf(tr[[col]])
  list(train = e(tr[[col]]), test = e(te[[col]]))
}
rk_uncon <- rank_on_train(train_r, test_r, "custom_index")
rk_con   <- rank_on_train(train_r, test_r, "custom_index_pos")

train_r$custom_rank_uncon <- rk_uncon$train
test_r$custom_rank_uncon  <- rk_uncon$test
train_r$custom_rank_con   <- rk_con$train
test_r$custom_rank_con    <- rk_con$test

cat("Cor (constrained rank vs SVI):", round(cor(train_r$custom_rank_con, train_r$RPL_THEMES), 3), "\n")
cat("Cor (unconstrained rank vs SVI):", round(cor(train_r$custom_rank_uncon, train_r$RPL_THEMES), 3), "\n")

eval_one <- function(df, pred_col) {
  resid <- df$diabetes - df[[pred_col]]
  data.frame(
    predictor = pred_col,
    rmse = sqrt(mean(resid^2)),
    mae  = mean(abs(resid)),
    r2   = 1 - sum(resid^2) / sum((df$diabetes - mean(df$diabetes))^2),
    cor  = cor(df[[pred_col]], df$diabetes)
  )
}

rescale_train <- function(df_train, df_test, pred_col) {
  tr <- data.frame(y = df_train$diabetes, pred = df_train[[pred_col]])
  te <- data.frame(pred = df_test[[pred_col]])
  fit <- lm(y ~ pred, data = tr)
  list(fit = fit, test_pred = predict(fit, newdata = te))
}

# fit each predictor on train, apply to test
svi_rs <- rescale_train(train_r, test_r, "RPL_THEMES")
test_r$svi_fit <- svi_rs$test_pred

con_rs <- rescale_train(train_r, test_r, "custom_rank_con")
test_r$con_fit <- con_rs$test_pred

unc_rs <- rescale_train(train_r, test_r, "custom_rank_uncon")
test_r$unc_fit <- unc_rs$test_pred

# sanity check
cat("Lengths:", length(test_r$svi_fit), length(test_r$con_fit), length(test_r$unc_fit), "\n")
cat("SVI pred range:", round(range(test_r$svi_fit), 2), "\n")
cat("Con pred range:", round(range(test_r$con_fit), 2), "\n")
cat("Unc pred range:", round(range(test_r$unc_fit), 2), "\n")

# headline table
results <- rbind(
  eval_one(test_r, "svi_fit"),
  eval_one(test_r, "con_fit"),
  eval_one(test_r, "unc_fit")
)
print(results)

# deltas vs SVI (positive = custom better)
baseline <- results[results$predictor == "svi_fit", ]
deltas <- results %>%
  filter(predictor != "svi_fit") %>%
  mutate(
    d_r2   = r2   - baseline$r2,
    d_rmse = rmse - baseline$rmse,   # negative = custom better
    d_mae  = mae  - baseline$mae,
    d_cor  = cor  - baseline$cor
  )
cat("\nDeltas vs SVI (d_r2 > 0, d_rmse < 0 = custom better):\n")
print(deltas)


# calibration plot

library(ggplot2)

plot_df <- bind_rows(
  data.frame(observed  = test_r$diabetes,
             predicted = test_r$svi_fit,
             model     = "CDC SVI"),
  data.frame(observed  = test_r$diabetes,
             predicted = test_r$con_fit,
             model     = "Custom (constrained)"),
  data.frame(observed  = test_r$diabetes,
             predicted = test_r$unc_fit,
             model     = "Custom (unconstrained)")
)

p_calib <- ggplot(plot_df, aes(x = predicted, y = observed, color = model)) +
  geom_point(alpha = 0.18, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 40)) +
  labs(
    x = "Predicted diabetes prevalence (%)",
    y = "Observed diabetes prevalence (%)",
    title = "Test-set calibration: SVI vs custom SDOH indices",
    subtitle = "Dashed line = perfect calibration",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_calib)
ggsave("calibration_test.png", p_calib, width = 7, height = 5.5, dpi = 150)

binned <- plot_df %>%
  group_by(model) %>%
  mutate(bin = ntile(predicted, 10)) %>%
  group_by(model, bin) %>%
  summarise(
    mean_pred = mean(predicted),
    mean_obs  = mean(observed),
    n         = n(),
    .groups   = "drop"
  )

print(binned, n = 30)

p_bin <- ggplot(binned, aes(x = mean_pred, y = mean_obs, color = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "Mean predicted (%)", y = "Mean observed (%)",
    title = "Decile-binned calibration",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_bin)
ggsave("calibration_binned.png", p_bin, width = 7, height = 5, dpi = 150)

# Give SVI the same functional flexibility as the custom index
library(splines)

svi_spline <- lm(diabetes ~ ns(RPL_THEMES, df = 4), data = train_r)
test_r$svi_spline_fit <- predict(svi_spline, newdata = test_r)

# also try decile-factor version (SVI's native use)
svi_decile <- lm(diabetes ~ factor(ntile(RPL_THEMES, 10)), data = train_r)
test_r$svi_decile_fit <- predict(svi_decile, newdata = test_r)

rbind(
  eval_one(test_r, "svi_fit"),          # linear
  eval_one(test_r, "svi_spline_fit"),   # spline
  eval_one(test_r, "svi_decile_fit"),   # decile factor
  eval_one(test_r, "con_fit"),          # custom constrained
  eval_one(test_r, "unc_fit")           # custom unconstrained
)

# refit the whole pipeline on F_TOTAL == 0 subset only
sub <- svi_pct %>% filter(F_TOTAL == 0)

# quick version: just re-evaluate the existing fits on the F_TOTAL==0 test tracts
test_sub <- test_r %>% filter(F_TOTAL == 0)
cat("Test tracts with F_TOTAL==0:", nrow(test_sub), "\n")

rbind(
  eval_one(test_sub, "svi_fit"),
  eval_one(test_sub, "svi_spline_fit"),
  eval_one(test_sub, "con_fit"),
  eval_one(test_sub, "unc_fit")
)