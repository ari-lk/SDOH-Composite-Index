library(dplyr)
library(tidyr)

svi_pct <- merged %>%
  mutate(
    p_pov150  = E_POV150  / E_TOTPOP * 100,
    p_unemp   = E_UNEMP   / E_TOTPOP * 100,
    p_nohsdp  = E_NOHSDP  / E_TOTPOP * 100,
    p_uninsur = E_UNINSUR / E_TOTPOP * 100,
    p_age65   = E_AGE65   / E_TOTPOP * 100,
    p_age17   = E_AGE17   / E_TOTPOP * 100,
    p_disabl  = E_DISABL  / E_TOTPOP * 100,
    p_sngpnt  = E_SNGPNT  / E_TOTPOP * 100,
    p_limeng  = E_LIMENG  / E_TOTPOP * 100,
    p_minrty  = E_MINRTY  / E_TOTPOP * 100,
    p_groupq  = E_GROUPQ  / E_TOTPOP * 100,
    p_hburd   = E_HBURD   / E_HH     * 100,
    p_crowd   = E_CROWD   / E_HH     * 100,
    p_noveh   = E_NOVEH   / E_HH     * 100,
    p_munit   = E_MUNIT   / E_HU     * 100,
    p_mobile  = E_MOBILE  / E_HU     * 100
  ) %>%
  mutate(across(starts_with("p_"), ~ ifelse(is.finite(.x), .x, NA_real_)))

cat("Rows with any NA predictor:",
    sum(!complete.cases(select(svi_pct, starts_with("p_")))), "\n")

# sanity check: compare your computed p_pov150 to the SVI's EP_POV150? 
# we didn't keep EP_* — skip; the validation in Phase 1 already confirmed denominators
summary(select(svi_pct, starts_with("p_")))

svi_pct <- svi_pct %>%
  filter(complete.cases(select(., starts_with("p_"))))
cat("Analysis tracts:", nrow(svi_pct), "\n")

set.seed(20261118)

svi_pct <- svi_pct %>% mutate(county = substr(GEOID, 1, 5))

train_idx <- svi_pct %>%
  group_by(county) %>%
  slice_sample(prop = 0.7) %>%
  ungroup() %>%
  pull(GEOID)

train <- svi_pct %>% filter(GEOID %in% train_idx)
test  <- svi_pct %>% filter(!GEOID %in% train_idx)

cat("Train:", nrow(train), "| Test:", nrow(test), "\n")

pct_vars <- c("p_pov150","p_unemp","p_hburd","p_nohsdp","p_uninsur",
              "p_age65","p_age17","p_disabl","p_sngpnt","p_limeng",
              "p_minrty","p_munit","p_mobile","p_crowd","p_noveh","p_groupq")

ecdf_train <- lapply(pct_vars, function(v) ecdf(train[[v]]))
names(ecdf_train) <- pct_vars

add_ranks <- function(df) {
  df %>%
    mutate(across(all_of(pct_vars),
                  ~ ecdf_train[[cur_column()]](.x),
                  .names = "r_{.col}"))
}

train_r <- add_ranks(train)
test_r  <- add_ranks(test)

# verify
cat("Train rank range (should be ~[0,1]):\n"); print(summary(train_r$r_p_pov150))
cat("Test rank range (may be slightly outside if distributions differ):\n"); print(summary(test_r$r_p_pov150))
cat("Train rows:", nrow(train_r), "| Test rows:", nrow(test_r), "\n")
cat("Any NA ranks in test:", sum(!complete.cases(select(test_r, starts_with("r_")))), "\n")

saveRDS(train_r, "train_r.rds")
saveRDS(test_r,  "test_r.rds")
