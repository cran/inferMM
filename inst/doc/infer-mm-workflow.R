## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(inferMM)

one_curve <- subset(sdl_demo, enzyme == "1111")

fit <- fit_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  variance = "sqrt"
)

summary(fit)

fit_auto <- fit_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  variance = "auto",
  power_selection = "quasi_aic"
)

summary(fit_auto)

## -----------------------------------------------------------------------------
coef(fit)
confint(fit)
set.seed(1)
confint(fit, method = "bootstrap", B = 99)
head(predict(fit, newdata = seq(0, 80, length.out = 6), interval = "prediction"))

## ----eval = FALSE-------------------------------------------------------------
# report_mm(fit, interval_type = "confidence")
# set.seed(1)
# report_mm(fit, method = "bootstrap", B = 99,
#           interval_type = "confidence")

## -----------------------------------------------------------------------------
screen <- screen_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  power_values = c(0.4, 0.6),
  include_auto = TRUE,
  quiet = TRUE
)

screen$table[, c("model", "selected_model", "quasi_aic", "quasi_bic", "rmse")]

## -----------------------------------------------------------------------------
h <- variance_function(function(S) 1 + S / 20,
                       label = "Var = gamma (1 + S/20)")
custom_fit <- fit_mm(one_curve$s_uM, one_curve$v_uM_per_min, variance = h)
summary(custom_fit)
confint(custom_fit)
predict(custom_fit, newdata = c(0, 10, 40), interval = "prediction")

custom_shapes <- list(
  linear_floor = function(S) 1 + S / 20,
  saturation = function(S) 1 + 9 * S / (20 + S)
)
custom_screen <- screen_mm(
  one_curve$s_uM, one_curve$v_uM_per_min,
  variance_models = c("constant", "sqrt"),
  variance_functions = custom_shapes, quiet = TRUE
)
custom_screen$table[, c("model", "quasi_aic", "quasi_bic", "success")]

## -----------------------------------------------------------------------------
set.seed(19)
skew_dat <- simulate_mm_data(n = 200, error = "skewed", skew_shape = 4)
assay_shape <- function(S) 1 + 9 * S / (20 + S)
fit_hm <- fit_mm_moments(skew_dat$x, skew_dat$y, variance = assay_shape)
summary(fit_hm)
confint(fit_hm)
fit_hm$vcov_joint

fit_gamma <- fit_mm_moments(
  skew_dat$x, skew_dat$y, variance = assay_shape,
  moments = "gamma", shape = 4
)
coef(fit_gamma)
predict(fit_gamma, c(10, 20, 80), interval = "confidence")

## ----eval = FALSE-------------------------------------------------------------
# set.seed(20)
# confint(fit_hm, method = "bootstrap", B = 399)  # pairs only
# report_mm(fit_hm, interval_type = "confidence")

## -----------------------------------------------------------------------------
grouped <- group_mm(
  data = sdl_demo,
  s = "s_uM",
  v = "v_uM_per_min",
  groups = "enzyme",
  variance_models = c("constant", "log", "sqrt", "cuberoot"),
  power_values = c(0.4, 0.6),
  include_auto = TRUE,
  variance_functions = custom_shapes,
  quiet = TRUE
)

head(grouped$comparison$best_by_group[, c("group_label", "model", "selected_model", "quasi_aic", "quasi_bic", "rmse")])

## ----eval = FALSE-------------------------------------------------------------
# plot(grouped, interval_type = "confidence")

## -----------------------------------------------------------------------------
cluster_fit <- cluster_mm(
  data = subset(alves_demo, enzyme == "BG"),
  s = "substrate_conc",
  v = "activity",
  cluster = "core",
  variance = "sqrt"
)

summary(cluster_fit)
confint(cluster_fit)
head(predict(cluster_fit, newdata = seq(0, 700, length.out = 6),
             interval = "confidence"))

## -----------------------------------------------------------------------------
set.seed(100)
sim_dat <- simulate_mm_data(variance_shape = "hill", error = "skewed")
head(sim_dat)

