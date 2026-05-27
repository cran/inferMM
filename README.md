# inferMM

`inferMM` provides variance-aware Michaelis-Menten estimation and inference for
enzyme-kinetic data with concentration-dependent heteroscedasticity.

The package is designed around a compact workflow:

- fit one curve with `fit_mm()`
- compare working variance models with `screen_mm()`
- repeat that workflow across many enzymes with `group_mm()`
- fit repeated or clustered assays with `cluster_mm()`
- summarize and plot fitted objects with standard S3 methods plus `report_mm()`

## Installation

```r
# install.packages("remotes")
remotes::install_github("mijeong-kim/inferMM")
```

## Bundled data

The package ships with two demo datasets.

- `sdl_demo`: self-driving laboratory Michaelis-Menten panel
- `alves_demo`: filtered soil exoenzyme kinetics panel from Alves et al. (2021)

```r
library(inferMM)
data(sdl_demo)
data(alves_demo)
```

## Minimal example

```r
library(inferMM)

one_curve <- subset(sdl_demo, enzyme == "1111")

fit <- fit_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  variance = "sqrt"
)

summary(fit)
confint(fit)
```

## Screening working variance models

```r
screen <- screen_mm(
  x = one_curve$s_uM,
  y = one_curve$v_uM_per_min,
  quiet = TRUE
)

screen$table[, c("model", "selected_model", "quasi_aic", "quasi_bic", "rmse")]
```

## Grouped analyses

```r
grouped <- group_mm(
  data = sdl_demo,
  s = "s_uM",
  v = "v_uM_per_min",
  groups = "enzyme",
  variance_models = c("constant", "log", "sqrt", "cuberoot"),
  quiet = TRUE
)

grouped$comparison$best_by_group[
  , c("group_label", "model", "selected_model", "quasi_aic", "quasi_bic", "rmse")
]
```

## Clustered analyses

```r
cluster_fit <- cluster_mm(
  data = subset(alves_demo, enzyme == "BG"),
  s = "substrate_conc",
  v = "activity",
  cluster = "core",
  variance = "sqrt"
)

summary(cluster_fit)
confint(cluster_fit)
```

For sparse clustered fits, default interval reporting is intentionally cautious:
printed summaries may suppress intervals, and bootstrap intervals should be read
as sensitivity analyses rather than routine default inference.

## Reporting and plotting

```r
report_mm(fit, interval_type = "confidence")
plot(grouped, interval_type = "confidence")
predict(fit, newdata = seq(0, 80, length.out = 6), interval = "prediction")
```

## Repository contents

- `R/`: package source code
- `man/`: function documentation
- `data/`: bundled `.rda` demo data
- `inst/extdata/`: raw CSV mirrors of the demo data
- `vignettes/`: end-to-end workflow vignette
- `tests/`: unit tests

For manuscript-oriented simulation code and saved paper outputs, see the
separate repository `inferMM-cils-repro`.
