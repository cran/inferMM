normalize_variance_model <- function(model) {
  aliases <- c(
    x_half = "sqrt",
    x_third = "cuberoot",
    cube_root = "cuberoot",
    power_auto = "auto",
    auto_power = "auto"
  )

  if (model %in% names(aliases)) {
    return(aliases[[model]])
  }

  model
}

default_power_grid <- function() {
  sort(unique(c(seq(0.2, 1.0, by = 0.05), 1 / 3, 1 / 2)))
}

validate_power_grid <- function(power_grid) {
  out <- sort(unique(as.numeric(power_grid)))
  keep <- is.finite(out) & out > 0
  out <- out[keep]
  if (!length(out)) {
    stop("power_grid must contain at least one positive finite value.")
  }
  out
}

format_power_value <- function(power, digits = 2L) {
  txt <- formatC(as.numeric(power), format = "f", digits = digits)
  txt <- sub("0+$", "", txt)
  txt <- sub("\\.$", "", txt)
  txt
}

selection_label <- function(selection) {
  switch(
    selection,
    quasi_aic = "quasi-AIC",
    quasi_bic = "quasi-BIC",
    selection
  )
}

power_model_name <- function(power) {
  paste0("power(", format_power_value(power), ")")
}

auto_model_name <- function(selection) {
  paste0("auto(", selection, ")")
}

power_label <- function(power) {
  paste0("Var = gamma S^", format_power_value(power))
}

selected_variance_model_name <- function(model, power = NULL) {
  if (identical(model, "power")) {
    return(power_model_name(power))
  }
  model
}

selected_variance_label <- function(model, power = NULL) {
  if (identical(model, "power")) {
    return(power_label(power))
  }
  model_display_label(model)
}

auto_selected_label <- function(label, selection) {
  paste0(label, " (selected by ", selection_label(selection), " from log/power candidates)")
}

resolve_variance_spec <- function(variance = "log",
                                  power = NULL) {
  variance <- normalize_variance_model(variance)
  fun <- variance_function(model = variance, power = power)

  list(
    model = attr(fun, "variance_model", exact = TRUE),
    label = attr(fun, "variance_label", exact = TRUE),
    fun = fun
  )
}

prepare_variance_candidates <- function(variance_models = default_variance_models(),
                                        power_values = NULL,
                                        include_auto = FALSE,
                                        power_selection = c("quasi_aic", "quasi_bic"),
                                        power_grid = default_power_grid()) {
  power_selection <- match.arg(power_selection)
  power_grid <- validate_power_grid(power_grid)

  models <- unique(vapply(variance_models, normalize_variance_model, character(1)))
  if ("power" %in% models && is.null(power_values)) {
    stop("Use 'power_values' to add fixed power models to screen_mm() or group_mm().")
  }
  if ("auto" %in% models) {
    include_auto <- TRUE
  }
  models <- setdiff(models, c("power", "auto"))

  bad <- setdiff(models, c("constant", "log", "sqrt", "cuberoot"))
  if (length(bad)) {
    stop("Unsupported variance model(s): ", paste(bad, collapse = ", "))
  }

  candidates <- lapply(models, function(model) {
    list(
      id = model,
      label = model_display_label(model),
      fit_args = list(variance = model)
    )
  })

  if (!is.null(power_values)) {
    power_values <- validate_power_grid(power_values)
    candidates <- c(
      candidates,
      lapply(power_values, function(power) {
        list(
          id = power_model_name(power),
          label = power_label(power),
          fit_args = list(variance = "power", power = power)
        )
      })
    )
  }

  if (isTRUE(include_auto)) {
    candidates <- c(
      candidates,
      list(
        list(
          id = auto_model_name(power_selection),
          label = paste0("Auto: best of log/power by ", selection_label(power_selection)),
          fit_args = list(
            variance = "auto",
            power_selection = power_selection,
            power_grid = power_grid
          )
        )
      )
    )
  }

  keep <- !duplicated(vapply(candidates, `[[`, character(1), "id"))
  candidates[keep]
}

#' Create Working Variance Functions
#'
#' Construct a working variance function for Michaelis-Menten models.
#'
#' @param model Character string naming the working variance specification.
#' @param power Optional exponent used when `model = "power"`.
#'
#' @return A function of substrate concentration `x`.
#' @export
variance_function <- function(model = c("constant", "log", "sqrt", "cuberoot", "power"),
                              power = NULL) {
  model <- normalize_variance_model(match.arg(model))

  out <- switch(
    model,
    constant = function(x) rep(1, length(x)),
    log = function(x) log(x + 1),
    sqrt = function(x) x^(1 / 2),
    cuberoot = function(x) x^(1 / 3),
    power = {
      if (is.null(power) || !is.finite(power) || length(power) != 1L || power <= 0) {
        stop("A positive finite scalar 'power' must be supplied when model = 'power'.")
      }
      force(power)
      function(x) x^power
    }
  )

  attr(out, "variance_model") <- model
  attr(out, "variance_label") <- if (identical(model, "power")) {
    power_label(power)
  } else {
    model_display_label(model)
  }
  attr(out, "variance_power") <- if (identical(model, "power")) power else NULL
  out
}
