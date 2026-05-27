if (!requireNamespace("testthat", quietly = TRUE)) {
  quit(save = "no", status = 0)
}

library(testthat)
library(inferMM)

test_check("inferMM")
