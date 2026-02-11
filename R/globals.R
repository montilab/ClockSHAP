# Silence R CMD check NOTES for non-standard evaluation symbols used in
# dplyr/ggplot2 pipelines inside plotting code.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "Feature", "Eff", "Sign", "IsOther",
      "Start", "End", "y", "fill_key", "Mid"
    )
  )
}
