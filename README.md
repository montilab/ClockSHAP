# ClockSHAP

ClockSHAP is an interpretation framework for **linear-additive aging clocks** (e.g., transcriptomic, epigenetic, proteomic clocks) that explains *why* a sample is predicted to be biologically older or younger than expected given its chronological age.

The method was developed in the context of tissue-anchored transcriptomic aging clocks, but is designed to support comparative aging analyses across diverse clock contexts.

Instead of asking:

> “Which features increase or decrease predicted age?”

ClockSHAP asks:

> “Which features make this sample *older or younger than expected for its age*?”

It does this by decomposing **age deviation** (also called “delta age” or “age acceleration”; predicted − expected) into **age-adjusted, per-feature contributions** that add up exactly to the deviation.

## What ClockSHAP returns

For each sample, ClockSHAP computes:

- **Predicted age**: the clock’s predicted biological age.
- **Expected age**: the clock’s prediction for a *reference* individual at the same chronological age.
- **Deviation**: `predicted - expected`.
- **Per-feature contributions** (`phi`): a matrix of feature attributions whose row-sums equal the deviation.

## Installation

ClockSHAP is currently distributed via GitHub.

```r
# install.packages("pak")
pak::pak("<org-or-user>/ClockSHAP")
```

## Quick start

```r
library(ClockSHAP)

set.seed(1)

# Synthetic data: samples x features
features <- matrix(rnorm(20), nrow = 5, ncol = 4)
colnames(features) <- paste0("F", 1:4)
age <- c(40, 50, 60, 70, 80)

# A linear clock specification:
#   predicted_age = alpha + sum_k beta_k * z_k
# where z_k = (x_k - mu_k) / sigma_k
clock <- linear_clock(
  alpha = 10,
  beta  = setNames(runif(4), colnames(features)),
  mu    = setNames(rep(0, 4), colnames(features)),
  sigma = setNames(rep(1, 4), colnames(features))
)

# Fit an age-feature reference profile from a cohort (per-feature linear models)
ref <- fit_reference_profile(
  features = as.data.frame(features),
  age = age
)

# Decompose deviation into feature contributions
cs <- clockshap(features = features, age = age, clock = clock, reference = ref)

cs
summary(cs)

# Extract components
phi       <- clockshap_phi(cs)
predicted <- clockshap_predicted(cs)
expected  <- clockshap_expected(cs)
deviation <- clockshap_deviation(cs)

# Visualize one sample
plot_clockshap_waterfall(cs, sample = 1, top_n = 10)
```

## The reference profile (why “expected” is age-adjusted)

Aging-clock input features often change systematically with chronological age.
If we use a single, global baseline (e.g., the unconditional mean feature vector), then the “baseline” can be inappropriate for older vs. younger individuals.

ClockSHAP addresses this by defining a **reference profile**:

$$E[X_k \mid age] = \gamma_{0k} + \gamma_{1k}\cdot age$$

The “expected age” is the clock’s prediction evaluated at the age-conditioned expected feature values. The resulting **deviation** and **feature attributions** are therefore interpretable as:

> effects that explain *relative* biological aging (older/younger than expected for chronological age),
> not just absolute predicted age.

## Why we call it “SHAP”

ClockSHAP is called “SHAP” because, for **linear-additive models**, SHAP values have a **closed-form solution**.

Consider a standardized linear clock:

$$\hat{y}(x) = \alpha + \sum_k \beta_k\, z_k$$
where \(z_k = (x_k - \mu_k)/\sigma_k\).

If we choose a baseline feature vector \(z^*\), then the SHAP value for feature \(k\) in a linear model (under the common *interventional / additive* formulation) reduces to:

$$\phi_k = \beta_k\,(z_k - z_k^*)$$

and the attributions add up exactly:

$$\sum_k \phi_k = \hat{y}(x) - \hat{y}(z^*)$$

### The ClockSHAP twist: an age-conditioned baseline

In aging clocks, the most meaningful baseline depends on chronological age.
ClockSHAP therefore uses an **age-conditioned baseline** derived from the reference profile:

$$x_k^*(age) = \gamma_{0k} + \gamma_{1k}\cdot age$$
$$z_k^*(age) = (x_k^*(age) - \mu_k)/\sigma_k$$

ClockSHAP then defines:

- **Expected age**: \(\hat{y}(z^*(age))\)
- **Deviation**: \(\hat{y}(z) - \hat{y}(z^*(age))\)
- **Feature contributions**: \(\phi_k = \beta_k\,(z_k - z_k^*(age))\)

This preserves the desirable SHAP property (exact additivity) while aligning the baseline with the biology of aging clocks.

## Status

This package is under active development.
The public API may change prior to a stable release.

## Citation

Manuscript in preparation.

For manuscript submission, we recommend citing a **tagged GitHub release** and (optionally) an archived DOI (e.g., Zenodo) corresponding to the exact version used in the paper.
