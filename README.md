# ClockSHAP

ClockSHAP is an interpretability framework for **linear-additive aging clocks** (e.g., transcriptomic, epigenetic, proteomic clocks fit using linear regression, ridge regression, elastic net, etc.).

It answers the question:

> **Why is this sample predicted to be biologically older or younger than expected for its chronological age?**

The method was developed in the context of tissue-anchored transcriptomic aging clocks, but is designed to support comparative aging analyses across diverse clock contexts.

## The key point (terminology matters)

Aging-clock papers often use terms like **delta age** and **age acceleration**, sometimes with different meanings.
ClockSHAP is built around a specific, age-conditioned target:

- **Predicted (clock) age**: `predicted`
- **Expected (clock) age for the same chronological age**: `expected`
- **Age-matched deviation** (package output: `deviation`):

`deviation = predicted - expected`

ClockSHAP explains `deviation` by decomposing it into per-feature contributions (`phi`) that add up exactly:

`rowSums(phi) = deviation`

> ClockSHAP explains **predicted − expected**, not **predicted − chronological age**.

If you want common alternatives, see **“Terminology”** below.

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
# predicted = alpha + sum_k beta[k] * z[k]
# where z[k] = (x[k] - mu[k]) / sigma[k]
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

# Decompose age-matched deviation into feature contributions
cs <- clockshap(features = features, age = age, clock = clock, reference = ref)

# Extract components
phi       <- clockshap_phi(cs)
predicted <- clockshap_predicted(cs)
expected  <- clockshap_expected(cs)
deviation <- clockshap_deviation(cs)

# Visualize one sample
plot_clockshap_waterfall(cs, sample = 1, top_n = 10)
```

## How “expected” is defined (the reference profile)

Aging-clock input features often change systematically with chronological age.
A single global baseline (e.g., the unconditional mean feature vector) can be inappropriate for comparing younger vs. older samples.

ClockSHAP therefore uses an **age-conditioned reference profile**. For each feature `k`, we estimate an age trend:

```text
x_ref[k](age) = gamma0[k] + gamma1[k] * age
```

Then we compute:

```text
expected = clock( x_ref(age) )
deviation = predicted - expected
```

Interpretation:
- `expected` is the clock’s **typical output** for reference individuals of the **same chronological age**.
- `deviation` is **how far above/below that age-matched expectation** a sample lies.

## Terminology (delta age, age acceleration, and what ClockSHAP explains)

Below are three quantities that are often discussed together but are *not the same*.

### 1) Delta age (predicted − chronological age)

```text
delta_age = predicted - age
```

This is simple, but it can be **age-biased** when the clock is miscalibrated (e.g., slope < 1).

### 2) Residual age acceleration / “relative age acceleration” (RAA)

A common fix is to regress clock output on age in some cohort and use residuals:

```r
age_accel <- residuals(stats::lm(predicted ~ age))
```

This is a **post hoc calibration** step intended to remove systematic age trends in prediction error.

### 3) ClockSHAP deviation (predicted − expected)

ClockSHAP uses:

```text
deviation = predicted - expected
```

where `expected` is computed by applying the **same clock** to an **age-conditioned reference profile**.

Why this often removes age bias without an extra regression step:
- if the clock tends to over-predict at younger ages and under-predict at older ages, that pattern appears in both `predicted` and `expected`, assuming upstream feature harmonization/ data normalization.
- subtracting them removes this global age trend **by construction**

Relationship to residual-based acceleration:
- If your reference profile is estimated from a cohort similar to what you would use for `lm(predicted ~ age)`, then `expected(age)` will often approximate `E[predicted | age]`.
- In that case, `deviation` and residual-based acceleration can be similar numerically.
- They are not guaranteed to be identical, because `expected(age)` is defined through a **feature-level** reference profile, not directly through a regression on `predicted`.

## Why we call it “SHAP”

For **linear-additive models**, SHAP-style additive attributions have a closed-form expression.

Write the clock on standardized features:

```text
predicted = alpha + sum_k beta[k] * z[k]
where z[k] = (x[k] - mu[k]) / sigma[k]
```

Given a baseline standardized feature vector `z_ref`, linear-model SHAP values reduce to:

```text
phi[k] = beta[k] * ( z[k] - z_ref[k] )
```

and additivity holds exactly:

```text
sum_k phi[k] = predicted(x) - predicted(z_ref)
```

### The ClockSHAP twist: the baseline depends on age

In aging clocks, the most meaningful baseline depends on chronological age.
ClockSHAP sets the baseline using the reference profile:

```text
x_ref[k](age) = gamma0[k] + gamma1[k] * age
z_ref[k](age) = ( x_ref[k](age) - mu[k] ) / sigma[k]
```

So ClockSHAP uses:

```text
expected(age)  = alpha + sum_k beta[k] * z_ref[k](age)
phi[k]         = beta[k] * ( z[k] - z_ref[k](age) )
deviation      = predicted - expected(age)
sum_k phi[k]   = deviation
```

This preserves exact additivity while making the baseline biologically appropriate for aging-clock interpretation.

## What ClockSHAP returns

For each sample:

- `predicted`: clock predicted age
- `expected`: clock expected age at the same chronological age (from the reference profile)
- `deviation`: `predicted - expected` (age-matched deviation)
- `phi`: per-feature contributions that sum to `deviation`

Accessors:

- `clockshap_predicted(cs)`
- `clockshap_expected(cs)`
- `clockshap_deviation(cs)`
- `clockshap_phi(cs)`

## Practical notes

- **Deviation is reference-dependent.** Your choice of reference cohort/profile defines what “expected for age” means.
- **“Years” are clock-relative.** A deviation of +5 under one clock is not guaranteed comparable to +5 under a different clock, training set, or feature space.
- **Upstream preprocessing matters.** ClockSHAP assumes the feature matrix is already analysis-ready and consistent with the clock (transformations, normalization, feature matching, batch correction, QC).
The *importance of this step can not be overstated*!! Mismatching feature spaces, normalization steps, or batch effects can lead to technical rather than biological deviations. 

For broader discussion of calibration, age bias, and computational challenges in clock analysis, see:

- Epigenetic ageing clocks: statistical methods and emerging computational challenges (Nature Reviews Genetics). https://doi.org/10.1038/s41576-024-00807-w

## Status

This package is under active development.
