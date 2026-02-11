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
# Option 1 (recommended): pak
install.packages("pak")
pak::pak("montilab/ClockSHAP")

# Option 2: remotes
install.packages("remotes")
remotes::install_github("montilab/ClockSHAP")

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

Aging-clock input features often change systematically with chronological age. A single global baseline (e.g., the unconditional mean feature vector) can be inappropriate for comparing younger vs. older samples.

ClockSHAP therefore uses an **age-conditioned reference profile**. For each feature `k`, we estimate an age trend in a reference cohort:

```text
x_ref[k](age) = gamma0[k] + gamma1[k] * age
```

where `gamma0[k]` and `gamma1[k]` are stored in the fitted `reference_profile` object.

Given a clock that maps features `x` to a predicted age `y_hat = f(x)`, ClockSHAP defines:

```text
y_exp(age)  = f( x_ref(age) )
deviation   = y_hat - y_exp(age)
```

Interpretation:
- `expected` is the clock’s **typical output** for reference individuals of the **same chronological age**.
- `deviation` is **how far above/below that age-matched expectation** a sample lies.

In this package, `reference_profile` is fit as per-feature linear models by default (the `gamma0/gamma1` parameters above). You can treat `y_exp(age)` as an age-matched baseline defined through the feature space rather than a direct regression on predicted values.


## Terminology (delta age, age acceleration, and what ClockSHAP explains)

Below are three quantities that are often discussed together but are **not the same**. ClockSHAP is built to explain **deviation from an age-matched reference expectation**, not raw prediction error.

### 1) Absolute age acceleration (predicted − chronological age) (sometimes refered to as delta age)

```text
AAA = y_hat - age
```

This is simple, but it can be **age-biased** when the clock has regression-to-the-mean (e.g., slope < 1), leading to younger samples being over-predicted and older samples under-predicted.

### 2) Residual age acceleration (RAA; regression residuals)

A common fix is to regress clock output on age in a chosen cohort and use residuals:

```r
RAA <- residuals(stats::lm(predicted ~ age))
```

Equivalently, if the fitted line in that cohort is:

```text
y_hat ≈ beta0 + beta1 * age
```

then:

```text
RAA = y_hat - (beta0 + beta1 * age)
```

This is a **post hoc calibration** intended to remove systematic age trends in prediction error.

### 3) ClockSHAP deviation (age-matched deviation; predicted − expected)

ClockSHAP uses the **age-conditioned** target:

```text
deviation = y_hat - y_exp(age)
```

where `y_exp(age)` (stored as `expected`) is computed by applying the *same clock* to an **age-conditioned reference profile** (see “How expected is defined” below).

Why this often reduces age bias without an extra regression step:
- If the clock tends to over-predict at younger ages and under-predict at older ages, that pattern can appear in both `y_hat` and `y_exp(age)` (assuming upstream harmonization / normalization).
- Subtracting them targets deviation from an **age-matched expectation** by construction.

Relationship to residual-based acceleration (RAA):
- If `y_exp(age)` approximates `E[y_hat | age]` for a cohort similar to what you would use for `lm(predicted ~ age)`, then `deviation` and RAA can be similar numerically.
- They are not guaranteed to be identical in general, because `y_exp(age)` is defined through a **feature-level** reference profile, not directly through a regression on predicted values.
- In the common case where (i) the clock is linear-additive and (ii) the reference profile is linear in age per feature (the default in this package), `y_exp(age)` becomes a linear function of age and `deviation` can coincide with reference-cohort residualization.

**Naming note:** we recommend calling this quantity **“deviation”** or **“age-matched deviation”** (and using ΔAge / RAA only when you explicitly mean those definitions).


## Why we call it “SHAP”

For **linear-additive models**, SHAP-style additive attributions have a closed-form expression.

Write the clock on standardized features:

```text
y_hat = alpha + sum_k beta[k] * z[k]
z[k]  = (x[k] - mu[k]) / sigma[k]
```

Given a baseline standardized feature vector `z_ref`, linear-model SHAP values reduce to:

```text
phi[k] = beta[k] * ( z[k] - z_ref[k] )
```

and additivity holds exactly:

```text
sum_k phi[k] = y_hat(x) - y_hat(z_ref)
```

### The ClockSHAP twist: the baseline depends on age

In aging clocks, the most meaningful baseline depends on chronological age. ClockSHAP sets the baseline using the reference profile:

```text
x_ref[k](age) = gamma0[k] + gamma1[k] * age
z_ref[k](age) = ( x_ref[k](age) - mu[k] ) / sigma[k]
```

So ClockSHAP uses:

```text
y_exp(age)   = alpha + sum_k beta[k] * z_ref[k](age)
phi[k](age)  = beta[k] * ( z[k] - z_ref[k](age) )
deviation    = y_hat - y_exp(age)
sum_k phi[k](age) = deviation
```

This preserves exact additivity while making the baseline biologically appropriate for aging-clock interpretation.

**API mapping:** `phi` is returned by `clockshap_phi(cs)`, and `gamma0/gamma1` are stored in the fitted `reference_profile` returned by `fit_reference_profile()`.


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
