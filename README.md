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
- Because `expected` comes from a chosen reference (e.g., GTEx tissue-of-origin), `deviation` should be interpreted **relative to that reference**.


In this package, `reference_profile` is fit as per-feature linear models by default (the `gamma0/gamma1` parameters above). You can treat `y_exp(age)` as an age-matched baseline defined through the feature space rather than a direct regression on predicted values.


## Terminology (delta age, age acceleration, and what ClockSHAP explains)

Below are three quantities that are often discussed together but are **not the same**. ClockSHAP is built to explain **deviation from an age-matched reference expectation**, not raw prediction error.

### 1) Absolute age acceleration (AAA; predicted − chronological age)

```text
AAA = y_hat - age
```

This is simple, but it can be **age-biased** when the clock has regression-to-the-mean (e.g., slope < 1), leading to younger samples being over-predicted and older samples under-predicted.

NOTE: AAA is sometimes referred to as delta(Δ) Age; we use AAA here for parallelism with RAA.


### 2) Relative age acceleration (RAA; regression residuals)

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

Naming map (ClockSHAP object fields):
```text
predicted = y_hat
expected  = y_exp(age)
deviation = predicted - expected
phi       = per-feature contributions (summing to deviation)
```

Key invariant (per sample i):
```text
sum_k phi_i[k] = deviation_i
```

where `y_exp(age)` (stored as `expected`) is computed by applying the *same clock* to an **age-conditioned reference profile** (see “How expected is defined”).

#### Why this can reduce age bias without an extra regression step

- `deviation = y_hat − y_exp(age)` subtracts an **age-matched baseline** instead of subtracting chronological age directly. This avoids the classic regression-to-the-mean artifact that makes `AAA = y_hat − age` systematically positive at young ages and negative at old ages when the clock slope is < 1.
- In practice, this helps when `y_exp(age)` tracks the **typical clock output at that age** in the reference setting (i.e., it behaves like an age-conditional expectation for the reference). If cross-cohort effects or biology cause the clock’s age-trend to differ strongly between the reference and the target, `deviation` may still show residual age trends (which should be interpreted as reference mismatch and/or biological decoupling, depending on context).

#### Relationship to relative age acceleration (RAA)

- RAA is often defined as residuals from `lm(y_hat ~ age)` in some cohort. ClockSHAP’s `deviation` is similar in spirit: it measures how far a sample’s clock output lies above/below the **age-conditional expectation**, except that ClockSHAP defines that expectation through a **feature-level reference profile** rather than a regression on `y_hat`.
- Because of this difference in how the expectation is constructed, `deviation` and regression residuals are not guaranteed to match in all settings. They will be numerically close when `y_exp(age)` closely approximates the cohort’s age-conditional mean clock output.
- In the default setup used here (linear-additive clock; reference profile fit as per-feature linear trends `gamma0 + gamma1 * age`), `y_exp(age)` is linear in age, so `deviation` often matches the residual from a reference-cohort fit of `y_hat ~ age` (as we observed in our GTEx→TCGA application).


**Reference note:** `expected` (and thus `deviation`) is **reference-anchored** (e.g., a specific GTEx tissue-of-origin). This quantity is not a universal biological-age score.

**Interpretation note:** `deviation` is a cross-sectional **state shift** relative to an age-matched reference expectation; it does not by itself imply a within-individual aging *rate*.

**Naming note:** we recommend calling this quantity **“deviation”** or **“age-matched deviation”** (and using AAA / RAA only when you explicitly mean those definitions).



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

Note: `phi` values are defined relative to the age-conditioned baseline (via `gamma0/gamma1`); changing the reference changes the decomposition target.


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
