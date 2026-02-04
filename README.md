# ClockSHAP

ClockSHAP is an interpretation framework for *linear-additive aging clocks* (e.g., transcriptomic, epigenetic, or proteomic clocks fit using linear regression, ridge regression, elastic net, etc.) 
that explains why a sample is predicted to be biologically older or younger than expected given its chronological age.
The method was developed in the context of tissue-anchored transcriptomic aging clocks, but is designed to support comparative aging analyses across diverse clock contexts.

Traditional aging clocks typically answer the question:

> “Is this sample biologically older or younger than expected?”

ClockSHAP goes one step further and asks:

> “*Why* is this sample predicted to be biologically older or younger than expected for its age?””

It does this by decomposing **age deviation** (also called “delta age” or “age acceleration”; predicted age − expected age) into **age-adjusted, per-feature contributions** that add up exactly to the deviation.

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

## Why ClockSHAP does not use “relative age acceleration” (RAA)

Many aging-clock analyses define **age acceleration** as the difference between a clock’s predicted age and chronological age:

\[
\text{Age acceleration} = \hat{y} - \text{age}
\]

In practice, most clocks exhibit **regression-to-the-mean** / miscalibration (e.g., slope < 1), so this quantity can be **age-biased**:
- younger individuals tend to be predicted *too old*
- older individuals tend to be predicted *too young*

To address this, it is common to compute **relative age acceleration (RAA)** by removing systematic age trends in predicted age, e.g. by fitting:

\[
\hat{y} = a + b\cdot \text{age} + \varepsilon
\]

and using the residuals (or an equivalent bias-corrected form) as the “relative” acceleration:

\[
\text{RAA} = \hat{\varepsilon}
\]

RAA is therefore a **post hoc calibration step** that makes \(\hat{y} - \text{age}\) more comparable across ages by correcting age-dependent prediction bias.

### ClockSHAP’s alternative: deviation relative to an age-conditioned expectation

ClockSHAP uses a different estimand. Rather than comparing predicted age to chronological age, ClockSHAP defines **deviation** as:

\[
\text{Deviation} = \hat{y} - E[\hat{y}\mid \text{age}]
\]

where \(E[\hat{y}\mid \text{age}]\) is the clock’s **average prediction** for reference individuals at the same chronological age (estimated from a user-supplied reference cohort).

Because both **predicted** and **expected** ages are produced by the **same clock**, they inherit the same global age-dependent miscalibration (e.g., slope compression). As a result, systematic age bias cancels **by construction**:

- if the clock tends to predict younger samples too old and older samples too young, that behavior is present in both \(\hat{y}\) and \(E[\hat{y}\mid \text{age}]\)
- subtracting them removes this global age trend without requiring an additional RAA step

In short: **ClockSHAP does not require RAA** because deviation is already an internally age-adjusted quantity.

### Caveats and practical requirements

ClockSHAP deviations should be interpreted as:

> **relative deviation within a specific clock and reference context**, not an absolute measure of biological age.

Key implications:

- **“Years” are clock-relative, not absolute.**  
  A deviation of +5 years under one clock is not guaranteed to be comparable to +5 years under a different clock, a different training procedure, or a different feature space.

- **Reference population matters.**  
  Expected age \(E[\hat{y}\mid \text{age}]\) is defined by your chosen reference cohort. For stable, interpretable deviations, the reference should:
  - have a **sufficient age range** (and ideally substantial overlap with the age range of the samples of interest)
  - be large enough to estimate a smooth age-conditioned expectation
  - be biologically relevant to the intended comparison (e.g., tissue/context matched when appropriate)

- **Preprocessing and harmonization must be done upstream.**  
  ClockSHAP assumes the feature matrix is already “analysis-ready.” Users should handle *before* running ClockSHAP:
  - feature preprocessing/normalization consistent with the clock (e.g., transformations, scaling/standardization inputs)
  - feature-space harmonization (same features, same units/definitions, consistent ordering)
  - batch effect correction / technical confounding control, where applicable
  - missingness and QC filtering consistent across cohorts

ClockSHAP is an interpretability framework: it explains deviations given a clock and a reference, but it does **not** attempt to fix calibration, batch, or distribution-shift issues internally. 
Failure to address these upstream may lead to misleading results.

For more discussion on the idea of relative-age-acceleration, potential issues of mismatched feature space/batch effects, and related topics, see [Epigenetic ageing clocks: statistical methods and emerging computational challenges](https://doi.org/10.1038/s41576-024-00807-w) and references therein.

## Status

This package is under active development.
The public API may change prior to a stable release.

## Citation

Manuscript in preparation.

For manuscript submission, we recommend citing a **tagged GitHub release** and (optionally) an archived DOI (e.g., Zenodo) corresponding to the exact version used in the paper.
