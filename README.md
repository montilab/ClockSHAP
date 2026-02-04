# ClockSHAP

ClockSHAP is an interpretation framework for **linear-additive aging clocks** (e.g., transcriptomic, epigenetic, proteomic clocks) that explains *why* a sample is predicted to be biologically older or younger **relative to what is expected for its chronological age**.

The method was developed in the context of tissue-anchored transcriptomic aging clocks, but is designed to support comparative aging analyses across diverse clock contexts.

Instead of asking:

> “Which features increase or decrease predicted age?”

ClockSHAP asks:

> “Which features make this sample *older or younger than expected for its age*?”

**Key point:** ClockSHAP explains **predicted − expected**, **not** predicted − chronological.

---

## What ClockSHAP explains (terminology you can cite)

Aging-clock papers often use terms like **“delta age”** and **“age acceleration”** inconsistently.  
To avoid ambiguity, this package uses the following definitions:

- **Chronological age** (\(\tau\)): the observed/known age of the sample (input `age`).
- **Predicted age** (\(\hat{y}\)): the clock’s output for that sample (output `predicted`).

Two common derived quantities in the literature:

- **Delta age** (sometimes “age difference”):  
  $$\Delta = \hat{y} - \tau$$  
  i.e., **predicted − chronological**.

- **Age acceleration** (common operational definition):  
  $$AA = \hat{y} - E[\hat{y} \mid \tau]$$  
  i.e., **predicted − age-conditioned expected prediction**, where the expectation is often estimated by regressing predicted age on chronological age in a reference cohort.

### ClockSHAP’s target: *age-matched deviation* (predicted − expected)

ClockSHAP does **not** explain \(\hat{y}-\tau\) directly.  
Instead, ClockSHAP explains an age-conditioned deviation that we call **age-matched deviation**:

$$D(\tau) = \hat{y} - \hat{y}^{\*}(\tau)$$

where \(\hat{y}^{\*}(\tau)\) (output `expected`) is the clock’s prediction for a **reference individual at the same chronological age**.  
That reference prediction is computed using an *age-conditioned baseline in feature space* (see “The reference profile” below).

**In the package,** `deviation` **means** this *age-matched deviation* \(D(\tau)\).

> Practical interpretation: ClockSHAP tells you which features push a sample **older/younger than expected among peers of the same age** (under the chosen reference profile).

**Relationship to “age acceleration”:** \(D(\tau)\) is conceptually closest to an age-acceleration residual, but is defined **mechanistically** via a reference profile in *feature space* (rather than only via a regression of \(\hat{y}\) on \(\tau\)).

---

## What ClockSHAP returns

For each sample, ClockSHAP computes:

- **Predicted age** (`predicted`): the clock’s predicted biological age, \(\hat{y}\).
- **Expected clock age** (`expected`): the clock’s prediction for a *reference* individual at the same chronological age, \(\hat{y}^{\*}(\tau)\).
- **Age-matched deviation** (`deviation`): \(\hat{y}-\hat{y}^{\*}(\tau)\).
- **Per-feature contributions** (`phi`): a matrix of feature attributions whose row-sums equal the deviation.

---

## Related quantities (and how they differ)

ClockSHAP is designed to explain **age-matched deviation** \(D(\tau)=\hat{y}-\hat{y}^{\*}(\tau)\).  
If your downstream analysis uses another definition, you can still compute it from the same outputs:

```r
# Predicted age from the clock
predicted <- clockshap_predicted(cs)

# 1) Delta age (predicted - chronological)
delta_age <- predicted - age

# 2) Regression-style age acceleration (residual from predicted ~ chronological)
age_accel <- residuals(stats::lm(predicted ~ age))
```

These quantities can be useful, but they answer **different questions** and generally will not equal `clockshap_deviation(cs)` unless your chosen reference makes \(\hat{y}^{\*}(\tau)\approx \tau\) (or matches the specific regression you use).

---
## Installation

ClockSHAP is currently distributed via GitHub.

```r
# install.packages("pak")
pak::pak("<org-or-user>/ClockSHAP")
```

---

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

# Decompose age-matched deviation into feature contributions
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

---

## The reference profile (why “expected” is age-adjusted)

Aging-clock input features often change systematically with chronological age.
If we use a single, global baseline (e.g., the unconditional mean feature vector), then the “baseline” can be inappropriate for older vs. younger individuals.

ClockSHAP addresses this by defining an **age-conditioned reference profile** for each feature:

$$E[X_k \mid \tau] = \gamma_{0k} + \gamma_{1k}\cdot \tau$$

Given a clock and a reference profile, ClockSHAP computes the **expected clock age** by evaluating the clock at the expected (age-conditioned) feature values:

$$\hat{y}^{\*}(\tau) = \hat{y}\big(E[X \mid \tau]\big)$$

The resulting **deviation** and **feature attributions** are therefore interpretable as effects that explain:

> relative biological aging (**older/younger than expected for chronological age**),
> not merely absolute predicted age.

### A note on “age bias” and calibration

Some analyses “correct” age-related bias by fitting a regression of \(\hat{y}\) on \(\tau\) and working with residuals, or by re-scaling predictions to match a diagonal \(\hat{y}\approx\tau\).

ClockSHAP takes a different approach:

- it keeps the clock as-is (no forced calibration step inside the method), and
- it defines “expected” through an **explicit, age-conditioned baseline in feature space**.

If you want your deviations to behave like a particular form of age acceleration, you should **choose and report** the reference cohort/profile used to compute \(\hat{y}^{\*}(\tau)\).

---

## Why we call it “SHAP”

ClockSHAP is called “SHAP” because, for **linear-additive models**, SHAP values have a **closed-form solution**.

Consider a standardized linear clock:

$$\hat{y}(x) = \alpha + \sum_k \beta_k\, z_k$$
where \(z_k = (x_k - \mu_k)/\sigma_k\).

If we choose a baseline feature vector \(z^*\), then the SHAP value for feature \(k\) in a linear model (under the common *interventional / additive* formulation) reduces to:

$$\phi_k = \beta_k\,(z_k - z_k^*)$$

and the attributions add up exactly:

$$\sum_k \phi_k = \hat{y}(z) - \hat{y}(z^*)$$

### The ClockSHAP twist: an age-conditioned baseline

In aging clocks, the most meaningful baseline depends on chronological age.
ClockSHAP therefore uses an **age-conditioned baseline** derived from the reference profile:

$$x_k^*(\tau) = \gamma_{0k} + \gamma_{1k}\cdot \tau$$
$$z_k^*(\tau) = (x_k^*(\tau) - \mu_k)/\sigma_k$$

ClockSHAP then defines:

- **Expected clock age**: \(\hat{y}^{\*}(\tau) = \hat{y}(z^*(\tau))\)
- **Age-matched deviation**: \(D(\tau) = \hat{y}(z) - \hat{y}(z^*(\tau))\)
- **Feature contributions**: \(\phi_k = \beta_k\,(z_k - z_k^*(\tau))\)

This preserves the desirable SHAP property (exact additivity) while aligning the baseline with the biology of aging clocks.

---

## Status

This package is under active development.
The public API may change prior to a stable release.

---

## Citation

Manuscript in preparation.

For manuscript submission, we recommend citing a **tagged GitHub release** and (optionally) an archived DOI (e.g., Zenodo) corresponding to the exact version used in the paper.
