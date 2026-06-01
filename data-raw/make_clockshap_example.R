## Generates the bundled synthetic example dataset `clockshap_example`.
##
## Run with the package loaded:
##   devtools::load_all(".")
##   source("data-raw/make_clockshap_example.R")
##
## Produces data/clockshap_example.rda.
## The data are entirely synthetic; gene symbols are used only for legibility.

set.seed(42)

genes <- c("CDKN2A", "GDF15", "EDA2R", "IL6", "LMNB1", "ELN", "COL1A1", "FOXO3")
n     <- 120
age   <- round(runif(n, 30, 85), 1)

age_slope <- c(
  CDKN2A = 0.040, GDF15 = 0.038, EDA2R = 0.042, IL6 = 0.036,
  LMNB1 = -0.039, ELN = -0.037, COL1A1 = 0.035, FOXO3 = -0.034
)
baseline <- c(
  CDKN2A = 3, GDF15 = 4, EDA2R = 2, IL6 = 2.5,
  LMNB1 = 8, ELN = 7, COL1A1 = 5, FOXO3 = 6
)
noise_sd <- c(
  CDKN2A = 0.80, GDF15 = 0.85, EDA2R = 0.80, IL6 = 0.90,
  LMNB1 = 0.80, ELN = 0.85, COL1A1 = 0.85, FOXO3 = 0.80
)

X <- sapply(
  genes,
  function(g) baseline[g] + age_slope[g] * age + rnorm(n, 0, noise_sd[g])
)
rownames(X) <- sprintf("S%03d", seq_len(n))
colnames(X) <- genes
names(age) <- rownames(X)

ref_trend <- vapply(genes, function(g) coef(lm(X[, g] ~ age)), numeric(2))
age_expected <- function(a) ref_trend[1, ] + ref_trend[2, ] * a

mu     <- colMeans(X)
sigma  <- apply(X, 2, sd)
Z      <- scale(X, center = mu, scale = sigma)
lambda <- 40
beta   <- as.numeric(solve(
  crossprod(Z) + lambda * diag(length(genes)),
  crossprod(Z, age - mean(age))
))
names(beta) <- genes
alpha <- mean(age)

clock <- linear_clock(alpha = alpha, beta = beta, mu = mu, sigma = sigma)

shift_genes <- c("CDKN2A", "GDF15", "IL6")
shift_z     <- c(0.7, 0.7, 0.6)

X["S001", ] <- age_expected(age["S001"])
X["S001", shift_genes] <- X["S001", shift_genes] + shift_z * sigma[shift_genes]
X["S002", ] <- age_expected(age["S002"])
X["S002", shift_genes] <- X["S002", shift_genes] - shift_z * sigma[shift_genes]

clockshap_example <- list(features = X, age = age, clock = clock)

if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(clockshap_example, overwrite = TRUE)
} else {
  dir.create("data", showWarnings = FALSE)
  save(clockshap_example,
       file = file.path("data", "clockshap_example.rda"),
       compress = "xz",
       version = 2)
}

ref <- fit_reference_profile(features = as.data.frame(X), age = age)
cs  <- clockshap(features = X, age = age, clock = clock, reference = ref)

cat(sprintf(
  "predicted %.1f-%.1f | cor(pred,age)=%.2f | slope=%.2f | dev sd=%.1f\n",
  min(cs$predicted), max(cs$predicted), cor(cs$predicted, age),
  coef(lm(cs$predicted ~ age))[2], sd(cs$deviation)
))
cat(sprintf(
  "S001 deviation = %+.1f | S002 deviation = %+.1f\n",
  cs$deviation[1], cs$deviation[2]
))
cat(sprintf(
  "additivity max|rowSums(phi) - deviation| = %.1e\n",
  max(abs(rowSums(cs$phi) - cs$deviation))
))
