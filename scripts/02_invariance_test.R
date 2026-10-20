# purpose: test measurement invariance of four-class LCA solution across Wave 1 and Wave 2
#   using multi-group LCA (glca package)
#
# approach:
#   configural model: item-response probabilities free by wave (allows classes to differ)
#   metric model: item-response probabilities constrained equal across waves (classes must be identical)
#   LR test: compares configural vs metric fit; non-significance indicates invariance
#
# outputs:
#   - Chi-square test of metric vs configural invariance
#   - Interpretation of class replication and non-invariance implications

source("scripts/00_packages.R")
library(glca)
library(poLCA)

set.seed(123)

pool.clean <- import("data/pool.clean.rds", trust = TRUE)

# listwise deletion
pool_complete <- pool.clean |>
  filter(complete.cases(across(all_of(truleo_items))))

pool_glca <- pool_complete |>
  mutate(
    across(all_of(truleo_items), ~ factor(., levels = 1:5)),
    wave = factor(post, levels = c(0, 1), labels = c("Wave 1", "Wave 2"))
  )

cat("N Wave 1:", sum(pool_glca$wave == "Wave 1"), "\n")
cat("N Wave 2:", sum(pool_glca$wave == "Wave 2"), "\n\n")

glca_formula <- as.formula(
  paste0("item(", paste(truleo_items, collapse = ", "), ") ~ 1")
)
polca_formula <- as.formula(
  paste0("cbind(", paste(truleo_items, collapse = ", "), ") ~ 1")
)


# step 1: configural model ------

cat("Fitting configural model (measure.inv = FALSE) ...\n")
fit_configural <- glca(glca_formula,
                       group    = wave,
                       data     = pool_glca,
                       nclass   = 4,
                       measure.inv = FALSE,
                       coeff.inv   = FALSE,
                       n.init   = 50,
                       seed     = 123,
                       verbose  = FALSE)

cat("Configural LL:", round(fit_configural$gof$loglik, 2), "\n")
cat("Configural npar:", fit_configural$model$npar, "\n")


# step 2: metric model -----

cat("\nFitting metric model (measure.inv = TRUE) ...\n")
fit_metric <- glca(glca_formula,
                   group    = wave,
                   data     = pool_glca,
                   nclass   = 4,
                   measure.inv = TRUE,
                   coeff.inv   = FALSE,
                   n.init   = 50,
                   seed     = 123,
                   verbose  = FALSE)

cat("Metric LL:", round(fit_metric$gof$loglik, 2), "\n")
cat("Metric npar:", fit_metric$model$npar, "\n")


# step 3: formal LRT -------------------------------------------------------

cat("\n--- Formal Measurement Invariance LRT (gofglca) ---\n")
lrt_result <- gofglca(fit_configural, fit_metric, test = "chisq")
print(lrt_result)

# extract key values for reporting
lrt_stat <- lrt_result$dtable["1", "Deviance"]
lrt_df   <- lrt_result$dtable["1", "Df"]
lrt_p    <- lrt_result$dtable["1", "Pr(>Chi)"]

cat(sprintf("\nLRT χ²(%d) = %.2f, p = %.4f\n", lrt_df, lrt_stat, lrt_p))


# step 4: save -------------------------------------------------------------

save(fit_configural, fit_metric, lrt_result,
     lrt_stat, lrt_df, lrt_p,
     file = "data/lca_invariance.RData")

cat("\nSaved: data/lca_invariance.RData\n")

