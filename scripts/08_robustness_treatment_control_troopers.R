# 08_robustness_treatment_control_troopers.R
# purpose: robustness check—apply pre-fit Wave 2 model to treatment vs control arms (troopers only subsample)
#   Tests whether RCT assignment to treatment (Truleo access) affects class distribution for front-line officers
#   Restricts to rank1 == "Officer/trooper" to isolate front-line perspective
#   Uses best-fit Wave 2 model (k=4) without re-estimation
#
# outputs:
#   - Class prevalence comparison (treatment vs control) for Wave 2 troopers only
#   - Stacked bar chart with single legend
#   - Chi-square test comparing class distributions between arms (trooper subsample)

source("scripts/00_packages.R")

library(poLCA)

set.seed(123)

# Load pre-fit models from wave-specific LCA
load("data/wave_specific_lca_results.RData")

pool.clean <- import("data/pool.clean.rds", trust = TRUE)

truleo_items <- c(
  paste0("truleo1_", 1:6),
  paste0("truleo2_", 1:7)
)

# Listwise deletion, filter to Wave 2, and troopers only
pool_complete <- pool.clean |>
  filter(complete.cases(across(all_of(truleo_items))),
         post == 1,
         rank1 == "Officer/trooper")

cat("ROBUSTNESS CHECK: TREATMENT VS CONTROL (TROOPERS ONLY, WAVE 2)\n")
cat("Applying pre-fit Wave 2 model to treatment and control arms\n")
cat("Total Wave 2 Troopers N:", nrow(pool_complete), "\n\n")

# Split by treatment arm
treatment <- pool_complete |> filter(truleo_access== 1)
control <- pool_complete |> filter(truleo_access== 0)

cat("Treatment N:", nrow(treatment), "\n")
cat("Control N:", nrow(control), "\n\n")

# CLASSIFY ARMS USING PRE-FIT WAVE 2 MODEL ----


# Function to compute posterior probabilities using pre-fit model
compute_posterior_log <- function(model, y) {
  y      <- matrix(as.integer(y), nrow = nrow(y), ncol = ncol(y))
  k      <- as.integer(length(model$P))
  n      <- as.integer(nrow(y))
  log_pr <- log(model$P)
  log_ll <- matrix(0.0, nrow = n, ncol = k)
  for (j in seq_along(model$probs)) {
    item_mat     <- model$probs[[j]]
    for (ki in seq_len(k)) {
      log_ll[, ki] <- log_ll[, ki] +
        log(item_mat[ki, y[, j]] + 1e-300)
    }
  }
  log_post  <- sweep(log_ll, 2, log_pr, "+")
  max_lp    <- apply(log_post, 1, max)
  log_denom <- max_lp + log(rowSums(exp(log_post - max_lp)))
  exp(log_post - log_denom)
}

cat("Classifying treatment arm...\n")
treatment_post <- compute_posterior_log(best_w2, as.matrix(treatment[, truleo_items]))
treatment_classified <- treatment |>
  mutate(class = max.col(treatment_post), arm = "Treatment")

cat("Classifying control arm...\n")
control_post <- compute_posterior_log(best_w2, as.matrix(control[, truleo_items]))
control_classified <- control |>
  mutate(class = max.col(control_post), arm = "Control")

# COMPARE CLASS PREVALENCES ----


cat("\n--- CLASS PREVALENCES (TROOPERS ONLY) ---\n")
cat("Treatment:\n")
print(prop.table(table(treatment_classified$class)))
cat("\nControl:\n")
print(prop.table(table(control_classified$class)))

# Chi-square test for treatment vs control (troopers only)
tbl_arm <- table(c(treatment_classified$class, control_classified$class),
                 c(rep("Treatment", nrow(treatment_classified)), rep("Control", nrow(control_classified))))
chi_arm <- chisq.test(tbl_arm)
cat("\nChi-square test (Treatment vs Control, Troopers Only):\n")
cat("χ²(3) =", round(chi_arm$statistic, 2), ", p =", round(chi_arm$p.value, 3), "\n")

# STACKED BAR CHART ----


fig <- bind_rows(treatment_classified, control_classified) |>
  mutate(
    arm = factor(arm, levels = c("Treatment", "Control")),
    class = factor(class, levels = 1:4, labels = paste0("Class ", 1:4))
  ) |>
  count(arm, class) |>
  group_by(arm) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |>
  ggplot(aes(x = arm, y = prop, fill = class)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = c("Class 1" = "#F18F01", "Class 2" = "#2E86AB",
               "Class 3" = "#C73E1D", "Class 4" = "#777777"),
    name = NULL
  ) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, 
       y = "Proportion") +
       # title = "Wave 2 Troopers: Treatment vs Control") +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

fig

ggsave("output/appendix/fig_robustness_treatment_control_troopers.png", fig,
       dpi = 350, width = 8, height = 6)
cat("\nSaved: output/appendix/fig_robustness_treatment_control_troopers.png\n")
