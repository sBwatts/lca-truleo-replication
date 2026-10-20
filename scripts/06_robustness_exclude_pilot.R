# 06_robustness_exclude_pilot.R
# purpose: robustness check—apply pre-fit Wave 1 and Wave 2 models to main study vs pilot subsamples
#   Tests whether class structure and prevalences differ between pilot (early deployment) and main study
#   Uses best-fit models (k=4) without re-estimation
#
# outputs:
#   - Class prevalence comparison (main study vs pilot) for Wave 1 and Wave 2
#   - Stacked bar charts (separate, non-shared legends to emphasize independent Wave structures)
#   - Chi-square tests comparing class distributions between pilot and main study

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

# Listwise deletion
pool_complete <- pool.clean |>
  filter(complete.cases(across(all_of(truleo_items))))

cat("ROBUSTNESS CHECK: EXCLUDE PILOT STUDY\n")
cat("Applying pre-fit models from wave-specific LCA\n")
cat("Total sample N:", nrow(pool_complete), "\n\n")

# Split by pilot status
main_study <- pool_complete |> filter(pilot == 0)
pilot_study <- pool_complete |> filter(pilot == 1)

cat("Main study N:", nrow(main_study), "\n")
cat("Pilot study N:", nrow(pilot_study), "\n\n")

# CLASSIFY SUBSAMPLES USING PRE-FIT MODELS ----


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

# Wave 1: Apply best_w1 to both main and pilot
cat("Classifying Wave 1 respondents using pre-fit model...\n")

main_w1 <- main_study |> filter(post == 0)
pilot_w1 <- pilot_study |> filter(post == 0)

main_w1_post <- compute_posterior_log(best_w1, as.matrix(main_w1[, truleo_items]))
pilot_w1_post <- compute_posterior_log(best_w1, as.matrix(pilot_w1[, truleo_items]))

main_w1_classified <- main_w1 |>
  mutate(class = max.col(main_w1_post), sample_type = "Main Study")

pilot_w1_classified <- pilot_w1 |>
  mutate(class = max.col(pilot_w1_post), sample_type = "Pilot")

# Wave 2: Apply best_w2 to both main and pilot
cat("Classifying Wave 2 respondents using pre-fit model...\n")

main_w2 <- main_study |> filter(post == 1)
pilot_w2 <- pilot_study |> filter(post == 1)

main_w2_post <- compute_posterior_log(best_w2, as.matrix(main_w2[, truleo_items]))
pilot_w2_post <- compute_posterior_log(best_w2, as.matrix(pilot_w2[, truleo_items]))

main_w2_classified <- main_w2 |>
  mutate(class = max.col(main_w2_post), sample_type = "Main Study")

pilot_w2_classified <- pilot_w2 |>
  mutate(class = max.col(pilot_w2_post), sample_type = "Pilot")

# COMPARE CLASS PREVALENCES ----


cat("\n--- WAVE 1: CLASS PREVALENCES ---\n")
cat("Main Study:\n")
print(prop.table(table(main_w1_classified$class)))
cat("\nPilot:\n")
print(prop.table(table(pilot_w1_classified$class)))

# Chi-square test for Wave 1
tbl_w1 <- table(main_w1_classified$class, main_w1_classified$sample_type)
# Recalculate with pilot
tbl_w1 <- table(c(main_w1_classified$class, pilot_w1_classified$class),
                c(rep("Main Study", nrow(main_w1_classified)), rep("Pilot", nrow(pilot_w1_classified))))
chi_w1 <- chisq.test(tbl_w1)
cat("\nWave 1 Chi-square test (Main Study vs Pilot):\n")
cat("χ²(3) =", round(chi_w1$statistic, 2), ", p =", round(chi_w1$p.value, 3), "\n")

cat("\n--- WAVE 2: CLASS PREVALENCES ---\n")
cat("Main Study:\n")
print(prop.table(table(main_w2_classified$class)))
cat("\nPilot:\n")
print(prop.table(table(pilot_w2_classified$class)))

# Chi-square test for Wave 2
tbl_w2 <- table(c(main_w2_classified$class, pilot_w2_classified$class),
                c(rep("Main Study", nrow(main_w2_classified)), rep("Pilot", nrow(pilot_w2_classified))))
chi_w2 <- chisq.test(tbl_w2)
cat("\nWave 2 Chi-square test (Main Study vs Pilot):\n")
cat("χ²(3) =", round(chi_w2$statistic, 2), ", p =", round(chi_w2$p.value, 3), "\n")

# STACKED BAR CHARTS ----


library(patchwork)

# Wave 1 stacked bar
fig_w1 <- bind_rows(main_w1_classified, pilot_w1_classified) |>
  mutate(
    sample_type = factor(sample_type, levels = c("Main Study", "Pilot")),
    class = factor(class, levels = 1:4, labels = paste0("Class ", 1:4))
  ) |>
  count(sample_type, class) |>
  group_by(sample_type) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |>
  ggplot(aes(x = sample_type, y = prop, fill = class)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = c("Class 1" = "#777777", "Class 2" =  "#C73E1D",
               "Class 3" = "#2E86AB", "Class 4" = "#F18F01"),
    name = NULL
  ) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion") +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

# Wave 2 stacked bar
fig_w2 <- bind_rows(main_w2_classified, pilot_w2_classified) |>
  mutate(
    sample_type = factor(sample_type, levels = c("Main Study", "Pilot")),
    class = factor(class, levels = 1:4, labels = paste0("Class ", 1:4))
  ) |>
  count(sample_type, class) |>
  group_by(sample_type) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |>
  ggplot(aes(x = sample_type, y = prop, fill = class)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = c("Class 1" = "#F18F01", "Class 2" = "#2E86AB",
               "Class 3" = "#C73E1D", "Class 4" = "#777777"),
    name = NULL
  ) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion") +
  theme_bw() +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

fig_combined <- (fig_w1 | fig_w2) +
  plot_annotation(tag_levels = "A")

fig_combined

ggsave("output/appendix/fig_robustness_exclude_pilot.png", fig_combined,
       dpi = 350, width = 10, height = 8)
cat("\nSaved: output/appendix/fig_robustness_exclude_pilot.png\n")
