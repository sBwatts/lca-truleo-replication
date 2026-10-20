# purpose: robustness check—fit independent LCA models excluding comprehension items
#   Excludes: truleo2_1 ("Understand how it works"), truleo2_2 ("Understand dept. plan for use")
#   Tests whether class structure is robust to removal of understanding-focused items
#
# outputs:
#   - Combined LCA enumeration table (k=1:5, both waves) with LL, AIC, BIC, aBIC, entropy, min class
#   - Item-response probability figures (stacked horizontal bars) for Wave 1 and Wave 2 classes
#   - Modal response tables (most probable response by item and class)
#   - Wave 1 and Wave 2 class prevalence comparison
#   - Docx table: combined fit indices

source("scripts/00_packages.R")

library(poLCA)

set.seed(123)

pool.clean <- import("data/pool.clean.rds", trust = TRUE)

# Define attitude items (EXCLUDE comprehension items)
truleo_items <- c(
  paste0("truleo1_", 1:6),
  paste0("truleo2_", 3:7)  # Skip truleo2_1 and truleo2_2
)

item_labels <- c(
  truleo1_1 = "Fish for policy violations",
  truleo1_2 = "Improve my public interactions",
  truleo1_3 = "Aid in justifying actions",
  truleo1_4 = "Improve colleagues' interactions",
  truleo1_5 = "Make job more dangerous",
  truleo1_6 = "Valuable addition overall",
  truleo2_3 = "Excited to see performance captured",
  truleo2_4 = "Concerned about dept. use",
  truleo2_5 = "Improve perf. eval. process",
  truleo2_6 = "Publicize good work",
  truleo2_7 = "Will be used to punish officers"
)

# Listwise deletion
pool_complete <- pool.clean |>
  filter(complete.cases(across(all_of(truleo_items))))

cat("ROBUSTNESS CHECK 1: EXCLUDE COMPREHENSION ITEMS\n")
cat("Items included:", length(truleo_items), "\n")
cat("Analytic sample:", nrow(pool_complete), "\n\n")

# Split by wave
wave1_data <- pool_complete |> filter(post == 0)
wave2_data <- pool_complete |> filter(post == 1)

cat("Wave 1 N:", nrow(wave1_data), "\n")
cat("Wave 2 N:", nrow(wave2_data), "\n\n")

# LCA FORMULA ----


lca_formula <- as.formula(
  paste0("cbind(", paste(truleo_items, collapse = ", "), ") ~ 1")
)

# FIT MODELS: k = 1:5 FOR EACH WAVE ----


cat("=== FITTING WAVE 1 MODELS ===\n")
lca_wave1 <- list()

for (k in 1:5) {
  cat("Fitting Wave 1 k =", k, "...\n")
  lca_wave1[[k]] <- poLCA(
    lca_formula,
    data = wave1_data,
    nclass = k,
    maxiter = 5000,
    nrep = ifelse(k == 1, 1, 50),
    verbose = FALSE
  )
}

cat("\n=== FITTING WAVE 2 MODELS ===\n")
lca_wave2 <- list()

for (k in 1:5) {
  cat("Fitting Wave 2 k =", k, "...\n")
  lca_wave2[[k]] <- poLCA(
    lca_formula,
    data = wave2_data,
    nclass = k,
    maxiter = 5000,
    nrep = ifelse(k == 1, 1, 50),
    verbose = FALSE
  )
}

# CALCULATE ENTROPY ----


calc_entropy <- function(model) {
  p <- model$posterior
  n <- nrow(p)
  k <- ncol(p)
  1 - (-sum(p * log(p + 1e-10)) / (n * log(k)))
}

# ENUMERATION TABLES ----


n_w1 <- nrow(wave1_data)
fit_wave1 <- tibble(
  k = 1:5,
  log_lik = sapply(lca_wave1[1:5], `[[`, "llik"),
  n_params = sapply(lca_wave1[1:5], `[[`, "npar"),
  AIC = sapply(lca_wave1[1:5], `[[`, "aic"),
  BIC = sapply(lca_wave1[1:5], `[[`, "bic"),
  aBIC = -2 * log_lik + log((n_w1 + 2) / 24) * n_params,
  entropy = c(NA, sapply(lca_wave1[2:5], calc_entropy)),
  min_class = c(NA, sapply(lca_wave1[2:5],
                           function(m) min(table(m$predclass))))
)

cat("\n--- WAVE 1 ENUMERATION TABLE ---\n")
print(fit_wave1)

n_w2 <- nrow(wave2_data)
fit_wave2 <- tibble(
  k = 1:5,
  log_lik = sapply(lca_wave2[1:5], `[[`, "llik"),
  n_params = sapply(lca_wave2[1:5], `[[`, "npar"),
  AIC = sapply(lca_wave2[1:5], `[[`, "aic"),
  BIC = sapply(lca_wave2[1:5], `[[`, "bic"),
  aBIC = -2 * log_lik + log((n_w2 + 2) / 24) * n_params,
  entropy = c(NA, sapply(lca_wave2[2:5], calc_entropy)),
  min_class = c(NA, sapply(lca_wave2[2:5],
                           function(m) min(table(m$predclass))))
)

cat("\n--- WAVE 2 ENUMERATION TABLE ---\n")
print(fit_wave2)

# COMBINED FIT INDICES TABLE ----


combined_fit <- bind_rows(
  fit_wave1 |> mutate(wave = "Wave 1"),
  fit_wave2 |> mutate(wave = "Wave 2")
) |>
  dplyr::select(wave, k, log_lik, n_params, AIC, BIC, aBIC, entropy, min_class) |>
  mutate(
    log_lik = round(log_lik, 2),
    AIC = round(AIC, 2),
    BIC = round(BIC, 2),
    aBIC = round(aBIC, 2),
    entropy = round(entropy, 3),
    min_class = as.integer(min_class)
  )

cat("\n--- COMBINED LCA ENUMERATION TABLE (EXCLUDE COMPREHENSION) ---\n")
print(combined_fit)

# Export to docx
combined_fit_table <- combined_fit |>
  rename(
    "Wave" = wave,
    "k" = k,
    "Log-Likelihood" = log_lik,
    "Parameters" = n_params,
    "AIC" = AIC,
    "BIC" = BIC,
    "aBIC" = aBIC,
    "Entropy" = entropy,
    "Min Class Size" = min_class
  ) |>
  flextable::flextable() |>
  flextable::autofit()

flextable::save_as_docx(combined_fit_table, path = "output/appendix/table_lca_enumeration_exclude_comprehension.docx")
cat("Saved: output/appendix/table_lca_enumeration_exclude_comprehension.docx\n")

# SELECT k = 4 FOR BOTH WAVES ----


k_best <- 4
best_w1 <- lca_wave1[[k_best]]
best_w2 <- lca_wave2[[k_best]]

cat("\n--- WAVE 1 (k=4) ---\n")
cat("Entropy:", round(calc_entropy(best_w1), 3), "\n")
cat("Class prevalences:\n")
print(prop.table(table(best_w1$predclass)))

cat("\n--- WAVE 2 (k=4) ---\n")
cat("Entropy:", round(calc_entropy(best_w2), 3), "\n")
cat("Class prevalences:\n")
print(prop.table(table(best_w2$predclass)))

# EXTRACT PROBABILITIES ----


extract_probs <- function(model) {
  imap_dfr(model$probs, function(mat, item) {
    as.data.frame(mat) |>
      mutate(item = item, class = paste0("Class ", seq_len(nrow(mat)))) |>
      pivot_longer(-c(item, class), names_to = "response", values_to = "prob")
  })
}

probs_w1 <- extract_probs(best_w1)
probs_w2 <- extract_probs(best_w2)

# ITEM-RESPONSE PROBABILITY FIGURES ----


library(patchwork)

fig_w1 <- probs_w1 |>
  mutate(
    response_num = as.integer(str_extract(response, "\\d+")),
    class_label = factor(paste0("Class ", str_extract(class, "\\d+")),
                         levels = paste0("Class ", 1:4)),
    item_label = factor(item_labels[item],
                        levels = rev(unname(item_labels)))
  ) |>
  ggplot(aes(x = item_label, y = prob,
             fill = factor(response_num, levels = 5:1))) +
  geom_col(width = 0.8) +
  coord_flip() +
  facet_wrap(~class_label, nrow = 1) +
  scale_fill_manual(
    values = c("5" = "#1a5276", "4" = "#2e86c1",
               "3" = "#aab7b8", "2" = "#e59866", "1" = "#922b21"),
    labels = c("5" = "Str. Agree", "4" = "Agree",
               "3" = "Neutral", "2" = "Disagree", "1" = "Str. Disagree"),
    name = NULL,
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Probability", title = "Wave 1: Item-Response Probabilities (Excl. Comprehension Items)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 7.5),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = unit(1.2, "lines")
  )

fig_w2 <- probs_w2 |>
  mutate(
    response_num = as.integer(str_extract(response, "\\d+")),
    class_label = factor(paste0("Class ", str_extract(class, "\\d+")),
                         levels = paste0("Class ", 1:4)),
    item_label = factor(item_labels[item],
                        levels = rev(unname(item_labels)))
  ) |>
  ggplot(aes(x = item_label, y = prob,
             fill = factor(response_num, levels = 5:1))) +
  geom_col(width = 0.8) +
  coord_flip() +
  facet_wrap(~class_label, nrow = 1) +
  scale_fill_manual(
    values = c("5" = "#1a5276", "4" = "#2e86c1",
               "3" = "#aab7b8", "2" = "#e59866", "1" = "#922b21"),
    labels = c("5" = "Str. Agree", "4" = "Agree",
               "3" = "Neutral", "2" = "Disagree", "1" = "Str. Disagree"),
    name = NULL,
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Probability", title = "Wave 2: Item-Response Probabilities (Excl. Comprehension Items)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 7.5),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = unit(1.2, "lines")
  )

fig_combined <- (fig_w1 / fig_w2) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_combined

ggsave("output/appendix/fig_robustness_exclude_comprehension.png", fig_combined,
       dpi = 350, width = 10, height = 8)
cat("\nSaved: output/appendix/fig_robustness_exclude_comprehension.png\n")

# SAVE RESULTS ----


save(
  best_w1, best_w2,
  fit_wave1, fit_wave2,
  probs_w1, probs_w2,
  file = "data/robustness_exclude_comprehension.RData"
)

cat("\nSaved: data/robustness_exclude_comprehension.RData\n")
