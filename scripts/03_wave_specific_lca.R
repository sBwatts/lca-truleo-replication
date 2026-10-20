# purpose: fit independent LCA models to Wave 1 and Wave 2 separately
#   Wave 1 and Wave 2 classes estimated independently (measurement non-invariance)
#
# outputs:
#   - Combined LCA enumeration table (k=1:5, both waves) with LL, AIC, BIC, aBIC, entropy, min class
#   - Item-response probability figures (stacked horizontal bars) for Wave 1 and Wave 2 classes
#   - Modal response tables (most probable response by item and class)
#   - Demographic stratification tables by class (age, tenure, org justice, rank, gender, race, education, military service)
#   - Wave 1 and Wave 2 class prevalence comparison
#   - Item mean response convergence plot (Wave 1 vs Wave 2, 1-5 scale, with overlay t-test table)
#   - Item-level analysis: paired t-tests and Cohen's d effect sizes (sorted by effect size)
#   - RData file with all fit objects, classified data, and results

source("scripts/00_packages.R")

library(poLCA)
library(effsize)

set.seed(123)

pool.clean <- import("data/pool.clean.rds", trust = TRUE)

truleo_items <- c(
  paste0("truleo1_", 1:6),
  paste0("truleo2_", 1:7)
)

item_labels <- c(
  truleo1_1 = "Fish for policy violations",
  truleo1_2 = "Improve my public interactions",
  truleo1_3 = "Aid in justifying actions",
  truleo1_4 = "Improve colleagues' interactions",
  truleo1_5 = "Make job more dangerous",
  truleo1_6 = "Valuable addition overall",
  truleo2_1 = "Understand how it works",
  truleo2_2 = "Understand dept. plan for use",
  truleo2_3 = "Excited to see performance captured",
  truleo2_4 = "Concerned about dept. use",
  truleo2_5 = "Improve perf. eval. process",
  truleo2_6 = "Publicize good work",
  truleo2_7 = "Will be used to punish officers"
)

# Listwise deletion
pool_complete <- pool.clean |>
  filter(complete.cases(across(all_of(truleo_items))))

cat("Analytic sample after listwise deletion:", nrow(pool_complete), "\n\n")

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


# Wave 1
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

cat("\n--- WAVE 1 LCA ENUMERATION TABLE (k=1:5) ---\n")
print(fit_wave1)

# Wave 2
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

cat("\n--- WAVE 2 LCA ENUMERATION TABLE (k=1:5) ---\n")
print(fit_wave2)

# COMBINED FIT INDICES TABLE ----


# Create combined table for export
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

cat("\n--- COMBINED LCA ENUMERATION TABLE ---\n")
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

combined_fit_table 

flextable::save_as_docx(combined_fit_table, path = "output/main/table_lca_enumeration_combined.docx")
cat("Saved: output/main/table_lca_enumeration_combined.docx\n")

# SELECT k = 4 FOR BOTH WAVES ----


k_best <- 4
best_w1 <- lca_wave1[[k_best]]
best_w2 <- lca_wave2[[k_best]]

cat("\n--- SELECTED k = 4 FOR BOTH WAVES ---\n\n")

cat("WAVE 1 (k=4):\n")
cat("Log-Likelihood:", round(best_w1$llik, 2), "\n")
cat("AIC:", round(best_w1$aic, 2), "\n")
cat("BIC:", round(best_w1$bic, 2), "\n")
cat("Entropy:", round(calc_entropy(best_w1), 3), "\n")
cat("Minimum class size:", min(table(best_w1$predclass)), "\n")
cat("Class prevalences:\n")
print(prop.table(table(best_w1$predclass)))

cat("\nWAVE 2 (k=4):\n")
cat("Log-Likelihood:", round(best_w2$llik, 2), "\n")
cat("AIC:", round(best_w2$aic, 2), "\n")
cat("BIC:", round(best_w2$bic, 2), "\n")
cat("Entropy:", round(calc_entropy(best_w2), 3), "\n")
cat("Minimum class size:", min(table(best_w2$predclass)), "\n")
cat("Class prevalences:\n")
print(prop.table(table(best_w2$predclass)))

# EXTRACT PROBABILITIES FOR VISUALIZATION ----


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


# Wave 1 figure
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
  labs(x = NULL, 
       y = "Probability") +
       # title = "Wave 1: Item-Response Probabilities by Class") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 7.5),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = unit(1.2, "lines")
  )

# Wave 2 figure
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
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, 
       y = "Probability") +
       # title = "Wave 2: Item-Response Probabilities by Class") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 7.5),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.spacing.x = unit(1.2, "lines")
  )

# Combine into patchwork figure (A/B panels, one column, two rows)
library(patchwork)

fig_combined <- (fig_w1 / fig_w2) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_combined

ggsave("output/main/fig_wave_comparison_item_response_probs.png", fig_combined,
       dpi = 350, width = 10, height = 8)
cat("\nSaved: output/main/fig_wave_comparison_item_response_probs.png\n")

# MODAL RESPONSE TABLES ----


# Wave 1
modal_w1 <- probs_w1 |>
  group_by(item, class) |>
  slice_max(prob, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    response_num = as.integer(str_extract(response, "\\d+")),
    cell = sprintf("%d (%.0f%%)", response_num, prob * 100),
    class_num = as.numeric(str_extract(class, "\\d+"))
  ) |>
  dplyr::select(item, class_num, cell) |>
  pivot_wider(names_from = class_num, values_from = cell, names_prefix = "Class ") |>
  left_join(tibble(item = names(item_labels), Item = unname(item_labels)), by = "item") |>
  dplyr::select(Item, `Class 1`, `Class 2`, `Class 3`, `Class 4`)

cat("\n--- WAVE 1 MODAL RESPONSES ---\n")
print(modal_w1)

modal_w1 |>
  flextable::flextable() |>
  flextable::autofit() |>
  flextable::save_as_docx(path = "output/appendix/table_modal_response_wave1.docx")
cat("Saved: output/appendix/table_modal_response_wave1.docx\n")

# Wave 2
modal_w2 <- probs_w2 |>
  group_by(item, class) |>
  slice_max(prob, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    response_num = as.integer(str_extract(response, "\\d+")),
    cell = sprintf("%d (%.0f%%)", response_num, prob * 100),
    class_num = as.numeric(str_extract(class, "\\d+"))
  ) |>
  dplyr::select(item, class_num, cell) |>
  pivot_wider(names_from = class_num, values_from = cell, names_prefix = "Class ") |>
  left_join(tibble(item = names(item_labels), Item = unname(item_labels)), by = "item") |>
  dplyr::select(Item, `Class 1`, `Class 2`, `Class 3`, `Class 4`)

cat("\n--- WAVE 2 MODAL RESPONSES ---\n")
print(modal_w2)

modal_w2 |>
  flextable::flextable() |>
  flextable::autofit() |>
  flextable::save_as_docx(path = "output/appendix/table_modal_response_wave2.docx")
cat("Saved: output/appendix/table_modal_response_wave2.docx\n")

# DEMOGRAPHIC TABLES ----


# Classify respondents based on modal class assignment
w1_classified <- wave1_data |>
  mutate(
    lca_class_num = best_w1$predclass,
    lca_class = paste0("Class ", lca_class_num)
  ) |>
  bind_cols(
    as.data.frame(best_w1$posterior) |>
      setNames(paste0("prob_class_", 1:4))
  )

w2_classified <- wave2_data |>
  mutate(
    lca_class_num = best_w2$predclass,
    lca_class = paste0("Class ", lca_class_num)
  ) |>
  bind_cols(
    as.data.frame(best_w2$posterior) |>
      setNames(paste0("prob_class_", 1:4))
  )

# Combine for wave comparison
combined <- bind_rows(
  w1_classified |> mutate(wave = "Wave 1"),
  w2_classified |> mutate(wave = "Wave 2")
)

# DESCRIPTIVE TABLE BY WAVE ----


cat("\n--- DESCRIPTIVE TABLE BY WAVE ---\n")

desc_table_wave <- combined |>
  mutate(
    wave = factor(wave, levels = c("Wave 1", "Wave 2"))
  ) |>
  dplyr::select(wave, leo_years, age, oj, sup, rank1, school, female, race1, military) |>
  tbl_summary(
    by = wave,
    label = list(
      leo_years ~ "Years as LEO",
      age ~ "Age (years)",
      oj ~ "Perceived Organizational Justice (1-5 scale)",
      sup ~ "Perceived Supervisor Support (1-5 scale)",
      rank1 ~ "Rank",
      school ~ "Education",
      female ~ "Gender",
      race1 ~ "Race",
      military ~ "Military service"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd}) [Min: {min}, Max: {max}]",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) |>
  add_p() |>
  as_flex_table() |>
  flextable::set_table_properties(width = 1, layout = "autofit")

desc_table_wave <- set_table_properties(desc_table_wave, opts_word = list(repeat_headers = FALSE))

desc_table_wave

desc_table_wave |>
  flextable::save_as_docx(path = "output/main/desc_table_wave_comparison.docx")
cat("Saved: output/main/desc_table_wave_comparison.docx\n")

# DESCRIPTIVE TABLES BY CLASS (WAVE 1 AND WAVE 2) ----


vars <- c("leo_years", "age", "zoj", "zsup", "rank1", "female", "race1", "school", "military")

make_tbl <- function(dat, wave_label) {
  dat |>
    dplyr::select(lca_class, all_of(vars)) |>
    rename(Class = lca_class) |>
    tbl_summary(
      by = Class,
      label = list(
        age ~ "Age (years)",
        leo_years ~ "Years as LEO",
        rank1 ~ "Rank",
        female ~ "Gender",
        race1 ~ "Race",
        school ~ "Education",
        zoj ~ "Perceived Org. Justice (Z-score)",
        zsup ~ "Perceived Supervisor Support (Z-score)",
        military ~ "Military service"
      ),
      statistic = list(
        all_continuous() ~ "{mean} ({sd})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      missing = "no"
    ) |>
    add_ci(
      method = list(all_continuous() ~ "t.test",
                    all_categorical() ~ "wilson"),
      pattern = "{stat} [{ci}]"
    )
}

# Wave 1 descriptive table by class
cat("\n--- WAVE 1 DESCRIPTIVE TABLE BY CLASS ---\n")

desc_table_w1_class <- make_tbl(w1_classified, "Wave 1") |>
  as_flex_table() |>
  flextable::add_footer_lines("Note: Cell ns may not sum to class totals due to item-level missingness on demographic variables.") |>
  flextable::set_table_properties(width = 1, layout = "autofit")

desc_table_w1_class <- set_table_properties(desc_table_w1_class, opts_word = list(repeat_headers = FALSE))

desc_table_w1_class

desc_table_w1_class |>
  flextable::save_as_docx(path = "output/appendix/desc_table_class_wave1.docx")
cat("Saved: output/appendix/desc_table_class_wave1.docx\n")

# Wave 2 descriptive table by class
cat("\n--- WAVE 2 DESCRIPTIVE TABLE BY CLASS ---\n")

desc_table_w2_class <- make_tbl(w2_classified, "Wave 2") |>
  as_flex_table() |>
  flextable::add_footer_lines("Note: Cell ns may not sum to class totals due to item-level missingness on demographic variables.") |>
  flextable::set_table_properties(width = 1, layout = "autofit")

desc_table_w2_class <- set_table_properties(desc_table_w2_class, opts_word = list(repeat_headers = FALSE))

desc_table_w2_class

desc_table_w2_class |>
  flextable::save_as_docx(path = "output/appendix/desc_table_class_wave2.docx")
cat("Saved: output/appendix/desc_table_class_wave2.docx\n")

# CLASS PREVALENCES BY WAVE ----


cat("\n--- CLASS PREVALENCES ---\n")

w1_prev <- prop.table(table(w1_classified$lca_class_num))
w2_prev <- prop.table(table(w2_classified$lca_class_num))

prev_comparison <- tibble(
  class = paste0("Class ", 1:4),
  wave1_n = as.numeric(table(w1_classified$lca_class_num)),
  wave1_prop = round(w1_prev, 3),
  wave2_n = as.numeric(table(w2_classified$lca_class_num)),
  wave2_prop = round(w2_prev, 3),
)

print(prev_comparison)

# WAVE 1 vs WAVE 2 ITEM-RESPONSE PROBABILITY CONVERGENCE PLOT ----


cat("\n--- CREATING CONVERGENCE PLOT ---\n")

# Calculate mean item responses (1-5 scale, treated as continuous)
# Captures full distribution of each item across all respondents
convergence_data <- tibble(
  item = names(item_labels),
  item_label = unname(item_labels),
  w1_mean = NA_real_,
  w2_mean = NA_real_
)

# Split data by wave and convert to numeric for means
w1_data_means <- pool_complete |>
  filter(post == 0) |>
  dplyr::select(all_of(truleo_items)) |>
  mutate(across(everything(), ~as.numeric(as.character(.))))

w2_data_means <- pool_complete |>
  filter(post == 1) |>
  dplyr::select(all_of(truleo_items)) |>
  mutate(across(everything(), ~as.numeric(as.character(.))))

# Calculate mean response for each item in each wave
# Use w1_classified and w2_classified which retain original response data
for (i in seq_along(truleo_items)) {
  item_col <- truleo_items[i]

  # Wave 1: mean from w1_classified respondents
  w1_vals <- as.numeric(as.character(w1_classified[[item_col]]))
  convergence_data$w1_mean[i] <- mean(w1_vals, na.rm = TRUE)

  # Wave 2: mean from w2_classified respondents
  w2_vals <- as.numeric(as.character(w2_classified[[item_col]]))
  convergence_data$w2_mean[i] <- mean(w2_vals, na.rm = TRUE)
}

# Assign item categories
convergence_data <- convergence_data |>
  mutate(
    category = case_when(
      item %in% c("truleo1_2", "truleo1_3", "truleo1_4", "truleo1_6", "truleo2_3", "truleo2_5", "truleo2_6") ~ "Benefits",
      item %in% c("truleo2_1", "truleo2_2") ~ "Understanding",
      item %in% c("truleo1_1", "truleo1_5", "truleo2_4", "truleo2_7") ~ "Concerns"
    ),
    item_num = 1:13
  )

library(ggrepel)

# Create legend labels with item numbers, ordered by category
convergence_data <- convergence_data |>
  mutate(
    legend_label = paste0(item_num, ". ", item_label),
    # Order by category for legend grouping
    category_order = factor(category,
                           levels = c("Benefits", "Understanding", "Concerns")),
    legend_label = factor(legend_label,
                         levels = legend_label)
  )

# Create color mapping preserving category colors
color_map <- c(
  "Benefits" = "#2E86AB",
  "Understanding" = "#F18F01",
  "Concerns" = "#C73E1D"
)

legend_colors <- setNames(
  color_map[as.character(convergence_data$category)],
  convergence_data$legend_label
)

# Create convergence plot with leader lines from item numbers to dots
fig_convergence <- convergence_data |>
  ggplot(aes(x = w1_mean, y = w2_mean, color = legend_label)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_text_repel(
    aes(label = item_num),
    size = 2.8,
    fontface = "bold",
    segment.size = 0.3,
    segment.color = "black",
    box.padding = 0.3,
    point.padding = 0.1,
    max.overlaps = Inf,
    force = 3,
    force_pull = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = legend_colors,
    name = NULL,
    guide = guide_legend(ncol = 1)
  ) +
  scale_x_continuous(limits = c(0.8, 5.2), breaks = 1:5) +
  scale_y_continuous(limits = c(0.8, 5.2), breaks = 1:5) +
  labs(
    x = "Wave 1: Mean Response (1-5 scale)",
    y = "Wave 2: Mean Response (1-5 scale)",
    # title = "Item Mean Response Convergence Between Waves",
    # subtitle = "Points on diagonal = stable; above diagonal = increased agreement; below diagonal = decreased agreement"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.text = element_text(size = 8)
  )

fig_convergence

ggsave("output/main/fig_convergence_wave1_wave2.png", fig_convergence,
       dpi = 350, width = 9, height = 7)
cat("Saved: output/main/fig_convergence_wave1_wave2.png\n")

# Print convergence summary
cat("\n--- ITEM MEAN RESPONSE SHIFTS (1-5 scale) ---\n")
convergence_summary <- convergence_data |>
  mutate(
    shift = w2_mean - w1_mean
  ) |>
  arrange(desc(abs(shift))) |>
  dplyr::select(item_num, item_label, category, w1_mean, w2_mean, shift) |>
  mutate(
    w1_mean = round(w1_mean, 2),
    w2_mean = round(w2_mean, 2),
    shift = round(shift, 2)
  ) |>
  rename(
    "#" = item_num,
    Item = item_label,
    Category = category,
    "Wave 1 Mean" = w1_mean,
    "Wave 2 Mean" = w2_mean,
    "Shift (W2-W1)" = shift
  )

print(convergence_summary)

# ITEM-LEVEL ANALYSIS: T-TESTS AND EFFECT SIZES ----


cat("\n--- ITEM-LEVEL COMPARISON (WAVE 1 VS WAVE 2) ---\n\n")

# Run t-test and calculate Cohen's d for each item
# Use w1_data_means and w2_data_means from convergence plot section (already properly converted)
item_results <- tibble(
  item = character(),
  label = character(),
  w1_mean = numeric(),
  w1_sd = numeric(),
  w2_mean = numeric(),
  w2_sd = numeric(),
  diff = numeric(),
  t_stat = numeric(),
  p_value = numeric(),
  cohens_d = numeric(),
)

for (i in seq_along(truleo_items)) {
  item <- truleo_items[i]

  # Extract the column and convert to numeric vector
  w1_data <- as.numeric(as.character(pull(w1_data_means, all_of(item))))
  w2_data <- as.numeric(as.character(pull(w2_data_means, all_of(item))))

  # Remove NAs
  w1_data <- w1_data[!is.na(w1_data)]
  w2_data <- w2_data[!is.na(w2_data)]

  # T-test
  t_test <- t.test(w1_data, w2_data)

  # Cohen's d
  cohens_d <- cohen.d(w2_data, w1_data)$estimate

  item_results <- item_results |>
    add_row(
      item = item,
      label = item_labels[item],
      w1_mean = mean(w1_data, na.rm = TRUE),
      w1_sd = sd(w1_data, na.rm = TRUE),
      w2_mean = mean(w2_data, na.rm = TRUE),
      w2_sd = sd(w2_data, na.rm = TRUE),
      diff = mean(w2_data, na.rm = TRUE) - mean(w1_data, na.rm = TRUE),
      t_stat = t_test$statistic,
      p_value = t_test$p.value,
      cohens_d = cohens_d,
    )
}

# Sort by absolute effect size
item_results <- item_results |>
  arrange(desc(abs(cohens_d)))

cat("Item-Level Shifts (sorted by effect size):\n\n")
print(item_results |>
  dplyr::select(label, w1_mean, w2_mean, diff, t_stat, p_value, cohens_d) |>
  mutate(
    w1_mean = round(w1_mean, 2),
    w2_mean = round(w2_mean, 2),
    diff = round(diff, 2),
    t_stat = round(t_stat, 2),
    p_value = round(p_value, 4),
    cohens_d = round(cohens_d, 3),
  ))

# Export to docx
item_table <- item_results |>
  dplyr::select(label, w1_mean, w2_mean, diff, t_stat, p_value, cohens_d) |>
  mutate(
    w1_mean = round(w1_mean, 2),
    w2_mean = round(w2_mean, 2),
    diff = round(diff, 2),
    t_stat = round(t_stat, 2),
    p_value = format(round(p_value, 4), nsmall = 4),
    cohens_d = round(cohens_d, 3),
  ) |>
  rename(
    "Item" = label,
    "Wave 1 M" = w1_mean,
    "Wave 2 M" = w2_mean,
    "Difference" = diff,
    "t-statistic" = t_stat,
    "p-value" = p_value,
    "Cohen's d" = cohens_d,
  )

item_table |>
  flextable::flextable() |>
  flextable::autofit() |>
  flextable::save_as_docx(path = "output/main/table_item_level_shifts.docx")


cat("\nSaved: output/main/table_item_level_shifts.docx\n\n")

# Summary
cat("Largest positive shift (more agreement):\n")
print(item_results |> filter(diff == max(diff)) |> dplyr::select(label, diff, cohens_d))
cat("\nLargest negative shift (less agreement):\n")
print(item_results |> filter(diff == min(diff)) |> dplyr::select(label, diff, cohens_d))

# Save results
save(
  lca_wave1, lca_wave2, best_w1, best_w2, k_best,
  fit_wave1, fit_wave2,
  probs_w1, probs_w2,
  w1_classified, w2_classified, combined,
  prev_comparison,
  convergence_data, convergence_summary,
  item_results,
  truleo_items, item_labels,
  file = "data/wave_specific_lca_results.RData"
)

cat("\nSaved: data/wave_specific_lca_results.RData\n")
