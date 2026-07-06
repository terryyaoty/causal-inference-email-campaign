# ============================================================================
# MSIN0303 Group Assignment 
# ============================================================================
# Preparation

# ---- Setup ----
# install.packages(c("tidyverse", "fixest", "broom", "modelsummary"))
library(tidyverse)
library(fixest)
library(broom)
library(modelsummary)
library(grf)
library(scales)

# ---- Load data ----
getwd()
setwd("C:/Users/步步高点读机/Desktop/UCL MS/Term 3/Causal Inference/Group/")
df <- read.csv("hillstrom email.csv.csv")

set.seed(42)
# ============================================================
# Q0. Problem Framing + Experimental Validation
# MSIN0303 Group Project  |  Hillstrom Email Campaign
# ============================================================
# ============================================================
# 1. Data Overview
# ============================================================

# Check structure
glimpse(df)

# Group sizes
table(df$segment)

# Outcome means by group
df %>%
  group_by(segment) %>%
  summarise(
    n          = n(),
    avg_spend  = mean(spend),
    avg_visit  = mean(visit),
    avg_conv   = mean(conversion),
    avg_recency = mean(recency),
    avg_history = mean(history)
  )

# ============================================================
# 2. Balance Table
# Compare pre-treatment covariate means across three groups
# ============================================================

covariates <- c("recency", "history", "mens", "womens", "newbie")

balance_table <- df %>%
  group_by(segment) %>%
  summarise(across(all_of(covariates), mean)) %>%
  pivot_longer(-segment, names_to = "variable", values_to = "mean") %>%
  pivot_wider(names_from = segment, values_from = mean)

print(balance_table)

# ============================================================
# 3. SMD Calculation
# SMD = (mean_treated - mean_control) / sd_control
# Rule of thumb: |SMD| < 0.1 indicates good balance
# ============================================================

compute_smd <- function(data, var, treated_label, control_label = "No E-Mail") {
  treated <- data %>% filter(segment == treated_label) %>% pull(!!sym(var))
  control <- data %>% filter(segment == control_label) %>% pull(!!sym(var))
  smd <- (mean(treated) - mean(control)) / sd(control)
  return(smd)
}

# Compute SMD for all covariates against No E-Mail control
smd_results <- expand_grid(
  variable  = covariates,
  treatment = c("Mens E-Mail", "Womens E-Mail")
) %>%
  rowwise() %>%
  mutate(SMD = compute_smd(df, variable, treatment)) %>%
  ungroup()

print(smd_results)

# ============================================================
# 4. Love Plot (SMD visualisation)
# ============================================================

ggplot(smd_results, aes(x = SMD, y = variable, color = treatment, shape = treatment)) +
  geom_point(size = 3) +
  geom_vline(xintercept =  0,   linetype = "solid",  color = "black") +
  geom_vline(xintercept =  0.1, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = -0.1, linetype = "dashed", color = "gray50") +
  labs(
    title    = "Balance Check: Standardised Mean Differences",
    subtitle = "Compared to No E-Mail control | Dashed lines = ±0.1 threshold",
    x        = "Standardised Mean Difference (SMD)",
    y        = "Covariate",
    color    = "Treatment",
    shape    = "Treatment"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("balance_smd_plot.png", width = 8, height = 5, dpi = 150)

# ============================================================
# 5. Balance Summary
# ============================================================

smd_results %>%
  mutate(balanced = abs(SMD) < 0.1) %>%
  group_by(treatment) %>%
  summarise(
    n_covariates = n(),
    n_balanced   = sum(balanced),
    max_abs_SMD  = max(abs(SMD))
  )

# ============================================================
# Q1 — Average Treatment Effects (ATE)
# ============================================================
#
# Goal: estimate the causal effect of each email type on
# customer behaviour, with 95% confidence intervals.
#
# Identification: Part 1 verified randomisation was
# successful (max |SMD| < 0.01). With clean random assignment,
# the OLS coefficient is an unbiased estimate of the ATE.
#
# Three outcomes:
#   spend       — primary (in dollars, comparable to £0.30 cost)
#   visit       — supporting evidence (binary)
#   conversion  — supporting evidence (binary)
#
# Methods used:
#   feols()         OLS regression
#   modelsummary()  side-by-side regression tables
#   coef(), confint()  extract estimates and CIs
# ============================================================

# Sanity check vs Yui's group sizes
table(df$segment)

# Set "No E-Mail" as the reference level so the OLS coefficients
# read out directly as "Mens vs Control" and "Womens vs Control"
df <- df %>%
  mutate(segment = factor(segment,
                          levels = c("No E-Mail", "Mens E-Mail", "Womens E-Mail")))


# ============================================================
# 1. Three OLS regressions (one per outcome)
# ============================================================
# In a clean RCT, OLS of the outcome on the treatment indicator
# recovers the ATE directly. No covariates needed .

ate_spend <- feols(spend      ~ segment, data = df)
ate_visit <- feols(visit      ~ segment, data = df)
ate_conv  <- feols(conversion ~ segment, data = df)


# ============================================================
# 2. Compare all three side-by-side
# ============================================================

modelsummary(
  list(
    "Visit (prob.)"      = ate_visit,
    "Conversion (prob.)" = ate_conv,
    "Spend ($)"          = ate_spend
  ),
  stars    = TRUE,
  gof_omit = "IC|Log|Adj|Pseudo|Within"
)


# ============================================================
# 3. Build a tidy table of ATEs with 95% CIs (for plotting)
# ============================================================

extract_ate <- function(model, outcome_name) {
  est <- coef(model)
  ci  <- confint(model)
  tibble(
    outcome   = outcome_name,
    treatment = names(est),
    estimate  = as.numeric(est),
    conf.low  = ci[, 1],
    conf.high = ci[, 2]
  ) %>%
    filter(treatment != "(Intercept)") %>%
    mutate(treatment = recode(treatment,
                              "segmentMens E-Mail"   = "Mens E-Mail",
                              "segmentWomens E-Mail" = "Womens E-Mail"))
}

ate_table <- bind_rows(
  extract_ate(ate_visit, "visit (prob)"),
  extract_ate(ate_conv,  "conversion (prob)"),
  extract_ate(ate_spend, "spend ($)")
) %>%
  # Funnel order: visit -> conversion -> spend
  mutate(outcome = factor(outcome,
                          levels = c("visit (prob)",
                                     "conversion (prob)",
                                     "spend ($)")))

print(ate_table)
write_csv(ate_table, "q1_ate_table.csv")


# ============================================================
# 4. Main chart — ATE on spend (hero shot)
# ============================================================
# Both ATEs sit clearly above 0 and above the £0.30 send cost.

ggplot(ate_table %>% filter(outcome == "spend ($)"),
       aes(x = treatment, y = estimate, color = treatment)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.15, linewidth = 1) +
  geom_hline(yintercept = 0,    linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0.30, linetype = "dotted", color = "red") +
  annotate("text", x = 2.45, y = 0.33, label = "send cost £0.30",
           color = "red", size = 3.5, hjust = 1) +
  labs(
    title    = "Causal effect of email on customer spend",
    subtitle = "2-week post-campaign window | OLS with 95% CIs",
    x = NULL,
    y = "Estimated lift in spend ($)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("q1_ate_spend_plot.png", width = 8, height = 5, dpi = 150)


# ============================================================
# 5. Robustness chart — all three outcomes
# ============================================================
# Funnel logic: email -> visit -> conversion -> spend.
# If the effect is real, all three should move in the same
# direction.

ggplot(ate_table,
       aes(x = treatment, y = estimate, color = treatment)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.15, linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ factor(outcome,
                      levels = c("visit (prob)",
                                 "conversion (prob)",
                                 "spend ($)")),
             scales = "free_y") +
  labs(
    title    = "ATEs across all three outcomes",
    subtitle = "Consistent positive lift across the funnel | 95% CIs",
    x = NULL, y = "Estimated effect"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))

ggsave("q1_ate_all_outcomes_plot.png", width = 11, height = 5, dpi = 150)


# ============================================================
# Q1 Done. Outputs saved to working directory:
#   q1_ate_table.csv             ATE estimates with 95% CIs
#   q1_ate_spend_plot.png        Hero shot (spend, primary)
#   q1_ate_all_outcomes_plot.png Robustness (all 3 outcomes)
# ============================================================


# ============================================================
# Q2 — Heterogeneous Treatment Effects & Targeting Policy
# ============================================================
#
# Builds on Q1's ATE results. Q1 told us emails work on average.
# Q2 asks: who should we email, with what creative, given 30% budget?
#
# Method (Lecture 10 cookbook):
#   - Train two binary causal forests (grf package):
#       Forest M: Mens   vs No E-Mail
#       Forest W: Womens vs No E-Mail
#   - Predict τ̂_M(x), τ̂_W(x) for every customer
#   - For each customer pick the better creative
#   - Rank by net_value (= best τ̂ − £0.30 cost)
#   - Target top 30% with net_value > 0
#
# Assumes Q1 script has been run, so `df` is in memory with
# segment as a factor (No E-Mail as reference level).
# ============================================================
# ============================================================
# 1. Build feature matrix
# ============================================================
# grf needs a numeric matrix. We dummify the categorical
# variables (zip_code, channel, history_segment) using
# model.matrix(). The "- 1" drops the intercept.

X_full <- model.matrix(
  ~ recency + history + mens + womens + newbie +
    zip_code + channel + history_segment - 1,
  data = df
)

cat("Feature matrix dimensions:", dim(X_full), "\n")


# ============================================================
# 2. Build the two binary sub-samples
# ============================================================
# causal_forest() only handles binary treatment, so we split:
#   Forest M: keep Mens + No E-Mail (drop Womens)
#   Forest W: keep Womens + No E-Mail (drop Mens)

idx_M <- df$segment %in% c("Mens E-Mail",   "No E-Mail")
idx_W <- df$segment %in% c("Womens E-Mail", "No E-Mail")

X_M <- X_full[idx_M, ]
Y_M <- df$spend[idx_M]
W_M <- as.integer(df$segment[idx_M] == "Mens E-Mail")  # 1 = treated, 0 = control

X_W <- X_full[idx_W, ]
Y_W <- df$spend[idx_W]
W_W <- as.integer(df$segment[idx_W] == "Womens E-Mail")


# ============================================================
# 3. Fit the two causal forests
# ============================================================
# Each forest takes 2-5 minutes. 2000 trees follows Lecture 10
# cookbook defaults (more trees = more stable τ̂ estimates).

cat("Fitting Forest M...\n")
forest_M <- causal_forest(X_M, Y_M, W_M, num.trees = 2000, seed = 42)

cat("Fitting Forest W...\n")
forest_W <- causal_forest(X_W, Y_W, W_W, num.trees = 2000, seed = 42)


# ============================================================
# 4. Sanity check — forest ATE should ≈ Q1 OLS ATE
# ============================================================
# If these match, the forest is well-calibrated.
# Recall from Q1: Mens ATE ≈ $0.77, Womens ATE ≈ $0.42.

ate_forest_M <- average_treatment_effect(forest_M)
ate_forest_W <- average_treatment_effect(forest_W)

cat("\n--- Sanity check ---\n")
cat(sprintf("Forest M ATE: %.3f (SE %.3f)  | Q1 OLS was ~0.77\n",
            ate_forest_M[1], ate_forest_M[2]))
cat(sprintf("Forest W ATE: %.3f (SE %.3f)  | Q1 OLS was ~0.42\n",
            ate_forest_W[1], ate_forest_W[2]))


# ============================================================
# 5. Predict τ̂ for ALL 64,000 customers
# ============================================================
# Each forest was trained on 2/3 of the data, but the learned
# τ(x) function can be applied to any customer. We need τ̂ for
# everyone (including the people who got the *other* email
# arm) to build a policy.

tau_M_all <- predict(forest_M, newdata = X_full)$predictions
tau_W_all <- predict(forest_W, newdata = X_full)$predictions

df_policy <- df %>%
  mutate(tau_M = tau_M_all,
         tau_W = tau_W_all)


# ============================================================
# 6. Diagnostics — heterogeneity in τ̂
# ============================================================

# (a) Histogram: how much do treatment effects vary across customers?
df_tau_long <- df_policy %>%
  select(tau_M, tau_W) %>%
  pivot_longer(everything(), names_to = "email_type", values_to = "tau") %>%
  mutate(email_type = factor(recode(email_type,
                                    tau_M = "Mens email",
                                    tau_W = "Womens email"),
                             levels = c("Mens email", "Womens email")))

ggplot(df_tau_long, aes(x = tau, fill = email_type)) +
  geom_histogram(bins = 60, alpha = 0.85) +
  facet_wrap(~ email_type, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = 0,    linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0.30, linetype = "solid",  color = "red") +
  annotate("text", x = 0.32, y = Inf, vjust = 1.5,
           label = "£0.30 = send cost", color = "red", size = 3.5, hjust = 0) +
  labs(
    title    = "Distribution of estimated treatment effects τ̂(x)",
    subtitle = "Customers to the right of the red line are profitable to target",
    x = "Estimated uplift in spend ($)",
    y = "Number of customers"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("q2_tau_histograms.png", width = 9, height = 6, dpi = 150)


# (b) Variable importance — what drives heterogeneity?
vi_df <- tibble(
  feature      = colnames(X_full),
  importance_M = as.numeric(variable_importance(forest_M)),
  importance_W = as.numeric(variable_importance(forest_W))
) %>%
  arrange(desc(importance_M))

print(vi_df)
write_csv(vi_df, "q2_variable_importance.csv")

# Variable importance bar chart (top 8 features)
vi_plot_data <- vi_df %>%
  slice_max(importance_M, n = 8) %>%
  pivot_longer(c(importance_M, importance_W),
               names_to = "email_type",
               values_to = "importance") %>%
  mutate(
    email_type = recode(email_type,
                        importance_M = "Mens email",
                        importance_W = "Womens email"),
    feature = factor(feature, levels = rev(unique(feature)))
  )

ggplot(vi_plot_data,
       aes(x = importance, y = feature, fill = email_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Mens email"   = "#F8766D",
                               "Womens email" = "#00BFC4")) +
  labs(
    title    = "What drives heterogeneity in email response?",
    subtitle = "Variable importance from both causal forests",
    x = "Importance (share of splits)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank())

ggsave("q2_variable_importance_plot.png", width = 9, height = 6, dpi = 150)
# (c) Best linear projection — human-readable HTE
# "Each extra unit of X changes τ by β" — useful for the slide.
blp_M <- best_linear_projection(forest_M, X_M)
blp_W <- best_linear_projection(forest_W, X_W)

cat("\n--- Best Linear Projection: Mens email ---\n");  print(blp_M)
cat("\n--- Best Linear Projection: Womens email ---\n"); print(blp_W)

# Save to CSV (grf returns a matrix object — convert explicitly)
blp_to_df <- function(blp) {
  data.frame(
    term      = rownames(blp),
    estimate  = blp[, "Estimate"],
    std_error = blp[, "Std. Error"],
    t_value   = blp[, "t value"],
    p_value   = blp[, "Pr(>|t|)"],
    row.names = NULL
  )
}
write_csv(blp_to_df(blp_M), "q2_blp_mens.csv")
write_csv(blp_to_df(blp_W), "q2_blp_womens.csv")


# ============================================================
# 7. Build targeting policy (the Step 5 logic in code)
# ============================================================

COST       <- 0.30                       # send cost per email
BUDGET_PCT <- 0.30                       # max share to target
N_CAP      <- ceiling(nrow(df) * BUDGET_PCT)

df_policy <- df_policy %>%
  mutate(
    # For each customer pick the better creative
    best_email = if_else(tau_M >= tau_W, "Mens E-Mail", "Womens E-Mail"),
    best_tau   = pmax(tau_M, tau_W),
    net_value  = best_tau - COST          # expected $ lift net of cost
  ) %>%
  arrange(desc(net_value)) %>%
  mutate(rank = row_number())

# Target if (i) ranked in top 30% AND (ii) net_value > 0
df_policy <- df_policy %>%
  mutate(
    targeted       = (rank <= N_CAP) & (net_value > 0),
    assigned_email = if_else(targeted, best_email, "No E-Mail")
  )

# Summary
policy_summary <- df_policy %>%
  count(assigned_email) %>%
  mutate(pct = n / sum(n))
print(policy_summary)

cat(sprintf("\nPolicy targets %d customers (%.1f%% of base)\n",
            sum(df_policy$targeted),
            100 * mean(df_policy$targeted)))


# ============================================================
# 8. Lift curve — visualises why targeting beats blanket
# ============================================================

df_lift <- df_policy %>%
  arrange(desc(net_value)) %>%
  mutate(
    pct_targeted = row_number() / n(),
    cum_net      = cumsum(net_value)
  )

# Reference points: peak (unconstrained optimum) and 30% cap
peak_idx <- which.max(df_lift$cum_net)
peak_pct <- df_lift$pct_targeted[peak_idx]
peak_net <- df_lift$cum_net[peak_idx]

ggplot(df_lift, aes(x = pct_targeted, y = cum_net)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_vline(xintercept = 0.30,     linetype = "dashed", color = "red") +
  geom_vline(xintercept = peak_pct, linetype = "dotted", color = "darkgreen") +
  annotate("text", x = 0.31, y = peak_net * 0.4, hjust = 0,
           label = "30% budget cap", color = "red", size = 3.8) +
  annotate("text", x = peak_pct + 0.01, y = peak_net * 0.8, hjust = 0,
           label = sprintf("Unconstrained optimum: %.0f%%", 100 * peak_pct),
           color = "darkgreen", size = 3.8) +
  scale_x_continuous(labels = percent_format()) +
  scale_y_continuous(labels = dollar_format()) +
  labs(
    title    = "Cumulative net revenue as we target more customers",
    subtitle = "Customers ordered by personalised net uplift (τ̂ − £0.30)",
    x = "Share of customer base targeted",
    y = "Cumulative expected net revenue"
  ) +
  theme_minimal(base_size = 13)

ggsave("q2_lift_curve.png", width = 9, height = 5, dpi = 150)


# ============================================================
# 9. Export policy table
# ============================================================

write_csv(df_policy, "q2_policy_table.csv")

cat("\n============================================================\n")
cat("Done. Outputs saved to working directory:\n")
cat("  q2_tau_histograms.png            — heterogeneity visualised\n")
cat("  q2_variable_importance.csv       — what drives heterogeneity (table)\n")
cat("  q2_variable_importance_plot.png  — what drives heterogeneity (chart)\n")
cat("  q2_blp_mens.csv / q2_blp_womens.csv — BLP (interpretable HTE)\n")
cat("  q2_lift_curve.png                — targeting vs blanket\n")
cat("  q2_policy_table.csv              — per-customer policy → hand to Member D\n")
cat("============================================================\n")
