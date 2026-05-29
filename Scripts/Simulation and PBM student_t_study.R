###############################################################################
#  SIMULATION: MLE for Two-Parameter Student's t-Distribution
#  Under SRS, RSS (c=1), and RSS2 (c=2)
#  Sample sizes : n = 4, 7, 10, 15, 20   
#  Parameters   : nu = 10 (fixed/known), 
#  Replications : N = 10,000
###############################################################################

#  1. EM-MLE for two-parameter t 
# E-step: w_i = (nu+1) / (nu + ((x_i - mu)/sigma)^2)
# M-step: mu = sum(w_i*x_i)/sum(w_i),  sigma^2 = sum(w_i*(x_i-mu)^2)/n

mle_t_em <- function(x, nu, tol = 1e-8, max_iter = 500) {
  n   <- length(x)
  mu  <- median(x)
  sig <- mad(x, constant = 1) + 1e-6
  for (iter in seq_len(max_iter)) {
    mu_old <- mu;  sig_old <- sig
    w   <- (nu + 1) / (nu + ((x - mu) / sig)^2)
    mu  <- sum(w * x) / sum(w)
    sig <- sqrt(pmax(sum(w * (x - mu)^2) / n, 1e-12))
    if (abs(mu - mu_old) < tol && abs(sig - sig_old) < tol) break
  }
  c(mu = mu, sigma = sig)
}

#  2. RSS / RSS2 generation
# r      = set size (= n throughout this study)
# c_cyc  = number of cycles (1 for RSS, 2 for RSS2)
# Total sample size = r * c_cyc  (= n for RSS, = 2n for RSS2)
#
# Procedure per cycle:
#   Draw r^2 units → form r sets of size r → sort each set
#   → take the i-th order statistic from the i-th set

gen_rss <- function(r, c_cyc, nu, mu_true, sigma_true) {
  out <- numeric(r * c_cyc)
  for (cyc in seq_len(c_cyc)) {
    for (i in seq_len(r)) {
      draw <- sort(mu_true + sigma_true * rt(r, df = nu))
      out[(cyc - 1L) * r + i] <- draw[i]
    }
  }
  out
}

#  3. Simulation for one (mu, sigma) scenario 
run_scenario <- function(mu_true, sigma_true, nu,
                         n_vec = c(4, 7, 10, 15, 20),
                         N     = 10000,
                         seed  = 1) {
  set.seed(seed)
  results <- data.frame()
  
  for (n in n_vec) {
    # SRS and RSS use n observations; RSS2 uses 2n observations (r=n, c=2)
    b_mu_srs  <- b_sig_srs  <- numeric(N)
    b_mu_rss  <- b_sig_rss  <- numeric(N)   # RSS  c=1, total = n
    b_mu_rss2 <- b_sig_rss2 <- numeric(N)   # RSS2 c=2, total = 2n
    
    for (rep in seq_len(N)) {
      
      # SRS (n observations)
      xs   <- mu_true + sigma_true * rt(n, df = nu)
      e_s  <- tryCatch(mle_t_em(xs, nu),
                       error = function(e) c(mu = NA, sigma = NA))
      
      # RSS: r = n, c = 1  → n observations
      x1   <- gen_rss(n, 1L, nu, mu_true, sigma_true)
      e_1  <- tryCatch(mle_t_em(x1, nu),
                       error = function(e) c(mu = NA, sigma = NA))
      
      # RSS2: r = n, c = 2 → 2n observations
      x2   <- gen_rss(n, 2L, nu, mu_true, sigma_true)
      e_2  <- tryCatch(mle_t_em(x2, nu),
                       error = function(e) c(mu = NA, sigma = NA))
      
      b_mu_srs[rep]  <- e_s["mu"]    - mu_true
      b_sig_srs[rep] <- e_s["sigma"] - sigma_true
      b_mu_rss[rep]  <- e_1["mu"]    - mu_true
      b_sig_rss[rep] <- e_1["sigma"] - sigma_true
      b_mu_rss2[rep] <- e_2["mu"]    - mu_true
      b_sig_rss2[rep]<- e_2["sigma"] - sigma_true
    }
    
    ok <- !is.na(b_mu_srs) & !is.na(b_mu_rss) & !is.na(b_mu_rss2)
    
    MSE <- function(v) mean(v[ok]^2)
    
    mse_mu_srs  <- MSE(b_mu_srs);  mse_mu_rss  <- MSE(b_mu_rss);  mse_mu_rss2  <- MSE(b_mu_rss2)
    mse_sig_srs <- MSE(b_sig_srs); mse_sig_rss <- MSE(b_sig_rss); mse_sig_rss2 <- MSE(b_sig_rss2)
    
    results <- rbind(results, data.frame(
      n         = n,
      param     = c("mu", "sigma"),
      Bias_SRS  = round(c(mean(b_mu_srs[ok]),  mean(b_sig_srs[ok])),  5),
      MSE_SRS   = round(c(mse_mu_srs,  mse_sig_srs),  5),
      Bias_RSS  = round(c(mean(b_mu_rss[ok]),  mean(b_sig_rss[ok])),  5),
      MSE_RSS   = round(c(mse_mu_rss,  mse_sig_rss),  5),
      Bias_RSS2 = round(c(mean(b_mu_rss2[ok]), mean(b_sig_rss2[ok])), 5),
      MSE_RSS2  = round(c(mse_mu_rss2, mse_sig_rss2), 5),
      RE1       = round(c(mse_mu_srs / mse_mu_rss,  mse_sig_srs / mse_sig_rss),  5),
      RE2       = round(c(mse_mu_srs / mse_mu_rss2, mse_sig_srs / mse_sig_rss2), 5)
    ))
  }
  results
}

NU   <- 10
NVEC <- c(4, 7, 10, 15, 20)
N    <- 10000

cat("  MLE Simulation: Two-Parameter Student's t (nu = 10)\n")
cat("  SRS: n obs | RSS (c=1): n obs, r=n | RSS2 (c=2): 2n obs, r=n\n")
cat("  N = 10,000 replications | n = 4, 7, 10, 15, 20\n")
cat("================================================================\n\n")

cat("Running Table 1: mu = 0,    sigma = 1    ...\n"); flush.console()
T5 <- run_scenario(0,    1,    NU, NVEC, N, seed = 101)

cat("Running Table 2: mu = 0.25, sigma = 0.75 ...\n"); flush.console()
T6 <- run_scenario(0.25, 0.75, NU, NVEC, N, seed = 202)

cat("Running Table 3: mu = 0.5,  sigma = 0.5  ...\n"); flush.console()
T7 <- run_scenario(0.5,  0.5,  NU, NVEC, N, seed = 303)

cat("Running Table 4: mu = 0.75, sigma = 1    ...\n"); flush.console()
T8 <- run_scenario(0.75, 1,    NU, NVEC, N, seed = 404)

#tables
print_table <- function(res, tbl, header) {
  cat("\n", strrep("=", 92), "\n", sep = "")
  cat(sprintf("  Table %s | %s\n", tbl, header))
  cat(strrep("-", 92), "\n")
  cat(sprintf("%-4s %-6s  %10s %10s  %10s %10s  %11s %10s  %8s %8s\n",
              "n", "Param", "Bias(SRS)", "MSE(SRS)",
              "Bias(RSS)", "MSE(RSS)", "Bias(RSS2)", "MSE(RSS2)", "RE1", "RE2"))
  cat(strrep("-", 92), "\n")
  
  prev <- -1
  for (i in seq_len(nrow(res))) {
    r <- res[i, ]
    if (r$n != prev && prev != -1) cat(strrep(".", 92), "\n")
    cat(sprintf("%-4d %-6s  %10.5f %10.5f  %10.5f %10.5f  %11.5f %10.5f  %8.5f %8.5f\n",
                r$n, r$param,
                r$Bias_SRS, r$MSE_SRS, r$Bias_RSS, r$MSE_RSS,
                r$Bias_RSS2, r$MSE_RSS2, r$RE1, r$RE2))
    prev <- r$n
  }
  cat(strrep("=", 92), "\n")
  pass <- all(res$RE1 > 1) && all(res$RE2 > 1)
  if (pass) cat("  CHECK: All RE1 and RE2 > 1 for both mu and sigma.\n\n")
  else {
    bad <- res[res$RE1 <= 1 | res$RE2 <= 1, c("n","param","RE1","RE2")]
    cat("  WARNING: RE <= 1 found:\n"); print(bad); cat("\n")
  }
}

print_table(T1, "1", "nu = 10,  mu = 0,    sigma = 1")
print_table(T2, "2", "nu = 10,  mu = 0.25, sigma = 0.75")
print_table(T3, "3", "nu = 10,  mu = 0.5,  sigma = 0.5")
print_table(T4, "4", "nu = 10,  mu = 0.75, sigma = 1")

# summary
T5$Scenario <- "Table1: mu=0,    sigma=1"
T6$Scenario <- "Table2: mu=0.25, sigma=0.75"
T7$Scenario <- "Table3: mu=0.5,  sigma=0.5"
T8$Scenario <- "Table4: mu=0.75, sigma=1"
all_res <- rbind(T1, T2, T3, T4)

write.csv(all_res, "student_t_final_results.csv", row.names = FALSE)

cat(strrep("=", 55), "\n")
cat("  OVERALL SUMMARY\n")
cat(strrep("-", 55), "\n")
cat(sprintf("  Total (n, param, scenario) rows : %d\n", nrow(all_res)))
cat(sprintf("  Rows with RE1 < 1              : %d\n", sum(all_res$RE1 < 1)))
cat(sprintf("  Rows with RE2 < 1              : %d\n", sum(all_res$RE2 < 1)))
cat(sprintf("  Minimum RE1                    : %.5f\n", min(all_res$RE1)))
cat(sprintf("  Minimum RE2                    : %.5f\n", min(all_res$RE2)))
cat(strrep("=", 55), "\n")
cat("  Saved: student_t_final_results.csv\n\n")

# PDFs and HRFs of the Two-Parameter Student's t Distribution
#  nu = 10 fixed throughout (degrees of freedom).
##################################################

# Packages
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("gridExtra",quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("scales",   quietly = TRUE)) install.packages("scales")

library(ggplot2)
library(gridExtra)
library(scales)

#  Distribution functions

# PDF of t(nu, mu, sigma)
pdf_t <- function(x, nu, mu, sigma) {
  dt((x - mu) / sigma, df = nu) / sigma
}

# CDF of t(nu, mu, sigma)
cdf_t <- function(x, nu, mu, sigma) {
  pt((x - mu) / sigma, df = nu)
}

# Hazard Rate Function: h(x) = f(x) / (1 - F(x))
hrf_t <- function(x, nu, mu, sigma) {
  f <- pdf_t(x, nu, mu, sigma)
  S <- 1 - cdf_t(x, nu, mu, sigma)
  ifelse(S < 1e-10, NA_real_, f / S)
}

#  2. Parameter combinations and plot settings 
nu <- 10

params <- list(
  list(mu = 0,    sigma = 1,    label = expression(paste(mu, " = 0,  ",    sigma, " = 1"))),
  list(mu = 0.25, sigma = 0.75, label = expression(paste(mu, " = 0.25, ", sigma, " = 0.75"))),
  list(mu = 0.5,  sigma = 0.5,  label = expression(paste(mu, " = 0.5,  ", sigma, " = 0.5"))),
  list(mu = 0.75, sigma = 1,    label = expression(paste(mu, " = 0.75, ", sigma, " = 1")))
)

colours <- c("#1b7837", "#762a83", "#d6604d", "#4393c3")
ltypes  <- c("solid", "dashed", "dotdash", "longdash")
x_min <- min(sapply(params, function(p) p$mu - 4 * p$sigma))
x_max <- max(sapply(params, function(p) p$mu + 4 * p$sigma))
x_seq <- seq(x_min, x_max, length.out = 1000)

# 3.data frames
df_pdf <- data.frame()
df_hrf <- data.frame()

for (p in params) {
  lbl <- paste0("mu=", p$mu, ", sigma=", p$sigma)
  
  pdf_vals <- pdf_t(x_seq, nu, p$mu, p$sigma)
  hrf_vals <- hrf_t(x_seq, nu, p$mu, p$sigma)
  
  df_pdf <- rbind(df_pdf,
                  data.frame(x = x_seq, y = pdf_vals, param = lbl))
  df_hrf <- rbind(df_hrf,
                  data.frame(x = x_seq, y = hrf_vals, param = lbl))
}

param_levels <- sapply(params, function(p)
  paste0("mu=", p$mu, ", sigma=", p$sigma))

df_pdf$param <- factor(df_pdf$param, levels = param_levels)
df_hrf$param <- factor(df_hrf$param, levels = param_levels)
legend_labels <- c(
  expression(paste(mu, " = 0,    ", sigma, " = 1")),
  expression(paste(mu, " = 0.25, ", sigma, " = 0.75")),
  expression(paste(mu, " = 0.5,  ", sigma, " = 0.5")),
  expression(paste(mu, " = 0.75, ", sigma, " = 1"))
)

base_theme <- theme_classic(base_size = 13) +
  theme(
    legend.position      = "bottom",
    legend.title         = element_blank(),
    legend.text          = element_text(size = 10),
    legend.key.width     = unit(1.8, "cm"),
    legend.background    = element_rect(fill = "white", colour = NA),
    axis.title           = element_text(size = 12, face = "plain"),
    axis.text            = element_text(size = 10),
    plot.title           = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.grid.minor     = element_blank(),
    panel.grid.major     = element_line(colour = "grey92", linewidth = 0.4)
  )

# 4. PDF plot
p_pdf <- ggplot(df_pdf, aes(x = x, y = y,
                            colour = param, linetype = param)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = colours, labels = legend_labels) +
  scale_linetype_manual(values = ltypes, labels = legend_labels) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "PDF",
    x     = "x",
    y     = "f(x)"
  ) +
  base_theme +
  guides(colour   = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

#5. HRF plot
hrf_x_min <- x_min + 0.5
hrf_x_max <- x_max - 0.5
df_hrf_plot <- subset(df_hrf, x >= hrf_x_min & x <= hrf_x_max & !is.na(y))
hrf_ymax <- quantile(df_hrf_plot$y, 0.995, na.rm = TRUE)

p_hrf <- ggplot(df_hrf_plot, aes(x = x, y = y,
                                 colour = param, linetype = param)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = colours, labels = legend_labels) +
  scale_linetype_manual(values = ltypes, labels = legend_labels) +
  scale_y_continuous(limits = c(0, hrf_ymax),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(hrf_x_min, hrf_x_max)) +
  labs(
    title = "HRF",
    x     = "x",
    y     = "h(x)"
  ) +
  base_theme +
  guides(colour   = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

#  6. Combine into one figure 
combined <- gridExtra::grid.arrange(
  p_pdf, p_hrf,
  ncol   = 2,
  top    = grid::textGrob(
    expression(paste("PDFs and HRFs of the Two-Parameter Student's ", italic(t),
                     " Distribution  (", nu, " = 10)")),
    gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
)

#7. Save
ggsave("student_t_pdf_hrf.png",
       plot   = combined,
       width  = 10,
       height = 5,
       dpi    = 300,
       bg     = "white")

ggsave("student_t_pdf_hrf.pdf",
       plot   = combined,
       width  = 10,
       height = 5)

cat("Saved: student_t_pdf_hrf.png  (300 dpi)\n")
cat("Saved: student_t_pdf_hrf.pdf  (vector)\n")


# =============================================================================
# Maximum Likelihood Estimation of the Two-Parameter Student's t-Distribution
# Under Simple Random Sampling (SRS) and Ranked Set Sampling (RSS)
# Dataset: Physical Body Measurements (Kaggle), anthropometric body measurements, n = 507
# Degrees of freedom: nu = 10 (fixed, as per simulation study)
# =============================================================================

pkgs <- c("MASS", "fitdistrplus", "ggplot2", "gridExtra",
          "dplyr", "knitr", "kableExtra", "moments", "nortest")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

library(MASS)
library(fitdistrplus)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(moments)
library(nortest)


# LOAD DATA

df <- read.csv("C:/Desktop/SEMESTER 2/ST/Dataset/dataset-310405444.csv")

# Full variable names
col_labels <- c(
  bia_di = "Biacromial Diameter",
  bii_di = "Biiliac Diameter",
  bit_di = "Bitrochanteric Diameter",
  che_de = "Chest Depth",
  che_di = "Chest Diameter",
  elb_di = "Elbow Diameter",
  wri_di = "Wrist Diameter",
  kne_di = "Knee Diameter",
  ank_di = "Ankle Diameter",
  sho_gi = "Shoulder Girth",
  che_gi = "Chest Girth",
  wai_gi = "Waist Girth",
  nav_gi = "Navel Girth",
  hip_gi = "Hip Girth",
  thi_gi = "Thigh Girth",
  bic_gi = "Bicep Girth",
  for_gi = "Forearm Girth",
  kne_gi = "Knee Girth",
  cal_gi = "Calf Girth",
  ank_gi = "Ankle Girth",
  wri_gi = "Wrist Girth",
  age    = "Age (years)",
  wgt    = "Weight (kg)",
  hgt    = "Height (cm)",
  sex    = "Sex (0=F, 1=M)"
)

cat("Dataset dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Sex distribution: Female =", sum(df$sex == 0), "| Male =", sum(df$sex == 1), "\n")
cat("Age range:", min(df$age), "-", max(df$age), "years\n\n")


# EM ALGORITHM FOR 2-PARAMETER t-DISTRIBUTION
#
# Treats t(nu, mu, sigma) as a scale mixture of normals (Liu & Rubin 1995)
# E-step: w_i = (nu + 1) / (nu + ((x_i - mu)/sigma)^2)
# M-step: mu  = sum(w_i * x_i) / sum(w_i)
#          sigma^2 = (1/n) * sum(w_i * (x_i - mu)^2)
# Fixed nu = 10 (as in the simulation study)
# ─────────────────────────────────────────────────────────────────────────────

em_t <- function(x, nu = 10, max_iter = 1000, tol = 1e-8) {
  n   <- length(x)
  mu  <- mean(x)
  sig <- sd(x)

  for (iter in seq_len(max_iter)) {
    mu_old  <- mu
    sig_old <- sig

    # E-step
    w <- (nu + 1) / (nu + ((x - mu) / sig)^2)

    # M-step
    mu  <- sum(w * x) / sum(w)
    sig <- sqrt(sum(w * (x - mu)^2) / n)

    if (abs(mu - mu_old) < tol && abs(sig - sig_old) < tol) break
  }
  list(mu = mu, sigma = sig, iterations = iter)
}


# KS TEST ACROSS ALL VARIABLES
# Objective 1: Identify variables that fit the 2-parameter t-distribution
#at the full population level with p-value >= 0.05

nu_fixed <- 10
cont_vars <- setdiff(names(df), "sex")

screening_results <- lapply(cont_vars, function(var) {
  x       <- df[[var]]
  fit     <- em_t(x, nu = nu_fixed)
  mu_hat  <- fit$mu
  sig_hat <- fit$sigma

  cdf_t <- function(q) pt((q - mu_hat) / sig_hat, df = nu_fixed)
  ks_res <- ks.test(x, cdf_t)

  data.frame(
    Variable  = var,
    Full_Name = col_labels[var],
    n         = length(x),
    mu_hat    = round(mu_hat, 4),
    sigma_hat = round(sig_hat, 4),
    Skewness  = round(skewness(x), 4),
    Kurtosis  = round(kurtosis(x) - 3, 4),
    KS_stat   = round(ks_res$statistic, 4),
    p_value   = round(ks_res$p.value, 4),
    Decision  = ifelse(ks_res$p.value > 0.05, "Accept H0",
                       ifelse(ks_res$p.value > 0.01, "Borderline", "Reject H0")),
    stringsAsFactors = FALSE
  )
})

screen_df <- do.call(rbind, screening_results)
screen_df <- screen_df[order(-screen_df$p_value), ]

cat("KS TEST RESULTS: ALL VARIABLES (Full Population, SRS)\n")
print(screen_df, row.names = FALSE)

cat("\n TOP CANDIDATES (p > 0.30) \n")
top_vars <- screen_df[screen_df$p_value > 0.30, "Variable"]
print(screen_df[screen_df$p_value > 0.30, ], row.names = FALSE)


#PRIMARY VARIABLES SELECTED FOR STUDY

primary_vars   <- c("cal_gi")
secondary_vars <- c("kne_gi")

#APPLICATION TO REAL DATA
# Calf Girth (cal_gi, p=0.66)
# Three sample sizes: n=6, n=10, n=17
# Three schemes per n: SRS | RSS (c=1) | RSS2 (c=2)



# Standard errors via observed Fisher information for location-scale t
# SE(mu)    = sigma / sqrt( n * nu/(nu+3) )
# SE(sigma) = sigma / sqrt( 2n * nu^2 / ((nu+1)(nu+3)) )
se_t <- function(sigma, n, nu = 10) {
  se_mu  <- sigma / sqrt(n * (nu / (nu + 3)))
  se_sig <- sigma / sqrt(2 * n * nu^2 / ((nu + 1) * (nu + 3)))
  list(se_mu = se_mu, se_sig = se_sig)
}

# RSS sample generator (one or two cycles)
# For each cycle: draw r^2 units -> r sets of size r -> select i-th order stat
generate_rss <- function(pop, r, cycles = 1) {
  out <- c()
  for (cyc in seq_len(cycles)) {
    for (i in seq_len(r)) {
      pool <- sample(pop, r * r, replace = TRUE)
      sets <- matrix(pool, nrow = r, ncol = r)
      out  <- c(out, apply(sets, 1, function(row) sort(row)[i])[i])
    }
  }
  out
}

build_bilal_table <- function(srs_obs, rss_obs, rss2_obs,
                               var_label, n_label, nu = 10) {
  fit_srs  <- em_t(srs_obs,  nu)
  fit_rss  <- em_t(rss_obs,  nu)
  fit_rss2 <- em_t(rss2_obs, nu)

  z95 <- qnorm(0.975)

  make_row <- function(fit, obs) {
    se   <- se_t(fit$sigma, length(obs), nu)
    cdf  <- function(q) pt((q - fit$mu) / fit$sigma, df = nu)
    ks   <- ks.test(obs, cdf)
    list(
      mu     = fit$mu,
      sigma  = fit$sigma,
      se_mu  = se$se_mu,
      se_sig = se$se_sig,
      lo_mu  = fit$mu    - z95 * se$se_mu,
      hi_mu  = fit$mu    + z95 * se$se_mu,
      lo_sig = fit$sigma - z95 * se$se_sig,
      hi_sig = fit$sigma + z95 * se$se_sig,
      ks_stat = ks$statistic,
      ks_pval = ks$p.value
    )
  }

  r_srs  <- make_row(fit_srs,  srs_obs)
  r_rss  <- make_row(fit_rss,  rss_obs)
  r_rss2 <- make_row(fit_rss2, rss2_obs)


  cat("\n")
  cat(paste(rep("-", 72), collapse = ""), "\n")
  cat(sprintf("Table: MLEs of mu and sigma of the t(nu=10) distribution\n"))
  cat(sprintf("Variable: %-30s  Sample size: %s\n", var_label, n_label))
  cat(paste(rep("-", 72), collapse = ""), "\n")
  cat(sprintf("%-14s %12s %12s %12s %12s %12s %12s\n",
              "", "SRS", "", "RSS (c=1)", "", "RSS2 (c=2)", ""))
  cat(sprintf("%-14s %12s %12s %12s %12s %12s %12s\n",
              "", "mu", "sigma", "mu", "sigma", "mu", "sigma"))
  cat(paste(rep("-", 72), collapse = ""), "\n")
  cat(sprintf("%-14s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
              "MLEs",
              r_srs$mu,   r_srs$sigma,
              r_rss$mu,   r_rss$sigma,
              r_rss2$mu,  r_rss2$sigma))
  cat(sprintf("%-14s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
              "SEs",
              r_srs$se_mu,   r_srs$se_sig,
              r_rss$se_mu,   r_rss$se_sig,
              r_rss2$se_mu,  r_rss2$se_sig))
  cat(sprintf("%-14s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
              "Lower (95%)",
              r_srs$lo_mu,   r_srs$lo_sig,
              r_rss$lo_mu,   r_rss$lo_sig,
              r_rss2$lo_mu,  r_rss2$lo_sig))
  cat(sprintf("%-14s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
              "Upper (95%)",
              r_srs$hi_mu,   r_srs$hi_sig,
              r_rss$hi_mu,   r_rss$hi_sig,
              r_rss2$hi_mu,  r_rss2$hi_sig))
  cat(paste(rep("-", 72), collapse = ""), "\n")
  cat(sprintf("%-14s %12s %12s %12s %12s %12s %12s\n",
              "", "SRS", "", "RSS (c=1)", "", "RSS2 (c=2)", ""))
  cat(sprintf("%-14s %24.4f %24.4f %24.4f\n",
              "KS (stat)",
              r_srs$ks_stat, r_rss$ks_stat, r_rss2$ks_stat))
  cat(sprintf("%-14s %24.4f %24.4f %24.4f\n",
              "KS (p-value)",
              r_srs$ks_pval, r_rss$ks_pval, r_rss2$ks_pval))
  cat(paste(rep("-", 72), collapse = ""), "\n")

  invisible(list(srs = r_srs, rss = r_rss, rss2 = r_rss2))
}

print_rss_matrix <- function(obs, r, label) {
  cat(sprintf("\n%s (r = %d, n = %d)\n", label, r, r))
  cat(paste(rep("-", r * 9), collapse = ""), "\n")
  m <- matrix(obs, nrow = r, ncol = r, byrow = FALSE)
  for (row in seq_len(r)) {
    cat(do.call(sprintf, c(list(paste(rep("%8.3f", r), collapse = "")), as.list(m[row, ]))), "\n")
  }
}


cat("\n\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("SECTION 5: FULL POPULATION FIT — SRS BASELINE\n")
cat(paste(rep("=", 72), collapse = ""), "\n")

for (v in c(primary_vars, secondary_vars)) {
  x    <- df[[v]]
  fit  <- em_t(x, nu_fixed)
  se   <- se_t(fit$sigma, length(x), nu_fixed)
  cdf  <- function(q) pt((q - fit$mu) / fit$sigma, df = nu_fixed)
  ks   <- ks.test(x, cdf)
  z95  <- qnorm(0.975)

  cat(sprintf("\nVariable : %s  (n = %d)\n", col_labels[v], length(x)))
  cat(sprintf("  mu_hat = %.4f   SE = %.4f   95%% CI: [%.4f, %.4f]\n",
              fit$mu, se$se_mu,
              fit$mu - z95*se$se_mu, fit$mu + z95*se$se_mu))
  cat(sprintf("  sig_hat= %.4f   SE = %.4f   95%% CI: [%.4f, %.4f]\n",
              fit$sigma, se$se_sig,
              fit$sigma - z95*se$se_sig, fit$sigma + z95*se$se_sig))
  cat(sprintf("  KS stat= %.4f   p-value = %.4f\n",
              ks$statistic, ks$p.value))
}

#   Table A  = RSS scheme (c=1) observed matrix    
#   Table B  = RSS2 scheme (c=2) observed matrix   
#   Table C  = MLE results: SRS | RSS | RSS2       


set.seed(42) 

sample_sizes <- c(6, 10, 17)   

all_results <- list()  

for (v in primary_vars) {

  pop <- df[[v]] 

  cat("\n\n")
  cat(paste(rep("=", 72), collapse = ""), "\n")
  cat(sprintf("VARIABLE: %s\n", toupper(col_labels[v])))
  cat(paste(rep("=", 72), collapse = ""), "\n")

  var_results <- list()

  for (n in sample_sizes) {

    cat(sprintf("\n\n>>> n = %d <<<\n", n))

    # SRS sample
    srs_obs <- sample(pop, n, replace = FALSE)
    cat(sprintf("\nSRS observed values (n=%d):\n", n))
    cat(do.call(sprintf, c(list(paste(rep("%8.3f", n), collapse = "")), as.list(srs_obs))), "\n")

    # RSS (c=1): r=n, one cycle 
    rss_obs <- generate_rss(pop, r = n, cycles = 1)

    # r x r ranked matrix
    cat(sprintf("\nTable (RSS c=1, r=%d): ranked set scheme from t(nu=10, mu, sigma)\n", n))
    cat(sprintf("One Cycle — diagonal entries are the selected observations\n"))
    set.seed(42 + n)  
    rss_matrix <- matrix(NA, nrow = n, ncol = n)
    for (i in seq_len(n)) {
      pool <- sample(pop, n * n, replace = TRUE)
      rss_matrix[i, ] <- sort(pool[((i-1)*n + 1):(i*n)])
    }
    rss_obs <- diag(rss_matrix)  # i-th order stat from i-th set
    for (row in seq_len(n)) {
      cat(do.call(sprintf, c(list(paste(rep("%8.3f", n), collapse = "")), as.list(rss_matrix[row, ]))))
      cat(sprintf("  <- selected: %.3f\n", rss_matrix[row, row]))
    }

    # RSS2 (c=2): two cycles
    # Unequal cycles: cycle1 has r=ceiling(n/2), cycle2 has r=floor(n/2)
    r1 <- ceiling(n / 2)
    r2 <- floor(n / 2)

    set.seed(42 + n + 100)
    rss2_obs <- c()

    cat(sprintf("\nTable (RSS2 c=2): two-cycle ranked set scheme\n"))
    cat(sprintf("Cycle 1 (r=%d):\n", r1))
    for (i in seq_len(r1)) {
      pool <- sample(pop, r1 * r1, replace = TRUE)
      row_sorted <- sort(pool[((i-1)*r1 + 1):(i*r1)])
      cat(do.call(sprintf, c(list(paste(rep("%8.3f", r1), collapse = "")), as.list(row_sorted))))
      cat(sprintf("  <- selected: %.3f\n", row_sorted[i]))
      rss2_obs <- c(rss2_obs, row_sorted[i])
    }
    cat(sprintf("Cycle 2 (r=%d):\n", r2))
    for (i in seq_len(r2)) {
      pool <- sample(pop, r2 * r2, replace = TRUE)
      row_sorted <- sort(pool[((i-1)*r2 + 1):(i*r2)])
      cat(do.call(sprintf, c(list(paste(rep("%8.3f", r2), collapse = "")), as.list(row_sorted))))
      cat(sprintf("  <- selected: %.3f\n", row_sorted[i]))
      rss2_obs <- c(rss2_obs, row_sorted[i])
    }

    # MLEs + SEs + CI + KS 
    res <- build_bilal_table(
      srs_obs  = srs_obs,
      rss_obs  = rss_obs,
      rss2_obs = rss2_obs,
      var_label = col_labels[v],
      n_label   = paste0("n = ", n)
    )

    var_results[[paste0("n", n)]] <- res
  }

  all_results[[v]] <- var_results
}

# Empirical CDF, PDF, PP-plot 


cat("\n\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

plot_bilal_figure <- function(x, var_label, nu = 10) {
  fit   <- em_t(x, nu)
  mu    <- fit$mu
  sig   <- fit$sigma
  cdf_t <- function(q) pt((q - mu) / sig, df = nu)
  ks    <- ks.test(x, cdf_t)

  x_sort  <- sort(x)
  n       <- length(x)
  emp_cdf <- (1:n) / n
  fit_cdf <- pt((x_sort - mu) / sig, df = nu)

  xg   <- seq(min(x) - 1.5*sig, max(x) + 1.5*sig, length.out = 500)
  dens <- dt((xg - mu) / sig, df = nu) / sig

  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0,0,2,0))

  # Panel 1: Empirical CDF vs fitted CDF
  plot(x_sort, emp_cdf, type = "s", lwd = 2, col = "black",
       xlab = "x", ylab = "F(x)",
       main = paste0("Empirical vs Fitted CDF"),
       ylim = c(0, 1))
  lines(xg, pt((xg - mu) / sig, df = nu),
        col = "steelblue", lwd = 2, lty = 2)
  legend("bottomright", bty = "n", cex = 0.85,
         legend = c("Empirical", paste0("t(\u03bd=10)")),
         col = c("black", "steelblue"), lwd = 2, lty = c(1, 2))


  # Panel 2: Histogram + fitted PDF
  hist(x, breaks = 20, freq = FALSE,
       col = "#DBEAFE", border = "white",
       xlab = "x", ylab = "f(x)",
       main = "Histogram + Fitted PDF",
       ylim = c(0, max(dens) * 1.1))
  lines(xg, dens, col = "steelblue", lwd = 2)
  legend("topright", bty = "n", cex = 0.85,
         legend = c(sprintf("\u03bc\u0302=%.4f", mu),
                    sprintf("\u03c3\u0302=%.4f", sig)),
         col = c("steelblue", "steelblue"), lwd = 2)

  # Panel 3: P-P plot
  plot(fit_cdf, emp_cdf, pch = 16, cex = 0.6, col = "steelblue",
       xlab = "Theoretical probability",
       ylab = "Empirical probability",
       main = paste0("P-P Plot for t(\u03bd=10)"),
       xlim = c(0,1), ylim = c(0,1))
  abline(0, 1, lty = 2, lwd = 1.5, col = "black")
  legend("topleft", bty = "n", cex = 0.85,
         legend = sprintf("KS = %.4f\np = %.4f", ks$statistic, ks$p.value))

  mtext(paste("Estimated CDF, PDF and PP-plot"),
        outer = TRUE, cex = 1, font = 2)

  par(mfrow = c(1,1), oma = c(0,0,0,0))
}

for (v in primary_vars) {
  plot_bilal_figure(df[[v]], col_labels[v])
  dev.copy(png,
           filename = paste0("figure2_", v, ".png"),
           width = 1400, height = 480, res = 130)
  dev.off()
  cat("Saved: figure2_", v, ".png\n", sep = "")
}





