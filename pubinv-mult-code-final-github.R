# Code for the paper "Achieving geoeconomic goals by boosting the economy without raising the public debt ratio?"
rm(list = ls())

## -- minimal packages
suppressPackageStartupMessages({
  library(data.table)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(here)
  library(lpirfs)
})

## Load & prepare the data
dt <- fread(here("data_pubinv_final.csv"), check.names = FALSE)

# calculate growth in private investment ratio
dt <- ddply(dt, "ccode", transform, INVGDP_diff = c(NA, diff(INVGDP))) # calculate real growth

# Functions
# cumulative IRFs
cum_irf <- function(obj, scale = 1) {
  mu <- as.numeric(obj$irf_panel_mean) * scale
  lo <- as.numeric(obj$irf_panel_low)  * scale
  hi <- as.numeric(obj$irf_panel_up)   * scale
  tibble(
    h   = 0:(length(mu) - 1),
    cum = cumsum(mu),
    lo  = cumsum(lo),
    hi  = cumsum(hi)
  )
}

# standard error bands
se_from_bands <- function(mid, lo, hi) {
  0.5 * ((hi - mid) + (mid - lo))
}

# multipliers
mult_from_ratio <- function(gdp_lp, ratio_lp, r_share) {
  gdp_c   <- cum_irf(gdp_lp,   scale = 100) %>% dplyr::rename(X = cum, X_lo = lo, X_hi = hi)
  ratio_c <- cum_irf(ratio_lp, scale =   1) %>% dplyr::rename(R = cum, R_lo = lo, R_hi = hi)
  out <- dplyr::left_join(gdp_c, ratio_c, by = "h")
  
  se_X <- se_from_bands(out$X, out$X_lo, out$X_hi)
  se_R <- se_from_bands(out$R, out$R_lo, out$R_hi)
  
  # adjust for impact of change in GDP on investment ratio
  den_pp <- out$R + r_share * out$X
  Dsafe  <- pmax(den_pp, 1e-12)
  
  mult <- out$X / Dsafe
  dMdX <- out$R / (Dsafe^2)
  dMdR <- -out$X / (Dsafe^2)
  #Delta method
  se_M <- sqrt((dMdX^2) * (se_X^2) + (dMdR^2) * (se_R^2))
  
  tibble(h = out$h, multiplier = mult, lo_1se = mult - se_M, hi_1se = mult + se_M)
}

plot_mult <- function(tbl, ttl = "Real GDP: cumulative investment multiplier") {
  ggplot(tbl, aes(h, multiplier)) +
    geom_ribbon(aes(ymin = lo_1se, ymax = hi_1se), fill = "grey70", alpha = 0.4) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.8, color = "#E41A1C") +
    labs(title = ttl, x = "Years after the shock", y = "multiplier") +
    theme_minimal(base_size = 12)
}

plot_cum <- function(obj, scale = 1, ttl = "", ylab = "") {
  cc <- cum_irf(obj, scale = scale)
  ggplot(cc, aes(h, cum)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.4) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.8, color = "#E41A1C") +
    labs(title = ttl, x = "Years after the shock", y = ylab) +
    theme_minimal(base_size = 12)
}

## Baseline LPs and the multiplier
# GDP LP
lp_gdp <- lp_lin_panel(
  data_set       = dt,
  endog_data     = "log_RGDP",
  shock          = "forecasterror",
  diff_shock     = FALSE,
  panel_model    = "within",
  panel_effect   = "twoways",
  robust_cov     = "vcovSCC",
  c_exog_data    = c("REER"),
  l_exog_data    = c("growth_RGDP","PDEBT","forecasterror","NOMLRATE"),
  lags_exog_data = 2,
  confint        = 1,
  hor            = 4
)

# Public investment LP
lp_ratio <- lp_lin_panel(
  data_set       = dt,
  endog_data     = "PUBINVRATIO",
  shock          = "forecasterror",
  diff_shock     = FALSE,
  panel_model    = "within",
  panel_effect   = "twoways",
  robust_cov     = "vcovSCC",
  c_exog_data    = c("REER"),
  l_exog_data    = c("growth_RGDP","PDEBT","forecasterror","NOMLRATE"),
  lags_exog_data = 2,
  confint        = 1,
  hor            = 4
)

# Average share rbar
rbar <- mean(dt$PUBINVRATIO, na.rm = TRUE)
if (!is.na(rbar) && rbar > 1) rbar <- rbar / 100

mult_base <- mult_from_ratio(lp_gdp, lp_ratio, rbar)
print(mult_base)

# plot multiplier
g_mult <- plot_mult(mult_base)

## Other responses (private inv., debt, unemployment)
lp_priv <- lp_lin_panel(
  data_set       = dt,
  endog_data     = "INVGDP",
  shock          = "forecasterror",
  diff_shock     = FALSE,
  panel_model    = "within",
  panel_effect   = "twoways",
  robust_cov     = "vcovSCC",
  c_exog_data    = c("REER"),
  l_exog_data    = c("growth_RGDP","PDEBT","forecasterror","NOMLRATE","INVGDP_diff"),
  lags_exog_data = 2,
  confint        = 1,
  hor            = 4
)

lp_debt <- lp_lin_panel(
  data_set       = dt,
  endog_data     = "PDEBT",
  shock          = "forecasterror",
  diff_shock     = FALSE,
  panel_model    = "within",
  panel_effect   = "twoways",
  robust_cov     = "vcovSCC",
  c_exog_data    = c("REER"),
  l_exog_data    = c("growth_RGDP","PDEBT","forecasterror","NOMLRATE"),
  lags_exog_data = 2,
  confint        = 1,
  hor            = 4
)

lp_unemp <- lp_lin_panel(
  data_set       = dt,
  endog_data     = "UNRATE",
  shock          = "forecasterror",
  diff_shock     = FALSE,
  panel_model    = "within",
  panel_effect   = "twoways",
  robust_cov     = "vcovSCC",
  c_exog_data    = c("REER"),
  l_exog_data    = c("growth_RGDP","PDEBT","forecasterror","NOMLRATE","UNRATE"),
  lags_exog_data = 2,
  confint        = 1,
  hor            = 4
)

# plots
g_priv <- plot_cum(lp_priv,  scale = 1, ttl = "Private investment ratio", ylab = "percentage points")
g_debt <- plot_cum(lp_debt,  scale = 1, ttl = "Public debt ratio",        ylab = "percentage points")
g_unem <- plot_cum(lp_unemp, scale = 1, ttl = "Unemployment rate",        ylab = "percentage points")

layout_base <- ggarrange(g_mult, g_unem, g_priv, g_debt,
                         labels = c("A)", "B)", "C)", "D)"),
                         ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave("econometric_results_BMF.jpg", plot = layout_base, width = 10, height = 6, bg = "white")

## Robustness checks
# Updated helper to allow custom dataset and rbar (defaults preserve prior behavior)
run_mult <- function(ctrls, lags, data_set = dt, r_share = rbar) {
  gdp_lp <- lp_lin_panel(
    data_set       = data_set,
    endog_data     = "log_RGDP",
    shock          = "forecasterror",
    diff_shock     = FALSE,
    panel_model    = "within",
    panel_effect   = "twoways",
    robust_cov     = "vcovSCC",
    c_exog_data    = c("REER"),
    l_exog_data    = ctrls,
    lags_exog_data = lags,
    confint        = 1,
    hor            = 4
  )
  ratio_lp <- lp_lin_panel(
    data_set       = data_set,
    endog_data     = "PUBINVRATIO",
    shock          = "forecasterror",
    diff_shock     = FALSE,
    panel_model    = "within",
    panel_effect   = "twoways",
    robust_cov     = "vcovSCC",
    c_exog_data    = c("REER"),
    l_exog_data    = ctrls,
    lags_exog_data = lags,
    confint        = 1,
    hor            = 4
  )
  mult_from_ratio(gdp_lp, ratio_lp, r_share)
}

## Robustness 1: replace growth_RGDP with OUTPUTGAP
ctrl_rob1 <- c("OUTPUTGAP","PDEBT","forecasterror","NOMLRATE")
mult_rob1 <- run_mult(ctrl_rob1, lags = 2)
print(mult_rob1)

## Robustness 2: increase lag length to 3 (baseline controls)
ctrl_base <- c("growth_RGDP","PDEBT","forecasterror","NOMLRATE")
mult_rob2 <- run_mult(ctrl_base, lags = 3)
print(mult_rob2)

## Robustness 3: add PRIMARYBAL to baseline controls
ctrl_rob3 <- c("growth_RGDP","PDEBT","forecasterror","NOMLRATE","PRIMARYBAL")
mult_rob3 <- run_mult(ctrl_rob3, lags = 2)
print(mult_rob3)

## Robustness 4: exclude Ireland
dt_no_irl <- subset(dt, ccode != "IRL")
rbar_no_irl <- mean(dt_no_irl$PUBINVRATIO, na.rm = TRUE)
if (!is.na(rbar_no_irl) && rbar_no_irl > 1) rbar_no_irl <- rbar_no_irl / 100
mult_rob4 <- run_mult(ctrl_base, lags = 2, data_set = dt_no_irl, r_share = rbar_no_irl)
print(mult_rob4)
