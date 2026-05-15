# ====================================================================
# SECTION 3:  PC — Portfolio Choice  (G&L ch. 4.5)
# ====================================================================
# Adds government bills (B) alongside cash (H). Central bank sets
# the interest rate r exogenously (HORIZONTALIST closure) and passively
# absorbs whatever bills households don't want to hold.
#
# Wealth:        V = H + B^h
# Tobin block:   B^h / V = lam0 + lam1 * r - lam2 * (YD / V)
#                (cash share is the complement, by accounting)
#
# IMPORTANT: in code we use the LINEARISED form
#   B^h = lam0 * V + lam1 * r * V - lam2 * YD
# which is algebraically identical but doesn't divide by V — so the
# model is well-defined even when V starts at 0.
#
# State dim.: 2  (V, Bh; lagged H, Bs, Bcb are linked by identities)
# --------------------------------------------------------------------

pc_model <- function(cur, prev, p, t) {
  r <- p$r
  G <- p$G
  
  # --- Goods market and household disposable income ---
  Y  <- cur[["C"]] + G
  Tx <- p$theta * (Y + r * prev[["Bh"]])              # tax on Y + coupon income
  YD <- Y - Tx + r * prev[["Bh"]]                     # YD includes coupons
  V  <- prev[["V"]] + YD - cur[["C"]]                 # wealth identity
  C  <- p$alpha1 * YD + p$alpha2 * prev[["V"]]
  
  # --- Tobin portfolio block (linearised, see header) ---
  Bh <- p$lam0 * V + p$lam1 * r * V - p$lam2 * YD
  H  <- V - Bh                                        # cash residual
  
  # --- Government / Central bank ---
  Bs  <- prev[["Bs"]] + G + r * prev[["Bs"]] -
    Tx - r * prev[["Bcb"]]                       # gov. budget on bills
  Bcb <- Bs - Bh                                      # CB absorbs slack
  Hs  <- Bcb                                          # CB liability = its assets
  
  c(Y = Y, Tx = Tx, YD = YD, C = C, V = V,
    Bh = Bh, H = H, Bs = Bs, Bcb = Bcb, Hs = Hs)
}

params_pc <- list(alpha1 = 0.6, alpha2 = 0.4, theta = 0.2,
                  G = 20, r = 0.025,
                  lam0 = 0.635, lam1 = 5, lam2 = 0.01)

# Clean zero start — no more "1e-6 hack" needed thanks to linearised Tobin
init_pc <- setNames(rep(0, 10),
                    c("Y","Tx","YD","C","V","Bh","H","Bs","Bcb","Hs"))

pc <- simulate_sfc(pc_model, as.list(init_pc), T_periods = 80, params_pc)

# Three consistency checks — all should be ~ 1e-15
cat("\nPC checks (machine precision expected):\n")
cat("  V = H + Bh    :", with(pc, max(abs(V - (H + Bh)))), "\n")
cat("  H = Hs        :", with(pc, max(abs(H - Hs))), "\n")
cat("  Bs = Bh + Bcb :", with(pc, max(abs(Bs - (Bh + Bcb)))), "\n")

# --- Plot --------------------------------------------------------
plot(pc$Bh, type = "l", lwd = 2, col = "steelblue",
     ylab = "Stocks", xlab = "period",
     main = "PC: bills and cash converge to steady-state portfolio")
lines(pc$H, lwd = 2, col = "firebrick")
legend("right", c("Bh (bills)", "H (cash)"),
       col = c("steelblue", "firebrick"), lwd = 2)

# --- Causal DAG (household side only, for readability) -----------
labels_pc <- c("Y","T","YD","C","V","Bh","H","V_lag","Bh_lag","r","G")
M_pc <- matrix(c(
  # Y T YD C V Bh H Vl Bhl r  G
  0,0, 0,1,0, 0,0, 0, 0, 0,1,    # Y  <- C, G
  1,0, 0,0,0, 0,0, 0, 1, 1,0,    # T  <- Y, r*Bh_lag
  1,1, 0,0,0, 0,0, 0, 1, 1,0,    # YD <- Y, T, r*Bh_lag
  0,0, 1,0,0, 0,0, 1, 0, 0,0,    # C  <- YD, V_lag
  0,0, 1,1,0, 0,0, 1, 0, 0,0,    # V  <- YD, C, V_lag
  0,0, 1,0,1, 0,0, 0, 0, 1,0,    # Bh <- V, YD, r
  0,0, 0,0,1, 1,0, 0, 0, 0,0,    # H  <- V, Bh
  0,0, 0,0,0, 0,0, 0, 0, 0,0,
  0,0, 0,0,0, 0,0, 0, 0, 0,0,
  0,0, 0,0,0, 0,0, 0, 0, 0,0,
  0,0, 0,0,0, 0,0, 0, 0, 0,0),
  nrow = 11, byrow = TRUE)
plot_causal(M_pc, labels_pc, "PC: r enters via portfolio AND interest income")

# --- Jacobian (independent states: V and Bh) ---------------------
steadystate_pc <- unlist(tail(pc, 1))
J_pc  <- compute_jacobian(pc_model, steadystate_pc, params_pc,
                          state_vars = c("V", "Bh"))
ev_pc <- check_stability(J_pc, "PC")


# ====================================================================
# SECTION 4:  PCEX — PC with adaptive expectations  (G&L ch. 4.7)
# ====================================================================
# Bills B^h are PLANNED at the start of the period based on EXPECTED
# wealth V^e and EXPECTED income YDe. Realised wealth differs ->
# cash H absorbs the gap.
#
# Pedagogical point: B^h is the intentional (chosen) asset,
# H is the residual buffer. This asymmetry comes from the timing of
# information (planning happens before YD is observed), not from
# anything in the Tobin equations themselves.
#
# State dim.: 3  (V, Bh, YD via YDe = YD_lag)
# --------------------------------------------------------------------

pcex_model <- function(cur, prev, p, t) {
  r <- p$r
  G <- p$G
  
  YDe <- prev[["YD"]]                                  # naive expectation
  C   <- p$alpha1 * YDe + p$alpha2 * prev[["V"]]
  Y   <- C + G
  Tx  <- p$theta * (Y + r * prev[["Bh"]])
  YD  <- Y - Tx + r * prev[["Bh"]]                     # realised
  Ve  <- prev[["V"]] + YDe - C                         # EXPECTED end-of-period wealth
  V   <- prev[["V"]] + YD  - C                         # REALISED wealth
  
  # Tobin block on EXPECTED wealth and EXPECTED income
  # (Bh is planned; H is residual)
  Bh  <- p$lam0 * Ve + p$lam1 * r * Ve - p$lam2 * YDe
  H   <- V - Bh                                        # buffer absorbs V - Ve
  
  Bs  <- prev[["Bs"]] + G + r * prev[["Bs"]] -
    Tx - r * prev[["Bcb"]]
  Bcb <- Bs - Bh
  Hs  <- Bcb
  
  c(Y = Y, Tx = Tx, YDe = YDe, YD = YD, C = C,
    Ve = Ve, V = V, Bh = Bh, H = H, Bs = Bs, Bcb = Bcb, Hs = Hs)
}

init_pcex <- setNames(rep(0, 12),
                      c("Y","Tx","YDe","YD","C","Ve","V","Bh","H","Bs","Bcb","Hs"))

pcex <- simulate_sfc(pcex_model, as.list(init_pcex), 80, params_pc)

cat("\nPCEX checks:\n")
cat("  V = H + Bh :", with(pcex, max(abs(V - (H + Bh)))), "\n")
cat("  H = Hs     :", with(pcex, max(abs(H - Hs))), "\n")

# --- Plot: forecast error and its decay --------------------------
plot(pcex$V - pcex$Ve, type = "h", col = "darkorange", lwd = 2,
     ylab = "V - Ve  (forecast error)", xlab = "period",
     main = "PCEX: forecast errors absorbed by H, then decay")
abline(h = 0, col = "grey")

# --- Jacobian (3 states: V, Bh, YD) ------------------------------
steadystate_pcex <- unlist(tail(pcex, 1))
J_pcex  <- compute_jacobian(pcex_model, steadystate_pcex, params_pc,
                            state_vars = c("V", "Bh", "YD"))
ev_pcex <- check_stability(J_pcex, "PCEX")


# Stitch a baseline path with a post-shock path so we can plot both.
# --------------------------------------------------------------------

stitch_shock <- function(baseline, shocked, shock_t) {
  rbind(baseline[1:shock_t, ], shocked[-1, ])
}

params_pc2 <- params_pc; params_pc2$r <- 0.035
pcex_shock <- simulate_sfc(pcex_model, as.list(pcex[40, ]),
                           T_periods = 40, params_pc2)
pcex_full  <- stitch_shock(pcex, pcex_shock, shock_t = 40)

op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(pcex_full$H,  type = "l", lwd = 2, col = "firebrick",
     xlab = "period", ylab = "H",  main = "Cash (buffer)")
abline(v = 40, lty = 3)
plot(pcex_full$Bh, type = "l", lwd = 2, col = "steelblue",
     xlab = "period", ylab = "Bh", main = "Bills (planned)")
abline(v = 40, lty = 3)
plot(pcex_full$Y,  type = "l", lwd = 2,
     xlab = "period", ylab = "Y",  main = "Income")
abline(v = 40, lty = 3)
plot(pcex_full$V - pcex_full$Ve, type = "l", lwd = 2, col = "darkorange",
     xlab = "period", ylab = "V - Ve", main = "Forecast error")
abline(h = 0, col = "grey"); abline(v = 40, lty = 3)
par(op)