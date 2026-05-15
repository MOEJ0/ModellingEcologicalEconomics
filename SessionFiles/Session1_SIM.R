
# ====================================================================
# SECTION 1:  SIM — Simplest Integrated Model  (G&L ch. 3)
# ====================================================================
# Two sectors: households and government. One asset: cash (H).
# Government pays households with cash; households consume, pay
# taxes, hold the rest.
#
# Key identity:   Delta H_t = G_t - T_t
# Steady state:   Y* = G / theta
# State dim.:     1  (the only lagged stock is H)
# --------------------------------------------------------------------

sim_model <- function(cur, prev, p, t) {
  G  <- p$G
  
  Y  <- cur[["C"]] + G                            # (1) goods market
  Tx <- p$theta * Y                               # (2) tax rule
  YD <- Y - Tx                                    # (3) disposable income
  C  <- p$alpha1 * YD + p$alpha2 * prev[["H"]]    # (4) consumption
  H  <- prev[["H"]] + G - Tx                      # (5) money stock identity
  
  c(Y = Y, Tx = Tx, YD = YD, C = C, H = H)
}

# Parameters and initial state
params_sim <- list(alpha1 = 0.6, alpha2 = 0.4, theta = 0.2, G = 20)
init_sim   <- list(Y = 0, Tx = 0, YD = 0, C = 0, H = 0)

# Run
sim <- simulate_sfc(sim_model, init_sim, T_periods = 60, params_sim)

# Accounting check (machine precision when the model is correct)
cat("\nSIM check  Delta H = G - T :",
    with(sim, max(abs(diff(H) - (params_sim$G - Tx[-1])))), "\n")

# --- Plot --------------------------------------------------------
plot(sim$Y, type = "l", lwd = 2, ylab = "Y, C, H", xlab = "period",
     main = "SIM: convergence to stationary state (Y* = G/theta = 100)")
lines(sim$C, col = "steelblue", lwd = 2)
lines(sim$H, col = "firebrick", lwd = 2)
legend("bottomright", c("Y", "C", "H"),
       col = c("black", "steelblue", "firebrick"), lwd = 2)

# --- Causal DAG --------------------------------------------------
labels_sim <- c("Y", "T", "YD", "C", "H", "G", "theta", "H_lag")
M_sim <- matrix(c(
  # Y T YD C H G th Hl
  0,0, 0,1,0,1,0, 0,    # Y  <- C, G
  1,0, 0,0,0,0,1, 0,    # T  <- Y, theta
  1,1, 0,0,0,0,0, 0,    # YD <- Y, T
  0,0, 1,0,0,0,0, 1,    # C  <- YD, H_lag
  0,1, 0,0,0,1,0, 1,    # H  <- T, G, H_lag
  0,0, 0,0,0,0,0, 0,    # G  exogenous
  0,0, 0,0,0,0,0, 0,    # theta parameter
  0,0, 0,0,0,0,0, 0),   # H_lag predetermined
  nrow = 8, byrow = TRUE)
plot_causal(M_sim, labels_sim, "SIM: causal structure")

# --- Jacobian / stability ----------------------------------------
steadystate_sim <- unlist(tail(sim, 1))
J_sim  <- compute_jacobian(sim_model, steadystate_sim, params_sim, state_vars = "H")
ev_sim <- check_stability(J_sim, "SIM")
# Analytical check: lambda = 1 - theta*alpha2 / (1 - alpha1*(1-theta))
# With alpha1=0.6, alpha2=0.4, theta=0.2  ->  lambda ~= 0.846
