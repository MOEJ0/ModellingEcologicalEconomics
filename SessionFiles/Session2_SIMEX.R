# ====================================================================
# SECTION 2:  SIMEX — SIM with adaptive expectations  (G&L ch. 3.7)
# ====================================================================
# Same accounting as SIM. Households now decide C_t BEFORE observing
# YD_t. They form an expectation YDe_t = YD_{t-1} (naive adaptive)
# and let cash H absorb the forecast error (buffer mechanism).
#
# Steady state is identical to SIM (in steadystate the expectation is correct).
# Transition path differs because expectations lag reality.
# State dim.: 2  (H, YD)
# --------------------------------------------------------------------

simex_model <- function(cur, prev, p, t) {
  G   <- p$G
  
  YDe <- prev[["YD"]]                              # naive expectation
  C   <- p$alpha1 * YDe + p$alpha2 * prev[["H"]]   # consumption planned
  Y   <- C + G                                     # implies Y from C, G
  Tx  <- p$theta * Y
  YD  <- Y - Tx                                    # realised
  H   <- prev[["H"]] + YD - C                      # H absorbs YD - YDe
  
  c(Y = Y, Tx = Tx, YD = YD, YDe = YDe, C = C, H = H)
}

params_simex <- params_sim                         # same parameters as SIM
init_simex   <- list(Y = 0, Tx = 0, YD = 0, YDe = 0, C = 0, H = 0)

simex <- simulate_sfc(simex_model, init_simex, 60, params_simex)

# Accounting check (household budget identity)
cat("SIMEX check Delta H = YD - C :",
    with(simex, max(abs(diff(H) - (YD - C)[-1]))), "\n")

# --- Plot --------------------------------------------------------
plot(sim$Y, type = "l", lwd = 2, ylim = c(0, 110),
     ylab = "Y", xlab = "period",
     main = "SIM vs SIMEX — same steady state, different path")
lines(simex$Y, col = "darkorange", lwd = 2, lty = 2)
legend("bottomright", c("SIM (perfect foresight)", "SIMEX (adaptive)"),
       col = c("black", "darkorange"), lwd = 2, lty = c(1, 2))

# --- Causal DAG (YDe closes a longer loop) -----------------------
labels_simex <- c("Y","T","YD","C","H","YDe","H_lag","YD_lag")
M_simex <- matrix(c(
  0,0,0,1,0,0,0,0,    # Y    <- C
  1,0,0,0,0,0,0,0,    # T    <- Y
  1,1,0,0,0,0,0,0,    # YD   <- Y, T
  0,0,0,0,0,1,1,0,    # C    <- YDe, H_lag
  0,0,1,1,0,0,1,0,    # H    <- YD, C, H_lag
  0,0,0,0,0,0,0,1,    # YDe  <- YD_lag
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0),
  nrow = 8, byrow = TRUE)
plot_causal(M_simex, labels_simex, "SIMEX: YDe enters the loop")

# --- Jacobian: state is now 2-D (H AND YD via YDe = YD_lag) ------
steadystate_simex <- unlist(tail(simex, 1))
J_simex  <- compute_jacobian(simex_model, steadystate_simex, params_simex,
                             state_vars = c("H", "YD"))
ev_simex <- check_stability(J_simex, "SIMEX")