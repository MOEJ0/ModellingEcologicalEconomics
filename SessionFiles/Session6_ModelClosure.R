
# ====================================================================
# SECTION 6a:  SIM-VERTICALIST — closure exercise
# ====================================================================
# SAME equations as SIM, but SWAP exogenous and endogenous variables:
#   Exogenous  : H-bar (central bank sets the money stock)
#   Endogenous : G    (government spending adjusts so that
#                      Delta H = G - T holds with H predetermined)
#
# Steady state:  Y* = alpha2 * H-bar / [ (1 - theta)(1 - alpha1) ]
# With our numbers, H-bar = 80 gives Y* = 100  (same as SIM with G = 20).
#
# PEDAGOGICAL POINT: same equations + same parameters, but a different
# closure rule gives a different CAUSAL story. The DAG arrows literally
# flip; the Jacobian eigenvalue collapses to ~ 0 because H is no longer
# endogenous and the only dynamic state vanishes.
# --------------------------------------------------------------------

sim_vert_model <- function(cur, prev, p, t) {
  H_t <- p$H_target                                    # exogenous
  dH  <- H_t - prev[["H"]]
  Y   <- (cur[["C"]] + dH) / (1 - p$theta)             # from Y(1-theta) = C + dH
  Tx  <- p$theta * Y
  G   <- dH + Tx                                       # endogenous!
  YD  <- Y - Tx
  C   <- p$alpha1 * YD + p$alpha2 * prev[["H"]]
  c(Y = Y, Tx = Tx, YD = YD, C = C, H = H_t, G = G)
}

params_vert <- list(alpha1 = 0.6, alpha2 = 0.4, theta = 0.2, H_target = 80)
# Start H AT the target so we don't see a one-shot jump at t = 1.
init_vert   <- list(Y = 0, Tx = 0, YD = 0, C = 0, H = 80, G = 0)

sim_vert <- simulate_sfc(sim_vert_model, init_vert, 60, params_vert)

# Cross-closure comparison: same Y*, different exogeneity story
cat("\nSteady-state comparison across closures:\n")
cat(sprintf("  SIM (fiscal-led) : Y* = %.2f, G  = %.2f, H  = %.2f\n",
            tail(sim$Y, 1), params_sim$G, tail(sim$H, 1)))
cat(sprintf("  SIM-verticalist  : Y* = %.2f, G* = %.2f, H-bar = %.2f\n",
            tail(sim_vert$Y, 1), tail(sim_vert$G, 1), params_vert$H_target))

steadystate_vert <- unlist(tail(sim_vert, 1))
J_vert  <- compute_jacobian(sim_vert_model, steadystate_vert, params_vert,
                            state_vars = "H")
ev_vert <- check_stability(J_vert, "SIM-vert")




# ====================================================================
# SECTION 6b:  CROSS-MODEL COMPARISON
# ====================================================================

results <- data.frame(
  model      = c("SIM","SIMEX","PC","PCEX","SIM-IO","SIM-vert"),
  state_dim  = c(length(ev_sim),  length(ev_simex), length(ev_pc),
                 length(ev_pcex), length(ev_io),    length(ev_vert)),
  max_abs_ev = c(max(Mod(ev_sim)),   max(Mod(ev_simex)),
                 max(Mod(ev_pc)),    max(Mod(ev_pcex)),
                 max(Mod(ev_io)),    max(Mod(ev_vert))),
  stable     = c(max(Mod(ev_sim))    < 1,
                 max(Mod(ev_simex))  < 1,
                 max(Mod(ev_pc))     < 1,
                 max(Mod(ev_pcex))   < 1,
                 max(Mod(ev_io))     < 1,
                 max(Mod(ev_vert))   < 1)
)
cat("\n----- Cross-model summary -----\n")
print(results, row.names = FALSE)