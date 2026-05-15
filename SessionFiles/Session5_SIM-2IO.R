# ====================================================================
# SECTION 5:  SIM-IO — SIM with Leontief production block
# ====================================================================
# Same macro SIM, but Y is decomposed into n sectoral outputs linked
# by a matrix A of technical (input-output) coefficients.
#
#   x_t = A x_t + f_t       =>   x_t = (I - A)^(-1) f_t = L f_t
#   f_t = beta * C_t + gamma * G_t        (sectoral final demand)
#   va  = (1 - colSums(A)) * x            (value added per sector)
#
# The IO block is STATIC within a period — it decodes Y into
# sectors but adds no new state variable. So the Jacobian is the
# same as SIM's: the IO extension changes composition, not dynamics.
# --------------------------------------------------------------------

sim_io_model <- function(cur, prev, p, t) {
  
  # Macro block (identical to SIM)
  G  <- p$G
  Y  <- cur[["C"]] + G
  Tx <- p$theta * Y
  YD <- Y - Tx
  C  <- p$alpha1 * YD + p$alpha2 * prev[["H"]]
  H  <- prev[["H"]] + G - Tx
  
  # IO block: decompose aggregate demand into sectoral outputs
  f  <- p$beta * C + p$gamma * G
  x  <- as.numeric(p$L %*% f)
  va <- (1 - colSums(p$A)) * x
  
  c(Y = Y, Tx = Tx, YD = YD, C = C, H = H,
    setNames(x,  paste0("x",  seq_along(x))),
    setNames(va, paste0("va", seq_along(va))))
}

# Two-sector economy: 1 = goods, 2 = services
A_mat <- matrix(c(0.30, 0.20,
                  0.15, 0.10), nrow = 2, byrow = TRUE)
L_mat <- solve(diag(2) - A_mat)                  # Leontief inverse

params_io <- list(
  alpha1 = 0.6, alpha2 = 0.4, theta = 0.2, G = 20,
  A = A_mat, L = L_mat,
  beta  = c(0.4, 0.6),                           # household consumption shares
  gamma = c(0.7, 0.3)                            # government skewed to goods
)

init_io <- setNames(rep(0, 9),
                    c("Y","Tx","YD","C","H","x1","x2","va1","va2"))

io <- simulate_sfc(sim_io_model, as.list(init_io), 60, params_io)

cat("\nSIM-IO checks:\n")
cat("  Y = va1 + va2 :", with(io, max(abs(Y - (va1 + va2)))), "\n")
cat("  Y = C + G     :", with(io, max(abs(Y - (C + params_io$G)))), "\n")

# --- Plot sectoral outputs ---------------------------------------
matplot(io[, c("x1", "x2")], type = "l", lwd = 2, lty = 1,
        col = c("steelblue", "darkorange"),
        xlab = "period", ylab = "Gross sectoral output",
        main = "SIM-IO: sectoral outputs (n = 2)")
legend("bottomright", c("x1 goods", "x2 services"),
       col = c("steelblue", "darkorange"), lwd = 2)

# --- Jacobian: same as SIM (IO adds no state) --------------------
steadystate_io <- unlist(tail(io, 1))
J_io  <- compute_jacobian(sim_io_model, steadystate_io, params_io,
                          state_vars = "H")
ev_io <- check_stability(J_io, "SIM-IO")