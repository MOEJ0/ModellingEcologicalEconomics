# ====================================================================
# SECTION 0:  GENERIC INFRASTRUCTURE
# ====================================================================

# --------------------------------------------------------------------
# 0.1 Gauss-Seidel solver
# ----------------------w----------------------------------------------
# OUTER loop runs over periods t = 1, 2, ...
# INNER loop solves the SIMULTANEOUS system within each period by
# iterating cur <- model_fn(cur, prev, params, t) until consecutive
# guesses agree to tolerance `tol`. Predetermined values (lagged
# stocks, parameters) live in `prev`; current-period unknowns being
# solved live in `cur`.
# --------------------------------------------------------------------

simulate_sfc <- function(model_fn, init, T_periods, params,
                         tol = 1e-9, max_iter = 500) {
  
  vars <- names(init)
  out  <- matrix(NA_real_, nrow = T_periods + 1, ncol = length(vars),
                 dimnames = list(NULL, vars))
  out[1, ] <- unlist(init)
  
  for (t in 2:(T_periods + 1)) {
    prev <- out[t - 1, ]
    cur  <- prev                                  # initial guess for period t
    
    for (k in seq_len(max_iter)) {
      new <- model_fn(cur, prev, params, t)
      
      # Defensive guard: surface NaN/Inf instead of cryptic if() errors
      if (any(!is.finite(new))) {
        stop(sprintf("Non-finite at t=%d, iter=%d in: %s",
                     t, k,
                     paste(names(new)[!is.finite(new)], collapse = ", ")))
      }
      
      if (max(abs(new - cur)) < tol) { cur <- new; break }
      cur <- new
    }
    out[t, ] <- cur
  }
  as.data.frame(out)
}


# --------------------------------------------------------------------
# 0.2 Numerical Jacobian at steady state
# --------------------------------------------------------------------
# Treats the model as a first-order difference equation
#       x_t = f(x_{t-1}; parameters)
# and computes J = d f / d x evaluated at the steady state x*
# via forward finite differences.
#
# STABILITY of the steady state  <=>  all eigenvalues of J have
# absolute values strictly less than 1 (i.e. lie inside the unit circle).
#
# Arguments
#   Steady State : named vector with steady-state values of ALL variables
#                (use unlist(tail(sim, 1)) after a converged simulation)
#                 i.e. the last value of this convergence must be an equilibrium
#   state_vars : character vector of names that are TRUE dynamic states
#                — independent lagged stocks. Other variables enter
#                only contemporaneously and don't add to state dimension.
# --------------------------------------------------------------------

compute_jacobian <- function(model_fn, steadystate, params, state_vars,
                             eps = 1e-6) {
  
  # Solve ONE period given a hypothetical "previous" state.
  # Used both to recover the steady state numerically and to evaluate
  # the perturbed transitions.
  solve_step <- function(prev) {
    cur <- prev
    for (k in 1:500) {
      new <- model_fn(cur, prev, params, t = 0)
      if (max(abs(new - cur)) < 1e-12) return(new)
      cur <- new
    }
    cur
  }
  
  base <- solve_step(steadystate)[state_vars]
  J <- matrix(0, length(state_vars), length(state_vars),
              dimnames = list(state_vars, state_vars))
  
  for (j in seq_along(state_vars)) {
    p <- steadystate
    p[state_vars[j]] <- p[state_vars[j]] + eps
    J[, j] <- (solve_step(p)[state_vars] - base) / eps
  }
  J
}


# --------------------------------------------------------------------
# 0.3 Stability check helper
# --------------------------------------------------------------------

check_stability <- function(J, label = "") {
  ev <- eigen(J)$values
  cat(sprintf("[%s] eigenvalues: %s | max|lambda|=%.4f -> %s\n",
              label,
              paste(round(ev, 4), collapse = ", "),
              max(Mod(ev)),
              ifelse(max(Mod(ev)) < 1, "STABLE", "UNSTABLE")))
  invisible(ev)
}


# --------------------------------------------------------------------
# 0.4 Causal DAG plotter
# --------------------------------------------------------------------
# M[i, j] = 1  <=>  variable i is computed FROM variable j.
# That is, there is a directed edge j -> i in the causal graph.
# Useful for seeing which variables form the simultaneous loop and
# which ones are predetermined.
# --------------------------------------------------------------------

plot_causal <- function(M, labels, title = "") {
  dg <- graph_from_adjacency_matrix(t(M), mode = "directed")
  V(dg)$name <- labels
  plot(dg, layout = layout_with_fr(dg),
       vertex.size = 30, vertex.color = "lightblue",
       vertex.label.color = "black", vertex.label.cex = 0.9,
       vertex.frame.color = NA, edge.arrow.size = 0.35,
       edge.color = "gray30", main = title, margin = -0.05)
}

