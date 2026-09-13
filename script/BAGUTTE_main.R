# ============================================================================
# Red-Black Gibbs Sampling for Spatially Correlated Change Point Detection
# ============================================================================

library(doParallel)
library(foreach)

#' Extract edges and uncentered HPD credible intervals as candidate edges
#' @param cp_prob Vector of marginal posterior probabilities along the line
#' @param cp_locations Local indices of the chosen point estimates (MAP/Medoid). 
#'        Defined as the first point of the second segment.
#' @param idx Global spatial indices corresponding to the line
#' @param window_size How far to look around the CP to calculate the local mass
#' @param target_mass The credible level (e.g., 0.95 for 95% CI)
extract_edge_and_ci <- function(cp_prob, cp_locations, idx, window_size = 10, target_mass = 0.95) {
  n_line <- length(cp_prob)
  results <- list()
  
  if (length(cp_locations) == 0) return(results)
  
  for (loc_cp in cp_locations) {
    # Skip if it's the very first or last point. 
    # Must be >= 2 because an edge requires a previous point (loc_cp - 1).
    if (loc_cp <= 1 || loc_cp > n_line) next
    
    # 1. The Point Estimate Edge (Optimal Boundary)
    # [Last point of previous segment, First point of next segment]
    edge_pair <- c(idx[loc_cp - 1], idx[loc_cp])
    
    # 2. The HPD Credible Interval: Finding all high-probability "candidate edges"
    # Left bound must be at least 2 to ensure (k-1) > 0 for edge creation
    left_bound <- max(2, loc_cp - window_size)
    right_bound <- min(n_line, loc_cp + window_size)
    local_indices <- left_bound:right_bound
    
    local_probs <- cp_prob[local_indices]
    total_local_mass <- sum(local_probs)
    
    # Sort by probability (descending) to get the Highest Posterior Density (HPD)
    sort_order <- order(local_probs, decreasing = TRUE)
    sorted_probs <- local_probs[sort_order]
    sorted_indices <- local_indices[sort_order]
    
    # Accumulate mass until target is hit
    cum_mass <- cumsum(sorted_probs) / total_local_mass
    keep_n <- which(cum_mass >= target_mass)[1]
    if (is.na(keep_n)) keep_n <- length(sorted_probs) # Fallback
    
    # The credible set of local indices k (re-sorted spatially for continuity)
    ci_local <- sort(sorted_indices[1:keep_n])
    
    # ==========================================
    # CORE MODIFICATION: CI as a set of Edges
    # For every valid k in the CI, create an edge [k-1, k]
    # ==========================================
    ci_edges_global <- cbind(idx[ci_local - 1], idx[ci_local])
    colnames(ci_edges_global) <- c("upstream_node", "downstream_node")
    
    # Extract the posterior probabilities for these specific candidate edges
    ci_probs <- cp_prob[ci_local]
    
    results[[length(results) + 1]] <- list(
      cp_local = loc_cp,
      edge_global = edge_pair,           # Vector of length 2: The optimal boundary
      ci_edges_global = ci_edges_global, # N x 2 Matrix: All possible boundaries in 95% CI
      ci_probs = ci_probs                # Vector of length N: Probabilities for each boundary
    )
  }
  
  return(results)
}

#' Red-Black Block Gibbs Sampler for Multiple Lines with MRF Prior
#'
#' @param count Expression count matrix
#' @param loc Spatial coordinates
#' @param gene_number Gene indices to analyze
#' @param direction Direction for line extraction
#' @param n_chains Number of MCMC chains
#' @param n_iter Number of iterations per chain
#' @param burnin Burnin period
#' @param n_rb_cycles Number of red-black Gibbs cycles
#' @param cp_distance Minimum distance between change points
#' @param beta_sd_scale Scaling for beta proposal
#' @param lambda Poisson prior rate for number of CPs
#' @param K_int Maximum number of change points (0 = auto)
#' @param mrf_alpha MRF prior strength (encourages spatial correlation)
#' @param mrf_beta MRF prior penalty (discourages isolated CPs)
#' @param mrf_bandwidth Spatial bandwidth for MRF (in index units)
#' @param mrf_gamma Neighbor compounding exponent (gamma>1 boosts support when both neighbors align)
#' @param cluster Gene cluster information
#' @param n_cores Number of cores for parallel processing
#'
test_coefficient_multiple_mrf <- function(
    count, loc,
    gene_number = 1,
    direction   = c("row", "col", "main diag", "minor diag"),
    cpp_file = "BAGUTTE_cpt_MRF_patched_mrf_gamma_warmstart_v3.cpp",
    cp_decision = c("medoid", "map", "ppi"),
    cp_medoid_tol = NULL,
    medoid_loss = c("avg", "hausdorff"),
    medoid_max_samples = 2000,
    n_chains = 4,
    n_iter   = 3000,
    burnin   = 1500,
    n_rb_cycles = 10,  # Number of red-black alternating cycles
    cp_distance = 5,
    beta_sd_scale = 5,
    beta_slope_prior_sd = NULL,  # slope prior SD for beta1; if NULL computed from beta_sd_scale and line length
    lambda = 1e-4,
    K_int = 0,
    mrf_alpha = 2.0,
    mrf_beta = 0.5,
    mrf_bandwidth = 5.0,
    mrf_gamma = 1.0,
    mrf_method = "binary",
    manual_s = NULL,
    cluster,
    n_cores = parallel::detectCores(logical = FALSE))
{
  direction <- match.arg(direction)
  cp_decision <- match.arg(cp_decision)
  medoid_loss <- match.arg(medoid_loss)
  if (is.null(cp_medoid_tol)) cp_medoid_tol <- cp_distance
  
  # Load C++ functions
  Rcpp::sourceCpp(cpp_file)
  
  # ---- Helper functions ------------------------------------------
  rot <- function(th) matrix(c(cos(th), -sin(th), sin(th), cos(th)), 2, 2, byrow = TRUE)
  
  R <- switch(direction,
              "col"        = rot(0),
              "row"        = rot(pi/2),
              "main diag"  = rot(pi/3 + pi/2),
              "minor diag" = rot(-pi/3 + pi/2))
  
  n_spots <- nrow(loc)
  G       <- length(gene_number)
  if (!is.null(manual_s)) {
    if (length(manual_s) == 1) s <- rep(manual_s, n_spots)
    else s <- manual_s
    message("Using user-provided size factors.")
  } else {
    s <- rowSums(count) / exp(mean(log(rowSums(count))))
  }
  
  # ---- Identify lines & Assign Colors (Inheritance Strategy) ----
  loc_rot <- loc %*% t(R)
  
  # 1. First DBSCAN: Identify "Physical Rows" (Primary Rows)
  db <- dbscan::dbscan(matrix(loc_rot[, 1], ncol = 1), eps = 0.1, minPts = 5)
  
  # Get spot indices for all physical rows (remove noise points labeled 0)
  primary_groups <- split(seq_len(n_spots), db$cluster)
  primary_groups <- primary_groups[names(primary_groups) != "0"]
  
  # 2. Spatially sort physical rows and pre-assign colors
  # Calculate the average position of each row
  primary_pos <- sapply(primary_groups, function(idx) mean(loc_rot[idx, 1]))
  
  # Sort from top to bottom (Spatial Sorting)
  ord <- order(primary_pos)
  primary_groups_sorted <- primary_groups[ord]
  
  # Assign colors directly here: Odd rows Red, Even rows Black
  # This ensures that all fragments of the same row will inherit this color later
  primary_colors <- rep(c("red", "black"), length.out = length(primary_groups_sorted))
  
  # 3. Orthogonal Splitting (Orthogonal Clustering) & Color Inheritance
  R_orth <- rot(pi/2) %*% R
  
  sub_slice_index <- list()
  line_colors <- character() # To store the final color of each segment
  group_id <- 1
  
  # Iterate through sorted physical rows
  for (i in seq_along(primary_groups_sorted)) {
    group <- primary_groups_sorted[[i]]
    p_color <- primary_colors[i] # [INHERITANCE] Get the color of the current physical row
    
    if (length(group) < 3) next
    
    loc_sub <- loc[group, , drop = FALSE]
    loc_rot_orth <- loc_sub %*% t(R_orth)
    
    # Second DBSCAN: Split broken lines (handle gaps)
    db2 <- dbscan::dbscan(matrix(loc_rot_orth[, 1], ncol = 1), eps = 1.1, minPts = 3)
    sub_groups <- split(group, db2$cluster)
    
    # Save split sub-segments
    for (sg_name in names(sub_groups)) {
      if (sg_name == "0") next # Skip noise
      sg <- sub_groups[[sg_name]]
      if (length(sg) == 0) next
      
      sub_slice_index[[group_id]] <- sg
      
      # [INHERITANCE] Sub-segments directly inherit the physical row's color
      # This ensures segments 10a and 10b both become Red, as they come from the same loop
      line_colors[group_id] <- p_color 
      
      group_id <- group_id + 1
    }
  }
  
  slice_index <- sub_slice_index
  n_lines <- length(slice_index)
  
  message(sprintf("Found %d lines (segments) in direction '%s'", n_lines, direction))
  
  # ---- Build neighbor structure -----------------------------------
  # (This part remains unchanged; calculate positions to build neighbor relationships)
  line_positions <- matrix(0, nrow = n_lines, ncol = 2)
  for (i in seq_along(slice_index)) {
    idx <- slice_index[[i]]
    line_positions[i, ] <- colMeans(loc_rot[idx, , drop = FALSE])
  }
  
  calculated_dist <- median(diff(sort(unique(primary_pos)))) # Use physical row spacing for robustness
  
  neighbor_graph <- build_neighbor_structure(
    line_ids = seq_len(n_lines),
    line_positions = line_positions,
    neighbor_distance = calculated_dist * 0.9 
  )
  
  # ---- Determine red-black coloring --------------------
  
  red_lines <- which(line_colors == "red")
  black_lines <- which(line_colors == "black")
  
  message(sprintf("Red-Black Assignment (Inherited): %d Red, %d Black segments.", 
                  length(red_lines), length(black_lines)))
  
  # ---- Initialize storage ------------------------------------------
  output_dir <- file.path(getwd(), "fitted_plots_HCC1_MRF", direction)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Storage for CP states across all lines
  n_clusters <- max(cluster$cluster_labels_origin)
  all_line_results <- vector("list", n_lines)
  
  # Initialize CP states for all lines
  cp_states <- vector("list", n_lines)
  for (i in seq_len(n_lines)) {
    cp_states[[i]] <- list()
    for (c in seq_len(n_clusters)) {
      cp_states[[i]][[c]] <- integer(0)  # Empty initial CPs
    }
  }
  
  # ---- Red-Black Gibbs Sampling ------------------------------------
  for (cycle in seq_len(n_rb_cycles)) {
    message(sprintf("\n========== Red-Black Cycle %d/%d ==========", cycle, n_rb_cycles))
    
    # ---- UPDATE RED LINES (using black neighbors) ------------------
    message(sprintf("Updating RED lines (%d lines)...", length(red_lines)))
    
    red_results <- update_line_group(
      line_indices = red_lines,
      slice_index = slice_index,
      neighbor_graph = neighbor_graph,
      cp_states = cp_states,
      count = count,
      loc = loc,
      loc_rot = loc_rot,
      gene_number = gene_number,
      cluster = cluster,
      s = s,
      n_chains = n_chains,
      n_iter = n_iter,
      burnin = burnin,
      cp_distance = cp_distance,
      beta_sd_scale = beta_sd_scale,
      beta_slope_prior_sd = beta_slope_prior_sd,
      lambda = lambda,
      K_int = K_int,
      mrf_alpha = mrf_alpha,
      mrf_beta = mrf_beta,
      mrf_bandwidth = mrf_bandwidth,
      mrf_gamma = mrf_gamma,
      mrf_method = mrf_method,
      n_cores = n_cores,
      direction = direction,
      output_dir = output_dir,
      cpp_file = cpp_file,
      cp_decision = cp_decision,
      cp_medoid_tol = cp_medoid_tol,
      medoid_loss = medoid_loss,
      medoid_max_samples = medoid_max_samples
    )
    
    # Update CP states for red lines
    for (i in seq_along(red_lines)) {
      line_id <- red_lines[i]
      res <- red_results[[i]]
      
      if (!is.null(res)) {
        cp_states[[line_id]] <- res$cp_states
        all_line_results[[line_id]] <- res
      }
    }
    
    # ---- UPDATE BLACK LINES (using red neighbors) ------------------
    message(sprintf("Updating BLACK lines (%d lines)...", length(black_lines)))
    
    black_results <- update_line_group(
      line_indices = black_lines,
      slice_index = slice_index,
      neighbor_graph = neighbor_graph,
      cp_states = cp_states,
      count = count,
      loc = loc,
      loc_rot = loc_rot,
      gene_number = gene_number,
      cluster = cluster,
      s = s,
      n_chains = n_chains,
      n_iter = n_iter,
      burnin = burnin,
      cp_distance = cp_distance,
      beta_sd_scale = beta_sd_scale,
      beta_slope_prior_sd = beta_slope_prior_sd,
      lambda = lambda,
      K_int = K_int,
      mrf_alpha = mrf_alpha,
      mrf_beta = mrf_beta,
      mrf_bandwidth = mrf_bandwidth,
      mrf_gamma = mrf_gamma,
      mrf_method = mrf_method,
      n_cores = n_cores,
      direction = direction,
      output_dir = output_dir,
      cpp_file = cpp_file,
      cp_decision = cp_decision,
      cp_medoid_tol = cp_medoid_tol,
      medoid_loss = medoid_loss,
      medoid_max_samples = medoid_max_samples
    )
    
    # Update CP states for black lines
    for (i in seq_along(black_lines)) {
      line_id <- black_lines[i]
      res <- black_results[[i]]
      
      if (!is.null(res)) {
        cp_states[[line_id]] <- res$cp_states
        all_line_results[[line_id]] <- res
      }
    }
    
    message(sprintf("Completed cycle %d", cycle))
  }
  
  # ---- Aggregate results -------------------------------------------
  message("\nAggregating final results...")
  
  slope_loc <- cbind(
    loc,
    matrix(0, n_spots, G,
           dimnames = list(NULL, paste0("slope_g", seq_len(G))))
  )

  mu_loc <- cbind(       
    loc,
    matrix(0, n_spots, G, dimnames = list(NULL, paste0("mu_g", seq_len(G))))
  )
                        
  cp_indicator_mat <- matrix(0, nrow = nrow(count), ncol = n_clusters)
  
  for (line_id in seq_len(n_lines)) {
    res <- all_line_results[[line_id]]
    if (is.null(res)) next
    
    # Fill slope matrix
    slope_loc[res$idx, 3:(2 + G)] <- res$upd
    mu_loc[res$idx, 3:(2 + G)] <- res$upd_mu
      
    # Fill CP indicator matrix
    for (k in seq_along(res$final_cp_index)) {
      cp_idx_k <- res$final_cp_index[[k]]
      if (length(cp_idx_k) == 0) next
      
      global_idx_k <- res$idx[cp_idx_k]
      cp_indicator_mat[global_idx_k, k] <- 1
    }
  }
  
  list(
    direction   = direction,
    slope_loc   = slope_loc,
    mu_loc      = mu_loc,
    cp_matrix   = cp_indicator_mat,
    line_results = all_line_results,
    cp_states_final = cp_states,
    neighbor_graph = neighbor_graph
  )
}

# ---- Posterior CP summaries: posterior medoid of CP sets -----------------
# A "posterior medoid" is a decision-theoretic point summary: pick a sampled CP
# configuration that minimizes posterior expected loss under a set-distance.
# Here we use a symmetric *truncated screening distance* (avg or Hausdorff-like).

cp_make_key <- function(idx) {
  idx <- as.integer(idx)
  if (length(idx) == 0) "" else paste(idx, collapse = ",")
}

cp_parse_key <- function(key) {
  if (is.null(key) || length(key) == 0 || is.na(key) || key == "") return(integer(0))
  as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
}

cp_screening_avg_trunc <- function(A, B, w) {
  A <- as.integer(A); B <- as.integer(B)
  if (length(A) == 0) return(0.0)
  if (length(B) == 0) return(as.double(w))
  d <- vapply(A, function(a) min(abs(a - B)), 0.0)
  mean(pmin(d, w))
}

cp_screening_max_trunc <- function(A, B, w) {
  A <- as.integer(A); B <- as.integer(B)
  if (length(A) == 0) return(0.0)
  if (length(B) == 0) return(as.double(w))
  d <- vapply(A, function(a) min(abs(a - B)), 0.0)
  max(pmin(d, w))
}

cp_set_distance <- function(A, B, w, loss = c("avg", "hausdorff")) {
  loss <- match.arg(loss)
  if (loss == "avg") {
    max(cp_screening_avg_trunc(A, B, w), cp_screening_avg_trunc(B, A, w))
  } else {
    max(cp_screening_max_trunc(A, B, w), cp_screening_max_trunc(B, A, w))
  }
}

posterior_medoid_cp <- function(cp_indicator, w, loss = c("avg", "hausdorff"),
                                max_samples = 2000, seed = 1L) {
  loss <- match.arg(loss)
  if (is.null(dim(cp_indicator))) cp_indicator <- matrix(cp_indicator, ncol = 1)

  S <- ncol(cp_indicator)
  if (S > max_samples) {
    # deterministic thinning to keep it stable across runs
    keep <- unique(round(seq(1, S, length.out = max_samples)))
    cp_indicator <- cp_indicator[, keep, drop = FALSE]
    S <- ncol(cp_indicator)
  }

  keys <- character(S)
  for (s in seq_len(S)) {
    idx <- which(cp_indicator[, s] > 0.5)
    keys[s] <- cp_make_key(idx)
  }

  tab <- sort(table(keys), decreasing = TRUE)
  uniq_keys <- names(tab)
  weights <- as.numeric(tab)
  uniq_sets <- lapply(uniq_keys, cp_parse_key)

  U <- length(uniq_sets)
  total_w <- sum(weights)
  if (U == 1) {
    return(list(
      cps = uniq_sets[[1]],
      key = uniq_keys[[1]],
      post_mass = 1.0,
      risk = 0.0,
      unique = tab,
      loss = loss,
      w = w
    ))
  }

  # weighted posterior expected loss for each unique configuration
  risk <- numeric(U)
  for (i in seq_len(U)) {
    Ai <- uniq_sets[[i]]
    ssum <- 0.0
    for (j in seq_len(U)) {
      if (j == i) next
      d <- cp_set_distance(Ai, uniq_sets[[j]], w = w, loss = loss)
      ssum <- ssum + weights[j] * d
    }
    risk[i] <- ssum / total_w
  }

  i_star <- which.min(risk)
  list(
    cps = uniq_sets[[i_star]],
    key = uniq_keys[[i_star]],
    post_mass = weights[i_star] / total_w,
    risk = risk[i_star],
    unique = tab,
    loss = loss,
    w = w
  )
}

# ---- Posterior CP summaries: configuration-level MAP -----------------------
# MAP here means the *most probable CP configuration* under the posterior, i.e.
# the sampled CP set (after burn-in) with the highest posterior probability.
# We approximate this by the empirical mode among sampled configurations.

posterior_map_cp <- function(cp_indicator, max_samples = 2000) {
  if (is.null(dim(cp_indicator))) cp_indicator <- matrix(cp_indicator, ncol = 1)

  S <- ncol(cp_indicator)
  if (S > max_samples) {
    keep <- unique(round(seq(1, S, length.out = max_samples)))
    cp_indicator <- cp_indicator[, keep, drop = FALSE]
    S <- ncol(cp_indicator)
  }

  keys <- character(S)
  for (s in seq_len(S)) {
    idx <- which(cp_indicator[, s] > 0.5)
    keys[s] <- cp_make_key(idx)
  }

  tab <- sort(table(keys), decreasing = TRUE)
  key_star <- names(tab)[1]
  cps_star <- cp_parse_key(key_star)

  list(
    cps = cps_star,
    key = key_star,
    post_mass = as.numeric(tab[1]) / sum(tab),
    unique = tab
  )
}

#' Update a group of lines in parallel (helper for red-black Gibbs)
update_line_group <- function(
    line_indices,
    slice_index,
    neighbor_graph,
    cp_states,
    count, loc, loc_rot,
    gene_number,
    cluster,
    s,
    n_chains,
    n_iter,
    burnin,
    cp_distance,
    beta_sd_scale,
    beta_slope_prior_sd,
    lambda,
    K_int,
    mrf_alpha,
    mrf_beta,
    mrf_bandwidth,
    mrf_gamma,
    mrf_method,
    cp_decision,
    cp_medoid_tol,
    medoid_loss,
    medoid_max_samples,
    n_cores,
    direction,
    output_dir,
    cpp_file)
{
  # Setup parallel cluster
  cl <- makeCluster(min(n_cores, length(line_indices)))
  on.exit(stopCluster(cl), add = TRUE)
  registerDoParallel(cl)

  cpp_file_abs <- normalizePath(cpp_file, mustWork = TRUE)
  
  # Export functions and data
  clusterExport(cl, varlist = c("cpp_file_abs"), envir = environment())
  clusterEvalQ(cl, {
    library(Rcpp)
    Rcpp::sourceCpp(cpp_file_abs)
  })
  
  clusterExport(cl,
                varlist = c("count", "loc", "loc_rot", "gene_number", "s",
                            "cp_distance", "beta_sd_scale", "beta_slope_prior_sd", "lambda", 
                            "K_int", "mrf_alpha", "mrf_beta", "mrf_bandwidth", "mrf_gamma", "mrf_method",
                            "cp_decision", "cp_medoid_tol", "medoid_loss", "medoid_max_samples",
                            "cp_make_key", "cp_parse_key", "cp_screening_avg_trunc", "cp_screening_max_trunc",
                            "cp_set_distance", "posterior_medoid_cp", "posterior_map_cp",
                            "extract_edge_and_ci", # Ensure the new edge-based CI helper is exported to parallel workers
                            "cluster", "slice_index", "neighbor_graph", "cp_states",
                            "n_chains", "n_iter", "burnin", "direction", "output_dir",
                           "compute_posterior_summaries_multi_gene_cube"),
                envir = environment())
  
  # Process lines in parallel
  results <- foreach(line_id = line_indices,
                     .packages = c("Rcpp", "RcppArmadillo")) %dopar% {
    
    idx <- slice_index[[line_id]]
    
    if (length(idx) < 10) {
      return(NULL)
    }
    
    # Get neighbor CP states
    neighbors <- neighbor_graph[[line_id]]
    neighbor_cp_list <- list()
    
    for (neighbor_type in c("upper_left", "upper_right", "lower_left", "lower_right")) {
      neighbor_ids <- neighbors[[neighbor_type]]
      
      if (length(neighbor_ids) > 0) {
        for (nb_id in neighbor_ids) {
          if (nb_id <= length(cp_states) && !is.null(cp_states[[nb_id]])) {
            # Collect CPs for this neighbor line as a per-cluster list (to avoid cross-cluster crosstalk)
            neighbor_cp_list <- c(neighbor_cp_list, list(cp_states[[nb_id]]))
          }
        }
      }
    }
    
    # Order spots along the line
    ord <- loc_rot[idx, 2]
    ord_rank <- order(ord)
    ord <- ord[ord_rank]
    idx <- idx[ord_rank]
    
    # Extract data
    pattern_all <- count[idx, gene_number, drop = FALSE]
    G_all <- length(gene_number)
    s_sub <- s[idx]
    
    keep_genes <- cluster$cluster_labels_origin != 0
    if (any(!keep_genes)) {
      active_idx <- which(keep_genes)
    } else {
      active_idx <- seq_len(G_all)
    }
    
    gene_active <- gene_number[active_idx]
    pattern_active <- pattern_all[, active_idx, drop = FALSE]
    cluster_active <- cluster$cluster_labels
    G_active <- length(gene_active)
    
    # Initialize CPs (warm start from previous outer iteration / previous RB cycle)
    C_max <- max(cluster_active)
    cp_init <- vector("list", C_max)

    prev_state_list <- NULL
    if (line_id <= length(cp_states)) prev_state_list <- cp_states[[line_id]]

    n_line <- length(ord)

    if (is.list(prev_state_list) && length(prev_state_list) > 0) {
      for (cc in seq_len(C_max)) {
        v <- integer(0)
        if (cc <= length(prev_state_list) && !is.null(prev_state_list[[cc]])) {
          v <- as.integer(prev_state_list[[cc]])
        }
        v <- sort(unique(v))
        # Keep within valid range (exclude endpoints)
        v <- v[v > 1 & v < n_line]
        # Thin to satisfy cp_distance (greedy)
        if (length(v) > 1 && cp_distance > 1) {
          vv <- integer(0)
          for (t in v) {
            if (length(vv) == 0 || (t - tail(vv, 1)) >= cp_distance) {
              vv <- c(vv, t)
            }
          }
          v <- vv
        }
        cp_init[[cc]] <- v
      }
    } else {
      for (cc in seq_len(C_max)) cp_init[[cc]] <- integer(0)
    }
    
    # ---- beta slope prior SD (controls prior on beta1 in Laplace marginal) ----
    if (is.null(beta_slope_prior_sd) || !is.finite(beta_slope_prior_sd)) {
      x_range <- max(ord) - min(ord)
      if (x_range <= 0) x_range <- 1
      beta_slope_prior_sd <- beta_sd_scale / (2 * x_range)
    }
    beta_slope_prior_sd <- as.numeric(beta_slope_prior_sd)
    chains <- vector("list", n_chains)
    
    # Run RJMCMC with MRF prior
    for (chain_id in seq_len(n_chains)) {
      chains[[chain_id]] <- rjmcmc_single_line_with_neighbors(
        x = ord,
        y = pattern_active,
        gene_clusters = cluster_active,
        s_ = s_sub,
        cp_init = cp_init,
        neighbor_cp_states = neighbor_cp_list,
        n_iter = n_iter,
        burnin = burnin,
        cp_distance = cp_distance,
        beta_proposal_sd = c(1, beta_slope_prior_sd),
        lambda = lambda,
        K_int = K_int,
        tau_logphi = 1,
        mrf_alpha = mrf_alpha,
        mrf_beta = mrf_beta,
        mrf_bandwidth = mrf_bandwidth,
        mrf_gamma = mrf_gamma,
        mrf_method = mrf_method
      )
    }
    
    # Merge chains and summarize
    cp_all <- abind::abind(lapply(chains, `[[`, "cp_indicator"), along = 3)
    beta_all <- do.call(abind::abind, c(lapply(chains, `[[`, "beta"), along = 3))
    intercept_all <- do.call(abind::abind, c(lapply(chains, `[[`, "intercept"), along = 3))
    mu_all <- do.call(abind::abind, c(lapply(chains, `[[`, "mu_loc"), along = 3))
      
    # Summarize change points per cluster
    n_clusters <- dim(cp_all)[2]
    final_cp_x_list <- vector("list", n_clusters)
    final_cp_idx <- vector("list", n_clusters)
    
    for (c in seq_len(n_clusters)) {
      mat_c <- cp_all[, c, ]
      if (is.matrix(mat_c)) {
        # Posterior probability of CP at each location (marginal PPI)
        cp_prob <- rowMeans(mat_c)
        
        if (cp_decision == "medoid") {
          w <- as.integer(round(cp_medoid_tol))
          w <- max(1L, w)
          med <- posterior_medoid_cp(
            cp_indicator = mat_c,
            w = w,
            loss = medoid_loss,
            max_samples = medoid_max_samples
          )
          cp_locations <- med$cps
          
          final_cp_idx[[c]] <- cp_locations
          final_cp_x_list[[c]] <- list(
            idx = cp_locations,
            x = ord[cp_locations],
            medoid_key = med$key,
            medoid_post_mass = med$post_mass,
            medoid_risk = med$risk
          )
        } else if (cp_decision == "map") {
          mp <- posterior_map_cp(
            cp_indicator = mat_c,
            max_samples = medoid_max_samples
          )
          cp_locations <- mp$cps

          final_cp_idx[[c]] <- cp_locations
          final_cp_x_list[[c]] <- list(
            idx = cp_locations,
            x = ord[cp_locations],
            map_key = mp$key,
            map_post_mass = mp$post_mass
          )
        } else {
          # Simple marginal thresholding (legacy)
          cp_locations <- which(cp_prob > 0.5)
          final_cp_idx[[c]] <- cp_locations
          final_cp_x_list[[c]] <- list(idx = cp_locations, x = ord[cp_locations])
        }
        
        # ---- Add candidate edges and spatial HPD credible intervals ----
        cp_details <- extract_edge_and_ci(
          cp_prob = cp_prob,
          cp_locations = cp_locations,
          idx = idx, # Pass global indices array for this line
          window_size = cp_distance * 2, # Utilize cp_distance to scale the search window naturally
          target_mass = 0.95
        )
        final_cp_x_list[[c]]$details <- cp_details
      }
    }
    
    # Compute posterior summaries
    summary <- compute_posterior_summaries_multi_gene_cube(beta_all)
    intercept_sum <- compute_posterior_summaries_multi_gene_cube(intercept_all)
    mu_sum <- compute_posterior_summaries_multi_gene_cube(mu_all)
      
    # Create update matrix
    upd <- matrix(0, nrow = length(idx), ncol = G_active)
    upd_mu <- matrix(0, nrow = length(idx), ncol = G_active)
      
    for (g in seq_len(G_active)) {
      cl <- cluster_active[g]
      cp_x <- final_cp_x_list[[cl]]$x
      
      if (length(cp_x) > 0) {
        seg_g <- cut(ord, breaks = c(-Inf, cp_x, Inf), labels = FALSE)
        seg_groups <- split(seq_along(seg_g), seg_g)
      } else {
        seg_groups <- list(`1` = seq_along(ord))
      }
      
      slope <- summary$mean[, g]
      lower <- summary$lower[, g]
      upper <- summary$upper[, g]
      mu_val <- mu_sum$mean[, g]
        
      for (sel in seg_groups) {
        if (length(sel) > 3) {
          if (mean(lower[sel] > 0) > 0.95 || mean(upper[sel] < 0) > 0.95) {
            upd[sel, g] <- slope[sel]
          }
        }
          upd_mu[sel, g] <- mu_val[sel]
      }
    }
    
    # Map back to full gene matrix
    upd_full <- matrix(0, nrow = nrow(upd), ncol = G_all)
    upd_full[, active_idx] <- upd

    upd_mu_full <- matrix(0, nrow = nrow(upd_mu), ncol = G_all) 
    upd_mu_full[, active_idx] <- upd_mu                         
      
    # Extract CP states for this line
    line_cp_states <- vector("list", n_clusters)
    for (c in seq_len(n_clusters)) {
      line_cp_states[[c]] <- final_cp_idx[[c]]
    }
    
    list(
      idx = idx,
      upd = upd_full,
      upd_mu = upd_mu_full,
      final_cp_index = final_cp_idx,
      final_cp_details = final_cp_x_list, # Include edges and CI info here!
      ord = ord,
      pattern = pattern_all,
      cp_states = line_cp_states
    )
  }
  
  return(results)
}

#' Helper function to compute posterior summaries from 3D array
compute_posterior_summaries_multi_gene_cube <- function(beta_cube) {
  # beta_cube: [G, n, S] where G=genes, n=spots, S=samples
  G <- dim(beta_cube)[1]
  n <- dim(beta_cube)[2]
  
  mean_mat <- matrix(0, n, G)
  lower_mat <- matrix(0, n, G)
  upper_mat <- matrix(0, n, G)
  
  for (g in seq_len(G)) {
    for (i in seq_len(n)) {
      samples <- beta_cube[g, i, ]
      mean_mat[i, g] <- mean(samples)
      lower_mat[i, g] <- quantile(samples, 0.025)
      upper_mat[i, g] <- quantile(samples, 0.975)
    }
  }
  
  list(
    mean = mean_mat,
    lower = lower_mat,
    upper = upper_mat
  )
}

# Helper function to correct Lattice locations
find_optimal_tilt <- function(coords, range_deg = seq(-5, 5, by = 0.1)) {
  scores <- c()
  range_rad <- range_deg * (pi / 180)
  mat <- as.matrix(coords)
  
  for (theta in range_rad) {
    R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2, byrow = TRUE)
    rot_coords <- mat %*% t(R)
    y_vals <- rot_coords[, 2]
    d <- density(y_vals, adjust = 0.5)
    scores <- c(scores, sum(d$y^2)) 
  }
  best_idx <- which.max(scores)
  return(list(angle = range_deg[best_idx], scores = scores))
}