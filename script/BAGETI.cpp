#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ============================================================================
// PART 1: Core helper functions (from original)
// ============================================================================

inline double prob_add(int K, int Kmax) {
  if (K == 0) return 1.0;
  if (K < Kmax) return 0.3;
  return 0.0;
}

inline double prob_remove(int K, int Kmax) {
  if (K == 0) return 0.0;
  if (K < Kmax) return 0.3;
  return 0.6;
}

inline double prob_move(int K, int /*Kmax*/) {
  if (K == 0) return 0.0;
  return 0.4;
}

inline double slog(double v) {
  return (v <= 0.0) ? -INFINITY : std::log(v);
}

double log_sum_exp(double a, double b) {
  return (a > b) ? a + log1p(exp(b - a)) : b + log1p(exp(a - b));
}

double rtruncnorm(double mu, double sigma, double a, double b) {
  double u = R::runif(0, 1);
  double Fa = R::pnorm(a, mu, sigma, 1, 0);
  double Fb = R::pnorm(b, mu, sigma, 1, 0);
  return R::qnorm(Fa + u*(Fb - Fa), mu, sigma, 1, 0);
}

double dtruncnorm(double x, double mean, double sd, double a, double b, bool logd = false) {
  double pdf = R::dnorm(x, mean, sd, false);
  double cdf_a = R::pnorm(a, mean, sd, true, false);
  double cdf_b = R::pnorm(b, mean, sd, true, false);
  double trunc_const = cdf_b - cdf_a;
  double result = pdf / trunc_const;
  return logd ? log(result) : result;
}

inline double dtruncnorm_log(double x, double mean, double sd, double a, double b) {
  return dtruncnorm(x, mean, sd, a, b, true);
}

double dtpois(int k, int K_max, double lambda) {
  if (k < 0 || k > K_max) {
    return -std::numeric_limits<double>::infinity();
  }
  double log_poisson = k * std::log(lambda) - lambda - std::lgamma(k + 1);
  double log_normalization = R::ppois(K_max, lambda, true, true);
  return log_poisson - log_normalization;
}

inline double clamp(double v, double lo, double hi) {
  return std::max(lo, std::min(hi, v));
}

inline double safe_log(double v) {
  return std::log(std::max(v, 1e-12));
}

// ============================================================================
// PART 2: Laplace approximation functions (from original)
// ============================================================================

// [[Rcpp::export]]
Rcpp::List laplace_marginal_NB(const arma::vec& y,
                               const arma::vec& x,
                               const arma::vec& s,
                               double            phi,
                               const arma::vec& mu0,
                               const arma::mat& Sigma0)
{
  const double LOG2PI = std::log(2.0 * M_PI);
  const int d = 2;

  if (y.n_elem == 0) {
    return Rcpp::List::create(
      Rcpp::_["log_marginal"] = 0.0,
      Rcpp::_["beta"] = mu0
    );
  }

  arma::mat Prec = arma::inv_sympd(Sigma0);
  arma::vec beta = mu0;
  arma::mat H(2,2);
  arma::vec g(2);

  for (int it = 0; it < 60; ++it) {
    H.zeros(); g.zeros();

    for (arma::uword i = 0; i < y.n_elem; ++i) {
      double eta = beta[0] + beta[1]*x[i] + safe_log(s[i]);
      eta = clamp(eta, -50.0, 50.0);
      double mu = std::exp(eta);
      double denom = phi + mu;

      double g_eta = phi * (y[i] - mu) / denom;
      double w = (y[i] + phi) * phi * mu / (denom*denom);

      g[0] += g_eta;
      g[1] += g_eta * x[i];
      H(0,0) += w;
      H(0,1) += w * x[i];
      H(1,1) += w * x[i] * x[i];
    }
    H(1,0) = H(0,1);

    arma::vec diff = beta - mu0;
    g -= Prec * diff;
    H += Prec;

    arma::vec step;
    bool ok = arma::solve(step, H, g, arma::solve_opts::fast);
    if (!ok || !step.is_finite()) step = 1e-3 * g;

    double nrm = arma::norm(step, 2);
    if (nrm > 2.0) step *= (2.0 / nrm);

    beta += step;
    if (arma::norm(step, 2) < 1e-6) break;
  }

  double loglik = 0.0;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    double eta = beta[0] + beta[1]*x[i] + safe_log(s[i]);
    eta = clamp(eta, -50.0, 50.0);
    double mu = std::exp(eta);
    loglik += R::dnbinom_mu(y[i], phi, mu, true);
  }

  arma::vec diff = beta - mu0;
  double logprior = -0.5 * arma::as_scalar(diff.t() * Prec * diff)
                  - 0.5 * d * LOG2PI
                  - 0.5 * std::log(std::max(arma::det(Sigma0), 1e-300));

  double val=0.0, sign=1.0;
  bool det_ok = arma::log_det(val, sign, H);
  double logdetH = (det_ok && sign > 0) ? val : std::log(std::max(arma::det(H + arma::eye(2,2)*1e-8), 1e-300));

  double logmarg = loglik + logprior + 0.5*d*LOG2PI - 0.5*logdetH;

  return Rcpp::List::create(
    Rcpp::_["log_marginal"] = logmarg,
    Rcpp::_["beta"] = beta
  );
}

inline Rcpp::List
laplace_marginal_NB_mask(const arma::vec& y,
                         const arma::vec& x,
                         const arma::vec& s,
                         const arma::uvec& mask,
                         double           phi,
                         const arma::vec& mu0,
                         const arma::mat& Sigma0)
{
  arma::uvec keep = find(mask == 0);
  if (keep.is_empty())
    return laplace_marginal_NB(arma::vec(), arma::vec(), arma::vec(), phi, mu0, Sigma0);

  return laplace_marginal_NB(y.elem(keep), x.elem(keep), s.elem(keep), phi, mu0, Sigma0);
}

// ============================================================================
// PART 3: MRF PRIOR - Spatial correlation for change points
// ============================================================================

// Structure to store neighbor information for each line
struct LineNeighbors {
  int line_id;
  std::vector<int> upper_left;   // Line IDs of upper-left neighbors
  std::vector<int> upper_right;  // Line IDs of upper-right neighbors
  std::vector<int> lower_left;   // Line IDs of lower-left neighbors
  std::vector<int> lower_right;  // Line IDs of lower-right neighbors
};

// [[Rcpp::export]]
List build_neighbor_structure(const IntegerVector& line_ids,
                               const NumericMatrix& line_positions,
                               double neighbor_distance = 0.1) {
  int n_lines = line_ids.size();
  List neighbor_list(n_lines);
  
  for (int i = 0; i < n_lines; ++i) {
    IntegerVector ul, ur, ll, lr;
    
    double pos_i = line_positions(i, 0);  // position along perpendicular direction
    
    for (int j = 0; j < n_lines; ++j) {
      if (i == j) continue;
      
      double pos_j = line_positions(j, 0);
      double dist = std::abs(pos_i - pos_j);
      
        if (dist > (neighbor_distance * 0.5) && dist < (neighbor_distance * 1.5)) {        // Determine direction based on second coordinate or ordering
        double offset = line_positions(j, 1) - line_positions(i, 1);
        
        if (pos_j < pos_i) {  // Upper neighbors
          if (offset < 0) ul.push_back(line_ids[j]);
          else ur.push_back(line_ids[j]);
        } else {  // Lower neighbors
          if (offset < 0) ll.push_back(line_ids[j]);
          else lr.push_back(line_ids[j]);
        }
      }
    }
    
    neighbor_list[i] = List::create(
      Named("line_id") = line_ids[i],
      Named("upper_left") = ul,
      Named("upper_right") = ur,
      Named("lower_left") = ll,
      Named("lower_right") = lr
    );
  }
  
  return neighbor_list;
}

// Compute MRF prior contribution for change points

// Compute MRF prior contribution for change points (cluster-specific neighborhood list)
// Robust "soft-OR" matching: for each CP, reward if ANY neighbor line has a nearby CP.
// Uses Gaussian-like kernel (exp(-d^2/(2h^2))) and aggregates across neighbor lines as 1 - prod(1-s_l).
// - alpha: strength of attraction to neighbor-aligned CPs
// - beta : penalty for isolated CPs (including "vacuum" case where neighbors have no CPs)
// - spatial_bandwidth: h (in index units)
// [[Rcpp::export]]
double compute_mrf_prior_log(const std::vector<int>& cp_current,
                             const List& neighbor_cps,
                             double alpha = 2.0,
                             double beta  = 0.5,
                             double spatial_bandwidth = 5.0,
                             double gamma = 1.0,
                             std::string method = "binary") {
  
  if (cp_current.empty()) return 0.0;

  // h: hard radius in binary mode, sigma in gaussian mode
  const double h = std::max(1e-8, spatial_bandwidth);
  const double gam = std::max(1e-8, gamma);
  
  double log_prior = 0.0;
  int n_neighbors = neighbor_cps.size();

  // Pre-process neighbor data
  std::vector< std::vector<int> > nb_vecs;
  nb_vecs.reserve(n_neighbors);
  int total_nb_cps = 0;
  for (int nb = 0; nb < n_neighbors; ++nb) {
    std::vector<int> v = as<std::vector<int>>(neighbor_cps[nb]);
    total_nb_cps += (int)v.size();
    nb_vecs.push_back(std::move(v));
  }

  // Vacuum Penalty
  if (total_nb_cps == 0) {
    log_prior -= beta * (double)cp_current.size();
    return log_prior;
  }

  for (int cp : cp_current) {
    double match_score = 0.0;

    if (method == "binary") {
      // === BINARY MODE (default) ===
      // Logic: If at least one neighbor is within range, it is considered fully supported (Score = 1.0)
      // "One or multiple neighbors have the same effect, unaffected by distance (as long as within range)"
      bool is_supported = false;
      for (int nb = 0; nb < n_neighbors; ++nb) {
        const auto& v = nb_vecs[nb];
        for (int nb_cp : v) {
          double d = std::abs((double)cp - (double)nb_cp);
          if (d <= h) { // Hard threshold decision
            is_supported = true;
            break; 
          }
        }
        if (is_supported) break; // If one support is found, no need to check other neighbors
      }
      match_score = is_supported ? 1.0 : 0.0;

    } else {
      // === GAUSSIAN MODE (original) ===
      // Logic: Gaussian decay + Gamma superposition
      double prod_no_match = 1.0; 
      for (int nb = 0; nb < n_neighbors; ++nb) {
        const auto& v = nb_vecs[nb];
        double best = 0.0;
        for (int nb_cp : v) {
          double d = std::abs((double)cp - (double)nb_cp);
          double s = std::exp(-(d * d) / (2.0 * h * h));
          if (s > best) best = s;
          if (best > 0.999) break; 
        }
        double one_minus = 1.0 - std::min(1.0, std::max(0.0, best));
        prod_no_match *= std::pow(one_minus, gam);
      }
      match_score = 1.0 - prod_no_match;
    }

    // Unified scoring formula:
    // Supported (Score=1) -> +alpha
    // Isolated  (Score=0) -> -beta
    log_prior += alpha * match_score - beta * (1.0 - match_score);
  }

  return log_prior;
}


// Extract cluster-specific neighbor CPs.
// Supports two input formats:
//  (1) neighbor_cp_states is a List of integer vectors (already filtered).
//  (2) neighbor_cp_states is a List over neighbor lines, each element is a List over clusters,
//      where element[[c-1]] is the CP vector for cluster c (1-indexed in our loops).
static inline List neighbor_cps_for_cluster(const List& neighbor_cp_states, int c_one_indexed) {
  int nnb = neighbor_cp_states.size();
  if (nnb == 0) return List();

  List out(nnb);
  int c0 = c_one_indexed - 1;

  for (int nb = 0; nb < nnb; ++nb) {
    SEXP el = neighbor_cp_states[nb];

    if (TYPEOF(el) == VECSXP) {
      List nb_all(el);
      if (c0 >= 0 && c0 < nb_all.size()) {
        out[nb] = nb_all[c0];
      } else {
        out[nb] = IntegerVector(0);
      }
    } else {
      // assume integer vector (flat format)
      out[nb] = el;
    }
  }
  return out;
}



// More sophisticated MRF prior using Ising-like model
// [[Rcpp::export]]
double compute_mrf_prior_ising(const IntegerVector& cp_indicator,
                               const List& neighbor_indicators,
                               double theta = 1.0) {
  int n = cp_indicator.size();
  double log_prior = 0.0;
  
  for (int i = 0; i < n; ++i) {
    if (cp_indicator[i] == 0) continue;
    
    // Check neighbors
    int sum_neighbors = 0;
    int n_neighbors = neighbor_indicators.size();
    
    for (int nb = 0; nb < n_neighbors; ++nb) {
      IntegerVector nb_ind = as<IntegerVector>(neighbor_indicators[nb]);
      if (i < nb_ind.size() && nb_ind[i] > 0) {
        sum_neighbors += nb_ind[i];
      }
    }
    
    // Ising-like: positive theta encourages agreement
    log_prior += theta * sum_neighbors;
  }
  
  return log_prior;
}

// ============================================================================
// PART 4: ZINB likelihood and cluster functions (from original, condensed)
// ============================================================================

inline double zinb_loglik_single_condN(
    const arma::vec& y,
    const arma::uvec& n_mat,
    const arma::vec& x,
    const std::vector<int>& cp,
    const std::vector<arma::vec>& beta_seg,
    double phi,
    const arma::vec& s)
{
  int n = (int)y.n_elem;
  std::vector<int> brk = cp;
  std::sort(brk.begin(), brk.end());
  brk.insert(brk.begin(), 0);
  brk.push_back(n);

  double ll = 0.0;
  for (size_t k = 0; k < brk.size() - 1; ++k) {
    int a = brk[k], b = brk[k+1];
    const arma::vec& beta = beta_seg[k];

    for (int i = a; i < b; ++i) {
      if (n_mat[i] == 1) continue;
      double eta = beta[0] + beta[1]*x[i];
      eta = clamp(eta, -50.0, 50.0);
      double mu = std::exp(eta) * s[i];
      mu = std::max(mu, 1e-12);
      ll += R::dnbinom_mu(y[i], phi, mu, true);
    }
  }
  return ll;
}

inline void clean_cp(std::vector<int>& cp) {
  std::sort(cp.begin(), cp.end());
  cp.erase(std::unique(cp.begin(), cp.end()), cp.end());
}


inline int candidate_count(int n, int cp_distance, const std::vector<int>& cp) {
  int A = 0;
  for (int pos = cp_distance; pos < n - cp_distance; ++pos) {
    bool ok = true;
    for (int cpi : cp) {
      if (std::abs(pos - cpi) < cp_distance) { ok = false; break; }
    }
    if (ok) ++A;
  }
  return A;
}

// For a given sorted CP vector and an index of CP to move, return the number of admissible
// positions in its local window (excluding the current CP position).
inline int move_window_count(const std::vector<int>& cp_sorted, int idx_cp, int n, int cp_distance) {
  int K = (int)cp_sorted.size();
  if (K == 0 || idx_cp < 0 || idx_cp >= K) return 0;
  int left_bd  = (idx_cp == 0 ? 0 : cp_sorted[idx_cp - 1]) + cp_distance;
  int right_bd = (idx_cp == K - 1 ? n : cp_sorted[idx_cp + 1]) - cp_distance;
  int curr = cp_sorted[idx_cp];
  int W = 0;
  for (int pos = left_bd; pos < right_bd; ++pos) {
    if (pos == curr) continue;
    ++W;
  }
  return W;
}

// ============================================================================
// PART 5: RED-BLACK GIBBS SAMPLER (single iteration for one color)
// ============================================================================

// [[Rcpp::export]]
List rjmcmc_single_line_with_neighbors(
    const arma::vec& x,
    const arma::mat& y,
    const IntegerVector& gene_clusters,
    const arma::vec& s_,
    SEXP cp_init,
    const List& neighbor_cp_states,  // Current CP states from neighbor lines
    int n_iter,
    int burnin = 1000,
    int cp_distance = 5,
    const NumericVector& beta_proposal_sd = NumericVector::create(1.0, 0.05),
    double lambda = 1e-4,
    int K_int = 0,
    double tau_logphi = 1.0,
    double mrf_alpha = 2.0,
    double mrf_beta = 0.5,
    double mrf_bandwidth = 5.0,
    double mrf_gamma = 1.0,
    const NumericVector& phi_prior_vec = NumericVector::create(2.0, 0.1),
    const NumericVector& p_prior_vec = NumericVector::create(1.0, 9.0),
    String p_mode = "global",
    std::string mrf_method = "binary")
{
  int n = x.n_elem;
  int G = y.n_cols;
  
  arma::vec s_use = s_;
  if (s_use.n_elem != (arma::uword)n) {
    s_use = arma::ones<arma::vec>(n);
  }

  // Initialize clusters
  IntegerVector cluster_labels = clone(gene_clusters);
  int C = max(cluster_labels);
  
  std::map<int, std::vector<int>> cluster_genes;
  for (int g = 0; g < G; ++g) {
    int c = cluster_labels[g];
    cluster_genes[c].push_back(g);
  }

  // Initialize change points per cluster
  std::map<int, std::vector<int>> cp_map;

  // Warm start: allow cp_init to be either an integer vector (used for all clusters)
  // or a list of integer vectors (one per cluster, length C).
  if (Rf_isNull(cp_init) || TYPEOF(cp_init) == INTSXP || TYPEOF(cp_init) == REALSXP) {
    std::vector<int> cp_init_vec;
    if (!Rf_isNull(cp_init)) {
      cp_init_vec = as<std::vector<int>>(cp_init);
    }
    for (int c = 1; c <= C; ++c) {
      cp_map[c] = cp_init_vec;
    }
  } else if (TYPEOF(cp_init) == VECSXP) {
    List cp_init_list(cp_init);
    for (int c = 1; c <= C; ++c) {
      if ((c - 1) < cp_init_list.size() && !Rf_isNull(cp_init_list[c - 1])) {
        cp_map[c] = as<std::vector<int>>(cp_init_list[c - 1]);
      } else {
        cp_map[c] = std::vector<int>();
      }
    }
  } else {
    for (int c = 1; c <= C; ++c) {
      cp_map[c] = std::vector<int>();
    }
  }

  // sanitize: sort+unique and drop invalid endpoints
  for (int c = 1; c <= C; ++c) {
    auto &v = cp_map[c];
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    v.erase(std::remove_if(v.begin(), v.end(), [&](int t){ return t <= 0 || t >= n; }), v.end());
  }

  // Initialize parameters
  arma::vec phis(G);
  for (int g = 0; g < G; ++g) {
    phis[g] = 10.0;
  }

  double p_glb = 0.1;
  arma::vec p_spot(n, arma::fill::value(p_glb));
  arma::umat n_mat(n, G, arma::fill::zeros);

  double alpha_p = p_prior_vec[0];
  double beta_p = p_prior_vec[1];

  // Beta lists
  std::vector<std::vector<arma::vec>> beta_lists(G);
  
  // Prior for beta
  arma::vec mu0 = {0.0, 0.0};
  arma::mat Sigma0 = arma::eye(2, 2);
  Sigma0(0,0) = 25.0;
  Sigma0(1,1) = beta_proposal_sd[1] * beta_proposal_sd[1];

  int K_max = (K_int > 0) ? K_int : std::max(3, (int)(n / 20));

  // Storage
  int n_save = n_iter - burnin;
  arma::cube cp_indicator_store(n, C, n_save, arma::fill::zeros);
  arma::cube beta_store_pointwise(G, n, n_save);
  arma::cube intercept_store_pointwise(G, n, n_save);
  arma::cube mu_store_pointwise(G, n, n_save);
  arma::mat p_store(n, n_save);
  arma::mat phi_store(G, n_save);
  arma::vec ll_store(n_save);
  arma::mat k_store(C, n_save);
  arma::vec accept(n_iter, arma::fill::zeros);

  // Lambda function for cluster marginal
  auto cluster_marginal_collapsed_cluster = [&](
      const std::vector<int>& cp,
      const std::vector<int>& g_ids,
      std::vector<std::vector<arma::vec>>& beta_out) -> double
  {
    std::vector<int> brk = cp;
    std::sort(brk.begin(), brk.end());
    brk.insert(brk.begin(), 0);
    brk.push_back(n);

    int K = (int)(brk.size() - 1);
    beta_out.resize(g_ids.size());

    double total_marg = 0.0;
    for (size_t j = 0; j < g_ids.size(); ++j) {
      int g = g_ids[j];
      beta_out[j].resize(K);

      // Work around Armadillo subview_col lacking .elem(): copy the column once per gene
      const arma::vec  y_g = y.col(g);
      const arma::uvec n_g = arma::conv_to<arma::uvec>::from(n_mat.col(g));


      for (int k = 0; k < K; ++k) {
        int a = brk[k], b = brk[k+1];
        arma::uvec idx = arma::regspace<arma::uvec>(a, b-1);
        arma::uvec mask = n_g.elem(idx);

        Rcpp::List lap = laplace_marginal_NB_mask(
          y_g.elem(idx), x.elem(idx), s_use.elem(idx),
          mask, phis[g], mu0, Sigma0
        );

        total_marg += Rcpp::as<double>(lap["log_marginal"]);
        beta_out[j][k] = Rcpp::as<arma::vec>(lap["beta"]);
      }
    }
    return total_marg;
  };

  auto cluster_ll = [&](
      const std::vector<int>& cp,
      const std::vector<std::vector<arma::vec>>& beta_lists_in,
      const std::vector<int>& g_ids,
      const arma::mat& y_in,
      const arma::umat& n_mat_in,
      const arma::vec& x_in,
      const arma::vec& phis_in,
      const arma::vec& s_in) -> double
  {
    double ll = 0.0;
    for (int g : g_ids) {
      ll += zinb_loglik_single_condN(y_in.col(g), n_mat_in.col(g),
                                     x_in, cp, beta_lists_in[g], phis_in[g], s_in);
    }
    return ll;
  };

  // Main MCMC loop
  for (int iter = 0; iter < n_iter; ++iter) {
    if ((iter + 1) % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // Update each cluster
    for (int c = 1; c <= C; ++c) {
      std::vector<int> g_ids = cluster_genes[c];
      if (g_ids.empty()) continue;

      std::vector<int> cp = cp_map[c];
      clean_cp(cp);
      int K_old = (int)cp.size();

      // Current marginal
      std::vector<std::vector<arma::vec>> beta_curr_cluster;
      double ll_current = cluster_marginal_collapsed_cluster(cp, g_ids, beta_curr_cluster);

      // Add MRF prior contribution
      double log_prior_current = dtpois(K_old, K_max, lambda);
      if (neighbor_cp_states.size() > 0) {
        log_prior_current += compute_mrf_prior_log(cp, neighbor_cps_for_cluster(neighbor_cp_states, c), 
                                                   mrf_alpha, mrf_beta, mrf_bandwidth, mrf_gamma, mrf_method);
      }

      // Proposal
      std::vector<int> cp_prop = cp;
      std::string move_type;
      double log_proposal_ratio = 0.0;
      bool proposal_valid = true;

      // Normalize move-type probabilities in current state (depends only on K)
      double p_add_old = prob_add(K_old, K_max);
      double p_rem_old = prob_remove(K_old, K_max);
      double p_mov_old = prob_move(K_old, K_max);
      double total_old = p_add_old + p_rem_old + p_mov_old;
      p_add_old /= total_old; p_rem_old /= total_old; p_mov_old /= total_old;

      double u = R::runif(0, 1);

      // Keep proposal bookkeeping for Hastings ratio
      int A_old = -1, A_new = -1;
      int K_new = K_old;
      int W_old = -1, W_new = -1;
      int moved_old_cp = -1, moved_new_cp = -1;

      if (u < p_add_old) {
        move_type = "add";
        std::vector<int> candidates;
        for (int pos = cp_distance; pos < n - cp_distance; ++pos) {
          bool ok = true;
          for (int cpi : cp) {
            if (std::abs(pos - cpi) < cp_distance) { ok = false; break; }
          }
          if (ok) candidates.push_back(pos);
        }

        A_old = (int)candidates.size();
        if (A_old > 0) {
          int new_cp = candidates[(int)(R::runif(0, 1) * A_old)];
          cp_prop.push_back(new_cp);
          clean_cp(cp_prop);
          K_new = (int)cp_prop.size();

          // Normalize move-type probabilities in proposed state
          double p_add_new = prob_add(K_new, K_max);
          double p_rem_new = prob_remove(K_new, K_max);
          double p_mov_new = prob_move(K_new, K_max);
          double total_new = p_add_new + p_rem_new + p_mov_new;
          p_add_new /= total_new; p_rem_new /= total_new; p_mov_new /= total_new;

          // log(q_reverse/q_forward) for add:
          // q_fwd = p_add_old * 1/A_old
          // q_rev = p_rem_new * 1/K_new
          log_proposal_ratio = slog(p_rem_new) - std::log((double)K_new)
                             - slog(p_add_old) + std::log((double)A_old);
        } else {
          proposal_valid = false;
        }
      }
      else if (u < p_add_old + p_rem_old) {
        move_type = "remove";
        if (K_old > 0) {
          int idx = (int)(R::runif(0, 1) * K_old);
          if (idx >= K_old) idx = K_old - 1;
          cp_prop.erase(cp_prop.begin() + idx);
          clean_cp(cp_prop);
          K_new = (int)cp_prop.size();

          A_new = candidate_count(n, cp_distance, cp_prop);
          if (A_new <= 0) {
            proposal_valid = false;
          } else {
            // Normalize move-type probabilities in proposed state
            double p_add_new = prob_add(K_new, K_max);
            double p_rem_new = prob_remove(K_new, K_max);
            double p_mov_new = prob_move(K_new, K_max);
            double total_new = p_add_new + p_rem_new + p_mov_new;
            p_add_new /= total_new; p_rem_new /= total_new; p_mov_new /= total_new;

            // log(q_reverse/q_forward) for remove:
            // q_fwd = p_rem_old * 1/K_old
            // q_rev = p_add_new * 1/A_new
            log_proposal_ratio = slog(p_add_new) - std::log((double)A_new)
                               - slog(p_rem_old) + std::log((double)K_old);
          }
        } else {
          proposal_valid = false;
        }
      }
      else {
        move_type = "move";
        if (K_old > 0) {
          int idx_cp = (int)(R::runif(0, 1) * K_old);
          if (idx_cp >= K_old) idx_cp = K_old - 1;
          moved_old_cp = cp[idx_cp];

          W_old = move_window_count(cp, idx_cp, n, cp_distance);
          if (W_old <= 0) {
            proposal_valid = false;
          } else {
            // Sample new CP uniformly from its admissible window
            int left_bd  = (idx_cp == 0 ? 0 : cp[idx_cp - 1]) + cp_distance;
            int right_bd = (idx_cp == K_old - 1 ? n : cp[idx_cp + 1]) - cp_distance;

            int offset = (int)(R::runif(0, 1) * W_old);
            int new_cp = -1;
            for (int pos = left_bd; pos < right_bd; ++pos) {
              if (pos == moved_old_cp) continue;
              if (offset == 0) { new_cp = pos; break; }
              --offset;
            }
            if (new_cp < 0) {
              proposal_valid = false;
            } else {
              moved_new_cp = new_cp;
              cp_prop[idx_cp] = new_cp;
              clean_cp(cp_prop);
              K_new = (int)cp_prop.size(); // should equal K_old

              // Find index of moved_new_cp in sorted cp_prop
              int idx_new = -1;
              for (int ii = 0; ii < (int)cp_prop.size(); ++ii) {
                if (cp_prop[ii] == moved_new_cp) { idx_new = ii; break; }
              }
              W_new = move_window_count(cp_prop, idx_new, n, cp_distance);
              if (W_new <= 0) {
                proposal_valid = false;
              } else {
                // Normalize move-type probabilities in proposed state
                double p_add_new = prob_add(K_new, K_max);
                double p_rem_new = prob_remove(K_new, K_max);
                double p_mov_new = prob_move(K_new, K_max);
                double total_new = p_add_new + p_rem_new + p_mov_new;
                p_add_new /= total_new; p_rem_new /= total_new; p_mov_new /= total_new;

                // log(q_reverse/q_forward) for move:
                // q_fwd = p_mov_old * 1/K_old * 1/W_old
                // q_rev = p_mov_new * 1/K_new * 1/W_new
                log_proposal_ratio = slog(p_mov_new) - std::log((double)K_new) - std::log((double)W_new)
                                   - slog(p_mov_old) + std::log((double)K_old) + std::log((double)W_old);
              }
            }
          }
        } else {
          proposal_valid = false;
        }
      }

// MH accept/reject
      bool accepted = false;
      if (proposal_valid) {
        int K_new = (int)cp_prop.size();

        std::vector<std::vector<arma::vec>> beta_prop_cluster;
        double ll_proposed = cluster_marginal_collapsed_cluster(cp_prop, g_ids, beta_prop_cluster);

        double log_prior_proposed = dtpois(K_new, K_max, lambda);
        if (neighbor_cp_states.size() > 0) {
          log_prior_proposed += compute_mrf_prior_log(cp_prop, neighbor_cps_for_cluster(neighbor_cp_states, c),
                                                      mrf_alpha, mrf_beta, mrf_bandwidth, mrf_gamma, mrf_method);
        }

        double log_accept = (ll_proposed - ll_current) 
                          + (log_prior_proposed - log_prior_current)
                          + log_proposal_ratio;

        if (std::isfinite(log_accept) && std::log(R::runif(0,1)) < log_accept) {
          accepted = true;
          cp = cp_prop;
          cp_map[c] = cp;
          accept[iter] = (move_type=="add") ? 1 : (move_type=="remove") ? 2 : 3;
          
          for (size_t j = 0; j < g_ids.size(); ++j) {
            int g = g_ids[j];
            beta_lists[g] = beta_prop_cluster[j];
          }
        } else {
          for (size_t j = 0; j < g_ids.size(); ++j) {
            int g = g_ids[j];
            beta_lists[g] = beta_curr_cluster[j];
          }
        }
      }

      // Phi updates (same as original)
      double logphi_lo = std::log(1.0);
      double logphi_hi = std::log(1000000.0);

      
      // Gene-level collapsed marginal (Laplace) for phi updates: integrates out beta segment-wise
      auto gene_marginal_collapsed = [&](const std::vector<int>& cp_in,
                                         int g,
                                         double phi_val,
                                         std::vector<arma::vec>& beta_out) -> double
      {
        std::vector<int> brk = cp_in;
        std::sort(brk.begin(), brk.end());
        brk.insert(brk.begin(), 0);
        brk.push_back(n);
        int K = (int)(brk.size() - 1);
        beta_out.resize(K);

        // Work around Armadillo subview_col lacking .elem(): copy the column once per call
        const arma::vec  y_g = y.col(g);
        const arma::uvec n_g = arma::conv_to<arma::uvec>::from(n_mat.col(g));

        double total_marg = 0.0;

        for (int k = 0; k < K; ++k) {
          int a = brk[k], b = brk[k+1];
          arma::uvec idx = arma::regspace<arma::uvec>(a, b-1);
          arma::uvec mask = n_g.elem(idx);

          Rcpp::List lap = laplace_marginal_NB_mask(
            y_g.elem(idx), x.elem(idx), s_use.elem(idx),
            mask, phi_val, mu0, Sigma0
          );
          total_marg += Rcpp::as<double>(lap["log_marginal"]);
          beta_out[k] = Rcpp::as<arma::vec>(lap["beta"]);
        }
        return total_marg;
      };

      for (int g : g_ids) {
        double log_phi_curr = std::log(std::max(phis[g], 1e-12));
        double log_phi_prop = rtruncnorm(log_phi_curr, tau_logphi, logphi_lo, logphi_hi);
        double phi_prop = std::exp(log_phi_prop);

        std::vector<arma::vec> beta_curr, beta_prop;
        double ll_curr = gene_marginal_collapsed(cp, g, phis[g], beta_curr);
        double ll_prop = gene_marginal_collapsed(cp, g, phi_prop, beta_prop);

        double log_prior_curr = R::dgamma(phis[g], phi_prior_vec[0], 1.0/phi_prior_vec[1], true);
        double log_prior_prop = R::dgamma(phi_prop, phi_prior_vec[0], 1.0/phi_prior_vec[1], true);

        double log_q_fwd = dtruncnorm_log(log_phi_prop, log_phi_curr, tau_logphi, logphi_lo, logphi_hi);
        double log_q_back = dtruncnorm_log(log_phi_curr, log_phi_prop, tau_logphi, logphi_lo, logphi_hi);

        double log_acc = (ll_prop - ll_curr) + (log_prior_prop - log_prior_curr)
                       + (log_q_back - log_q_fwd) + (log_phi_prop - log_phi_curr);

        bool accepted_phi = false;
        if (std::isfinite(log_acc) && std::log(R::runif(0,1)) < log_acc) {
          phis[g] = phi_prop;
          beta_lists[g] = beta_prop;
          accepted_phi = true;
        }
        if (!accepted_phi) {
          beta_lists[g] = beta_curr;
        }
      }
    }

    // Update p and n_mat (same as original)
    arma::mat mu(n, G, arma::fill::zeros);
    for (int g = 0; g < G; ++g) {
      int c = cluster_labels[g];
      std::vector<int> seg = cp_map[c];
      std::sort(seg.begin(), seg.end());
      seg.insert(seg.begin(), 0);
      seg.push_back(n);

      for (size_t k = 0; k < seg.size() - 1; ++k) {
        int a = seg[k], b = seg[k+1];
        const arma::vec& beta = beta_lists[g][k];
        for (int i = a; i < b; ++i) {
          double eta = beta[0] + beta[1] * x[i];
          eta = clamp(eta, -50.0, 50.0);
          mu(i,g) = std::exp(eta) * s_use[i];
          mu(i,g) = std::max(mu(i,g), 1e-12);
        }
      }
    }

    if (p_mode == "global") {
      unsigned n_zero = 0, s_zero = 0;
      for (int i = 0; i < n; ++i) {
        for (int g = 0; g < G; ++g) {
          if (y(i,g) == 0) {
            ++n_zero;
            double nb0 = R::dnbinom_mu(0, phis[g], mu(i,g), false);
            double denom = p_glb + (1.0 - p_glb) * nb0;
            denom = std::max(denom, 1e-300);
            double pi_z1 = p_glb / denom;
            n_mat(i,g) = (R::runif(0,1) < pi_z1);
            if (n_mat(i,g)) ++s_zero;
          } else {
            n_mat(i,g) = 0;
          }
        }
      }
      unsigned N_total = (unsigned)n * (unsigned)G;  // all (i,g)
        double a = alpha_p + (double)s_zero;
        double b = beta_p  + (double)(N_total - s_zero);
        p_glb = R::rbeta(a, b);
        p_spot.fill(p_glb);
    }

    // Store results after burnin
    if (iter >= burnin) {
      int s_idx = iter - burnin;

      for (int c = 1; c <= C; ++c) {
        for (int idx : cp_map[c]) {
          if (idx >= 0 && idx < n) cp_indicator_store(idx, c-1, s_idx) = 1.0;
        }
      }

      for (int g = 0; g < G; ++g) {
        int c = cluster_labels[g];
        std::vector<int> brk = cp_map[c];
        std::sort(brk.begin(), brk.end());
        brk.insert(brk.begin(), 0);
        brk.push_back(n);

        for (size_t k = 0; k < beta_lists[g].size(); ++k) {
          int a = brk[k], b = brk[k+1];
          for (int i = a; i < b; ++i) {
            intercept_store_pointwise(g, i, s_idx) = beta_lists[g][k][0];
            beta_store_pointwise(g, i, s_idx) = beta_lists[g][k][1];

            double eta = beta_lists[g][k][0] + beta_lists[g][k][1] * x[i];
            eta = clamp(eta, -50.0, 50.0);
            mu_store_pointwise(g, i, s_idx) = std::exp(eta) * s_use[i];
          }
        }
      }

      p_store.col(s_idx) = p_spot;
      phi_store.col(s_idx) = phis;

      double total_ll = 0.0;
      for (int c = 1; c <= C; ++c) {
        total_ll += cluster_ll(cp_map[c], beta_lists, cluster_genes[c], y, n_mat, x, phis, s_use);
      }
      ll_store[s_idx] = total_ll;

      for (int c = 1; c <= C; ++c) {
        k_store(c-1, s_idx) = (double)cp_map[c].size();
      }
    }
  }

  return List::create(
    Named("beta") = beta_store_pointwise,
    Named("intercept") = intercept_store_pointwise,
    Named("mu_loc") = mu_store_pointwise,
    Named("cp_indicator") = cp_indicator_store,
    Named("p") = p_store,
    Named("phi") = phi_store,
    Named("log_likelihoods") = ll_store,
    Named("K") = k_store,
    Named("accept") = accept
  );
}