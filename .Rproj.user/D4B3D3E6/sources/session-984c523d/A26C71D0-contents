#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdio>

using namespace Rcpp;
using namespace arma;

struct NeighborCache {
  std::vector<uvec> neighbors;
  uvec cell_indices;
};

// Safe sparse matrix multiplication
sp_mat safe_spmat_mult(const sp_mat& A, const sp_mat& B) {
  if (A.n_cols != B.n_rows) {
    char msg[256];
    snprintf(msg, sizeof(msg),
             "Matrix dimension mismatch: %llux%llu vs %llux%llu",
             (unsigned long long)A.n_rows, (unsigned long long)A.n_cols,
             (unsigned long long)B.n_rows, (unsigned long long)B.n_cols);
    throw std::runtime_error(msg);
  }

  try {
    return A * B;
  } catch (const std::bad_alloc& e) {
    throw std::runtime_error("Memory allocation failed during matrix multiplication");
  } catch (...) {
    throw std::runtime_error("Unknown error during matrix multiplication");
  }
}

// Custom conversion function
sp_mat safe_spmat_conversion(SEXP mat) {
  S4 dgC(mat);
  IntegerVector dims = dgC.slot("Dim");
  IntegerVector i = dgC.slot("i");
  IntegerVector p = dgC.slot("p");
  NumericVector x = dgC.slot("x");

  const uword n_rows = dims[0];
  const uword n_cols = dims[1];

  sp_mat result(n_rows, n_cols);

  for (int col = 0; col < n_cols; col++) {
    int start = p[col];
    int end = p[col+1];
    for (int j = start; j < end; j++) {
      uword row_idx = i[j];
      if (row_idx >= n_rows) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "Row index %llu out of bounds (max %llu) at col %d, position %d",
                 (unsigned long long)row_idx, (unsigned long long)n_rows, col, j);
        throw std::runtime_error(msg);
      }
      result(row_idx, col) = x[j];
    }
  }
  return result;
}

// [[Rcpp::export]]
List sdm_sim_rcpp(SEXP A, SEXP M_orig, SEXP g0_input,
                  const Rcpp::NumericVector& suit_values,
                  const Rcpp::List& adj_list,
                  int nsteps,
                  bool stochastic_dispersal,
                  bool disp_prop2_suitability,
                  double disper_prop,
                  bool progress_bar) {
  // Set random seed if needed
  if(stochastic_dispersal){
    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(0);
  }

  // Convert inputs with dimension checks
  sp_mat A_mat = Rcpp::as<arma::sp_mat>(A);
  sp_mat M_mat = Rcpp::as<arma::sp_mat>(M_orig);
  sp_mat g0 = Rcpp::as<arma::sp_mat>(g0_input);
  const uword n = A_mat.n_rows;

  if (g0.n_rows != n || g0.n_cols != 1) {
    stop("Initial state vector must be Nx1");
  }

  // Binary suitability
  vec binary_suit(n, fill::zeros);
  for(uword i = 0; i < n; ++i) {
    binary_suit(i) = (A_mat(i, i) > 0) ? 1.0 : 0.0;
  }

  // Suitability probabilities
  vec suit_probs = as<vec>(suit_values);
  if(disp_prop2_suitability) {
    suit_probs *= disper_prop;
  }

  // Build neighbor cache
  std::vector<std::vector<uword>> neighbors(n);
  CharacterVector list_names = adj_list.names();

  for (R_xlen_t i = 0; i < adj_list.size(); i++) {
    if (adj_list[i] == R_NilValue) continue;

    int src_cell_1based = std::stoi(as<std::string>(list_names[i]));
    uword src_cell = static_cast<uword>(src_cell_1based - 1);

    if (src_cell >= n) continue;

    DataFrame df = as<DataFrame>(adj_list[i]);
    if (!df.containsElementNamed("ToNonNaCell")) continue;

    NumericVector to_non_na = df["ToNonNaCell"];
    neighbors[src_cell].reserve(to_non_na.size());

    for (R_xlen_t j = 0; j < to_non_na.size(); j++) {
      uword nb_index = static_cast<uword>(to_non_na[j]-1);
      if (nb_index < n) {
        neighbors[src_cell].push_back(nb_index);
      }
    }
  }

  // Pre-allocate output
  List sdm(nsteps + 1);
  sdm[0] = g0;

  // Deterministic setup
  sp_mat AMA_base;
  if (!stochastic_dispersal) {
    AMA_base = A_mat * M_mat * A_mat;
  }

  // Progress setup
  const int bar_width = 50;
  int progress_interval = std::max(1, nsteps / bar_width);
  if (progress_bar) {
    Rprintf("Simulation progress:\n[");
    for (int i = 0; i < bar_width; i++) Rprintf(" ");
    Rprintf("]\r[");
  }
  // Main simulation loop
  sp_mat g0_next;
  for (int step = 1; step <= nsteps; ++step) {
    if (progress_bar && (step % 10 == 0)) Rcpp::checkUserInterrupt();

    g0_next = g0;

    if (stochastic_dispersal) {
      // Generate random numbers in bulk
      vec rand_values(n, fill::randu);
      const uvec& occ_locs = find(g0 > 0);
      const uword n_occ = occ_locs.n_elem;

      for (uword i = 0; i < n_occ; ++i) {
        const uword cell = occ_locs(i);
        if (cell >= neighbors.size()) continue;
        const std::vector<uword>& cell_neighbors = neighbors[cell];
        const uword n_nb = cell_neighbors.size();

        for (uword j = 0; j < n_nb; ++j) {
          const uword nb_index = cell_neighbors[j];

          if (binary_suit(nb_index) < 0.5) continue;
          const double prob = disp_prop2_suitability ?
          suit_probs(nb_index) : disper_prop;

          if (rand_values(nb_index) < prob) {
            // Use rand_values[j] for THIS specific dispersal event
            g0_next(nb_index, 0) = 1.0;
          }
        }
      }
    } else {
      // Deterministic dispersal
      g0_next = AMA_base * g0 + g0;
      g0_next.transform([](double val) { return val > 0 ? 1.0 : 0.0; });
    }

    g0 = g0_next;
    sdm[step] = g0;

    // Progress update
    if (progress_bar && (step % progress_interval == 0 || step == nsteps)) {
      int pos = static_cast<int>(bar_width * static_cast<double>(step) / nsteps);
      Rprintf("\r[");
      for (int i = 0; i < pos; i++) Rprintf("=");
      if (pos < bar_width) Rprintf(">");
      for (int i = pos+1; i < bar_width; i++) Rprintf(" ");
      Rprintf("] %3d%%", static_cast<int>(100.0 * step / nsteps));
    }
  }

  if (progress_bar) Rprintf("\n");
  return sdm;
}
