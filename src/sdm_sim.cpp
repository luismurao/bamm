#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdio>

using namespace Rcpp;
using namespace arma;

// Safe sparse matrix multiplication
sp_mat safe_spmat_mult(const sp_mat& A, const sp_mat& B) {
  if (A.n_cols != B.n_rows) {
    Rf_error("Matrix dimension mismatch: %llux%llu vs %llux%llu",
             (unsigned long long)A.n_rows, (unsigned long long)A.n_cols,
             (unsigned long long)B.n_rows, (unsigned long long)B.n_cols);
  }
  try {
    return A * B;
  } catch (const std::bad_alloc&) {
    Rf_error("Memory allocation failed during matrix multiplication.");
  } catch (...) {
    Rf_error("Unknown error during matrix multiplication.");
  }
  return sp_mat(); // unreachable, silences compiler warning
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

  // ── Convert inputs ──────────────────────────────────────────────────────────
  sp_mat A_mat = Rcpp::as<arma::sp_mat>(A);
  sp_mat M_mat = Rcpp::as<arma::sp_mat>(M_orig);
  sp_mat g0    = Rcpp::as<arma::sp_mat>(g0_input);
  const uword n = A_mat.n_rows;

  // ── Input validation (Rf_error para compatibilidad con UBSAN) ───────────────
  if (n == 0) {
    Rf_error("set_A sparse matrix has zero rows. "
               "Check the threshold in model2sparse(): it may be too high, "
               "leaving no valid cells in the model.");
  }
  if (A_mat.n_rows != A_mat.n_cols) {
    Rf_error("set_A must be a square matrix (n x n). "
               "Got %llux%llu.",
               (unsigned long long)A_mat.n_rows,
               (unsigned long long)A_mat.n_cols);
  }
  if (M_mat.n_rows == 0 || M_mat.n_cols == 0) {
    Rf_error("set_M sparse matrix has zero dimensions. "
               "Check the output of adj_mat().");
  }
  if (M_mat.n_rows != n || M_mat.n_cols != n) {
    Rf_error("set_M dimensions (%llux%llu) must match set_A (%llux%llu).",
             (unsigned long long)M_mat.n_rows,
             (unsigned long long)M_mat.n_cols,
             (unsigned long long)n,
             (unsigned long long)n);
  }
  if (g0.n_rows != n || g0.n_cols != 1) {
    Rf_error("initial_points must be a %llux1 sparse vector, got %llux%llu.",
             (unsigned long long)n,
             (unsigned long long)g0.n_rows,
             (unsigned long long)g0.n_cols);
  }
  if ((uword)suit_values.size() != n) {
    Rf_error("suit_values length (%llu) must match matrix dimension (%llu).",
             (unsigned long long)suit_values.size(),
             (unsigned long long)n);
  }
  if (nsteps <= 0) {
    Rf_error("nsteps must be a positive integer (got %d).", nsteps);
  }
  if (disper_prop < 0.0 || disper_prop > 1.0) {
    Rf_error("disper_prop must be between 0 and 1 (got %f).", disper_prop);
  }
  if (g0.n_nonzero == 0) {
    Rf_error("initial_points has no occupied cells. "
               "Check that occurrence coordinates fall within the model extent.");
  }
  // ────────────────────────────────────────────────────────────────────────────

  // Random seed for stochastic mode
  if (stochastic_dispersal) {
    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(0);
  }

  // Binary suitability vector (diagonal of A indicates suitable cells)
  vec binary_suit(n, fill::zeros);
  for (uword i = 0; i < n; ++i) {
    binary_suit(i) = (A_mat(i, i) > 0) ? 1.0 : 0.0;
  }

  // Suitability probabilities
  vec suit_probs = as<vec>(suit_values);
  if (disp_prop2_suitability) {
    suit_probs *= disper_prop;
  }

  // Build neighbor lookup from adj_list
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
      uword nb_index = static_cast<uword>(to_non_na[j] - 1);
      if (nb_index < n) {
        neighbors[src_cell].push_back(nb_index);
      }
    }
  }

  // Pre-allocate output list (step 0 = initial state)
  List sdm(nsteps + 1);
  sdm[0] = g0;

  // Precompute AMA for deterministic mode
  sp_mat AMA_base;
  if (!stochastic_dispersal) {
    AMA_base = A_mat * M_mat * A_mat;
  }

  // Progress bar setup
  const int bar_width = 50;
  int progress_interval = std::max(1, nsteps / bar_width);
  if (progress_bar) {
    Rprintf("Simulation progress:\n[");
    for (int i = 0; i < bar_width; i++) Rprintf(" ");
    Rprintf("]\r[");
  }

  // ── Main simulation loop ────────────────────────────────────────────────────
  sp_mat g0_next;
  for (int step = 1; step <= nsteps; ++step) {

    if (progress_bar && (step % 10 == 0)) Rcpp::checkUserInterrupt();

    g0_next = g0;

    if (stochastic_dispersal) {
      const uvec occ_locs = find(g0 > 0);
      const uword n_occ = occ_locs.n_elem;

      if (n_occ > 0) {
        vec rand_values(n, fill::randu);  // safe: n > 0 guaranteed above

        for (uword i = 0; i < n_occ; ++i) {
          const uword cell = occ_locs(i);
          if (cell >= neighbors.size()) continue;

          const std::vector<uword>& cell_neighbors = neighbors[cell];
          const uword n_nb = cell_neighbors.size();

          for (uword j = 0; j < n_nb; ++j) {
            const uword nb_index = cell_neighbors[j];
            if (binary_suit(nb_index) < 0.5) continue;

            const double prob = disp_prop2_suitability
            ? suit_probs(nb_index)
              : disper_prop;

            if (rand_values(nb_index) < prob) {
              g0_next(nb_index, 0) = 1.0;
            }
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

    // Update progress bar
    if (progress_bar && (step % progress_interval == 0 || step == nsteps)) {
      int pos = static_cast<int>(bar_width * static_cast<double>(step) / nsteps);
      Rprintf("\r[");
      for (int i = 0; i < pos; i++) Rprintf("=");
      if (pos < bar_width) Rprintf(">");
      for (int i = pos + 1; i < bar_width; i++) Rprintf(" ");
      Rprintf("] %3d%%", static_cast<int>(100.0 * step / nsteps));
    }
  }
  // ────────────────────────────────────────────────────────────────────────────

  if (progress_bar) Rprintf("\n");
  return sdm;
}
