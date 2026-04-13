#define ARMA_64BIT_WORD
#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <string>
#include <cstdio>

using namespace Rcpp;
using namespace arma;

// ── Conversión segura de dgCMatrix → arma::sp_mat ───────────────────────────
// Usa el constructor CSC nativo de Armadillo con punteros a memoria ya
// inicializada por std::vector — evita fill_zeros sobre puntero nulo (UBSAN).
sp_mat safe_spmat_conversion(SEXP mat) {
  if (Rf_isNull(mat)) return sp_mat(0, 0);

  S4 dgC(mat);
  IntegerVector dims  = dgC.slot("Dim");
  IntegerVector i_vec = dgC.slot("i");
  IntegerVector p_vec = dgC.slot("p");
  NumericVector x_vec = dgC.slot("x");

  const uword n_rows = (uword)dims[0];
  const uword n_cols = (uword)dims[1];
  const uword nnz    = (uword)x_vec.size();

  if (n_rows == 0 || n_cols == 0 || nnz == 0) {
    return sp_mat(n_rows, n_cols);
  }

  // Construir directamente usando arma::Col (más eficiente)
  arma::Col<uword> row_ind(nnz);
  arma::Col<uword> col_ptr(n_cols + 1);
  arma::Col<double> vals(nnz);

  for (uword k = 0; k < nnz; k++) {
    row_ind[k] = (uword)i_vec[k];
    vals[k]    = x_vec[k];
  }
  for (uword k = 0; k <= n_cols; k++) {
    col_ptr[k] = (uword)p_vec[k];
  }

  return sp_mat(row_ind, col_ptr, vals, n_rows, n_cols);
}
// ── Verificación de matriz sparse válida ─────────────────────────────────────
bool is_valid_sparse_matrix(const sp_mat& mat) {
  return mat.n_rows > 0 && mat.n_cols > 0 &&
    mat.n_nonzero <= mat.n_rows * mat.n_cols;
}

// ── Multiplicación segura de sparse matrices ─────────────────────────────────
sp_mat safe_matrix_multiply(const sp_mat& A, const sp_mat& B) {
  if (A.n_cols != B.n_rows) {
    Rf_error("Matrix dimensions incompatible: %llux%llu * %llux%llu",
             (unsigned long long)A.n_rows, (unsigned long long)A.n_cols,
             (unsigned long long)B.n_rows, (unsigned long long)B.n_cols);
  }
  if (A.n_nonzero == 0 || B.n_nonzero == 0) {
    return sp_mat(A.n_rows, B.n_cols);
  }
  return A * B;
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

  // ── Conversión segura de inputs ──────────────────────────────────────────
  sp_mat A_mat = safe_spmat_conversion(A);
  sp_mat M_mat = safe_spmat_conversion(M_orig);
  sp_mat g0    = safe_spmat_conversion(g0_input);
  const uword n = A_mat.n_rows;

  // ── Validaciones ─────────────────────────────────────────────────────────
  if (n == 0) {
    Rf_error("set_A sparse matrix has zero rows. "
               "Check the threshold in model2sparse(): "
               "it may be too high, leaving no valid cells in the model.");
  }
  if (A_mat.n_rows != A_mat.n_cols) {
    Rf_error("set_A must be a square matrix (n x n).");
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
    Rf_error("nsteps must be a positive integer.");
  }
  if (disper_prop < 0.0 || disper_prop > 1.0) {
    Rf_error("disper_prop must be a value between 0 and 1.");
  }
  if (g0.n_nonzero == 0) {
    Rf_error("initial_points has no occupied cells. "
               "Check that occurrence coordinates fall within the model extent.");
  }

  // ── Semilla aleatoria ────────────────────────────────────────────────────
  if (stochastic_dispersal) {
    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(123);
  }

  // ── Suitability binaria — std::vector evita fill_zeros de Armadillo ──────
  std::vector<double> binary_suit_v(n, 0.0);
  for (uword i = 0; i < n; ++i) {
    binary_suit_v[i] = (A_mat(i, i) > 0) ? 1.0 : 0.0;
  }

  std::vector<double> suit_probs_v(n);
  for (uword i = 0; i < n; ++i) {
    suit_probs_v[i] = disp_prop2_suitability ?
    (double)suit_values[i] * disper_prop : disper_prop;
  }

  // ── Cache de vecinos ─────────────────────────────────────────────────────
  std::vector<std::vector<uword>> neighbors(n);

  if (adj_list.size() > 0) {
    CharacterVector list_names = adj_list.names();

    for (R_xlen_t i = 0; i < adj_list.size(); i++) {
      if (adj_list[i] == R_NilValue) continue;
      if (i >= list_names.size()) continue;

      std::string name = as<std::string>(list_names[i]);
      int src_cell_1based = 0;
      try { src_cell_1based = std::stoi(name); } catch (...) { continue; }

      uword src_cell = static_cast<uword>(src_cell_1based - 1);
      if (src_cell >= n) continue;

      try {
        DataFrame df = as<DataFrame>(adj_list[i]);
        if (!df.containsElementNamed("ToNonNaCell")) continue;

        NumericVector to_non_na = df["ToNonNaCell"];
        for (R_xlen_t j = 0; j < to_non_na.size(); j++) {
          double nb_val = to_non_na[j];
          if (Rcpp::NumericVector::is_na(nb_val)) continue;
          uword nb_index = static_cast<uword>(nb_val - 1);
          if (nb_index < n) neighbors[src_cell].push_back(nb_index);
        }
      } catch (...) { continue; }
    }
  }

  // ── Pre-asignación de salida ─────────────────────────────────────────────
  List sdm(nsteps + 1);
  sdm[0] = g0;

  // ── Setup determinístico ─────────────────────────────────────────────────
  sp_mat AMA_base;
  if (!stochastic_dispersal) {
    sp_mat AM = safe_matrix_multiply(A_mat, M_mat);
    AMA_base  = safe_matrix_multiply(AM, A_mat);
  }

  // ── Barra de progreso ────────────────────────────────────────────────────
  const int bar_width = 50;
  int progress_interval = std::max(1, nsteps / bar_width);
  if (progress_bar) {
    Rprintf("Simulation progress:\n[");
    for (int i = 0; i < bar_width; i++) Rprintf(" ");
    Rprintf("]\r[");
    R_FlushConsole();
  }

  // ── Loop principal ───────────────────────────────────────────────────────
  sp_mat g0_next;
  for (int step = 1; step <= nsteps; ++step) {
    if (step % 10 == 0) Rcpp::checkUserInterrupt();

    g0_next = g0;

    if (stochastic_dispersal) {
      const uvec occ_locs = find(g0 > 0);
      const uword n_occ = occ_locs.n_elem;

      if (n_occ > 0) {
        // Números aleatorios vía RNG de R — no usa randu<vec> de Armadillo
        Rcpp::NumericVector rv = Rcpp::runif((int)n);
        const double* rand_ptr = rv.begin();

        for (uword i = 0; i < n_occ; ++i) {
          const uword cell = occ_locs(i);
          if (cell >= neighbors.size()) continue;

          const std::vector<uword>& cell_neighbors = neighbors[cell];
          for (uword j = 0; j < cell_neighbors.size(); ++j) {
            const uword nb = cell_neighbors[j];
            if (nb >= n) continue;
            if (binary_suit_v[nb] < 0.5) continue;

            const double prob = suit_probs_v[nb];
            if (prob >= 1.0 || rand_ptr[nb] < prob) {
              g0_next(nb, 0) = 1.0;
            }
          }
        }
      }

    } else {
      if (AMA_base.n_rows > 0) {
        sp_mat temp = safe_matrix_multiply(AMA_base, g0);
        g0_next = temp + g0;
        // Binarizar iterando solo elementos no-cero
        for (sp_mat::iterator it = g0_next.begin(); it != g0_next.end(); ++it) {
          if (*it > 0) *it = 1.0;
        }
      }
    }

    g0 = g0_next;
    sdm[step] = g0;

    if (progress_bar && (step % progress_interval == 0 || step == nsteps)) {
      int pos = static_cast<int>(bar_width * static_cast<double>(step) / nsteps);
      if (pos > bar_width) pos = bar_width;
      Rprintf("\r[");
      for (int i = 0; i < pos; i++) Rprintf("=");
      if (pos < bar_width) Rprintf(">");
      for (int i = pos + 1; i < bar_width; i++) Rprintf(" ");
      Rprintf("] %3d%%", static_cast<int>(100.0 * step / nsteps));
      R_FlushConsole();
    }
  }

  if (progress_bar) { Rprintf("\n"); R_FlushConsole(); }
  return sdm;
}
