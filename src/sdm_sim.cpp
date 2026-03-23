#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdio>

using namespace Rcpp;
using namespace arma;

// ── Conversión segura de dgCMatrix → arma::sp_mat ───────────────────────────
sp_mat safe_spmat_conversion(SEXP mat) {
  try {
    S4 dgC(mat);
    IntegerVector dims  = dgC.slot("Dim");
    IntegerVector i_vec = dgC.slot("i");
    IntegerVector p_vec = dgC.slot("p");
    NumericVector x_vec = dgC.slot("x");

    const uword n_rows = (uword)dims[0];
    const uword n_cols = (uword)dims[1];
    const uword nnz    = (uword)x_vec.size();

    // Retorno para matrices vacías
    if (n_rows == 0 || n_cols == 0) {
      return sp_mat(0, 0);
    }

    if (nnz == 0) {
      return sp_mat(n_rows, n_cols);
    }

    // Validar tamaños
    if (p_vec.size() != (size_t)n_cols + 1) {
      Rcpp::warning("p_vec size mismatch: expected %llu, got %llu",
                    (unsigned long long)n_cols + 1,
                    (unsigned long long)p_vec.size());
      return sp_mat(n_rows, n_cols);
    }

    // Construcción por batch insert — segura bajo UBSAN
    arma::urowvec row_idx(nnz);
    arma::urowvec col_idx(nnz);
    arma::vec     values(nnz);

    uword k = 0;
    const int n_cols_int = (int)n_cols;

    for (int col = 0; col < n_cols_int; col++) {
      int start = p_vec[col];
      int end = (col + 1 < n_cols_int) ? p_vec[col + 1] : (int)n_cols;

      // Validar límites
      if (start < 0 || end > (int)i_vec.size()) {
        Rcpp::warning("Invalid index range at column %d: [%d, %d)", col, start, end);
        continue;
      }

      for (int j = start; j < end && k < nnz; j++, k++) {
        int row_idx_val = i_vec[j];
        if (row_idx_val < 0 || row_idx_val >= (int)n_rows) {
          Rcpp::warning("Row index %d out of range at column %d", row_idx_val, col);
          continue;
        }

        row_idx(k) = (uword)row_idx_val;
        col_idx(k) = (uword)col;
        values(k)  = x_vec[j];
      }
    }

    if (k == 0) {
      return sp_mat(n_rows, n_cols);
    }

    // Redimensionar si no se usaron todos los elementos
    if (k < nnz) {
      row_idx.resize(k);
      col_idx.resize(k);
      values.resize(k);
    }

    arma::umat locations(2, k);
    locations.row(0) = row_idx;
    locations.row(1) = col_idx;

    return sp_mat(locations, values, n_rows, n_cols);

  } catch (std::exception& e) {
    Rcpp::stop("Error in sparse matrix conversion: %s", e.what());
  }
  return sp_mat(0, 0);
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

  // ── Conversión segura de inputs ─────────────────────────────────────────────
  sp_mat A_mat = safe_spmat_conversion(A);
  sp_mat M_mat = safe_spmat_conversion(M_orig);
  sp_mat g0    = safe_spmat_conversion(g0_input);
  const uword n = A_mat.n_rows;

  // ── Validaciones ────────────────────────────────────────────────────────────
  if (n == 0) {
    Rf_error("set_A sparse matrix has zero rows. "
               "Check the threshold in model2sparse(): it may be too high, "
               "leaving no valid cells in the model.");
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
  // ────────────────────────────────────────────────────────────────────────────

  // Semilla aleatoria
  if (stochastic_dispersal) {
    Rcpp::Environment baseEnv("package:base");
    Rcpp::Function setSeed = baseEnv["set.seed"];
    setSeed(0);
  }

  // Suitability binaria
  vec binary_suit(n, fill::zeros);
  for (uword i = 0; i < n; ++i) {
    binary_suit(i) = (A_mat(i, i) > 0) ? 1.0 : 0.0;
  }

  // Probabilidades de suitability
  vec suit_probs = as<vec>(suit_values);
  if (disp_prop2_suitability) {
    suit_probs *= disper_prop;
  }

  // Cache de vecinos
  std::vector<std::vector<uword>> neighbors(n);
  CharacterVector list_names = adj_list.names();

  for (R_xlen_t i = 0; i < adj_list.size(); i++) {
    if (adj_list[i] == R_NilValue) continue;

    // Validar que el nombre existe y es convertible
    if (i >= list_names.size()) continue;

    std::string name = as<std::string>(list_names[i]);
    int src_cell_1based;
    try {
      src_cell_1based = std::stoi(name);
    } catch (...) {
      continue;
    }

    uword src_cell = static_cast<uword>(src_cell_1based - 1);
    if (src_cell >= n) continue;

    DataFrame df = as<DataFrame>(adj_list[i]);
    if (!df.containsElementNamed("ToNonNaCell")) continue;

    NumericVector to_non_na = df["ToNonNaCell"];

    if (to_non_na.size() > 0) {
      neighbors[src_cell].reserve(to_non_na.size());

      for (R_xlen_t j = 0; j < to_non_na.size(); j++) {
        double nb_val = to_non_na[j];
        if (NumericVector::is_na(nb_val)) continue;

        uword nb_index = static_cast<uword>(nb_val - 1);
        if (nb_index < n) {
          neighbors[src_cell].push_back(nb_index);
        }
      }
    }
  }

  // Pre-asignación de salida
  List sdm(nsteps + 1);
  sdm[0] = g0;

  // Setup determinístico
  sp_mat AMA_base;
  if (!stochastic_dispersal) {
    // Validar que las multiplicaciones son posibles
    if (A_mat.n_cols != M_mat.n_rows || M_mat.n_cols != A_mat.n_rows) {
      Rf_error("Matrix dimensions incompatible for AMA multiplication");
    }
    AMA_base = A_mat * M_mat * A_mat;
  }

  // Setup barra de progreso
  const int bar_width = 50;
  int progress_interval = std::max(1, nsteps / bar_width);
  if (progress_bar) {
    Rprintf("Simulation progress:\n[");
    for (int i = 0; i < bar_width; i++) Rprintf(" ");
    Rprintf("]\r[");
  }

  // ── Loop principal de simulación ────────────────────────────────────────────
  sp_mat g0_next;
  for (int step = 1; step <= nsteps; ++step) {
    if (progress_bar && (step % 10 == 0)) Rcpp::checkUserInterrupt();

    g0_next = g0;

    if (stochastic_dispersal) {
      const uvec occ_locs = find(g0 > 0);
      const uword n_occ = occ_locs.n_elem;

      if (n_occ > 0) {
        vec rand_values = randu<vec>(n);

        for (uword i = 0; i < n_occ; ++i) {
          const uword cell = occ_locs(i);
          if (cell >= neighbors.size()) continue;

          const std::vector<uword>& cell_neighbors = neighbors[cell];
          const uword n_nb = cell_neighbors.size();

          for (uword j = 0; j < n_nb; ++j) {
            const uword nb_index = cell_neighbors[j];

            // Validar que nb_index está en rango
            if (nb_index >= n) continue;

            // Verificar suitability
            if (binary_suit(nb_index) < 0.5) continue;

            const double prob = disp_prop2_suitability ?
            suit_probs(nb_index) : disper_prop;

            // Verificar que prob es válida
            if (prob <= 0.0) continue;
            if (prob >= 1.0 || rand_values(nb_index) < prob) {
              g0_next(nb_index, 0) = 1.0;
            }
          }
        }
      }

    } else {
      // Dispersal determinístico
      sp_mat temp = AMA_base * g0;
      g0_next = temp + g0;
      // Transformación segura
      for (uword i = 0; i < g0_next.n_rows; ++i) {
        if (g0_next(i, 0) > 0) {
          g0_next(i, 0) = 1.0;
        }
      }
    }

    g0 = g0_next;
    sdm[step] = g0;

    // Actualización barra de progreso
    if (progress_bar && (step % progress_interval == 0 || step == nsteps)) {
      int pos = static_cast<int>(bar_width * static_cast<double>(step) / nsteps);
      if (pos > bar_width) pos = bar_width;
      Rprintf("\r[");
      for (int i = 0; i < pos; i++) Rprintf("=");
      if (pos < bar_width) Rprintf(">");
      for (int i = pos + 1; i < bar_width; i++) Rprintf(" ");
      Rprintf("] %3d%%", static_cast<int>(100.0 * step / nsteps));
    }
  }

  if (progress_bar) Rprintf("\n");
  return sdm;
}
