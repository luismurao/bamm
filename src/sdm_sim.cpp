#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdio>

using namespace Rcpp;
using namespace arma;

// Función auxiliar para verificar si una matriz es válida
bool is_valid_sparse_matrix(const sp_mat& mat) {
  return mat.n_rows > 0 && mat.n_cols > 0 && mat.n_nonzero <= mat.n_rows * mat.n_cols;
}

// ── Conversión segura de dgCMatrix → arma::sp_mat ───────────────────────────
sp_mat safe_spmat_conversion(SEXP mat) {
  // Verificar que no sea NULL
  if (Rf_isNull(mat)) {
    Rf_warning("NULL matrix encountered, returning empty matrix");
    return sp_mat(0, 0);
  }

  try {
    S4 dgC(mat);

    // Verificar slots
    if (!dgC.hasSlot("Dim") || !dgC.hasSlot("i") || !dgC.hasSlot("p") || !dgC.hasSlot("x")) {
      Rf_error("Invalid dgCMatrix: missing required slots");
    }

    IntegerVector dims  = dgC.slot("Dim");
    IntegerVector i_vec = dgC.slot("i");
    IntegerVector p_vec = dgC.slot("p");
    NumericVector x_vec = dgC.slot("x");

    if (dims.size() < 2) {
      Rf_error("Invalid dimensions vector");
    }

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

    // Validar tamaños de vectores
    if ((size_t)p_vec.size() != (size_t)n_cols + 1) {
      Rf_warning("p_vec size mismatch: expected %llu, got %llu",
                 (unsigned long long)n_cols + 1,
                 (unsigned long long)p_vec.size());
      return sp_mat(n_rows, n_cols);
    }

    if ((size_t)i_vec.size() != nnz) {
      Rf_warning("i_vec size mismatch: expected %llu, got %llu",
                 (unsigned long long)nnz,
                 (unsigned long long)i_vec.size());
      return sp_mat(n_rows, n_cols);
    }

    // Construcción segura
    std::vector<uword> row_idx;
    std::vector<uword> col_idx;
    std::vector<double> values;

    row_idx.reserve(nnz);
    col_idx.reserve(nnz);
    values.reserve(nnz);

    const int n_cols_int = (int)n_cols;

    for (int col = 0; col < n_cols_int; col++) {
      int start = p_vec[col];
      int end = (col + 1 < n_cols_int) ? p_vec[col + 1] : (int)nnz;

      // Validar límites
      if (start < 0 || start > end || end > (int)nnz) {
        Rf_warning("Invalid index range at column %d: [%d, %d)", col, start, end);
        continue;
      }

      for (int j = start; j < end; j++) {
        if (j >= (int)i_vec.size()) break;

        int row_idx_val = i_vec[j];
        if (row_idx_val < 0 || row_idx_val >= (int)n_rows) {
          Rf_warning("Row index %d out of range at column %d", row_idx_val, col);
          continue;
        }

        double val = x_vec[j];
        if (!Rcpp::NumericVector::is_na(val) && val != 0.0) {
          row_idx.push_back((uword)row_idx_val);
          col_idx.push_back((uword)col);
          values.push_back(val);
        }
      }
    }

    const uword final_nnz = values.size();
    if (final_nnz == 0) {
      return sp_mat(n_rows, n_cols);
    }

    arma::umat locations(2, final_nnz);
    for (uword k = 0; k < final_nnz; k++) {
      locations(0, k) = row_idx[k];
      locations(1, k) = col_idx[k];
    }

    arma::vec vals(values);

    return sp_mat(locations, vals, n_rows, n_cols);

  } catch (std::exception& e) {
    Rf_error("Error in sparse matrix conversion: %s", e.what());
  } catch (...) {
    Rf_error("Unknown error in sparse matrix conversion");
  }

  return sp_mat(0, 0);
}

// Versión segura de multiplicación de matrices
sp_mat safe_matrix_multiply(const sp_mat& A, const sp_mat& B) {
  if (A.n_cols != B.n_rows) {
    Rf_error("Matrix dimensions incompatible for multiplication: %llu x %llu and %llu x %llu",
             (unsigned long long)A.n_rows, (unsigned long long)A.n_cols,
             (unsigned long long)B.n_rows, (unsigned long long)B.n_cols);
  }

  if (A.n_nonzero == 0 || B.n_nonzero == 0) {
    return sp_mat(A.n_rows, B.n_cols);
  }

  try {
    return A * B;
  } catch (std::exception& e) {
    Rf_error("Matrix multiplication failed: %s", e.what());
  }

  return sp_mat(A.n_rows, B.n_cols);
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

  try {
    // ── Conversión segura de inputs ─────────────────────────────────────────────
    sp_mat A_mat = safe_spmat_conversion(A);
    sp_mat M_mat = safe_spmat_conversion(M_orig);
    sp_mat g0    = safe_spmat_conversion(g0_input);

    // Verificar que las conversiones fueron exitosas
    if (!is_valid_sparse_matrix(A_mat)) {
      Rf_error("Failed to convert A matrix");
    }

    const uword n = A_mat.n_rows;

    // ── Validaciones ────────────────────────────────────────────────────────────
    if (n == 0) {
      Rf_error("set_A sparse matrix has zero rows");
    }

    if (A_mat.n_rows != A_mat.n_cols) {
      Rf_error("set_A must be a square matrix (n x n)");
    }

    if (M_mat.n_rows == 0 || M_mat.n_cols == 0) {
      Rf_error("set_M sparse matrix has zero dimensions");
    }

    if (M_mat.n_rows != n || M_mat.n_cols != n) {
      Rf_error("set_M dimensions (%llux%llu) must match set_A (%llux%llu)",
               (unsigned long long)M_mat.n_rows,
               (unsigned long long)M_mat.n_cols,
               (unsigned long long)n,
               (unsigned long long)n);
    }

    if (g0.n_rows != n || g0.n_cols != 1) {
      Rf_error("initial_points must be a %llux1 sparse vector, got %llux%llu",
               (unsigned long long)n,
               (unsigned long long)g0.n_rows,
               (unsigned long long)g0.n_cols);
    }

    if ((uword)suit_values.size() != n) {
      Rf_error("suit_values length (%llu) must match matrix dimension (%llu)",
               (unsigned long long)suit_values.size(),
               (unsigned long long)n);
    }

    if (nsteps <= 0) {
      Rf_error("nsteps must be a positive integer");
    }

    if (disper_prop < 0.0 || disper_prop > 1.0) {
      Rf_error("disper_prop must be a value between 0 and 1");
    }

    if (g0.n_nonzero == 0) {
      Rf_error("initial_points has no occupied cells. "
                 "Check that occurrence coordinates fall within the model extent.");
    }
    // ────────────────────────────────────────────────────────────────────────────

    // Semilla aleatoria (solo si es estocástico)
    if (stochastic_dispersal) {
      Rcpp::Environment baseEnv("package:base");
      Rcpp::Function setSeed = baseEnv["set.seed"];
      setSeed(123);  // Usar una semilla fija para reproducibilidad
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

    if (adj_list.size() > 0) {
      CharacterVector list_names = adj_list.names();

      for (R_xlen_t i = 0; i < adj_list.size(); i++) {
        if (adj_list[i] == R_NilValue) continue;

        // Obtener índice de celda fuente
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

        // Extraer vecinos
        try {
          DataFrame df = as<DataFrame>(adj_list[i]);
          if (!df.containsElementNamed("ToNonNaCell")) continue;

          NumericVector to_non_na = df["ToNonNaCell"];

          for (R_xlen_t j = 0; j < to_non_na.size(); j++) {
            double nb_val = to_non_na[j];
            if (Rcpp::NumericVector::is_na(nb_val)) continue;

            uword nb_index = static_cast<uword>(nb_val - 1);
            if (nb_index < n) {
              neighbors[src_cell].push_back(nb_index);
            }
          }
        } catch (...) {
          // Ignorar errores en la extracción de vecinos
          continue;
        }
      }
    }

    // Pre-asignación de salida
    List sdm(nsteps + 1);
    sdm[0] = g0;

    // Setup determinístico
    sp_mat AMA_base;
    if (!stochastic_dispersal) {
      sp_mat AM = safe_matrix_multiply(A_mat, M_mat);
      AMA_base = safe_matrix_multiply(AM, A_mat);
    }

    // Setup barra de progreso
    const int bar_width = 50;
    int progress_interval = std::max(1, nsteps / bar_width);
    if (progress_bar) {
      Rprintf("Simulation progress:\n[");
      for (int i = 0; i < bar_width; i++) Rprintf(" ");
      Rprintf("]\r[");
      R_FlushConsole();
    }

    // ── Loop principal de simulación ────────────────────────────────────────────
    sp_mat g0_next;

    for (int step = 1; step <= nsteps; ++step) {
      // Verificar interrupción del usuario cada 10 pasos
      if (step % 10 == 0) {
        Rcpp::checkUserInterrupt();
      }

      if (stochastic_dispersal) {
        // Simulación estocástica
        g0_next = g0;

        const uvec occ_locs = find(g0 > 0);
        const uword n_occ = occ_locs.n_elem;

        if (n_occ > 0) {
          vec rand_values = randu<vec>(n);

          for (uword i = 0; i < n_occ; ++i) {
            const uword cell = occ_locs(i);
            if (cell >= neighbors.size()) continue;

            const std::vector<uword>& cell_neighbors = neighbors[cell];

            for (uword j = 0; j < cell_neighbors.size(); ++j) {
              const uword nb_index = cell_neighbors[j];
              if (nb_index >= n) continue;
              if (binary_suit(nb_index) < 0.5) continue;

              double prob = disp_prop2_suitability ? suit_probs(nb_index) : disper_prop;
              if (prob > 0 && (prob >= 1.0 || rand_values(nb_index) < prob)) {
                g0_next(nb_index, 0) = 1.0;
              }
            }
          }
        }
      } else {
        // Dispersal determinístico con verificaciones
        if (AMA_base.n_rows > 0 && AMA_base.n_cols > 0) {
          sp_mat temp = safe_matrix_multiply(AMA_base, g0);
          g0_next = temp + g0;

          // Binarizar de forma segura
          for (uword i = 0; i < g0_next.n_rows; ++i) {
            if (g0_next(i, 0) > 0) {
              g0_next(i, 0) = 1.0;
            }
          }
        } else {
          g0_next = g0;
        }
      }

      g0 = g0_next;
      sdm[step] = g0;

      // Actualización barra de progreso
      if (progress_bar && (step % progress_interval == 0 || step == nsteps)) {
        int pos = static_cast<int>(bar_width * static_cast<double>(step) / nsteps);
        if (pos > bar_width) pos = bar_width;

        Rprintf("\r[");
        for (int i = 0; i < pos; i++) {
          Rprintf("=");
        }
        if (pos < bar_width) {
          Rprintf(">");
        }
        for (int i = pos + 1; i < bar_width; i++) {
          Rprintf(" ");
        }
        Rprintf("] %3d%%", static_cast<int>(100.0 * step / nsteps));
        R_FlushConsole();
      }
    }

    if (progress_bar) {
      Rprintf("\n");
      R_FlushConsole();
    }

    return sdm;

  } catch (std::exception& e) {
    Rf_error("Exception in sdm_sim_rcpp: %s", e.what());
  } catch (...) {
    Rf_error("Unknown exception in sdm_sim_rcpp");
  }

  return R_NilValue;
}
