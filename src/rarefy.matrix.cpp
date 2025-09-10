#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rarefy_matrix(NumericMatrix mat, int N) {
  int n_samples = mat.nrow();
  int n_taxa = mat.ncol();
  NumericMatrix result(n_samples, n_taxa);

  for (int i = 0; i < n_samples; i++) {
    double row_sum_d = 0;
    for (int j = 0; j < n_taxa; j++) {
      row_sum_d += mat(i, j);
    }

    int row_sum = static_cast<int>(row_sum_d);

    if (row_sum < N) {
      // 深度不足，整行设为 NA
      for (int j = 0; j < n_taxa; j++) {
        result(i, j) = NA_REAL;
      }
    } else {
      int remainingN = N;
      int remainingTotal = row_sum;

      for (int j = 0; j < n_taxa; j++) {
        int count = static_cast<int>(mat(i, j));
        if (count == 0) {
          result(i, j) = 0;
        } else {
          // Rf_rhyper(white, black, draws) 返回 double
          int sampled = static_cast<int>(Rf_rhyper(count, remainingTotal - count, remainingN));
          result(i, j) = sampled;
          remainingN -= sampled;
          remainingTotal -= count;
        }
      }
    }
  }

  return result;
}
