#include <Rcpp.h>
#include <unistd.h>
#include <chrono>
#include <cmath>
using namespace Rcpp;

//' Main function of mapprox() in R
//'
//' \code{mapprox} is an wrapper of \code{linear_interpolation} and
//' \code{linear_interpolation_cpp}.
//'
//' This is the main function of \code{mapprox(..., use_cpp = TRUE)}.
//' Provided data are cleaned beforehand in \code{mapprox} and
//' then passed to this function.
//' \code{x} must be completely gridded data,
//' and duplicated values are not allowed.
//' Also, no NA in \code{x}, \code{y}, and \code{xout} is allowed.
//'
//' @param x Reference data of explanatory variables.
//'   This must be a data frame in which elements are all numeric.
//'   Also, this should be ordered by \code{order} function.
//' @param y Reference data of a target variable.
//'   This must be a numeric vector with length \code{nrow(x)}.
//' @param xout Values of explanatory variables where interpolation is to
//'   take place. This must be a data frame and
//'   \code{ncol(xout)} must be the same as \code{ncol(x)}.
//' @param rule A scalar integer determining how interpolation outside the
//'   given explanatory variables works. This must be eigher
//'   \code{1}, \code{2}, or \code{3}.
//'   If \code{rule = 1}, \code{NA}s are returned for such points.
//'   If \code{rule = 2}, the values at the closest points in the region
//'   given as the reference data are returned.
//'   If \code{rule = 3}, linear extension of the interpolated result is applied.
//' @param verbose Logical. Should progress of the process printed?
//'
//' @return A numeric vector of interpolated result.
//'
//' @seealso [mapprox()]
// [[Rcpp::export]]
NumericVector linear_interpolation_cpp(List x_r,
                                       NumericVector y,
                                       List xout_r,
                                       int rule,
                                       bool verbose) {

  // Pre-processing ------------------------------------------------------------

  // convert x and xout to NumericVector
  int n_var = x_r.size();    // n. of explanatory variables
  NumericVector x   [n_var];
  NumericVector xout[n_var];
  for (int j = 0; j < n_var; j++) { // j not i for consistency with latter part
    x   [j] = x_r   [j];
    xout[j] = xout_r[j];
  }
  int n_out = xout[0].size(); // size of output data
  // int n_ref = x   [0].size(); // size of reference data
  // Rprintf("n_var: %d\n", n_var); // check
  // Rprintf("n_ref: %d\n", n_ref); // check
  // Rprintf("n_out: %d\n", n_out); // check

  // check if x and xout have reshaped correctly
  // for (int i = 0; i < n_var; i++) {
  //   Rprintf("x[%d]: ", i);
  //   for (int j = 0; j < n_ref; j++) {
  //     Rprintf("%.1f ", x[i][j]);
  //   }
  //   Rprintf("\n");
  // }
  // for (int i = 0; i < n_var; i++) {
  //   Rprintf("xout[%d]: ", i);
  //   for (int j = 0; j < n_out; j++) {
  //     Rprintf("%.1f ", xout[i][j]);
  //   }
  //   Rprintf("\n");
  // }

  // extract unique set of x
  NumericVector x_set[n_var];
  for (int j = 0; j < n_var; j++) {
    NumericVector x_j = clone(x[j]);
    std::sort(x_j.begin(), x_j.end());
    x_j.erase(std::unique(x_j.begin(), x_j.end()), x_j.end());
    x_set[j] = x_j;
    // check
    // Rprintf("x_set[%d]: ", j);
    // for (int i = 0; i < x_j.size(); i++) Rprintf("%.1f ", x_j[i]);
    // Rprintf("\n");
  }

  // length of grids
  NumericVector grid_len[n_var];
  for (int j = 0; j < n_var; j++) {
    int grid_len_size_j = x_set[j].size() - 1;
    NumericVector grid_len_j(grid_len_size_j);
    for (int i = 0; i < grid_len_size_j; i++) {
      grid_len_j[i] = x_set[j][i + 1] - x_set[j][i];
    }
    grid_len[j] = grid_len_j;
  }
  // check
  // for (int j = 0; j < n_var; j++) {
  //   Rprintf("grid_len[%d]: ", j);
  //   for (int i = 0; i < grid_len[j].size(); i++) Rprintf("%.1f ", grid_len[j][i]);
  //   Rprintf("\n");
  // }
  // Rprintf("\n");


  // Interpolation -------------------------------------------------------------
  auto t1 = std::chrono::system_clock::now();

  NumericVector yout(n_out, NA_REAL);
  bool print_progress = false;
  for (int i = 0; i < n_out; i++) {

    // i_th data point for output
    double xout_i[n_var];
    for (int j = 0; j < n_var; j++) xout_i[j] = xout[j][i];

    // check
    // Rprintf("%d: x1 ", i);
    // for (int j = 0; j < n_var; j++) Rprintf("%.1f ", xout_i[j]);

    // If rule = 1, the returned values should be NA
    // for xout placed out of range of x.
    if (rule == 1) {
      bool is_out_i = false;
      for (int j = 0; j < n_var; j++) {
        int maxind = x_set[j].size() - 1;
        is_out_i = is_out_i || xout_i[j] < x_set[j][0];
        is_out_i = is_out_i || xout_i[j] > x_set[j][maxind];
      }
      if (is_out_i) continue;
    }

    // When rule = 2, values of xout is replaced with closest value of x
    // if they are out of range of x.
    if (rule == 2) {
      for (int j = 0; j < n_var; j++) {
        int maxind = x_set[j].size() - 1;
        double xout_ij = xout_i[j];
        xout_ij = (xout_ij > x_set[j][maxind])?x_set[j][maxind]:xout_ij;
        xout_ij = (xout_ij < x_set[j][0]     )?x_set[j][0]     :xout_ij;
        xout_i[j] = xout_ij;
      }
    }

    // check
    // Rprintf("x2 ");
    // for (int j = 0; j < n_var; j++) Rprintf("%.1f ", xout_i[j]);

    // Where does each data of xout place? Which grid expanded by x?
    // An element of `smaller_index` (k) indicate
    // that the xout_i[j] is between k_th and (k + 1)_th values of x
    // (x_set[j][k] and x_set[j][k + 1]).
    int smaller_ind_i[n_var];
    double grid_len_i[n_var]; // lengths of grids where the data point placed
    for (int j = 0; j < n_var; j++) {
      int n_val = x_set[j].size();
      int k = 0;
      while (true) {
        if (xout_i[j] < x_set[j][k]) {
          smaller_ind_i[j] = (k == 0)?0:(k - 1);
          break;
        }
        if (k == n_val - 2) {
          smaller_ind_i[j] = k;
          break;
        }
        k++;
      }
      grid_len_i[j] = grid_len[j][smaller_ind_i[j]];
    }

    // check
    // Rprintf("si ");
    // for (int j = 0; j < n_var; j++) Rprintf("%d ", smaller_ind_i[j]);
    // Rprintf("gl ");
    // for (int j = 0; j < n_var; j++) Rprintf("%.1f ", grid_len_i[j]);
    // Rprintf("\n");

    // All 0/1 patterns corresponding to each explanatory variable
    double yout_i = 0;
    for (int k = 0; k < std::pow(2, n_var); k++) {

      // index of k_th grid point surrounding xout_i
      int added_ind_k[n_var];
      int ind_k[n_var];
      int k_rest = k;
      for (int j = 0; j < n_var - 1; j++) {
        int divider_j = std::pow(2, n_var - j - 1);
        added_ind_k[j] = k_rest / divider_j;
        ind_k[j] = smaller_ind_i[j] + added_ind_k[j];
        k_rest = k_rest % divider_j;
      }
      added_ind_k[n_var - 1] = k_rest;
      ind_k[n_var - 1] = smaller_ind_i[n_var - 1] + k_rest;

      // check
      // Rprintf("  %d: ai ", k);
      // for (int j = 0; j < n_var; j++) Rprintf("%d ", added_ind_k[j]);
      // Rprintf(" k ", k);
      // for (int j = 0; j < n_var; j++) Rprintf("%d ", ind_k[j]);

      // x of the k_th grid point
      // double x_k[n_var];
      // for (int j = 0; j < n_var; j++) x_k[j] = x_set[j][ind_k[j]];

      // check
      // Rprintf("xk ");
      // for (int j = 0; j < n_var; j++) Rprintf("%.1f ", x_k[j]);

      // y of the k_th grid point
      int ind_y_k = 0;
      for (int j = 0; j < n_var; j++) {
        ind_y_k += ind_k[j];
        if (j == n_var - 1) break;
        ind_y_k *= x_set[j + 1].size();
      }
      double y_k = y[ind_y_k];

      // check
      // Rprintf("l %d yk %.1f ", l, y_k);

      // weight of y_k
      double weight_k = 1;
      for (int l = 0; l < n_var; l++) {
        double x_nei_l = x_set[l][smaller_ind_i[l] - added_ind_k[l] + 1];
        double weight_kl = (xout_i[l] - x_nei_l) / grid_len_i[l];
        if (added_ind_k[l] == 0) weight_kl *= -1;
        weight_k *= weight_kl;
      }
      yout_i += y_k * weight_k;

      Rcpp::checkUserInterrupt();

      // check
      // Rprintf("w %.1f ", weight_k);
      // Rprintf("\n");

    }

    // show rest time
    auto t2 = std::chrono::system_clock::now();
    double t_elap = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    double t_rest = t_elap / (i + 1) * (n_out - i - 1);
    if (!print_progress && verbose && t_elap > 10) {
      print_progress = true;
      Rprintf("\n");
    }
    if (print_progress && i % 1000 == 0) {
      Rprintf("\r  Loop %i/%i: rest time %.0fs", i + 1, n_out, t_rest);
    }

    yout[i] = yout_i;

  }

  if (print_progress) Rprintf("\n                              ");
  // if (print_progress) Rprintf("\n");
  return yout;

}
