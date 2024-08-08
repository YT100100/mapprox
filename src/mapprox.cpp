#include <Rcpp.h>
#include <chrono>
#include <cmath>
// #include <time.h>
using namespace Rcpp;

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
  int n_ref = x   [0].size(); // size of reference data
  int n_out = xout[0].size(); // size of output data
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
  // Rprintf("Step 3/3: Preparation...     \n");

  NumericVector yout(n_out, NA_REAL);
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
          smaller_ind_i[j] = (k == 0)?0:k - 1;
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

      // one of the grid points surrounding xout_i
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

      // x of the focusing grid point
      double x_k[n_var];
      for (int j = 0; j < n_var; j++) x_k[j] = x_set[j][ind_k[j]];

      // check
      // Rprintf("xk ");
      // for (int j = 0; j < n_var; j++) Rprintf("%.1f ", x_k[j]);

      // y of the focusing grid point
      int l = 0;
      while (true) {
        bool go_next = false;
        for (int j = 0; j < n_var; j++) {
          if (x[j][l] != x_k[j]) go_next = true;
          if (go_next) break;
        }
        if (!go_next) break;
        l++;
      }
      double y_k = y[l];

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

      // check
      // Rprintf("w %.1f ", weight_k);
      // Rprintf("\n");

    }

    yout[i] = yout_i;
    Rcpp::checkUserInterrupt();

  }

  // 時刻の表示
  // std::time_t now_c = std::chrono::system_clock::to_time_t(t1);
  // std::cout << "現在時刻: " << std::ctime(&now_c);

  return yout;

}
