#include <Rcpp.h>
#include <chrono>
// #include <time.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector linear_interpolation_cpp(List x_r,
                                       NumericVector y,
                                       List xout_r,
                                       int rule,
                                       bool verbose) {

  // Pre-processing ------------------------------------------------------------
  // fetch n. of explanatory variables
  int n_var = x_r.size();
  Rprintf("n_var: %i\n", n_var);

  // reshape x and xout
  NumericVector x[n_var];
  NumericVector xout[n_var];
  for (int i = 0; i < n_var; i++) {
    x   [i] = x_r   [i];
    xout[i] = xout_r[i];
  }

  // fetch data size
  int n_ref = x   [0].size();
  int n_out = xout[0].size();
  Rprintf("n_ref: %i\n", n_ref);
  Rprintf("n_out: %i\n", n_out);

  // check if x and xout have reshaped correctly
  for (int i = 0; i < n_var; i++) {
    Rprintf("x[%i]: ", i);
    for (int j = 0; j < n_ref; j++) {
      Rprintf("%.1f ", x[i][j]);
    }
    Rprintf("\n");
  }
  for (int i = 0; i < n_var; i++) {
    Rprintf("xout[%i]: ", i);
    for (int j = 0; j < n_out; j++) {
      Rprintf("%.1f ", xout[i][j]);
    }
    Rprintf("\n");
  }

  // extract unique set of x
  NumericVector x_set[n_var];
  for (int i = 0; i < n_var; i++) {

    NumericVector x_i = x[i];
    std::sort(x_i.begin(), x_i.end());
    x_i.erase(std::unique(x_i.begin(), x_i.end()), x_i.end());
    x_set[i] = x_i;

    Rprintf("x_set[%i]: ", i);
    for (int j = 0; j < x_i.size(); j++) Rprintf("%.1f ", x_i[j]);
    Rprintf("\n");

  }


  // Interpolation -------------------------------------------------------------
  auto t1 = std::chrono::system_clock::now();
  Rprintf("Step 3/3: Preparation...     \n");

  float xout_i[n_var];
  float x_i_min, x_i_max, xout_ij;
  for (int i = 0; i < n_out; i++) {

    for (int j = 0; j < n_var; j++) xout_i[j] = xout[j][i];

    // When rule = 2, values of  xout is replaced with closest value of x
    // if they are out of range of x.
    if (rule == 2) {
      for (int j = 0; j < n_var; j++) {
      // for (int j = 0; j < 1; j++) {
        x_i_min = x_set[j][0];
        x_i_max = x_set[j][x_set[j].size() - 1];
        xout_ij = xout_i[j];
        xout_ij = (xout_ij > x_i_max)?x_i_max:xout_ij;
        xout_ij = (xout_ij < x_i_min)?x_i_min:xout_ij;
        xout_i[j] = xout_ij;
      }
    }

    // check if xout have converted correctly
    Rprintf("xout[%i] rule2 adj: ", i);
    for (int j = 0; j < n_var; j++) {
      Rprintf("%.1f ", xout_i[j]);
    }
    Rprintf("\n");

    // Where does each data of xout place? Which grid expanded by x?
    // An element of `smaller_indices` (i_rc in row r and column c) indicate
    // that the xout[r, c] is between i_rc_th and (i_rc + 1)_th values of x
    // (x_list[[c]][i_rc + (0:1)]).
    // smaller_indices <- mapply(function(x, vals) {
    //
    //   i <- sapply(vals, function(val) sum(val - x >= 0)) # faster
    //   // i <- sapply(vals, function(val) which(val - x < 0)[1] - 1) # slower
    //
    //   // When rule = 3, xout placed out of the region of x are treated
    //   // as if they are in the closest grid.
    //   pmax(pmin(i, length(x) - 1), 1)
    //
    // }, x = x_list, vals = xout)


  }




  // 時刻の表示
  // std::time_t now_c = std::chrono::system_clock::to_time_t(t1);
  // std::cout << "現在時刻: " << std::ctime(&now_c);



  NumericVector yout(n_out);
  for (int i = 0; i < n_out; i++) {
    yout[i] = i;
  }
  return yout;


//   // Interpolation (1): Preparation --------------------------------------------
//   t1 <- proc.time()[3]
//   if (verbose) cat('Step 3/4: Preparation...     ')
//

}
