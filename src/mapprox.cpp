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


  // Interpolation (1): Preparation --------------------------------------------
  auto t1 = std::chrono::system_clock::now();
  Rprintf("Step 3/4: Preparation...     \n");

  // When rule = 2, values of  xout is replaced with closest value of x
  // if they are out of range of x.
  if (rule == 2) {

    float x_i_min, x_i_max, xout_ij;
    for (int i = 0; i < n_var; i++) {

      x_i_min = x_set[i][0];
      x_i_max = x_set[i][x_set[i].size() - 1];

      for (int j = 0; j < n_out; j++) {
        xout_ij = xout[i][j];
        xout_ij = (xout_ij > x_i_max)?x_i_max:xout_ij;
        xout_ij = (xout_ij < x_i_min)?x_i_min:xout_ij;
        xout[i][j] = xout_ij;
      }

    }

    // check if xout have converted correctly
    Rprintf("Processing of data out of region (rule 2)\n");
    for (int i = 0; i < n_var; i++) {
      Rprintf("xout[%i]: ", i);
      for (int j = 0; j < n_out; j++) {
        Rprintf("%.1f ", xout[i][j]);
      }
      Rprintf("\n");
    }

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
