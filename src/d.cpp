#include <RcppArmadillo.h>
#include <iostream>
#include <utility>
// [[Rcpp::depends(RcppArmadillo)]]

typedef struct {
  double mean;
  double se;
  double z;
  uint n;
  uint n_block;
} JKStat ;

namespace Rcpp {
  template <>
  SEXP wrap(const JKStat& x) {
    return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("mean") = x.mean,
      Rcpp::Named("se") = x.se,
      Rcpp::Named("z") = x.z,
      Rcpp::Named("n") = x.n,
      Rcpp::Named("n_block") = x.n_block));
  };
}

arma::mat blocked_col_sums(const arma::mat m, uint n, uint s) {
  const uint limit = std::min(m.n_rows, s * n);
  
  arma::rowvec aux;
  arma::mat blocked(n, m.n_cols, arma::fill::zeros);
  for (uint i = 0; i < limit; i++) {
    aux = m.row(i);
    aux(arma::find_nonfinite(aux)).zeros();
    blocked.row(i/s) += aux;
  }
  return(blocked);
}

arma::vec blocked_sums(const arma::vec v, uint n, uint s) {
  const uint limit = std::min(v.size(), s * n);
  
  arma::vec aux = arma::vec(n, arma::fill::zeros);
  for (uint i = 0; i < limit; i++) {
    if (!arma::is_finite(v(i)))
      continue;
    aux(i/s) += v(i);
  }
  return(aux);
}

//[[Rcpp::export]]
arma::mat d_statistic(const arma::mat& fq, const arma::uvec cols) {
  arma::mat result = arma::mat(fq.n_rows, 2);
  arma::vec w, x, y, z;
  w = fq.col(cols(0) - 1);
  x = fq.col(cols(1) - 1);
  y = fq.col(cols(2) - 1);
  z = fq.col(cols(3) - 1);
  result.col(0) = (w - x) % (y - z);
  result.col(1) = (w + x - 2 * (w % x)) % (y + z - 2 * (y % z));
  return(result);
}

JKStat jackknife_d_(
    const arma::mat& fq, 
    const arma::uvec cols, 
    uint size) {
  uint n = ceil((double) fq.n_rows / (double) size);
  double val, mn, se;
  arma::vec jkv, snv, wei, hj, jv;
  
  arma::mat dstat = d_statistic(fq, cols);
  arma::mat bstat = blocked_col_sums(dstat, n, size);
  
  val = arma::sum(bstat.col(0)) / arma::sum(bstat.col(1));
  jkv = (arma::sum(bstat.col(0)) - bstat.col(0)) /
    (arma::sum(bstat.col(1)) - bstat.col(1));
  snv = arma::vec(dstat.n_rows, arma::fill::zeros);
  snv.elem(arma::find_finite(dstat.col(0))).fill(1);
  snv = blocked_sums(snv, n, size);
  wei = 1 - snv / arma::sum(snv);
  mn = n * val - arma::sum(wei % jkv);
  hj = arma::sum(snv) / snv;
  jv = arma::pow((hj * val - (hj - 1) % jkv) - mn, 2) / (hj - 1);
  se = sqrt(arma::mean(jv.elem(arma::find_finite(jv))));

  JKStat result = { mn, se, mn/se, (uint) arma::sum(snv), n };
  return result;
}

JKStat jackknife_f4_(
    const arma::mat& fq, 
    const arma::uvec cols, 
    uint size) {
  uint n = ceil((double) fq.n_rows / (double) size);
  double val, mn, se;
  arma::vec jkv, snv, wei, hj, jv;
  
  arma::mat dstat = d_statistic(fq, cols);
  arma::mat bstat = blocked_col_sums(dstat, n, size);
  
  val = arma::sum(bstat.col(0));
  jkv = (arma::sum(bstat.col(0)) - bstat.col(0));
  snv = arma::vec(fq.n_rows, arma::fill::zeros);
  snv.elem(arma::find_finite(dstat.col(0))).fill(1);
  snv = blocked_sums(snv, n, size);
  wei = 1 - snv / arma::sum(snv);
  mn = n * val - arma::sum(wei % jkv);
  hj = arma::sum(snv) / snv;
  jv = arma::pow((hj * val - (hj - 1) % jkv) - mn, 2) / (hj - 1);
  se = sqrt(arma::mean(jv.elem(arma::find_finite(jv))));

  JKStat result = {mn, se, mn / se, (uint) arma::sum(snv), n};
  return result;
}

//[[Rcpp::export]]
JKStat jackknife_d(
    const arma::mat& fq, 
    const arma::uvec cols, 
    uint size) {
  return jackknife_d_(fq, cols, size);
}

//[[Rcpp::export]]
Rcpp::List jackknife_ds(
    const arma::mat& fq, 
    const arma::umat& cols, 
    uint size) {
  uint last_check = 0;
  Rcpp::List result(cols.n_rows);
  for (uint i = 0; i < cols.n_rows; i++) {
    result[i] = jackknife_d_(fq, arma::trans(cols.row(i)), size);
    if (i * 80 / cols.n_rows >= last_check) {
      Rcpp::checkUserInterrupt();
      std::cout << ".";
      std::cout.flush();
      last_check++;
    }
  }
  std::cout << std::endl;
  return result;
}

//[[Rcpp::export]]
Rcpp::List jackknife_fs(
    const arma::mat& fq, 
    const arma::umat& cols, 
    uint size) {
  uint last_check = 0;
  Rcpp::List result(cols.n_rows);
  for (uint i = 0; i < cols.n_rows; i++) {
    result[i] = jackknife_f4_(fq, arma::trans(cols.row(i)), size);
    if (i * 80 / cols.n_rows >= last_check) {
      Rcpp::checkUserInterrupt();
      std::cout << ".";
      std::cout.flush();
      last_check++;
    }
  }
  std::cout << std::endl;
  return result;
}

JKStat jackknife_f3_(
  const arma::mat& ac,
  const arma::mat& an,
  const arma::uvec cols,
  uint size) {
  // Compute hc correction
  arma::vec c0, n0, hc;
  c0 = ac.col(cols(0)-1);
  n0 = an.col(cols(0)-1);
  hc = c0 % (n0-c0) / (arma::pow(n0, 3) - arma::pow(n0, 2));
  // Compute F3
  arma::vec f3, fi;
  arma::mat fq = ac.cols(cols-1) / an.cols(cols-1);
  f3 = (fq.col(0) - fq.col(1)) % (fq.col(0) - fq.col(2)) - hc;
  fi = arma::vec(fq.n_rows, arma::fill::zeros);
  fi.elem(arma::find_finite(f3)).ones();
  // Block computation
  uint n = arma::sum(fi);
  uint nb = ceil((double) fq.n_rows / (double) size);
  f3 = blocked_sums(f3, nb, size);
  fi = blocked_sums(fi, nb, size);
  // Jackknife computation
  arma::vec jkv, wei, jh, jv;
  double val, mn, se;
  jkv = (arma::sum(f3) - f3) / (arma::sum(fi) - fi);

  val = arma::sum(f3) / arma::sum(fi);
  wei = 1 - fi / arma::sum(fi);
  mn  = nb * val - arma::sum(wei % jkv);

  jh = arma::sum(fi) / fi;
  jv = arma::pow((jh * val - (jh-1) % jkv) - mn, 2)/(jh - 1);
  se = sqrt(arma::mean(jv.elem(arma::find_finite(jv))));
  // Return results
  JKStat result = {mn, se, mn / se, n, nb};
  return result;
}

//[[Rcpp::export]]
Rcpp::List jackknife_f3s(
    const arma::mat& ac, 
    const arma::mat& an, 
    const arma::umat& cols, 
    uint size) {
  uint last_check = 0;
  Rcpp::List result(cols.n_rows);
  for (uint i = 0; i < cols.n_rows; i++) {
    result[i] = jackknife_f3_(ac, an, arma::trans(cols.row(i)), size);
    if (i * 80 / cols.n_rows >= last_check) {
      Rcpp::checkUserInterrupt();
      std::cout << ".";
      std::cout.flush();
      last_check++;
    }
  }
  std::cout << std::endl;
  return result;
}
