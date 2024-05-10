#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//double non0_min(vec &p);
vec my_pava(vec &values, vec &weight, bool decreasing);
vec colfdr(vec &xis, vec &f1, vec &f2, vec &f3);
vec fdr(vec &lfdr);

// [[Rcpp::export]]
SEXP em_pava(SEXP pa_in, SEXP pb_in, SEXP x_in, SEXP pi0a_in, SEXP pi0b_in, SEXP pi0x_in) {
  try{
    const int maxIter = 200;
    const double tol = 1e-3;
    const double pvalue_cutoff = 1e-15;
    const double f_cutoff = 1e-15;

    vec pa = as<arma::vec>(pa_in), pb = as<arma::vec>(pb_in), x = as<arma::vec>(x_in);
    const double pi0_pa = Rcpp::as<double>(pi0a_in), pi0_pb = Rcpp::as<double>(pi0b_in), pi0_x = Rcpp::as<double>(pi0x_in);

    int m = pa.size();

    double min_a = pa.elem(find(pa>0)).min(), min_b = pb.elem(find(pb>0)).min(), min_x = x.elem(find(x>0)).min();
    pa.elem(find(pa==0)).fill(pvalue_cutoff < min_a ? pvalue_cutoff : min_a);
    pb.elem(find(pb==0)).fill(pvalue_cutoff < min_b ? pvalue_cutoff : min_b);
    x.elem(find(x==0)).fill(pvalue_cutoff < min_x ? pvalue_cutoff : min_x);

    vec p1 = sort(pa), p2 = sort(pb), p3 = sort(x);
    uvec ix1 = sort_index(pa), ix2 = sort_index(pb), ix3 = sort_index(x);

    vec p1_diff(m), p2_diff(m), p3_diff(m);
    p1_diff(0) = p1(0);
    p2_diff(0) = p2(0);
    p3_diff(0) = p3(0);
    for(int i = 1; i<m; ++i){
      p1_diff(i) = p1(i) - p1(i-1);
      p2_diff(i) = p2(i) - p2(i-1);
      p3_diff(i) = p3(i) - p3(i-1);
    }

    // Initialization
    double xi000=pi0_pa*pi0_pb*pi0_x, xi001=pi0_pa*pi0_pb*(1-pi0_x), xi010=pi0_pa*(1-pi0_pb)*pi0_x, xi011=pi0_pa*(1-pi0_pb)*(1-pi0_x);
    double xi100=(1-pi0_pa)*pi0_pb*pi0_x, xi101=(1-pi0_pa)*pi0_pb*(1-pi0_x), xi110=(1-pi0_pa)*(1-pi0_pb)*pi0_x, xi111=(1-pi0_pa)*(1-pi0_pb)*(1-pi0_x);

    vec f0 = ones(m,1), f1(m), f2(m), f3(m);
    f1 = 1 - pa;
    f2 = 1 - pb;
    f3 = 1 - x;

    vec loglik(maxIter);
    loglik(0) = -datum::inf;

    vec f(m), f000(m), f001(m), f010(m), f011(m), f100(m), f101(m), f110(m), f111(m);
    vec gamma000(m), gamma001(m), gamma010(m), gamma011(m), gamma100(m), gamma101(m), gamma110(m), gamma111(m);
    vec Q1(m), Q2(m), Q3(m), y1(m), y2(m), y3(m), res1(m), res2(m), res3(m);
    double loglik_delta;

    // std::cout<<"EM begins:"<<std::endl;

    for (int i = 1; i < maxIter; i++){
      // E-step
      f000 = xi000 * f0 % f0 % f0;
      f001 = xi001 * f0 % f0 % f3;
      f010 = xi010 * f0 % f2 % f0;
      f011 = xi011 * f0 % f2 % f3;
      f100 = xi100 * f1 % f0 % f0;
      f101 = xi101 * f1 % f0 % f3;
      f110 = xi110 * f1 % f2 % f0;
      f111 = xi111 * f1 % f2 % f3;
      f = f000 + f001 + f010 + f011 + f100 + f101 + f110 + f111;
      gamma000 = f000 / f;
      gamma001 = f001 / f;
      gamma010 = f010 / f;
      gamma011 = f011 / f;
      gamma100 = f100 / f;
      gamma101 = f101 / f;
      gamma110 = f110 / f;
      gamma111 = f111 / f;

      // M-step
      // update f1 and f2
      Q1 = gamma100 + gamma101 + gamma110 + gamma111;
      Q2 = gamma010 + gamma011 + gamma110 + gamma111;
      Q3 = gamma001 + gamma011 + gamma101 + gamma111;
      Q1 = Q1(ix1);
      Q2 = Q2(ix2);
      Q3 = Q3(ix3);

      y1 = - p1_diff * sum(Q1) / Q1;
      y2 = - p2_diff * sum(Q2) / Q2;
      y3 = - p3_diff * sum(Q3) / Q3;

      y1.elem(find_nonfinite(y1)).fill(y1.elem(find_finite(y1)).min());
      y2.elem(find_nonfinite(y2)).fill(y2.elem(find_finite(y2)).min());
      y3.elem(find_nonfinite(y3)).fill(y3.elem(find_finite(y3)).min());

      res1 = my_pava(y1, Q1, true);
      res2 = my_pava(y2, Q2, true);
      res3 = my_pava(y3, Q3, true);

      f1 = -1 / res1;
      f1 = f1 / sum(f1 % p1_diff);
      f1(ix1) = f1;
      f1.elem(find_nan(f1)).fill(f1.min());

      f2 = -1 / res2;
      f2 = f2 / sum(f2 % p2_diff);
      f2(ix2) = f2;
      f2.elem(find_nan(f2)).fill(f2.min());

      f3= -1 / res3;
      f3 = f3 / sum(f3 % p3_diff);
      f3(ix3) = f3;
      f3.elem(find_nan(f3)).fill(f3.min());

      double min_f1 = f1.elem(find(f1>0)).min(), min_f2 = f2.elem(find(f2>0)).min(), min_f3 = f3.elem(find(f3>0)).min();
      f1.elem(find(f1<=0)).fill(f_cutoff < min_f1 ? f_cutoff : min_f1);
      f2.elem(find(f2<=0)).fill(f_cutoff < min_f2 ? f_cutoff : min_f2);
      f3.elem(find(f3<=0)).fill(f_cutoff < min_f3 ? f_cutoff : min_f3);

      // update xi's
      xi000 = mean(gamma000);
      xi001 = mean(gamma001);
      xi010 = mean(gamma010);
      xi011 = mean(gamma011);
      xi100 = mean(gamma100);
      xi101 = mean(gamma101);
      xi110 = mean(gamma110);
      xi111 = mean(gamma111);

      // calculate the updated log-likelihood
      loglik(i) = sum(Q1%log(f1)+Q2%log(f2)+Q3%log(f3))+sum(gamma000*log(xi000)
                  +gamma001*log(xi001)+gamma010*log(xi010)+gamma011*log(xi011)
                  +gamma100*log(xi100)+gamma101*log(xi101)+gamma110*log(xi110)
                  +gamma111*log(xi111));
      loglik_delta = abs((loglik(i) - loglik(i-1))/loglik(i-1));

      // std::cout<<i<<". "<< loglik(i) << ", delta = " << loglik_delta << std::endl;

      if(loglik_delta <= tol){
        break;
      }
    }
    vec xis = {xi000, xi001, xi010, xi011, xi100, xi101, xi110, xi111};
    vec Lfdr = colfdr(xis, f1, f2, f3);
    vec radj = fdr(Lfdr);

    return Rcpp::List::create(Rcpp::Named("Lfdr") = Lfdr,
                              Rcpp::Named("radj") = radj,
                              Rcpp::Named("loglik") = loglik,
                              Rcpp::Named("xi") = xis,
                              Rcpp::Named("f1") = f1,
                              Rcpp::Named("f2") = f2,
                              Rcpp::Named("f3") = f3);
  } catch( std::exception &ex ) {
    forward_exception_to_r(ex);
    //std::cout<< "aaaaa=========================" << std::endl;
    return Rcpp::List::create();
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
    //std::cout<< "bbbbb========================" << std::endl;
    return Rcpp::List::create();
  }
}

// double non0_min(vec &p){
//   double _min = std::numeric_limits<double>::max();
//
//   p.for_each([&_min](double &val) { if(val > 0 && val < _min) _min = val; });
//
//   return _min;
// }

vec colfdr(vec &xis, vec &f1, vec &f2, vec &f3){
  int m = f1.size();
  vec f0 = ones(m,1);

  vec f000 = xis(0) * f0 % f0 % f0;
  vec f001 = xis(1) * f0 % f0 % f3;
  vec f010 = xis(2) * f0 % f2 % f0;
  vec f011 = xis(3) * f0 % f2 % f3;
  vec f100 = xis(4) * f1 % f0 % f0;
  vec f101 = xis(5) * f1 % f0 % f3;
  vec f110 = xis(6) * f1 % f2 % f0;
  vec f111 = xis(7) * f1 % f2 % f3;
  vec f = f000 + f001 + f010 + f011 + f100 + f101 + f110 + f111;

  vec Lfdr = (f000+f001+f010+f011+f100+f101)/f;

  return Lfdr;
}

vec fdr(vec &lfdr){
  int m = lfdr.size();

  vec ordered_lfdr = sort(lfdr), s = linspace(1,m,m);
  uvec ix_lfdr = sort_index(lfdr);

  vec radj = cumsum(ordered_lfdr)/s;
  radj(ix_lfdr) = radj;

  return radj;
}

// inline bool compare(double x, double y, Ordering ordering)
// {
//   return ordering == Ordering::Increasing ? x > y : x < y;
// }

vec my_pava(vec &values, vec &weight, bool decreasing)
{
  if(decreasing){
    values = reverse(values);
    weight = reverse(weight);
  }
  vec w(values.size(), fill::zeros);
  vec x(values.size(), fill::zeros);
  x[0] = values[0];
  w[0] = weight[0];
  unsigned j = 0;
  vec s(values.size(), fill::zeros);

  for (unsigned i = 1; i < values.size(); i++) {
    j += 1;
    x[j] = values[i];
    w[j] = weight[i];
    while (j > 0 && x[j - 1]>x[j]) {
      x[j - 1] = (w[j] * x[j] + w[j - 1] * x[j - 1]) / (w[j] + w[j - 1]);
      w[j - 1] += w[j];
      j -= 1;
    }
    s[j + 1] = i + 1;
  }

  vec ww(values.size(), fill::zeros);
  vec xx(values.size(), fill::zeros);
  for (unsigned k = 0; k < j + 1; k++) {
    for (unsigned i = s[k]; i < s[k + 1]; i++) {
      ww[i] = w[k];
      xx[i] = x[k];
    }
  }

  if(decreasing){
    xx = reverse(xx);
  }

  return xx;
}
