
// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif


// user includes
#include <stdlib.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <vector>

// declarations

// [[Rcpp::export]] 
extern "C" {
  SEXP file193b1b67af14( SEXP x, SEXP y, SEXP N, SEXP w, SEXP p) ;
}

// definition

SEXP file193b1b67af14( SEXP x, SEXP y, SEXP N, SEXP w, SEXP p){
  BEGIN_RCPP
  
// using namespace Rcpp;
  RNGScope scope;
  
  NumericVector xa(x);
  NumericVector ya(y);
  int Na = as<int>(N);
  int n_xa = xa.size(), n_ya = ya.size();
  
  // arma::mat wa = Rcpp::as<arma::mat>(w);
  arma::sp_mat wa = Rcpp::as<arma::sp_mat>(w);
  
  Rcpp::NumericMatrix n(p); //one column to each x[i]/y
  NumericVector limite(n_xa);
  arma::mat z(n_ya, n_xa);
  arma::mat beta(n_ya, n_xa);
  arma::colvec one(xa.begin(), xa.size());
  arma::colvec betaM(ya.begin(), ya.size());
  arma::colvec betaP(ya.begin(), ya.size());
  arma::colvec Prob(ya.begin(), ya.size());
  arma::colvec betaC(ya.begin(), ya.size());
  z.col(0) = arma::zeros<arma::mat>(n_ya, 1);
  NumericVector nT(n_ya);
  //  NumericVector nV(n_ya);
  Rcpp::NumericMatrix nV(n_ya, n_xa);
  
  NumericVector init(n_xa);
  int i, j, randLR;
  
  for (i=0; i<n_xa; i++) {
    for (j=0; j<n_ya; j++) {
      z(j,i) = 0; 
    }
    randLR = rand() % n_ya;
    // randonly initialize z
    z(randLR,i)=1;
    init(i)=randLR+1;
    limite(i)  = i;
    one(i) = 1;
  }
  // the number of ways c can be produced
  
  beta = wa * z;  
  one = beta * one;
  
  int k, l,l2, pos;
  NumericMatrix prob_table(n_xa, Na);
  NumericMatrix class_table(n_xa, Na);
  
  NumericVector likelihood(Na);
 // system("mkdir beta");
  
  // normalize input Likelihood to (0,1)
  for(int m=0;m<n_xa;m++){
      NumericVector nV0(n_ya);
      nV0 = n(_, m);
      std::transform(nV0.begin(), nV0.end(), nV0.begin(), std::bind2nd(std::divides<double>(),std::accumulate(nV0.begin(), nV0.end(), 0.0)));
      nV(_,m )=nV0;
  }
  
  int oldval;
  for(k=0; k<Na; k++){          // foreach sample iteration
  
    // limite = f(limite, n_xa);
    std::random_shuffle(limite.begin(),limite.end());
    for (l=0; l<n_xa; l++) {             // foreach mass
      l2=limite(l);     
      // trick, where the z.col came from	
      betaM = one - wa * z.col(l2);
      
//      String file = "beta/beta_"+toString(k+1)+"_"+toString(l+1)+".txt";
//      betaM.save(file, arma::raw_ascii);
    
      // sum delta and normalize
      Prob = betaM;
      std::transform(Prob.begin(), Prob.end(), Prob.begin(), std::bind2nd(std::plus<double>(),1));
      std::transform(Prob.begin(), Prob.end(), Prob.begin(), std::bind2nd(std::divides<double>(),std::accumulate(Prob.begin(), Prob.end(), 0.0)));
      
     // multipling the the likelihoods
      int nP;
      for (nP=0; nP<betaP.size(); nP++) {
        betaP(nP) = Prob(nP)*nV(nP,l2);
        if(z(nP, l2)==1){
          oldval = nP;
        }  
        
      }
      
      // continuing control
      if (std::accumulate(betaP.begin(), betaP.end(), 0.0) > 0) {
        std::transform(betaP.begin(), betaP.end(), betaP.begin(), std::bind2nd(std::divides<double>(),std::accumulate(betaP.begin(), betaP.end(), 0.0)));
      }
      else {
        // std::transform(betaP.begin(), betaP.end(), betaP.begin(), std::bind2nd(std::plus<double>(),1));
      }
      // a way to sample according a distribution, based on rogers (2009)
      double base = ::Rf_runif(0,1);
      
      std::partial_sum (betaP.begin(), betaP.end(), betaC.begin());
      int sz = betaC.size();
      int idx;
      for (idx=0; idx<sz; idx++) {
        if(betaC[idx] > base) {
          pos = idx;
          break;
        }
      }
      
      if(pos!=oldval){
// redundant
//        one = betaM - wa * z.col(l2);  
// redundant

// short
//        int nP;
//        for (nP=0; nP<betaP.size(); nP++) {
//          z(nP, l2) = 0;
//        }
        z(oldval,l2) = 0; // others shoulde be 0 
// short

        z(pos, l2) = 1;
        one = betaM + wa * z.col(l2);
      }
            
      prob_table(l2,k) = betaP(pos);
      pos = pos + 1;
      class_table(l2,k) = pos;
    } //end for l
    
  
    // calculate Likelihood
    double log_sum=0.0;
    for(int m=0;m<n_xa;m++){
      
      betaM = one - wa * z.col(m);
      // sum delta and normalize
      Prob = betaM;
      std::transform(Prob.begin(), Prob.end(), Prob.begin(), std::bind2nd(std::plus<double>(),1));
      std::transform(Prob.begin(), Prob.end(), Prob.begin(), std::bind2nd(std::divides<double>(),std::accumulate(Prob.begin(), Prob.end(), 0.0)));
      
      double sum=0;
      for(int c=0;c<n_ya;c++){
        sum = sum + z(c,m) * nV(c,m) * Prob(c);
      }
      log_sum = log_sum +log(sum);
    }
    // save Likelihood
    likelihood(k)=log_sum;


} // end for k
  

  
  //		betaP = z.col(l2);
  return Rcpp::List::create( Rcpp::Named("beta")=beta, Rcpp::Named("z")=z, Rcpp::Named("log_likelihoods")=likelihood, Rcpp::Named("pos")=pos,  Rcpp::Named("n")=nV, Rcpp::Named("betaM")=betaM, Rcpp::Named("betaP")=betaP, Rcpp::Named("Prob")=Prob,Rcpp::Named("prob_table")=prob_table, Rcpp::Named("class_table")=class_table, Rcpp::Named("init")=init );
  
  END_RCPP
}
  