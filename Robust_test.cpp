// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include<math.h>
#include <cmath> 
#include <algorithm>  
//#include "kendallc.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

// Function returns the rank vector of the set of observations
vec rankify(vec X) {
  int N = X.n_elem;
  vec Rank_X(N);// Rank Vector
  
  for(int i = 0; i < N; i++)
  {
    int r = 1, s = 1;
    // Count no of smaller elements
    // in 0 to i-1
    for(int j = 0; j < i; j++) {
      if (X[j] < X[i] ) r++;
      if (X[j] == X[i] ) s++;
    }
    // Count no of smaller elements
    // in i+1 to N-1
    for (int j = i+1; j < N; j++) {
      if (X[j] < X[i] ) r++;
      if (X[j] == X[i] ) s++;
    }
    // Use Fractional Rank formula fractional_rank = r + (n-1)/2
    Rank_X[i] = r + (s-1) * 0.5;       
  }
  // Return Rank Vector
  return Rank_X;
}

// [[Rcpp::export]]
// function that returns Pearson correlation coefficient.
float spearman(vec X,vec Y)
{
  X = rankify(X);
  Y = rankify(Y);
  int n =  X.n_elem;
  float sum_X = 0, sum_Y = 0,sum_XY = 0,squareSum_X = 0,squareSum_Y = 0;
  
  for (int i = 0; i < n; i++)
  {
    // sum of elements of array X.
    sum_X = sum_X + X[i];
    // sum of elements of array Y.
    sum_Y = sum_Y + Y[i];
    // sum of X[i] * Y[i].
    sum_XY = sum_XY + X[i] * Y[i];
    // sum of square of array elements.
    squareSum_X = squareSum_X + X[i] * X[i];
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
  }
  
  // use formula for calculating correlation coefficient.
  float corr = (float)(n * sum_XY -sum_X * sum_Y) /
                  sqrt((n * squareSum_X -sum_X * sum_X) *
                    (n * squareSum_Y -sum_Y * sum_Y));
  
  return corr;
}



// [[Rcpp::export]]
mat Spear_mat(mat Z){
  int m=Z.n_cols;
  int i,j;
  mat corr(m, m, fill::none);
  for(i = 0;i < m;i++){
    for(j = 0;j < m;j++){
      corr(i,j) = spearman(Z.col(i),Z.col(j));
    }
  }
  return corr;
}



// [[Rcpp::export]]
vec T_rho(mat Z){
  mat corr = Spear_mat(Z);
  mat cor2 = corr%corr;
  int N = Z.n_rows;
  int m = Z.n_cols;
  vec input(m,fill::zeros);
  vec output(2);
  cor2.diag() = input;
  double sum_rhopq = accu(cor2)/2.0;
  double C = m*(m-1)/2.0/(N-1);
  double T_rho= sum_rhopq - 1.0*C;
  double  gamma = 1.0*m/N;
  output(0) = T_rho;
  double pvalue;
  pvalue = 2.0*(1.0- normcdf(abs(T_rho),0.0,gamma));
  output(1) = pvalue;
  return output;
}

// [[Rcpp::export]]
vec S_Nm(mat Z){
  mat corr = Spear_mat(Z);
  mat cor2 = corr%corr;
  int N = Z.n_rows;
  int m = Z.n_cols;
  vec input(m,fill::zeros);
  vec output(2);
  cor2.diag() = input;
  double sum_rhopq = accu(cor2)/2.0;
  double C = m*(m-1)/2.0/(N-1);
  double UP,DOWN,Sigma_2_Nm,S_Nm;
  UP = m*(m-1)*(25*pow(N,3.0)-57*pow(N,2.0)-40*N+108);
  DOWN = 25*pow((N-1),3.0)*N*(N+1);
  Sigma_2_Nm = 1.0*UP/DOWN;
  S_Nm = pow(Sigma_2_Nm,-0.5)*(sum_rhopq-C);
  output(0) = S_Nm;
  double pvalue;
  pvalue = 2.0*(1.0 - normcdf(abs(S_Nm)));
  output(1) = pvalue;
  return output;
}

// [[Rcpp::export]]
vec t_nm(mat Z){
  int N = Z.n_rows;
  int n = N-1;
  int m = Z.n_cols;
  mat R1 = cor(Z);
  R1.diag() = zeros(R1.n_cols); 
  double Tn,esd,pvalue;
  vec out(2);
  
  Tn = accu(R1%R1)/2.0-m*(m-1)/2.0/n;
  esd = sqrt(m*(m-1)*(n-1.0)/(n+2.0)/n/n);
  Tn = Tn*1.0/esd;
  pvalue = 2*(1-normcdf(abs(Tn)));
  out(0) = Tn;
  out(1) = pvalue;
  
  return out;
}

// [[Rcpp::export]]
vec T3(mat Z){
  mat corr = cor(Z);
  mat cor2 = corr*corr;
  mat corE = corr%corr;
  mat cor4 = corE%corE;
  int N = Z.n_rows;
  int m = Z.n_cols;
  vec output(2);
  float lambda_3,a_20,a_40;
  lambda_3 = (N-1)*(trace(cor2)-pow(trace(corr),2)/(N-1))/(N-2)/(trace(corE));
  a_20 = (N-1)*(trace(corE))/m/(N+1);
  a_40 = trace(cor4)/m;
  float T3=(N-1)*(lambda_3-1)/(2*sqrt(1-a_40/(m*pow(a_20,2))));
  output(0) = T3;
  double pvalue;
  pvalue = 2.0*(1 - normcdf(abs(T3)));
  output(1) = pvalue;
  return output;
}

// [[Rcpp::export]]
//gendata multivariate normal distribution
mat GenData(int n, int m,float rho = 0){
  vec meanv(m, fill::zeros);
  mat A(m, m, fill::eye),B(m, m, fill::ones);
  mat covm = (1-rho)*A + rho*B;
  mat X = mvnrnd(meanv, covm, n);
  return X.t();
}//not t and cauchy

vec func(mat Z){
  vec X(2,fill::none);
  return X;
}

// [[Rcpp::export]]
mat RESULT(vec N,vec M,int index,float rho = 0,int nsim = 5000,float alpha = 0.05){
  int length_n = N.n_elem;
  int length_m = M.n_elem;
  mat RES(length_n,length_m,fill::none);
  mat X;
  vec res(nsim);
  int i,j,k;
  vec (*funcPtr)(mat);//funcPtr is short for 'function pointer'
  vec fit(1);
  vec test;
  
  if(index == 1){
    funcPtr = &S_Nm;
  }
  else if (index== 2){
    funcPtr = &T_rho;
  }
  else if (index== 3){
    funcPtr = &t_nm;
  }
  else if (index== 4){
    funcPtr = &T3;
  }
  else{
    funcPtr = &func;
  }
  
  for(i = 0;i < length_n;i++){
    for(j = 0;j < length_m;j++){
      for(k = 0;k < nsim;k++){
        X = GenData(N(i),M(j),rho);
        fit = funcPtr(X);
        res(k) = 1*( fit(1) < alpha);
      }
      RES(i,j) = accu(res)/nsim/1.0;
    }
  }
  return RES;
}
