// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
NumericMatrix extractFeatureMatrix(Rcpp::List massList,Rcpp::List intList, NumericVector Ref ,int length_ref, int length_spectra, double Bin_range) {
  double right;
  double left;
  arma::uvec bin_index;
  arma::vec temp_col(length_ref);
  arma::mat store;
// Bin_range=0.002
  arma::vec temp_mas;
  arma::vec temp_int;
  arma::vec ref=as<arma::vec>(Ref);
  
   for(int i=0;i<length_spectra;i++){

   temp_mas=as<arma::vec>(massList(i));
   temp_int=as<arma::vec>(intList(i));
   
   for(int j=0;j<length_ref;j++){
    right=ref(j)*(1+Bin_range);
    left=ref(j)*(1-Bin_range);
bin_index=find(temp_mas<=right && temp_mas>=left);
temp_col(j)=sum(temp_int.elem(bin_index));
}
   store.insert_cols(store.n_cols, temp_col);
   
 }
 
  return (Rcpp::wrap(store));
}


// [[Rcpp::export]]
NumericMatrix extractFeatureMatrix_V2(Rcpp::List massList,Rcpp::List intList, NumericVector Ref ,int length_ref, int length_spectra, double Bin_range) {
  double right;
  double left;
  arma::uvec bin_index;
  arma::vec temp_col(length_ref);
  arma::mat store;
  // Bin_range=0.002
  arma::vec temp_mas;
  arma::vec temp_int;
  arma::vec ref=as<arma::vec>(Ref);
  
  for(int i=0;i<length_spectra;i++){
    
    temp_mas=as<arma::vec>(massList(i));
    temp_int=as<arma::vec>(intList(i));
    
    for(int j=0;j<length_ref;j++){

      
      if(j==0){
        left=ref(j);
        right=(ref(j)+ref(j+1))/2;
      } else if(j==(length_ref-1)){
        right=ref(j);
        left=(ref(j)+ref(j-1))/2;
      } else{
        
        right=(ref(j)+ref(j+1))/2;
        left=(ref(j)+ref(j-1))/2;
      }
      
 
      
      bin_index=find(temp_mas<=right && temp_mas>=left);
      temp_col(j)=sum(temp_int.elem(bin_index));
    }
    store.insert_cols(store.n_cols, temp_col);
    
  }
  
  return (Rcpp::wrap(store));
}


// [[Rcpp::export]]
NumericMatrix extractFeatureMatrix2(Rcpp::List massList,Rcpp::List intList, NumericVector Ref ,int length_ref, int length_spectra, double halfwindowsize) {
  double right;
  double left;
  arma::uvec bin_index;
  arma::vec temp_col(length_ref);
  arma::mat store;
  // Bin_range=0.002
  arma::vec temp_mas;
  arma::vec temp_int;
  arma::vec ref=as<arma::vec>(Ref);
  
  for(int i=0;i<length_spectra;i++){
    
    temp_mas=as<arma::vec>(massList(i));
    temp_int=as<arma::vec>(intList(i));
    
    for(int j=0;j<length_ref;j++){
      right=ref(j)+halfwindowsize;
      left=ref(j)-halfwindowsize;
      bin_index=find(temp_mas<=right && temp_mas>=left);
      temp_col(j)=sum(temp_int.elem(bin_index));
    }
    store.insert_cols(store.n_cols, temp_col);
    
  }
  
  
  
  return (Rcpp::wrap(store));
}

/*
// [[Rcpp::export]]
NumericMatrix extractFeatureMatrix3(Rcpp::List massList,Rcpp::List intList, NumericVector Ref ,int length_ref, int length_spectra, double halfwindowsize_points) {
  double right;
  double left;
  arma::uvec bin_index;
  arma::vec temp_col(length_ref);
  arma::mat store;
  // Bin_range=0.002
  arma::vec temp_mas;
  arma::vec temp_int;
  arma::vec ref=as<arma::vec>(Ref);
  
  for(int i=0;i<length_spectra;i++){
    
    temp_mas=as<arma::vec>(massList(i));
    temp_int=as<arma::vec>(intList(i));
    
    for(int j=0;j<length_ref;j++){
      if((j!=0) && (j!=(length_ref-1))){
        right=ref(j+halfwindowsize_points);
        left=ref(j-halfwindowsize_points);
      } else if(j==0){
        right=ref(j+halfwindowsize_points);
        left=ref(j);
        
      } else if (j==(length_ref-1)){
        right=ref(j);
        left=ref(j-halfwindowsize_points);
      }
      

 
      bin_index=find(temp_mas<=right && temp_mas>=left);
      temp_col(j)=sum(temp_int.elem(bin_index));
    }
    store.insert_cols(store.n_cols, temp_col);
    
  }
  
  
  
  return (Rcpp::wrap(store));
}
*/
