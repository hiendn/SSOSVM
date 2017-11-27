// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]   

//'@importFrom Rcpp sourceCpp
//'@useDynLib SSOSVM

#include "RcppArmadillo.h"

Rcpp::NumericVector export_vec(arma::vec y)
{
  Rcpp::NumericVector tmp = Rcpp::wrap(y);
  tmp.attr("dim") = R_NilValue;
  return tmp;
}

double psiFun(arma::vec THETA,  arma::rowvec Yrow, double EPSILON ){
  double psi_temp = std::pow(std::sqrt(std::pow(1-arma::dot(THETA,Yrow),2.0)+EPSILON)+1-arma::dot(THETA,Yrow),2.0)/(2*std::sqrt(std::pow(1-arma::dot(THETA,Yrow),2.0)+EPSILON));
  return(psi_temp);
}

double chiFun(arma::vec THETA,  arma::rowvec Yrow, double EPSILON ){
  double chi_temp = std::exp(-arma::dot(THETA,Yrow))/(1+std::exp(-arma::dot(THETA,Yrow)));
  return(chi_temp);
}

double omegaFun(arma::vec THETA,  arma::rowvec Yrow, double EPSILON ){
  double omega_temp = std::sqrt(std::pow(1-arma::dot(THETA,Yrow),2.0)+EPSILON);
  return(omega_temp);
}



//'@export
// [[Rcpp::export]]
Rcpp::List SquareHingeC(arma::mat& YMAT,  int DIM = 2, double EPSILON = 0.00001, bool returnAll = false) {

    int NN = YMAT.n_rows;
    double LAMBDA = 1.0/NN;
    arma::vec THETA =  arma::zeros(DIM+1,1);
    arma::vec THETA_OLD = THETA;
    arma::mat THETA_list;
    
    // Make I_BAR_p
    arma::mat IBAR (DIM+1,DIM+1,arma::fill::eye);
    IBAR(0,0) = 0.0;
    

    if(returnAll){
      THETA_list = arma::zeros(NN, DIM+1);
    }else{
      THETA_list = arma::zeros(1, DIM+1);
    }
    
    //Initialize little psi vector
    double psi_temp = psiFun(THETA,  YMAT.row(1), EPSILON); 
    arma::vec psi = arma::zeros(NN,1);
    psi(0) = psi_temp;

      //Main loop
    for(int ii = 1; ii<NN; ii++) {
      
      //Update little psi vector
      psi(ii) =  psiFun(THETA,  YMAT.row(ii), EPSILON);
      
      // Put old value of theta in a vector
      THETA_OLD = THETA;
      
      // Turn psi into a column matrix
      arma::vec psi_mat = psi.rows(0, ii);
          
      arma::mat THSQ =  YMAT.rows(0,ii).t() * YMAT.rows(0,ii);
      arma::mat THSQ2 = YMAT.rows(0,ii).t() * ((YMAT.rows(0,ii)*THETA_OLD)+0.5*psi_mat);
      
      // Update Theta
      THETA = arma::pinv(THSQ+LAMBDA*NN*IBAR)*THSQ2;
    
      if(returnAll){
        THETA_list.row(ii)=THETA.t();
      }

    }
    
    Rcpp::List retList = Rcpp::List::create(
      Rcpp::Named("THETA")= export_vec(THETA),
      Rcpp::Named("NN")= NN,
      Rcpp::Named("DIM")= DIM,
      Rcpp::Named("THETA_list")= THETA_list
  );
  
  
  return(retList);
}



//'@export
// [[Rcpp::export]]
Rcpp::List HingeC(arma::mat& YMAT,  int DIM = 2, double EPSILON = 0.00001, bool returnAll = false) {
  
  int NN = YMAT.n_rows;
  double LAMBDA = 1.0/NN;
  arma::vec THETA =  arma::zeros(DIM+1,1);
  arma::vec THETA_OLD = THETA;
  arma::mat THETA_list;
 
  if(returnAll){
    THETA_list = arma::zeros(NN, DIM+1);
  }else{
    THETA_list = arma::zeros(1, DIM+1);
  }
  
  // Make I_BAR_p
  arma::mat IBAR(DIM+1,DIM+1,arma::fill::eye);
  IBAR(0,0) = 0.0;
  
  //Initialize omega vector
  double omega_temp = omegaFun(THETA,  YMAT.row(1), EPSILON); 
  arma::vec omega = arma::zeros(NN,1);
  omega(0) = omega_temp;
  
  //Main loop
  for(int ii = 1; ii<NN; ii++) {
      
    //Update omega
   
    omega(ii) =  omegaFun(THETA,  YMAT.row(ii), EPSILON);
      
    //Make inverse OMEGA matrix
    arma::vec omega_mat = omega.rows(0, ii);
    arma::mat OMEGA_INV =  arma::diagmat(1.0/omega_mat);
    arma::vec ONES(ii+1, arma::fill::ones);

    //Compute theta update
    THETA = arma::pinv(YMAT.rows(0,ii).t()*OMEGA_INV*YMAT.rows(0,ii)+4*NN*IBAR*LAMBDA)*(YMAT.rows(0,ii).t()*(OMEGA_INV)*(ONES+omega_mat));
    
    
    if(returnAll){
      THETA_list.row(ii)=THETA.t();
    }
    
  }
  
  Rcpp::List retList = Rcpp::List::create(
    Rcpp::Named("THETA")= export_vec(THETA),
    Rcpp::Named("NN")= NN,
    Rcpp::Named("DIM")= DIM,
    Rcpp::Named("THETA_list")= THETA_list
  );
  
  
  return(retList);
}


//'@export
// [[Rcpp::export]]
Rcpp::List LogisticC(arma::mat& YMAT, int DIM = 2, double EPSILON = 0.00001, bool returnAll = false) {
  
  int NN = YMAT.n_rows;
  double LAMBDA = 1.0/NN;
  arma::vec THETA =  arma::zeros(DIM+1);
  arma::vec THETA_OLD = THETA;
  arma::mat THETA_list;
  
  if(returnAll){
    THETA_list = arma::zeros(NN, DIM+1);
  }else{
    THETA_list = arma::zeros(1, DIM+1);
  }
  
  // Make I_BAR_p
  arma::mat IBAR (DIM+1,DIM+1,arma::fill::eye);
  IBAR(0,0) = 0.0;
  
  //Initialize omega vector
  double chi_temp = chiFun(THETA,  YMAT.row(1), EPSILON); 
  arma::vec chi = arma::zeros(NN);
  chi(0) = chi_temp;

  //Main loop
  for(int ii = 1; ii<NN; ii++) {
    
    //Update chi vector
    chi(ii) = chiFun(THETA,  YMAT.row(ii), EPSILON); 
      
    // Store old theta away
    THETA_OLD = THETA ;
    
    //Turn chi vector into column vector
    arma::vec chi_mat = chi.rows(0, ii);
    //Compute Theta
    THETA = arma::pinv(YMAT.rows(0,ii).t()*YMAT.rows(0,ii)+8*LAMBDA*NN*IBAR)*YMAT.rows(0,ii).t()*(YMAT.rows(0,ii)*THETA_OLD+4*chi_mat);

    if(returnAll){
      THETA_list.row(ii)=THETA.t();
    }
    
  }
  
  Rcpp::List retList = Rcpp::List::create(
    Rcpp::Named("THETA")= export_vec(THETA),
    Rcpp::Named("NN")= NN,
    Rcpp::Named("DIM")= DIM,
    Rcpp::Named("THETA_list")= THETA_list
  );
  
  
  return(retList);
}
