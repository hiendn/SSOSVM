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
Rcpp::List SquareHingeC(arma::mat& YMAT,  int DIM = 2, double EPSILON = 0.00001, bool returnAll = false, double rho=1.0) {

    int NN = YMAT.n_rows;
    double LAMBDA;
    LAMBDA = 1.0/NN;
    arma::vec THETA =  arma::zeros(DIM+1,1);
    
    arma::vec THETA_OLD = THETA;
    arma::mat THETA_list;
    
    // Make I_BAR_p
    arma::mat IBAR (DIM+1,DIM+1,arma::fill::eye);
    IBAR(0,0) = 0.0;
    
    arma::mat store(DIM+1,DIM+1,arma::fill::zeros);
    arma::mat THSQ2(DIM+1,1,arma::fill::zeros); 
    
    arma::mat Part2a(DIM+1,1,arma::fill::zeros); 
    arma::mat Part2b(DIM+1,1,arma::fill::zeros);
    arma::mat Part2c(DIM+1,1,arma::fill::zeros); 

    arma::mat Inter(DIM+1,DIM+1,arma::fill::zeros);
    
    if(returnAll){
      THETA_list = arma::zeros(NN, DIM+1);
    }else{
      THETA_list = arma::zeros(1, DIM+1);
    }
    
    //Initialize little psi vector
    double psi_temp = psiFun(THETA,  YMAT.row(1), EPSILON); 
    arma::vec psi = arma::zeros(NN,1);
    psi(0) = psi_temp;

    
    store =  YMAT.row(0).t()*YMAT.row(0);
    Part2b = YMAT.row(0).t()*((0.5/rho)*psi(0));
    Part2c = store*THETA;
    
      //Main loop
    for(int ii = 1; ii<NN; ii++) {
      //Update little psi vector
      psi(ii) =  psiFun(THETA,  YMAT.row(ii), EPSILON);
      
      // Turn psi into a column matrix
      //arma::vec psi_mat = psi.rows(0, ii);
      
      Inter = YMAT.row(ii).t()*YMAT.row(ii);
          
      store += Inter;
      
      Part2a  = arma::inv_sympd(store+(LAMBDA*NN*IBAR/rho));
      
      Part2b += YMAT.row(ii).t()*((0.5/rho)*psi(ii));
      
      Part2c += Inter*THETA_OLD;
      
      //THSQ2 = THSQ2a+THSQ2b;
      
      // Update Theta
      THETA = Part2a*Part2c+Part2a*Part2b;
      
      // Put old value of theta in a vector
      THETA_OLD = THETA;
    
      if(returnAll){
        THETA_list.row(ii)=THETA.t();
      }

    }
    
    Rcpp::List retList = Rcpp::List::create(
      Rcpp::Named("THETA")= export_vec(THETA),
      Rcpp::Named("NN")= NN,
      Rcpp::Named("DIM")= DIM,
      Rcpp::Named("THETA_list")= THETA_list,
      Rcpp::Named("PSI")=export_vec(psi)
  );
  
  
  return(retList);
}




//'@export
// [[Rcpp::export]]
Rcpp::List HingeC(arma::mat& YMAT,  int DIM = 2, double EPSILON = 0.00001, bool returnAll = false, double rho=1.0) {
  
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
  
  
  arma::mat  store(DIM+1,DIM+1,arma::fill::zeros);
  arma::mat  store2(DIM+1,DIM+1,arma::fill::zeros);
  store = YMAT.row(0).t()*YMAT.row(0)+4*LAMBDA*NN*IBAR;
  store2 = (YMAT.row(0)*(1.0/omega(0))).t()*(1.0+omega.row(0)) ;
  
 
  
  //Main loop
  for(int ii = 1; ii<NN; ii++) {
      
    //Update omega
   
    omega(ii) =  omegaFun(THETA,  YMAT.row(ii), EPSILON);
      
    store += YMAT.row(ii).t()*(1.0/omega(ii))*YMAT.row(ii);
    store2 += (YMAT.row(ii)*(1.0/omega(ii))).t()*(1.0+omega.row(ii)) ;
    
    //Compute theta update
    THETA = arma::inv_sympd(store)*(store2);
                                                                                                                                                                                
    
    if(returnAll){
      THETA_list.row(ii)=THETA.t();
    }
    
  }
  
  Rcpp::List retList = Rcpp::List::create(
    Rcpp::Named("THETA")= export_vec(THETA),
    Rcpp::Named("NN")= NN,
    Rcpp::Named("DIM")= DIM,
    Rcpp::Named("THETA_list")= THETA_list,
    Rcpp::Named("OMEGA")=export_vec(omega)
  );
  
  
  return(retList);
}

//'@export
// [[Rcpp::export]]
Rcpp::List LogisticC(arma::mat& YMAT, int DIM = 2, double EPSILON = 0.00001, bool returnAll = false, double rho=1.0) {
  
  int NN = YMAT.n_rows;
  double LAMBDA;
  LAMBDA = 1.0/NN;
  
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
  
  
  arma::mat store (DIM+1,DIM+1,arma::fill::zeros);
  arma::mat Part2(DIM+1,1,arma::fill::zeros);
  arma::mat Part2a(DIM+1,1,arma::fill::zeros);
  arma::mat Part2b(DIM+1,1,arma::fill::zeros);
  arma::mat Part2c(DIM+1,1,arma::fill::zeros); 
  
  arma::mat Inter(DIM+1,DIM+1,arma::fill::zeros);
  
  //Initialize omega vector
  double chi_temp = chiFun(THETA,  YMAT.row(1), EPSILON); 
  arma::vec chi = arma::zeros(NN);
  chi(0) = chi_temp;
    
    
  
  store = YMAT.row(0).t()*YMAT.row(0); // Only compute the inner product here
  Part2b = YMAT.row(0).t()*((4/rho)*chi(0));
  Part2c = store*THETA;
  
  //Main loop
  for(int ii = 1; ii<NN; ii++) {
    
    //Rcpp::Rcout << "Doing "<< ii << " of  " << NN << std::endl;
    
    //Update chi vector
    chi(ii) = chiFun(THETA,  YMAT.row(ii), EPSILON); 
    
    //Turn chi vector into column vector
    //arma::vec chi_mat = chi.rows(0, ii);
    
    Inter = YMAT.row(ii).t()*YMAT.row(ii);
    
    store += Inter;
    
    Part2a = arma::inv_sympd(store+(8*LAMBDA*NN*IBAR/rho));
    
    Part2b += YMAT.row(ii).t()*((4/rho)*chi(ii));
    
    Part2c += Inter*THETA_OLD;
    
    //Compute Theta
    THETA = Part2a*Part2c+Part2a*Part2b;
    
    // Store old theta away
    THETA_OLD = THETA ;
      
      if(returnAll){
        THETA_list.row(ii)=THETA.t();
      }
      
  }
  
  Rcpp::List retList = Rcpp::List::create(
    Rcpp::Named("THETA")= export_vec(THETA),
    Rcpp::Named("NN")= NN,
    Rcpp::Named("DIM")= DIM,
    Rcpp::Named("THETA_list")= THETA_list,
    Rcpp::Named("CHI")= export_vec(chi)
  );


return(retList);
}