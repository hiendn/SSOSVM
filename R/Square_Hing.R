
generateSim <- function(NN = 10^4, seed=NULL) {
  if(!is.null(seed)){set.seed(seed)}
  # Simulation
  XX <- mvtnorm::rmvnorm(NN/2,c(0,-2)+4)
  XX <- rbind(XX,mvtnorm::rmvnorm(NN/2,c(0,2)+4))
  # XX <- scale(XX)
  SAMPLE <- sample(1:NN)
  XX <- XX[SAMPLE,]
  YY <- c(rep(-1,NN/2),rep(1,NN/2))
  YY <- YY[SAMPLE]
  YMAT <- cbind(YY,YY*XX)

  return(list(YMAT=YMAT, XX=XX, YY=YY))
}

SquareHinge<-function(YMAT,DIM=2, EPSILON = 10^-5, returnAll = F, verbose=F){
  NN<-nrow(YMAT)
  THETA <- matrix(0,DIM+1,1)
  LAMBDA <- 1/NN
  THETA_OLD <- THETA

  # Initialize little psi vector
  psi <- (sqrt((1-sum(THETA*YMAT[1,]))^2+EPSILON)+1-sum(THETA*YMAT[1,]))^2/(2*sqrt((1-sum(THETA*YMAT[1,]))^2+EPSILON))
 
                                                                              
  if(returnAll){
    THETA_list<-array(0,c(NN,DIM+1))
  }else{
    THETA_list<-NA
  }

  # Main loop
  for (ii in 2:NN) {
    # Update little psi vector
    psi <- c(psi,(sqrt((1-sum(THETA*YMAT[ii,]))^2+EPSILON)+1-sum(THETA*YMAT[ii,]))^2/
               (2*sqrt((1-sum(THETA*YMAT[ii,]))^2+EPSILON)))
    # Put old value of theta in a vector
    THETA_OLD <- THETA
    # Turn psi into a column matrix
    psi_mat <- matrix(psi,ii,1)
    # Make I_BAR_p
    IBAR <- diag(c(0,rep(1,DIM)))
    # Update Theta
    THETA <- MASS::ginv(t(YMAT[1:ii,])%*%YMAT[1:ii,]+LAMBDA*NN*IBAR)%*%
              t(YMAT[1:ii,])%*%(YMAT[1:ii,]%*%THETA_OLD+0.5*psi_mat)
    
    

    if(verbose){
      message(paste(ii,THETA))
    }

    if(returnAll){
      THETA_list[ii,]<-THETA
    }

    
  }
  return(list(THETA=THETA, NN=NN, DIM=DIM, THETA_list=THETA_list))

}

Hinge<-function(YMAT,DIM=2, EPSILON = 10^-5, returnAll = F, verbose=F){

  NN<-nrow(YMAT)
  THETA <- matrix(0,DIM+1,1)
  LAMBDA <- 1/NN
  THETA_OLD <- THETA

  # Initialize little omega vector
  omega <- sqrt((1-sum(YMAT[1,]*THETA))^2+EPSILON)


  if(returnAll){
    THETA_list<-array(0,c(NN,DIM+1))
  }else{
    THETA_list<-NA
  }

  # Main loop
  for  (ii in 2:NN) {
      # Update omega
      omega <- c(omega,sqrt((1-sum(YMAT[ii,]*THETA))^2+EPSILON))
      # Make inverse OMEGA matrix
      OMEGA_INV <- diag(1/omega)
      # Make IBAR_p
      IBAR <- diag(c(0,rep(1,DIM)))
      # Make identity matrix
      EYE <- diag(ii)
      # Make vector of ones
      ONES <- as.matrix(rep(1,ii),ii,1)
      # turn omega vector into a column vector
      omega_mat <- as.matrix(omega,ii,1)
      # Compute theta update
      THETA <- MASS::ginv(t(YMAT[1:ii,])%*%(OMEGA_INV)%*%YMAT[1:ii,]+4*NN*IBAR*LAMBDA)%*%
        t(YMAT[1:ii,])%*%(OMEGA_INV)%*%(ONES+omega_mat)

      if(verbose){
        message(paste(ii,THETA))
      }

      if(returnAll){
        THETA_list[ii,]<-THETA
      }
    }

    return(list(THETA=THETA, NN=NN, DIM=DIM, THETA_list=THETA_list))
}


Logistic<-function(YMAT,DIM=2, EPSILON = 10^-5, returnAll = F, verbose=F){

  NN<-nrow(YMAT)
  THETA <- matrix(0,DIM+1,1)
  LAMBDA <- 1/NN
  THETA_OLD <- THETA

    #
  # Initialize chi vector
  chi <- exp(-sum(THETA*YMAT[1,]))/(1+exp(-sum(THETA*YMAT[1,])))


  if(returnAll){
    THETA_list<-array(0,c(NN,DIM+1))
  }else{
    THETA_list<-NA
  }

  # Main loop
  for (ii in 2:NN) {
    # Update chi vector
    chi <- c(chi,exp(-sum(THETA*YMAT[ii,]))/(1+exp(-sum(THETA*YMAT[ii,]))))

    # Store old theta away
    THETA_OLD <- THETA

    # Turn chi vector into column vector
    chi_mat <- matrix(chi,ii,1)

    # Make I_P_BAR
    IBAR <- diag(c(0,rep(1,DIM)))

    # Compute Theta
    THETA <- MASS::ginv(t(YMAT[1:ii,])%*%YMAT[1:ii,]+8*LAMBDA*NN*IBAR)%*%
      t(YMAT[1:ii,])%*%(YMAT[1:ii,]%*%THETA_OLD+4*chi_mat)


    if(verbose){
      message(paste(ii,THETA))
    }

    if(returnAll){
      THETA_list[ii,]<-THETA
    }
  }

  return(list(THETA=THETA, NN=NN, DIM=DIM, THETA_list=THETA_list))

}



# # Initialization (YMAT is BIG_Y_BAR in text)

# 
# # Plots only make sense for 2D XX
# plot(YMAT$XX,pch=(YMAT$YY+2),col=heat.colors(sq1$NN)[1:sq1$NN])
# curve(-sq1$THETA[1]/sq1$THETA[3]-sq1$THETA[2]/sq1$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='black',add=T)
# curve(-h1$THETA[1]/h1$THETA[3]-h1$THETA[2]/h1$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='blue',add=T)
# curve(-l1$THETA[1]/sq1$THETA[3]-l1$THETA[2]/l1$THETA[3]*x,from=min(YMAT$XX),to=max(YMAT$XX),col='green',add=T)
# 


