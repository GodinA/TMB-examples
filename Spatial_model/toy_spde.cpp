// R-INLA Spatial example from Thorson and Coilin, ICES 2014
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
  using namespace R_inla; 
  using namespace density; 
  using namespace Eigen; 
  
  // data 
  DATA_VECTOR(counts); // y
  DATA_MATRIX(X); // matrix of fixed effect - Intercept + eastings
  DATA_IVECTOR(meshidxloc); // mesh loc
  DATA_STRUCT(spde,spde_t); // spde list object 
  
  // parameters list
  PARAMETER_VECTOR(beta); // betas
  PARAMETER(log_tau); // spatial variance
  PARAMETER(log_kappa); // spatial scale 
  PARAMETER_VECTOR(x);  // spatial vector 
  
  // transformations
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  
  // initialize negative log likelihood
  Type nll = 0.0; 
  
  // define precision matrix for the Gaussian weights w - eqn (10) in Lindgren et al. (2011)
  SparseMatrix<Type> Q = Q_spde(spde,kappa); 
  
  nll += GMRF(Q)(x); // Evaluate density of GMRF
  
  vector<Type> Xbeta = X*beta;  
  for(int i=0; i< counts.size(); i++){   
    // process model
    Type eta = Xbeta(i) + x(meshidxloc(i))/tau; // linear predictor
    Type lambda = exp(eta); // mean & variance poisson distribution 
    // observation model
    nll -= dpois(counts(i),lambda, true); 
  }
  
  return nll;
}

