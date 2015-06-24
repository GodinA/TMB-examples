# R-INLA Spatial example from Thorson, J. T., & Coilin, M. (2014). Mixed
# effects: a unifying framework for statistical modelling in fisheries biology.
# ICES Journal of Marine Science. doi:10.1093/icesjms/fst048

# https://github.com/James-Thorson/mixed-effects/tree/master/Spatial_model

## Data: simulate 1000 data points from a Poisson distribution within a square 
## sampling domain, where the logarithm of the expected value increases from 
## West to East and has a residual pattern along the diagonal. Model: spatial 
## model using a Poisson distribution for measurement errors, including eastings
## as a covariate, and assuming spatial residuals follow a Gaussian Random Field
## (GRF).

rm(list = ls(all = TRUE))
setwd("~/Dropbox/TMB-examples/Spatial_model") 
#install_github("kaskr/adcomp/TMB") 

# load library
library(TMB)
require(splancs)
require(rgl)
require(INLA)
library(lattice)
library(gridExtra)

# replicate simulated data 
n = 1000
s1 = runif(n, min=0, max=1)
s2 = runif(n, min=0, max=1)
a = 1
ln_y_hat = 2.0*exp(-((s1-s2)/0.2)^2) + a*s1   # spatial process     
y = rpois(n, lambda=exp(ln_y_hat))                 
coords <- cbind(s1, s2)
data = data.frame(y=y,coords=coords,eastings=coords[,1])
data = na.omit(data)
coords = as.matrix(data[,2:3])
head(data)

# quick plot
with(data,plot(coords.s1,coords.s2,cex=sqrt(y+0.1)/2))

# boundr
pl.dom <- cbind( c(min(coords[,1]),max(coords[,1]),max(coords[,1]),min(coords[,1])), c(min(coords[,2]),min(coords[,2]),max(coords[,2]),max(coords[,2])))

# mesh
mesh <- inla.mesh.2d(coords, pl.dom, max.edge=1, offset=-0.5)
# plot(mesh,main='',asp=1)

# Fit model with R-INLA ========================
spde <- inla.spde2.matern(mesh, alpha=2)
A <- inla.spde.make.A(mesh, loc=coords)  

stk <- inla.stack(data=list(resp=data$y), 
  A=list(A,1,1), 
  effects=list(i=1:spde$n.spde, intercept=rep(1,nrow(data)), x=data$eastings), 
  tag='est')

formula = resp ~ 0 + intercept + x + f(i, model=spde)

result = inla(formula, family="poisson", 
  data=inla.stack.data(stk), 
  control.predictor=list(A=inla.stack.A(stk), compute=TRUE), 
  control.family=list(hyper = list(theta = list(initial = log(1/0.01^2), fixed = FALSE))), 
  keep=FALSE)##verbose=TRUE

summary(result)

# range of GRF 
rf <- inla.spde2.result(result, 'i', spde) 
c(mean=inla.emarginal(function(x) x, rf$marginals.range[[1]]), 
  q=inla.hpdmarginal(0.95, rf$marginals.range[[1]]))[c(1,2,3)]

# all parameters
round(tables <- t(sapply(c(result$marginals.fix, 
  result$marginals.hy[1:2], 
  result$marginals.range[1]), function(m) 
    c(mean=inla.emarginal(function(x) x, m), 
      lim=inla.hpdmarginal(0.95, m)))),4)

# Calculate true expected surface 
s1_proj = outer( seq(0,1,length=100), rep(1,100) )
s2_proj = outer( rep(1,100), seq(0,1,length=100) )
ln_y_exp = 2.0*exp(-((s1_proj-s2_proj)/0.2)^2) + a*s1_proj

# Caclualte estimated expected surface - INLA
mesh_proj = inla.mesh.projector(mesh,xlim=c(0,1), ylim=c(0,1), dims=c(100,100))

ln_y_hat = inla.mesh.project(mesh_proj, result$summary.ran$i$mean) + 
  result$summary.fix['intercept','mean'] + 
  result$summary.fix['x','mean']*s1_proj

Col = colorRampPalette(colors=c("blue","white","red"))
At = seq( min(ln_y_exp,ln_y_hat),max(ln_y_exp,ln_y_hat), length=100)

# quick plots
grid.arrange(
  levelplot(ln_y_exp, xlab='', ylab='', main='ln(True value)', col.regions=rev(topo.colors(99)), at=At, scales=list(draw=FALSE)),
  levelplot(ln_y_hat, xlab='', ylab='', main='ln(Exp. value)', col.regions=rev(topo.colors(99)), at=At, scales=list(draw=FALSE)), 
  nrow=1
)

# Fit model with TMB ========================
# Fixed effects 
X = model.matrix(~1+eastings,data=data);head(X)
data_tmb <- list(counts=data$y,meshidxloc=mesh$idx$loc-1,X=as.matrix(X))
str(data_tmb)

# SPDE part
spde_param<-inla.spde2.matern(mesh, alpha=2)$param.inla
# need to get all SparseMatrix elements of the SPDE object - to calculate the
# inverse-covariance matrix of GRF c("M0","M1","M2") is G0, G1, G2 in eqn (10)
# in Lindgren 2011
data_tmb$spde <- spde_param[c("M0","M1","M2")]
n_s = nrow(data_tmb$spde$M0);n_s # idem spde_param$n
str(data_tmb)

# Giving starting values for the parameters 
# Tested a few 
par <- list(beta=c(0,0),log_tau=-2.0,log_kappa=2.5,x=rep(0.0,n_s))
# par2 <- list(beta=c(0,0),log_tau=0,log_kappa=0,x=rep(0.0,n_s))
# par3 <- list(beta=c(0,0),log_tau=2.0,log_kappa=-2.5,x=rep(0.0,n_s))

# compile and load C++
compile("toy_spde.cpp")
dyn.load(dynlib("toy_spde"))

obj <- MakeADFun(data_tmb,par,random="x",DLL="toy_spde")
# obj2 <- MakeADFun(data_tmb,par2,random="x",DLL="spde_acg")
# obj3 <- MakeADFun(data_tmb,par3,random="x",DLL="spde_acg")

opt <- nlminb(obj$par,obj$fn,obj$gr)
# opt2 <- nlminb(obj2$par,obj2$fn,obj2$gr)
# opt3 <- nlminb(obj3$par,obj3$fn,obj3$gr)

# c(opt$objective,
#   opt2$objective,
#   opt3$objective)

opt$par

pl <- obj$env$parList()
ln_y_hat_tmb= inla.mesh.project(mesh_proj, pl$x/exp(opt$par["log_tau"])) + 
  opt$par[1] + 
  opt$par[2]*s1_proj

# plot all
grid.arrange(
  levelplot(ln_y_exp, xlab='', ylab='', main='ln(True value)', col.regions=rev(topo.colors(99)), at=At, scales=list(draw=FALSE)),
  levelplot(ln_y_hat, xlab='', ylab='', main='ln(Exp. value) (INLA)', col.regions=rev(topo.colors(99)), at=At, scales=list(draw=FALSE)), 
  levelplot(ln_y_hat_tmb, xlab='', ylab='', main='ln(Exp. value) (TMB)', col.regions=rev(topo.colors(99)), at=At, scales=list(draw=FALSE)), 
  nrow=3
)
# head(summary(sdreport(obj)),n=4)
