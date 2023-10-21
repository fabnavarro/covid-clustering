require(Clustfun)
require(fda)
require(mixtools)
require(HDclassif)
require(parallel)
require(icamix)
require(VGAM)
library(rwavelet)
library(far)


rsimu4_1 <- function(n,
                     bdelta=0,
                     covariateeffect=0,
                     model=1){
  delta <- sample(0:bdelta, n, replace=TRUE)
  size <- 512
  z <- sort(sample(1:3, n, replace=TRUE, prob = c(.5,.25,.25)))
  obs <- matrix(NA, n, size)
  if (model==1){
    effect <- .5
    fbase <- sin(2.5 * pi * (1:(size))/size)
    fbasebis <- sin((2.5-effect) * pi * (1:(size))/size)
    fbase <- fbase * (fbase>0)
    fbasebis <- fbasebis * (fbasebis>0)
    r <- rbind(fbase, (1+effect) * fbase, fbasebis)
    for (i in 1:n) obs[i,] <- c(rep(0, delta[i]), r[z[i], (1:(size-delta[i]))] )
  }else if (model==2){
    x <- seq(0,1,length.out = size)
    f1 <- sin(5*pi*x)+ sin(2*pi*x)
    fg1 <- t(replicate(sum(z==1),f1+rnorm(size)))
    gamma1 <- gamma2 <-  matrix(0, nrow=4, ncol = 4, byrow = T)
    diag(gamma1) <- c(0.672, 0.228, 0.9, 0.34)
    gamma1[1,2] <- -0.134
    gamma1[2,1] <- 0.364
    diag(gamma2) <- diag(gamma1)
    fg2 <- simul.far(m=sum(z==2), n=size, base=base.simul.far(24, 5),
                     d.rho=diag(c(0.45, 0.9, 0.34, 0.45)), alpha=gamma1, cst1=1)
    fg3 <- simul.far(m=sum(z==3), n=size, base=base.simul.far(24, 5), d.rho=diag(c(0.45, 0.9, 0.34, 0.45)),
                     alpha=gamma2, cst1=1)
    fall <- rbind(fg1,fg2$var,fg3$var)
    for (i in 1:n) obs[i,] <- c(rep(0, delta[i]), fall[i, (1:(size-delta[i]))] )
  }
  covariates <- matrix(rnorm(n * 2), n, 2)
  obsnorm <- obs
  mu <- rep(NA, n)
  for (i in 1:n){
    obs[i,] <- obs[i,] + sqrt(0.5) * rnorm(size)
    obsnorm[i,] <- obs[i,]
    obs[i,] <- obs[i,] *  (covariateeffect * ((sum(covariates[i,]/sqrt(2)))**2) + 1) / (1+covariateeffect)
    mu[i] <- (covariateeffect * ((sum(covariates[i,]/sqrt(2)))**2) + 1) / (1+covariateeffect)
  }
  list(obs=obs, obsnorm=obsnorm, z=z,  covariates=covariates, bdelta=bdelta, covariateeffect=covariateeffect, mu = mu, model=model)
}

oneBasisComparison <- function(ech){
  # Estimatio of the basis coefficients
  zTInp <- clustfun(ech$obs, 3, type = "TI", covariates = NULL)$clustering[[1]]$zhat
  coeffsplines <- t((smooth.basis(1:ncol(ech$obs), t(ech$obs), create.bspline.basis(c(1, ncol(ech$obs)), nbasis=9, norder=4), fdnames=list("Time", "Obs", "Deg C"))$fd)$coefs)
  coeffexpo <- t((smooth.basis((1:ncol(ech$obs))/ncol(ech$obs), t(ech$obs), create.power.basis(c(0,1), nbasis=9), fdnames=list("Time", "Obs", "Deg C"))$fd)$coefs)
  coeffpolygonal <- t((smooth.basis(1:ncol(ech$obs), t(ech$obs), create.polygonal.basis(seq(1, ncol(ech$obs), length.out = 9)), fdnames=list("Time", "Obs", "Deg C"))$fd)$coefs)
  zsplinesnp <- apply(npMSL(scale(coeffsplines), mu0=3,  bw=nrow(ech$obs)**(-1/5), verb = FALSE)$posteriors, 1, which.max)
  zexponp <- apply(npMSL(scale(coeffexpo), mu0=3,  bw=nrow(ech$obs)**(-1/5), verb = FALSE)$posteriors, 1, which.max)
  zpolygonalnp <- apply(npMSL(scale(coeffpolygonal), mu0=3,  bw=nrow(ech$obs)**(-1/5), verb = FALSE)$posteriors, 1, which.max)
  # output
  data.frame(bdelta=rep(ech$bdelta,4),
             covariateeffect=rep(ech$covariateeffect, 4),
             model=rep(ech$model,4),
             size=rep(nrow(ech$obs), 4),
             basis=c("TI", "Expo", "Bsplines", "polyg"),
             value=c(ARI(ech$z, zTInp), ARI(ech$z, zexponp), ARI(ech$z, zsplinesnp), ARI(ech$z, zpolygonalnp))
  )
}

locfct <- function(eval, index, Ystar){
  u <- dnorm(index, eval, sd = length(Ystar)**(-1/5))
  sum(u*Ystar)/sum(u)
}

oneCovariateComparison <- function(ech){
  # Estimatio of the basis coefficients
  zTInp <- clustfun(ech$obs, 3, type = "TI", covariates = NULL)$clustering[[1]]$zhat
  zTInpcov <- clustfun(ech$obs, 3, type = "TI", covariates = ech$covariates)$clustering[[1]]$zhat
  oracle <- clustfun(ech$obsnorm, 3, type = "TI", covariates = NULL)$clustering[[1]]$zhat
  data.frame(bdelta=rep(ech$bdelta, 3),
             covariateeffect=rep(ech$covariateeffect, 3),
             model=rep(ech$model, 3),
             size=rep(nrow(ech$obs), 3),
             cov=c("no-cov","cov","oracle"),
             value=c(ARI(ech$z, zTInp), ARI(ech$z, zTInpcov), ARI(ech$z, oracle))
  )
}

oneCovariateInvestigate <- function(ech){
  # Estimatio of the basis coefficients
  rescov <- clustfun(ech$obs, 3, type = "TI", covariates = ech$covariates)
  oracle <- clustfun(ech$obsnorm, 3, type = "TI", covariates = NULL)
  data.frame(bdelta=rep(ech$bdelta, 3),
             covariateeffect=rep(ech$covariateeffect, 3),
             model=rep(ech$model, 3),
             size=rep(nrow(ech$obs), 3),
             cov=c("mu","param","obstocluster"),
             value=c(mean((exp(rescov$reg$m) - ech$mu)**2), mean((rescov$reg$siparam - 1/sqrt(2))**2), mean(rowSums((rescov$obstocluster - oracle$obstocluster)**2)))
  )
}

oneClusteringComparison <- function(ech){
  # Estimatio of the basis coefficients
  decompo <- clustfun(ech$obs, 3, type = "TI", covariates = ech$covariates)
  ARIzNPdiago <- try(ARI(decompo$clustering[[1]]$zhat, ech$z))
  if (class(ARIzNPdiago)=="try-error") ARIzNPdiago <- 0
  ARIzHddc <- try(ARI(hddc(decompo$obstocluster, 3, model="ALL")$class, ech$z))
  if (class(ARIzHddc)=="try-error") ARIzHddc <- 0
  ARIzNPfull <- try(ARI(ESTIMATEDMEMBER(EMFASTICAALG(decompo$obstocluster, 3, h=decompo$clustering[[1]]$band, tol=1e-8)), ech$z), silent = T)
  loc <- 0
  while (class(ARIzNPfull)=="try-error"){
    loc <- loc + 1
    ARIzNPfull <- try(ARI(ESTIMATEDMEMBER(EMFASTICAALG(decompo$obstocluster, 3, h=decompo$clustering[[1]]$band, tol=1e-8)), ech$z), silent = T)
    if (loc>20){
      ARIzNPfull <- 0
      break;
    }
  }
  data.frame(bdelta=rep(ech$bdelta,3),
             covariateeffect=rep(ech$covariateeffect, 3),
             model=rep(ech$model,3),
             size=rep(nrow(ech$obs), 3),
             clustering=c("npdiago", "npfull", "GMMsparse"),
             value=c(ARIzNPdiago, ARIzNPfull,   ARIzHddc)
  )
}

