library(glmnet)

get_bestfit <- function(cvfit, type.measure, exclude.zero=FALSE) {
  bestouts <- cvfit$measures[[type.measure]]
  extreme  <- ifelse(type.measure == "rmse", min, max)
  
  if (exclude.zero) {
    val <- extreme(bestouts[cvfit$nzero>0])
  } else {
    val <- extreme(bestouts)
  }
  
  ind <- which(bestouts == val)
  
  bestfit  <- list(
    measure = type.measure, 
    val     = val, 
    ind     = ind, 
    lam     = cvfit$lambda[ind], 
    preval  = cvfit$fit.preval[,ind], 
    nzero   = cvfit$nzero[ind], 
    coef    = coef(cvfit, s=cvfit$lambda[ind])
  )
  bestfit
}

run_cvglmnet <- function(X, y, foldid=NULL, keep=T, parallel=T, type.measure="rsq", exclude.zero=FALSE, ...) 
{
  if (!(type.measure %in% c("rsq", "r", "rmse"))) stop("unknown type.measure: ", type.measure)
  
  rmse <- function(x1, x2) sqrt(mean((x1-x2)^2))
  
  #cvfit  <- cv.glmnet(X, y, keep=keep, parallel=parallel)
  cvfit  <- cv.glmnet(X, y, keep=keep, parallel=parallel, foldid=foldid, ...)
  
  rs     <- sapply(1:length(cvfit$lambda), function(i) cor(cvfit$fit.preval[,i], y))
  rsqs   <- rs^2
  rmses  <- sapply(1:length(cvfit$lambda), function(i) rmse(cvfit$fit.preval[,i], y))
  cvfit$measures <- list(
    r = rs, 
    rsq = rsqs, 
    rmse = rmses
  )
  
  cvfit$bestfit <- get_bestfit(cvfit, type.measure, exclude.zero)
  
  return(cvfit)
}



# Make function to run a repeated CV

run_repeated_cvglmnet <- function(X, y, family="gaussian", parallel=F, ezero=F, nreps=10, cfolds=NULL, ...)
{
  # ezero: means if any fits with no coefficients should be excluded
  
  #nreps <- 10; ezero <- F; parallel <- F; family <- "gaussian"
  
  library(cvTools)
  if (is.null(cfolds)) {
    cfolds <- cvFolds(nrow(X), K=10, R=nreps)
  }
  
  # Predict y with model X repeated nreps
  cvrep.fits <- llply(1:nreps, function(i) {
    foldid <- cfolds$which[cfolds$subsets[,i]]
    run_cvglmnet(X, y, foldid=foldid, parallel=F, exclude.zero=ezero, family=family, keep=T, ...)
  }, .parallel=parallel)
  
  # Calculate the average model fits across the 10 repeats
  r2s <- sapply(cvrep.fits, function(cvfit) cvfit$bestfit$val)
  mean.r2s <- mean(r2s)
  
  # Combine the predicted values
  predvals <- sapply(cvrep.fits, function(x) x$bestfit$preval)
  predvals <- rowMeans(predvals)
  
  # Get the residuals
  resids <- lm(y ~ predvals)$residuals
  
  # Combine the r2 and # non-zero values across the lambdas
  lambda <- cvrep.fits[[1]]$lambda
  cvr2   <- matrix(cvrep.fits[[1]]$measures$rsq, nrow=1)
  cvnz   <- matrix(cvrep.fits[[1]]$nzero, nrow=1)
  for (i in 2:nreps) {
    lambda0 <- cvrep.fits[[i]]$lambda
    rsq0    <- cvrep.fits[[i]]$measures$rsq
    nzeros0 <- cvrep.fits[[i]]$nzero
    
    newlambda <- intersect(lambda, lambda0)
    cvr2      <- rbind(cvr2[,is.element(lambda, newlambda)], 
                       rsq0[is.element(lambda0, newlambda)])
    cvnz      <- rbind(cvnz[,is.element(lambda, newlambda)], 
                       nzeros0[is.element(lambda0, newlambda)])
    lambda    <- newlambda
  }
  # average across the repeats
  df.perlambda <- cbind(lambda=lambda, rsq=colMeans(cvr2), nzeros=colMeans(cvnz))
  
  # Another way to get the max r^2 instead of mean.r2s
  if (ezero) {
    max.rsq <- max(df.perlambda[df.perlambda[,3] > 0, 2])
  } else {
    max.rsq <- max(df.perlambda[,2])
  }
  wmax.rsq<- df.perlambda[,2] == max.rsq
  lam <- df.perlambda[wmax.rsq,1] # selected lambda
  max.nz <- df.perlambda[wmax.rsq,3]
  
  # Determine the coefficients from the full model
  # using the lambda with the best r2
  ## get full model fit first
  full.fit <- glmnet(X, y, family=family, lambda=lam)#, ...)
  ## get the coefficients
  coefs <- as.vector(full.fit$beta)
  
  list(mean.max.r2s=mean.r2s, fitted=predvals, resids=resids, 
       df=df.perlambda, max.r2=max.rsq, max.nonzero=max.nz, 
       lambda=lam, coefs=coefs, 
       rep.fits=cvrep.fits, full.fit=full.fit)
}
