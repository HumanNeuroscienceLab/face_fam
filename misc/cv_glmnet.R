library(glmnet)

get_bestfit <- function(cvfit, type.measure, exclude.zero=FALSE) {
  bestouts <- cvfit$measures[[type.measure]]
  extreme  <- ifelse(type.measure == "rmse", min, max)
  
  if (exclude.zero) {
    val <- extreme(bestouts[cvfit$nzero>0])
  } else {
    val <- extreme(bestouts)
  }
  
  ind <- which(bestouts == val)[1] # make the less complex model if more than one
  
  bestfit  <- list(
    measure = type.measure, 
    val     = val, 
    ind     = ind, 
    lam     = cvfit$lambda[ind], 
    nzero   = cvfit$nzero[ind], 
    coef    = coef(cvfit, s=cvfit$lambda[ind])
  )
  
  if (cvfit$family == "binomial") {
    bestfit$prob <- cvfit$probs[,ind]
    bestfit$class <- cvfit$classes[,ind]
    bestfit$confMat <- caret::confusionMatrix(data=cvfit$classes[,ind], ref=cvfit$y)
  } else if (cvfit$family == "multinomial") {
    bestfit$prob <- cvfit$probs[,,ind]
    bestfit$class <- cvfit$classes[,ind]
    bestfit$confMat <- caret::confusionMatrix(data=cvfit$classes[,ind], ref=cvfit$y)
  } else {
    bestfit$preval <- cvfit$fit.preval[,ind]
  }
  
  bestfit
}

run_cvglmnet <- function(X, y, keep=T, family="gaussian", type.measure="rsq", exclude.zero=FALSE, ...) 
{
  if (!(type.measure %in% c("rsq", "r", "rmse", "Accuracy", "Kappa"))) stop("unknown type.measure: ", type.measure)
  
  rmse <- function(x1, x2) sqrt(mean((x1-x2)^2))
  
  cvfit  <- cv.glmnet(X, y, keep=keep, family=family, ...)
  
  if (family=="gaussian") {
    rs     <- sapply(1:length(cvfit$lambda), function(i) cor(cvfit$fit.preval[,i], y))
    rsqs   <- rs^2
    rmses  <- sapply(1:length(cvfit$lambda), function(i) rmse(cvfit$fit.preval[,i], y))
    cvfit$measures <- data.frame(
      lambda = cvfit$lambda, 
      nzero = cvfit$nzero, 
      r = rs, 
      rsq = rsqs, 
      rmse = rmses
    )
  } else if (family=="multinomial") {
    lev_y <- levels(y)
    nlam <- length(cvfit$lambda)
    # the left-out probabilities for each fold
    cvprobs <- cvfit$fit.preval[,,1:nlam]
    # get the dominant class for each lambda and observation
    cvclasses <- apply(cvprobs, c(1,3), function(x) lev_y[which.max(x)])
    # get the summary scores or measures
    cvsummary <- adply(cvclasses, .(2), function(x) {
      dat <- data.frame(obs=y, pred=factor(x, levels=lev_y))
      ret <- caret::multiClassSummary(dat, lev=lev_y)
      ret
    })
    # save
    cvfit$probs <- cvprobs
    cvfit$classes <- cvclasses
    cvfit$measures <- cbind(lambda=cvfit$lambda, nzero=cvfit$nzero, cvsummary[,-1])
  } else if (family=="binomial") {
    lev_y <- levels(y)
    nlam <- length(cvfit$lambda)
    # the left-out probabilities for each fold
    cvprobs <- cvfit$fit.preval[,1:nlam]
    # get the dominant class for each lambda and observation
    cvclasses <- apply((cvprobs>0.5)*1, 2, factor, levels=c(0,1), labels=lev_y)
    # get the summary scores or measures
    cvsummary <- t(apply(cvclasses, 2, caret::postResample, y))
    # save
    cvfit$probs <- cvprobs
    cvfit$classes <- cvclasses
    cvfit$measures <- cbind(lambda=cvfit$lambda, nzero=cvfit$nzero, cvsummary)
    cvfit$measures <- as.data.frame(cvfit$measures)
  }
  
  cvfit$family <- family
  cvfit$y <- y
  cvfit$bestfit <- get_bestfit(cvfit, type.measure, exclude.zero)
  
  return(cvfit)
}



# Make function to run a repeated CV

run_repeated_cvglmnet <- function(X, y, k=10, family="gaussian", 
                                  type.measure='rsq', parallel=F, ezero=F, 
                                  nreps=10, cfolds=NULL, ...)
{
  # ezero: means if any fits with no coefficients should be excluded
  
  #nreps <- 10; ezero <- F; parallel <- F; family <- "gaussian"
  
  #library(cvTools)
  if (is.null(cfolds)) {
    folds <- caret::createMultiFolds(y, k=k, times=nreps)
    cfolds <- sapply(1:nreps, function(i) {
      if (k < 10) {
        if (nreps > 9) {
          names <- sprintf("Fold%i.Rep%02i", 1:k, i)
        } else {
          names <- sprintf("Fold%i.Rep%i", 1:k, i)
        }
      } else {
        if (nreps > 9) {
          names <- sprintf("Fold%02i.Rep%02i", 1:k, i)
        } else {
          name <- sprintf("Fold%02i.Rep%i", 1:k, i)
        }
      }
      vec <- vector("numeric", length(y))
      for (fi in 1:length(names)) {
        vec[-folds[[names[fi]]]] <- fi
      }
      vec
    })
  }
  
  # Predict y with model X repeated nreps
  cvrep.fits <- llply(1:nreps, function(i) {
    #cat(i,"\n")
    run_cvglmnet(X, y, foldid=cfolds[,i], exclude.zero=ezero, family=family, 
                 type.measure = type.measure, keep=T)#, ...)
  }, .parallel=parallel)
  
  # Calculate the average model fits across the 10 repeats
  res <- sapply(cvrep.fits, function(cvfit) cvfit$bestfit$val)
  mean.res <- mean(res)
  
  # Combine the predicted values
  if (family == "gaussian") {
    # Combine the predicted values
    predvals <- sapply(cvrep.fits, function(x) x$bestfit$preval)
    predvals <- rowMeans(predvals)
    # Get the residuals
    resids <- lm(y ~ predvals)$residuals
    # No predicted class
    predclass <- NULL
  } else if (family == "binomial") {
    # Combine the predicted values
    predvals <- sapply(cvrep.fits, function(x) x$bestfit$prob)
    predvals <- rowMeans(predvals)
    # Get the residuals
    resids <- glm(y ~ predvals, family=binomial(link='logit'))$residuals
    # Get the predicted class
    predclass <- factor((predvals>0.5)*1, levels=c(0,1), labels=levels(y))
  } else if (family == "multinomial") {
    # Combine the predicted values
    predvals <- laply(cvrep.fits, function(x) x$bestfit$prob)
    predvals <- apply(predvals, 2:3, mean)
    # Get the residuals
    resids <- nnet::multinom(y ~ predvals, maxit=500)$residuals
    # Get the predicted class
    predclass <- apply(predvals, 1, function(x) levels(y)[which.max(x)])
    predclass <- factor(predclass, levels=levels(y))
  }
  
  # Combine the r2 and # non-zero values across the lambdas
  lambda <- cvrep.fits[[1]]$lambda
  cvres  <- matrix(cvrep.fits[[1]]$measures[[type.measure]], nrow=1)
  cvnz   <- matrix(cvrep.fits[[1]]$nzero, nrow=1)
  for (i in 2:nreps) {
    lambda0 <- cvrep.fits[[i]]$lambda
    res0    <- cvrep.fits[[i]]$measures[[type.measure]]
    nzeros0 <- cvrep.fits[[i]]$nzero
    
    newlambda <- intersect(lambda, lambda0)
    cvres     <- rbind(cvres[,is.element(lambda, newlambda)], 
                       res0[is.element(lambda0, newlambda)])
    cvnz      <- rbind(cvnz[,is.element(lambda, newlambda)], 
                       nzeros0[is.element(lambda0, newlambda)])
    lambda    <- newlambda
  }
  # average across the repeats
  df.perlambda <- cbind(lambda=lambda, res=colMeans(cvres), nzeros=colMeans(cvnz))
  
  # Another way to get the max r^2 instead of mean.r2s
  if (ezero) {
    max.res <- max(df.perlambda[df.perlambda[,3] > 0, 2])
  } else {
    max.res <- max(df.perlambda[,2])
  }
  wmax.res <- which(df.perlambda[,2] == max.res)[1]
  if (df.perlambda[wmax.res,3] == 0) wmax.res <- min(wmax.res+1, nrow(df.perlambda))
  lam <- df.perlambda[wmax.res,1] # selected lambda
  max.nz <- df.perlambda[wmax.res,3]
  
  # Determine the coefficients from the full model
  # using the lambda with the best r2
  ## get full model fit first
  full.fit <- glmnet(X, y, family=family)#, ...)
  ## get the coefficients
  coefs0 <- coef(full.fit, s=lam)
  if (family=="multinomial") {
    intercept <- sapply(coefs0, function(x) x[1,])
    coefs <- sapply(coefs0, function(x) x[-1,])
  } else {
    intercept <- coefs0[1,]
    coefs <- coefs0[-1,]
  }
  
  list(mean.max.res=mean.res, max.res=res, 
       fitted=predvals, resids=resids, fitted.class=predclass, 
       df=df.perlambda, max.res=max.res, max.nonzero=max.nz, 
       lambda=lam, intercept=intercept, coefs=coefs, 
       rep.fits=cvrep.fits, full.fit=full.fit)
}
