#' @export
#' @importFrom stats rbinom runif

simConf <- function(n,nconf,dconf = c('beta',2,2),
                    nexpred, dexpred = c('beta',2,2),
                    noutpred, doutpred = c('beta',2,2),
                    effect = 2, escale = 'OR',
                    brate,erate){
  set.seed(12345)
  # Define random variable generators
  ncov <- nconf + nexpred + noutpred
  genConf <- paste('r',dconf[1],"(n,",dconf[2],",",dconf[3],")",sep='')
  genExPred <- paste('r',dexpred[1],"(n,",dexpred[2],",",dexpred[3],")",sep='')
  genOutPred <- paste('r',doutpred[1],"(n,",doutpred[2],",",doutpred[3],")",sep='')

  dataConf <- as.data.frame(replicate(nconf, eval(parse(text = genConf ))))
  dataExPred <- as.data.frame(replicate(nexpred, eval(parse(text = genExPred ))))
  dataOutPred <- as.data.frame(replicate(noutpred, eval(parse(text = genOutPred ))))

  covData <- data.frame(dataConf,dataExPred,dataOutPred)
  #X <- rbind(covData[c(1:(n/2)),]-0.5,covData[c(((n/2)+1):n),]+0.5);dim(X)
  X <- covData
  colnames(X) <- paste("x", 1:ncov, sep = "")
  #covCoef <- runif(ncov,-.5,.5)
  #covCoef <- round(runif(ncov,.2,.7),2)*sample(c(-1,1),ncov,replace=T)
  nconfpart <- round(nconf/5,0)
  nconfremain <- nconf-nconfpart

  covCoef <- c(round(runif(nconfpart,.6,.7),2),
               round(runif(nconfremain,.1,.2),2)*sample(c(-1,1),nconfremain,replace=T),
               round(runif(nexpred,.1,.2)*(-1),2),
               round(runif(noutpred,.5,.6),2))
  expredCols <- (nconf+1):(nconf+nexpred)
  outpredCols <- (nconf+nexpred+1):(ncov)

  ## Generate potential outcomes ##
  Xout <- X[,-expredCols]
  const.out <- log(brate/(1-brate))
  #logOR1 <- const.out + as.numeric(effect)*1 + t(t(as.vector(covCoef[-expredCols])) %*% t(Xout)  )
  logOR1 <- const.out + as.numeric(log(effect))*1 + t(t(as.vector(covCoef[-expredCols])) %*% t(Xout))
  prob.logOR1 = exp(logOR1)/(1+exp(logOR1))
  outcome1 = rbinom(n,1,prob.logOR1)
  #prop.table(table(outcome1))

  logOR0 = logOR1 - as.numeric(effect)*1
  prob.logOR0 = exp(logOR0)/(1+exp(logOR0))
  outcome0 = rbinom(n,1,prob.logOR0)
  #prop.table(table(outcome0))


  ## Generate treatment ##
  Xtrt <- X[,-outpredCols]
  const.ex <- log(erate/(1-erate))
  logExp = const.ex + t(t(as.vector(covCoef[-outpredCols])) %*% t(Xtrt)  )
  prob.logExp = exp(logExp)/(1+exp(logExp)) # probability to get exposed
  exposure = rbinom(n,1,prob.logExp); ## Generate exposure

  obs.outcome = exposure*outcome1 + (1-exposure)*outcome0

  ## Generate data
  dat = cbind(outcome1, outcome0, obs.outcome, exposure, X)
  dat
}
