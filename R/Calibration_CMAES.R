#' write the configuration for CMA-ES method
#' \code{write.cmaes} 
#' @param x Configuration of CMA-ES calibration method.
#' @param file full path of file
#' @param backup TRUE/FALSE
#' @return Parameters for CMA-ES configuration
#' @export
write.cmaes <- function(x=NULL, file, backup = TRUE){
  if(is.null(x)){
    cn = c('lambda',
           'ncores',
           'maxgen',
           'stopfitness',
           'sigma',
           'updateic',
           'nspingup',
           'walltime',
           'tsd_obs',
           'path_out'
    )
    x=cbind(48, 48, 100, 0.3, 0.8, 0, 0, 86400)
    x = data.frame(x)
    x=cbind(x, 'xxx', 'cmaes_out')
    names(x)=cn
    write.cmaes(x, file)
  }else{
    write.config(x, file=file, backup = backup)
  }
  return(x)
}
#' Automatically calibrate model, that requires user-defined model and objective function.
#' \code{CMAES}
#' @param CV  Input data, list of data.
#' @param cmd Command to run the model.
#' @param objfunc  User-defined objective function which return the objective values.
#' @param Call_Model Function that calls model simulation.
#' @param lambda Number of children in each generation
#' @param maxstep Maximum generations
#' @param ncores Number of cores to simulate. 1 = one thread.
#' @param sigma Sigma Value to sample (0, 1)
#' @param stopfitness The optimal value. When the objective value is smaller than stopfitness, calibration success.
#' @param debug Whether debug Model. 
#' @param ... More options passing to objfunc.
#' @importFrom doParallel registerDoParallel
#' @export
CMAES<- function (CV, cmd, objfunc,  Call_Model, 
                  lambda = CV$method$LAMBDA,
                  maxstep = CV$method$MAXGEN,
                  # ncores = max(CV$method$NCORES, 1),
                  sigma = CV$method$SIGMA,
                  stopfitness=CV$method$STOPFITNESS,
                  debug=FALSE,
                  ...){
  #This function is modified from CRAN package 'cmaes'.
  # stopifnot(is.numeric(par))
  stopifnot(!is.null(objfunc))
  # if( !file.exists(cmd) ){
  # if(system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)!=0){
  #   stop(paste(cmd, "does not exist.\n"))
  # }
  dir.out = as.character(type.convert(CV$PATH_OUT))
  range = CV$range
  message('CMAES::Output path: ', dir.out)
  dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)
  # if(ncores > 1){
  #     doParallel::registerDoParallel(ncores)
  # }
  para.name = names(CV$range)
  para.id = which(CV$range[1, ]!=0)
  N = sum(length(para.id)) # number of parameters to calibration.
  stopifnot( N > 0)
  
  message('CMAES::Calibration on parameters: ')
  print(CV$range[,para.id])
  
  if(is.null(lambda) ){ lambda = 4 + floor(3 * log(N))}
  if (sigma < 0.1 || sigma > 0.9){
    sigma <- max(min(sigma, 0.9), 0.1) }
  sigma = rep(sigma, N)
  # lower = as.numeric(range['min', para.id])
  # upper = as.numeric(range['max', para.id] )
  lower = rep(-1, N)
  upper = rep(1, N)
  
  xmean <- rep(0, N)
  
  mu <- lambda/2
  weights <- log(mu + 0.5) - log(1:mu)
  mu <- floor(mu)
  weights <- weights/sum(weights)
  mueff <- sum(weights)^2/sum(weights^2)
  cc <- (4 + mueff/N)/(N + 4 + 2 * mueff/N)
  cs <- (mueff + 2)/(N + mueff + 5)
  c1 <- 2/((N + 1.3)^2 + mueff)
  cmu <- min(1 - c1, 2 * (mueff - 2 + 1/mueff)/((N + 2)^2 +  mueff))
  damps <- 1 + 2 * max(0, sqrt((mueff - 1)/(N + 1)) - 1) + cs
  pc <- ps <- numeric(N)
  B <- diag(N)
  D <- rep(1, N)
  C <- as.matrix(B) %*% diag(D^2) %*% B
  invsqrtC <- as.matrix(B) %*% diag(D^-1) %*% B
  eigeneval <- 0
  chiN <- N^0.5 * (1 - 1/(4 * N) + 1/(21 * N^2))
  ml.triu <- function(M, k = 0) {
    if (k == 0)
      M[lower.tri(M, diag = FALSE)] <- 0
    else M[col(M) <= row(M) + k - 1] <- 0
    return(M)
  }
  iGen <- 1
  arr = array(0, dim=c(N, lambda, maxstep))
  BestOBJ = numeric(maxstep)
  arx <- matrix(0, nrow = N, ncol = lambda)
  colnames(arx) = paste0('sample', 1:lambda)
  rownames(arx) = rownames(para.name[para.id])
  arfitness <- numeric(lambda)
  
  write.config(CV$range, file=file.path(dir.out, paste0(CV$prjname,'.range.txt')), backup = FALSE)
  for(iGen in 1:maxstep) {
    message('\n\n')
    message('\t========================')
    message('CMAES::', iGen, '/', maxstep)
    message('\t========================')
    arx = arx * 0;
    for (k in 1:lambda) {
      arxk <- xmean + sigma * B %*% (D * rnorm(N, sd=sigma))
      arxk <- ifelse(arxk > lower, ifelse(arxk < upper, arxk, upper),lower)
      arx[, k] <- arxk
    }
    message('CMAES::',' sigma = ') ; print(sigma); message('')
    
    # Range + arx to .calib files
    # run the model  and get GOF back.
    # debug(Call_Model)
    xout = Call_Model(iGen=iGen, pop = arx, #ncores=ncores, 
                      CV=CV, CMD.EXE = cmd, objfunc=objfunc, 
                      debug=debug,
                      ...)
    
    # print(xout$fitness)
    arfitness = abs(CV$method$BESTGOF - xout$fitness)
    message('CMAES::Fitness at Geneartion ', iGen, '/', maxstep, ':')
    print(arfitness)
    
    arr[,,iGen] = arx;
    SortID <- order(arfitness, decreasing = FALSE)
    arfitness <- arfitness[SortID]
    iBest = SortID[1]
    calibmat = xout$varlist$calibmat
    bestcalib = calibmat[SortID[1], ]
    
    gof.tab = cbind(iGen, SortID[1:5], arfitness[1:5]);
    colnames(gof.tab) = c('Gen', 'Job', 'GOF')
    if(iGen ==1){
      write.table(gof.tab, file = file.path(CV$PATH_OUT, 'gof.csv'),
                  append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }else{
      write.table(gof.tab, file = file.path(CV$PATH_OUT, 'gof.csv'),
                  append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    message('================================')
    message('CMAES::','\n\nOBJECTIVE VALUES in Genearation ', iGen)
    print(gof.tab)
    message('================================')
    # CV$calib = bestcalib # UPDATE THE CALIB into the CV
    write.config(bestcalib, file=file.path(dir.out, paste0('calib_Gen.', iGen, '.calib') ), backup = FALSE)
    # vlist = pre.files(iGen = iGen, pop = pop, CV=CV);
    # Obj.Func(jobid = SortID[1], CV = CV, vlist = )
    # if( arfitness < 1) { #if the simulation is good enough.
    #   fun.updateinit(CV, iBest)
    # }
    xold <- xmean
    xmean <- arx[, SortID[1:mu]] %*% weights
    ps <- (1 - cs) * ps + sqrt(cs * (2 - cs) * mueff) * invsqrtC %*% (xmean - xold)/sigma
    hsig <- norm(ps, "F")/sqrt(1 - (1 - cs)^(2 * iGen/lambda))/chiN < 1.4 + 2/(N + 1)
    pc <- (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mueff) * (xmean - xold)/sigma
    artmp <- (1/sigma) * (arx[, SortID[1:mu]] - matrix(1, 1, mu) %x% xold)
    C <- (1 - c1 - cmu) * C + c1 * (pc %*% t(pc) + (1 - hsig) * cc * (2 - cc) * C) + cmu * artmp %*% diag(weights) %*% t(artmp)
    sigma <- sigma * exp((cs/damps) * (norm(ps, "F")/chiN - 1))
    if (iGen - eigeneval > lambda/(c1 + cmu)/N/10) {
      eigeneval <- iGen
      C <- ml.triu(C) + t(ml.triu(C, 1))
      if (any(is.nan(C)))
        break
      EigVal <- eigen(C, symmetric = TRUE)
      B <- EigVal$vectors
      D <- sqrt(EigVal$values)
      invsqrtC <- B %*% diag(D^-1) %*% t(B)
    }
    BestOBJ[iGen] = arfitness[1]
    if (arfitness[1] <= stopfitness || max(D) > 1e+07 * min(D)){
      message('CMAES::','The best fitness (', arfitness[1], ') is less than threshold (', stopfitness, ').')
      message('CMAES::','Current Generation = ', iGen)
      message('CMAES::','Best Calib:')
      print(bestcalib)
      break
    }
    if( diff(range(arfitness, na.rm=TRUE)) / mean(arfitness, na.rm=TRUE) < 0.01){
      message('CMAES::','The POSSIBLE best fitness is reached: ', arfitness[1], '.')
      message('CMAES::','Current Generation = ', iGen)
      message('CMAES::','Best Calib:')
      print(bestcalib)
      break
    }
  }
  xmin <- arx[, SortID[1]]
  ymin <- calibmat[, SortID[1]]
  return(list(xmin = xmin,
              ymin = ymin,
              BestOBJ  = BestOBJ[1:iGen])
  )
}
