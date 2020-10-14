#' Fit the log-linear model for bias factors
#'
#' This function fits the log-linear model to the data using the least squares.
#' @param otu.tab Feature table of the observed abundances.Each row is one sample and each column is one taxon. It doesn't have to be normalized. However, data will be transformed into relative abundances for analysis.
#' @param x The covariate matrix that you want to test for bias.
#' @param offset The true relative abundances matrix that you want to test as an offset. It should be of the same dimension of the otu.tab. If NULL, it is assumed that all taxa that are present in each sample have the same abundance levels within the sample.
#' @param delta The indicator matrix of whether the corresponding taxon is present in each sample. Should be of the same dimension as otu.tab. If NULL, it is taken as all taxa that have nonzero reads in a sample is designed to be present in the sample.
#' @param tol Tolerance level for matrix inversion. Default is 10^-8. Should be a small number.
#' @return coefficient estimate of the log linear model (beta.matrix); and some other variables that need to be used for inference.
#' @export
#'
LLBias.fit = function(otu.tab, x, offset=NULL, delta=NULL,
                     cols.are.otus=TRUE, tol=10^-8) {

  otu.tab=as.matrix(otu.tab)
  n.otu = ncol(otu.tab)
  n.obs = nrow(otu.tab)

  if (!is.null(delta)){
    if (nrow(delta) != n.obs | ncol(delta) != n.otu){
      stop("Delta matrix is not of the same dimension as the feature table (otu.tab).")
    }
  }

  if (!is.null(offset)){
    if (nrow(offset) != n.obs | ncol(offset) != n.otu){
      stop("True abundance (offset) matrix is not of the same dimension as the feature table (otu.tab).")
    }
  }

  if (NROW(x) != n.obs){
      stop("n.obs don't match between otu table and delta matrix")
  }
  if (NCOL(x) == 1){
    x = as.matrix(x, nrow = NROW(x), ncol = 1)
  }


  if (is.null(delta)) {
    delta=mat.or.vec(n.obs,n.otu)
    delta[ otu.tab>0 ]=1
  }else{
    otu.tab=otu.tab*delta
  }

  # Remove all samples that have only a single taxon present
  if (any(rowSums(delta) <= 1)){
    flag = which(rowSums(delta) <= 1)
    flag1 = which(rowSums(delta) == 1)
    if (length(flag1) > 0){
      print(paste0(length(flag1), "samples have only one taxon present by design."))
    }
    flag0 = which(rowSums(delta) == 0)
    if (length(flag0) > 0){
      print(paste0(length(flag0), "samples have no taxon present by design."))
    }

    otu.tab = otu.tab[-flag, ]
    delta = delta[-flag, ]
    offset = offset[-flag, ]
    x = x[-flag, ]
    n.obs = nrow(otu.tab)
    print(paste("These samples are removed from analysis as they contribute no
          information. We are left with", o.obs, "samples for analysis"))
  }

  n.beta = NCOL(x)
  n.var=n.beta*n.otu
  n.y=n.obs*n.otu
  y=rep(0,n.obs*n.otu)

  k.s=rowSums(delta)
  theta = 1/rowSums(otu.tab)* otu.tab
  theta[ delta==1 ]= log(theta[delta==1] )
  if (!is.null(offset)){
    theta[ delta==1 ] = theta[ delta==1 ] - log( offset[ delta==1 ] )
  }
  theta =theta - rowSums(theta*delta )/k.s
  theta=theta*delta

  #   setup linear system
  y=as.vector(t(theta))
  x.big=matrix(0, nrow=n.y, ncol=n.var)
  for (i in 1:n.obs) {
    p.i=diag( delta[i,] ) - ( delta[i,] %o% delta[i,] )/k.s[i]
    low=(i-1)*n.otu+1
    up=i*n.otu
    x.big[low:up,]=kronecker( p.i, x[i,,drop=FALSE] )
  }
  #
  #   solve linear system by obtaining x.big^-
  #
  x.big.ginv=g.inv(x.big, tol)
  beta=x.big.ginv %*% y
  beta.matrix = matrix(beta, ncol = n.otu)

  re = list(beta.matrix = beta.matrix)  ## End of estimation


  mu=matrix(x.big %*% beta, byrow=TRUE, nrow=n.obs, ncol=n.otu)
  resid=theta - mu

  # Calcualte the Sigma for permutation.

  sigma.resid =comp.var(clr=resid, delta=delta, k.s=k.s)

  res=list( x.big=x.big, y=y, beta.matrix = beta.matrix,
            resid=resid, mu=mu, sigma.resid = sigma.resid,
            n.otu = n.otu, n.obs = n.obs, x.big.ginv = x.big.ginv,
            data = list(otu.tab = otu.tab, x = x, offset = offset, delta = delta),
            theta = theta)
  return(res)
}


#' Hypothesis testing for the log-linear model
#'
#' This function tests a contrast of the beta coefficients in the log-linear model
#' @param mod The fitted log-linear model from LLBias.fit
#' @param C.matrix The contrast matrix to be tested. The number of columns should be n.otu \eqn{\times} n.variable (number of variable in the design matrix \eqn{x} for estimation). The algorithm will assess the testablity of the contrasts.
#' @param C.list Another way to specify the contrast matrix. It consists a list of matrices each with the same dimensionality of the \eqn{\beta} coefficient matrix. The function will decipher this list and make it into a C.matrix. Only one of C.matrix and C.list needs to be specified.
#' @param n.perm Number of permutations to be conducted. Default = 1000.
#' @param tol Tolerance level for matrix inversion. Default is 10^-8. Should be a small number.
#' @return p-value from the hypothesis testing
#' @export
#'
LLBias.test = function(mod,  C.matrix = NULL,  n.perm = 1000, tol = 1e-8, C.list = NULL){

  if ((!is.null(C.list))&(!is.null(C.matrix))){
    stop("Both C.list and C.matrix specify the contrast matrix. Only one of them needs to be specified.")
  }
  if ((is.null(C.list))&(is.null(C.matrix))){
    stop("A contrast matrix is needed for hypothesis testing")
  }


  n.otu = mod$n.otu
  n.obs = mod$n.obs
  n.beta = NCOL(mod$data$x) # x has been changed into matrix format
  n.var=n.beta*n.otu
  n.y=n.obs*n.otu
  x.big = mod$x.big
  y = mod$y
  k.s = rowSums(mod$data$delta)
  delta = mod$data$delta
  beta.est = c(mod$beta.matrix)
  x.big.ginv = mod$x.big.ginv
  mu = mod$mu
  theta = mod$theta
  resid = mod$resid
  sigma.resid = mod$sigma.resid
  beta.matrix = mod$beta.matrix
 if ((!is.null(C.matrix))){
   if (is.vector(C.matrix)){
     if (length(C.matrix) != length(beta.est)){
       stop("dimension of Contrast do not match the data. Please redefine contrast")
     }
     C.matrix = matrix(C.matrix, nrow = 1)
   }
 }

  if (!is.null(C.list)){
    if (!all(unlist(lapply(C.list, dim)) == rep(c(NROW(beta.matrix), NCOL(beta.matrix)), length(C.list)))){
      stop("Some elements in C.list has wrong dimensionality. Each of them should be of the same dimension as the dimension as the beta.matrix")
    }
    C.matrix = matrix(NA, length(C.list), NCOL(C.list[[1]])*NROW(C.list[[1]]))
    for (i in 1:length(C.list)){
      C.matrix[i, ] = as.vector(C.list[[i]])
    }
  }

  resid.null.decor= matrix(NA, n.obs, n.otu)
  sigma.null.i.half= array(data=NA, dim=c(n.obs,n.otu,n.otu))
  sigma.resid.null=  matrix(data=NA, n.otu,n.otu)

  mu.null=array(data=NA,dim=c(n.obs,n.otu))
  beta.k=rep(0, n.var)  # matrix of the vector of beta

  ## Check the testablity of the contrast matrix.

  svd.x = svd(x.big)
  R0 = svd.x$v[, abs(svd.x$d) < 1e-10]
  if (abs(norm(C.matrix %*% R0, "F")) > 1e-5){
    stop("The hypothesis is not testable. Check your contrast matrix.")
  }

  n.constraint=dim(C.matrix)[1]
  beta.k =limSolve::lsei(A=x.big, B=y, E=C.matrix, F=rep(0,n.constraint))$X

  mu.null=matrix( x.big %*% beta.k, byrow=TRUE, nrow=n.obs, ncol=n.otu)
  resid.null = theta - mu.null
  sigma.resid.null = comp.var(clr=resid.null, delta=delta, k.s=k.s)

  for (i in 1:n.obs) {
    p.i=diag(delta[i,]) - 1/k.s[i] * (delta[i,] %o% delta[i,])
    sigma.i=p.i %*% sigma.resid.null %*% p.i *(k.s[i] -1)/k.s[i]
    sigma.i.minus.half=g.inv(sigma.i,tol=tol,power=-0.5)
    resid.null.decor[i,]=resid.null[i,] %*% sigma.i.minus.half
    sigma.null.i.half[i,,]=g.inv(sigma.i,tol=tol,power=0.5)
  }
  resid.null.decor=delta*resid.null.decor


  beta.matrix=matrix(beta.est,ncol=n.otu)
  Fstat = sum(resid.null^2)/sum(resid^2) - 1

  p.F = 0

  for (np in 1:n.perm){

    # initiate the size of the matrix
    # resid.decor.perm=matrix(0, n.obs, n.otu)

    resid.decor.perm = resid.null.decor

    for (j in 1:n.obs) {
      use=which( delta[j,]>0 )
      if (length(use) == 1){
        resid.decor.perm[j,use]=resid.null.decor[j, use]
      }else{
        iperm.j=sample(use,replace= F)
        resid.decor.perm[j,use]=resid.null.decor[j, iperm.j]
      }
    }

    resid.perm =  matrix(NA, n.obs, n.otu)
    for (i in 1:n.obs) {
      resid.perm[i,]=resid.decor.perm[i,] %*% sigma.null.i.half[i,,]
    }

    theta.perm=mu.null +resid.perm
    y.perm=as.vector(t(theta.perm))
    beta.perm = x.big.ginv %*% y.perm
    mu.perm = x.big %*% beta.perm
    resid.perm = y.perm - mu.perm  ## Full model in this case.
    beta.k.perm = lsei(A = x.big, B = y.perm, E = C.matrix, F=rep(0,n.constraint))$X

    # mu.perm.null=x.reduced %*% beta.k.perm
    mu.perm.null =x.big  %*% beta.k.perm
    resid.perm.null = y.perm - mu.perm.null

    Fstat.perm = sum(resid.perm.null^2)/sum(resid.perm^2) - 1
    p.F= p.F + ifelse(Fstat >= Fstat.perm,0,1) + ifelse(Fstat == Fstat.perm,0.5,0)
  }

  p.F=p.F/n.perm

  return(p.F)
}

