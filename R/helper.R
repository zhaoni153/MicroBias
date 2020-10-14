g.inv = function(x, tol, power=-1) {
  # calcuate the genralized inverse of a matrix.
  # power: gives the
  svd.x=svd(x)
  use=(svd.x$d>tol)
  dinv.ut=(svd.x$d[use])^power *t( svd.x$u[,use] )
  x.ginv=svd.x$v[,use] %*% dinv.ut
  return(x.ginv)
}

comp.var = function( clr, delta, k.s, tol=10^-8, shrink=FALSE ) {
  # calculate the compositional variance \Sigma for use in the permutation algorithm
  n.obs=dim(clr)[1]
  n.otu=dim(clr)[2]
  clr=clr - rowSums( clr*delta )/k.s
  clr=clr*delta
  #
  #   calculate compositional mean mu
  #
  mu=rep(0,n.otu)
  p.m.cum=matrix(0,nrow=n.otu, ncol=n.otu)
  p.big=matrix(0,nrow=n.otu^2, ncol=n.otu^2)
  for (i in 1:n.obs) {
    mu=mu + clr[i,]
    p.m.i=diag( delta[i,] ) - 1/k.s[i]* delta[i,] %o% delta[i,]
    p.m.cum=p.m.cum+p.m.i
    p.big=p.big + kronecker(p.m.i, p.m.i)
  }
  p.m.ginv=g.inv( p.m.cum, tol=tol )
  mu = p.m.ginv %*% mu
  #
  #   calculate compositional variance
  #
  v.mat=matrix(0,nrow=n.otu, ncol=n.otu)
  clr.ctr=matrix(0,nrow=n.obs, ncol=n.otu)
  for (i in 1:n.obs) {
    p.m.i=diag( delta[i,] ) - 1/k.s[i]* delta[i,] %o% delta[i,]
    clr.ctr[i,]=clr[i,] - p.m.i %*% mu
    #        v.mat=v.mat + clr.ctr[i,] %o% clr.ctr[i,]
  }

  if (shrink==TRUE) {
    sigma=sqrt( colMeans( clr.ctr^2 ) - colMeans(clr.ctr)^2 )
    target=sigma*p.m.cum/n.obs
    target=t( sigma*t(target) )
    v.mat=n.obs*CovEst.2003LW(clr.ctr, target=target)$S
  }
  else {
    v.mat=t( clr.ctr ) %*% clr.ctr
  }
  v.vec=as.vector(v.mat)
  p.big.ginv=g.inv( p.big, tol=tol )
  sigma.vector=p.big.ginv %*% v.vec
  sigma=matrix(sigma.vector, nrow=n.otu)
  #    res=list(mu=mu, p.m.cum=p.m.cum, p.big=p.big, v.mat=v.mat, sigma=sigma, clr.ctr=clr.ctr, clr=clr, delta=delta)
  return(sigma)
}




























