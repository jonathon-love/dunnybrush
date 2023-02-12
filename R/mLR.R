
#' @importFrom MASS ginv
mLR <- function(data=NULL, weights=NULL, design=NULL, nstep=Inf, parallel=NULL, init=NULL, tol=NULL, input=NULL) {
  # solves the monotonic Linear Regression problem
  # data is a n-vector of observations
  # w is an n x n matrix of weights (default=identity matrix)
  # a is a design matrix (e.g., for main effect, additive model, etc.)
  # (default=intercept model)
  # nstep is number of topes to search (default=exhaustive search)
  # parallel is optional list of parallelism classes of the difference matrix of a
  # init is optional initial tope for search
  # tol is tolerance (absolute values less than tol set to zero)
  # input is output from a previous call to mLR
  # returns optimal fit, fitted values (predictions), and tope (sign vector).
  #
  # load packages if required
  #if(require("lpSolveAPI")==F){install.packages ("lpSolveAPI"); library(lpSolveAPI)}
  #if(require("nloptr")==F){install.packages ("nloptr"); library(nloptr)}
  #if(require("MASS")==F){install.packages ("MASS"); library(MASS)}
  #if(require("quadprog")==F){install.packages ("quadprog"); library(quadprog)}

  if (is.null(nstep)){nstep <- Inf}
  if(!is.null(input)){
    y <- input$data; pred <- input$predictions; w <- input$weights; a <- input$design;
    xhat <- sign(mat2diff(input$permutation)); open <- input$open; closed <- input$closed; itr <- input$nstep; minitr <- input$solstep;
    type <- input$type; minfit <- input$fit; parallel <- input$parallel; nstep <- itr+nstep; tol <- input$tol
  } else {
    y <- data; w <- weights; a <- design
    if (is.data.frame(y)) {y = y[,]} # convert to vector
    n <- length(y)
    if (is.null(w)){w <- diag(rep(1,n))} # default identity matrix
    if (is.data.frame(w)) {w=data.matrix(w)} # convert to matrix
    if (is.vector(w)) {w <- diag(w)} # convert vector to matrix form
    if (is.null(a)){a <- matrix(rep(1,n),n,1)} # default intercept model
    if (is.data.frame(a)) {a=data.matrix(a)} # convert to matrix
    if (is.null(tol)) {tol <- 1e-5} # set tolerance for values nearly zero
    open <- NULL; closed <- NULL; parallel <- NULL
    itr <- 0; minitr <- 0 # initialize step count and solution step
  }
# initialize some matrices etc.
    dy = mat2diff(y) # difference vector of y
    da <- mat2diff(a,1) # difference matrix of a
    if(is.vector(a)){a<-as.matrix(a)}
# deal with special cases
    r <- qr(da)$rank
    if(r==0){
      type <- 'Data is in DA'
      xhat <- rep(0,nrow(da)) # xhat is null vector
      mr <- lsqisotonic1(y,w,xhat)
      pred <- mr$pred; minfit <- mr$fit;
    }
    else if(r==1){
      type <- 'DA rank 1'
      k <- which.max(colSums(abs(da))); t <- da[,k]
      mr1 <- lsqisotonic1(y,w,t)
      mr2 <- lsqisotonic1(y,w,-t)
      if(mr2$fit < mr1$fit){
        minfit <- mr2$fit; xhat <- -t; pred <- mr2$pred
      } else {
        minfit <- mr1$fit; xhat <- t; pred <- mr1$pred
      }
    }
    else if(!is.null(solveLP(dy,da))){
      type <- 'Data is in DA'
      xhat <- sign(dy); minfit <- 0; pred <- y;
    }
    else { # search topes for a solution
      type <- 'Search'
      cond <- condense(a) # condense a
      dac <- mat2diff(cond$matrix,1) # difference matrix of condensation of a
      dc <- mat2diff(diag(nrow(cond$matrix))) # difference matrix of order nrow(ac)
      if (is.null(parallel)){parallel <- parallelclass(dac)} # calculate parallelism classes of dac
      rankdc <- rep(0,length(parallel))
      for (i in 1:length(parallel)){rankdc[i]<-qr(dc[parallel[[i]],])$rank} # store dc sub-matrix ranks
      if(is.null(input)){
        if (is.null(init)){ # calculate initial tope
          yp <- a%*%ginv(t(a)%*%w%*%a)%*%t(a)%*%w%*%y # weighted projection of y onto col(a)
          dyp <- mat2diff(yp[cond$id]) # difference matrix of condensation of yp
          dyp[abs(dyp)<tol] <- 0 # set small differences to zero
          init <- as.vector(sign(dyp)) # initial tope
        }
        tc <- init
        # check if tc contains zeros and is therefore not a tope
        if (any(tc==0)){
          uc <- maxsignaugment(tc) # find best initial tope
          bestfit <-Inf
          for (i in nrow(uc)){
            x <- solveLP(uc[i,],dac)
            if (!is.null(x)){
              u <- decondense(uc[i,],cond); mr <- lsqisotonic1(y,w,u)
              if (mr$fit < bestfit){
                bestfit <- mr$fit; tc <- as.vector(uc[i,])
              }
            }
          }
        }
        t <- decondense(tc, cond); mr <- lsqisotonic1(y,w,t) # # decondense initial tope and fit
        open <- list() # create open list
        open[[1]] <- list("fit"=mr$fit, "pred"=mr$pred, "tope"=tc) # initialize open
        closed <- t(as.matrix(tope2perm(tc))) # add initial tope to closed (as permutation)
        minfit <- Inf; xhat <- c() # initialize optimal fit and solution tope
      }
      H <- dac%*%ginv(t(dac)%*%dac)%*%t(dac) # projection matrix on to col(dac)
      n <- nrow(cond$matrix); nc <- n*(n+1)*(2*n+1)/6 # sum of permutation elements squared

      # start search
      while (length(open)>0 && itr < nstep && minfit > 0) {
        itr <- itr + 1 # increment step count
        tc <- open[[1]]$tope; f <- open[[1]]$fit; yp <- open[[1]]$pred # retrieve values from open
        open[[1]] <- NULL # remove first element from open
        if (minfit-f>tol){minfit <- f; xhat <- decondense(t,cond); pred <- yp; minitr <- itr} # update if better fit
        else if(minfit-f==tol){xhat <- rbind(xhat,decondense(t,cond))}
        # find adjacent topes
        ds <- abs(dc%*%t(dc)%*%tc)/2 # indicator of adjacent topes in dc
        for (i in 1:length(parallel)) {
          s <- parallel[[i]]; j <- ds[s]==1
          if (sum(j)==rankdc[i]){ # parallelism class defines a face of the cone of tc in col(d)
            xc <- tc; xc[s] <- -xc[s]; # create candidate adjacent tope
            pxc <- tope2perm(xc)
            closed_xc = closed%*%pxc # check if in closed and therefore visited before
            if (max(closed_xc) < nc) { # not in closed so continue
              q <- as.vector(xc%*%H); q[abs(q)<tol] <- 0 # projection test
              if (all(xc==sign(q))) { # xc passes projection test
                tope <- 1
              } else {
                if (is.null(solveLP(xc,dac))) {tope <- 0} else {tope <- 1} # solve LP problem
              }
              if (tope==1) { # xc is an adjacent tope
                x <- decondense(xc,cond); mr <- lsqisotonic1(y,w,x) # calculate fit
                open[[length(open)+1]] <- list("fit"=mr$fit, "pred"=mr$pred, "tope"=xc); # add adjacent tope to open
              }
              closed <- rbind(closed,pxc) # add xc to closed as permutation
            }
          }
        }
        if (length(open) > 0) {
          f <- sapply(open, function(x) {x[[1]]}) # create vector of fit values in open
          open <- open[order(f)] # re-order open by fit
        }
      }
    }

  # arrange output
    if (length(pred)==0) {
      mr <- lsqisotonic1(y,w,xhat)
      pred <- mr$pred; minfit <- mr$fit
    }
  if (qr(da)$rank > 0){x <- solveLP(xhat,da)} else {x <- rep(0,ncol(a))}
  if (!is.null(x)){xp <- as.vector(a%*%x)} else {xp <- NULL}
  xhat <- as.vector(xhat)
  rownames(closed) <- NULL
  output <- list("fit"=minfit, "data"=y, "predicted"=pred, "permutation"=tope2perm(xhat), "solution"=x, "estimate"=xp,
                 "nstep"=itr, "solstep"=minitr, "weights"=w, "design"=a, "parallel"=parallel, "tol"=tol,
                 "open"=open, "closed"=closed, "type"=type)
  return(output)
}

factdesign <- function(levels=NULL, addflag=0){
  # returns design matrix for factorial combination defined by levels,
  # a vector consisting of the number of levels of each factor
  # if addflag is 0 then it returns additive model
  # otherwise it returns all the interaction terms as well
  nfact <- length(levels)
  c <- 1:prod(levels)
  d <- c(rep(1,nfact))
  for (i in 1:nfact-1) {
    d[i] <- prod(levels[(i+1):length(levels)])
  }
  b <- matrix(0,length(c),nfact)
  for (i in 1:length(c)) {
    j <- i
    for (k in 1:nfact) {
      z <- floor((j-1)/d[k]); b[i,k] <- z+1
      j <- j - z*d[k]
    }
  }
  m <- vector(mode = "list", length = nfact)
  for (k in 1:nfact) {
    z <- diag(levels[k])
    m[[k]] <- z[,-1]
    if (!is.matrix(m[[k]])) {m[[k]] = matrix(m[[k]],length(m[[k]]),1)}
  }
  A <- matrix(); A <- A[-1]
  for (i in 1:nrow(b)) {
    aa <- c()
    for (k in 1:ncol(b)) {
      aa <- c(aa, m[[k]][b[i,k],]);
    }
    A <- rbind(A, aa)
  }
  rownames(A) <- NULL # remove annoying row names

  if (addflag != 0) {
    # add interaction terms if asked for
    for (ifact in 2:nfact) {
      k = combn(nfact,ifact); k = t(k)
      for (j in 1:nrow(k)) {
        b <- vector(mode = "list", length = ncol(k))
        for (i in 1:ncol(k)) {
          jj <- k[j,i]
          if (jj==1) {
            i1 <- 1
          }
          else {
            i1 <- sum(levels[1:jj-1]-1)+1
          }
          i2 <- i1 + levels[jj] - 2 # delimits main effect block
          aa = A[,i1:i2]
          if (is.matrix(aa)) {b[[i]] <- aa}
          else {b[[i]] <- matrix(aa,length(aa),1) }
        }
        B <- b[[1]]
        for (i in 2:length(b)) {
          bb <- matrix(); bb <- bb[-1]
          for (ib in 1:ncol(b[[i]])) {
            bb <- cbind(bb,B*b[[i]][,ib])
          }
          B <- bb
        }
        A <- cbind(A,B)
      }
    }
  }
  A <- cbind(c(rep(1,nrow(A))),A) # add intercept
  return(A)
}

mat2diff <- function(A=NULL, flag=0) {
  # returns difference matrix DA = D*A
  # where A is a matrix - usually a design matrix
  # if flag == 0 then all null columns are removed from DA
  A <- as.matrix(A)
  n <- nrow(A);
  if (n > 1) {
    DA <- matrix(0,choose(n,2),ncol(A))
    k <- 0
    for (i in 1:(n-1)) {
      for (j in (i+1):n){
        k <- k + 1
        DA[k,] <- A[i,]-A[j,]
      }
    }
    # remove null columns
    if (flag == 0) {
      if (ncol(DA) > 1){
        s <- apply(DA != 0, 2, sum); s <- which(s==0)
        if (length(s) > 0) {
          DA <- DA[,-s]
          if (!is.matrix(DA)) {DA=matrix(DA,length(DA),1)}
        }
      }
    }
  }
  else {DA <- matrix(0,1,length(A)) }
  return (DA)
}

lsqisotonic1 <- function(y=NULL, w=NULL, t=NULL){
  # calculates isotonic regression on vector y with weights w
  # t is a sign vector specifying the required order of fitted values
  # output is fitted values ($pred) and weighted least-squares fit ($fit)
  #
  # based on Matlab function lsqisotonic
  # Copyright 2003-2004 The MathWorks, Inc.
  # Revision: 1.1.6.3 Date: 2004/02/01 22:10:40

  n <- length(y)
  if (is.null(w)){w <- diag(rep(1,n))}
  if (!is.null(t)){
    d <- mat2diff(diag(1,n))
    x <- as.vector(t)%*%d

    # Sort points ascending in x, break ties with y.
    # force secondary approach to ties
    iord = 1:n
    xy = cbind(t(x),-y,iord)
    xyord <- xy[order(xy[,1],xy[,2]),]
    ord = xyord[,3]; iord[ord] <- 1:n
    xyord[,2]=-xyord[,2]

    # Initialize fitted values to the given values.
    yhat <- xyord[,2]
    # convert vector of weights to matrix form or matrix of weights to vector
    # form
    if (is.null(w)){w <- diag(rep(1,length(y)))} # default w = identity matrix
    if (is.vector(w)){
      W <- diag(w)
    } else {
      W <- w; w <- diag(W)
    }
    block <- 1:n
    w <- w[ord]# reorder w as a column

    # Merge zero-weight points with preceding pos-weighted point (or
    # with the following pos-weighted point if at start).

    posWgts <- (w > 0)
    if (any(!posWgts)){
      idx <- cumsum(posWgts); idx[idx == 0] <- 1
      w <- w[posWgts]
      yhat <- yhat[posWgts]
      block <- idx[block]
    }

    diffs <- diff(yhat)
    while (any(diffs < 0)){
      # If all blocks are monotonic, then we're done.
      # Otherwise, merge blocks of non-increasing fitted values, and set the
      # fitted value within each block equal to a constant, the weighted mean
      # of values in that block.
      idx <- cumsum(c(1,diffs>0))
      sumyhat <- tapply(w*yhat,idx,sum)
      w <- tapply(w,idx,sum)
      yhat <- sumyhat / w
      block <- idx[block]
      diffs <- diff(yhat)
    }

    # Broadcast merged blocks out to original points, and put back in original order.
    yhat <- yhat[block]; yhat <- yhat[iord]
    # calculate fit
    v <- as.vector(y - yhat); names(yhat) <- NULL
    fit <- v%*%W%*%v; fit <- as.numeric(fit)
  } else {
    out <- NULL
  }
  # store in out
  out <- list("pred" = yhat, "fit" = fit)
  return (out)
}

d2t <- function(x=NULL, n=NULL){
  # Converts a decimal number x into a trinary array of (optional) length n
  # Based on the Matlab function d2b.m
  t <- 3
  if (is.null(n)){n <- ceiling(log(max(x))/log(t)) + 1}
  y <- matrix(0,length(x),n)
  for (j in 1:length(x)){
    z <- x[j]
    if (z < 0){
      y[j,n] <- -1
    } else if (z ==0){
      y[j,n] <- 0
    } else if (z ==1){
      y[j,n] <- 1
    } else {
      c <- ceiling(log(z)/log(t)) + 1 # Number of divisions necessary ( rounding up the log2(x) )
      yy <- rep(0,c) # Initialize output array
      for (i in 1:c){
        r <- floor(z / t)
        yy[c+1-i] <- z - t*r
        z <- r
      }
      y[j,(n-c+1):n] <- yy
    }
  }
  if (all(y[,1]==0)){y <- y[,-1]} # delete leading zeros if not required
  return(y)
}

maxsignaugment <- function(v){
  # v is a m x n matrix (usually of signs)
  # for each row X of v, find the set of topes w
  # such that for each Y in w, Y o X = X
  # return the set of unique topes so formed
if (is.vector(v)) {v <- t(as.matrix(v))}
if (ncol(v)==1) {v <- t(v)}
s <- matrix(); s <- s[-1]
for (i in 1:nrow(v)){
    u <- matrix(v[i,],1,ncol(v))
    k <- which(u==0); n <- length(k)
    if (n > 0){
        r <- d2t(0:(3^n-1)); r[r==2] <- -1
        p <- kronecker(matrix(1,nrow(as.matrix(r))),u)
        p[,k] <- r
        s <- rbind(s,p)
    } else {
        s <- rbind(s,u)
    }
}
m <- apply(s!=0,1,sum)==ncol(s) # find sign vectors that are topes
s <- s[m,] # keep them
s <- s[!duplicated(s), ] # remove duplicated sign vectors
return(s)
}

parallelclass <- function(a){
  # returns list of parallel classes of matrix a
  out <- simplematrix(a)
  p = list()
  for (i in 1:length(out$id)){
    p = append(p,list(which(out$index==out$id[i])))
  }
  return(p)
}

simplematrix <- function(a){
  # simplifies matrix a by removing null and parallel rows
  # s is the simple matrix
  # ia is a list of numbers of rows of a included in s
  # ic is a list of numbers of rows of s corresponding to each row of a
  # if any ic==0 then the corresponding row of a is null
  #
  ia <- matrix(); ia <- ia[-1]; ic <- rep(0,nrow(a))
  j <- apply(a!=0,1,sum); ic[j==0] <- -1
  # find parallel rows
  for (i in 1:nrow(a)){
    if (ic[i]==0){
      ic[i] <- i; ia <- c(ia,i)
      if (i < nrow(a)){
        for (j in (i+1):nrow(a)){
          if (qr(a[c(i,j),])$rank==1 && ic[j]==0){
            ic[j] <- i
          }
        }
      }
    }
  }
  s <- a[ia,]
  ic[ic<0] <- 0
  out = list("matrix"=s,"id"=ia,"index"=ic)
  return(out)
}

condense <- function(a){
  # returns condensed version of a (redundant rows removed)
  # id is vector of retained rows of a
  # index is vector of rows of a indexed by id
  #
  a <- as.matrix(a)
  ic <- rep(0,nrow(a)); id <- ic; k <- 0
  # find identical rows
  for (i in 1:nrow(a)){
    if (ic[i]==0){
      k <- k + 1; ic[i] <- k; id[k] <- i
      if (i < nrow(a)){
        for (j in (i+1):nrow(a)){
          if (all(a[i,]==a[j,]) & ic[j]==0){ic[j] <- k}
        }
      }
    }
  }
  if (k < length(id)){j = (k+1):length(id); id <- id[-j]}
  ac <- a[id,]
  out <- list('matrix'= ac,'id'=id,'index'=ic)
  return (out)
}

decondense <- function(x=NULL, cond=NULL){
  # decondenses a condensed sign difference vector according to cond
  # cond is output from condense(a) where a is the design matrix
  if(is.null(x)){cat('Error: Tope not specified.')} else {
    ic <- cond$index; u <- unique(ic)
    if(length(ic)>length(u)){
      if(length(ic)==1){y<-rep(x,length(ic))}
      else {
        p <- tope2perm(x)
        p <- p[ic]
        y <- sign(mat2diff(p))
      }
    } else {y <- x}
    y = as.vector(y)
    return(y)
  }
}

#' @importFrom lpSolveAPI make.lp set.objfn get.solutioncount get.variables add.constraint set.bounds
solveLP <- function(y=NULL, a=NULL){
  # determines if sign(y) is a covector of matrix a
  # solves LP problem
  # if no solution then returns NULL

  n <- nrow(a); m <- ncol(a)
  lprec <- make.lp(n, m)
  set.objfn(lprec, rep(0,m))
  C <- a
  for (i in 1:nrow(C)) {
    if (y[i]>0) {add.constraint(lprec, C[i,], ">=", 1)}
    else if (y[i]<0) {add.constraint(lprec, C[i,], "<=", -1)}
    else {add.constraint(lprec, C[i,], "=", 0)}
  }
  set.bounds(lprec, lower = rep(-Inf,m), upper = rep(Inf,m),columns = 1:m)
  solve(lprec)
  exitflag <- 1; if(get.solutioncount(lprec)==0){exitflag <- 0}
  if(exitflag==1){x <- get.variables(lprec)}else{x<-NULL}
  return(x)
}

tope2perm <- function(x=NULL){
  # converts tope x to corresponding permutation
  n <- (1+sqrt(1+8*length(x)))/2
  d <- mat2diff(diag(n))
  p <- t(d)%*%x; p[] <- (p+n+1)/2
  p <- as.vector(p)
  return(p)
}
