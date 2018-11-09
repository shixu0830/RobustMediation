################################################################
### Necessary functions for GENIUS method:
### Isabel Fulcher (isabelfulcher@g.harvard.edu): updated functions in the
### GENIUS package to incorporate C for single IV models only 
### (and only for continuous cases for A,Y -- can easily update though)
### Xu Shi (xushi@hsph.harvard.edu)
################################################################
genius_addY_withC <- function(Y,A,G,C,formula=A~G+C,alpha=0.05,lower=-10,upper=10) {
  if (is.data.frame(Y)) {
    Y=as.vector(data.matrix(Y));
  }
  if (is.data.frame(A)) {
    A=as.vector(data.matrix(A));
  }
  if (is.data.frame(G)) {
    G=data.matrix(G)
  }
  if (is.data.frame(C)) {
    C=data.matrix(C)
  }
  if (class(C) == "matrix"){
    #number of covariates
    nC =dim(C)[2];
    if (nC==1) {C=as.vector(C);}
  } else {
    nC =1;
  }
  if (class(G) == "matrix") {
    #number of IVs
    nIV =dim(G)[2];
    #sample size
    N = dim(G)[1];
    if (nIV==1) {G=as.vector(G);}
  } else {
    nIV =1;
    N = length(G);
  }
  A.binary = all(A %in% 0:1);
  
  if (nIV>1) {
    
    if (A.binary) {
      glm.AG = stats::glm(formula(formula), family=binomial, x=T, data=cbind(A,as.data.frame(G),C));  ###XS:add C
    } else {
      glm.AG = stats::lm(formula(formula), x=T, data=cbind(A,as.data.frame(G),C)); ###XS:add C
    }
    #Multiple-IV GMM
    gmm.fun <- function(beta, X) {
      sweep(X,2,apply(X,2,mean),"-")*(A-glm.AG$fit)*(Y-beta*A)
    }
    
    gmm.out = gmm::gmm(gmm.fun,x=G,t0=c(lower,upper),optfct="optimize",
                       type="iterative",wmatrix="optimal",vcov="iid",centeredVcov=TRUE);
    
    beta.est=gmm.out$coef;
    
    if (A.binary) {
      mm= mm.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      M = M.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      B = B.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
    } else {
      mm= mm.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      M = M.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      B = B.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
    }
    
    beta.var= diag((1/N)*(solve(B)%*%M)%*%mm%*%t((solve(B)%*%M)))[nIV+length(glm.AG$coef)+1];
  } else {
    if (A.binary) {
      glm.AG = stats::glm(formula(formula), family=binomial, x=T, data=cbind(A,as.data.frame(G),C)); ###XS:add C
    } else {
      glm.AG = stats::lm(formula(formula), x=T, data=cbind(A,as.data.frame(G),C)); ###XS:add C
    }
    
    if (nC > 0){
      #call the formula (possibility to make more flexible)
      formula_GC=G~C
      
      #single-IV estimator with covariates
      #glm.GC = stats::glm(formula(formula_GC), family=binomial, x=T, data=cbind(A,as.data.frame(G)));  ###XS: change data
      glm.GC = stats::glm(formula(formula_GC), family=binomial, x=T, data=cbind(G,as.data.frame(C)));
      
      #single-IV estimator with covariates
      beta.est = mean((G-glm.GC$fit)* (A-glm.AG$fit)*Y)/mean((G-glm.GC$fit)*(A-glm.AG$fit)*A);
      
      #variance calculation
      mm = mm.linear.withC(c(glm.GC$coef,glm.AG$coef,beta.est),
                           glm.GC=glm.GC,glm.AG=glm.AG,A=A,G=G,Y=Y,C=C,N=N)
      B = (1/N)*numDeriv::jacobian(B.deriv,c(glm.GC$coef,glm.AG$coef,beta.est),
                                   glm.GC=glm.GC,glm.AG=glm.AG,A=A,G=G,Y=Y,C=C,N=N)
      
      var = solve(B)%*%mm%*%t(solve(B));
      beta.var = (1/N)*var[length(glm.GC$coef)+length(glm.AG$coef)+1,length(glm.GC$coef)+length(glm.AG$coef)+1]
    } else {
      #single-IV estimator with no covariates
      beta.est=mean((G-mean(G))* (A-glm.AG$fit) * Y)/mean((G-mean(G))* (A-glm.AG$fit) * A);
      
      #variance calculation
      if (A.binary) {
        mm= mm.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
        B = B.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      } else {
        mm= mm.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
        B = B.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      }
      beta.var=diag((1/N)*(solve(B))%*%mm%*%t((solve(B))))[nIV+length(glm.AG$coef)+1];
    }
    
  }
  ci = beta.est + c(-1,1)*stats::qnorm(1-alpha/2)*sqrt(beta.var);
  pval = 2*pnorm(-abs(beta.est/sqrt(beta.var)));
  object <- list(beta.est=unname(beta.est), beta.var = beta.var, ci=ci, pval=unname(pval))
  class(object) <- "genius"
  return(object)
}

# "Meat" part of variance allowing for covariates in single IV model (FIX)
mm.linear.withC <- function (b, glm.GC, glm.AG, A, G, Y, C, N) {
  
  design.mat.GC <- as.matrix(glm.GC$x);
  design.mat.AG <- as.matrix(glm.AG$x);
  
  dim_gc = ncol(design.mat.GC);
  dim_ag = ncol(design.mat.AG);
  P <- dim_gc + dim_ag + 1 
  
  mean_GC = plogis(design.mat.GC%*%b[1:dim_gc])
  mean_AG = design.mat.AG%*%b[(dim_gc+1):(dim_gc+dim_ag)]
  
  #normal equations
  g.e = design.mat.GC*as.vector(G-mean_GC);
  a.e = design.mat.AG*as.vector(A-mean_AG);
  
  #IV equations
  iv.e = as.vector(G-mean_GC)*as.vector(A-mean_AG)*(Y-b[P]*A);
  
  #centering the IV moment conditions
  #iv.e = iv.e-mean(iv.e)
  
  tildem =cbind(g.e, a.e, iv.e);
  
  meat = (1/N)*(t(tildem)%*%tildem)
  
  return(meat)
}

#using this function instead because I'm too lazy to take derivs by hand
B.deriv <- function (b, glm.GC, glm.AG, A, G, Y, C, N) {
  
  design.mat.GC <- as.matrix(glm.GC$x);
  design.mat.AG <- as.matrix(glm.AG$x);
  
  dim_gc = ncol(design.mat.GC);
  dim_ag = ncol(design.mat.AG);
  P <- dim_gc + dim_ag + 1 
  
  mean_GC = plogis(design.mat.GC%*%b[1:dim_gc])
  mean_AG = design.mat.AG%*%b[(dim_gc+1):(dim_gc+dim_ag)]
  
  #normal equations
  g.e = design.mat.GC*as.vector(G-mean_GC);
  a.e = design.mat.AG*as.vector(A-mean_AG);
  
  #IV equations
  iv.e = as.vector(G-mean_GC)*as.vector(A-mean_AG)*(Y-b[P]*A);
  #iv.e = iv.e-mean(iv.e)
  
  tildem =cbind(g.e, a.e, iv.e);
  
  deriv = matrix(1,1,N)%*%as.matrix(tildem)
  
  return(deriv)
}


### Three necessary functions from the GENIUS package
genius_addY <- function(Y,A,G,formula=A~G,alpha=0.05,lower=-10,upper=10) {
  if (is.data.frame(Y)) {
    Y=as.vector(data.matrix(Y));
  }
  if (is.data.frame(A)) {
    A=as.vector(data.matrix(A));
  }
  if (is.data.frame(G)) {
    G=data.matrix(G)
  }
  if (class(G) == "matrix") {
    #number of IVs
    nIV =dim(G)[2];
    #sample size
    N = dim(G)[1];
    if (nIV==1) {G=as.vector(G);}
  } else {
    nIV =1;
    N = length(G);
  }
  A.binary = all(A %in% 0:1);
  
  if (nIV>1) {
    
    if (A.binary) {
      glm.AG = stats::glm(formula(formula), family=binomial, x=T, data=cbind(A,as.data.frame(G)));
    } else {
      glm.AG = stats::lm(formula(formula), x=T, data=cbind(A,as.data.frame(G)));
    }
    #Multiple-IV GMM
    gmm.fun <- function(beta, X) {
      sweep(X,2,apply(X,2,mean),"-")*(A-glm.AG$fit)*(Y-beta*A)
    }
    
    gmm.out = gmm::gmm(gmm.fun,x=G,t0=c(lower,upper),optfct="optimize",
                       type="iterative",wmatrix="optimal",vcov="iid",centeredVcov=TRUE);
    
    beta.est=gmm.out$coef;
    
    if (A.binary) {
      mm= mm.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      M = M.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      B = B.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
    } else {
      mm= mm.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      M = M.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      B = B.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
    }
    
    beta.var= diag((1/N)*(solve(B)%*%M)%*%mm%*%t((solve(B)%*%M)))[nIV+length(glm.AG$coef)+1];
  }
  else {
    if (A.binary) {
      glm.AG = stats::glm(formula(formula), family=binomial, x=T, data=cbind(A,as.data.frame(G)));
    } else {
      glm.AG = stats::lm(formula(formula), x=T, data=cbind(A,as.data.frame(G)));
    }
    #single-IV estimator
    beta.est=mean((G-mean(G))* (A-glm.AG$fit) * Y)/mean((G-mean(G))* (A-glm.AG$fit) * A);
    
    if (A.binary) {
      mm= mm.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      B = B.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
    } else {
      mm= mm.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
      B = B.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
    }
    beta.var=diag((1/N)*(solve(B))%*%mm%*%t((solve(B))))[nIV+length(glm.AG$coef)+1];
  }
  ci = beta.est + c(-1,1)*stats::qnorm(1-alpha/2)*sqrt(beta.var);
  pval = 2*pnorm(-abs(beta.est/sqrt(beta.var)));
  object <- list(beta.est=unname(beta.est), beta.var = beta.var, ci=ci, pval=unname(pval))
  class(object) <- "genius"
  return(object)
}

B.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
  dim_a = length(b)-1;
  B = matrix(0,(dim_a+n_iv+1),(dim_a+n_iv+1));
  B[1:n_iv,1:n_iv] = -diag(n_iv);
  B[(n_iv+1):(dim_a+n_iv),(n_iv+1):(dim_a+n_iv)] = -(1/n_sample)*(t(design.mat)%*%(design.mat));
  if (n_iv>1) {
    A_prime = dudbeta.linear(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
    Omega_inv = W.opt.linear(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
    dudmu = -diag(n_iv)*mean(as.vector(exposure-(design.mat%*%b[1:dim_a]))*(outcome-b[dim_a+1]*exposure));
    dudpsi = -(1/n_sample)*t(sweep(iv,2,apply(iv,2,mean),"-")*(outcome-b[dim_a+1]*exposure))%*%design.mat;
    B[(dim_a+n_iv+1),] = (A_prime%*%Omega_inv)%*%cbind(dudmu, dudpsi,t(A_prime));
  }
  else {
    dudmu = -1*mean(as.vector(exposure-(design.mat%*%b[1:dim_a]))*(outcome-b[dim_a+1]*exposure));
    # dudpsi= -1*apply( (iv-mean(iv))*(outcome-b[dim_a+1]*exposure)*(design.mat) ,2, mean);##RS:20180429change
    dudpsi= -1*apply( c((iv-mean(iv))*(outcome-b[dim_a+1]*exposure))*(design.mat) ,2, mean);
    dudb  = -mean((iv-mean(iv))*as.vector(exposure-(design.mat%*%b[1:dim_a]))*exposure);
    B[(dim_a+n_iv+1),] = c(dudmu, dudpsi, dudb);
  }
  B
}


mm.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
  
  dim_a = length(b)-1;
  if (n_iv>1) {
    #G-/mu
    mu.e = sweep(iv,2,apply(iv,2,mean),"-");
    #normal equations
    nr.e = design.mat*as.vector(exposure-(design.mat%*%b[1:dim_a]));
    #IV equations
    iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-(design.mat%*%b[1:dim_a]))*(outcome-b[dim_a+1]*exposure);
    #centering the IV moment conditions
    iv.e = sweep(iv.e,2, apply(iv.e,2,mean),"-");
    tildem =cbind(mu.e, nr.e, iv.e);
  }
  else {
    #G-/mu
    mu.e = iv-mean(iv);
    #normal equations
    nr.e = design.mat*as.vector(exposure-(design.mat%*%b[1:dim_a]));
    #IV equations
    iv.e = (iv-mean(iv))*as.vector(exposure-(design.mat%*%b[1:dim_a]))*(outcome-b[dim_a+1]*exposure);
    #centering the IV moment conditions
    iv.e = iv.e-mean(iv.e);
    tildem =cbind(mu.e, nr.e, iv.e);
  }
  (1/n_sample)*(t(tildem)%*%tildem);
}
