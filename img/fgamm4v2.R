# Mathew McLean
# June 13, 2013
# test FLM vs FGAM using RLRsim
# in this version attempt to reduce test to testing one component
#     by setting sigma2=alpha*sigma3 where alpha= hat(sigma2)/hat(sigma3)
# to do: implement oracle tests with sigma_FLM and alpha fixed at true values

require(mgcv)
require(lme4)
#require(gamm4)
require(fda)
require(RLRsim)
#require(caret)

testfunRatio <- function(y,X,tvals=seq(0,1,l=ncol(X)),family=gaussian(),
                   splinepars=list(k=c(8,8),m=list(c(2,2),c(2,2),extendXrange=.01)),
                   REML=TRUE,oracle=0,sig2FLM=NULL,ratio=NULL,twoVC=FALSE ){
  
  N <- length(y)
  numT <- length(tvals)
  Kx <- splinepars$k[1]
  Kt <- splinepars$k[2]
  dx <- 2
  dt <- 2
  extendXrange <- ifelse(is.null(splinepars$extendXrange),.01,
                         splinepars$extendXrange)
  
  DDx <- crossprod({
    DD <- diag(Kx)
    if(dx>0) for(i in 1:dx) DD <- diff(DD)
    DD
  }) 
  DDt <- crossprod({
    DD <- diag(Kt)
    if(dt>0) for(i in 1:dt) DD <- diff(DD)
    DD
  })
  
  
  # eigen decomp of marginal penalties:
  eDDx <- eigen(DDx)
  nullx <- (Kx-dx+1):Kx
  if(dx>0) eDDx$values[nullx] <- 0
  
  eDDt <- eigen(DDt)
  nullt <- (Kt-dt+1):Kt
  if(dt>0) eDDt$values[nullt] <- 0
  
  bbt <- create.bspline.basis(range(tvals), nbasis=Kt)
  bbtevals <- eval.basis(tvals, bbt)
  
  
  Xrange <-range(X, na.rm=TRUE)
  extXrange <- Xrange + c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx <- create.bspline.basis(extXrange, nbasis=Kx)
  L <- rep(1/numT, numT)
  
  UDinvmat_t <- eDDt$vectors[,1:(Kt-2)]%*%diag(1/sqrt(eDDt$values[1:(Kt-2)]))
  Zt <- bbtevals%*%UDinvmat_t
  UDinvmat_x <- eDDx$vectors[,1:(Kx-2)]%*%diag(1/sqrt(eDDx$values[1:(Kx-2)]))
  
  ###########################################################################
  # 1)
  # fit full model to obtain estimate of alpha1=sigma2/sigma3 and alpha2=sigma3/sigma2
  if(!twoVC){
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                    Zt*xvals, Zx,Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) )  
    
    nupp <- dx*dt-1
    upind <- 1:nupp
    dims <- if(FALSE) c(nupp,Kt-2,Kx-2,Kx-2,(Kx-2)*(Kt-2)) else c(nupp,Kt-2,2*(Kx-2),(Kx-2)*(Kt-2))
    dnames <- if(FALSE) c('para','x.Zt','Zx','Zx.t','Zx.Zt') else c('para','x.Zt','onet.Zx','Zx.Zt')
    npars <- sum(dims)
    npp <- npars-nupp
    parinds <- 1:npars
    names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
    # browser()
    colnames(Design) <- names(parinds)
    
    dat <- list()
    dat$y <- y
    dat$X <- Design[,1:nupp]
    # browser()
    form <- 'y~0+X'
    for(i in 2:length(dims)){
      dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
      form <- paste(form,'+ (1|',dnames[i],')')
    }
    
    form <-as.formula(form)
    tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
    # browser()
  
    names(tempfit$fr$fixef) <- c('1','x','x*t')
  
    tn <- names(tempfit$FL$fl)
    tn <- tn[attr(tempfit$FL$fl,'assign')]
    ind <- 1:length(tn)
    sn <- dnames[2:length(dnames)]
    snames <- if(FALSE) c('x*f(t)','f(x)','t*f(x)','f(x,t)') else c('x*f(t)','f(x)+t*f(x)','f(x,t)')
    for(i in 1:(length(dnames)-1)){
      k <- ind[sn[i]==tn]
      tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
        as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
      attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
    }
    
    fullfit <- do.call(lme4:::lmer_finalize,tempfit)
   # browser()
    temp <- lme4:::VarCorr(fullfit)
    if(temp[[1]]==0 | temp[[2]]==0)
      stop('Variance component estimated to be zero')
    alpha1 <- attr(temp[[1]],'stddev')/attr(temp[[2]],'stddev')
    alpha2 <- 1/alpha1
  }
  
#  can't use RLRTSim directly because need log-likelihood values from lme4 fits under null and alternative  
#   qrX <- qr(Design[,1:nupp])
#   Z <- cbind(Design[,(sum(dims[1:2])+1):sum(dims[1:3])],
#              Design[,(sum(dims[1:3])+1):ncol(Design)])
#   test1 <- RLRTSim(Design[,1:nupp],Z=Z,qrX=qrX,
#                    sqrt.Sigma=diag(ncol(Z)))
  
  
  ##########################################################################
  # 2)
  # fit no nuisance parm model with only one var comp. m to exactRLRT or exactLRT for test 1
  # alpha = sigma2/sigma3
  if(!twoVC){
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                    alpha1*Zx,alpha1*Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) )  
  }else{
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                    Zx,Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) )  
  }
  
  nupp <- dx*dt-1
  
  upind <- 1:nupp
  dims <- c(nupp,2*(Kx-2)+(Kx-2)*(Kt-2))
  dnames <- c('para','onet.Zx')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )

  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
    names(tempfit$fr$fixef) <- c('1','x','x*t')
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2]
  snames <- 'f_1(x)+t*f_2(x)'
  for(i in 1:1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  #browser()
  m1 <- do.call(lme4:::lmer_finalize,tempfit)
 
#   ##################################################################################################
#   # 3)
#   # first fit no nuisance parm model with only one var comp. m arg to exactRLRT or exactLRT for test 2
#   
#   Design <- t( apply(X,1,function(xvals){
#     Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
#     Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
#                   Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
#     crossprod(L,Bxi)
#   }) )  
#   
#   nupp <- dx*dt-1
#   
#   upind <- 1:nupp
#   dims <- c(nupp,(Kx-2)*(Kt-2))
#   dnames <- c('para','Zx.Zt')
#   npars <- sum(dims)
#   npp <- npars-nupp
#   parinds <- 1:npars
#   names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
#   
#   colnames(Design) <- names(parinds)
#   
#   dat <- list()
#   dat$y <- y
#   dat$X <- Design[,1:nupp]
#   # browser()
#   form <- 'y~0+X'
#   for(i in 2:length(dims)){
#     dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
#     form <- paste(form,'+ (1|',dnames[i],')')
#   }
#   
#   form <-as.formula(form)
#   tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
#   # browser()
#   names(tempfit$fr$fixef) <- c('1','x','x*t')
#   
#   tn <- names(tempfit$FL$fl)
#   tn <- tn[attr(tempfit$FL$fl,'assign')]
#   ind <- 1:length(tn)
#   sn <- dnames[2]
#   snames <- 'f(x,t)'
#   for(i in 1){
#     k <- ind[sn[i]==tn]
#     tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
#       as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
#     attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
#   }
#   
#   m2 <- do.call(lme4:::lmer_finalize,tempfit)
  
  ###########################################################################################
  # 3)
  # fit model under alternative with alpha1 (mA arg for test one)
  if(!twoVC){
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                    Zt*xvals,alpha1*Zx, alpha1*Zx*tvals,Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) ) 
  }else{
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                    Zt*xvals,Zx, Zx*tvals,Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) ) 
  }
  
  nupp <- dx*dt-1
  
  upind <- 1:nupp
  dims <- c(nupp,Kt-2,2*(Kx-2)+(Kx-2)*(Kt-2))
  dnames <- c('para','x.Zt','onet.Zx')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )

  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  names(tempfit$fr$fixef) <- c('1','x','x*t')
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:3]
  snames <- c('x*f(t)','f_1(x)+t*f_2(x)')
  for(i in 1:2){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  mA1 <- do.call(lme4:::lmer_finalize,tempfit)
  
#   ###########################################################################################
#   # 4)
#   # fit model under alternative with [1 t].Zx component zero (mA arg for test two)
#   
#   Design <- t( apply(X,1,function(xvals){
#     Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
#     Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
#                   Zt*xvals,Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
#     crossprod(L,Bxi)
#   }) ) 
#   
#   nupp <- dx*dt-1
#   
#   upind <- 1:nupp
#   dims <- c(nupp,Kt-2,(Kx-2)*(Kt-2))
#   dnames <- c('para','x.Zt','Zx.Zt')
#   npars <- sum(dims)
#   npp <- npars-nupp
#   parinds <- 1:npars
#   names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
#   
#   colnames(Design) <- names(parinds)
#   
#   dat <- list()
#   dat$y <- y
#   dat$X <- Design[,1:nupp]
#   # browser()
#   form <- 'y~0+X'
#   for(i in 2:length(dims)){
#     dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
#     form <- paste(form,'+ (1|',dnames[i],')')
#   }
#   
#   
#   
#   form <-as.formula(form)
#   tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
#   # browser()
#   names(tempfit$fr$fixef) <- c('1','x','x*t')
#   
#   tn <- names(tempfit$FL$fl)
#   tn <- tn[attr(tempfit$FL$fl,'assign')]
#   ind <- 1:length(tn)
#   sn <- dnames[2:3]
#   snames <- c('x*f(t)','f(x,t)')
#   for(i in 1:2){
#     k <- ind[sn[i]==tn]
#     tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
#       as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
#     attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
#   }
#   
#   mA2 <- do.call(lme4:::lmer_finalize,tempfit)
  
  ###########################################################################################
  # 4)
  # fit model under null (FLM) with one variance component (m0 arg to exactRLRT) same for test 1 and 2
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zt*xvals) #=Z_p
    crossprod(L,Bxi)
  }) ) 
  
  nupp <- dx*dt-1
  
  upind <- 1:nupp
  dims <- c(nupp,Kt-2)
  dnames <- c('para','x.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  names(tempfit$fr$fixef) <- c('1','x','x*t')
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:2]
  snames <- 'x*f(t)'
  for(i in 1:1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  m0 <- do.call(lme4:::lmer_finalize,tempfit)
  
  #browser()
  #test1 <- exactLRT(m=m1,mA=mA1,m0=m0)
  test1r <- exactRLRT(m1,mA1,m0)
  #test2 <- exactLRT(m1,mA1,m0)
  #test2r <- exactRLRT(m2,mA2,m0)
  
  pval1r <- test1r$p
#  pval2r <- test2r$p
  
  reject90 <- pval1r < .05
  reject95 <- pval1r < .025 
  reject99 <- pval1r < .005 
  
  ret <- c(
    pval_np=pval1r,reject90=reject90,reject95=reject95,reject99=reject99
    )
  
  return(ret)
}

testfun <- function(y,X,tvals=seq(0,1,l=ncol(X)), fixed.effects, family=gaussian(),
                    splinepars=list(k=c(8,8), m=list(c(2,2),c(2,2)), extendXrange=.01),
                    REML=TRUE, oracle=0, formratio=FALSE, sig2FLM=NULL, ratio=NULL,
                    psplines=TRUE){
  
  N <- length(y)
  numT <- length(tvals)
  if (missing(fixed.effects)){
    num.fe <- 0
  }else if (is.null(dim(fixed.effects))){
    num.fe <- 1
  }else{
    num.fe <- ncol(fixed.effects)
  }
  Kx <- splinepars$k[1]
  Kt <- splinepars$k[2]
  dx <- 2
  dt <- 2
  extendXrange <- ifelse(is.null(splinepars$extendXrange),.01,
                         splinepars$extendXrange)

  bbt <- create.bspline.basis(range(tvals), nbasis=Kt)
  bbtevals <- eval.basis(tvals, bbt)
  
  
  Xrange <-range(X, na.rm=TRUE)
  extXrange <- Xrange + c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx <- create.bspline.basis(extXrange, nbasis=Kx)
  L <- rep(1/numT, numT)
  
  if(psplines){
    DDx <- crossprod ({
      DD <- diag(Kx)
      if (dx>0) 
        for (i in 1:dx) 
          DD <- diff(DD)
      DD
    }) 
    DDt <- crossprod({
      DD <- diag(Kt)
      if (dt>0) 
        for (i in 1:dt) 
          DD <- diff(DD)
      DD
    })
  }else{
    DDx <- getbasispenalty(bbx, int2Lfd(2))
    DDt <- getbasispenalty(bbt, int2Lfd(2))
  }
  
  # eigen decomp of marginal penalties:
  eDDx <- eigen(DDx)
  nullx <- (Kx-dx+1):Kx
  if(dx>0) eDDx$values[nullx] <- 0
  
  eDDt <- eigen(DDt)
  nullt <- (Kt-dt+1):Kt
  if(dt>0) eDDt$values[nullt] <- 0
  
  UDinvmat_t <- eDDt$vectors[,1:(Kt-2)]%*%diag(1/sqrt(eDDt$values[1:(Kt-2)]))
  Zt <- bbtevals%*%UDinvmat_t
  UDinvmat_x <- eDDx$vectors[,1:(Kx-2)]%*%diag(1/sqrt(eDDx$values[1:(Kx-2)]))
  
  ##########################################################################
  # 1)
  # first fit no nuisance parm model with only one var comp. m to exactRLRT or exactLRT for test 1
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zx, Zx*tvals) #=Z_p
    crossprod(L,Bxi)
  }) )  
  
  nupp <- dx*dt-1
  if (num.fe)
    Design <- cbind(fixed.effects,Design)
  
  nupp <- nupp + num.fe
  
  upind <- 1:nupp
  dims <- c(nupp,2*(Kx-2))
  dnames <- c('para','onet.Zx')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  if (num.fe){
    names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
  }else{
    names(tempfit$fr$fixef) <- c('1','x','x*t')
  }
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2]
  snames <- 'f_1(x)+t*f_2(x)'
  for(i in 1:1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  #browser()
  m1 <- do.call(lme4:::lmer_finalize,tempfit)
  
  ##################################################################################################
  # 2)
  # first fit no nuisance parm model with only one var comp. m arg to exactRLRT or exactLRT for test 2
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
    crossprod(L,Bxi)
  }) )  
  
  nupp <- dx*dt-1
  
  upind <- 1:nupp
  dims <- c(nupp,(Kx-2)*(Kt-2))
  dnames <- c('para','Zx.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  names(tempfit$fr$fixef) <- c('1','x','x*t')
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2]
  snames <- 'f(x,t)'
  for(i in 1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  m2 <- do.call(lme4:::lmer_finalize,tempfit)
  
  ###########################################################################################
  # 3)
  # fit model under alternative with Zx.Zt component zero (mA arg for test one)
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zt*xvals,Zx, Zx*tvals) #=Z_p
    crossprod(L,Bxi)
  }) ) 
  
  nupp <- dx*dt-1
  
  upind <- 1:nupp
  dims <- c(nupp,Kt-2,2*(Kx-2))
  dnames <- c('para','x.Zt','onet.Zx')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  names(tempfit$fr$fixef) <- c('1','x','x*t')
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:3]
  snames <- c('x*f(t)','f_1(x)+t*f_2(x)')
  for(i in 1:2){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  mA1 <- do.call(lme4:::lmer_finalize,tempfit)
  
  ###########################################################################################
  # 4)
  # fit model under alternative with [1 t].Zx component zero (mA arg for test two)
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zt*xvals,Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
    crossprod(L,Bxi)
  }) ) 
  
  nupp <- dx*dt-1
  
  upind <- 1:nupp
  dims <- c(nupp,Kt-2,(Kx-2)*(Kt-2))
  dnames <- c('para','x.Zt','Zx.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  names(tempfit$fr$fixef) <- c('1','x','x*t')
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:3]
  snames <- c('x*f(t)','f(x,t)')
  for(i in 1:2){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  mA2 <- do.call(lme4:::lmer_finalize,tempfit)
  
  ###########################################################################################
  # 5)
  # fit model under null (FLM) with one variance component (m0 arg to exactRLRT) same for test 1 and 2
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zt*xvals) #=Z_p
    crossprod(L,Bxi)
  }) ) 
  
  nupp <- dx*dt-1
  
  upind <- 1:nupp
  dims <- c(nupp,Kt-2)
  dnames <- c('para','x.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  names(tempfit$fr$fixef) <- c('1','x','x*t')
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:2]
  snames <- 'x*f(t)'
  for(i in 1:1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  m0 <- do.call(lme4:::lmer_finalize,tempfit)
 # browser()
  var1 <- lme4:::VarCorr(mA1)[[1]]  #onet.Zx
  var1.0 <- lme4:::VarCorr(m1)[[1]]
  var2 <- lme4:::VarCorr(mA2)[[1]]  #Zt.Zx
  var2.0 <- lme4:::VarCorr(m2)[[1]]
  var0 <- lme4:::VarCorr(m0)[[1]]
  if(var0==0)
    return(rep(NA,5))
  if((var1==0 | var1.0==0) & (var2==0 | var2.0==0)){
    return(c(NA,NA,0,0,0))
  }else if(var1>0 & (var2==0 | var2.0==0)){
    test1r <- exactRLRT(m1,mA1,m0)
    pval1r <- test1r$p
    pval2r <- NA
    reject90 <- pval1r < .05 
    reject95 <- pval1r < .025
    reject99 <- pval1r < .005 
  }else if((var1==0 | var1.0==0) & var2 >0){
    test2r <- exactRLRT(m2,mA2,m0)
    pval1r <- NA
    pval2r <- test2r$p
    reject90 <- pval2r < .05
    reject95 <- pval2r < .025
    reject99 <- pval2r < .005
  }else{
    test1r <- exactRLRT(m1,mA1,m0)
    test2r <- exactRLRT(m2,mA2,m0)  
    pval1r <- test1r$p
    pval2r <- test2r$p    
    reject90 <- pval1r < .05 & pval2r < .05
    reject95 <- pval1r < .025 & pval2r < .025
    reject99 <- pval1r < .005 & pval2r < .005
  }
  
  ret <- c(
    pval_np=pval1r,pval_pp=pval2r,reject90=reject90,reject95=reject95,reject99=reject99
  )
  
  return(ret)
}

FGAMM4 <- function(y, X, tvals=seq(0,1,l=ncol(X)),fixed.effects, family=gaussian(),
                   splinepars=list (k=c(8,8), m=list (c(2,2),c(2,2)), extendXrange=.01),
                   checkRankDeficiency=0, REML=TRUE, nullBasis=0, full=TRUE, 
                   psplines=TRUE){
  
  N <- length(y)
  numT <- length(tvals)
  if (missing(fixed.effects)){
    num.fe <- 0
  }else if (is.null(dim(fixed.effects))){
    num.fe <- 1
  }else{
    num.fe <- ncol(fixed.effects)
  }
  Kx <- splinepars$k[1]
  Kt <- splinepars$k[2]
  dx <- 2
  dt <- 2
  extendXrange <- ifelse (is.null(splinepars$extendXrange), .00,
                         splinepars$extendXrange)
  
  bbt <- create.bspline.basis(range(tvals), nbasis=Kt)
  bbtevals <- eval.basis(tvals, bbt)
  
  
  Xrange <-range(X, na.rm=TRUE)
  extXrange <- Xrange + c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx <- create.bspline.basis(extXrange, nbasis=Kx)
  L <- rep(1/numT, numT)
  
  if(psplines){
    DDx <- crossprod ({
      DD <- diag(Kx)
      if (dx>0) 
        for (i in 1:dx) 
          DD <- diff(DD)
      DD
    }) 
    DDt <- crossprod({
      DD <- diag(Kt)
      if (dt>0) 
        for (i in 1:dt) 
          DD <- diff(DD)
      DD
    })
  }else{
    DDx <- getbasispenalty(bbx, int2Lfd(2))
    DDt <- getbasispenalty(bbt, int2Lfd(2))
  }
  
  
  # eigen decomp of marginal penalties:
  eDDx <- eigen(DDx)
  nullx <- (Kx-dx+1):Kx
  if(dx>0) 
    eDDx$values[nullx] <- 0
  
  eDDt <- eigen(DDt)
  nullt <- (Kt-dt+1):Kt
  if(dt>0) 
    eDDt$values[nullt] <- 0
  
  if(TRUE){
    UDinvmat_t <- eDDt$vectors[, 1:(Kt-2)]%*%diag(1/sqrt(eDDt$values[1:(Kt-2)]))
    Zt <- bbtevals%*%UDinvmat_t
    UDinvmat_x <- eDDx$vectors[, 1:(Kx-2)]%*%diag(1/sqrt(eDDx$values[1:(Kx-2)]))
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      Bxi<-	cbind( 1+numeric(numT), xvals, xvals*tvals,  # =Z_0
                   Zt*xvals, Zx,Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) )  
    
    nupp <- dx*dt-1
  }
  if (num.fe)
    Design <- cbind(fixed.effects,Design)
  
  nupp <- nupp + num.fe
  upind <- 1:nupp
  dims <- if (full) c(nupp, Kt-2, Kx-2, Kx-2, (Kx-2)*(Kt-2)) else c(nupp,Kt-2,2*(Kx-2),(Kx-2)*(Kt-2))
  dnames <- if (full) c('para','x.Zt','Zx','Zx.t','Zx.Zt') else c('para','x.Zt','onet.Zx','Zx.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  # browser()
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  if (TRUE){
    if (num.fe){
      names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
    }else{
      names(tempfit$fr$fixef) <- c('1','x','x*t')
    }
  }
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:length(dnames)]
  snames <- if(full) c('x*f(t)','f(x)','t*f(x)','f(x,t)') else c('x*f(t)','f(x)+t*f(x)','f(x,t)')
  for (i in 1:(length(dnames)-1)){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  ret <- list()
  ret$lme4fit <- do.call(lme4:::lmer_finalize,tempfit)
  ret$bbx <- bbx
  ret$bbt <- bbt
  ret$UDinvmat_x <- UDinvmat_x
  ret$UDinvmat_t <- UDinvmat_t
  ret$parinds <- parinds
  ret$dnames <- dnames
  ret$full <- full
  
  vc <- lme4:::VarCorr(ret$lme4fit)
  ret$sigma2 <- as.numeric(attr(vc,'sc'))^2
  sp <- numeric(4)
  for (i in 1:(length(dnames)-1)){
    if (as.numeric(vc[[i]]) >0 ){
      sp[i] <- ret$sigma2/as.numeric(vc[[i]])
    }else{
      sp[i] <- 1e10
    }
  }
  names(sp) <- names(vc)
  sp <- sp[dnames[-1]]
  names(sp) <- snames
  ret$sp <- sp
  blups <- unlist(as(lme4::ranef(ret$lme4fit), 'vector'))
  ind <- unlist(sapply(snames,function(x){grep(x,
                                               substr(names(blups), 1, nchar(x)), fixed=T)}))
  blups <- blups[ind]
  names(blups) <- unlist( mapply(FUN=function(y,x)rep(y,e=x), y=snames, x=dims[-1]) )  
  ret$snames <- snames  
  ret$blups <- blups
  ret$coefficients <- c(lme4::fixef(ret$lme4fit),blups)
  ret$residuals <- residuals(ret$lme4fit)
  ret$fitted.values <- dat$y - ret$residuals
  ret$nullBasis <- nullBasis
  ret$psplines <- psplines
  
  return(ret)
}

flmm4 <- function(y,X,tvals=seq(0,1,l=ncol(X)),family=gaussian(),
                   splinepars=list(k=c(8,8),m=list(c(2,2),c(2,2)),extendXrange=.01),
                   checkRankDeficiency=0,REML=TRUE,nullBasis=0 ){
  
  N <- length(y)
  numT <- length(tvals)
  Kx <- splinepars$k[1]
  Kt <- splinepars$k[2]
  dx <- 2
  dt <- 2
  extendXrange <- ifelse(is.null(splinepars$extendXrange),.01,
                         splinepars$extendXrange)
  
  DDx <- crossprod({
    DD <- diag(Kx)
    if(dx>0) for(i in 1:dx) DD <- diff(DD)
    DD
  }) 
  DDt <- crossprod({
    DD <- diag(Kt)
    if(dt>0) for(i in 1:dt) DD <- diff(DD)
    DD
  })
  
  
  # eigen decomp of marginal penalties:
  eDDx <- eigen(DDx)
  nullx <- (Kx-dx+1):Kx
  if(dx>0) eDDx$values[nullx] <- 0
  
  eDDt <- eigen(DDt)
  nullt <- (Kt-dt+1):Kt
  if(dt>0) eDDt$values[nullt] <- 0
  
  bbt <- create.bspline.basis(range(tvals), nbasis=Kt)
  bbtevals <- eval.basis(tvals, bbt)
  
  
  Xrange <-range(X, na.rm=TRUE)
  extXrange <- Xrange + c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx <- create.bspline.basis(extXrange, nbasis=Kx)
  L <- rep(1/numT, numT)
  
  if(!nullBasis){
    UDinvmat_t <- eDDt$vectors[,1:(Kt-2)]%*%diag(1/sqrt(eDDt$values[1:(Kt-2)]))
    Zt <- bbtevals%*%UDinvmat_t
    UDinvmat_x <- eDDx$vectors[,1:(Kx-2)]%*%diag(1/sqrt(eDDx$values[1:(Kx-2)]))
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    #  Bxi<-	cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
    #               Zx, Zt*xvals, Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, Zt*xvals) #=Z_p
      crossprod(L,Bxi)
    }) )  
    
    nupp <- dx*dt-1
  }else if(nullBasis==2){
    stop('null basis==2 not implemented')
    UDinvmat_t <- cbind(eDDt$vectors[,1:(Kt-2)]/eDDt$values[1:(Kt-2)],eDDt$vectors[,(Kt-1):Kt])
    Zt <- bbtevals%*%UDinvmat_t
    UDinvmat_x <- cbind(eDDx$vectors[,1:(Kx-2)]/eDDx$values[1:(Kx-2)],eDDx$vectors[,(Kx-1):Kx])
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      Bxi<-	cbind( Zx[,(Kx-1):Kx],Zt[,(Kt-1):Kt], #=Z_0
                   Zx[,1:(Kt-2)], Zt[,1:(Kt-2)]*xvals, Zx[,1:(Kt-2)]*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) )  
    
    nupp <- dx*dt    
  }else{
    stop('null basis!=0 not implemented')
    UDinvmat_t <- cbind(eDDt$vectors[,1:(Kt-2)]/eDDt$values[1:(Kt-2)],eDDt$vectors[,(Kt-1):Kt])
    Zt <- bbtevals%*%UDinvmat_t
    UDinvmat_x <- cbind(eDDx$vectors[,1:(Kx-2)]/eDDx$values[1:(Kx-2)],eDDx$vectors[,(Kx-1):Kx])
    nB2 <- Zt[,(Kt-1):Kt]
    nB2 <- (nB2%*%eigen(crossprod(scale(nB2,scale=F)))$vec)[,2:1]
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      nB1 <- Zx[,(Kx-1):Kx]
      nB1 <- (nB1%*%eigen(crossprod(scale(nB1,scale=F)))$vec)[,2:1]
      nB <- nB1[,c(1,2,2)]
      nB[,3] <- nB[,3]*nB2[,2]
      Bxi <- cbind( nB, 
                    Zx[,1:(Kx-2)]*nB2[,1],Zt[,1:(Kt-2)]*nB1[,2],Zx[,1:(Kx-2)]*nB2[,2], 
                    Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) )  
    nupp <- 3
  }
  upind <- 1:nupp
  dims <- c(nupp,Kt-2) #c(nupp,Kx-2,Kt-2,Kx-2,(Kx-2)*(Kt-2))
  dnames <- c('para','x.Zt') #c('para','Zx','x.Zt','Zx.t','Zx.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  #browser()
 # print(findLinearCombos(Design)$remove)
  
  # find numerical rank deficiency and remove these columns from fitting
  if(checkRankDeficiency){
    stop('checking rank deficiency not implemented')
    removeind <- NULL
    if(checkRankDeficiency==1){
      for( i in 1:length(dims) ){
        curind <- which(names(parinds)==dnames[i])
        eigvals <- eigen(crossprod(Design[,curind]) ,TRUE,TRUE,FALSE)$val
        temp <- curind[eigvals < eigvals[1]*sqrt(.Machine$double.eps)]
        dims[i] <- dims[i]-length(temp)
        removeind <- c(removeind,temp)
      }
    }else if(checkRankDeficiency==2){
      eigvals <- eigen(crossprod(Design) ,TRUE,TRUE,FALSE)$val
      removeind <- parinds[eigvals < eigvals[1]*sqrt(.Machine$double.eps)]
      mt <- table(names(removeind))
      if(length(mt)>0){
        for(i in 1:length(mt)){
          dims[dnames==names(mt)[i]] <- dims[dnames==names(mt)[i]] - mt[[i]]
        }
      }
      #browser()
    }else{
      removeind <- findLinearCombos(Design)$remove
      mt <- table(names(parinds[removeind]))
      if(length(mt)>0){
        for(i in 1:length(mt)){
          dims[dnames==names(mt)[i]] <- dims[dnames==names(mt)[i]] - mt[[i]]
        }
      }
    }
    if(length(removeind)>0){
      Design <- Design[,-removeind]
      parinds <- parinds[-removeind]
      npars <- npars-length(removeind)
      npp <- sum(dims[-1])
      nupp <- dims[1]
    }
  }
  colnames(Design) <- names(parinds)
  parinds <- 1:ncol(Design)
  names(parinds) <- colnames(Design)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  if(!nullBasis){
    names(tempfit$fr$fixef) <- c('1','x','x*t')
  }else if(nullBasis==2){
    if(checkRankDeficiency){
      if(!any(removeind<5))
        names(tempfit$fr$fixef) <- c('nullX1','nullX2','nullT1','nullT2')
    }else{
      names(tempfit$fr$fixef) <- c('nullX1','nullX2','nullT1','nullT2')
    }
  }else{
    if(checkRankDeficiency){
      if(!any(removeind<5))
        names(tempfit$fr$fixef) <- c('1','nX2','nX2.nT2')
    }else{
      names(tempfit$fr$fixef) <- c('1','nX2','nX2.nT2')
    }   
  }
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2] #dnames[2:5]
  snames <- 'x*f(t)' #c('f(x)','x*f(t)','t*f(x)','f(x,t)')
  for(i in 1:1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  ret <- list()
  ret$lme4fit <- do.call(lme4:::lmer_finalize,tempfit)
  ret$bbx <- bbx
  ret$bbt <- bbt
  ret$UDinvmat_x <- UDinvmat_x
  ret$UDinvmat_t <- UDinvmat_t
  ret$parinds <- parinds
  ret$dnames <- dnames
  if(checkRankDeficiency)  ret$removeind <- removeind
  vc <- lme4:::VarCorr(ret$lme4fit)
  ret$sigma2 <- as.numeric(attr(vc,'sc'))^2
  sp <- numeric(1)
  for(i in 1:1){
    if( as.numeric(vc[[i]]) >0 ){
      sp[i] <- ret$sigma2/as.numeric(vc[[i]])
    }else{
      sp[i] <- 1e10
    }
  }
  names(sp) <- names(vc)
  sp <- sp[dnames[-1]]
  names(sp) <- snames
  ret$sp <- sp
  blups <- unlist(as(lme4::ranef(ret$lme4fit),'vector'))
  ind <- unlist(sapply(snames,function(x){grep(x,
                                               substr(names(blups),1,nchar(x)),fixed=T)}))
  blups <- blups[ind]
  names(blups) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=snames,x=dims[-1]) )  
  ret$snames <- snames  
  ret$blups <- blups
  ret$coefficients <- c(lme4::fixef(ret$lme4fit),blups)
  ret$residuals <- residuals(ret$lme4fit)
  ret$fitted.values <- dat$y - ret$residuals
  ret$nullBasis <- 0
  
  return(ret)
}

# constructBasis=function(xvals,tvals=seq(0,1,l=length(xvals)),bbx,bbt,
# 				UDinvmat_t,UDinvmat_x){
#   Kx <- bbx$nbasis
#   Kt <- bbt$nbasis
#   numT <- length(xvals)
#   stopifnot(length(xvals)==length(tvals))
#   bbtevals <- eval.basis(tvals, bbt)
#   Zt <- bbtevals%*%UDinvmat_t
# 
#   Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
#   Design <- cbind( 1+numeric(numT), xvals, xvals*tvals, Zx, Zt*xvals,
#          Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)])
#   
#   return(Design)
# }

plotF <- function(fit,xvals,tvals=seq(0,1,l=length(xvals)),trueF=NULL,sx=NULL,
	st=NULL,difference=FALSE,components=FALSE,Xobs=NULL,tobs=NULL,insideOnly=FALSE,
                  pretty=TRUE, contour=FALSE, use.lattice=TRUE, ...){
  
  if(insideOnly) stopifnot(!is.null(Xobs) & !is.null(tobs))
  gridvals=expand.grid(tvals,xvals)
  btvals <-  eval.basis(gridvals$Var1,fit$bbt)
  bxvals <- eval.basis(gridvals$Var2,fit$bbx)
  Zx <- bxvals%*%fit$UDinvmat_x
  Zt <- btvals%*%fit$UDinvmat_t
  ncx <- ncol(Zx)
  nct <- ncol(Zt)
  #browser()
  if(insideOnly){
    Xrng <- cbind(rep(tobs, e=nrow(Xobs)), as.vector(Xobs))[chull(cbind(rep(tobs, e=nrow(Xobs)), as.vector(Xobs))),]
    inside <- insidePoly(cbind(gridvals$Var1,gridvals$Var2),Xrng)
  }
 
  if(!is.null(Xobs) & !is.null(tobs)){
    dat <- data.frame(x=as.vector(Xobs),y=rep(tobs,nrow(Xobs)))
  }else{
    dat <- data.frame(x=numeric(0),y=numeric(0),z=numeric(0))
  }
  if(!fit$nullBasis){
    basisevals <- cbind(1+numeric(nrow(btvals)),gridvals$Var2,gridvals$Var2*gridvals$Var1,
                      Zt*gridvals$Var2, Zx, Zx*gridvals$Var1, 
                      Zx[,rep(1:ncx,e=nct)]*Zt[,rep(1:nct,t=ncx)])
  }else if(fit$nullBasis==2){
     basisevals <- cbind(Zx[,(ncx-1):ncx],Zt[,(nct-1):nct],Zx[,1:(ncx-2)],
                 Zt[,1:(nct-2)]*gridvals$Var2, Zx[,1:(ncx-2)]*gridvals$Var1, 
                 Zx[,rep(1:(ncx-2),e=nct-2)]*Zt[,rep(1:(nct-2),t=ncx-2)])
  }else{
    nB1 <- Zx[,(ncx-1):ncx]
    nB1 <- (nB1%*%eigen(crossprod(scale(nB1,scale=F)))$vec)[,2:1]
    nB2 <- Zt[,(nct-1):nct]
    nB2 <- (nB2%*%eigen(crossprod(scale(nB2,scale=F)))$vec)[,2:1]
    nB <- nB1[,c(1,2,2)]
    nB[,3] <- nB[,3]*nB2[,2]
    basisevals <- cbind( nB, 
                  Zx[,1:(ncx-2)]*nB2[,1],Zt[,1:(nct-2)]*nB1[,2],Zx[,1:(ncx-2)]*nB2[,2], 
                  Zx[,rep(1:(ncx-2),e=nct-2)]*Zt[,rep(1:(nct-2),t=ncx-2)])
  }
#basisevals[,3] <- -basisevals[,3]
#basisevals[,2] <- .1*basisevals[,2]
#browser()
  #if(!is.null(fit$removeind)) basisevals <- basisevals[,-fit$removeind]
  blups <- fit$coef[names(fit$parinds)!='para']
  if (num.fe <- length(grep('fe',names(fit$coef)))){  # check if any fixed effects not part of smooth
    fit$coef <- fit$coef[-grep('fe',names(fit$coef))]
    fit$parinds <- fit$parinds[-c(1:num.fe)]
  }
  fsurfvals <- basisevals%*%fit$coef
  
  if(insideOnly){
    allres <- fsurfvals
    fsurfvals[!inside] <- NA
  }
  
  if(!is.null(Xobs) & !is.null(tobs)) 
    dat$z <- rep(min(fsurfvals,na.rm=T), length(dat$x))
  a <- wireframe(fsurfvals~gridvals$Var2*gridvals$Var1, scales=list(arrows=FALSE),
              main=list(expression(paste(phantom(000), hat(F)(x,t)))),
		  xlab='x', ylab='t',zlab='', pts=dat,panel.3d.wireframe =
          function(x, y, z,xlim, ylim, zlim,
                   xlim.scaled, ylim.scaled, zlim.scaled,pts,...) {
              panel.3dwire(x = x, y = y, z = z,xlim = xlim,ylim = ylim,
			            zlim = zlim,xlim.scaled = xlim.scaled,
			            ylim.scaled = ylim.scaled,zlim.scaled = zlim.scaled,...)
              xx <- xlim.scaled[1] + diff(xlim.scaled) *
                      (pts$x - xlim[1]) / diff(xlim)
              yy <- ylim.scaled[1] + diff(ylim.scaled) *
                      (pts$y - ylim[1]) / diff(ylim)
              zz <- zlim.scaled[1] + diff(zlim.scaled) *
                      (pts$z - zlim[1]) / diff(zlim)
              panel.3dscatter(x = xx,y = yy,z = zz,xlim = xlim,ylim = ylim,
                  zlim = zlim,xlim.scaled = xlim.scaled,ylim.scaled = ylim.scaled,
			            zlim.scaled = zlim.scaled,cex=.3,...)
          }, ...)
  if (is.null(trueF)){
    if(!components){
     print(a)
    }else{
      #browser()
      if (fit$full){
        dnames <- c('para','x.Zt','Zx','Zx.t','Zx.Zt') 
        plotlabs <- c(expression(paste(phantom(000),hat(beta)[0]+hat(beta)[1]*x+hat(beta)[2]*x%.%t)),
                     expression(paste(phantom(000),x%.%hat(f)(t))),
                     expression(paste(phantom(000),hat(g)[1](x))),
                     expression(paste(phantom(000)*t%.%hat(g)[2](x))),
                     expression(paste(phantom(000),hat(h)(x,t))))
      }else{ 
        dnames <- c('para','x.Zt','onet.Zx','Zx.Zt')
#         plotlabs <- c(expression(paste(phantom(000),hat(beta)[0]+hat(beta)[1]*x+hat(beta)[2]*x%.%t)),
#                     expression(paste(phantom(000),x%.%hat(g)(t))),
#                     expression(paste(phantom(000),hat(h)[1](x)+t%.%hat(h)[2](x))),
#                     #expression(paste(phantom(000)*t%.%hat(h)[2](x))),
#                     expression(paste(phantom(000),hat(f)(x,t))))
        plotlabs <- c(expression(paste(hat(beta)[0]+hat(beta)[1]*x+hat(beta)[2]*x%.%t)),
                      expression(paste(x%.%hat(f)(t))),
                      expression(paste(hat(g)[1](x)+t%.%hat(g)[2](x))),
                      #expression(paste(phantom(000)*t%.%hat(h)[2](x))),
                      expression(paste(hat(h)(x,t))))
      }
      
      coefs <- fit$coef
      if(!pretty){
       # browser()
        plots <- vector('list',length(dnames))
        for(i in 1:length(dnames)){
          if(any(names(fit$parinds)==dnames[i])){
            surfvals <- basisevals[,names(fit$parinds)==dnames[i]]%*%
              coefs[names(fit$parinds)==dnames[i]]
          }else{
            surfvals <- numeric(length(gridvals$Var1))
          }
          if(insideOnly) surfvals[!inside] <- NA
          plots[[i]]=wireframe(surfvals~gridvals$Var2*gridvals$Var1,
                               scales=list(arrows=FALSE),xlab='x',ylab='t',zlab='',
                               main=list(plotlabs[i]),...)
        }
        print(a,split=c(1,1,3,2),more=TRUE)
        print(plots[[1]],split=c(1,2,3,2),more=TRUE)
        print(plots[[2]],split=c(2,1,3,2),more=TRUE)
        print(plots[[3]],split=c(2,2,3,2),more=TRUE)
        if(fit$full){
          print(plots[[4]],split=c(3,1,3,2),more=TRUE)
          print(plots[[5]],split=c(3,2,3,2),more=FALSE)
        }else{
          print(plots[[4]],split=c(3,1,3,2),more=FALSE)
        }
      }else{  # pretty=TRUE
        if(insideOnly) fsurfvals[!inside] <- NA
        dat <- data.frame(t=gridvals$Var1,x=gridvals$Var2,surf=fsurfvals,
                          surf.name=rep(1,length(fsurfvals)))
        for(i in 1:length(dnames)){
          if(any(names(fit$parinds)==dnames[i])){
            surfvals <- basisevals[,names(fit$parinds)==dnames[i]]%*%
              coefs[names(fit$parinds)==dnames[i]]
            if(insideOnly) surfvals[!inside] <- NA
            tdat <- data.frame(t=gridvals$Var1,x=gridvals$Var2,surf=surfvals,
                              surf.name=rep(i+1,length(fsurfvals)))
            dat <- rbind(dat,tdat)
          }else{
            stop('oops')
            surfvals <- numeric(length(gridvals$Var1))
          }

        }  # end for loop
        snames <- c(expression(paste(phantom(000),hat(F)(x,t))), plotlabs)
        #dat$surf.names <- as.factor(dat$surf.name)
        #levels(dat$surf.names) <- snames
        #browser()
        if(use.lattice){
        #  browser()
          dots <- list(...)
          args <- c(list(x = surf~x*t | surf.name, data=dat,par.strip.text=list(cex=2),
                         scales=list(arrows=FALSE, x=list(cex=1)), xlab='x', ylab='t', zlab='',
                         strip=strip.custom(factor.levels=snames, style=1, strip.levels=c(F,T),
                                            strip.names=c(F,F)),
                         snames=snames, region=TRUE,vars=1:6), dots)
          if(contour){
            print(do.call(contourplot, args))
          }else{
            print(do.call(wireframe, args))
          }

        }else{
          snames <- as.list(snames)
           MyLabeller <- function(variable, value){
             return(bquote(.(snames[value])))
           }
          pg <- ggplot(dat, aes(x, t, z = surf, fill=surf), geom='polygon') +  #   # 
                 #geom_tile() + geom_contour(size=1.2,binwidth=1) +
            #scale_fill_continuous(low='lightpink1', high = 'purple4') +
            #scale_colour_gradient(low='thistle1', high = 'darkslateblue') +
            #scale_fill_gradientn(colours=c('thistle1','thistle2', 'plum1', 'plum4', 'lightpink1', 'lightpink3',  'purple1', 'purple4', 'slateblue4')) +
            scale_fill_gradientn(colours=c('plum', 'plum2','orchid','orchid2',
                                           'orchid3','purple1','purple2','purple3','purple4','slateblue4'))+
            scale_colour_gradientn(colours=c('thistle','thistle1', 'thistle2', 'thistle3', 'thistle4',
                                           'mediumpurple',  'mediumpurple1', 'mediumpurple2',
                                            'mediumpurple3', 'mediumpurple4')) +
            #scale_colour_brewer(palette = 'Purples') +
            geom_tile(aes(fill=surf)) + 
            stat_contour(aes(colour=..level..),  bins=10.5, size=1.1)+  # ..level..
            facet_wrap(~surf.name)+
            facet_grid(~surf.name, labeller=MyLabeller) +
           # geom_text(data = lab_dat, aes(x= x, y= y, label=labs))+
          theme(text=element_text(size=I(20)))  # legend.position = "none",
          #scale_color_brewer(...,type='seq', palette='Purples') +
        #  pg <- pg + facet_wrap(~surf.name)  # facet_wrap has no labeller. open issue with ggplot2
          print(pg)
        } #fill_contour + geom_tile(aes(fill = value)) + stat_contour(bins=15.5)
      }  # end check for pretty
#       function(which.panel,factor.levels,...){
#         lab <- 1:5 #snames
#         strip.default(which.panel,factor.levels=lab,...)
#       }
    }
  }else{
    truesurfvals <- GenSurface(xval=gridvals$Var2,tval=gridvals$Var1,type=trueF,sx=sx,st=st) 
    if(!is.null(Xobs) & !is.null(tobs)) dat$z <- rep(min(truesurfvals-fsurfvals),length(dat$x))
    if(insideOnly){
      allres <- sqrt(mean((truesurfvals-allres)^2))
      truesurfvals[!inside] <- NA 
    }
    if(difference){
      a=wireframe(truesurfvals-fsurfvals~gridvals$Var2*gridvals$Var1,scales=list(arrows=FALSE),
                main=expression(paste(F(x,t)-hat(F)(x,t))),xlab='x',ylab='t',
		   zlab='',pts=dat,panel.3d.wireframe =
          function(x, y, z,xlim, ylim, zlim,
                   xlim.scaled, ylim.scaled, zlim.scaled,pts,...) {
              panel.3dwire(x = x, y = y, z = z,xlim = xlim,ylim = ylim,
			zlim = zlim,xlim.scaled = xlim.scaled,
			ylim.scaled = ylim.scaled,zlim.scaled = zlim.scaled,...)
              xx <- xlim.scaled[1] + diff(xlim.scaled) *
                      (pts$x - xlim[1]) / diff(xlim)
              yy <- ylim.scaled[1] + diff(ylim.scaled) *
                      (pts$y - ylim[1]) / diff(ylim)
              zz <- zlim.scaled[1] + diff(zlim.scaled) *
                      (pts$z - zlim[1]) / diff(zlim)
              panel.3dscatter(x = xx,y = yy,z = zz,xlim = xlim,ylim = ylim,
                  zlim = zlim,xlim.scaled = xlim.scaled,ylim.scaled = ylim.scaled,
			zlim.scaled = zlim.scaled,cex=.3,...)
      },...)
    #  print(a)

    } # end check for difference=TRUE  if difference FALSE use 'a' formed earlier
    if(!components){
      b=wireframe(truesurfvals~gridvals$Var2*gridvals$Var1,scales=list(arrows=FALSE),
                  main=list('F(x,t)'),xlab='x',ylab='t',zlab='',...)
      print(a,split=c(1,1,1,2),more=TRUE)
      print(b,split=c(1,2,1,2),more=FALSE)          
    }else{
       dnames <- if(fit$full) c('para','x.Zt','Zx','Zx.t','Zx.Zt') else c('para','x.Zt','[1 t]*Zx','Zx.Zt')
       plots <- vector('list',length(dnames))
       plotlabs <- c(expression(paste(phantom(000),hat(beta)[0]+hat(beta)[1]*x+hat(beta)[2]*x%.%t)),
          expression(paste(phantom(000),hat(g)[1](x))),
	    expression(paste(phantom(000),x%.%hat(f)(t))),
          expression(paste(phantom(000)*t%.%hat(g)[2](x))),
	    expression(paste(phantom(000),hat(h)(x,t))))
       coefs <- fit$coef
       for(i in 1:length(dnames)){
         if(any(names(fit$parinds)==dnames[i])){
    	      surfvals <- basisevals[,names(fit$parinds)==dnames[i]]%*%
				    coefs[names(fit$parinds)==dnames[i]]
         }else{
           surfvals <- numeric(length(gridvals$Var1))
         }
         if(insideOnly) surfvals[!inside] <- NA
	       plots[[i]]=wireframe(surfvals~gridvals$Var2*gridvals$Var1,
		     scales=list(arrows=FALSE),xlab='x',ylab='t',zlab='',
            main=list(plotlabs[i]),...)
	     }

    	print(plots[[1]],split=c(1,1,2,3),more=TRUE)
    	print(plots[[2]],split=c(1,2,2,3),more=TRUE)
    	print(plots[[3]],split=c(1,3,2,3),more=TRUE)
    	print(plots[[4]],split=c(2,1,2,3),more=TRUE)
           if(fit$full)
    	        print(plots[[5]],split=c(2,2,2,3),more=TRUE)
          print(a,split=c(2,3,2,3),more=FALSE)
    }
  }  # end else for check for null trueF
  if(!is.null(trueF)){
    if(!insideOnly){
	    return(sqrt(mean( (truesurfvals-fsurfvals)^2 )))
    }else{
      return(c( all=allres,inside=sqrt(mean( (truesurfvals-fsurfvals)^2,na.rm=TRUE )) ))
    }
  }
}

CreateDataFabian <- function(n=100,ngrid=30,sigmae=1,seed=sample(.Machine$integer.max,1)){
  

  ## true function f(x(t), t)
  test1 <- function(x,t, sx=0.3, st=0.4)
  { (pi^sx*st)*(1.2*exp(-(x-0.2)^2/sx^2-(t-0.3)^2/st^2)+
                  0.8*exp(-(x-0.7)^2/sx^2-(t-0.8)^2/st^2))
  }
  
  
  #vector of timepoints (same for all obs.)
  t <- seq(0, 1, l=ngrid)
  
  #generate random trajectories x(t)
  fxt <- function(t, k=40){
    eigenfct <- function(T, k){
      Tc <- 2*pi*(T - min(T))/(max(T)-min(T))
      if (k%%2 == 1)
        result <- sin(ceiling(k/2)*Tc)
      else
        result <- cos((k/2)*Tc)
      
      return(result/sqrt((max(T)-min(T))/2))
    }
    
    phiMat <- sapply(1:k, function(i) eigenfct(t, i)) 
    eigenvals <- exp(-(seq(0, (k-1)/2, by=.5)))
    f <- phiMat %*% (rt(k, df=5)* sqrt(eigenvals))
    #scale each function to range=[0,1] for applying test function
    f <- f - min(f)
    return(drop(f/max(f)))
  }  
  #matrix of functional predictors x_i(t)
  xtmat <- t(replicate(n, fxt(t)))
  
  #vis:  matplot(t, t(xtmat), type="l", lty=1, col=8)
  #     sapply(1:5, function(i) lines(t, fxt(t), col=i, lwd=2))	
  # # quantile transform 
  # ztmat <- apply(xtmat, 2, function(x){
  #   (rank(x)-1)/(length(x)-1)
  # })  
  #apply(ztmat, 2, range)
  
  
  ## matrix of evaluation points
  tmat <- matrix(t, nrow=n, nc=ngrid, byrow=TRUE)
  
  Fxt <- test1(xtmat, tmat)    ## evaluate F(x(t),t) on grid
  
  intTrapez <- function(ft, t){
    ## Trapezoidal rule integration of f(t) on equally spaced points t
    diff(range(t))/(2*length(t)) * (2*sum(ft) - ft[1]  - ft[length(T)]) 
  }
  
  ## /int F(x(t),t) dt
  intFxt <- apply(Fxt, 1, intTrapez, t=t) 
  
  # noisy observations
#   snr <- 10
#   y <- intFxt + sd(intFxt)/snr * scale(rnorm(n))
  y <- intFxt+sigmae*rnorm(n)
  
  # define linear operator Li that does trapezoidal rule integration:
  # int_a^b f(t) = (b-a)/2*n  (f(t_1) + 2*f(t_2)+ ... + 2*f(t_n-1) + f(t_n))
  Li <- diff(range(t))/(2*ngrid) * c(1, rep(2, ngrid-2), 1)
  L <- matrix(Li, nr=n, nc=ngrid, byrow=T)

  return(list(y=y,xtmat=xtmat,tmat=tmat,L=L,Fxt=Fxt))
}



FGAMM4mp <- function(y, X1, X2, tvals1=seq(0,1,l=ncol(X)), tvals2=seq(0,1,l=ncol(X)),
                     fixed.effects, qr.con, family=gaussian(),
                  splinepars1=list(k=c(8,8),m=list(c(2,2),c(2,2)),extendXrange=.01),
                  splinepars2=list(k=c(8,8),m=list(c(2,2),c(2,2)),extendXrange=.01),
                  REML=TRUE,psplines=TRUE){
  
  N <- length(y)
  num.t1 <- length(tvals1)
  num.t2 <- length(tvals2)
  if (missing(fixed.effects)){
    num.fe <- 0
  }else if (is.null(dim(fixed.effects))){
    num.fe <- 1
  }else{
    num.fe <- ncol(fixed.effects)
  }
  if(!all.equal(tvals1,tvals2))
    stop('different grid values for the functional predictors not implemented yet')
  Kx1 <- splinepars1$k[1]
  Kt1 <- splinepars1$k[2]
  Kx2 <- splinepars2$k[1]
  Kt2 <- splinepars2$k[2]
  dx <- 2
  dt <- 2
  extendXrange <- ifelse(is.null(splinepars1$extendXrange),.00,
                         splinepars1$extendXrange)
  
  bbt1 <- create.bspline.basis(range(tvals1), nbasis=Kt1)
  bbtevals1 <- eval.basis(tvals1, bbt1)
  bbt2 <- create.bspline.basis(range(tvals2), nbasis=Kt2)
  bbtevals2 <- eval.basis(tvals1, bbt2)
  
  
  Xrange1 <-range(X1, na.rm=TRUE)
  extXrange1 <- Xrange1 #+ c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  Xrange2 <-range(X2, na.rm=TRUE)
  extXrange2 <- Xrange2 #+ c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx1 <- create.bspline.basis(extXrange1, nbasis=Kx1)
  L1 <- rep(1/num.t1, num.t1)
  bbx2 <- create.bspline.basis(extXrange2, nbasis=Kx2)
  L2 <- rep(1/num.t2, num.t2)
  
  if(psplines){
    DDx1 <- crossprod ({
      DD <- diag(Kx1)
      if (dx>0) 
        for (i in 1:dx) 
          DD <- diff(DD)
      DD
    }) 
    DDt1 <- crossprod({
      DD <- diag(Kt1)
      if (dt>0) 
        for (i in 1:dt) 
          DD <- diff(DD)
      DD
    })
    DDx2 <- crossprod ({
      DD <- diag(Kx2)
      if (dx>0) 
        for (i in 1:dx) 
          DD <- diff(DD)
      DD
    }) 
    DDt2 <- crossprod({
      DD <- diag(Kt2)
      if (dt>0) 
        for (i in 1:dt) 
          DD <- diff(DD)
      DD
    })
  }else{
    DDx1 <- getbasispenalty(bbx1, int2Lfd(2))
    DDt1 <- getbasispenalty(bbt1, int2Lfd(2))
    DDx2 <- getbasispenalty(bbx2, int2Lfd(2))
    DDt2 <- getbasispenalty(bbt2, int2Lfd(2))
  }
  
  
  # eigen decomp of marginal penalties:
  eDDx1 <- eigen(DDx1)
  nullx1 <- (Kx1-dx+1):Kx1
  if (dx>0) eDDx1$values[nullx1] <- 0
  
  eDDt1 <- eigen(DDt1)
  nullt1 <- (Kt1-dt+1):Kt1
  if (dt>0) eDDt1$values[nullt1] <- 0
  
  DDx2 <- crossprod({
    DD <- diag(Kx2)
    if(dx>0) for(i in 1:dx) DD <- diff(DD)
    DD
  }) 
  DDt2 <- crossprod({
    DD <- diag(Kt2)
    if(dt>0) for(i in 1:dt) DD <- diff(DD)
    DD
  })
  
  
  # eigen decomp of marginal penalties:
  eDDx2 <- eigen(DDx2)
  nullx2 <- (Kx2-dx+1):Kx2
  if (dx>0) eDDx2$values[nullx2] <- 0
  
  eDDt2 <- eigen(DDt2)
  nullt2 <- (Kt2-dt+1):Kt2
  if (dt>0) eDDt2$values[nullt2] <- 0
  
  if(TRUE){
    UDinvmat_t1 <- eDDt1$vectors[,1:(Kt1-2)]%*%diag(1/sqrt(eDDt1$values[1:(Kt1-2)]))
    Zt1 <- bbtevals1%*%UDinvmat_t1
    UDinvmat_x1 <- eDDx1$vectors[,1:(Kx1-2)]%*%diag(1/sqrt(eDDx1$values[1:(Kx1-2)]))
    UDinvmat_t2 <- eDDt2$vectors[,1:(Kt2-2)]%*%diag(1/sqrt(eDDt2$values[1:(Kt2-2)]))
    Zt2 <- bbtevals2%*%UDinvmat_t2
    nb1.t <- bbtevals1%*%eDDt1$vectors[,(Kt1-1):Kt1]
    nb2.t <- bbtevals2%*%eDDt2$vectors[,(Kt2-1):Kt2]
    UDinvmat_x2 <- eDDx2$vectors[,1:(Kx2-2)]%*%diag(1/sqrt(eDDx2$values[1:(Kx2-2)]))
    Design <- matrix(unlist(lapply(1:N,function(i){
      bxvals1 <- eval.basis(X1[i, ], bbx1)
      Zx1 <- bxvals1%*%UDinvmat_x1
      bxvals2 <- eval.basis(X2[i, ], bbx2)
      Zx2 <- bxvals2%*%UDinvmat_x2
      nb1.x <- bxvals1%*%eDDx1$vectors[, (Kx1-1):Kx1]
      nb2.x <- bxvals2%*%eDDx2$vectors[, (Kx2-1):Kx2]
      #  Bxi<-	cbind( 1+numeric(num.t1), xvals, xvals*tvals, #=Z_0
      #               Zx, Zt*xvals, Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      #Bxi<-  cbind( 1+numeric(num.t1), xvals, xvals*tvals, Zt*xvals) #=Z_p
      Bxi<-  cbind( nb1.t[, 1]*nb1.x[, 1], nb1.t[, 2]*nb1.x[, 1], nb1.t[, 1]*nb1.x[, 2],  # =Z1_0
                    nb2.t[, 1]*nb2.x[, 1], nb2.t[, 2]*nb2.x[, 1], nb2.t[, 1]*nb2.x[, 2],
                    Zt1*X1[i, ], Zx1,Zx1*tvals1, 
                    Zx1[,rep(1:(Kx1-2),e=Kt1-2)]*Zt1[,rep(1:(Kt1-2),t=Kx1-2)],
                    Zt2*X2[i, ], Zx2, Zx2*tvals2,
                    Zx2[,rep(1:(Kx2-2),e=Kt2-2)]*Zt2[,rep(1:(Kt2-2),t=Kx2-2)]) #=Z_p
      #browser()
      crossprod(L1, Bxi)
    })), nr=N, nc=Kx1*Kt1+Kx2*Kt2-Kt1-Kt2+2, byrow=TRUE)
    nupp <- 2*dx*dt-2
  }
  if(num.fe)
    Design <- cbind(fixed.effects, Design)
  nupp <- nupp + num.fe
  upind <- 1:nupp
  dims <- c(nupp, Kt1-2, 2*(Kx1-2), (Kx1-2)*(Kt1-2), Kt2-2, 2*(Kx2-2), (Kx2-2)*(Kt2-2)) #c(nupp,Kx-2,Kt-2,Kx-2,(Kx-2)*(Kt-2))
  dnames <- c('para', 'x1.Zt1', 'onet1.Zx1', 'Zx1.Zt1', 'x2.Zt2', 'onet2.Zx2', 'Zx2.Zt2') #c('para','Zx','x.Zt','Zx.t','Zx.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  if (!missing(qr.con)){
    dat$X <- t(qr.qty(qr.con, t(Design[,1:nupp])))[,-(1:2)]
  }else{
    dat$X <- Design[,1:nupp]
  }
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)

#   if (TRUE){
#     if (num.fe){
#       names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), 
#                                    't11.x11', 't12.x11', 't11.x12',
#                                    't21.x21', 't22.x21', 't21.x22')
#     }else{
#       names(tempfit$fr$fixef) <- c('t11.x11', 't12.x11', 't11.x12',
#                                    't21.x21', 't22.x21', 't21.x22')
#     }
#   }
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:length(dnames)]
  snames <- c('x1*f(t1)', 'g1(x1)+t1*g2(x1)', 'h(x1,t1)', 
              'x2*f(t2)', 'g1(x2)+t2*g2(x2)', 'h(x2,t2)')
  for (i in 1:(length(dnames)-1)){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  ret <- list()
  ret$lme4fit <- do.call(lme4:::lmer_finalize,tempfit)
  ret$bbx1 <- bbx1
  ret$bbt1 <- bbt1
  ret$UDinvmat_x1 <- UDinvmat_x1
  ret$UDinvmat_t1 <- UDinvmat_t1
  ret$bbx2 <- bbx2
  ret$bbt2 <- bbt2
  ret$UDinvmat_x2 <- UDinvmat_x2
  ret$UDinvmat_t2 <- UDinvmat_t2
  ret$null.ev.t1 <- eDDt1$vectors[,(Kt1-1):Kt1]
  ret$null.ev.t2 <- eDDt2$vectors[,(Kt2-1):Kt2]
  ret$null.ev.x1 <- eDDx1$vectors[,(Kx1-1):Kx1]
  ret$null.ev.x2 <- eDDx2$vectors[,(Kx2-1):Kx2]
  ret$parinds <- parinds
  ret$dnames <- dnames
  
  vc <- lme4:::VarCorr(ret$lme4fit)
  ret$sigma2 <- as.numeric(attr(vc,'sc'))^2
  sp <- numeric(6)
  for (i in 1:(length(dnames)-1)){
    if (as.numeric(vc[[i]]) >0 ){
      sp[i] <- ret$sigma2/as.numeric(vc[[i]])
    }else{
      sp[i] <- 1e15
    }
  }
  names(sp) <- names(vc)
  sp <- sp[dnames[-1]]
  names(sp) <- snames
  ret$sp <- sp
  blups <- unlist(as(lme4::ranef(ret$lme4fit), 'vector'))
  ind <- unlist(sapply(snames,function(x){grep(x,
                                               substr(names(blups), 1, nchar(x)), fixed=T)}))
  blups <- blups[ind]
  names(blups) <- unlist( mapply(FUN=function(y,x)rep(y,e=x), y=snames, x=dims[-1]) )  
  ret$snames <- snames  
  ret$blups <- blups
  ret$fixef <- qr.qy(qr.con,c(0,0,lme4::fixef(ret$lme4fit)))
  names(ret$fixef) <- if (num.fe){
              c(outer('fe.', 1:num.fe, paste, sep=''), 'nX1.1', 'nX1.2', 'nX1.3', 'nX2.1', 
                     'nX2.2', 'nX2.3')
              }else{
                c('nX1.1', 'nX1.2', 'nX1.3', 'nX2.1', 'nX2.2', 'nX2.3')
              }
  ret$coefficients <- c(ret$fixef,blups)
  ret$residuals <- residuals(ret$lme4fit)
  ret$fitted.values <- dat$y - ret$residuals
  
  return(ret)
}

predict.fgamm <- function(fit, X, fe, tvals){
  numT <- length(tvals)
  btvals <-  eval.basis(tvals,fit$bbt)
  Kt <- fit$bbt$nb
  Kx <- fit$bbx$nb
  X[X > fit$bbx$range[2]] <- fit$bbx$range[2]
  X[X < fit$bbx$range[1]] <- fit$bbx$range[1]
  
  Zt <- btvals%*%fit$UDinvmat_t
  L <- matrix(1/numT,numT,1)
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,fit$bbx)%*%fit$UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals,  # =Z_0
                  Zt*xvals, Zx,Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
    crossprod(L,Bxi)
  }) )  
  if(!missing(fe))
    Design <- cbind(fe,Design)
  
  return(drop(Design%*%fit$coef))
}

predict.fgamm2mp <- function(fit, X1, X2, fe, tvals1, tvals2){
  N <- nrow(X1)
  num.t1 <- length(tvals1)
  num.t2 <- length(tvals2)
  btvals1 <-  eval.basis(tvals1, fit$bbt1)
  btvals2 <-  eval.basis(tvals2, fit$bbt2)
  Zt1 <- btvals1%*%fit$UDinvmat_t1
  Zt2 <- btvals1%*%fit$UDinvmat_t2
  Kt1 <- fit$bbt1$nb
  Kx1 <- fit$bbx1$nb
  Kt2 <- fit$bbt2$nb
  Kx2 <- fit$bbx2$nb
  X1[X1 > fit$bbx1$range[2]] <- fit$bbx1$range[2]
  X1[X1 < fit$bbx1$range[1]] <- fit$bbx1$range[1]
  X2[X2 > fit$bbx2$range[2]] <- fit$bbx2$range[2]
  X2[X2 < fit$bbx2$range[1]] <- fit$bbx2$range[1]
  
  nb1.t <- btvals1%*%fit$null.ev.t1
  nb2.t <- btvals2%*%fit$null.ev.t2
  L1 <- matrix(1/num.t1,num.t1,1)
  
  Design <- matrix(unlist(lapply(1:N,function(i){
    bxvals1 <- eval.basis(X1[i, ], fit$bbx1)
    Zx1 <- bxvals1%*%fit$UDinvmat_x1
    bxvals2 <- eval.basis(X2[i, ], fit$bbx2)
    Zx2 <- bxvals2%*%fit$UDinvmat_x2
    nb1.x <- bxvals1%*%fit$null.ev.x1
    nb2.x <- bxvals2%*%fit$null.ev.x2
    #  Bxi<-  cbind( 1+numeric(num.t1), xvals, xvals*tvals, #=Z_0
    #               Zx, Zt*xvals, Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
    #Bxi<-  cbind( 1+numeric(num.t1), xvals, xvals*tvals, Zt*xvals) #=Z_p
    Bxi<-  cbind( nb1.t[, 1]*nb1.x[, 1], nb1.t[, 2]*nb1.x[, 1], nb1.t[, 1]*nb1.x[, 2],  # =Z1_0
                  nb2.t[, 1]*nb2.x[, 1], nb2.t[, 2]*nb2.x[, 1], nb2.t[, 1]*nb2.x[, 2],
                  Zt1*X1[i, ], Zx1,Zx1*tvals1, 
                  Zx1[,rep(1:(Kx1-2),e=Kt1-2)]*Zt1[,rep(1:(Kt1-2),t=Kx1-2)],
                  Zt2*X2[i, ], Zx2, Zx2*tvals2,
                  Zx2[,rep(1:(Kx2-2),e=Kt2-2)]*Zt2[,rep(1:(Kt2-2),t=Kx2-2)]) #=Z_p
    #browser()
    crossprod(L1, Bxi)
  })), nr=N, nc=Kx1*Kt1+Kx2*Kt2-Kt1-Kt2+2, byrow=TRUE)

  if(!missing(fe))
    Design <- cbind(fe, Design)
  
  return(drop(Design%*%fit$coef))

}

# function for testing using standard tensor product spline construction

TestLinearitySimple <- function(y,X,tvals=seq(0,1,l=ncol(X)), family=gaussian(),
                    splinepars=list(k=c(8,8), m=list(c(2,2),c(2,2)), extendXrange=.01, bs='ps'),
                    REML=TRUE, removeCon=FALSE, cond.S=FALSE, fixed.effects, lambda.t){
  
  N <- length(y)
  numT <- length(tvals)
  if (missing(fixed.effects)){
    num.fe <- 1
  }else if (is.null(dim(fixed.effects))){
    num.fe <- 2
  }else{
    num.fe <- ncol(fixed.effects)+1
  }
  Kx <- splinepars$k[1]
  Kt <- splinepars$k[2]
  bs <- splinepars$bs
  dx <- 2
  dt <- 2
 # extendXrange <- ifelse(is.null(splinepars$extendXrange),.01,
#                         splinepars$extendXrange)
  L <- matrix(1/numT, N, numT)
  tmat <- matrix(tvals, N, numT, byrow=TRUE)
#   if(missing(lambda.t)){
#     if(missing(fixed.effects)){
#       fit <- gam(y~te(X, tmat, by=L, bs=bs, k=c(Kx, Kt)), method='REML') 
#     }else{
#       fit <- gam(y~fixed.effects+te(X, tmat, by=L, bs=bs, k=c(Kx, Kt)), method='REML')
#     }
#   }else{
#     
#   }
#  browser()
  args <- list(method='REML')
  args$fit <- switch((missing(lambda.t) || is.null(lambda.t))+1, TRUE, TRUE)
  args$formula <- switch((missing(fixed.effects) || is.null(fixed.effects))+1,
                         y~fixed.effects+te(X, tmat, by=L, bs=bs, k=c(Kx, Kt)),  # FALSE
                         y~te(X, tmat, by=L, bs=bs, k=c(Kx, Kt)))  # TRUE
  fit <- do.call(gam, args)
  
  sm <- fit$smooth[[1]]
  
  P1 <- sm$S[[1]]
  P2 <- sm$S[[2]]
  
  model.mat <- model.matrix(fit)
  if(removeCon){
    qrC <- attr(sm,'qrc')
    model.mat <- t(qr.qty(qrC, t(cbind(rep(0, nrow(model.mat)), model.mat[, -1]))))  # remove int.
    model.mat <- cbind(rep(1,nrow(model.mat)), model.mat)  # add back in int.
    P1 <- t(qr.qty(qrC, t(as.matrix(bdiag(0, P1)))))
    P1 <- qr.qy(qrC, P1)
    P2 <- t(qr.qty(qrC, t(as.matrix(bdiag(0, P2)))))
    P2 <- qr.qy(qrC, P2)
  #  S.mat.no.constr <- fit$sp[1]*P1.no.constr+fit$sp[2]*P2.no.constr
  #  eigen(S.mat.no.constr)$val
  }
  
  P1 <- as.matrix(bdiag(rep(0, num.fe), P1))
  P2 <- as.matrix(bdiag(rep(0, num.fe), P2))
  
#  S.mat <- fit$sp[1]*sm$S[[1]] + fit$sp[2]*sm$S[[2]]
  if(cond.S){
    eps <- .Machine$double.eps^.7*eigen(fit$sp[1]*P1 + fit$sp[2]*P2, 
                                      TRUE, TRUE)$val[1]
    precond.eigen <- eigen(P1/norm(P1, 'F') + P2/norm(P1, 'F'))
#   P1 <- crossprod(precond.eigen$vec[, precond.eigen$val > sqrt(.Machine$double.eps)], P1)%*%
#     precond.eigen$vec[, precond.eigen$val > sqrt(.Machine$double.eps)]
#   P2 <- crossprod(precond.eigen$vec[, precond.eigen$val > sqrt(.Machine$double.eps)], P2)%*%
#     precond.eigen$vec[, precond.eigen$val > sqrt(.Machine$double.eps)]
    P1 <- crossprod(precond.eigen$vec[, precond.eigen$val > eps], P1)%*%
      precond.eigen$vec[, precond.eigen$val > eps]
    P2 <- crossprod(precond.eigen$vec[, precond.eigen$val > eps], P2)%*%
      precond.eigen$vec[, precond.eigen$val > eps]
    model.mat <- model.mat%*%precond.eigen$vec[, precond.eigen$val > eps]
  }
  # eigen(fit$sp[1]*P1+fit$sp[2]*P2, TRUE, TRUE)$val # full rank
  # eigen(P1, TRUE, TRUE)$val  # not even close to full rank
####### stuff from Wood (2011) that I don't think you need
#   P.mat1 <- diag(1/sqrt(diag(P1)))
#   P.mat2 <- diag(1/sqrt(diag(P2)))
#   P1.precond <- P.mat1%*%P1%*%P.mat1
#   P2.precond <- P.mat2%*%P2%*%P.mat2  
#  P1.chol <- chol(P1.precond, pivot=TRUE)
#  piv.ord1 <- order(attr(P1.chol, 'pivot'))  
#  sqrtP1 <- P1.chol[, piv.ord1]%*%diag(1/diag(P.mat1))
  
  P2.chol <- suppressWarnings(chol(P2, pivot=TRUE))
  piv.ord2 <- order(attr(P2.chol, 'pivot'))
  sqrtP2 <- P2.chol[, piv.ord2]
  
  # all.equal(P1, crossprod(sqrt.P1))
  
 # P1 <- as.matrix(bdiag(rep(0, num.fe), P1))  # add intercept
  #P2 <- cbind(matrix(0, nrow(P1), num.fe), P1)  # as.matrix(bdiag(rep(0, num.fe), P2))
#  browser()
#   eigP2 <- eigen(P2)
#   pos.eval <- eigP2$val*(!eigP2$val < sqrt(.Machine$double.eps))  # truncate near-zero eig. vals
#   sqrtP2 <- switch((missing(lambda.t) || is.null(lambda.t))+1,
#                    sqrt(lambda.t)*eigP2$vec%*%diag(sqrt(pos.eval)),  # FALSE
#                    sqrt(fit$sp[2])*eigP2$vec%*%diag(sqrt(pos.eval)))  # TRUE
#   sqrtP2 <-  cbind(matrix(0, nrow(sqrtP2), num.fe), sqrtP2) # add zero columns for fixed effects
  augy <- c(y, rep(0,nrow(P2)))
  augX <- rbind(model.mat, sqrtP2)
  
  eigP1 <- eigen(P1)
 # pos.eval <- eigP1$val*(!eigP1$val < sqrt(.Machine$double.eps))  # truncate near-zero eig. vals
 #  sqrtP1 <- sqrt(fit$sp[1])*eigP1$vec%*%diag(sqrt(pos.eval))
 # augy2 <- c(y, rep(0,nrow(P2)+nrow(P1)))
 # augX2 <- rbind(model.mat, sqrtP2, sqrtP1)
 # logLik(lm(augy2~-1+augX2))
  zero.ev.ind <- eigP1$val < sqrt(.Machine$double.eps)
  newX <- augX%*%eigP1$vec[, zero.ev.ind]
  Z <- augX%*%eigP1$vec[, !zero.ev.ind]%*%diag(1/sqrt(eigP1$val[!zero.ev.ind]))
  #Z <- Z/norm(Z, 'F')
  #newX <- newX/norm(newX, 'F')
  # need constraint on Z
  
  # args <- list(formula=y ~ s(tmat, by=I(L*X), k=Kt, bs=bs), method='REML')
  if(!(missing(lambda.t) || is.null(lambda.t)))
    args$sp <- lambda.t
  # flm.fit <- do.call(gam, args)
  # rlrt.obs <- max(0, 2 * (logLik(fit)[1] - logLik(flm.fit)[1]))
  rownames(newX) <- NULL
  rownames(Z) <- NULL
  dat <- list(y=augy, X=newX)
 #  browser()
  nupp <- ncol(newX)
  
  upind <- 1:nupp
  dims <- c(nupp,ncol(Z))
  dnames <- c('para','np')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  augN <- length(augy)
  # browser()
  form <- 'y~0+X'
  dat[[dnames[2]]] <- factor(rep(1:dims[2],length=augN))
  form <- paste(form,'+ (1|',dnames[2],')')
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2]
  snames <- 'np'
    k <- ind[sn[1]==tn]
    tempfit$FL$trms[[1]]$A <- tempfit$FL$trms[[1]]$Zt <- 
      as(t(Z),'dgCMatrix')
    attr(tempfit$FL$trms[[1]]$ST,'dimnames') <- list(snames[1],snames[1])

  #browser()
  m <- do.call(lme4:::lmer_finalize,tempfit)

  #sim.dist <- RLRTSim(newX, Z, qr(newX), diag(ncol(Z)))
  #p.value <- ifelse(rlrt.obs==0, 1, mean(rlrt.obs < sim.dist))
  p.value <- exactRLRT(m)$p
  ret <- c(pval=p.value, reject90=p.value<.1, reject95=p.value<.05, reject99=p.value<.01)
  return(ret)
}
#     
#   bbt <- create.bspline.basis(range(xind), Kt)
#   bbx <- create.bspline.basis(range(X), Kx)  
#   btevals <- eval.basis(xind, bbt)  
#   L <- matrix(1/J, J, 1)
#   Design <- t(apply(X, 1, function(x){
#     bxevals <- eval.basis(x, bbx)
#       crossprod(L, bxevals[, rep(1:Kx, t=Kt)]*btevals[, rep(1:Kt, e=Kx)])
#   }))
#   
#   C.mat <- c(0, colSums(Design))
#   Design <- cbind(rep(1,N), Design)

testfunOracle <- function(y,X,tvals=seq(0,1,l=ncol(X)), family=gaussian(),
                    splinepars=list(k=c(8,8), m=list(c(2,2),c(2,2)), extendXrange=.01),
                    REML=TRUE, oracle=0, formratio=FALSE, sig2FLM=NULL, ratio=NULL,
                    psplines=TRUE, unknown.fe=FALSE, fixed.effects){
  
  N <- length(y)
  numT <- length(tvals)
  if (missing(fixed.effects)){
    num.fe <- 0
  }else if (is.null(dim(fixed.effects))){
    num.fe <- 1
  }else{
    num.fe <- ncol(fixed.effects)
  }
  Kx <- splinepars$k[1]
  Kt <- splinepars$k[2]
  dx <- 2
  dt <- 2
  extendXrange <- ifelse(is.null(splinepars$extendXrange),.01,
                         splinepars$extendXrange)
  
  bbt <- create.bspline.basis(range(tvals), nbasis=Kt)
  bbtevals <- eval.basis(tvals, bbt)
  
  
  Xrange <-range(X, na.rm=TRUE)
  extXrange <- Xrange + c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx <- create.bspline.basis(extXrange, nbasis=Kx)
  L <- rep(1/numT, numT)
  
  if(psplines){
    DDx <- crossprod ({
      DD <- diag(Kx)
      if (dx>0) 
        for (i in 1:dx) 
          DD <- diff(DD)
      DD
    }) 
    DDt <- crossprod({
      DD <- diag(Kt)
      if (dt>0) 
        for (i in 1:dt) 
          DD <- diff(DD)
      DD
    })
  }else{
    DDx <- getbasispenalty(bbx, int2Lfd(2))
    DDt <- getbasispenalty(bbt, int2Lfd(2))
  }
  
  # eigen decomp of marginal penalties:
  eDDx <- eigen(DDx)
  nullx <- (Kx-dx+1):Kx
  if(dx>0) eDDx$values[nullx] <- 0
  
  eDDt <- eigen(DDt)
  nullt <- (Kt-dt+1):Kt
  if(dt>0) eDDt$values[nullt] <- 0
  
  UDinvmat_t <- eDDt$vectors[,1:(Kt-2)]%*%diag(1/sqrt(eDDt$values[1:(Kt-2)]))
  Zt <- bbtevals%*%UDinvmat_t
  UDinvmat_x <- eDDx$vectors[,1:(Kx-2)]%*%diag(1/sqrt(eDDx$values[1:(Kx-2)]))
  
  ##########################################################################
  # 1)
  # first fit no nuisance parm model with only one var comp. m to exactRLRT or exactLRT for test 1
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi <-  cbind(Zx, Zx*tvals) #=Z_p
    if(unknown.fe)
      Bxi <- cbind(1+numeric(numT), xvals, xvals*tvals, Bxi)
    crossprod(L,Bxi)
  }) )  
  
  nupp <- dx*dt-1
  if (num.fe)
    Design <- cbind(fixed.effects,Design)
  
  nupp <- nupp + num.fe
  
  # upind <- 1:nupp
  if(unknown.fe){
    dims <- c(nupp, 2*(Kx-2))
    dnames <- c('para', 'onet.Zx')
  }else{
    dims <- c(2*(Kx-2))
    dnames <- c('onet.Zx')
  }
  n.dims <- unknown.fe + 1
  npars <- sum(dims)
  # npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  #browser()
  dat <- list()
  dat$y <- y
  if(unknown.fe){
    dat$X <- Design[, 1:nupp]
    form <- 'y~0+X'
  }else{
    dat$poop <- rnorm(length(y))
    form <- 'y~0+poop'
  }
  i <- n.dims
  dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
  form <- paste(form,'+ (1|',dnames[i],')')
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
#   if (num.fe){
#     names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
#   }else{
#     names(tempfit$fr$fixef) <- c('1','x','x*t')
#   }
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[n.dims]
  snames <- 'f_1(x)+t*f_2(x)'
  i <- 1
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  
  #browser()
  m1 <- do.call(lme4:::lmer_finalize,tempfit)
  
  ##################################################################################################
  # 2)
  # first fit no nuisance parm model with only one var comp. m arg to exactRLRT or exactLRT for test 2
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind(Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
    if(unknown.fe)
      Bxi <- cbind(1+numeric(numT), xvals, xvals*tvals, Bxi)
    crossprod(L,Bxi)
  }) )  
  
  # nupp <- dx*dt-1
  
 #  upind <- 1:nupp
  if(unknown.fe){
    dims <- c(nupp, (Kx-2)*(Kt-2))
    dnames <- c('para','Zx.Zt')
  }else{
    dims <- c((Kx-2)*(Kt-2))
    dnames <- c('Zx.Zt')
  }
  npars <- sum(dims)
 # npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
 # dat$X <- Design[,1:nupp]
  # browser()
  if(unknown.fe){
    dat$X <- Design[, 1:nupp]
    form <- 'y~0+X'
  }else{
    dat$poop <- rnorm(length(y))
    form <- 'y~0+poop'
  }
  
  i <- n.dims
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
 # names(tempfit$fr$fixef) <- c('1','x','x*t')
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[n.dims]
  snames <- 'f(x,t)'
  i <- 1
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  
  m2 <- do.call(lme4:::lmer_finalize,tempfit)
  
 #  browser()
  var1.0 <- lme4:::VarCorr(m1)[[1]]  # onet.Zx
  var2.0 <- lme4:::VarCorr(m2)[[1]]  # Zt.Zx

  if(var1.0==0 && var2.0==0){
    return(c(NA,NA,0,0,0))
  }else if(var1.0>0 && var2.0==0){
    test1r <- exactRLRT(m1)
    pval1r <- test1r$p
    pval2r <- NA
    reject90 <- pval1r < .05 
    reject95 <- pval1r < .025
    reject99 <- pval1r < .005 
  }else if(var1.0==0 & var2.0 >0){
    test2r <- exactRLRT(m2)
    pval1r <- NA
    pval2r <- test2r$p
    reject90 <- pval2r < .05
    reject95 <- pval2r < .025
    reject99 <- pval2r < .005
  }else{
    test1r <- exactRLRT(m1)
    pval1r <- test1r$p
    test2r <- exactRLRT(m2)
    pval2r <- test2r$p
    reject90 <- pval2r < .05  && pval1r < .05
    reject95 <- pval2r < .025 && pval1r < .025
    reject99 <- pval2r < .005 && pval1r < .005    
  }
  
  ret <- c(
    pval_np=pval1r,pval_pp=pval2r,reject90=reject90,reject95=reject95,reject99=reject99
  )
  
  return(ret)
}

TestKnowEqualVC <- function(y,X,tvals=seq(0,1,l=ncol(X)), family=gaussian(),
                          splinepars=list(k=c(8,8), m=list(c(2,2),c(2,2)), extendXrange=.01),
                          REML=TRUE, oracle=0, formratio=FALSE, sig2FLM=NULL, ratio=NULL,
                          psplines=TRUE, removeCon = FALSE, fixed.effects){
  
  N <- length(y)
  numT <- length(tvals)
  if (missing(fixed.effects)){
    num.fe <- 0
  }else if (is.null(dim(fixed.effects))){
    num.fe <- 1
  }else{
    num.fe <- ncol(fixed.effects)
  }
  Kx <- splinepars$k[1]
  Kt <- splinepars$k[2]
  dx <- 2
  dt <- 2
  extendXrange <- ifelse(is.null(splinepars$extendXrange),.01,
                         splinepars$extendXrange)
  
  bbt <- create.bspline.basis(range(tvals), nbasis=Kt)
  bbtevals <- eval.basis(tvals, bbt)
  
  
  Xrange <-range(X, na.rm=TRUE)
  extXrange <- Xrange + c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx <- create.bspline.basis(extXrange, nbasis=Kx)
  L <- rep(1/numT, numT)
  
  if(psplines){
    DDx <- crossprod ({
      DD <- diag(Kx)
      if (dx>0) 
        for (i in 1:dx) 
          DD <- diff(DD)
      DD
    }) 
    DDt <- crossprod({
      DD <- diag(Kt)
      if (dt>0) 
        for (i in 1:dt) 
          DD <- diff(DD)
      DD
    })
  }else{
    DDx <- getbasispenalty(bbx, int2Lfd(2))
    DDt <- getbasispenalty(bbt, int2Lfd(2))
  }
  
  # eigen decomp of marginal penalties:
  eDDx <- eigen(DDx)
  nullx <- (Kx-dx+1):Kx
  if(dx>0) eDDx$values[nullx] <- 0
  
  eDDt <- eigen(DDt)
  nullt <- (Kt-dt+1):Kt
  if(dt>0) eDDt$values[nullt] <- 0
  
  UDinvmat_t <- eDDt$vectors[,1:(Kt-2)]%*%diag(1/sqrt(eDDt$values[1:(Kt-2)]))
  Zt <- bbtevals%*%UDinvmat_t
  UDinvmat_x <- eDDx$vectors[,1:(Kx-2)]%*%diag(1/sqrt(eDDx$values[1:(Kx-2)]))
  
  ##########################################################################
  # 1)
  # first fit no nuisance parm model with only one var comp. m to exactRLRT or exactLRT for test 1
  # known that sig2=sig3
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind(1+numeric(numT), xvals, xvals*tvals, 
                 Zx, Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
    crossprod(L,Bxi)
  }) )  
  
  nupp <- dx*dt-1
  if (num.fe)
    Design <- cbind(fixed.effects,Design)
  
  nupp <- nupp + num.fe
  
  upind <- 1:nupp
  dims <- c(nupp, 2*(Kx-2)+(Kx-2)*(Kt-2))
  dnames <- c('para', 'nonFLM')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  # dat$poop <- rnorm(length(y))
  dat$X <- Design[,1:nupp]
   #browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  #browser()
    if (num.fe){
      names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
    }else{
      names(tempfit$fr$fixef) <- c('1','x','x*t')
    }
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2]
  snames <- 'nonFLM'
  for(i in 1:1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  #browser()
  m <- do.call(lme4:::lmer_finalize,tempfit)
  
  ###########################################################################################
  # 2)
  # fit full model assuming sig2=sig3 (mA arg for test two)
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zt*xvals, Zx, Zx*tvals,
                  Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
    crossprod(L,Bxi)
  }) ) 
  
  nupp <- dx*dt-1
  
  if (num.fe)
    Design <- cbind(fixed.effects,Design)
  
  nupp <- nupp + num.fe
  
  upind <- 1:nupp
  dims <- c(nupp, Kt-2, 2*(Kx-2) + (Kx-2)*(Kt-2))
  dnames <- c('para','FLM','nonFLM')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  if (num.fe){
    names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
  }else{
    names(tempfit$fr$fixef) <- c('1','x','x*t')
  }
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:3]
  snames <- c('FLM','nonFLM')
  for(i in 1:2){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  mA <- do.call(lme4:::lmer_finalize,tempfit)
  
  ###########################################################################################
  # 3)
  # fit model under null (FLM) with one variance component (m0 arg to exactRLRT) same for test 1 and 2
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zt*xvals) #=Z_p
    crossprod(L,Bxi)
  }) ) 
  
  nupp <- dx*dt-1
  
  if (num.fe)
    Design <- cbind(fixed.effects,Design)
  
  nupp <- nupp + num.fe
  
  upind <- 1:nupp
  dims <- c(nupp,Kt-2)
  dnames <- c('para','x.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  if (num.fe){
    names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
  }else{
    names(tempfit$fr$fixef) <- c('1','x','x*t')
  }
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:2]
  snames <- 'FLM'
  for(i in 1:1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  m0 <- do.call(lme4:::lmer_finalize,tempfit)
  # browser()
  var1 <- lme4:::VarCorr(m)[[1]]  #onet.Zx
  var2 <- lme4:::VarCorr(mA)[[1]]  #Zt.Zx
  var0 <- lme4:::VarCorr(m0)[[1]]
  if(var0==0)
    return(rep(NA,4))
  if(var1==0 || var2==0){
    return(c(NA,0,0,0))
  }else{
    test <- exactRLRT(m,mA,m0)
    pval <- test$p
    reject90 <- pval < .1 
    reject95 <- pval < .05
    reject99 <- pval < .01 
  }
  
  ret <- c(
    pval=pval, reject90=reject90, reject95=reject95, reject99=reject99
  )
  
  return(ret)
}

TestTwoVCBOOT <- function(y,X,tvals=seq(0,1,l=ncol(X)), family=gaussian(),
                            splinepars=list(k=c(8,8), m=list(c(2,2),c(2,2)), extendXrange=.01),
                            REML=TRUE, oracle=0, formratio=FALSE, sig2FLM=NULL, ratio=NULL,
                            psplines=TRUE, removeCon = FALSE, BOOT=FALSE, n.sims.boot=500,
                            fixed.effects){
  
  N <- length(y)
  numT <- length(tvals)
  if (missing(fixed.effects)){
    num.fe <- 0
  }else if (is.null(dim(fixed.effects))){
    num.fe <- 1
  }else{
    num.fe <- ncol(fixed.effects)
  }
  Kx <- splinepars$k[1]
  Kt <- splinepars$k[2]
  dx <- 2
  dt <- 2
  extendXrange <- ifelse(is.null(splinepars$extendXrange),.01,
                         splinepars$extendXrange)
  
  bbt <- create.bspline.basis(range(tvals), nbasis=Kt)
  bbtevals <- eval.basis(tvals, bbt)
  
  
  Xrange <-range(X, na.rm=TRUE)
  extXrange <- Xrange + c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx <- create.bspline.basis(extXrange, nbasis=Kx)
  L <- rep(1/numT, numT)
  
  if(psplines){
    DDx <- crossprod ({
      DD <- diag(Kx)
      if (dx>0) 
        for (i in 1:dx) 
          DD <- diff(DD)
      DD
    }) 
    DDt <- crossprod({
      DD <- diag(Kt)
      if (dt>0) 
        for (i in 1:dt) 
          DD <- diff(DD)
      DD
    })
  }else{
    DDx <- getbasispenalty(bbx, int2Lfd(2))
    DDt <- getbasispenalty(bbt, int2Lfd(2))
  }
  
  # eigen decomp of marginal penalties:
  eDDx <- eigen(DDx)
  nullx <- (Kx-dx+1):Kx
  if(dx>0) eDDx$values[nullx] <- 0
  
  eDDt <- eigen(DDt)
  nullt <- (Kt-dt+1):Kt
  if(dt>0) eDDt$values[nullt] <- 0
  
  UDinvmat_t <- eDDt$vectors[,1:(Kt-2)]%*%diag(1/sqrt(eDDt$values[1:(Kt-2)]))
  Zt <- bbtevals%*%UDinvmat_t
  UDinvmat_x <- eDDx$vectors[,1:(Kx-2)]%*%diag(1/sqrt(eDDx$values[1:(Kx-2)]))
  
  ##########################################################################
  # 1)
  # first fit no nuisance parm model with only one var comp. m to exactRLRT or exactLRT for test 1
  # known that sig2=sig3
  if(!BOOT){
    Design <- t( apply(X,1,function(xvals){
      Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
      Bxi<-  cbind(1+numeric(numT), xvals, xvals*tvals, 
                   Zx, Zx*tvals, Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
      crossprod(L,Bxi)
    }) )  
    
    nupp <- dx*dt-1
    if (num.fe)
      Design <- cbind(fixed.effects,Design)
    
    nupp <- nupp + num.fe
    
    upind <- 1:nupp
    dims <- c(nupp, 2*(Kx-2)+(Kx-2)*(Kt-2))
    dnames <- c('para', 'nonFLM')
    npars <- sum(dims)
    npp <- npars-nupp
    parinds <- 1:npars
    names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
    
    colnames(Design) <- names(parinds)
    
    dat <- list()
    dat$y <- y
    # dat$poop <- rnorm(length(y))
    dat$X <- Design[,1:nupp]
    #browser()
    form <- 'y~0+X'
    for(i in 2:length(dims)){
      dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
      form <- paste(form,'+ (1|',dnames[i],')')
    }
    
    form <-as.formula(form)
    tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
    #browser()
    if (num.fe){
      names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
    }else{
      names(tempfit$fr$fixef) <- c('1','x','x*t')
    }
    
    tn <- names(tempfit$FL$fl)
    tn <- tn[attr(tempfit$FL$fl,'assign')]
    ind <- 1:length(tn)
    sn <- dnames[2]
    snames <- 'nonFLM'
    for(i in 1:1){
      k <- ind[sn[i]==tn]
      tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
        as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
      attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
    }
    #browser()
    m <- do.call(lme4:::lmer_finalize,tempfit)
  }
  
  ###########################################################################################
  # 2)
  # fit full model assuming sig2=sig3 (mA arg for test two)
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zt*xvals, Zx, Zx*tvals,
                  Zx[,rep(1:(Kx-2),e=Kt-2)]*Zt[,rep(1:(Kt-2),t=Kx-2)]) #=Z_p
    crossprod(L,Bxi)
  }) ) 
  
  nupp <- dx*dt-1
  
  if (num.fe)
    Design <- cbind(fixed.effects,Design)
  
  nupp <- nupp + num.fe
  
  upind <- 1:nupp
  dims <- c(nupp, Kt-2, 2*(Kx-2) + (Kx-2)*(Kt-2))
  dnames <- c('para','FLM','nonFLM')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  if (num.fe){
    names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
  }else{
    names(tempfit$fr$fixef) <- c('1','x','x*t')
  }
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:3]
  snames <- c('FLM','nonFLM')
  for(i in 1:2){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  mA <- do.call(lme4:::lmer_finalize,tempfit)
  
  ###########################################################################################
  # 3)
  # fit model under null (FLM) with one variance component (m0 arg to exactRLRT) same for test 1 and 2
  
  Design <- t( apply(X,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tvals, #=Z_0
                  Zt*xvals) #=Z_p
    crossprod(L,Bxi)
  }) ) 
  
  nupp <- dx*dt-1
  
  if (num.fe)
    Design <- cbind(fixed.effects,Design)
  
  nupp <- nupp + num.fe
  
  upind <- 1:nupp
  dims <- c(nupp,Kt-2)
  dnames <- c('para','x.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  
  colnames(Design) <- names(parinds)
  
  dat <- list()
  dat$y <- y
  dat$X <- Design[,1:nupp]
  # browser()
  form <- 'y~0+X'
  for(i in 2:length(dims)){
    dat[[dnames[i]]] <- factor(rep(1:dims[i],length=N))
    form <- paste(form,'+ (1|',dnames[i],')')
  }
  
  form <-as.formula(form)
  tempfit <- lmer(form,data=dat,REML=REML,doFit=FALSE)
  # browser()
  if (num.fe){
    names(tempfit$fr$fixef) <- c(outer('fe.', 1:num.fe, paste, sep=''), '1', 'x', 'x*t')
  }else{
    names(tempfit$fr$fixef) <- c('1','x','x*t')
  }
  
  tn <- names(tempfit$FL$fl)
  tn <- tn[attr(tempfit$FL$fl,'assign')]
  ind <- 1:length(tn)
  sn <- dnames[2:2]
  snames <- 'FLM'
  for(i in 1:1){
    k <- ind[sn[i]==tn]
    tempfit$FL$trms[[k]]$A <- tempfit$FL$trms[[k]]$Zt <- 
      as(t(Design[,names(parinds)==sn[i]]),'dgCMatrix')
    attr(tempfit$FL$trms[[k]]$ST,'dimnames') <- list(snames[i],snames[i])
  }
  
  m0 <- do.call(lme4:::lmer_finalize,tempfit)
  # browser()
  if(!BOOT){
    var1 <- lme4:::VarCorr(m)[[1]]  #onet.Zx
    var2 <- lme4:::VarCorr(mA)[[1]]  #Zt.Zx
    var0 <- lme4:::VarCorr(m0)[[1]]
    if(var0==0)
      return(rep(NA,4))
    if(var1==0 || var2==0){
      return(c(NA,0,0,0))
    }else{
      test <- exactRLRT(m,mA,m0)
      pval <- test$p
      reject90 <- pval < .1 
      reject95 <- pval < .05
      reject99 <- pval < .01 
    }
    
    ret <- c(
      pval=pval, reject90=reject90, reject95=reject95, reject99=reject99
    )
    
    return(ret)
  }else{
   # browser()
    # from ?simulate-mer in lme4
    pboot <- function(m0,m1) {
      s <- simulate(m0)
      L0 <- logLik(refit(m0,s))
      L1 <- logLik(refit(m1,s))
      2*(L1-L0)
    }
    
    obsdev <- c(2*(logLik(mA)-logLik(m0)))
    pboot.res <- try(replicate(n.sims.boot, pboot(m0, mA)),FALSE)
    if (class(pboot.res)=='try-error'){
      return(rep(NA, 4))
    }else{
      p.val <- mean(obsdev < pboot.res)
      return(c(pval=p.val, reject90=p.val < .1, reject95=p.val < .05, reject99=p.val < .01))
    }
  }
}