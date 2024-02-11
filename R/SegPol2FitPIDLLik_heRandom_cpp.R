ggconst =1e30;
tttt <- environment();

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
GMOptimSegPolDF = function(par,fn, tol=1e-4,maxit=500)
{
  x2 = optim(par=par, fn, control=list(maxit=20000))
  try({x2 <- optim(par=x2$par, fn, method="BFGS")}, silent=T)
  x1 = x2$par;
  fx1 = x2$value
  new_combgrad = rep(NA,length(par))
  list(par=x1,fn=fx1,gr=new_combgrad, counts=c(nfunction=x2$counts[1],ngratient=NA),convergence=x2$convergence)
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
logit = function(x) .Call("logit",as.numeric(x))

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
expit <- function(x)
{
   res <- NULL
   res[x<0] <- exp(x[x<0])/(1+exp(x[x<0]))
   res[x>=0] <- exp(x[x>=0])/(1+exp(x[x>=0]))
   res
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
InciSegPol2PID <- function(nx, nEst_XPars, frestrknots)
{
  if((length(nx))==0)
  {
    return(NULL)
  } else
  {
    .Call("InciSegPol2PID_Rv",as.numeric(nx),as.numeric(nEst_XPars), as.numeric(frestrknots))
  }
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
ProjSegPol2PID <- function(XPars0,nVectX, grecov, uscmax,frestrknots,scfrac)
{
  if((length(XPars0))==0|(length(nVectX))==0|(length(grecov))==0|(length(uscmax))==0|(length(frestrknots))==0|(length(scfrac))==0)
  {
    stop("arguments should have non zero length")
  } else
  {
    if((length(XPars0)-3)%%2==0)
    {
      .Call("ProjSegPol2PID",as.numeric(XPars0),as.numeric(nVectX), as.numeric(grecov), as.numeric(uscmax),as.numeric(frestrknots), as.numeric(scfrac))
    } else
    {
      print(XPars0)
      stop("the length of XPars0 should be an odd number")
    }
  }
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
PrevRandEff <- function(prev,sigma)
{
  if((length(prev)==0)|(length(sigma)==0))
  {
    stop("arguments should have non zero length")
  } else
  {
    .Call("PrevRandEff",as.numeric(prev),as.numeric(sigma))
  }
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
nLBCheckCalc <- function(XPars, frestrknots)
{
  s_XPars=XPars;
  if(length(s_XPars)==1)
  {
    s_XPars = c(s_XPars,0,0)
  }
  if((length(s_XPars)-3)%%2==0)
  {
    .Call("LBCheckCalc",as.numeric(s_XPars),as.numeric(frestrknots))
  } else
  {
    stop("the length of XPars should be an odd number")
  }
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
nUBCheckCalc <- function(XPars, frestrknots, upb)
{
  s_XPars=XPars;
  if(length(s_XPars)==1)
  {
    s_XPars = c(s_XPars,0,0)
  }
  if((length(s_XPars)-3)%%2==0)
  {
    .Call("UBCheckCalc",as.numeric(s_XPars),as.numeric(frestrknots),as.numeric(upb))
  } else
  {
    stop("the length of XPars should be an odd number")
  }
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
LLikSegPol2PID_Random <- function(XPars, VectX,VectY,VectWeights, grecov, uscmax,frestrknots,scfrac)
{
  if(length(XPars)==2)
  {
    XPars = c(XPars[1],0,0,XPars[2])
  }

  XPars0 = XPars[-length(XPars)];
  g_sigma = exp(XPars[length(XPars)]);
  if(XPars[length(XPars)]>=10) g_sigma= exp(10)
  if(XPars[length(XPars)] <= -10) g_sigma= exp(-10)

  nVectX = c(0.001,VectX,1);
  nprojy = ProjSegPol2PID(XPars0,nVectX, grecov, uscmax,frestrknots,scfrac)

  nprojy[nprojy>10] = 10;
  nprojy = PrevRandEff(nprojy, g_sigma);

  projy = nprojy[-c(1,length(nprojy))]
  cprojy = 1-projy

  guyep=1e-10
  uu0 = projy<= guyep; uu1 = projy >= 1-guyep
  projy[uu0]= -expm1(-projy[uu0]^2)+1e-20
  projy[uu1]= 1+expm1(-(1-projy[uu1])^2)

  uu0 = cprojy<= guyep; uu1 = cprojy >= 1-guyep
  cprojy[uu0]= -expm1(-cprojy[uu0]^2)+1e-20
  cprojy[uu1]= 1+expm1(-(1-cprojy[uu1])^2)
  res = sum(VectWeights*(VectY*log(projy)+(1-VectY)*log(cprojy)));

  -res+tttt$ggconst*sum((max(nprojy[c(1,length(nprojy))]-0.2,0)^2))+100*(XPars[length(XPars)]< -10)*(XPars[length(XPars)]+10)^2 + 100*(XPars[length(XPars)]> 10)*(XPars[length(XPars)]-10)^2
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
SSSegPol2PID_Random <- function(XPars, VectX,VectY,VectWeights, grecov, uscmax,frestrknots,scfrac)
{
  if(length(XPars)==2)
  {
    XPars = c(XPars[1],0,0,XPars[2])
  }

  XPars0 = XPars[-length(XPars)];
  g_sigma = exp(XPars[-length(XPars)]);
  if(XPars[-length(XPars)]>=100) g_sigma= exp(100)
  nVectX = c(0.001,VectX,1);
  nprojy = ProjSegPol2PID(XPars0,nVectX, grecov, uscmax,frestrknots,scfrac)
  nprojy[nprojy>10] = 10;
  nprojy = PrevRandEff(nprojy, g_sigma);
  projy = nprojy[-c(1,length(nprojy))]

  cprojy = 1-projy
  res = sum(VectWeights*(VectY-projy)^2);
  res+tttt$ggconst*sum((max(nprojy[c(1,length(nprojy))]-0.2,0)^2))++100*(XPars[length(XPars)]< -10)*(XPars[length(XPars)]+10)^2 + 100*(XPars[length(XPars)]> 10)*(XPars[length(XPars)]-10)^2
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
GFitSegPol2PIDLLik_Random <- function (StartXPars, VectX, ObsY, VectWeights, incicheck = FALSE, recov0, uscmax, frestrknots,scfrac,upb=0.4)
{
  ScaledVectX = VectX;
  npar = length(StartXPars);
  old_XPars = StartXPars #runif(3) #

  #
  curvesmoothing <- function(x)
  {
    res <- 0
    res[(length(x)>=5)] <- 10000*(x[-c(1,2)][2]^2+10*x[length(x)]^2)#Extra constraint to avoid rebounds at the end or at the beginning of the trend
    res
  }

  fntomax = function(x)
  {
    x1 = x[-length(x)];
    res = LLikSegPol2PID_Random(x[-c(1,2)], ScaledVectX,ObsY,VectWeights,grecov=recov0, uscmax,frestrknots,scfrac)+
      tttt$ggconst*(x[1]^2+exp(x[1]))*nLBCheckCalc(x1[-c(1,2)],frestrknots)+
      tttt$ggconst*(x[2]^2+exp(x[2]))*nUBCheckCalc(x1[-c(1,2)],frestrknots,upb)+curvesmoothing(x1);
    res
  }
  condcheck = 1

  x0 = c(5,5,old_XPars)
  res = GMOptimSegPolDF(par=x0,fn=fntomax, tol=1e-8,maxit=15000) #
  x1 = res$par

  new_lossf = LLikSegPol2PID_Random(x1[-c(1,2)], ScaledVectX,ObsY,VectWeights,grecov=recov0, uscmax,frestrknots,scfrac)

  AIC = 2*new_lossf+2*npar
  list(par=x1,value=new_lossf, uscmax=uscmax, AIC=AIC, convergence=res$convergence)
}


#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
FitSegPol2PIDLLik_Random <- function (StartXPars, VectX, ObsY, VectWeights, incicheck = FALSE, recov0, uscmax, frestrknots,scfrac)
{
  ScaledVectX = VectX;
  npar = length(StartXPars);
  old_XPars = StartXPars #runif(3) #
  upb=0.4
  fn_gconst=1e30
  #Guyz
  curvesmoothing <- function(x)
  {
    res <- 0
    res[(length(x)>=5)] <- 10000*(x[-c(1,2)][2]^2+10*x[length(x)]^2)#Extra constraint to avoid rebounds at the end or at the beginning of the trend
    res
  }

  fntomax = function(x)
  {
    x1 = x[-length(x)];
    res = LLikSegPol2PID_Random(x[-c(1,2)], ScaledVectX,ObsY,VectWeights,grecov=recov0, uscmax,frestrknots,scfrac)+
      fn_gconst*(x[1]^2+exp(x[1]))*nLBCheckCalc(x1[-c(1,2)],frestrknots)+
      fn_gconst*(x[2]^2+exp(x[2]))*nUBCheckCalc(x1[-c(1,2)],frestrknots,upb)+curvesmoothing(x1)
    res
  }

  lstomax = function(x)
  {
    x1 = x[-length(x)];
    SSSegPol2PID_Random(x[-c(1,2)], ScaledVectX,ObsY,VectWeights,grecov=recov0, uscmax,frestrknots,scfrac)+
      fn_gconst*(x[1]^2+exp(x[1]))*nLBCheckCalc(x1[-c(1,2)],frestrknots)+
      fn_gconst*(x[2]^2+exp(x[2]))*nUBCheckCalc(x1[-c(1,2)],frestrknots,upb)+curvesmoothing(x1)
  }

  condcheck = 1
  olop0 = list(par=c(5,5,old_XPars), convergence=1)

  olop = GFitSegPol2PIDLLik_Random(StartXPars, VectX, ObsY, VectWeights, incicheck = FALSE, recov0, uscmax, frestrknots,scfrac);

  next_XPars = olop$par[-c(1:2)]
  condcheck = nLBCheckCalc(next_XPars[-length(next_XPars)],frestrknots)+nUBCheckCalc(next_XPars[-length(next_XPars)],frestrknots,upb)
  kk0= 0;
  cc = 1;

  fn_gconst = 1e30
  while(TRUE)
  {
    olop = optim(par=olop0$par, fn=fntomax, control=list(maxit=20000))
    try({olop = optim(par=olop$par, fn=fntomax, method="BFGS")}, silent = T)

    next_XPars = olop$par[-c(1:2)]
    condcheck = nLBCheckCalc(next_XPars[-length(next_XPars)],frestrknots)+nUBCheckCalc(next_XPars[-length(next_XPars)],frestrknots,upb)
    if((condcheck<=0) | (cc>10)) break;
    fn_gconst = max(1/condcheck,fn_gconst)*100
    cc = cc+1
  }

  olop$value = LLikSegPol2PID_Random(next_XPars, ScaledVectX,ObsY,VectWeights,grecov=recov0, uscmax,frestrknots,scfrac)
  AIC = 2*olop$value+2*npar
  list(par=next_XPars,value=olop$value, uscmax=uscmax, AIC=AIC)
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
BestFitPIDLLik_Random <- function(vX,vY,vW, maxknots=3,max_FProj=5,ficheck=FALSE,frecov)
{
  bestlres= lapply(1:(maxknots+1), function(ii) NA)
  gxmax = max_FProj;
  scvX = vX/gxmax;
  uscmax_in = gxmax+1;
  restrknots_in = (sort(vX)[length(vX)-1])/gxmax;
  if(length(unique(vX))>2)
  {
    suuX <- sort(unique(vX))
    if(length(suuX)==3) restrknots_in = sum(sort(suuX)[length(suuX)-c(1:2)])/2/gxmax else restrknots_in = sum(sort(suuX)[length(suuX)-c(2:3)])/2/gxmax;
  }
  scfrac_in = gxmax;
  vect_aic = numeric()
  ll = -1;
  while(TRUE)
  {
    XPars0 = c(rep(0,3+2*ll),0)
    XPars0[1] = logit(0.05)

    interm = FitSegPol2PIDLLik_Random(StartXPars=XPars0, VectX=scvX, ObsY=vY, VectWeights=vW, recov0=frecov, uscmax=gxmax+1,frestrknots=restrknots_in, scfrac=scfrac_in)
    XPars0 = interm$par

    XPars1 = c(rep(0,3+2*ll),0)
    XPars1[1] = logit(0.1)
    interm1 = FitSegPol2PIDLLik_Random(StartXPars=XPars1, VectX=scvX, ObsY=vY, VectWeights=vW, recov0=frecov, uscmax=gxmax+1,frestrknots=restrknots_in, scfrac=scfrac_in)

    if(interm1$AIC < interm$AIC)
    {
      interm = interm1;
      XPars1 = c(rep(0,3+2*ll),0)
      XPars1[1] = logit(0.15)
      interm1 = FitSegPol2PIDLLik_Random(StartXPars=XPars1, VectX=scvX, ObsY=vY, VectWeights=vW, recov0=frecov, uscmax=gxmax+1,frestrknots=restrknots_in, scfrac=scfrac_in)
      if(interm1$AIC < interm$AIC)
      {
        interm = interm1
      }
    }

    XPars2 = c(rep(0,3+2*ll),0)
    XPars2[1] = logit(0.01)
    interm2 = FitSegPol2PIDLLik_Random(StartXPars=XPars2, VectX=scvX, ObsY=vY, VectWeights=vW, recov0=frecov, uscmax=gxmax+1,frestrknots=restrknots_in, scfrac=scfrac_in)
    if(interm2$AIC < interm$AIC)
    {
      interm = interm2
    }

    XPars3 = c(rep(0,3+2*ll),0)
    XPars3[1] = logit(0.001)
    interm3 = FitSegPol2PIDLLik_Random(StartXPars=XPars2, VectX=scvX, ObsY=vY, VectWeights=vW, recov0=frecov, uscmax=gxmax,frestrknots=restrknots_in, scfrac=scfrac_in)
    if(interm3$AIC < interm$AIC)
    {
      interm = interm3;
    }

    vect_aic[ll+2] = interm$AIC
    bestlres[[ll+2]] = interm;

    ll=ll+1;
    if(ll>=maxknots)
    {
      break
    }
  }

  ind = which.min(vect_aic)[1]
  Est_XPars = bestlres[[ind]]$par
  funcestprojprev = function(x)
  {
    nx = x/gxmax
    nEst_XPars = Est_XPars[-length(Est_XPars)];
    if(length(nEst_XPars)==1)
    {
      nEst_XPars = c(nEst_XPars,0,0)
    }
    yproj = ProjSegPol2PID(nEst_XPars, nx,grecov=frecov, uscmax=gxmax+1, frestrknots=restrknots_in, scfrac=scfrac_in)
    yproj = PrevRandEff(yproj,exp(Est_XPars[length(Est_XPars)]))
    yproj
  }
  funcestprojinci = function(x)
  {
    nx = x/gxmax
    nEst_XPars = Est_XPars[-length(Est_XPars)];
    if(length(nEst_XPars)==1)
    {
      nEst_XPars = c(nEst_XPars,0,0)
    }
    iproj = InciSegPol2PID(nx, nEst_XPars, frestrknots=restrknots_in)
    iproj
  }
  list(All_Pars = bestlres[[ind]],AIC=vect_aic, best_AIC=vect_aic[ind],recov=frecov,gxmax=gxmax,funcestprojprev=funcestprojprev, scfrac=scfrac_in,funcestprojinci=funcestprojinci,frestrknots=restrknots_in,data=list(vX=vX,vY=vY,vW=vW))
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
BestFit2PIDLLik_Random <- function(vX=obs_time,vY=obs_prev,vW=obs_weig, xResBestFit)
{
  gxmax =  xResBestFit$gxmax;
  scvX = vX /gxmax;
  XPars0 = xResBestFit$All_Pars$par
  ficheck = xResBestFit$ficheck
  frecov = xResBestFit$recov
  restrknots_in = xResBestFit$frestrknots
  scfrac_in =xResBestFit$scfrac

  XPars1 = XPars0
  interm = FitSegPol2PIDLLik_Random(StartXPars=XPars1, VectX=scvX, ObsY=vY, VectWeights=vW, recov0=frecov,uscmax=gxmax+1,frestrknots=restrknots_in, scfrac=scfrac_in)

  Est_XPars = interm$par

  funcestprojprev <- function(x)
  {
    nx = x /gxmax
    nEst_XPars = Est_XPars[-length(Est_XPars)];
    if(length(nEst_XPars)==1)
    {
      nEst_XPars = c(nEst_XPars,0,0)
    }
    yproj = ProjSegPol2PID(nEst_XPars, nx, grecov=frecov, uscmax=gxmax+1, frestrknots=restrknots_in, scfrac=scfrac_in)
    yproj = PrevRandEff(yproj,exp(Est_XPars[length(Est_XPars)]))
    yproj
  }

  funcestprojinci <- function(x)
  {
    nx = x/gxmax
    nEst_XPars = Est_XPars[-length(Est_XPars)];
    if(length(nEst_XPars)==1)
    {
      nEst_XPars = c(nEst_XPars,0,0)
    }
    iproj = InciSegPol2PID(nx, nEst_XPars,frestrknots=restrknots_in)
    iproj
  }
  list(All_Pars = interm, recov=frecov,gxmax=gxmax,funcestprojprev=funcestprojprev, scfrac=scfrac_in, funcestprojinci=funcestprojinci,frestrknots=restrknots_in, data=list(vX=vX,vY=vY,vW=vW))
}

# GetLROR <- function(dat, L_surveytypes)
# {
#   fn <- function(par)
#   {
#     X <- exp(par)
#     projprev <- numeric()
#     projprev[dat$Data_type%in%L_surveytypes$PregWom] <- dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]*X[1]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]+dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]*X[1])
#     projprev[dat$Data_type%in%L_surveytypes$GeneWom] <- dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]*X[2]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]+dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]*X[2])
#     projprev[dat$Data_type%in%L_surveytypes$GeneMen] <- dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]*X[3]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]+dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]*X[3])
#
#     nnloss <-  -sum(dat$N_tested*(dat$p*log(pmax(projprev,1e-10)) + (1-dat$p)*log(1-projprev)), na.rm=T) + sum(par^2)/100
#     nnloss
#   }
#
#   par0 <- rep(0,3)
#   opt0 <- optim(par0,fn)
#   opt1 <- optim(opt0$par,fn, method="BFGS")
#
#   res <- exp(opt1$par)
#   names(res) <- c("PregWom","GeneWom","GeneMen")
#   res
# }

GetLROR <- function(dat, L_surveytypes)
{
  fn <- function(par)
  {
    X <- exp(par)
    projprev <- numeric()
    projprev[dat$Data_type%in%L_surveytypes$PregWom] <- dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]*X[1]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]+dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]*X[1])
    projprev[dat$Data_type%in%L_surveytypes$GeneWom] <- dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]*X[2]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]+dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]*X[2])
    projprev[dat$Data_type%in%L_surveytypes$GeneMen] <- dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]*X[3]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]+dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]*X[3])
    projprev[dat$Data_type%in%L_surveytypes$BloodDo] <- dat$estimprev[dat$Data_type%in%L_surveytypes$BloodDo]*X[3]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$BloodDo]+dat$estimprev[dat$Data_type%in%L_surveytypes$BloodDo]*X[3])

    nnloss <-  -sum(dat$N_tested*(dat$p*log(pmax(projprev,1e-10)) + (1-dat$p)*log(1-projprev)), na.rm=T) + sum(par^2)/1000
    nnloss
  }

  #Pregnant Women
  fm_pw <- function(par)
  {
    X <- exp(par)
    projprev <- numeric()
    projprev[dat$Data_type%in%L_surveytypes$PregWom] <- dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]*X[1]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]+dat$estimprev[dat$Data_type%in%L_surveytypes$PregWom]*X[1])
    nnloss <-  -sum(dat$N_tested*(dat$p*log(pmax(projprev,1e-10)) + (1-dat$p)*log(1-projprev)), na.rm=T) + sum(par^2)/100
    nnloss
  }

  par0_pw <- rep(0,2)
  opt0_pw <- optim(par0_pw,fm_pw)
  opt1_pw <- optim(opt0_pw$par,fm_pw, method="BFGS")$par

  #LR Women
  fm_gw <- function(par)
  {
    X <- exp(par)
    projprev <- numeric()
    projprev[dat$Data_type%in%L_surveytypes$GeneWom] <- dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]*X[1]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]+dat$estimprev[dat$Data_type%in%L_surveytypes$GeneWom]*X[1])
    nnloss <-  -sum(dat$N_tested*(dat$p*log(pmax(projprev,1e-10)) + (1-dat$p)*log(1-projprev)), na.rm=T) + sum(par^2)/100
    nnloss
  }

  par0_gw <- rep(0,2)
  opt0_gw <- optim(par0_gw,fm_gw)
  opt1_gw <- optim(opt0_gw$par,fm_gw, method="BFGS")$par

  #LR Men
  fm_gm <- function(par)
  {
    X <- exp(par)
    projprev <- numeric()
    projprev[dat$Data_type%in%L_surveytypes$GeneMen] <- dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]*X[1]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]+dat$estimprev[dat$Data_type%in%L_surveytypes$GeneMen]*X[1])
    nnloss <-  -sum(dat$N_tested*(dat$p*log(pmax(projprev,1e-10)) + (1-dat$p)*log(1-projprev)), na.rm=T) + sum(par^2)/100
    nnloss
  }

  par0_gm <- rep(0,2)
  opt0_gm <- optim(par0_gm,fm_gm)
  opt1_gm <- optim(opt0_gm$par,fm_gm, method="BFGS")$par

  #Blood Donors
  fm_bd <- function(par)
  {
    X <- exp(par)
    projprev <- numeric()
    projprev[dat$Data_type%in%L_surveytypes$BloodDo] <- dat$estimprev[dat$Data_type%in%L_surveytypes$BloodDo]*X[1]/(1-dat$estimprev[dat$Data_type%in%L_surveytypes$BloodDo]+dat$estimprev[dat$Data_type%in%L_surveytypes$BloodDo]*X[1])
    nnloss <-  -sum(dat$N_tested*(dat$p*log(pmax(projprev,1e-10)) + (1-dat$p)*log(1-projprev)), na.rm=T) + sum(par^2)/100
    nnloss
  }

  par0_bd <- rep(0,2)
  opt0_bd <- optim(par0_bd,fm_bd)
  opt1_bd <- optim(opt0_bd$par,fm_bd, method="BFGS")$par

  par1 <- c(par0_pw[1],opt1_gw[1], opt1_gm[1], opt1_bd[1])
  par1 <- optim(par1,fn)$par

  res <- exp(par1)
  names(res) <- c("PregWom","GeneWom","GeneMen")
  res
}

var_names_for_long <- c("Prevalence (%)", "Prevalence (%)","Prevalence (%)",
                        "Prevalence cases (#)", "Prevalence cases (#)", "Prevalence cases (#)",
                        "Incidence rate", "Incidence rate","Incidence rate",
                        "Incidence cases", "Incidence cases", "Incidence cases")

var_names_sexes <- c("Males", "Females","Both sexes",
                     "Males", "Females","Both sexes",
                     "Males", "Females","Both sexes",
                     "Males", "Females","Both sexes")


idxtoplot <- list(c("EstimatePrevF", "PrevLB_2.5%F", "PrevUB_97.5%F"),
                  c("EstimatePrevM", "PrevLB_2.5%M", "PrevUB_97.5%M"),
                  c("EstimatePrevM+F", "PrevLB_2.5%M+F", "PrevUB_97.5%M+F"),
                  c("CasePrev_EstF", "CasePrev_LB_2.5%F",  "CasePrev_UB_97.5%F"),
                  c("CasePrev_EstM", "CasePrev_LB_2.5%M",  "CasePrev_UB_97.5%M"),
                  c("CasePrev_EstM+F", "CasePrev_LB_2.5%M+F",  "CasePrev_UB_97.5%M+F"),
                  c("EstimateIncF", "IncLB_2.5%F", "IncUB_97.5%F"),
                  c("EstimateIncM", "IncLB_2.5%M", "IncUB_97.5%M"),
                  c("EstimateIncM+F", "IncLB_2.5%M+F", "IncUB_97.5%M+F"),
                  c("CaseInc_EstF", "CaseInc_LB_2.5%F", "CaseInc_UB_97.5%F"),
                  c("CaseInc_EstM", "CaseInc_LB_2.5%M", "CaseInc_UB_97.5%M"),
                  c("CaseInc_EstM+F", "CaseInc_LB_2.5%M+F", "CaseInc_UB_97.5%M+F")

)


Default_List_surveytypes <- vector("list",3)
names(Default_List_surveytypes) <- c("PregWom", "GeneWom", "GeneMen")
Default_List_surveytypes$PregWom <- c("ANC Routine screening","ANC Survey")
Default_List_surveytypes$GeneWom <- c("BloodDonor Screening Women","BloodDonor Screening Men + Women",
                              "Survey LowRisk Men+Women", "Survey LowRisk Men + Women","Survey LowRisk Women",
                              "ANC Routine screening","ANC Survey")
Default_List_surveytypes$GeneMen <- c("BloodDonor Screening Men","BloodDonor Screening Men + Women",
                              "Survey LowRisk Men+Women", "Survey LowRisk Men + Women", "Survey LowRisk Men")

Default_List_surveytypes$BloodDo <- c("BloodDonor Screening Men","BloodDonor Screening Men + Women",
                                      "BloodDonor Screening Women")

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
fCountryAnalysis_glob <- function(Nboots=1000,
                                  fname.data.file = name.data.file,
                                  Fmaxknots,
                                  FB_ProjMax_F,
                                  FB_ProjMax_B,
                                duration.syphilis,
                                data.syphilis,
                                PopSizeBoth,
                                PopSizeMen,
                                PopSizeWomen,
                                PopSizeMSM,
                                PopSizeFSW,
                                ISO3,
                                year_predict,
                                zerprev_adj=1/100,
                                MtoFRatio,
                                LRtoHRPOR,
                                fn_min_year_last_data,
                                fn_cores=4,
                                autosavefile=TRUE,
                                fn_filter_survey_LR = Default_List_surveytypes)
{
  PopSize = PopSizeBoth;
  PopSizeMales = PopSizeMen;
  PopSizeFemales = PopSizeWomen
  PopSizeMSM_in = PopSizeMSM;
  PopSizeFSW_in = PopSizeFSW;
  SyphData = data.syphilis;
  #Name of the output file
  name.out.file = fname.data.file
  if(grepl("/",name.out.file))
  {
    nn <- max(unlist(gregexpr('/', name.out.file)))
    name.out.file <- substr(name.out.file,nn+1, nchar(name.out.file))
  }
  name.out.file <- paste(paste(substr(name.out.file,1,nchar(name.out.file)-5),as.character(Sys.Date()),sep="_"),"_out.xlsx",sep="")
  if(nchar(name.out.file)>218)
  {
    name.out.file <- substr(name.out.file,nchar(name.out.file)-217,nchar(name.out.file))
  }

  #List of countries
  lcountry = unique(data.syphilis$Country) #lcountry = levels(data.syphilis$Country)

  res.syphilis = data.frame();
  #Loading the excel file where the results will be stored
  wb = openxlsx::loadWorkbook(fname.data.file);
  nnamessheetwb =  openxlsx::sheets(wb)#names(getSheets(wb));
  if(!is.element("SYPH_RBootstrap_All",nnamessheetwb))
  {
    Syphilis_Rbootstrap = "SYPH_RBootstrap_All"
    openxlsx::addWorksheet(wb,"SYPH_RBootstrap_All")#createSheet(wb,"SYPH_RBootstrap_All")#createSheet(wb,"SYPH_RBootstrap_All") # I am going to write out Syphilis results for all the countries here
  } else
  {
    Syphilis_Rbootstrap = "SYPH_RBootstrap_All"
  }

  if(!is.element("SYPH_YearCheck_Glob",nnamessheetwb))
  {
    Syphilis_YearCheck_Glob = "SYPH_YearCheck_Glob"
    openxlsx::addWorksheet(wb,"SYPH_YearCheck_Glob")
  } else
  {
    Syphilis_YearCheck_Glob = "SYPH_YearCheck_Glob"
  }

  if(!is.element("SYPH_AdjTestValCheck_Glob",nnamessheetwb))
  {
    Syphilis_AdjTestCheck_Glob = "SYPH_AdjTestValCheck_Glob"
    openxlsx::addWorksheet(wb,Syphilis_AdjTestCheck_Glob)
  } else
  {
    Syphilis_AdjTestCheck_Glob = "SYPH_AdjTestValCheck_Glob"
  }

  # Workbook Styles
  #New names
  ttt = c("Year", "Country","ISO3","WHO_Region","EstimatePrevF", "MedianPrevF", "PrevLB_2.5%F", "PrevUB_97.5%F",
          "EstimatePrevM", "MedianPrevM", "PrevLB_2.5%M", "PrevUB_97.5%M",
          "EstimatePrevM+F", "MedianPrevM+F", "PrevLB_2.5%M+F", "PrevUB_97.5%M+F")

  tttEl = c("CasePrev_EstF", "CasePrev_MedF", "CasePrev_LB_2.5%F", "CasePrev_UB_97.5%F",
            "CasePrev_EstM", "CasePrev_MedM", "CasePrev_LB_2.5%M", "CasePrev_UB_97.5%M",
            "CasePrev_EstM+F", "CasePrev_MedM+F", "CasePrev_LB_2.5%M+F", "CasePrev_UB_97.5%M+F")

  tttInc = c("EstimateIncF", "MedianIncF", "IncLB_2.5%F", "IncUB_97.5%F",
             "EstimateIncM", "MedianIncM", "IncLB_2.5%M", "IncUB_97.5%M",
             "EstimateIncM+F", "MedianIncM+F", "IncLB_2.5%M+F", "IncUB_97.5%M+F",
             "EstimateIncPregWom", "MedianIncPregWom", "IncLB_2.5%PregWom", "IncUB_97.5%PregWom")

  tttIncEl = c("CaseInc_EstF", "CaseInc_MedF", "CaseInc_LB_2.5%F", "CaseInc_UB_97.5%F",
               "CaseInc_EstM", "CaseInc_MedM", "CaseInc_LB_2.5%M", "CaseInc_UB_97.5%M",
               "CaseInc_EstM+F", "CaseInc_MedM+F", "CaseInc_LB_2.5%M+F", "CaseInc_UB_97.5%M+F",
               "NationalPop15-49yF","NationalPop15-49yM","NationalPop15-49yM+F", "Country_Curve_Fit","Date_Last_Run")

  syph_nn1 = c("year","Country","ISO3","WHO_Region","PrevEstF","PrevEstF_Med", "PrevEstF_LoB", "PrevEstF_UpB",
               "PrevEstM","PrevEstM_Med", "PrevEstM_LoB", "PrevEstM_UpB",
               "PrevEstMPlusF","PrevEstMPlusF_Med", "PrevEstMPlusF_LoB", "PrevEstMPlusF_UpB",
               "CasePrevEstF", "CasePrevMedF", "CasePrevF_LoB", "CasePrevF_UpB",
               "CasePrevEstM", "CasePrevMedM", "CasePrevM_LoB", "CasePrevM_UpB",
               "CasePrevEstMPlusF", "CasePrevMedMPlusF", "CasePrevMPlusF_LoB", "CasePrevMPlusF_UpB")

  syph_nninc = c("InciEstF","InciEstF_Med", "InciEstF_LoB", "InciEstF_UpB",
                 "InciEstM","InciEstM_Med", "InciEstM_LoB", "InciEstM_UpB",
                 "InciEstMPlusF","InciEstMPlusF_Med", "InciEstMPlusF_LoB", "InciEstMPlusF_UpB",
                 "CaseInciEstF", "CaseIncMedF", "CaseIncF_LoB", "CaseIncF_UpB",
                 "CaseInciEstM", "CaseIncMedM", "CaseIncM_LoB", "CaseIncM_UpB",
                 "CaseInciEstMPlusF", "CaseIncMedMPlusF", "CaseIncMPlusF_LoB", "CaseIncMPlusF_UpB",
                 "PopM","PopF","PopMPlusF")


  tttKPs = c("Year", "Country","ISO3","WHO_Region","EstimatePrevFSW", "MedianPrevFSW", "PrevLB_2.5%FSW", "PrevUB_97.5%FSW",
          "EstimatePrevMSM", "MedianPrevMSM", "PrevLB_2.5%MSM", "PrevUB_97.5%MSM",
          "EstimatePrevPregWom", "MedianPrevPregWom", "PrevLB_2.5%PregWom", "PrevUB_97.5%PregWom")

  tttElKPs = c("CasePrev_EstFSW", "CasePrev_MedFSW", "CasePrev_LB_2.5%FSW", "CasePrev_UB_97.5%FSW",
            "CasePrev_EstMSM", "CasePrev_MedMSM", "CasePrev_LB_2.5%MSM", "CasePrev_UB_97.5%MSM")

  syph_nn1KPs = c("year","Country","ISO3","WHO_Region","PrevEstFSW","PrevEstFSW_Med", "PrevEstFSW_LoB", "PrevEstFSW_UpB",
               "PrevEstMSM","PrevEstMSM_Med", "PrevEstMSM_LoB", "PrevEstMSM_UpB",
               "CasePrevEstFSW", "CasePrevMedFSW", "CasePrevFSW_LoB", "CasePrevFSW_UpB",
               "CasePrevEstMSM", "CasePrevMedMSM", "CasePrevMSM_LoB", "CasePrevMSM_UpB")

  nsep <- 0
  nligne_decal_titre <- 0
  nrownum <- 0

  est_all <- list()
  results_all <- list()
  countries_in <- list()

  PopSize_all <- list()
  PopSize_all2 <-  list()

  ii <- 0
  nncc <- 0
  nncc2 <- 0

  gfy0 <- min(year_predict)
  gly0 <- max(year_predict)+1
  gslen0 <- 30
  data.syphilis$popss <- NA
  WorldData <- data.syphilis[-(1:nrow(data.syphilis)),]

  ewho_regions <- unique(data.syphilis$WHO_region)
  WHO_regions_Data <-  lapply(1:length(ewho_regions), function(whor) WorldData)
  names(WHO_regions_Data) <- ewho_regions

  PopSize_all_WHO_regions <- lapply(1:length(ewho_regions), function(whor) list())
  PopSize_all2_WHO_regions <- lapply(1:length(ewho_regions), function(whor) list())

  names(PopSize_all_WHO_regions) <- ewho_regions
  names(PopSize_all2_WHO_regions) <- ewho_regions

  est_all_WHO_regions <- lapply(1:length(ewho_regions), function(whor) list())
  results_all_WHO_regions <- lapply(1:length(ewho_regions), function(whor) list())
  names(est_all_WHO_regions) <- ewho_regions
  names(results_all_WHO_regions) <- ewho_regions

  nncc_WHO_region <- rep(0,length(ewho_regions))
  nncc2_WHO_region <- rep(0,length(ewho_regions))

  DurSyph_all_WHO_regions <- rep(0,length(ewho_regions)) #For weighthed syphilis durations
  names(DurSyph_all_WHO_regions) <- ewho_regions

  SampSize_all_WHO_regions <- rep(0,length(ewho_regions)) #For weighthed syphilis durations
  names(SampSize_all_WHO_regions) <- ewho_regions

  All_CountryDataUse=data.frame(CountryName=factor(levels(lcountry)),Year=numeric(),WithinRange=factor(levels=c("Yes","No")),Flatlined=factor(levels=c("Yes","No")))
  ee0 <- environment()

  library(parallel)
  require(doParallel)
  NumberOfCluster <- fn_cores;
  cluster <- makeCluster(NumberOfCluster)
  registerDoParallel(cluster)

  all_res <- foreach(lc = lcountry)%dopar%{
    #all_res <- lapply(lcountry, function(lc){#foreach(lc = lcountry)%dopar%{
    require(SpectrumSyphilis)
    CountryModel <- "PID+Spline"

    cdata.syph = data.syphilis[data.syphilis$Country==lc & !is.na(data.syphilis$Country) & !is.na(data.syphilis$Prevalence),]
    #print(lc)
    nww <- cdata.syph$WghtNSpectrum
    cdata.syph <- cdata.syph[!is.na(nww),]
    cdata.syph$Weights <- nww[!is.na(nww)]

    cdata.syph <- cdata.syph[!is.na(cdata.syph$Year),]
    cdata.syph <- cdata.syph[cdata.syph$Weights>0,]
    # imputation of zero prevalence

    cdata.syph$Prevalence[cdata.syph$Prevalence==0] <- 1/cdata.syph$N_tested[cdata.syph$Prevalence==0]
    ##
    # Helps for sensitivity analysis
    cdata.syph$Prevalence[cdata.syph$Prevalence==0] <- cdata.syph$Prevalence[cdata.syph$Prevalence==0]*zerprev_adj

    ###***###
    check_maxyear = FALSE
    if(nrow(cdata.syph)>=1)
    {
      myy = max(cdata.syph$Year,na.rm=T)
      if(length(myy)==1)
      {
        if(myy>=fn_min_year_last_data)#Lastest data point should be later than fn_min_year_last_data
        {
          ww = cdata.syph$Weights[which(cdata.syph$Year==myy)]
          if(length(ww[ww>0])>=1) check_maxyear=TRUE # At least one survey conducted after 2011 with a weight >0
        }
      }
    }

    ###***###
    if(check_maxyear)
    {
      numeric_iso3 =  as.character(cdata.syph$ISO3[1]);
      alpha_iso3 =  ISO3$Alpha_ISO3[ISO3$Numeric_ISO3==numeric_iso3 & !is.na(ISO3$Numeric_ISO3)];
      alpha_iso3 <- alpha_iso3[!is.na(alpha_iso3)][1]

      vv = as.character(alpha_iso3)

      if(!is.na(vv))
      {
        cdata.syph$popss = sapply(round(cdata.syph$Year), function(xx){
          res=NA
          vxx = PopSize[PopSize$ISO3==vv & !is.na(PopSize$ISO3),-1]
          res0 = unlist(vxx[PopSize[1,-1]==xx])
          if(length(res0)) res = res0;
          res
        })
      } else
      {
        cdata.syph$popss = rep(NA,length(cdata.syph$Year))
      }

      cdata.syph = cdata.syph[cdata.syph$Year>=gfy0,]
      rr_whoreg = cdata.syph$WHO_region[1]
      amy = which(names(PopSize_all2_WHO_regions)==rr_whoreg)
      treat.zone = cdata.syph$WHO_region[1]
      if(treat.zone=="EURO" | treat.zone=="AMRO"| treat.zone=="AMR"| treat.zone=="EUR") treat.zone =  "A" else if (treat.zone=="AFRO" | treat.zone=="AFR") treat.zone = "C" else treat.zone = "B"
      dur_Syph = duration.syphilis$duration[which(duration.syphilis$t.zone==treat.zone)];
      tbadd = dur_Syph*PopSize[PopSize$ISO3==vv & !is.na(PopSize$ISO3),-1][1];
      nntest = ifelse(is.vector(tbadd),length(tbadd),nrow(tbadd))
    }

    if(!check_maxyear)
    {
      return(0)
    } else
    {
      numeric_iso3 =  cdata.syph$ISO3[1];
      alpha_iso3 =  ISO3$Alpha_ISO3[ISO3$Numeric_ISO3==as.character(numeric_iso3) & !is.na(ISO3$Numeric_ISO3)];
      alpha_iso3 <- alpha_iso3[!is.na(alpha_iso3)][1]

      resboot_Syphilis = data.frame(matrix(NA,ncol=length(year_predict),nrow=Nboots))
      names(resboot_Syphilis)=paste("year",year_predict,sep="_")
      resboot_SyphilisPrevF = resboot_Syphilis
      resboot_SyphilisPrevM = resboot_Syphilis
      resboot_SyphilisInciF = resboot_Syphilis
      resboot_SyphilisInciM = resboot_Syphilis

      resboot_SyphilisPrevFSW = resboot_Syphilis
      resboot_SyphilisPrevMSM = resboot_Syphilis
      resboot_SyphilisPrevPregWom = resboot_Syphilis
      resboot_SyphilisPrevGeneWom = resboot_Syphilis
      resboot_SyphilisPrevGeneMen = resboot_Syphilis

      #Blood donors
      resboot_SyphilisPrevBloodDo = resboot_Syphilis

      Out_syphilis = data.frame(year=year_predict)
      Out_syphilis$Country = lc;
      Out_syphilis$ISO3 = alpha_iso3

      Out_syphilis$WHO_Region = cdata.syph$WHO_region[1]
      Out_syphilisKPs <- Out_syphilis
      Out_syphilisPregWom <- Out_syphilis
      Out_syphilisGeneWom <- Out_syphilis
      Out_syphilisGeneMen <- Out_syphilis
      Out_syphilisBloodDo <- Out_syphilis
      #Incidence
      incboot_Syphilis = resboot_Syphilis

      #Syphilis Analysis
      tab = cdata.syph
      tab$Npos = tab$N_tested*tab$Prevalence/100

      treat.zone = tab$WHO_region[1]
      if(treat.zone=="EURO" | treat.zone=="AMRO" | treat.zone=="EUR" | treat.zone=="AMR") treat.zone =  "A" else if (treat.zone=="AFRO" | treat.zone=="AFR") treat.zone = "C" else treat.zone = "B"
      dur_Syph = duration.syphilis$duration[which(duration.syphilis$t.zone==treat.zone)]

      ee=environment()
      tabtab=tab; ee$tab = tab;
      lll =  sapply(1:length(tab$Weights), function(iii) {
        if(!is.na(tabtab$Weights[iii]))
          if(tabtab$Weights[iii]==0){
            ee$tab$Weights[iii]=mean(tabtab$Weights[tabtab$Weights!=0 & tabtab$Survey== tabtab$Survey[iii]],na.rm=T)
          }
      })

      rm(tabtab); rm(lll)
      tab = ee$tab
      tab=tab[!is.na(tab$Weights) & !is.na(tab$N_tested),]

      CountryDataUse=data.frame(CountryName=rep(lc,length(year_predict)),Year=year_predict,WithinRange=rep("No",length(year_predict)),
                                Flatlined=rep("No",length(year_predict)))
      levels(CountryDataUse$WithinRange)=c("No","Yes")
      levels(CountryDataUse$Flatlined)=c("No","Yes")
      CountryDataUse$WithinRange[CountryDataUse$Year>=min(tab$Year) & CountryDataUse$Year<=max(tab$Year)] ="Yes"
      CountryDataUse$Flatlined[CountryDataUse$Year<=(min(tab$Year)-FB_ProjMax_B) | CountryDataUse$Year>=(max(tab$Year)+FB_ProjMax_F)] ="Yes"

      vvprob = tab$Weights
      err_syph = lapply(1:nrow(tab),function(ll) which(tab$DX_Code==tab$DX_Code[ll] & tab$Data_type==tab$Data_type[ll]))
      wei_syph = lapply(1:nrow(tab),function(ll) {ww0 = tab$Weights[which(tab$DX_Code==tab$DX_Code[ll] & tab$Data_type==tab$Data_type[ll])]; if(length(ww0)>0) ww0=ww0/sum(ww0); ww0} )

      err_samp = function(ll) {vv=as.vector(err_syph[[ll]]); sprob= as.vector(wei_syph[[ll]]) ; if (length(vv)==1) res= vv else res=sample(vv, prob=sprob, replace=TRUE)[1]}
      ERR_SAMP = function(ll) as.vector(sapply(ll,err_samp))

      tab$Weights = tab$Weights/sum(tab$Weights) # Normalize the weights so that the sample size is not affected
      tab$Weights = tab$Weights*sum(tab$N_tested,na.rm=T)

      tab$variance = tab$Npos/tab$N_tested*(1-tab$Npos/tab$N_tested)/tab$N_tested;
      tab$p = tab$Npos/tab$N_tested
      tab$alpha = tab$p*(tab$p*(1-tab$p)/tab$variance-1)
      tab$beta = (1-tab$p)*(tab$p*(1-tab$p)/tab$variance-1)

      amin = floor(min(tab$Year))-(FB_ProjMax_B+0.5)
      bmax = ceiling(max(tab$Year)) + (FB_ProjMax_F+0.5)
      bmoinsa = bmax-amin

      aa <- 0.5
      bb <- bmoinsa-0.5

      Vect_Year_Pred = (Out_syphilis$year-amin);
      Vect_Year_Pred[Vect_Year_Pred>= bb] = bb
      Vect_Year_Pred[Vect_Year_Pred<= aa] = aa #Not needed?

      scyears <- (tab$Year-amin)

      Vect_Year_Pred_MSM = Vect_Year_Pred;
      Vect_Year_Pred_FSW = Vect_Year_Pred;
      Vect_Year_Pred_PregWom = Vect_Year_Pred;
      Vect_Year_Pred_LRWom = Vect_Year_Pred;
      Vect_Year_Pred_LRMen = Vect_Year_Pred;
      Vect_Year_Pred_BloodDonors = Vect_Year_Pred;
      scyears_PregWom <- scyears
      scyears_LRWom <- scyears
      scyears_LRMen <- scyears
      scyears_BloodDonor <- scyears

      typemodel=TRUE #

      icc = 0
      gmodestim0PID <- 0;

      tabb_all <- tab;
      tabb_all$estim <- tabb_all$error <- NA
      #if(!is.null(fn_filter_survey_LR))
      {
        tab <- subset(tabb_all,!(Data_type%in%c("MSM", "FSW"))) #SubPopForPrev

        if(nrow(tab)==0) return(0)

        vvprob = tab$Weights
        err_syph = lapply(1:nrow(tab),function(ll) which(tab$DX_Code==tab$DX_Code[ll] & tab$Data_type==tab$Data_type[ll]))
        wei_syph = lapply(1:nrow(tab),function(ll) {ww0 = tab$Weights[which(tab$DX_Code==tab$DX_Code[ll] & tab$Data_type==tab$Data_type[ll])]; if(length(ww0)>0) ww0=ww0/sum(ww0); ww0} )

        err_samp = function(ll) {vv=as.vector(err_syph[[ll]]); sprob= as.vector(wei_syph[[ll]]) ; if (length(vv)==1) res= vv else res=sample(vv, prob=sprob, replace=TRUE)[1]}
        ERR_SAMP = function(ll) as.vector(sapply(ll,err_samp))

        tab$Weights = tab$Weights/sum(tab$Weights) # Normalize the weights so that the sample size is not affected
        tab$Weights = tab$Weights*sum(tab$N_tested,na.rm=T)

        amin = floor(min(tab$Year))-(FB_ProjMax_B+0.5)
        bmax = ceiling(max(tab$Year)) + (FB_ProjMax_F+0.5)
        bmoinsa = bmax-amin

        aa <- 0.5
        bb <- bmoinsa-0.5

        Vect_Year_Pred = (Out_syphilis$year-amin);
        Vect_Year_Pred[Vect_Year_Pred>= bb] = bb
        Vect_Year_Pred[Vect_Year_Pred<= aa] = aa #Not needed?

        scyears <- (tab$Year-amin)

        Vect_Year_Pred_PregWom = Vect_Year_Pred;
        Vect_Year_Pred_LRWom = Vect_Year_Pred;
        Vect_Year_Pred_LRMen = Vect_Year_Pred;
        Vect_Year_Pred_BloodDonors = Vect_Year_Pred;
        scyears_PregWom <- scyears
        scyears_LRWom <- scyears
        scyears_LRMen <- scyears
        scyears_BloodDonor <- scyears
      }

      try(gmodestim0PID  <- BestFitPIDLLik_Random(vX=scyears, vY=tab$p, vW=tab$Weights, maxknots=Fmaxknots,max_FProj=bmoinsa,ficheck=FALSE, frecov=dur_Syph)
          ,silent=T)

      idxMSM <-  which(tabb_all$Data_type=="MSM")
      idxFSW <-  which(tabb_all$Data_type=="FSW")
      if(length(unique(tabb_all[idxMSM,]$Year))<2) idxMSM <- NULL
      if(length(unique(tabb_all[idxFSW,]$Year))<2) idxFSW <- NULL

      all_estMSM <- all_estFSW <- NULL

      if(is.list(gmodestim0PID))
      {
        tab$estim = gmodestim0PID$funcestprojprev(scyears) #pred_glmfit(Years=tab$Year,res)
        tab$estimprev = tab$estim

        tab$estim = log(tab$estim/(1-tab$estim))
        tab$error = log(tab$Npos/(tab$N_tested-tab$Npos))-tab$estim#res$residuals

        prevOR <- GetLROR(dat=tab, L_surveytypes=fn_filter_survey_LR)
        temp_prev = gmodestim0PID$funcestprojprev(Vect_Year_Pred) #pred_glmfit(Years=Out_syphilis$year,res) #

        Out_syphilisGeneWom$PrevEstF = Out_syphilis$PrevEstF = temp_prev*prevOR['GeneWom']/(1-temp_prev+temp_prev*prevOR['GeneWom'])
        Out_syphilisGeneMen$PrevEstM = Out_syphilis$PrevEstM = temp_prev*prevOR['GeneMen']/(1-temp_prev+temp_prev*prevOR['GeneMen'])
        Out_syphilisPregWom$PrevEstPregWom = temp_prev*prevOR['PregWom']/(1-temp_prev+temp_prev*prevOR['PregWom'])
        Out_syphilisBloodDo$PrevEstBloodDo = temp_prev*prevOR['BloodDo']/(1-temp_prev+temp_prev*prevOR['BloodDo'])

        Out_syphilisPregWom$CasePrevEstPregWom = Out_syphilisPregWom$PrevEstPregWom
        Out_syphilisBloodDo$CasePrevEstBloodDo = Out_syphilisBloodDo$PrevEstBloodDo
        #

        Out_syphilisGeneWom$CasePrevEstF = Out_syphilis$CasePrevEstF = Out_syphilis$PrevEstF
        Out_syphilisGeneMen$CasePrevEstM = Out_syphilis$CasePrevEstM = Out_syphilis$PrevEstM
        Out_syphilis$CasePrevEstMPlusF = Out_syphilis$PrevEstF

        Out_syphilis$InciEstF = gmodestim0PID$funcestprojinci(Vect_Year_Pred) #pred_glmfit(Years=Out_syphilis$year,res) #

        Out_syphilis$InciEstM = MtoFRatio*Out_syphilis$InciEstF
        Out_syphilis$CaseInciEstF = MtoFRatio*Out_syphilis$InciEstF
        Out_syphilis$CaseInciEstM = MtoFRatio*Out_syphilis$InciEstF
        Out_syphilis$CaseInciEstMPlusF = MtoFRatio*Out_syphilis$InciEstF

        ##
        if(length(idxMSM)>=2)
        {
          temp_MSM <- tabb_all[idxMSM,]
          modMSM <- glm(Prevalence/100~Year, data=temp_MSM, family = quasibinomial(link="logit"))

          aminMSM = floor(min(temp_MSM$Year))- FB_ProjMax_B
          bmaxMSM = ceiling(max(temp_MSM$Year)) + FB_ProjMax_F

          Vect_Year_Pred_MSM = Out_syphilis$year;
          Vect_Year_Pred_MSM[Vect_Year_Pred_MSM>= bmaxMSM] = bmaxMSM
          Vect_Year_Pred_MSM[Vect_Year_Pred_MSM<= aminMSM] = aminMSM #Not needed?

          all_estMSM$prev <- expit(coef(modMSM)[1]+coef(modMSM)[2]*(Vect_Year_Pred_MSM))
          Out_syphilisKPs$PrevEstMSM <- all_estMSM$prev

          tabb_all$estim[idxMSM] = coef(modMSM)[1] +coef(modMSM)[2]*tabb_all$Year[idxMSM] #pred_glmfit(Years=tab$Year,res)
          tabb_all$error[idxMSM] = log(tabb_all$Npos[idxMSM]/(tabb_all$N_tested[idxMSM]-tabb_all$Npos[idxMSM]))-tabb_all$estim[idxMSM]#res$residuals
        } else
        {
          Out_syphilisKPs$PrevEstMSM <- expit(-log(LRtoHRPOR$MtoMSM$POR) + logit(Out_syphilis$PrevEstM))
        }
        Out_syphilisKPs$CasePrevEstMSM <- Out_syphilisKPs$PrevEstMSM

        if(length(idxFSW)>=2)
        {
          temp_FSW <- tabb_all[idxFSW,]

          aminFSW = floor(min(temp_FSW$Year))- FB_ProjMax_B
          bmaxFSW = ceiling(max(temp_FSW$Year)) + FB_ProjMax_F

          Vect_Year_Pred_FSW = Out_syphilis$year;
          Vect_Year_Pred_FSW[Vect_Year_Pred_FSW>= bmaxFSW] = bmaxFSW
          Vect_Year_Pred_FSW[Vect_Year_Pred_FSW<= aminFSW] = aminFSW #Not needed?

          modFSW <- glm(Prevalence/100~Year, data=temp_FSW, family = quasibinomial(link="logit"))
          all_estFSW$prev <- expit(coef(modFSW)[1]+coef(modFSW)[2]*(Vect_Year_Pred_FSW))

          Out_syphilisKPs$PrevEstFSW <- all_estFSW$prev

          tabb_all$estim[idxFSW] = coef(modFSW)[1] +coef(modFSW)[2]*tabb_all$Year[idxFSW] #pred_glmfit(Years=tab$Year,res)
          tabb_all$error[idxFSW] = log(tabb_all$Npos[idxFSW]/(tabb_all$N_tested[idxFSW]-tabb_all$Npos[idxFSW]))-tabb_all$estim[idxFSW]#res$residuals
        } else
        {
          Out_syphilisKPs$PrevEstFSW <- expit(-log(LRtoHRPOR$FtoFSW$POR) + logit(Out_syphilis$PrevEstF))
        }
        Out_syphilisKPs$CasePrevEstFSW <- Out_syphilisKPs$PrevEstFSW

      } else
      {
        icc=icc+1;
        if(icc>=10)
        {
          Out_syphilisGeneWom$PrevEstF = Out_syphilis$PrevEstF[] = NA #
          Out_syphilisGeneMen$PrevEstM[] = Out_syphilis$PrevEstM[] = NA
          Out_syphilisGeneWom$CasePrevEstF[] = Out_syphilis$CasePrevEstF[] = NA #
          Out_syphilisGeneMen$CasePrevEstM[] = Out_syphilis$CasePrevEstM[] = NA

          Out_syphilis$InciEstM[] = NA #
          Out_syphilis$CaseInciEstF[] = NA #
          Out_syphilis$CaseInciEstM[] = NA #
          Out_syphilis$CaseInciEstMPlusF[] = NA #

          Out_syphilisKPs$PrevEstMSM[] <- NA
          Out_syphilisKPs$CasePrevEstMSM[] <- NA
          Out_syphilisKPs$PrevEstFSW[] <- NA
          Out_syphilisKPs$CasePrevEstFSW[] <- NA

          Out_syphilisPregWom$PrevEstPregWom[] <- NA
          Out_syphilisPregWom$CasePrevEstPregWom[] <- NA
        }
      }#End if(is.list(gmodestim0PID))

      #Bootstrap
      count=0
      if(is.list(gmodestim0PID))
      {
        while(1)
        {
          count=count+1
          #Syphilis
          u = runif(1)
          tab$p = tab$estim+tab$error[ERR_SAMP(1:length(tab$error))]
          tab$p = exp(tab$p)/(1+exp(tab$p))

          tab$alpha = tab$p*(tab$N_tested-1)
          tab$beta = (1-tab$p)*(tab$N_tested-1)
          tab$Npos = (qbeta(rep(u,length(tab$Npos)),tab$alpha,tab$beta)*tab$N_tested)

          tab$p = tab$Npos/tab$N_tested*runif(length(tab$Npos),0.75,1.25);
          tab.boot=tab
          dur_Syph_Boot = dur_Syph*runif(1,.5,1.5)

          gmodestim0PID$recov = dur_Syph_Boot
          icc = 0
          gmodestim0PID_Boot=0
          while(TRUE)
          {
            scyears_boot <- (tab.boot$Year-amin)
            try(gmodestim0PID_Boot <- BestFit2PIDLLik_Random(vX=scyears_boot,vY=tab.boot$p,vW=tab.boot$Weights, gmodestim0PID))
            if(is.list(gmodestim0PID_Boot))
            {
              tab.boot$estimprev = gmodestim0PID_Boot$funcestprojprev(scyears_boot)

              prevOR.boot <- GetLROR(dat=tab.boot, L_surveytypes=fn_filter_survey_LR)
              temp_prev.boot = gmodestim0PID_Boot$funcestprojprev(Vect_Year_Pred) #pred_glmfit(Years=Out_syphilis$year,res) #

              resboot_SyphilisPrevGeneWom[count,] = resboot_SyphilisPrevF[count,] = temp_prev.boot*prevOR.boot["GeneWom"]/(1-temp_prev.boot+temp_prev.boot*prevOR.boot["GeneWom"])#gmodestim0PID_Boot$funcestprojprev(Vect_Year_Pred)#pred_glmfit(Out_syphilis$year,res)
              resboot_SyphilisInciF[count,] = gmodestim0PID_Boot$funcestprojinci(Vect_Year_Pred)#pred_glmfit(Out_syphilis$year,res


              resboot_SyphilisPrevGeneMen[count,] = resboot_SyphilisPrevM[count,] = temp_prev.boot*prevOR.boot["GeneMen"]/(1-temp_prev.boot+temp_prev.boot*prevOR.boot["GeneMen"])
              resboot_SyphilisPrevPregWom[count,] = temp_prev.boot*prevOR.boot["PregWom"]/(1-temp_prev.boot+temp_prev.boot*prevOR.boot["PregWom"])

              resboot_SyphilisPrevBloodDo[count,] = temp_prev.boot*prevOR.boot["BloodDo"]/(1-temp_prev.boot+temp_prev.boot*prevOR.boot["BloodDo"])
              break;
            } else
            {
              icc=icc+1
            }
            if(icc>=10)
            {
              resboot_SyphilisPrevGeneWom[count,] = resboot_SyphilisPrevF[count,] = NA
              resboot_SyphilisInciF[count,] = NA
              resboot_SyphilisPrevGeneMen[count,] = resboot_SyphilisPrevM[count,] = NA
              resboot_SyphilisPrevPregWom[count,] <- NA
              resboot_SyphilisPrevBloodDo[count,] <- NA
              break;
            }
          }

          zzz <-  MtoFRatio*(1+(-0.3+0.6*expit(rnorm(1))))
          resboot_SyphilisInciM[count,] = zzz*resboot_SyphilisInciF[count,]

          numrep <- length(resboot_SyphilisInciM[count,])
          if(length(idxMSM)>=2)
          {
            tabb_all$p[idxMSM] = tabb_all$estim[idxMSM]+tabb_all$error[idxMSM][ERR_SAMP(1:length(tabb_all$error[idxMSM]))]
            tabb_all$p[idxMSM] = exp(tabb_all$p[idxMSM])/(1+exp(tabb_all$p[idxMSM]))

            tabb_all$alpha[idxMSM] = tabb_all$p[idxMSM]*(tabb_all$N_tested[idxMSM]-1)
            tabb_all$beta[idxMSM] = (1-tabb_all$p[idxMSM])*(tabb_all$N_tested[idxMSM]-1)
            tabb_all$Npos[idxMSM] = (qbeta(rep(u,length(tabb_all$Npos[idxMSM])),tabb_all$alpha[idxMSM],tabb_all$beta[idxMSM])*tabb_all$N_tested[idxMSM])

            tabb_all$p[idxMSM] = tabb_all$Npos[idxMSM]/tabb_all$N_tested[idxMSM]*runif(length(tabb_all$Npos[idxMSM]),0.75,1.25);

            temp_MSM <- tabb_all[idxMSM,] #tab.boot[idxMSM,]
            temp_MSM$Prevalence <- rbinom(length(idxMSM),round(temp_MSM$N_tested),temp_MSM$Prevalence/100)/round(temp_MSM$N_tested)
            modMSM <- glm(Prevalence~Year, data=temp_MSM, family = quasibinomial(link="logit"))

            #
            aminMSM = floor(min(temp_MSM$Year))- FB_ProjMax_B
            bmaxMSM = ceiling(max(temp_MSM$Year)) + FB_ProjMax_F

            Vect_Year_Pred_MSM = Out_syphilis$year;
            Vect_Year_Pred_MSM[Vect_Year_Pred_MSM>= bmaxMSM] = bmaxMSM
            Vect_Year_Pred_MSM[Vect_Year_Pred_MSM<= aminMSM] = aminMSM #Not needed?
            #

            all_estMSM$prev <- expit(coef(modMSM)[1]+coef(modMSM)[2]*(Vect_Year_Pred_MSM))

            resboot_SyphilisPrevMSM[count,] <- all_estMSM$prev
          } else
          {
            resboot_SyphilisPrevMSM[count,] <- expit(-log(LRtoHRPOR$MtoMSM$POR+rnorm(numrep,0,LRtoHRPOR$MtoMSM$sdlogPOR)) + logit(resboot_SyphilisPrevM[count,]))
          }

          if(length(idxFSW)>=2)
          {
            tabb_all$p[idxFSW] = tabb_all$estim[idxFSW]+tabb_all$error[idxFSW][ERR_SAMP(1:length(tabb_all$error[idxFSW]))]
            tabb_all$p[idxFSW] = exp(tabb_all$p[idxFSW])/(1+exp(tabb_all$p[idxFSW]))

            tabb_all$alpha[idxFSW] = tabb_all$p[idxFSW]*(tabb_all$N_tested[idxFSW]-1)
            tabb_all$beta[idxFSW] = (1-tabb_all$p[idxFSW])*(tabb_all$N_tested[idxFSW]-1)
            tabb_all$Npos[idxFSW] = (qbeta(rep(u,length(tabb_all$Npos[idxFSW])),tabb_all$alpha[idxFSW],tabb_all$beta[idxFSW])*tabb_all$N_tested[idxFSW])

            tabb_all$p[idxFSW] = tabb_all$Npos[idxFSW]/tabb_all$N_tested[idxFSW]*runif(length(tabb_all$Npos[idxFSW]),0.75,1.25);

            temp_FSW <- tabb_all[idxFSW,] #tab.boot[idxFSW,]
            temp_FSW$Prevalence <- rbinom(length(idxFSW),round(temp_FSW$N_tested),temp_FSW$Prevalence/100)/round(temp_FSW$N_tested)
            modFSW <- glm(Prevalence~Year, data=temp_FSW, family = quasibinomial(link="logit"))

            #
            aminFSW = floor(min(temp_FSW$Year))- FB_ProjMax_B
            bmaxFSW = ceiling(max(temp_FSW$Year)) + FB_ProjMax_F

            Vect_Year_Pred_FSW = Out_syphilis$year;
            Vect_Year_Pred_FSW[Vect_Year_Pred_FSW>= bmaxFSW] = bmaxFSW
            Vect_Year_Pred_FSW[Vect_Year_Pred_FSW<= aminFSW] = aminFSW #Not needed?
            #

            all_estFSW$prev <- expit(coef(modFSW)[1]+coef(modFSW)[2]*(Vect_Year_Pred_FSW))
            resboot_SyphilisPrevFSW[count,] <- all_estFSW$prev
          } else
          {
            resboot_SyphilisPrevFSW[count,] <- expit(-log(LRtoHRPOR$FtoFSW$POR+rnorm(numrep,0,LRtoHRPOR$FtoFSW$sdlogPOR)) + logit(resboot_SyphilisPrevF[count,]))
          }
          if(count>=Nboots) break
        } #End While

        CI_SyphilisPrevF = apply(resboot_SyphilisPrevF,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevF = CI_SyphilisPrevF

        CI_SyphilisPrevGeneWom = apply(resboot_SyphilisPrevGeneWom,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevGeneWom = CI_SyphilisPrevGeneWom

        CI_SyphilisPrevFSW = apply(resboot_SyphilisPrevFSW,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevFSW = CI_SyphilisPrevFSW

        CI_SyphilisInciF = apply(resboot_SyphilisInciF,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCaseInciF = apply((1-resboot_SyphilisPrevF)*resboot_SyphilisInciF,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))

        CI_SyphilisPrevM = apply(resboot_SyphilisPrevM,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevM = CI_SyphilisPrevM

        CI_SyphilisPrevGeneMen = apply(resboot_SyphilisPrevGeneMen,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevGeneMen = CI_SyphilisPrevGeneMen

        CI_SyphilisPrevMSM = apply(resboot_SyphilisPrevMSM,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevMSM = CI_SyphilisPrevMSM

        CI_SyphilisPrevPregWom = apply(resboot_SyphilisPrevPregWom,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevPregWom = CI_SyphilisPrevPregWom

        CI_SyphilisPrevBloodDo = apply(resboot_SyphilisPrevBloodDo,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevBloodDo = CI_SyphilisPrevBloodDo

        CI_SyphilisInciM = apply(resboot_SyphilisInciM,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCaseInciM = apply((1-resboot_SyphilisPrevM)*resboot_SyphilisInciM,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))

        popMPlusF <- unlist(PopSize[PopSize$ISO3==vv & !is.na(PopSize$ISO3),-1]);
        popM <- unlist(PopSizeMales[PopSizeMales$ISO3==vv & !is.na(PopSizeMales$ISO3),-1]);
        popF <- unlist(PopSizeFemales[PopSizeFemales$ISO3==vv & !is.na(PopSizeFemales$ISO3),-1]);

        popMSM <- unlist(PopSizeMSM_in[PopSizeMSM_in$ISO3==vv & !is.na(PopSizeMSM_in$ISO3),-1]);
        popFSW <- unlist(PopSizeFSW_in[PopSizeFSW_in$ISO3==vv & !is.na(PopSizeFSW_in$ISO3),-1]);

        if(is.na(vv))
        {
          popMPlusF <- popM <- popF <- rep(NA,nrow(Out_syphilis))
          popMSM <- popFSW <- rep(NA,nrow(Out_syphilis))
        }

        Out_syphilis$PopMPlusF = popMPlusF
        Out_syphilis$PopM = popM
        Out_syphilis$PopF = popF

        Out_syphilisKPs$PopMSM = popMSM
        Out_syphilisKPs$PopFSW = popFSW

        #Females
        Out_syphilisGeneWom$PrevEstF_Med <- Out_syphilis$PrevEstF_Med <- CI_SyphilisPrevF[2,]
        Out_syphilisGeneWom$PrevEstF_LoB <- Out_syphilis$PrevEstF_LoB <- CI_SyphilisPrevF[1,]
        Out_syphilisGeneWom$PrevEstF_UpB <- Out_syphilis$PrevEstF_UpB <- CI_SyphilisPrevF[3,]

        Out_syphilisGeneWom$CasePrevEstF <- Out_syphilis$CasePrevEstF <- Out_syphilis$CasePrevEstF*(popF-popFSW)
        Out_syphilisGeneWom$CasePrevMedF <- Out_syphilis$CasePrevMedF <- CI_SyphilisCasePrevF[2,]*(popF-popFSW)
        Out_syphilisGeneWom$CasePrevF_LoB <- Out_syphilis$CasePrevF_LoB <- CI_SyphilisCasePrevF[1,]*(popF-popFSW)
        Out_syphilisGeneWom$CasePrevF_UpB <- Out_syphilis$CasePrevF_UpB <- CI_SyphilisCasePrevF[3,]*(popF-popFSW)

        Out_syphilis$InciEstF_Med <- CI_SyphilisInciF[2,]
        Out_syphilis$InciEstF_LoB <- CI_SyphilisInciF[1,]
        Out_syphilis$InciEstF_UpB <- CI_SyphilisInciF[3,]

        Out_syphilis$CaseInciEstF <- Out_syphilis$CaseInciEstF*(popF)
        Out_syphilis$CaseIncMedF <- CI_SyphilisCaseInciF[2,]*(popF)
        Out_syphilis$CaseIncF_LoB <- CI_SyphilisCaseInciF[1,]*(popF)
        Out_syphilis$CaseIncF_UpB <- CI_SyphilisCaseInciF[3,]*(popF)

        #FSW
        Out_syphilisKPs$PrevEstFSW_Med <- CI_SyphilisPrevFSW[2,]
        Out_syphilisKPs$PrevEstFSW_LoB <- CI_SyphilisPrevFSW[1,]
        Out_syphilisKPs$PrevEstFSW_UpB <- CI_SyphilisPrevFSW[3,]

        Out_syphilisKPs$CasePrevEstFSW <- Out_syphilisKPs$CasePrevEstFSW*popFSW
        Out_syphilisKPs$CasePrevMedFSW <- CI_SyphilisCasePrevFSW[2,]*popFSW
        Out_syphilisKPs$CasePrevFSW_LoB <- CI_SyphilisCasePrevFSW[1,]*popFSW
        Out_syphilisKPs$CasePrevFSW_UpB <- CI_SyphilisCasePrevFSW[3,]*popFSW

        #Pregnant Women
        Out_syphilisPregWom$PrevEstPregWom_Med <- CI_SyphilisPrevPregWom[2,]
        Out_syphilisPregWom$PrevEstPregWom_LoB <- CI_SyphilisPrevPregWom[1,]
        Out_syphilisPregWom$PrevEstPregWom_UpB <- CI_SyphilisPrevPregWom[3,]

        #Boold Donors
        Out_syphilisBloodDo$PrevEstBloodDo_Med <- CI_SyphilisPrevBloodDo[2,]
        Out_syphilisBloodDo$PrevEstBloodDo_LoB <- CI_SyphilisPrevBloodDo[1,]
        Out_syphilisBloodDo$PrevEstBloodDo_UpB <- CI_SyphilisPrevBloodDo[3,]

        #General-Women
        Out_syphilisGeneWom$PrevEstGeneWom_Med <- CI_SyphilisPrevGeneWom[2,]
        Out_syphilisGeneWom$PrevEstGeneWom_LoB <- CI_SyphilisPrevGeneWom[1,]
        Out_syphilisGeneWom$PrevEstGeneWom_UpB <- CI_SyphilisPrevGeneWom[3,]
        ##***************************************************
        #Males
        Out_syphilisGeneMen$PrevEstM_Med <- Out_syphilis$PrevEstM_Med <- CI_SyphilisPrevM[2,]
        Out_syphilisGeneMen$PrevEstM_LoB <- Out_syphilis$PrevEstM_LoB <- CI_SyphilisPrevM[1,]
        Out_syphilisGeneMen$PrevEstM_UpB <- Out_syphilis$PrevEstM_UpB <- CI_SyphilisPrevM[3,]

        Out_syphilisGeneMen$CasePrevEstM <- Out_syphilis$CasePrevEstM <- Out_syphilis$CasePrevEstM*(popM-popMSM)
        Out_syphilisGeneMen$CasePrevMedM <- Out_syphilis$CasePrevMedM <- CI_SyphilisCasePrevM[2,]*(popM-popMSM)
        Out_syphilisGeneMen$CasePrevM_LoB <- Out_syphilis$CasePrevM_LoB <- CI_SyphilisCasePrevM[1,]*(popM-popMSM)
        Out_syphilisGeneMen$CasePrevM_UpB <- Out_syphilis$CasePrevM_UpB <- CI_SyphilisCasePrevM[3,]*(popM-popMSM)

        Out_syphilis$InciEstM_Med <- CI_SyphilisInciM[2,]
        Out_syphilis$InciEstM_LoB <- CI_SyphilisInciM[1,]
        Out_syphilis$InciEstM_UpB <- CI_SyphilisInciM[3,]

        Out_syphilis$CaseInciEstM <- Out_syphilis$CaseInciEstM*(popM)
        Out_syphilis$CaseIncMedM <- CI_SyphilisCaseInciM[2,]*(popM)
        Out_syphilis$CaseIncM_LoB <- CI_SyphilisCaseInciM[1,]*(popM)
        Out_syphilis$CaseIncM_UpB <- CI_SyphilisCaseInciM[3,]*(popM)

        #MSM
        Out_syphilisKPs$PrevEstMSM_Med <- CI_SyphilisPrevMSM[2,]
        Out_syphilisKPs$PrevEstMSM_LoB <- CI_SyphilisPrevMSM[1,]
        Out_syphilisKPs$PrevEstMSM_UpB <- CI_SyphilisPrevMSM[3,]

        Out_syphilisKPs$CasePrevEstMSM <- Out_syphilisKPs$CasePrevEstMSM*popMSM
        Out_syphilisKPs$CasePrevMedMSM <- CI_SyphilisCasePrevMSM[2,]*popMSM
        Out_syphilisKPs$CasePrevMSM_LoB <- CI_SyphilisCasePrevMSM[1,]*popMSM
        Out_syphilisKPs$CasePrevMSM_UpB <- CI_SyphilisCasePrevMSM[3,]*popMSM

        #General-men
        Out_syphilisGeneMen$PrevEstGeneMen_Med <- CI_SyphilisPrevGeneMen[2,]
        Out_syphilisGeneMen$PrevEstGeneMen_LoB <- CI_SyphilisPrevGeneMen[1,]
        Out_syphilisGeneMen$PrevEstGeneMen_UpB <- CI_SyphilisPrevGeneMen[3,]

        ##############################################################################
        #Males
        Out_syphilis$CasePrevEstM <- Out_syphilis$CasePrevEstM+Out_syphilisKPs$CasePrevEstMSM
        Out_syphilis$PrevEstM <- Out_syphilis$CasePrevEstM/popM

        resboot_SyphilisCasePrevM <- t(sapply(1:nrow(resboot_SyphilisPrevM), function(ii)
        {
          as.numeric(resboot_SyphilisPrevGeneMen[ii,])*(popM-popMSM)+as.numeric(resboot_SyphilisPrevMSM[ii,])*popMSM
        }))

        resboot_SyphilisPrevM <- t(sapply(1:nrow(resboot_SyphilisCasePrevM), function(ii)
        {
          as.numeric(resboot_SyphilisCasePrevM[ii,])/(popM)
        }))


        #Females
        Out_syphilis$CasePrevEstF <- Out_syphilis$CasePrevEstF+Out_syphilisKPs$CasePrevEstFSW
        Out_syphilis$PrevEstF <- Out_syphilis$CasePrevEstF/popF


        resboot_SyphilisCasePrevF <- t(sapply(1:nrow(resboot_SyphilisPrevF), function(ii)
        {
          as.numeric(resboot_SyphilisPrevGeneWom[ii,])*(popF-popFSW)+as.numeric(resboot_SyphilisPrevFSW[ii,])*popFSW
        }))

        resboot_SyphilisPrevF <- t(sapply(1:nrow(resboot_SyphilisCasePrevF), function(ii)
        {
          as.numeric(resboot_SyphilisCasePrevF[ii,])/(popF)
        }))


        #BothSexes
        Out_syphilis$CasePrevEstMPlusF <- Out_syphilis$CasePrevEstF + Out_syphilis$CasePrevEstM
        Out_syphilis$PrevEstMPlusF <- Out_syphilis$CasePrevEstMPlusF/(popM+popF)

        Out_syphilis$CaseInciEstMPlusF <- Out_syphilis$CaseInciEstF+Out_syphilis$CaseInciEstM
        Out_syphilis$InciEstMPlusF <- Out_syphilis$CaseInciEstMPlusF/(popM+popF)

        resboot_SyphilisCasePrevMPlusF <- t(sapply(1:nrow(resboot_SyphilisPrevM), function(ii)
        {
          as.numeric(resboot_SyphilisPrevGeneMen[ii,])*(popM-popMSM)+as.numeric(resboot_SyphilisPrevMSM[ii,])*popMSM+
            as.numeric(resboot_SyphilisPrevGeneWom[ii,])*(popF-popFSW)+as.numeric(resboot_SyphilisPrevFSW[ii,])*popFSW
        }))

        resboot_SyphilisPrevMPlusF <- t(sapply(1:nrow(resboot_SyphilisCasePrevMPlusF), function(ii)
        {
          as.numeric(resboot_SyphilisCasePrevMPlusF[ii,])/(popM+popF)
        }))

        ###############################
        resboot_SyphilisCaseInciMPlusF <- t(sapply(1:nrow(resboot_SyphilisInciM), function(ii)
        {
          (as.numeric((1-resboot_SyphilisPrevM[ii,])*resboot_SyphilisInciM[ii,])*popM+
             as.numeric((1-resboot_SyphilisPrevF[ii,])*resboot_SyphilisInciF[ii,])*popF)
        }))

        resboot_SyphilisInciMPlusF <- t(sapply(1:nrow(resboot_SyphilisCaseInciMPlusF), function(ii)
        {
          as.numeric(resboot_SyphilisCaseInciMPlusF[ii,])/(popM+popF)
        }))


        #########################################################################
        CI_SyphilisPrevMPlusF = apply(resboot_SyphilisPrevMPlusF,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCasePrevMPlusF = apply(resboot_SyphilisCasePrevMPlusF,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))

        CI_SyphilisInciMPlusF = apply(resboot_SyphilisInciMPlusF,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))
        CI_SyphilisCaseInciMPlusF = apply(resboot_SyphilisCaseInciMPlusF,2,function(x) quantile(x,c(.025,.5,.975),na.rm=T))

        Out_syphilis$PrevEstMPlusF_Med <- CI_SyphilisPrevMPlusF[2,]
        Out_syphilis$PrevEstMPlusF_LoB <- CI_SyphilisPrevMPlusF[1,]
        Out_syphilis$PrevEstMPlusF_UpB <- CI_SyphilisPrevMPlusF[3,]

        Out_syphilis$CasePrevMedMPlusF <- CI_SyphilisCasePrevMPlusF[2,]
        Out_syphilis$CasePrevMPlusF_LoB <- CI_SyphilisCasePrevMPlusF[1,]
        Out_syphilis$CasePrevMPlusF_UpB <- CI_SyphilisCasePrevMPlusF[3,]

        Out_syphilis$InciEstMPlusF_Med <- CI_SyphilisInciMPlusF[2,]
        Out_syphilis$InciEstMPlusF_LoB <- CI_SyphilisInciMPlusF[1,]
        Out_syphilis$InciEstMPlusF_UpB <- CI_SyphilisInciMPlusF[3,]

        Out_syphilis$CaseIncMedMPlusF <- CI_SyphilisCaseInciMPlusF[2,]
        Out_syphilis$CaseIncMPlusF_LoB <- CI_SyphilisCaseInciMPlusF[1,]
        Out_syphilis$CaseIncMPlusF_UpB <- CI_SyphilisCaseInciMPlusF[3,]

        Out_syphilis$Model = CountryModel
        col_idx <- as.numeric(sapply(c(syph_nn1[1:16],syph_nn1[17:28],syph_nninc[1:12],syph_nninc[13:27]), function(xx) which(names(Out_syphilis)==xx)))
        Out_syphilis <- Out_syphilis[,col_idx]

        Out_syphilisKPs <- Out_syphilisKPs[,c(1:4,7,11:13,5,17:19,8,14:16,6,20:22,9:10)]
        names(Out_syphilisKPs)[1:12] <- tttKPs[1:12]
        names(Out_syphilisKPs)[13:20] <- tttElKPs

        Out_syphilisPregWom <- Out_syphilisPregWom[,c(1:5, 7:9,6)]
        names(Out_syphilisPregWom)[1:4] <- tttKPs[1:4]
        names(Out_syphilisPregWom)[5:8] <- tttKPs[13:16]

        Out_syphilisGeneMen <- Out_syphilisGeneMen[,c(1:5,7:9,6, 10:15)]
        Out_syphilisGeneMen <- Out_syphilisGeneMen[,1:8]
        names(Out_syphilisGeneMen)[1:4] <- tttKPs[1:4]
        names(Out_syphilisGeneMen)[5:8] <- c("EstimatePrevGeneMen", "MedianPrevGeneMen", "PrevLB_2.5%GeneMen", "PrevUB_97.5%GeneMen")

        Out_syphilisGeneWom <- Out_syphilisGeneWom[,c(1:5,7:9,6, 10:15)]
        Out_syphilisGeneWom <- Out_syphilisGeneWom[,1:8]
        names(Out_syphilisGeneWom)[1:4] <- tttKPs[1:4]
        names(Out_syphilisGeneWom)[5:8] <- c("EstimatePrevGeneWom", "MedianPrevGeneWom", "PrevLB_2.5%GeneWom", "PrevUB_97.5%GeneWom")

        Out_syphilisBloodDo <- Out_syphilisBloodDo[,c(1:5, 7:9,6)]
        names(Out_syphilisBloodDo)[1:4] <- tttKPs[1:4]
        names(Out_syphilisBloodDo)[5:8] <- c("EstimatePrevBloodDo", "MedianPrevBloodDo", "PrevLB_2.5%BloodDo", "PrevUB_97.5%BloodDo")

        infoRun = data.frame(fitted=rep("Run", nrow(Out_syphilis)), DLastRun =as.character(rep(Sys.time(),nrow(Out_syphilis))))
        ctr_res <- list(Out_syphilis=Out_syphilis, infoRun=infoRun, CountryDataUse=CountryDataUse, CountryDataUse=CountryDataUse,
                        Out_syphilisKPs=Out_syphilisKPs, Out_syphilisPregWom=Out_syphilisPregWom, Out_syphilisGeneWom=Out_syphilisGeneWom,
                        Out_syphilisGeneMen=Out_syphilisGeneMen, Out_syphilisBloodDo=Out_syphilisBloodDo)
      } else
      {
        return(0)
      }#End if(is.list(gmodestim0PID))
    }#end if(!check_maxyear)
  }#End for loop
  #)

  stopCluster(cluster)

  #Printing results in file
  for(ctr_res in all_res)
  {
    if(length(ctr_res)<=1) next
    ################################################################################################################
    ########### End of the estimation procedure ####################################################################
    ################################################################################################################
    if(nrownum<=length(year_predict))
    {
      n_row2 = 1
      nligne_decal_titre = 1
      printed_names <- c(ttt,tttEl,tttInc[1:12],tttIncEl, names(ctr_res$Out_syphilisKPs)[-c(1:4)], names(ctr_res$Out_syphilisPregWom)[-c(1:4)],
                         names(ctr_res$Out_syphilisGeneWom)[-c(1:4)], names(ctr_res$Out_syphilisGeneMen)[-c(1:4)], names(ctr_res$Out_syphilisBloodDo)[-c(1:4)])
      openxlsx::writeData(wb,sheet=Syphilis_Rbootstrap,as.data.frame(matrix(printed_names,nrow=1)), startRow = n_row2, startCol = 1, colNames = FALSE, rowNames = FALSE)
    }

    infoRun <- data.frame(fitted=rep("Run", nrow(ctr_res$Out_syphilis)), DLastRun =as.character(rep(Sys.time(),nrow(ctr_res$Out_syphilis))))
    temp_data <- cbind(ctr_res$Out_syphilis, infoRun, ctr_res$Out_syphilisKPs[-c(1:4)], ctr_res$Out_syphilisPregWom[-c(1:4)],
                       ctr_res$Out_syphilisGeneWom[-c(1:4)], ctr_res$Out_syphilisGeneMen[-c(1:4)], ctr_res$Out_syphilisBloodDo[-c(1:4)])
    openxlsx::writeData(wb, sheet=Syphilis_Rbootstrap, temp_data,colNames = FALSE,rowNames = FALSE,
                        startRow = 1+nligne_decal_titre+nrownum, startCol = 1) #ctr_res$Out_syphilis
    nrownum = nrownum+nsep+nrow(ctr_res$Out_syphilis)
    All_CountryDataUse = rbind(All_CountryDataUse,ctr_res$CountryDataUse)
  }

  if(nrow(All_CountryDataUse)>=1)
  {
    openxlsx::writeData(wb, sheet=Syphilis_YearCheck_Glob, All_CountryDataUse, colNames = TRUE, rowNames = FALSE, startRow=1,
                        startCol=1)
  } else
  {
    print("No country fitted!!!!!")
  }

  if(nrow(SyphData[,is.element(names(SyphData),c("Country","ISO3","DiagTestAdjusteFactor"))])>=1)
  {
    openxlsx::writeData(wb, sheet = Syphilis_AdjTestCheck_Glob, SyphData[,is.element(names(SyphData), c("Country","ISO3","DiagTestAdjusteFactor"))],
                        colNames=TRUE, rowNames=FALSE, startRow=1, startCol = 1)
  }

  if(autosavefile) openxlsx::saveWorkbook(wb, name.out.file, overwrite = TRUE)
  all_res <- list(filename=name.out.file, SyphData=SyphData, wb=wb)
  invisible(all_res)
}

default_syph_durations <- list(Sypdur_A = 1.28, Sypdur_B = 2.42, Sypdur_C = 4.13)
charnumtransform <- function(x) as.numeric(gsub(",", "", as.character(x)))

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
RunFitSyphilis0 <- function(name.data.file,
                            name.popu.file = NULL,
                            num_boot=400,
                           in_Fmaxknots=2,
                           Fin_B_ProjMax_F=3,
                           Fin_B_ProjMax_B=3,
                           f_high_risk_adj = 1.1,
                           in_year_predict=1990:2023,
                           zero_prev_adj=1/100,
                           Sypdurs = default_syph_durations,
                           MtoFRatio=1,
                           f_LRtoHRPOR,
                           FIRST_YEAR=1990,
                           min_year_last_data=2011,
                           list_countries=NULL,
                           filter_survey = Default_List_surveytypes,
                           rautosavefile=TRUE,
                           appendtodefaultDB=FALSE)
{
  namesCol = c("Country","ISO3","ISO3_letters","WHO_region","Data_type","Data_type_code","Sex","Year"
               ,"Diagnostic_test","DX_Code","N_positive","N_tested","Prevalence",
               "Weight_for_Spectrum_fitting", "WghtNSpectrum")

  TFnamesCol = c("Country",	"ISO3",	"ISO3_letters",	"WHO region",	"Data type",	"Data type, code",
                 "Sex",	"Year",	"Diagnostic test",	"DX_Code",	 "N positive", 	 "N tested", 	"Prevalence",
                 "Weight_for_Spectrum_fitting",	"WghtNSpectrum")

  high_risk_adj = f_high_risk_adj

  SyphData <- openxlsx::read.xlsx(name.data.file,sheet="SyphData", startRow=1,
                                  cols= 1:(1 + length(namesCol)), check.names = FALSE, sep.names =" ")[, -c(7, 9)]

  names(SyphData)[names(SyphData)=="Weight for Spectrum fitting"] <- "Weight_for_Spectrum_fitting"
  names(SyphData)[names(SyphData)=="Midpoint study year"] <- "Year"
  names(SyphData)[names(SyphData)=="Population type"] <- "Data type"
  names(SyphData)[names(SyphData)=="Population code"] <- "Data type, code"

  suppressWarnings(SyphData$Prevalence <- charnumtransform(SyphData$Prevalence))#as.numeric(as.character(SyphData$Prevalence)))
  suppressWarnings(SyphData$`N positive` <- charnumtransform(SyphData$`N positive`))#as.numeric(as.character(SyphData$`N positive`)))
  suppressWarnings(SyphData$`N tested` <- charnumtransform(SyphData$`N tested`))#as.numeric(as.character(SyphData$`N tested`)))
  suppressWarnings(SyphData$Year <- charnumtransform(SyphData$Year))#as.numeric(as.character(SyphData$Year)))

  SyphData <- subset(SyphData, !is.na(Prevalence))
  SyphData$`N tested`[is.na(SyphData$`N tested`) & SyphData$Prevalence==0] <- 300
  SyphData$`N tested`[is.na(SyphData$`N tested`) & is.na(SyphData$`N positive`)] <- 300
  idx_zerprev <- which(SyphData$Prevalence==0)
  SyphData$Prevalence[idx_zerprev] = zero_prev_adj/SyphData$`N tested`[idx_zerprev]

  idx_np <-  which(!is.na(SyphData$`N positive`) & (is.na(SyphData$`N tested`)))
  SyphData$`N tested`[idx_np] <- SyphData$`N positive`[idx_np]/SyphData$Prevalence[idx_np]*100

  idx_nt <-  which(is.na(SyphData$`N positive`) & (!is.na(SyphData$`N tested`)))
  SyphData$`N positive`[idx_nt] <- SyphData$`N tested`[idx_nt]*SyphData$Prevalence[idx_nt]/100

  idx_ntp <-  which(is.na(SyphData$`N positive`) & (is.na(SyphData$`N tested`)))
  SyphData$`N positive`[idx_ntp] <- 300*SyphData$Prevalence[idx_ntp]/100
  SyphData$`N tested`[idx_ntp] <- 300

  SyphData$DX_Code[is.na(SyphData$DX_Code)] <- 1
  SyphData <- subset(SyphData,!is.na(SyphData$`N tested`) & SyphData$`N tested`>0)

  # Prevalence adjustment
  DiagnosticTest = Default_DiagnosticTest

  LowRisk <- c("ANC Routine screening","ANC Survey", "BloodDonor Screening Men", "Survey LowRisk Men", "BloodDonor Screening Men + Women",
               "Survey LowRisk Men+Women", "Survey LowRisk Women", "BloodDonor Screening Women")

  SyphData$Prevalence <- sapply(1:length(SyphData$Prevalence),function(ii){
    res = NA
    if(!is.na(SyphData$DX_Code[ii]))
    {
      hr_adj <- ifelse(any(LowRisk==SyphData$'Data type'[ii]),high_risk_adj,1);
      #adj = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==max(SyphData$DX_Code[ii],1))];
      adj = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==SyphData$DX_Code[ii])];
      res = SyphData$Prevalence[ii]*adj*hr_adj
    }
    res
  } )

  SyphData$DiagTestAdjusteFactor <- sapply(1:length(SyphData$Prevalence),function(ii){
    res = NA
    if(!is.na(SyphData$DX_Code[ii]))
    {
      #res = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==max(SyphData$DX_Code[ii],1))];
      res = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==SyphData$DX_Code[ii])];
    }
    res
  } )

  #End  Prevalence adjustment
  ISO3 <- SyphDatabase_ISO3
  namesCol = c(namesCol,"DiagTestAdjusteFactor")

  SyphData = SyphData[,is.element(names(SyphData),TFnamesCol)];
  genv = environment();
  ll = sapply(1:ncol(SyphData), function(jj) if(is.element(names(SyphData)[jj],TFnamesCol))
  {
    iiind = which(TFnamesCol==names(SyphData)[jj])
    names(genv$SyphData)[jj]= namesCol[iiind]
    jj
  })

  SyphData$WghtNSpectrum = SyphData$Weight_for_Spectrum_fitting
  SyphData <- subset(SyphData, Year>= FIRST_YEAR)
  data.syphilis = SyphData;

  if(appendtodefaultDB) data.syphilis <- rbind(data.syphilis,Default_SyphDataBase_20230902)
  #Syphilis
  Sypdur_A = Sypdurs$Sypdur_A # B
  Sypdur_B = Sypdurs$Sypdur_B
  Sypdur_C = Sypdurs$Sypdur_C

  duration.syphilis=data.frame(t.zone=c("A","B","C"),duration=c(Sypdur_A,Sypdur_B,Sypdur_C))
  if(!is.null(list_countries)) data.syphilis <- subset(data.syphilis, ISO3%in%list_countries)

  #if(!is.null(filter_survey)) data.syphilis <- subset(data.syphilis, Data_type%in%filter_survey)

  if(nrow(data.syphilis)==0)
  {
    warnings("No country seclected. Please check the countries names or Data type.")
  }

  fn_PopSizeBoth=Default_PopSize_20230902
  fn_PopSizeMen=Default_PopSizeMales_20230902
  fn_PopSizeWomen=Default_PopFemales_20230902

  if(!is.null(name.popu.file))
  {
    fn_PopSizeBoth = openxlsx::read.xlsx(name.popu.file,sheet="Population aged 15-49",startRow=2, cols=c(3,8:43),colNames=FALSE)
    colnames(fn_PopSizeBoth) =  c("ISO3", paste("Year_", 1990+(0:(ncol(fn_PopSizeBoth)-2)),sep=""))
    for(ii in 2:ncol(fn_PopSizeBoth))
    {
      fn_PopSizeBoth[,ii] <- as.numeric(as.character(fn_PopSizeBoth[,ii]))
    }

    fn_PopSizeMen = openxlsx::read.xlsx(name.popu.file,sheet="Males",startRow=2, cols=c(3,8:43),colNames=FALSE)
    colnames(fn_PopSizeMen) =  c("ISO3", paste("Year_", 1990+(0:(ncol(fn_PopSizeMen)-2)),sep=""))
    for(ii in 2:ncol(fn_PopSizeMen))
    {
      fn_PopSizeMen[,ii] <- as.numeric(as.character(fn_PopSizeMen[,ii]))
    }

    fn_PopSizeWomen= openxlsx::read.xlsx(name.popu.file,sheet="Females",startRow=2, cols=c(3,8:43),colNames=FALSE)
    colnames(fn_PopSizeWomen) =  c("ISO3", paste("Year_", 1990+(0:(ncol(fn_PopSizeWomen)-2)),sep=""))
    for(ii in 2:ncol(fn_PopSizeWomen))
    {
      fn_PopSizeWomen[,ii] <- as.numeric(as.character(fn_PopSizeWomen[,ii]))
    }
  }

  nnames_popsizes <- c("ISO3", paste("Year_", in_year_predict, sep=""))
  fn_PopSizeBoth <- fn_PopSizeBoth[,colnames(fn_PopSizeBoth)%in%nnames_popsizes]
  fn_PopSizeMen <- fn_PopSizeMen[,colnames(fn_PopSizeMen)%in%nnames_popsizes]
  fn_PopSizeWomen <- fn_PopSizeWomen[,colnames(fn_PopSizeWomen)%in%nnames_popsizes]

  if(ncol(fn_PopSizeBoth)!=length(nnames_popsizes)) stop("Population size data incompatible. The model is not allowed to predict before 1990 or beyond 2025")

  # Loading functions
  result = fCountryAnalysis_glob(Nboots=num_boot, fname.data.file = name.data.file, Fmaxknots=in_Fmaxknots,
                                 FB_ProjMax_F=Fin_B_ProjMax_F, FB_ProjMax_B=Fin_B_ProjMax_B,
                                 duration.syphilis, data.syphilis, PopSizeBoth=fn_PopSizeBoth,
                                 PopSizeMen=fn_PopSizeMen, PopSizeWomen=fn_PopSizeWomen, ISO3,
                                 year_predict=in_year_predict, zerprev_adj=zero_prev_adj, MtoFRatio,
                                 LRtoHRPOR = f_LRtoHRPOR,
                                 fn_min_year_last_data=min_year_last_data, autosavefile=rautosavefile,
                                 fn_filter_survey_LR = filter_survey) # Running the code
  if(!is.null(result)) class(result) <- "Syph-fit"
  invisible(result)
}

#Default prevalence Odds Ratios for MSM and FSW
OR_LRtoHR = list(MtoMSM = list(POR=0.9, sdlogPOR=0.1), FtoFSW=list(POR=0.9, sdlogPOR=0.1))

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
RunFitSyphilis1 <- function(name.data.file,
                            name.popu.file = NULL,
                            num_boot=400,
                            in_Fmaxknots=2,
                            Fin_B_ProjMax_F=3,
                            Fin_B_ProjMax_B=3,
                            f_high_risk_adj = 1.1,
                            in_year_predict=1990:2023,
                            zero_prev_adj=1/100,
                            Sypdurs = default_syph_durations,
                            DiagnosticTest = Default_DiagnosticTest,
                            MtoFRatio = 1,
                            f_LRtoHRPOR = OR_LRtoHR,
                            FIRST_YEAR=1990,
                            min_year_last_data=2011,
                            list_countries=NULL,
                            filter_survey = Default_List_surveytypes,
                            rautosavefile=TRUE,
                            appendtodefaultDB=FALSE)
{

  namesCol = c("Country","ISO3" ,"ISO3_letters", "WHO_region","Data_type", "Population_code", "Sex", "Year", "DX_Code",
               "N_positive", "N_tested", "Prevalence", "Weight", "WghtNSpectrum", "Weight_for_Spectrum_fitting")
  TFnamesCol = c("Country","ISO3", "ISO3_letters", "WHO region","Data type", "Population code", "Sex", "Year", "DX_Code",
                 "N positive", "N tested", "Prevalence", "Weight","WghtNSpectrum", "Weight_for_Spectrum_fitting","DiagTestAdjusteFactor")

  high_risk_adj = f_high_risk_adj;

  shna <- openxlsx::sheets(openxlsx::loadWorkbook(name.data.file))

  if(is.element("Data Entry",shna))
  {
    SyphData <- openxlsx::read.xlsx(name.data.file,sheet="Data Entry", startRow=1,
                                    cols= 1:(2 + length(namesCol)), check.names = FALSE, sep.names =" ")
  } else if (is.element("Combined",shna))
  {
    SyphData <- openxlsx::read.xlsx(name.data.file,sheet="Combined", startRow=1,
                                    cols= 1:(2 + length(namesCol)), check.names = FALSE, sep.names =" ")
  }else
  {
    stop("The data entry form does not seem have correct sheet names (either Data Entry or Combined) please check")
  }

  names(SyphData)[names(SyphData)=="Weight"] <- "Weight_for_Spectrum_fitting"
  names(SyphData)[names(SyphData)=="Midpoint study year"] <- "Year"
  names(SyphData)[names(SyphData)=="Population type"] <- "Data type"
  names(SyphData)[names(SyphData)=="Population code"] <- "Data type, code"
  names(SyphData)[names(SyphData)=="DX_Code - updated"] <- "DX_Code"
  names(SyphData)[names(SyphData)=="Duration estimate"] <- "Duration"

  SyphData$ISO3 <- SyphData$ISO3_letters

  SyphData$"Data type" <- sapply(SyphData$"Data type, code", function(xx){
    res <- "Other"
    #if(xx==1) res = "ANC Survey" else if(xx==2) res = "ANC Routine screening"
    if(xx==1){
      res = "ANC Survey"
    } else if(xx==2) {
      res = "ANC Routine screening"
    } else if(xx==3) {
      res = "Survey LowRisk Women"
    } else if(xx==4) {
      res = "Survey LowRisk Men"
    } else if(xx==5) {
      res = "Survey LowRisk Men+Women"
    } else if(xx==6){
      res = "BloodDonor Screening Women"
    } else if(xx==7){
      res = "BloodDonor Screening Men"
    } else if(xx==8){
      res = "BloodDonor Screening Men + Women"
    } else if (xx==9) {
      res = "FSW"
    } else if (xx%in%c(11,12)) {
      res = "MSM"
    }
    res
  })

  SyphData <- SyphData[SyphData$"Data type"%in%c("ANC Survey","ANC Routine screening","Survey LowRisk Women",
                                                 "Survey LowRisk Men", "Survey LowRisk Men+Women",
                                                 "Survey LowRisk Men + Women", "BloodDonor Screening Women",
                                                 "BloodDonor Screening Men", "BloodDonor Screening Men + Women",
                                                 "BloodDonor Screening Men+Women","FSW", "MSM"),]

  #SyphData$"Data type"[SyphData$"Data type"%in%c("ANC Routine screening","ANC Survey")] <- "Pregnant women/General-women"
  #SyphData$"Data type"[SyphData$"Data type"%in%c("BloodDonor Screening Men")] <- "Blood donors/General-men"
  #SyphData$"Data type"[SyphData$"Data type"%in%c("BloodDonor Screening Women")] <- "Blood donors/General-women"
  #SyphData$"Data type"[SyphData$"Data type"%in%c("BloodDonor Screening Men + Women")] <- "Blood donors/General-women/General-men"
  #SyphData$"Data type"[SyphData$"Data type"%in%c("Survey LowRisk Men")] <- "General-men"
  #SyphData$"Data type"[SyphData$"Data type"%in%c("Survey LowRisk Men+Women")] <- "General-women/General-men"
  #SyphData$"Data type"[SyphData$"Data type"%in%c("Survey LowRisk Men + Women")] <- "General-women/General-men"
  #SyphData$"Data type"[SyphData$"Data type"%in%c("Survey LowRisk Women")] <- "General-women"
  #SyphData <- SyphData[SyphData$"Data type"%in%c("Pregnant women/General-women","Blood donors/General-men",
  #                                                     "Blood donors/General-women","Blood donors/General-women/General-men",
  #                                                     "General-men","General-women/General-men","General-women", "FSW", "MSM"),]


  filter_survey_mod <- filter_survey
  #if(!is.null(filter_survey))
  #{
  #  filter_survey_mod[which(!(filter_survey_mod%in%c("ANC Survey","ANC Routine screening","MSM", "FSW")))] <- "Other"
  #}

  suppressWarnings(SyphData$Prevalence <- as.numeric(as.character(SyphData$Prevalence)))
  suppressWarnings(SyphData$`N positive` <- as.numeric(as.character(SyphData$`N positive`)))
  suppressWarnings(SyphData$`N tested` <- as.numeric(as.character(SyphData$`N tested`)))
  suppressWarnings(SyphData$`Duration` <- as.numeric(as.character(SyphData$`Duration`)))
  suppressWarnings(SyphData$Year <- as.numeric(as.character(SyphData$Year)))

  SyphData <- subset(SyphData, !is.na(Prevalence))
  SyphData$`N tested`[is.na(SyphData$`N tested`) & SyphData$Prevalence==0] <- 300
  SyphData$`N tested`[is.na(SyphData$`N tested`) & is.na(SyphData$`N positive`)] <- 300
  idx_zerprev <- which(SyphData$Prevalence==0)
  SyphData$Prevalence[idx_zerprev] = zero_prev_adj/SyphData$`N tested`[idx_zerprev]

  idx_np <-  which(!is.na(SyphData$`N positive`) & (is.na(SyphData$`N tested`)))
  SyphData$`N tested`[idx_np] <- SyphData$`N positive`[idx_np]/SyphData$Prevalence[idx_np]*100

  idx_nt <-  which(is.na(SyphData$`N positive`) & (!is.na(SyphData$`N tested`)))
  SyphData$`N positive`[idx_nt] <- SyphData$`N tested`[idx_nt]*SyphData$Prevalence[idx_nt]/100

  idx_ntp <-  which(is.na(SyphData$`N positive`) & (is.na(SyphData$`N tested`)))
  SyphData$`N positive`[idx_ntp] <- 300*SyphData$Prevalence[idx_ntp]/100
  SyphData$`N tested`[idx_ntp] <- 300

  SyphData$DX_Code[is.na(SyphData$DX_Code)] <- 1
  SyphData <- subset(SyphData,!is.na(SyphData$`N tested`) & SyphData$`N tested`>0)

  LowRisk <- c("ANC Routine screening","ANC Survey", "BloodDonor Screening Men", "Survey LowRisk Men", "BloodDonor Screening Men + Women",
               "Survey LowRisk Men+Women", "Survey LowRisk Women", "BloodDonor Screening Women")

  SyphData$Prevalence = sapply(1:length(SyphData$Prevalence),function(ii){
    res = NA;
    if(!is.na(SyphData$DX_Code[ii]))
    {
      hr_adj <- ifelse(any(LowRisk==SyphData$'Data type'[ii]),high_risk_adj,1);
      #adj = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==max(SyphData$DX_Code[ii],1))];
      adj = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==SyphData$DX_Code[ii])];
      res = SyphData$Prevalence[ii]*adj*hr_adj
    }
    res
  } )

  SyphData$DiagTestAdjusteFactor = sapply(1:length(SyphData$Prevalence),function(ii){
    res = NA;
    if(!is.na(SyphData$DX_Code[ii]))
    {
      #res = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==max(SyphData$DX_Code[ii],1))];
      res = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==SyphData$DX_Code[ii])];
    }
    res
  } )

  ISO3 <- SyphDatabase_ISO3
  namesCol = c(namesCol,"DiagTestAdjusteFactor")

  SyphData$ISO3 <- sapply(SyphData$ISO3, function(xx) {
    res <- NA
    if(sum(ISO3$Alpha_ISO3==xx & !is.na(xx) & !is.na(ISO3$Alpha_ISO3))>=1)
    {
      res <- ISO3$Numeric_ISO3[ISO3$Alpha_ISO3==xx & !is.na(xx) & !is.na(ISO3$Alpha_ISO3)][1]
    }
    res
  })

  SyphData = SyphData[,is.element(names(SyphData),TFnamesCol)];
  genv = environment();
  ll = sapply(1:ncol(SyphData), function(jj) if(is.element(names(SyphData)[jj],TFnamesCol))
  {
    iiind = which(TFnamesCol==names(SyphData)[jj])
    names(genv$SyphData)[jj]= namesCol[iiind]
    jj
  })

  SyphData$WghtNSpectrum = SyphData$Weight_for_Spectrum_fitting
  SyphData <- subset(SyphData, Year>= FIRST_YEAR)
  data.syphilis = SyphData;
  if(appendtodefaultDB) data.syphilis <- rbind(data.syphilis,Default_SyphDataBase_20230902)

  #Syphilis
  Sypdur_A = Sypdurs$Sypdur_A # B
  Sypdur_B = Sypdurs$Sypdur_B
  Sypdur_C = Sypdurs$Sypdur_C

  duration.syphilis=data.frame(t.zone=c("A","B","C"),duration=c(Sypdur_A,Sypdur_B,Sypdur_C))
  if(!is.null(list_countries)) data.syphilis <- subset(data.syphilis, ISO3_letters%in%list_countries)

  #if(!is.null(filter_survey)) data.syphilis <- subset(data.syphilis, Data_type%in%filter_survey)

  if(nrow(data.syphilis)==0)
  {
    warnings("No country seclected. Please check the countries names or Data type.")
    return(NULL)
  }

  fn_PopSizeBoth=Default_PopSize_20230902;
  fn_PopSizeMen=Default_PopSizeMales_20230902;
  fn_PopSizeWomen=Default_PopFemales_20230902

  fn_PopSizeMSM = Default_PopSizeMSM_20230902
  fn_PopSizeFSW = Default_PopFSW_20230902

  if(!is.null(name.popu.file))
  {
    fn_PopSizeBoth = openxlsx::read.xlsx(name.popu.file,sheet="Population aged 15-49",startRow=2, cols=c(3,8:43),colNames=FALSE)
    colnames(fn_PopSizeBoth) =  c("ISO3", paste("Year_", 1990+(0:(ncol(fn_PopSizeBoth)-2)),sep=""))
    for(ii in 2:ncol(fn_PopSizeBoth))
    {
      fn_PopSizeBoth[,ii] <- as.numeric(as.character(fn_PopSizeBoth[,ii]))
    }

    fn_PopSizeMen = openxlsx::read.xlsx(name.popu.file,sheet="Males",startRow=2, cols=c(3,8:43),colNames=FALSE)
    colnames(fn_PopSizeMen) =  c("ISO3", paste("Year_", 1990+(0:(ncol(fn_PopSizeMen)-2)),sep=""))
    for(ii in 2:ncol(fn_PopSizeMen))
    {
      fn_PopSizeMen[,ii] <- as.numeric(as.character(fn_PopSizeMen[,ii]))
    }

    fn_PopSizeWomen= openxlsx::read.xlsx(name.popu.file,sheet="Females",startRow=2, cols=c(3,8:43),colNames=FALSE)
    colnames(fn_PopSizeWomen) =  c("ISO3", paste("Year_", 1990+(0:(ncol(fn_PopSizeWomen)-2)),sep=""))
    for(ii in 2:ncol(fn_PopSizeWomen))
    {
      fn_PopSizeWomen[,ii] <- as.numeric(as.character(fn_PopSizeWomen[,ii]))
    }

    fn_PopSizeMSM= openxlsx::read.xlsx(name.popu.file,sheet="MSM",startRow=2, cols=c(3,8:43),colNames=FALSE)
    colnames(fn_PopSizeMSM) =  c("ISO3", paste("Year_", 1990+(0:(ncol(fn_PopSizeMSM)-2)),sep=""))
    for(ii in 2:ncol(fn_PopSizeMSM))
    {
      fn_PopSizeMSM[,ii] <- as.numeric(as.character(fn_PopSizeMSM[,ii]))
    }

    fn_PopSizeFSW= openxlsx::read.xlsx(name.popu.file,sheet="FSW",startRow=2, cols=c(3,8:43),colNames=FALSE)
    colnames(fn_PopSizeFSW) =  c("ISO3", paste("Year_", 1990+(0:(ncol(fn_PopSizeFSW)-2)),sep=""))
    for(ii in 2:ncol(fn_PopSizeFSW))
    {
      fn_PopSizeFSW[,ii] <- as.numeric(as.character(fn_PopSizeFSW[,ii]))
    }
  }

  nnames_popsizes <- c("ISO3", paste("Year_", in_year_predict, sep=""))
  fn_PopSizeBoth <- fn_PopSizeBoth[,colnames(fn_PopSizeBoth)%in%nnames_popsizes]
  fn_PopSizeMen <- fn_PopSizeMen[,colnames(fn_PopSizeMen)%in%nnames_popsizes]
  fn_PopSizeWomen <- fn_PopSizeWomen[,colnames(fn_PopSizeWomen)%in%nnames_popsizes]

  fn_PopSizeMSM <- fn_PopSizeMSM[,colnames(fn_PopSizeMSM)%in%nnames_popsizes]
  fn_PopSizeFSW <- fn_PopSizeFSW[,colnames(fn_PopSizeFSW)%in%nnames_popsizes]

  if(ncol(fn_PopSizeBoth)!=length(nnames_popsizes)) stop("Population size data incompatible. The model is not allowed to predict before 1990 or beyond 2025")

  # Loading functions
  result = fCountryAnalysis_glob(Nboots=num_boot, fname.data.file = name.data.file, Fmaxknots=in_Fmaxknots,
                                 FB_ProjMax_F=Fin_B_ProjMax_F, FB_ProjMax_B=Fin_B_ProjMax_B,
                                 duration.syphilis, data.syphilis, PopSizeBoth=fn_PopSizeBoth,
                                 PopSizeMen=fn_PopSizeMen, PopSizeWomen=fn_PopSizeWomen, PopSizeMSM=fn_PopSizeMSM,
                                 PopSizeFSW  =fn_PopSizeFSW, ISO3,
                                 year_predict=in_year_predict, zerprev_adj=zero_prev_adj, MtoFRatio,
                                 LRtoHRPOR = f_LRtoHRPOR,
                                 fn_min_year_last_data=min_year_last_data, autosavefile=rautosavefile,
                                 fn_filter_survey_LR = filter_survey_mod) # Running the code
  if(!is.null(result)){
    result$ListRiskGroupsComp <- filter_survey
    class(result) <- "Syph-fit"
  }
  invisible(result)
}

#' Syphilis prevalence fits and projection
#'
#' @param name.data.file name of the entry file. This is a data base containing Syphilis' survey and program data. This file needs a special format. Please contact the author (GMahiane@avenirhealth.org) if you have not been provided with such a file yet.
#' @param name.popu.file name of the population data entry file. This shall contain populations sizes of the 15-to-49 years old. By default, it is set to NULL, which makes the program use population sizes from Spectrum 2023 (WPP 2022) populations estimates.
#' @param num_boot integer, number of boostrap replications.
#' @param in_Fmaxknots integer, maximum number of knots allowed for the second order segmented polynomial (spline) used to model the incidence trend. This is set to 3 by default.
#' @param Fin_B_ProjMax_F integer, maximum year allowed for forward extrapolation for the projection before flat-lining. The default value is 3 years.
#' @param Fin_B_ProjMax_B integer, maximum year allowed for backward extrapolation for the projection before flat-lining. The default value is 3 years.
#' @param f_high_risk_adj non negative real number, adjustment factor for adjustment of prevalence among high risk individuals. The prevalence is deflated if this parameter is lower than 1 and inflated otherwise, so this is expected to be larger than 1. The default value is 1.1
#' @param in_year_predict sequence of years incidence and prevalence should be projected. This is set to 1990 to 2023 by default.
#' @param zero_prev_adj non-zero real number in the range 0 to 1, to adjust for zero prevalence in the data set. The default value is 1/100
#' @param Sypdurs list of three non negative real numbers representing Syphilis durations for three regions of the world. The default value for this is default_syph_durations = list(Sypdur_A=1.28, Sypdur_B=2.42, Sypdur_C=4.13)
#' @param f_DiagTestParam Data frame containing diagnostic parameters. The default can be obtained using the function get_DefaultDiagTests()
#' @param MtoFRatio non-zero real number representing an estimate of the male to female ratio. This is one by default.
#' @param LRtoHRPOR list containing Prevalence Odds Ratios for estimating prevalence among MSM and FSW. The default parameter can be obtained by typing get_DefaultPOR() in the console.
#' @param FIRST_YEAR non-zero real number indicating the earliest year for projections to start. This is 1990 by default.
#' @param min_year_last_data non-zero real number indicating the earliest year for data to be used. This is 2011 by default.
#' @param list_countries Vector with (alpha) ISO3 codes of countries estimates are to be obtained.
#' @param filter_survey Vector with names of data types that are to be included in the analysis.
#' @param f_autosavefile boolean indicating whether results should be automatically saved in to a file or not. This is set to FALSE by default
#' @return Creates and saves an xlsx file in the working directory and outputs a list of class "Syph-fit" with the name of that file and the prevalence data used. The created file has a name with the suffix "DDMMYY_zAdj_xx_out", where DDMMYY is the data when the fit was completed, and xx is the prevalence adjustment factor.
#' @examples Not available
RunFitSyphilis <- function(name.data.file,
                           name.popu.file = NULL,
                            num_boot=400,
                            in_Fmaxknots=2,
                           Fin_B_ProjMax_F=3,
                           Fin_B_ProjMax_B=3,
                            f_high_risk_adj = 1.1,
                           in_year_predict=1990:2023,
                            zero_prev_adj=1/100,
                            Sypdurs = default_syph_durations,
                            f_DiagnosticTest = Default_DiagnosticTest,
                            MtoFRatio=1,
                            LRtoHRPOR = OR_LRtoHR,
                            FIRST_YEAR=1990,
                           min_year_last_data=2011,
                           list_countries=NULL,
                            filter_survey = Default_List_surveytypes,
                           f_autosavefile=TRUE)
{
  tempsheets <- openxlsx::sheets(openxlsx::loadWorkbook(name.data.file))
  result <- NULL
  if((is.element("Data Entry",tempsheets)) | (is.element("Combined",tempsheets)))
  {
    result <- RunFitSyphilis1(name.data.file, name.popu.file, num_boot,
                             in_Fmaxknots, Fin_B_ProjMax_F, Fin_B_ProjMax_B,
                             f_high_risk_adj,in_year_predict,
                             zero_prev_adj, Sypdurs,
                             f_DiagnosticTest,
                             MtoFRatio,
                             f_LRtoHRPOR = LRtoHRPOR,
                             FIRST_YEAR, min_year_last_data, list_countries,
                             filter_survey, rautosavefile=f_autosavefile)
  } else if(is.element("SyphData",tempsheets))
  {
    result <- RunFitSyphilis0(name.data.file, name.popu.file, num_boot,
                             in_Fmaxknots, Fin_B_ProjMax_F,  Fin_B_ProjMax_B,
                             f_high_risk_adj,in_year_predict,
                             zero_prev_adj, Sypdurs, MtoFRatio,
                             FIRST_YEAR, min_year_last_data, list_countries,
                             filter_survey, rautosavefile=f_autosavefile)
  } else
  {
    stop("The data base is not appropriate")
  }
  if(!is.null(result)) class(result) <- "Syph-fit"
  invisible(result)
}

#' Save Syphilis' prevalence and incidence estimates to xlsx file
#'
#' @param xSyphfit object of class "Syph-fitj". This is a list Syphilis estimates.
#' @param fname name of the xlsx file to which the results are to be saved.
#' @return boolean, TRUE upon completion of the task
#' @examples Not available
saveSyphfit <- function(xSyphfit, fname=NULL)
{
  if(class(xSyphfit)!="Syph-fit") stop("class(xSyphfit) should be xCSProj")
  fnameout <- fname
  if(is.null(fname))
  {
    chardate <- Sys.Date();
    chardate <- gsub("-","", chardate)
    fnameout <- paste(paste("CSProj", chardate, sep="_"),".xlsx",sep="")
  }
  openxlsx::saveWorkbook(xSyphfit$wb, file=fnameout, overwrite = T)
  return(TRUE)
}

#' Generate Templates for Syphilis data entry forms
#' @param ctr_iso3 Optional, vector of alpha numeric ISO codes of countries for which the Samples should be obtained from. This is NULL by default, implying that the template does not have any data.
#' @return Creates and saves an xlsx file in the working directory
#' @examples GetTemplateSyphData()
GetTemplateSyphData <-  function(ctr_iso3=NULL)
{
  if(is.null(ctr_iso3))
  {
    tempplate <- Default_SyphDataBase_20230902[-(1:nrow(Default_SyphDataBase_20230902)),]
  } else
  {
    tempplate <- subset(Default_SyphDataBase_20230902, ISO3%in%ctr_iso3)
  }

  outfilename <- "Template.xlsx"
  if(!is.null(ctr_iso3)) outfilename <- paste(ctr_iso3,outfilename, sep="_")
  openxlsx::write.xlsx(tempplate,tempplate)
  res <- paste("Saved in", outfilename)
  return(res)
}


#' Congenital Syphilis estimation and projection using results from the procedure RunFitSyphilis
#'
#' @param syphfitfile name of a file created from the function RunFitSyphilis. This is a data base containing Syphilis estimates.
#' @param list_countries Vector containing the ISO codes of countries Congenital Syphilis estimates are to be estimated.
#' @param proj_years integer vector of years for which estimates and projections are to be made. By default, this is the sequence 1990:2025
#' @param min_year non-zero real number indicating the earliest year for data to be used. This is 2011 by default.
#' @param rSyph_preg Ratio of Syphilis prevalence among all women to prevalence among pregnant women. This is 1 by default.
#' @param CSinputfiles list of four file names for (1) prevalence: syphilis prevalence used for prevalence estimates, (2) screening: data for Syphilis screen among pregnant women, (3) csdb: congenital syphilis prevalence, and (4) demographics: all demographic data by country. By default, this is set to NULL, which means data 2023 default data bases are used.
#' @return A list of class "CSProj", of dataframes CongenDataOut, LongCongenDataOutForPlots, RegCSABO, LongRegCSABO, RegDataPrevInc, LongRegDataPrevIncForPlots containing national and regional congenital Syphilis estimates.
#' @examples Not available
#' @keywords internal
#' @export
CalcCS_p <- function(syphfitfile, list_countries=NULL, proj_years=1990:2025,min_year=2011, rSyph_preg=1, CSinputfiles=NULL, f_DiagnosticTest=Default_DiagnosticTest)#CalcCS <- function(xSyphfit, proj_years=1990:2025)
{
  if(!is.null(CSinputfiles))
  {
    pre_res_0 <- get_precalcsCS(CSinputfiles, f_DiagnosticTest)
    CongenDataIn <- pre_res_0$CongenDataIn
    SDGRegions <- pre_res_0$SDGRegions
    RiskABO_Asumptions <- pre_res_0$RiskABO_Asumptions
    LongEarlyANCData <- pre_res_0$LongEarlyANCData
    SE_Asumptions <- pre_res_0$SE_Asumptions
    TimingANC_Asumptions <- pre_res_0$TimingANC_Asumptions
    EarlyANC <- pre_res_0$EarlyANC
    fn_impute <- pre_res_0$fn_impute
    rm(pre_res_0)

    pre_res_1 <- get_rawCS(CSinputfiles, f_DiagnosticTest)
    CongenDataRaw <- pre_res_1$CongenDataRaw
    SyphDataRaw <- pre_res_1$SyphDataRaw
    rm(pre_res_1)
  }

  fittedprevIncfile <- syphfitfile
  IncideAndPrevalence <- openxlsx::read.xlsx(fittedprevIncfile, sheet="SYPH_RBootstrap_All")
  if(is.null(IncideAndPrevalence)) return(NULL)

  names(IncideAndPrevalence)[names(IncideAndPrevalence)=="WHO_Region"] <- "WHO Region"
  prevdata <- SyphDataRaw

  ctr_with_good_data <- c()
  for(ctr_iso in unique(prevdata$ISO3_letters))
  {
    temp_ctr_data <- subset(prevdata,ISO3_letters==ctr_iso & prevdata$Year>=min_year)
    if(nrow(temp_ctr_data)>=1)
    {
      ctr_with_good_data <- c(ctr_with_good_data,ctr_iso)
    }
  }

  IncideAndPrevalence <- subset(IncideAndPrevalence,ISO3%in%ctr_with_good_data)

  IncideAndPrevalence$SDG_Region <- sapply(IncideAndPrevalence$ISO3, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  all_countries_iso <- unique(CongenDataIn$`ISO code`)

  if(!is.null(list_countries))
  {
    #all_countries_iso <- unique(CongenDataIn$`ISO code`[is.element(CongenDataIn$Country,list_countries)])
    all_countries_iso <- unique(CongenDataIn$`ISO code`[is.element(CongenDataIn$`ISO code`,list_countries)])
  }
  results <- list()

  CongenDataOut <- data.frame();
  for(isoctr in all_countries_iso)
  {
    temp_ctr <- lapply(proj_years, function(xx){
      res <- CongenDataIn[1,];
      res[1,] <- NA
      idxyy <- which(CongenDataIn$`ISO code`==isoctr& CongenDataIn$Year==xx)
      if(length(idxyy)>=1)
      {
        res <- CongenDataIn[idxyy[1],]
      } else
      {
        idx_ctr <- which(CongenDataIn$`ISO code`==isoctr)
        res$Year = xx
        res$Country = CongenDataIn$Country[idx_ctr[1]];
        res$`ISO code` = CongenDataIn$`ISO code`[idx_ctr[1]];
        res$`ISO3nmb` = CongenDataIn$`ISO3nmb`[idx_ctr[1]];
        res$`WHO Region` = CongenDataIn$`WHO Region`[idx_ctr[1]];
      }
      res
    })

    nn <- names(CongenDataIn[1,])
    temp_ctr <- plyr::ldply (temp_ctr, data.frame)
    names(temp_ctr) <- nn
    #
    temp_ctr_prev <- lapply(proj_years, function(xx){
      idxyy <- which(IncideAndPrevalence$ISO3==isoctr& IncideAndPrevalence$Year==xx)
      tncol <- ncol(IncideAndPrevalence)-4
      res <- IncideAndPrevalence[1,-c(1:4)]
      res[1,] <- NA
      if(length(idxyy)>=1)
      {
        res <- IncideAndPrevalence[idxyy[1],-c(1:4)]
      }
      res
    })

    nn <- names(IncideAndPrevalence[1,-c(1:4)])
    temp_ctr_prev <- plyr::ldply (temp_ctr_prev, data.frame)
    names(temp_ctr_prev) <- nn

    if(nrow(temp_ctr_prev)>=1) temp_ctr <- cbind(temp_ctr,temp_ctr_prev)

    #SpectrumBirths
    temp_ctr$EstimatePrevF <- fn_impute('EstimatePrevPregWom',temp_ctr,IncideAndPrevalence)#fn_impute('EstimatePrevF',temp_ctr,IncideAndPrevalence)
    temp_ctr$'PrevLB_2.5%F' <- fn_impute('PrevLB_2.5%PregWom',temp_ctr,IncideAndPrevalence)#fn_impute('PrevLB_2.5%F',temp_ctr,IncideAndPrevalence)
    temp_ctr$'PrevUB_97.5%F' <- fn_impute('PrevUB_97.5%PregWom',temp_ctr,IncideAndPrevalence)#fn_impute('PrevUB_97.5%F',temp_ctr,IncideAndPrevalence)

    for(ii in seq_len(length(temp_ctr$'PrevUB_97.5%F')))
    {
      if(is.na(temp_ctr$'PrevUB_97.5%F'[ii]) | is.na(temp_ctr$EstimatePrevF[ii])) next
      if(temp_ctr$'PrevUB_97.5%F'[ii] < temp_ctr$EstimatePrevF[ii])
      {
        temp_ctr$'PrevUB_97.5%F'[ii] <- temp_ctr$EstimatePrevF[ii]+1.96*sqrt(temp_ctr$EstimatePrevF[ii]/4)
        temp_ctr$'PrevLB_2.5%F'[ii] <- max(0,temp_ctr$EstimatePrevF[ii]-1.96*sqrt(temp_ctr$EstimatePrevF[ii]/4))
      }
    }

    #Adult prevalence
    temp_ctr$'Adult prev. USED -- Spectrum or other' <- temp_ctr$EstimatePrevF#/rSyph_preg;
    temp_ctr$'LB, Adult prev. USED -- Spectrum or other' <- temp_ctr$`PrevLB_2.5%F`#/rSyph_preg;
    temp_ctr$'UB, Adult prev. USED -- Spectrum or other' <- temp_ctr$`PrevUB_97.5%F`#/rSyph_preg;

    #Prevalence among pregnant women
    temp_ctr$'Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)' <- temp_ctr$EstimatePrevF/rSyph_preg;
    temp_ctr$'LB on maternal prevalence, incl. Imputed' <- temp_ctr$`PrevLB_2.5%F`/rSyph_preg;
    temp_ctr$'UB on maternal prevalence, incl. Imputed' <- temp_ctr$`PrevUB_97.5%F`/rSyph_preg;

    temp_ctr$'Treated (%)' <- fn_impute('Treated (%)',temp_ctr,CongenDataIn)
    temp_ctr$'Syphilis-tested (1st ANC, %)' <- fn_impute('Syphilis-tested (1st ANC, %)',temp_ctr,CongenDataIn)
    temp_ctr$'Women with >= 1 ANC visit (%)' <- fn_impute('Women with >= 1 ANC visit (%)',temp_ctr,CongenDataIn)

    #Syphilis-infected pregnancies
    temp_ctr$'Syphilis-infected pregnancies' <- temp_ctr$'Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)'*temp_ctr$Pregnancies;
    temp_ctr$'Treated mothers' <- temp_ctr$`Treated (%)`/100*temp_ctr$'Syphilis-tested (1st ANC, %)'/100*temp_ctr$`Women with >= 1 ANC visit (%)`/100;
    temp_ctr$'Untreated mothers' <- 1-temp_ctr$'Treated mothers';

    idx_risk_untrea <- which(RiskABO_Asumptions[,1]=="Risk, from UNTREATED maternal syphilis")
    idx_risk_trea <- which(RiskABO_Asumptions[,1]=="Risk, from TREATED maternal syphilis")

    temp_ctr$'ABO cases (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$'Untreated mothers'*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]+temp_ctr$'Treated mothers'*RiskABO_Asumptions$`All outcomes`[idx_risk_trea])

    temp_ctr$'CS cases (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$'Untreated mothers'+temp_ctr$'Treated mothers'*RiskABO_Asumptions$`All outcomes`[idx_risk_trea])

    temp_ctr$'ABO, not seen in ANC (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (1-temp_ctr$`Women with >= 1 ANC visit (%)`/100)*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]

    temp_ctr$'ABO, ANC women not screened (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (1-temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]

    temp_ctr$'ABO, ANC-screened women not treated (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100*temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*(1-temp_ctr$`Treated (%)`/100)*
      RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]

    temp_ctr$'ABO = CS, ANC women treated (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100*temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*temp_ctr$`Treated (%)`/100*
      RiskABO_Asumptions$`All outcomes`[idx_risk_trea]

    temp_ctr$'CS cases, not seen in ANC (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (1-temp_ctr$`Women with >= 1 ANC visit (%)`/100)

    temp_ctr$'CS cases, ANC women not screened (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100)*(1-temp_ctr$`Syphilis-tested (1st ANC, %)`/100)

    temp_ctr$'CS cases, ANC-screened women not treated (2012 method)' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100)*temp_ctr$`Syphilis-tested (1st ANC, %)`/100*(1-temp_ctr$`Treated (%)`/100)

    dftt <- subset(LongEarlyANCData,ISO3==isoctr & LongEarlyANCData$`Start year`%in%proj_years)
    aborisktreated <- rep(0, nrow(temp_ctr))
    if(nrow(dftt)==0)
    {
      dftt_temp <- subset(LongEarlyANCData,SDG_Region== temp_ctr$SDG_Region[1])
      aborisktreated <- sapply(unique(temp_ctr$Year), function(xx){
        tresu <- 0;
        rr <- dftt_temp$`ABO risk, average, treated women`[dftt_temp$`Start year`==xx]
        if(length(rr)>=1) tresu <- mean(rr,na.rm=T)
        tresu
      } )
    }else
    {
      aborisktreated <- dftt$`ABO risk, average, treated women`
    }

    temp_ctr$'ABO cases' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Untreated mothers`*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]+temp_ctr$`Treated mothers`*aborisktreated)

    temp_ctr$'CS cases' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Untreated mothers`+temp_ctr$`Treated mothers`*aborisktreated)

    temp_ctr$'ABO, not seen in ANC' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (1-temp_ctr$`Women with >= 1 ANC visit (%)`/100)*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]

    temp_ctr$'ABO, ANC women not screened' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100)*(1-temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]

    temp_ctr$'ABO, ANC-screened women not treated' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100)*(temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*(1-temp_ctr$`Treated (%)`/100)*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]

    temp_ctr$'ABO = CS, ANC women treated' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100)*(temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*(temp_ctr$`Treated (%)`/100)*aborisktreated

    temp_ctr$'CS cases, not seen in ANC' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (1-temp_ctr$`Women with >= 1 ANC visit (%)`/100)

    temp_ctr$'CS cases, ANC women not screened' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100)*(1-temp_ctr$`Syphilis-tested (1st ANC, %)`/100)

    temp_ctr$'CS cases, ANC-screened women not treated' = temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*
      (temp_ctr$`Women with >= 1 ANC visit (%)`/100)*(temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*(1-temp_ctr$`Treated (%)`/100)

    temp_ctr$'CS / 100,000 live births' = temp_ctr$'CS cases'/temp_ctr$`Live Births`*100000;

    temp_ctr$'CS / 100,000, 2012, for RANK' = ifelse(temp_ctr$Year==2012,temp_ctr$'CS / 100,000 live births',NA)

    temp_ctr$'CS / 100,000, 2016, for RANK' = ifelse(temp_ctr$Year==2016,temp_ctr$'CS / 100,000 live births',NA)

    temp_ctr$'CS report completeness' = ifelse(!is.na(temp_ctr$`Congenital syphilis case REPORTS`),temp_ctr$`Congenital syphilis case REPORTS`/temp_ctr$`CS case report rate`,NA)

    temp_ctr$'ABO/100,000 live births' = temp_ctr$`ABO cases (2012 method)`/temp_ctr$`Live Births`*100000;

    temp_ctr$'ABO risk, treated mothers' = aborisktreated;
    temp_ctr$'Pregnancies, 2008' = ifelse(temp_ctr$Year==2008,temp_ctr$Pregnancies,NA);
    temp_ctr$'Pregnancies, 2012' = ifelse(temp_ctr$Year==2012,temp_ctr$Pregnancies,NA);
    temp_ctr$'Pregnancies, 2016' = ifelse(temp_ctr$Year==2016,temp_ctr$Pregnancies,NA);
    temp_ctr$'Pregnancies, 2020' = ifelse(temp_ctr$Year==2020,temp_ctr$Pregnancies,NA);

    temp_ctr$'ANC-1, imputed?' = NA;
    temp_ctr$'Screen cov, imputed?' = NA;
    temp_ctr$'Maternal prevalence, imputed?' = NA;

    temp_ctr$'3 ANC service coverages national?' = NA;
    temp_ctr$'Spectrum trend & national 3 ANC coverages Y=1' = NA;

    temp_ctr$'Liveborn with clinical CS' = temp_ctr$'ABO cases'*RiskABO_Asumptions$`CS, risk: liveborn`[idx_risk_untrea]/RiskABO_Asumptions$`treated mothers`[idx_risk_untrea];
    temp_ctr$'Prematurity or LBW due to CS' = temp_ctr$'ABO cases'*RiskABO_Asumptions$`CS risk: prematurity or LBW`[idx_risk_untrea]/RiskABO_Asumptions$`treated mothers`[idx_risk_untrea];
    temp_ctr$'Neonatal death due to CS' = temp_ctr$'ABO cases'*RiskABO_Asumptions$`CS risk: neonatal death`[idx_risk_untrea]/RiskABO_Asumptions$`treated mothers`[idx_risk_untrea];
    temp_ctr$'Stillbirth due to CS' = temp_ctr$'ABO cases'*RiskABO_Asumptions$`CS risk: stillbirth`[idx_risk_untrea]/RiskABO_Asumptions$`treated mothers`[idx_risk_untrea];
    temp_ctr$'Asymptomatic CS' = temp_ctr$'CS cases'-temp_ctr$'ABO cases';
    temp_ctr$'Asymptomatic, mother diagnosed' = NA;
    temp_ctr$'Asymptomatic, mother NOT diagnosed' = NA;

    temp_ctr$'EMTCT elimination certification country?' = NA;
    temp_ctr$'CS case report rate (incl. Imputations as 0.3 cases for reported 0s), Elimination countries' = NA;
    temp_ctr$'ABO = CS, ANC women treated, ELIM.CIES only' = NA;
    temp_ctr$'CS cases, not seen in ANC, ELIM.CIESonly' = NA;
    temp_ctr$'CS cases, ANC women not screened, ELIM.CIESonly' = NA;
    temp_ctr$'CS cases, ANC-screened women not treated, ELIM.CIESonly' = NA;

    temp_ctr$'Graph labels' = NA;
    temp_ctr$'Maternal duration of infection' = 2.4;
    temp_ctr$'Maternal incidence/person-year' = temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`/temp_ctr$'Maternal duration of infection';

    timefirstanc <- rep(0, nrow(temp_ctr))
    if(nrow(dftt)==0)
    {
      #dftt_temp <- subset(LongEarlyANCData,Region==paste(temp_ctr$`WHO Region`[1],"O",sep=""))
      dftt_temp <- subset(LongEarlyANCData,SDG_Region==temp_ctr$SDG_Region[1])
      timefirstanc <- sapply(unique(temp_ctr$Year), function(xx){
        tresu <- 0;
        rr <- dftt_temp$`Time of first ANC, national average`[dftt$`Start year`==xx]
        if(length(rr)>=1) tresu <- mean(rr,na.rm=T)
        tresu
      } )
    }else
    {
      timefirstanc <- dftt$`Time of first ANC, national average`
    }
    temp_ctr$'Additional maternal infections, from reinfection after treatment' = (40-timefirstanc)/52*temp_ctr$'Maternal incidence/person-year'*temp_ctr$'ABO = CS, ANC women treated';

    temp_ctr$'Relative increase in CS, from reinfection' = temp_ctr$'Additional maternal infections, from reinfection after treatment'/temp_ctr$'CS cases';
    temp_ctr$'Variance on adult (not only ANC) prevalence' = ((temp_ctr$`PrevUB_97.5%M+F`-temp_ctr$`PrevLB_2.5%M+F`)/2/1.96)^2;
    temp_ctr$'Variance on MATERNAL prevalence' = ((temp_ctr$`UB on maternal prevalence, incl. Imputed`-temp_ctr$`LB on maternal prevalence, incl. Imputed`)/2/1.96)^2;

    temp_ctr$'Variance on maternal prevalence * Pregnancies^2' = temp_ctr$'Variance on MATERNAL prevalence'*temp_ctr$Pregnancies^2
    temp_ctr$'Variance on # pregnancies' = 0

    temp_ctr$'Variance on ABO prob. For treated mothers' = ((temp_ctr$'ABO risk, treated mothers'*.5)/2/1.9)^2

    temp_ctr$'LB on ANC-1 coverage' = temp_ctr$`Women with >= 1 ANC visit (%)`/100-(1-temp_ctr$`Women with >= 1 ANC visit (%)`/100)*SE_Asumptions$`Lower-bound`[SE_Asumptions$X1=="ANC-1 coverage"]
    temp_ctr$'UB on ANC-1 coverage' = temp_ctr$`Women with >= 1 ANC visit (%)`/100-(1-temp_ctr$`Women with >= 1 ANC visit (%)`/100)*SE_Asumptions$`Upper-bound`[SE_Asumptions$X1=="ANC-1 coverage"]

    temp_ctr$'Variance on ANC-1 coverage' = ((temp_ctr$'UB on ANC-1 coverage'-temp_ctr$'LB on ANC-1 coverage')/2/1.96)^2

    temp_ctr$'LB on screen coverage' = temp_ctr$`Syphilis-tested (1st ANC, %)`/100-(1-temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*SE_Asumptions$`Lower-bound`[SE_Asumptions$X1=="Screening coverage"]
    temp_ctr$'UB on screen coverage' = temp_ctr$`Syphilis-tested (1st ANC, %)`/100-(1-temp_ctr$`Syphilis-tested (1st ANC, %)`/100)*SE_Asumptions$`Upper-bound`[SE_Asumptions$X1=="Screening coverage"]
    temp_ctr$'Variance on Screen coverage' = ((temp_ctr$'LB on screen coverage'-temp_ctr$'UB on screen coverage')/2/1.96)^2

    temp_ctr$'LB on Treat coverage' = temp_ctr$`Treated (%)`/100-(1-temp_ctr$`Treated (%)`/100)*SE_Asumptions$`Lower-bound`[SE_Asumptions$X1=="Treatment coverage"]
    temp_ctr$'UB on Treat coverage' = temp_ctr$`Treated (%)`/100-(1-temp_ctr$`Treated (%)`/100)*SE_Asumptions$`Upper-bound`[SE_Asumptions$X1=="Treatment coverage"]
    temp_ctr$'Variance on Treat coverage' = ((temp_ctr$'LB on Treat coverage'-temp_ctr$'UB on Treat coverage')/2/1.96)^2

    temp_ctr$'Variance on % of mothers treated, or untreated' = (temp_ctr$`Women with >= 1 ANC visit (%)`/100)^2*temp_ctr$'Variance on ANC-1 coverage' +
      (temp_ctr$`Syphilis-tested (1st ANC, %)`/100)^2*temp_ctr$'Variance on Screen coverage' + (temp_ctr$`Treated (%)`/100)^2*temp_ctr$'Variance on Treat coverage'

    temp_ctr$'Variance on CS case number' = ((temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*(temp_ctr$`Untreated mothers`+temp_ctr$`Treated mothers`*temp_ctr$`ABO risk, treated mothers`))^2)*
      temp_ctr$`Variance on # pregnancies`+
      ((temp_ctr$Pregnancies*(temp_ctr$`Untreated mothers`+temp_ctr$`Treated mothers`*temp_ctr$`ABO risk, treated mothers`))^2)*temp_ctr$`Variance on MATERNAL prevalence`+
      ((temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`)^2)*temp_ctr$`Variance on % of mothers treated, or untreated`+
      ((temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*temp_ctr$`ABO risk, treated mothers`)^2)*temp_ctr$`Variance on % of mothers treated, or untreated`+
      ((temp_ctr$`Untreated mothers`*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`)^2)*temp_ctr$`Variance on ABO prob. For treated mothers`

    temp_ctr$'Square of Variance on CS case number' = temp_ctr$'Variance on CS case number'^2

    temp_ctr$'SE on CS case number' = sqrt(temp_ctr$'Variance on CS case number')

    temp_ctr$'Variance on ABO case number' = ((temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*(temp_ctr$`Untreated mothers`*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]+temp_ctr$`Treated mothers`*temp_ctr$`ABO risk, treated mothers`))^2)*temp_ctr$`Variance on # pregnancies`+
      ((temp_ctr$Pregnancies*(temp_ctr$`Untreated mothers`*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea]+temp_ctr$`Treated mothers`*temp_ctr$`ABO risk, treated mothers`))^2)*temp_ctr$`Variance on MATERNAL prevalence`+
      ((temp_ctr$Pregnancies*temp_ctr$'Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)'*temp_ctr$'Untreated mothers')^2)*RiskABO_Asumptions$'Variance on ABO risk assumption'[idx_risk_untrea]+
      ((temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*RiskABO_Asumptions$`All outcomes`[idx_risk_untrea])^2)*temp_ctr$`Variance on % of mothers treated, or untreated`+
      ((temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*temp_ctr$`Treated mothers`)^2*temp_ctr$`Variance on ABO prob. For treated mothers`)+
      ((temp_ctr$Pregnancies*temp_ctr$`Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)`*temp_ctr$`ABO risk, treated mothers`)^2)*temp_ctr$`Variance on % of mothers treated, or untreated`

    temp_ctr$'Square of Variance on ABO case number' = temp_ctr$'Variance on ABO case number'^2
    temp_ctr$'SE on ABO case number' = sqrt(temp_ctr$'Variance on ABO case number')

    temp_ctr$'2.5% of CS case number' = temp_ctr$`CS cases`-1.96*temp_ctr$'SE on CS case number'
    temp_ctr$'2.5% of CS case number'[temp_ctr$'2.5% of CS case number'<0] <- 0
    temp_ctr$'97.5% of CS case number' = temp_ctr$`CS cases`+1.96*temp_ctr$'SE on CS case number'

    temp_ctr$'2.5% on ABO number' = temp_ctr$`ABO cases`-1.96*temp_ctr$'SE on ABO case number'
    temp_ctr$'2.5% on ABO number'[temp_ctr$'2.5% on ABO number'<0] <- 0
    temp_ctr$'97.5% on ABO number' = temp_ctr$`ABO cases`+1.96*temp_ctr$'SE on ABO case number'

    #More variances
    temp_ctr$'Variance on CS / 100,000 live births' = (temp_ctr$'CS cases'/temp_ctr$`Live Births`)^2*(temp_ctr$'Variance on CS case number'/(temp_ctr$'CS cases'^2)+
                                                                                                        0+0)*100000^2;
    temp_ctr$'Variance on ABO/100,000 live births' = (temp_ctr$`ABO cases (2012 method)`/temp_ctr$`Live Births`)^2*(temp_ctr$'Variance on ABO case number'/(temp_ctr$`ABO cases (2012 method)`^2)+
                                                                                                                      0+0)*100000^2;

    temp_ctr$'2.5% on CS / 100,000 live births' = temp_ctr$'CS / 100,000 live births' -1.96*sqrt(temp_ctr$'Variance on CS / 100,000 live births')
    temp_ctr$'2.5% on CS / 100,000 live births'[temp_ctr$'2.5% on CS / 100,000 live births'<0] = 0
    temp_ctr$'97.5% on CS / 100,000 live births' = temp_ctr$'CS / 100,000 live births' + 1.96*sqrt(temp_ctr$'Variance on CS / 100,000 live births')

    temp_ctr$'2.5% on ABO/100,000 live births' = temp_ctr$'ABO/100,000 live births' - 1.96*sqrt(temp_ctr$'Variance on ABO/100,000 live births')
    temp_ctr$'2.5% on ABO/100,000 live births'[temp_ctr$'2.5% on ABO/100,000 live births'<0] = 0
    temp_ctr$'97.5% on ABO/100,000 live births'= temp_ctr$'ABO/100,000 live births' + 1.96*sqrt(temp_ctr$'Variance on ABO/100,000 live births')

    CongenDataOut <- rbind(CongenDataOut,temp_ctr)
  }

  ################################################################################
  ################################################################################
  ################################################################################
  ################################################################################
  base_variables <- c( "Rank, 2012 ABO cases", "Rank, 2016 ABO cases","Rank 2012, CS case RATE",
                       "Rank 2016, CS case RATE", "Country","WHO Region","ISO code","ISO3nmb",
                       "Year", "Live Births", "Still births", "Stillbirths, Blencowe & Hogan 2016",
                       "Pregnancies", "Still/Live births", "Source of Live and/or Stillbirths",
                       "Women with >= 1 ANC visit (%)", "Source of ANC1", "N tested (1st ANC visit)",
                       "N, 1st visits", "Syphilis-tested (1st ANC, %)", "N tested (any visit ANC)",
                       "N any ANC visits", "Syphilis-tested (any ANC visit,%)", "Source of Test coverage",
                       "ANC women with syphilis, treated, N","Syphilis-infected ANC", "Treated (%)",
                       "Source of Treated", "Congenital syphilis case REPORTS", "CS case report rate",
                       "EstimatePrevF", "MedianPrevF", "PrevLB_2.5%F", "PrevUB_97.5%F", "EstimatePrevM",
                       "MedianPrevM", "PrevLB_2.5%M", "PrevUB_97.5%M", "EstimatePrevM+F","MedianPrevM+F",
                       "PrevLB_2.5%M+F", "PrevUB_97.5%M+F", "CasePrev_EstF", "CasePrev_MedF", "CasePrev_LB_2.5%F",
                       "CasePrev_UB_97.5%F","CasePrev_EstM", "CasePrev_MedM", "CasePrev_LB_2.5%M", "CasePrev_UB_97.5%M",
                       "CasePrev_EstM+F", "CasePrev_MedM+F", "CasePrev_LB_2.5%M+F", "CasePrev_UB_97.5%M+F", "EstimateIncF",
                       "MedianIncF", "IncLB_2.5%F", "IncUB_97.5%F", "EstimateIncM", "MedianIncM", "IncLB_2.5%M",
                       "IncUB_97.5%M","EstimateIncM+F", "MedianIncM+F", "IncLB_2.5%M+F", "IncUB_97.5%M+F",
                       "CaseInc_EstF", "CaseInc_MedF", "CaseInc_LB_2.5%F", "CaseInc_UB_97.5%F", "CaseInc_EstM",
                       "CaseInc_MedM", "CaseInc_LB_2.5%M", "CaseInc_UB_97.5%M", "CaseInc_EstM+F", "CaseInc_MedM+F",
                       "CaseInc_LB_2.5%M+F", "CaseInc_UB_97.5%M+F","NationalPop15-49yF", "NationalPop15-49yM",
                       "NationalPop15-49yM+F", "Country_Curve_Fit", "Date_Last_Run", "Adult prev. USED -- Spectrum or other",
                       "LB, Adult prev. USED -- Spectrum or other", "Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)",
                       "LB on maternal prevalence, incl. Imputed", "UB on maternal prevalence, incl. Imputed",
                       "Syphilis-infected pregnancies", "Treated mothers", "Untreated mothers", "ABO cases (2012 method)",
                       "CS cases (2012 method)", "ABO, not seen in ANC (2012 method)", "ABO, ANC women not screened (2012 method)",
                       "ABO, ANC-screened women not treated (2012 method)", "ABO = CS, ANC women treated (2012 method)",
                       "CS cases, not seen in ANC (2012 method)", "CS cases, ANC women not screened (2012 method)",
                       "CS cases, ANC-screened women not treated (2012 method)", "ABO cases", "CS cases",
                       "ABO, not seen in ANC", "ABO, ANC women not screened", "ABO, ANC-screened women not treated",
                       "ABO = CS, ANC women treated", "CS cases, not seen in ANC", "CS cases, ANC women not screened",
                       "CS cases, ANC-screened women not treated", "CS / 100,000 live births", "CS / 100,000, 2012, for RANK",
                       "CS / 100,000, 2016, for RANK",  "CS report completeness", "ABO/100,000 live births", "ABO risk, treated mothers",
                       "Pregnancies, 2008", "Pregnancies, 2012", "Pregnancies, 2016", "Pregnancies, 2020", "ANC-1, imputed?",
                       "Screen cov, imputed?", "Maternal prevalence, imputed?", "3 ANC service coverages national?",
                       "Spectrum trend & national 3 ANC coverages Y=1", "Liveborn with clinical CS", "Prematurity or LBW due to CS",
                       "Neonatal death due to CS", "Stillbirth due to CS", "Asymptomatic CS", "Asymptomatic, mother diagnosed",
                       "Asymptomatic, mother NOT diagnosed","Asymptomatic, mother NOT diagnosed", "EMTCT elimination certification country?",
                       "CS case report rate (incl. Imputations as 0.3 cases for reported 0s), Elimination countries",
                       "ABO = CS, ANC women treated, ELIM.CIES only", "CS cases, not seen in ANC, ELIM.CIESonly",
                       "CS cases, ANC women not screened, ELIM.CIESonly", "CS cases, ANC-screened women not treated, ELIM.CIESonly",
                       "Graph labels", "Maternal duration of infection", "Maternal incidence/person-year",
                       "Additional maternal infections, from reinfection after treatment", "Relative increase in CS, from reinfection",
                       "Variance on adult (not only ANC) prevalence","Variance on MATERNAL prevalence",
                       "Variance on maternal prevalence * Pregnancies^2", "Variance on # pregnancies", "Variance on ABO prob. For treated mothers",
                       "LB on ANC-1 coverage", "UB on ANC-1 coverage", "Variance on ANC-1 coverage", "LB on screen coverage",
                       "UB on screen coverage", "Variance on Screen coverage", "LB on Treat coverage", "UB on Treat coverage",
                       "Variance on Treat coverage", "Variance on % of mothers treated, or untreated", "Variance on CS case number",
                       "Square of Variance on CS case number", "SE on CS case number", "Variance on ABO case number",
                       "Square of Variance on ABO case number", "SE on ABO case number", "2.5% of CS case number",
                       "97.5% of CS case number", "2.5% on ABO number","97.5% on ABO number", "SDG Region")

  base_variables_desciption <- c( "Rank, 2012 ABO cases", "Rank, 2016 ABO cases","Rank 2012, CS case rate",
                                  "Rank 2016, CS case rate", "Country","WHO Region","ISO code","ISO3nmb",
                                  "Year", "Live Births", "Still births", "Stillbirths, Blencowe & Hogan 2016",
                                  "Pregnancies", "Still/Live births", "Source of Live and/or Stillbirths",
                                  "Women with >= 1 ANC visit (%)", "Source of ANC1", "N tested (1st ANC visit)",
                                  "N, 1st visits", "Syphilis-tested (1st ANC, %)", "N tested (any visit ANC)",
                                  "N any ANC visits", "Syphilis-tested (any ANC visit,%)", "Source of Test coverage",
                                  "ANC women with syphilis, treated, N","Syphilis-infected ANC", "Treated (%)",
                                  "Source of Treated", "Congenital syphilis case report", "CS case report rate",
                                  "Estimate, Prev. F.", "Median Prev. F.", "LB, Prev. F. (%)", "UB, Prev. F. (%)",
                                  "Estimate, Prev. M.", "Median, Prev. M.", "LB, Prev. M. (%)", "UB, Prev. M. (%)",
                                  "Estimate, Prev. M.+F.", "Median, Prev. M.+F.", "LB Prev. M.+F. (%)",
                                  "UB, Prev. M.+F. (%)", "Estimate, Case Prev. F.", "Median, Case Prev. F.", "LB, Case Prev. F.",
                                  "UB, Case Prev. F.","Estimate, Case Prev. M.", "Median, Case Prev. M.", "LB, Case Prev. M.",
                                  "UB, Case Prev. M.", "Estimate, Case Prev. M.+F.", "Median, Case Prev. M.+F.",
                                  "LB, Case Prev. M.+F.", "UB, Case Prev. M.+F.", "Estimate, Inc. F.",
                                  "Median, Inc. F.", "LB, Inc. F.", "UB, Inc. F.", "Estimate, Inc. M.", "Median, Inc. M.", "LB, Inc. M.",
                                  "UB, Inc. M.","Estimate, Inc. M.+F.", "Median, Inc. M.+F.", "LB, Inc. M.+F.", "UB, Inc. M.+F.",
                                  "Estimate, Case Inc. F.", "Medidan, Case Inc. F.", "LB, Case Inc. F.", "UB, Case Inc. F.", "Estimate, Case Inc. M.",
                                  "Median, Case Inc. M.", "LB, Case Inc. M.", "UB, Case Inc. M.", "Estimate, Case Inc. M.+F.", "Median, Case Inc. M.+F.",
                                  "LB, Case Inc. M.+F.", "UB, Case Inc. M.+F.","National Pop. 15-49y F.", "National Pop. 15-49y M.",
                                  "National Pop. 15-49y M.+F.", "Country Curve Fit", "Date last run", "Adult prev. used -- Spectrum or other",
                                  "LB, Adult prev. Used -- Spectrum or other", "Maternal syphilis prevalence",
                                  "LB, maternal prev.", "UB, maternal prev.",
                                  "Syphilis-infected pregnancies", "Treated mothers", "Untreated mothers", "ABO cases (2012 method)",
                                  "CS cases (2012 method)", "ABO, not seen in ANC (2012 method)", "ABO, ANC women not screened (2012 method)",
                                  "ABO, ANC-screened women not treated (2012 method)", "ABO = CS, ANC women treated (2012 method)",
                                  "CS cases, not seen in ANC (2012 method)", "CS cases, ANC women not screened (2012 method)",
                                  "CS cases, ANC-screened women not treated (2012 method)", "ABO cases", "CS cases",
                                  "ABO, not seen in ANC", "ABO, ANC women not screened", "ABO, ANC-screened women not treated",
                                  "ABO = CS, ANC women treated", "CS cases, not seen in ANC", "CS cases, ANC women not screened",
                                  "CS cases, ANC-screened women not treated", "CS / 100,000 live births", "CS / 100,000, 2012, for rank",
                                  "CS / 100,000, 2016, for rank",  "CS report completeness", "ABO/100,000 live births", "ABO risk, treated mothers",
                                  "Pregnancies, 2008", "Pregnancies, 2012", "Pregnancies, 2016", "Pregnancies, 2020", "ANC-1, imputed (Yes/No)",
                                  "Screen cov., imputed (Yes/No)", "Maternal prev., imputed (Yes/No)", "3 ANC service coverages national (Yes/No)",
                                  "Spectrum trend & national 3 ANC coverages", "Liveborn with clinical CS", "Prematurity or LBW due to CS",
                                  "Neonatal death due to CS", "Stillbirth due to CS", "Asymptomatic CS", "Asymptomatic, mother diagnosed",
                                  "Asymptomatic, mother not diagnosed","Asymptomatic, mother not diagnosed", "EMTCT elimination certification country?",
                                  "CS case report rate (incl. Imputations as 0.3 cases for reported 0s), Elimination countries",
                                  "ABO = CS, ANC women treated, ELIM.CIES only", "CS cases, not seen in ANC, ELIM.CIES only",
                                  "CS cases, ANC women not screened, ELIM.CIES only", "CS cases, ANC-screened women not treated, ELIM.CIES only",
                                  "Graph labels", "Maternal duration of infection", "Maternal incidence/person-year",
                                  "Additional maternal infections, from reinfection after treatment", "Relative increase in CS, from reinfection",
                                  "Variance on adult (not only ANC) prevalence","Variance on maternal prevalence",
                                  "Variance on maternal prevalence * Pregnancies^2", "Variance on # pregnancies", "Variance on ABO prob. For treated mothers",
                                  "LB on ANC-1 coverage", "UB on ANC-1 coverage", "Variance on ANC-1 coverage", "LB on screen coverage",
                                  "UB on screen coverage", "Variance on Screen coverage", "LB on Treat coverage", "UB on Treat coverage",
                                  "Variance on Treat coverage", "Variance on % of mothers treated, or untreated", "Variance on CS case number",
                                  "Square of Variance on CS case number", "SE on CS case number", "Variance on ABO case number",
                                  "Square of Variance on ABO case number", "SE on ABO case number", "LB, CS case number",
                                  "UB, CS case number", "LB, ABO number","UB, ABO number", "SDG Region")

  var_name_table <- data.frame(dbname=base_variables,chname=base_variables_desciption)
  fixidx <- c(1:9);
  idxtoplot <- list(c('Pregnancies','Pregnancies','Pregnancies'),
                    c("Treated (%)","Treated (%)","Treated (%)"),
                    c("Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)","LB on maternal prevalence, incl. Imputed","UB on maternal prevalence, incl. Imputed"),
                    c("ABO cases","2.5% on ABO number","97.5% on ABO number"),
                    c("CS cases","2.5% of CS case number","97.5% of CS case number"),
                    c("CS / 100,000 live births","2.5% on CS / 100,000 live births","97.5% on CS / 100,000 live births"),
                    c("ABO/100,000 live births","2.5% on ABO/100,000 live births","97.5% on ABO/100,000 live births"),
                    c("Syphilis-tested (1st ANC, %)","Syphilis-tested (1st ANC, %)","Syphilis-tested (1st ANC, %)"),
                    c("Women with >= 1 ANC visit (%)","Women with >= 1 ANC visit (%)","Women with >= 1 ANC visit (%)"))

  BaseData <- CongenDataOut[,which(is.element(names(CongenDataOut),var_name_table$dbname[fixidx]))]
  LongCongenDataOutForPlots <- data.frame()
  LongCongenDataOutForPlotsRaw <- data.frame()

  if(nrow(BaseData)>=1)
  {
    for(ii in seq_len(length(idxtoplot)))
    {
      v_name <- idxtoplot[[ii]][1];
      idxb <- which(names(CongenDataOut)==v_name)
      idxlb <- which(names(CongenDataOut)==idxtoplot[[ii]][2])
      idxub <- which(names(CongenDataOut)==idxtoplot[[ii]][3])

      temp_data <- BaseData;
      temp_data$indicator <- v_name;
      temp_data$value <- CongenDataOut[,idxb];
      temp_data$lower <- CongenDataOut[,idxlb];
      temp_data$upper <- CongenDataOut[,idxub];
      temp_data$datatype <- "Projected"
      LongCongenDataOutForPlots <- rbind(LongCongenDataOutForPlots,temp_data)
    }

    #Prevalence
    temp_data <- BaseData[1:nrow(SyphDataRaw),];
    temp_data[,] <- NA
    temp_data$Country <- SyphDataRaw$Country
    temp_data$`WHO Region` <- sapply(SyphDataRaw$WHO_region, function(xx) substr(xx,1,nchar(xx)-1))
    temp_data$`ISO code` <- SyphDataRaw$ISO3_letters
    temp_data$ISO3nmb <- SyphDataRaw$ISO3
    temp_data$Year <- SyphDataRaw$Year
    temp_data$indicator <- "Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)";
    temp_data$value <- SyphDataRaw$BestPrevalence;
    temp_data$lower <- SyphDataRaw$LowerPrevalence;
    temp_data$upper <- SyphDataRaw$UpperPrevalence;
    temp_data$datatype <- SyphDataRaw$Data_type
    LongCongenDataOutForPlotsRaw <- temp_data
    rm(temp_data)

  }

  #"Syphilis-tested (1st ANC, %)"
  temp_data <- CongenDataRaw[,1:9];
  temp_data$indicator <- "Syphilis-tested (1st ANC, %)";
  temp_data$value <- CongenDataRaw$`Syphilis-tested (1st ANC, %)`;
  temp_data$lower <- CongenDataRaw$`Syphilis-tested (1st ANC, %)`;
  temp_data$upper <- CongenDataRaw$`Syphilis-tested (1st ANC, %)`;
  if(nrow(temp_data)>=1) temp_data$datatype <- "Reported"
  if(length(names(LongCongenDataOutForPlotsRaw))==length(names(temp_data))) LongCongenDataOutForPlotsRaw <- rbind(LongCongenDataOutForPlotsRaw,temp_data)
  rm(temp_data)

  #"Women with >= 1 ANC visit (%)"
  temp_data <- CongenDataRaw[,1:9];
  if(nrow(temp_data)>=1)
  {
    temp_data$indicator <- "Women with >= 1 ANC visit (%)";
    temp_data$value <- CongenDataRaw$`Women with >= 1 ANC visit (%)`;
    temp_data$lower <- CongenDataRaw$`Women with >= 1 ANC visit (%)`;
    temp_data$upper <- CongenDataRaw$`Women with >= 1 ANC visit (%)`;
    temp_data$datatype <- "Reported"
    LongCongenDataOutForPlotsRaw <- rbind(LongCongenDataOutForPlotsRaw,temp_data)
  }

  LongCongenDataOutForPlotsRaw$SDG_Region <- sapply(LongCongenDataOutForPlotsRaw$`ISO code`, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  LongCongenDataOutForPlotsRaw <- subset(LongCongenDataOutForPlotsRaw,!is.na(SDG_Region))
  rm(temp_data)

  #Grabing data by regions
  if(nrow(  LongCongenDataOutForPlots)>=1)
  {
    LongCongenDataOutForPlots <- subset(LongCongenDataOutForPlots,!(Country%in%paste(unique(LongCongenDataOutForPlots$`WHO Region`),"total")))
    LongCongenDataOutForPlots$SDG_Region <- sapply(LongCongenDataOutForPlots$`ISO code`, function(x)
    {
      res <- NA
      idx <- which(SDGRegions$COUNTRY_CODE==x)
      if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
      res
    })
  }

  if(is.null(list_countries))
  {
    #SDG Regions
    for(reg in unique(LongCongenDataOutForPlots$SDG_Region))#for(reg in unique(LongCongenDataOutForPlots$`WHO Region`))
    {
      years <- unique(LongCongenDataOutForPlots$Year)
      country <- reg;
      rk <- NA
      isoCode <- paste(reg,"r",sep="_")
      temp_data <- subset(LongCongenDataOutForPlots,LongCongenDataOutForPlots$SDG_Region==reg)

      femalepop <- sapply(seq_len(nrow(temp_data)), function(ii){
        ctr <- temp_data$`ISO code`[ii];
        yy <- temp_data$Year[ii]
        idx <- which(CongenDataOut$`ISO code`==ctr & CongenDataOut$Year==yy)
        res <- NA
        if(length(idx)>=1)
        {
          vals <- CongenDataOut$`NationalPop15-49yF`[idx]
          if(any(!is.na(vals))) res <- mean(vals,na.rm=T)
        }
        res
      })

      pregnancies <- sapply(seq_len(nrow(temp_data)), function(ii){
        ctr <- temp_data$`ISO code`[ii];
        yy <- temp_data$Year[ii]
        idx <- which(CongenDataOut$`ISO code`==ctr & CongenDataOut$Year==yy)
        res <- NA
        if(length(idx)>=1)
        {
          vals <- CongenDataOut$Pregnancies[idx]
          if(any(!is.na(vals))) res <- mean(vals,na.rm=T)
        }
        res
      })

      livebirths <- sapply(seq_len(nrow(temp_data)), function(ii){
        ctr <- temp_data$`ISO code`[ii];
        yy <- temp_data$Year[ii]
        idx <- which(CongenDataOut$`ISO code`==ctr & CongenDataOut$Year==yy)
        res <- NA
        if(length(idx)>=1)
        {
          vals <- CongenDataOut$`Live Births`[idx]
          if(any(!is.na(vals))) res <- mean(vals,na.rm=T)
        }
        res
      })

      for(indic in c("Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)",
                     "Treated (%)","CS / 100,000 live births",
                     "ABO/100,000 live births", "Syphilis-tested (1st ANC, %)", "Women with >= 1 ANC visit (%)"))
      {
        value <- lower <- upper <- numeric()
        for(yii  in seq_len(length(years)))
        {
          #cat("--", reg,"; year=",years[yii], "\n" )
          value[yii] <- NA;
          lower[yii] <- NA;
          upper[yii] <- NA;

          av <- temp_data$value[temp_data$indicator==indic & temp_data$Year==years[yii]]
          indic_weights <- pregnancies[temp_data$indicator==indic & temp_data$Year==years[yii]];
          if(!indic%in%c("Treated (%)", "Syphilis-tested (1st ANC, %)", "Women with >= 1 ANC visit (%)")) indic_weights <- livebirths[temp_data$indicator==indic & temp_data$Year==years[yii]]

          if(!all(is.na(av)))
          {
            value[yii] <- weighted.mean(av,indic_weights,na.rm=T)
            vsd <- ((temp_data$upper[temp_data$indicator==indic & temp_data$Year==years[yii]]-temp_data$lower[temp_data$indicator==indic & temp_data$Year==years[yii]])/2/1.96)^2
            vsd <- sqrt(sum(vsd*indic_weights^2 ,na.rm=T)/sum(indic_weights,na.rm=T)^2)
            lower[yii] <- value[yii]-1.96*vsd;
            upper[yii] <- value[yii]+1.96*vsd;
          }#End if(!all(is.na(av)))
        }

        lower[lower<0] = 0;
        tt_dattomerge <- temp_data[1:length(years),]
        tt_dattomerge$'Rank, 2012 ABO cases' <- tt_dattomerge$'Rank, 2016 ABO cases' <- tt_dattomerge$'Rank 2012, CS case RATE' <-rk
        tt_dattomerge$'Rank 2016, CS case RATE' <- rk
        tt_dattomerge$Country <- country
        tt_dattomerge$'WHO Region' <- NA#reg #isoCode
        tt_dattomerge$'ISO code' <- isoCode
        tt_dattomerge$'ISO3nmb' <- isoCode
        tt_dattomerge$Year <- years
        tt_dattomerge$indicator <- indic
        tt_dattomerge$value <- value
        tt_dattomerge$upper <- upper
        tt_dattomerge$lower <- lower
        tt_dattomerge$datatype <- "Projected"
        tt_dattomerge$SDG_Region <- reg #isoCode
        LongCongenDataOutForPlots <- rbind(LongCongenDataOutForPlots,tt_dattomerge)
      }#End

      ###
      for(indic in c("Pregnancies","ABO cases",
                     "CS cases"))
      {
        value <- lower <- upper <- numeric()
        for(yii  in seq_len(length(years)))
        {
          value[yii] <- NA;
          lower[yii] <- NA;
          upper[yii] <- NA;

          av <- temp_data$value[temp_data$indicator==indic & temp_data$Year==years[yii]]
          if(!any(is.na(av)))
          {
            value[yii] <- sum(av)
            vsd <- ((temp_data$upper[temp_data$indicator==indic & temp_data$Year==years[yii]]-temp_data$lower[temp_data$indicator==indic & temp_data$Year==years[yii]])/2/1.96)^2
            vsd <- sqrt(sum(vsd, na.rm=T))
            lower[yii] <- value[yii]-1.96*vsd;
            upper[yii] <- value[yii]+1.96*vsd;
          }#End if(!all(is.na(av)))
        }

        lower[lower<0] = 0;
        tt_dattomerge <- temp_data[1:length(years),]
        tt_dattomerge$'Rank, 2012 ABO cases' <- tt_dattomerge$'Rank, 2016 ABO cases' <- tt_dattomerge$'Rank 2012, CS case RATE' <-rk
        tt_dattomerge$'Rank 2016, CS case RATE' <- rk
        tt_dattomerge$Country <- country
        tt_dattomerge$'WHO Region' <- NA #isoCode
        tt_dattomerge$'ISO code' <- isoCode
        tt_dattomerge$'ISO3nmb' <- isoCode
        tt_dattomerge$Year <- years
        tt_dattomerge$indicator <- indic
        tt_dattomerge$value <- value
        tt_dattomerge$upper <- upper
        tt_dattomerge$lower <- lower
        tt_dattomerge$datatype <- "Projected"
        tt_dattomerge$SDG_Region <- reg #isoCode
        LongCongenDataOutForPlots <- rbind(LongCongenDataOutForPlots,tt_dattomerge)
      }#End
    }

    ###Global
    years <- unique(LongCongenDataOutForPlots$Year)
    country <- "Global";
    rk <- NA
    isoCode <- paste(reg,"r",sep="_")

    ###
    temp_data <- subset(LongCongenDataOutForPlots, !(Country%in%unique(LongCongenDataOutForPlots$SDG_Region)))
    for(indic in c("Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)",
                   "Treated (%)","CS / 100,000 live births",
                   "ABO/100,000 live births", "Syphilis-tested (1st ANC, %)", "Women with >= 1 ANC visit (%)"))
    {
      temp_data <- subset(LongCongenDataOutForPlots, !(Country%in%unique(LongCongenDataOutForPlots$SDG_Region)))#LongCongenDataOutForPlots

      femalepop <- sapply(seq_len(nrow(temp_data)), function(ii){
        ctr <- temp_data$`ISO code`[ii];
        yy <- temp_data$Year[ii]
        idx <- which(CongenDataOut$`ISO code`==ctr & CongenDataOut$Year==yy)
        res <- NA
        if(length(idx)>=1)
        {
          vals <- CongenDataOut$`NationalPop15-49yF`[idx]
          if(any(!is.na(vals))) res <- mean(vals,na.rm=T)
        }
        res
      })

      pregnancies <- sapply(seq_len(nrow(temp_data)), function(ii){
        ctr <- temp_data$`ISO code`[ii];
        yy <- temp_data$Year[ii]
        idx <- which(CongenDataOut$`ISO code`==ctr & CongenDataOut$Year==yy)
        res <- NA
        if(length(idx)>=1)
        {
          vals <- CongenDataOut$Pregnancies[idx]
          if(any(!is.na(vals))) res <- mean(vals,na.rm=T)
        }
        res
      })

      livebirths <- sapply(seq_len(nrow(temp_data)), function(ii){
        ctr <- temp_data$`ISO code`[ii];
        yy <- temp_data$Year[ii]
        idx <- which(CongenDataOut$`ISO code`==ctr & CongenDataOut$Year==yy)
        res <- NA
        if(length(idx)>=1)
        {
          vals <- CongenDataOut$`Live Births`[idx]
          if(any(!is.na(vals))) res <- mean(vals,na.rm=T)
        }
        res
      })

      value <- lower <- upper <- numeric()
      for(yii  in seq_len(length(years)))
      {
        value[yii] <- NA;
        lower[yii] <- NA;
        upper[yii] <- NA;

        av <- temp_data$value[temp_data$indicator==indic & temp_data$Year==years[yii]]
        indic_weights <- pregnancies[temp_data$indicator==indic & temp_data$Year==years[yii]];
        if(!indic%in%c("Treated (%)", "Syphilis-tested (1st ANC, %)", "Women with >= 1 ANC visit (%)")) indic_weights <- livebirths[temp_data$indicator==indic & temp_data$Year==years[yii]]

        if(!all(is.na(av)))
        {
          value[yii] <- weighted.mean(av,indic_weights,na.rm=T)
          vsd <- ((temp_data$upper[temp_data$indicator==indic & temp_data$Year==years[yii]]-temp_data$lower[temp_data$indicator==indic & temp_data$Year==years[yii]])/2/1.96)^2
          vsd <- sqrt(sum(vsd*indic_weights^2,na.rm=T)/sum(indic_weights,na.rm=T)^2)
          lower[yii] <- value[yii]-1.96*vsd;
          upper[yii] <- value[yii]+1.96*vsd;
        }#End if(!all(is.na(av)))
      }

      lower[lower<0] = 0;
      tt_dattomerge <- temp_data[1:length(years),]
      tt_dattomerge$'Rank, 2012 ABO cases' <- tt_dattomerge$'Rank, 2016 ABO cases' <- tt_dattomerge$'Rank 2012, CS case RATE' <-rk
      tt_dattomerge$'Rank 2016, CS case RATE' <- rk
      tt_dattomerge$Country <- "Global"
      tt_dattomerge$'WHO Region' <- "Global" #isoCode
      tt_dattomerge$'ISO code' <- "Global"
      tt_dattomerge$'ISO3nmb' <- "Global"
      tt_dattomerge$Year <- years
      tt_dattomerge$indicator <- indic
      tt_dattomerge$value <- value
      tt_dattomerge$upper <- upper
      tt_dattomerge$lower <- lower
      tt_dattomerge$datatype <- "Projected"
      tt_dattomerge$SDG_Region <- "Global" #isoCode
      LongCongenDataOutForPlots <- rbind(LongCongenDataOutForPlots,tt_dattomerge)
    }#End

    ###
    for(indic in c("Pregnancies","ABO cases",
                   "CS cases"))
    {
      value <- lower <- upper <- numeric()
      for(yii  in seq_len(length(years)))
      {
        value[yii] <- NA;
        lower[yii] <- NA;
        upper[yii] <- NA;
        av <- temp_data$value[temp_data$indicator==indic & temp_data$Year==years[yii]]
        if(!any(is.na(av)))
        {
          value[yii] <- sum(av)
          vsd <- ((temp_data$upper[temp_data$indicator==indic & temp_data$Year==years[yii]]-temp_data$lower[temp_data$indicator==indic & temp_data$Year==years[yii]])/2/1.96)^2
          vsd <- sqrt(sum(vsd, na.rm=T))
          lower[yii] <- value[yii]-1.96*vsd;
          upper[yii] <- value[yii]+1.96*vsd;
        }#End if(!all(is.na(av)))
      }

      lower[lower<0] = 0;
      tt_dattomerge <- temp_data[1:length(years),]
      tt_dattomerge$'Rank, 2012 ABO cases' <- tt_dattomerge$'Rank, 2016 ABO cases' <- tt_dattomerge$'Rank 2012, CS case RATE' <-rk
      tt_dattomerge$'Rank 2016, CS case RATE' <- rk
      tt_dattomerge$Country <- "Global"
      tt_dattomerge$'WHO Region' <- "Global" #isoCode
      tt_dattomerge$'ISO code' <- "Global"
      tt_dattomerge$'ISO3nmb' <- "Global"
      tt_dattomerge$Year <- years
      tt_dattomerge$indicator <- indic
      tt_dattomerge$value <- value
      tt_dattomerge$upper <- upper
      tt_dattomerge$lower <- lower
      tt_dattomerge$datatype <- "Projected"
      tt_dattomerge$SDG_Region <- "Global" #isoCode
      LongCongenDataOutForPlots <- rbind(LongCongenDataOutForPlots,tt_dattomerge)
    }#End
  }
  LongCongenDataOutForPlots <- rbind(LongCongenDataOutForPlots,LongCongenDataOutForPlotsRaw)
  rm(LongCongenDataOutForPlotsRaw)
  ###############################################################################
  ###############################################################################
  ###############################################################################
  RegCSABO <- data.frame()
  if(is.null(list_countries))
  {
    for(reg in unique(CongenDataOut$SDG_Region))
    {
      all_reg_data <- subset(CongenDataOut,CongenDataOut$SDG_Region==reg)
      for(yy in unique(CongenDataOut$Year))
      {
        all_reg_data_yy <- subset(all_reg_data,Year==yy)
        livebirths <- sum(all_reg_data_yy$`Live Births`,na.rm=T)
        Livebornwithcliniccs <- sum(all_reg_data_yy$`Liveborn with clinical CS`,na.rm=T)/livebirths*100000
        PrematurityLBWduecs <- sum(all_reg_data_yy$`Prematurity or LBW due to CS`,na.rm=T)/livebirths*100000
        Neonataldeathduecs <- sum(all_reg_data_yy$`Stillbirth due to CS`,na.rm=T)/livebirths*100000
        StillbirthduetoCS <- sum(all_reg_data_yy$`Stillbirth due to CS`,na.rm=T)/livebirths*100000
        Asymptomaticcs <- sum(all_reg_data_yy$`Neonatal death due to CS`,na.rm=T)/livebirths*100000

        ABOnotseeninANC <- sum(all_reg_data_yy$`ABO, not seen in ANC`, na.rm=T)
        ABOANCwomennotscreened <- sum(all_reg_data_yy$`ABO, ANC women not screened`, na.rm=T)
        ABOANCscreenedwomennottreated <- sum(all_reg_data_yy$`ABO, ANC-screened women not treated`, na.rm=T)
        ABOANCtestedAndtreated <- sum(all_reg_data_yy$`ABO = CS, ANC women treated`, na.rm=T)

        temp <- data.matrix(matrix(c(Livebornwithcliniccs,PrematurityLBWduecs, Neonataldeathduecs, StillbirthduetoCS, Asymptomaticcs,livebirths,
                                     ABOnotseeninANC, ABOANCwomennotscreened, ABOANCscreenedwomennottreated, ABOANCtestedAndtreated), nrow=1))

        colnames(temp) <- c("Liveborn with clinical CS",
                            "Prematurity or LBW due to CS",
                            "Neonatal death due to CS",
                            "Stillbirth due to CS",
                            "Asymptomatic CS",
                            "Live births",
                            "ABO, not seen in ANC",
                            "ABO, ANC women not screened",
                            "ABO, ANC-screened women not treated",
                            "ABO, ANC women treated")

        temp <- as.data.frame(temp)

        alltemp <- data.frame(Country=NA,'WHO Region'=NA,'ISO code'=NA, SDG_Region=reg, Year=yy, check.names = F)
        alltemp <- cbind(alltemp,temp)
        RegCSABO <- rbind(RegCSABO,alltemp)
      }#End for(yy in unique(CongenDataOut$Year))
    }#for(reg in unique(CongenDataOut$`WHO Region`))

    #for(reg in unique(CongenDataOut$`WHO Region`))
    {
      reg <- "Global"
      all_reg_data <- CongenDataOut
      for(yy in unique(CongenDataOut$Year))
      {
        all_reg_data_yy <- subset(all_reg_data,Year==yy)
        livebirths <- sum(all_reg_data_yy$`Live Births`,na.rm=T)
        Livebornwithcliniccs <- sum(all_reg_data_yy$`Liveborn with clinical CS`,na.rm=T)/livebirths*100000
        PrematurityLBWduecs <- sum(all_reg_data_yy$`Prematurity or LBW due to CS`,na.rm=T)/livebirths*100000
        Neonataldeathduecs <- sum(all_reg_data_yy$`Stillbirth due to CS`,na.rm=T)/livebirths*100000
        StillbirthduetoCS <- sum(all_reg_data_yy$`Stillbirth due to CS`,na.rm=T)/livebirths*100000
        Asymptomaticcs <- sum(all_reg_data_yy$`Asymptomatic CS`,na.rm=T)/livebirths*100000

        ABOnotseeninANC <- sum(all_reg_data_yy$`ABO, not seen in ANC`, na.rm=T)
        ABOANCwomennotscreened <- sum(all_reg_data_yy$`ABO, ANC women not screened`, na.rm=T)
        ABOANCscreenedwomennottreated <- sum(all_reg_data_yy$`ABO, ANC-screened women not treated`, na.rm=T)
        ABOANCtestedAndtreated <- sum(all_reg_data_yy$`ABO = CS, ANC women treated`, na.rm=T)

        temp <- data.matrix(matrix(c(Livebornwithcliniccs,PrematurityLBWduecs, Neonataldeathduecs, StillbirthduetoCS, Asymptomaticcs,livebirths,
                                     ABOnotseeninANC, ABOANCwomennotscreened, ABOANCscreenedwomennottreated, ABOANCtestedAndtreated), nrow=1))

        colnames(temp) <- c("Liveborn with clinical CS",
                            "Prematurity or LBW due to CS",
                            "Neonatal death due to CS",
                            "Stillbirth due to CS",
                            "Asymptomatic CS",
                            "Live births",
                            "ABO, not seen in ANC",
                            "ABO, ANC women not screened",
                            "ABO, ANC-screened women not treated",
                            "ABO, ANC women treated")

        temp <- as.data.frame(temp)

        alltemp <- data.frame(Country=NA,'WHO Region'=NA,'ISO code'=NA, SDG_Region=reg,Year=yy, check.names = F)
        alltemp <- cbind(alltemp,temp)
        RegCSABO <- rbind(RegCSABO,alltemp)
      }#End for(yy in unique(CongenDataOut$Year))
    }#for(reg in unique(CongenDataOut$`WHO Region`))
  }

  LongRegCSABO <- data.frame()
  if(is.null(list_countries))
  {
    for(vn in c("Liveborn with clinical CS",
                "Prematurity or LBW due to CS",
                "Neonatal death due to CS",
                "Stillbirth due to CS",
                "Asymptomatic CS",
                "Live births",
                "ABO, not seen in ANC",
                "ABO, ANC women not screened",
                "ABO, ANC-screened women not treated",
                "ABO, ANC women treated"))
    {
      temp_data <- RegCSABO[,names(RegCSABO)%in%c("Country", "WHO Region", "ISO code", "SDG_Region", "Year", vn)]
      names(temp_data)[names(temp_data)==vn] <- "value"
      temp_data$indicator <- vn
      LongRegCSABO <- rbind(LongRegCSABO,temp_data)
    }
  }
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ################################################################################
  ###Global and regional Syphilis Prevalence estimates
  previnc_vars <- c("Country","WHO Region","ISO code", "Year", "EstimatePrevF", "MedianPrevF",
                    "PrevLB_2.5%F", "PrevUB_97.5%F", "EstimatePrevM",
                    "MedianPrevM", "PrevLB_2.5%M", "PrevUB_97.5%M", "EstimatePrevM+F","MedianPrevM+F",
                    "PrevLB_2.5%M+F", "PrevUB_97.5%M+F", "CasePrev_EstF", "CasePrev_MedF", "CasePrev_LB_2.5%F",
                    "CasePrev_UB_97.5%F","CasePrev_EstM", "CasePrev_MedM", "CasePrev_LB_2.5%M", "CasePrev_UB_97.5%M",
                    "CasePrev_EstM+F", "CasePrev_MedM+F", "CasePrev_LB_2.5%M+F", "CasePrev_UB_97.5%M+F", "EstimateIncF",
                    "MedianIncF", "IncLB_2.5%F", "IncUB_97.5%F", "EstimateIncM", "MedianIncM", "IncLB_2.5%M",
                    "IncUB_97.5%M","EstimateIncM+F", "MedianIncM+F", "IncLB_2.5%M+F", "IncUB_97.5%M+F",
                    "CaseInc_EstF", "CaseInc_MedF", "CaseInc_LB_2.5%F", "CaseInc_UB_97.5%F", "CaseInc_EstM",
                    "CaseInc_MedM", "CaseInc_LB_2.5%M", "CaseInc_UB_97.5%M", "CaseInc_EstM+F", "CaseInc_MedM+F",
                    "CaseInc_LB_2.5%M+F", "CaseInc_UB_97.5%M+F","NationalPop15-49yF", "NationalPop15-49yM",
                    "NationalPop15-49yM+F", "Country_Curve_Fit")

  RegDataPrevInc <- data.frame()
  if(is.null(list_countries))
  {
    for(reg in unique(CongenDataOut$SDG_Region))
    {
      all_reg_data <- subset(CongenDataOut,CongenDataOut$SDG_Region==reg)
      for(yy in unique(CongenDataOut$Year))
      {
        all_reg_data_yy <- subset(all_reg_data,Year==yy)
        EstimatePrevF <- weighted.mean(all_reg_data_yy$'EstimatePrevF',dplyr::coalesce(all_reg_data_yy$'NationalPop15-49yF',0),na.rm=T)
        MedianPrevF <- median(all_reg_data_yy$'EstimatePrevF',na.rm = T)

        vprevF <- ((all_reg_data_yy$'PrevUB_97.5%F'-all_reg_data_yy$'PrevLB_2.5%F')/1.96/2)^2
        sprevF <- sqrt(sum(vprevF*all_reg_data_yy$'NationalPop15-49yF'^2,na.rm=T)/sum(all_reg_data_yy$'NationalPop15-49yF',na.rm=T)^2)

        PrevLB_F <- max(EstimatePrevF-1.96*sprevF,0)
        PrevUB_F <- min(EstimatePrevF+1.96*sprevF,1)

        #
        EstimatePrevM <- weighted.mean(all_reg_data_yy$'EstimatePrevM',dplyr::coalesce(all_reg_data_yy$'NationalPop15-49yM',0),na.rm=T)
        MedianPrevM <- median(all_reg_data_yy$'EstimatePrevM',na.rm = T)

        vprevM <- ((all_reg_data_yy$'PrevUB_97.5%M'-all_reg_data_yy$'PrevLB_2.5%M')/1.96/2)^2
        sprevM <- sqrt(sum(vprevM*all_reg_data_yy$'NationalPop15-49yM'^2,na.rm=T)/sum(all_reg_data_yy$'NationalPop15-49yM',na.rm=T)^2)

        PrevLB_M <- max(EstimatePrevM-1.96*sprevM,0)
        PrevUB_M <- min(EstimatePrevM+1.96*sprevM,1)

        #
        EstimatePrevMF <- weighted.mean(all_reg_data_yy$'EstimatePrevM+F',dplyr::coalesce(all_reg_data_yy$'NationalPop15-49yM+F',0),na.rm=T)
        MedianPrevMF <- median(all_reg_data_yy$'EstimatePrevM+F',na.rm = T)

        vprevMF <- ((all_reg_data_yy$'PrevUB_97.5%M+F'-all_reg_data_yy$'PrevLB_2.5%M+F')/1.96/2)^2
        sprevMF <- sqrt(sum(vprevMF*all_reg_data_yy$'NationalPop15-49yM+F'^2,na.rm=T)/sum(all_reg_data_yy$'NationalPop15-49yM+F',na.rm=T)^2)

        PrevLB_MF <- max(EstimatePrevMF-1.96*sprevMF,0)
        PrevUB_MF <- min(EstimatePrevMF+1.96*sprevM,1)

        #
        CasePrev_EstF <- sum(all_reg_data_yy$CasePrev_EstF, na.rm=T)
        CasePrev_MedF <- sum(all_reg_data_yy$CasePrev_MedF, na.rm=T)

        vprevF <- ((all_reg_data_yy$'CasePrev_UB_97.5%F'-all_reg_data_yy$'CasePrev_LB_2.5%F')/1.96/2)^2
        sprevF <- sqrt(sum(vprevF,na.rm=T))

        CasePrevLB_F <- max(CasePrev_EstF - 1.96*sprevF,0)
        CasePrevUB_F <- CasePrev_EstF + 1.96*sprevF

        #
        CasePrev_EstM <- sum(all_reg_data_yy$CasePrev_EstM, na.rm=T)
        CasePrev_MedM <- sum(all_reg_data_yy$CasePrev_MedM, na.rm=T)

        vprevM <- ((all_reg_data_yy$'CasePrev_UB_97.5%M'-all_reg_data_yy$'CasePrev_LB_2.5%M')/1.96/2)^2
        sprevM <- sqrt(sum(vprevM,na.rm=T))

        CasePrevLB_M <- max(CasePrev_EstM - 1.96*sprevM,0)
        CasePrevUB_M <- CasePrev_EstM + 1.96*sprevM

        #
        CasePrev_EstMF <- sum(all_reg_data_yy$'CasePrev_EstM+F', na.rm=T)
        CasePrev_MedMF <- sum(all_reg_data_yy$'CasePrev_MedM+F', na.rm=T)

        vprevMF <- ((all_reg_data_yy$'CasePrev_UB_97.5%M+F'-all_reg_data_yy$'CasePrev_LB_2.5%M+F')/1.96/2)^2
        sprevMF <- sqrt(sum(vprevMF,na.rm=T))

        CasePrevLB_MF <- max(CasePrev_EstMF - 1.96*sprevMF,0)
        CasePrevUB_MF <- CasePrev_EstMF + 1.96*sprevMF

        #
        EstimateIncF <- weighted.mean(all_reg_data_yy$EstimateIncF,dplyr::coalesce(all_reg_data_yy$`NationalPop15-49yF`,0),na.rm=T)
        MedianIncF <- median(all_reg_data_yy$MedianIncF,na.rm = T)

        vprevF <- ((all_reg_data_yy$`IncUB_97.5%F`-all_reg_data_yy$`IncLB_2.5%F`)/1.96/2)^2
        sprevF <- sqrt(sum(vprevF*all_reg_data_yy$`NationalPop15-49yF`^2,na.rm=T)/sum(all_reg_data_yy$`NationalPop15-49yF`,na.rm=T)^2)

        IncLB_F <- max(EstimateIncF-1.96*sprevF,0)
        IncUB_F <- min(EstimateIncF+1.96*sprevF,1)

        #
        EstimateIncM <- weighted.mean(all_reg_data_yy$EstimateIncM,dplyr::coalesce(all_reg_data_yy$`NationalPop15-49yM`,0),na.rm=T)
        MedianIncM <- median(all_reg_data_yy$MedianIncM,na.rm = T)

        vprevM <- ((all_reg_data_yy$`IncUB_97.5%M`-all_reg_data_yy$`IncLB_2.5%M`)/1.96/2)^2
        sprevM <- sqrt(sum(vprevM*all_reg_data_yy$`NationalPop15-49yM`^2,na.rm=T)/sum(all_reg_data_yy$`NationalPop15-49yM`,na.rm=T)^2)

        IncLB_M <- max(EstimateIncM-1.96*sprevM,0)
        IncUB_M <- min(EstimateIncM+1.96*sprevM,1)

        #
        EstimateIncMF <- weighted.mean(all_reg_data_yy$'EstimateIncM+F',dplyr::coalesce(all_reg_data_yy$`NationalPop15-49yM+F`,0),na.rm=T)
        MedianIncMF <- median(all_reg_data_yy$'MedianIncM+F',na.rm = T)

        vprevMF <- ((all_reg_data_yy$`IncUB_97.5%M+F`-all_reg_data_yy$`IncLB_2.5%M+F`)/1.96/2)^2
        sprevMF <- sqrt(sum(vprevMF*all_reg_data_yy$`NationalPop15-49yM+F`^2,na.rm=T)/sum(all_reg_data_yy$`NationalPop15-49yM+F`,na.rm=T)^2)

        IncLB_MF <- max(EstimateIncMF-1.96*sprevMF,0)
        IncUB_MF <- min(EstimateIncMF+1.96*sprevMF,1)

        #
        CaseInc_EstF <- sum(all_reg_data_yy$CaseInc_EstF, na.rm=T)
        CaseInc_MedF <- sum(all_reg_data_yy$CaseInc_MedF, na.rm=T)

        vprevF <- ((all_reg_data_yy$`CaseInc_UB_97.5%F`-all_reg_data_yy$`CaseInc_LB_2.5%F`)/1.96/2)^2
        sprevF <- sqrt(sum(vprevF,na.rm=T))

        CaseInc_LB_F <- max(CaseInc_EstF - 1.96*sprevF,0)
        CaseInc_UB_F <- CaseInc_EstF + 1.96*sprevF

        #
        CaseInc_EstM <- sum(all_reg_data_yy$CaseInc_EstM, na.rm=T)
        CaseInc_MedM <- sum(all_reg_data_yy$CaseInc_MedM, na.rm=T)

        vprevM <- ((all_reg_data_yy$`CaseInc_UB_97.5%M`-all_reg_data_yy$`CaseInc_LB_2.5%M`)/1.96/2)^2
        sprevM <- sqrt(sum(vprevM,na.rm=T))

        CaseInc_LB_M <- max(CaseInc_EstM - 1.96*sprevM,0)
        CaseInc_UB_M <- CaseInc_EstM + 1.96*sprevM

        #
        CaseInc_EstMF <- sum(all_reg_data_yy$'CaseInc_EstM+F', na.rm=T)
        CaseInc_MedMF <- sum(all_reg_data_yy$'CaseInc_MedM+F', na.rm=T)

        vprevMF <- ((all_reg_data_yy$`CaseInc_UB_97.5%M+F`-all_reg_data_yy$`CaseInc_LB_2.5%M+F`)/1.96/2)^2
        sprevMF <- sqrt(sum(vprevMF,na.rm=T))

        CaseInc_LB_MF <- max(CaseInc_EstMF - 1.96*sprevMF,0)
        CaseInc_UB_MF <- CaseInc_EstMF + 1.96*sprevMF

        #
        NationalPop1549yF <- sum(all_reg_data_yy$'NationalPop15-49yF', na.rm=T)
        NationalPop1549yM <- sum(all_reg_data_yy$'NationalPop15-49yM', na.rm=T)
        NationalPop1549yMF <- sum(all_reg_data_yy$'NationalPop15-49yM+F', na.rm=T)

        temp <- data.matrix(matrix(c(EstimatePrevF, MedianPrevF, PrevLB_F, PrevUB_F, EstimatePrevM, MedianPrevM, PrevLB_M, PrevUB_M, EstimatePrevMF, MedianPrevMF,
                                     PrevLB_MF, PrevUB_MF, CasePrev_EstF, CasePrev_MedF, CasePrevLB_F, CasePrevUB_F, CasePrev_EstM, CasePrev_MedM, CasePrevLB_M,
                                     CasePrevUB_M, CasePrev_EstMF, CasePrev_MedMF, CasePrevLB_MF, CasePrevUB_MF, EstimateIncF, MedianIncF, IncLB_F, IncUB_F, EstimateIncM,
                                     MedianIncM, IncLB_M, IncUB_M, EstimateIncMF, MedianIncMF, IncLB_MF, IncUB_MF, CaseInc_EstF, CaseInc_MedF, CaseInc_LB_F, CaseInc_UB_F,
                                     CaseInc_EstM, CaseInc_MedM, CaseInc_LB_M, CaseInc_UB_M, CaseInc_EstMF, CaseInc_MedMF, CaseInc_LB_MF, CaseInc_UB_MF, NationalPop1549yF,
                                     NationalPop1549yM, NationalPop1549yMF), nrow=1))

        colnames(temp) <- c("EstimatePrevF", "MedianPrevF",
                            "PrevLB_2.5%F", "PrevUB_97.5%F", "EstimatePrevM",
                            "MedianPrevM", "PrevLB_2.5%M", "PrevUB_97.5%M", "EstimatePrevM+F","MedianPrevM+F",
                            "PrevLB_2.5%M+F", "PrevUB_97.5%M+F", "CasePrev_EstF", "CasePrev_MedF", "CasePrev_LB_2.5%F",
                            "CasePrev_UB_97.5%F","CasePrev_EstM", "CasePrev_MedM", "CasePrev_LB_2.5%M", "CasePrev_UB_97.5%M",
                            "CasePrev_EstM+F", "CasePrev_MedM+F", "CasePrev_LB_2.5%M+F", "CasePrev_UB_97.5%M+F", "EstimateIncF",
                            "MedianIncF", "IncLB_2.5%F", "IncUB_97.5%F", "EstimateIncM", "MedianIncM", "IncLB_2.5%M",
                            "IncUB_97.5%M","EstimateIncM+F", "MedianIncM+F", "IncLB_2.5%M+F", "IncUB_97.5%M+F",
                            "CaseInc_EstF", "CaseInc_MedF", "CaseInc_LB_2.5%F", "CaseInc_UB_97.5%F", "CaseInc_EstM",
                            "CaseInc_MedM", "CaseInc_LB_2.5%M", "CaseInc_UB_97.5%M", "CaseInc_EstM+F", "CaseInc_MedM+F",
                            "CaseInc_LB_2.5%M+F", "CaseInc_UB_97.5%M+F","NationalPop15-49yF", "NationalPop15-49yM",
                            "NationalPop15-49yM+F")

        temp <- as.data.frame(temp)
        alltemp <- data.frame(Country=NA,SDG_Region=reg,'ISO code'=NA, Year=yy, check.names = F)

        alltemp <- cbind(alltemp,temp)
        RegDataPrevInc <- rbind(RegDataPrevInc,alltemp)
      }
    }

    ##lobal
    #for(reg in unique(CongenDataOut$`WHO Region`))
    {
      reg <- "Global"
      all_reg_data <- CongenDataOut
      for(yy in unique(CongenDataOut$Year))
      {
        all_reg_data_yy <- subset(all_reg_data,Year==yy)
        EstimatePrevF <- weighted.mean(all_reg_data_yy$'EstimatePrevF',dplyr::coalesce(all_reg_data_yy$'NationalPop15-49yF',0),na.rm=T)
        MedianPrevF <- median(all_reg_data_yy$'EstimatePrevF',na.rm = T)

        vprevF <- ((all_reg_data_yy$'PrevUB_97.5%F'-all_reg_data_yy$'PrevLB_2.5%F')/1.96/2)^2
        sprevF <- sqrt(sum(vprevF*all_reg_data_yy$'NationalPop15-49yF'^2,na.rm=T)/sum(all_reg_data_yy$'NationalPop15-49yF',na.rm=T)^2)

        PrevLB_F <- max(EstimatePrevF-1.96*sprevF,0)
        PrevUB_F <- min(EstimatePrevF+1.96*sprevF,1)

        #
        EstimatePrevM <- weighted.mean(all_reg_data_yy$'EstimatePrevM',dplyr::coalesce(all_reg_data_yy$'NationalPop15-49yM',0),na.rm=T)
        MedianPrevM <- median(all_reg_data_yy$'EstimatePrevM',na.rm = T)

        vprevM <- ((all_reg_data_yy$'PrevUB_97.5%M'-all_reg_data_yy$'PrevLB_2.5%M')/1.96/2)^2
        sprevM <- sqrt(sum(vprevM*all_reg_data_yy$'NationalPop15-49yM'^2,na.rm=T)/sum(all_reg_data_yy$'NationalPop15-49yM',na.rm=T)^2)

        PrevLB_M <- max(EstimatePrevM-1.96*sprevM,0)
        PrevUB_M <- min(EstimatePrevM+1.96*sprevM,1)

        #
        EstimatePrevMF <- weighted.mean(all_reg_data_yy$'EstimatePrevM+F',dplyr::coalesce(all_reg_data_yy$'NationalPop15-49yM+F',0),na.rm=T)
        MedianPrevMF <- median(all_reg_data_yy$'EstimatePrevM+F',na.rm = T)

        vprevMF <- ((all_reg_data_yy$'PrevUB_97.5%M+F'-all_reg_data_yy$'PrevLB_2.5%M+F')/1.96/2)^2
        sprevMF <- sqrt(sum(vprevMF*all_reg_data_yy$'NationalPop15-49yM+F'^2,na.rm=T)/sum(all_reg_data_yy$'NationalPop15-49yM+F',na.rm=T)^2)

        PrevLB_MF <- max(EstimatePrevMF-1.96*sprevMF,0)
        PrevUB_MF <- min(EstimatePrevMF+1.96*sprevM,1)

        #
        CasePrev_EstF <- sum(all_reg_data_yy$CasePrev_EstF, na.rm=T)
        CasePrev_MedF <- sum(all_reg_data_yy$CasePrev_MedF, na.rm=T)

        vprevF <- ((all_reg_data_yy$'CasePrev_UB_97.5%F'-all_reg_data_yy$'CasePrev_LB_2.5%F')/1.96/2)^2
        sprevF <- sqrt(sum(vprevF,na.rm=T))

        CasePrevLB_F <- max(CasePrev_EstF - 1.96*sprevF,0)
        CasePrevUB_F <- CasePrev_EstF + 1.96*sprevF

        #
        CasePrev_EstM <- sum(all_reg_data_yy$CasePrev_EstM, na.rm=T)
        CasePrev_MedM <- sum(all_reg_data_yy$CasePrev_MedM, na.rm=T)

        vprevM <- ((all_reg_data_yy$'CasePrev_UB_97.5%M'-all_reg_data_yy$'CasePrev_LB_2.5%M')/1.96/2)^2
        sprevM <- sqrt(sum(vprevM,na.rm=T))

        CasePrevLB_M <- max(CasePrev_EstM - 1.96*sprevM,0)
        CasePrevUB_M <- CasePrev_EstM + 1.96*sprevM

        #
        CasePrev_EstMF <- sum(all_reg_data_yy$'CasePrev_EstM+F', na.rm=T)
        CasePrev_MedMF <- sum(all_reg_data_yy$'CasePrev_MedM+F', na.rm=T)

        vprevMF <- ((all_reg_data_yy$'CasePrev_UB_97.5%M+F'-all_reg_data_yy$'CasePrev_LB_2.5%M+F')/1.96/2)^2
        sprevMF <- sqrt(sum(vprevMF,na.rm=T))

        CasePrevLB_MF <- max(CasePrev_EstMF - 1.96*sprevMF,0)
        CasePrevUB_MF <- CasePrev_EstMF + 1.96*sprevMF

        #
        EstimateIncF <- weighted.mean(all_reg_data_yy$EstimateIncF,dplyr::coalesce(all_reg_data_yy$`NationalPop15-49yF`,0),na.rm=T)
        MedianIncF <- median(all_reg_data_yy$MedianIncF,na.rm = T)

        vprevF <- ((all_reg_data_yy$`IncUB_97.5%F`-all_reg_data_yy$`IncLB_2.5%F`)/1.96/2)^2
        sprevF <- sqrt(sum(vprevF*all_reg_data_yy$`NationalPop15-49yF`^2,na.rm=T)/sum(all_reg_data_yy$`NationalPop15-49yF`,na.rm=T)^2)

        IncLB_F <- max(EstimateIncF-1.96*sprevF,0)
        IncUB_F <- min(EstimateIncF+1.96*sprevF,1)

        #
        EstimateIncM <- weighted.mean(all_reg_data_yy$EstimateIncM,dplyr::coalesce(all_reg_data_yy$`NationalPop15-49yM`,0),na.rm=T)
        MedianIncM <- median(all_reg_data_yy$MedianIncM,na.rm = T)

        vprevM <- ((all_reg_data_yy$`IncUB_97.5%M`-all_reg_data_yy$`IncLB_2.5%M`)/1.96/2)^2
        sprevM <- sqrt(sum(vprevM*all_reg_data_yy$`NationalPop15-49yM`^2,na.rm=T)/sum(all_reg_data_yy$`NationalPop15-49yM`,na.rm=T)^2)

        IncLB_M <- max(EstimateIncM-1.96*sprevM,0)
        IncUB_M <- min(EstimateIncM+1.96*sprevM,1)

        #
        EstimateIncMF <- weighted.mean(all_reg_data_yy$'EstimateIncM+F',dplyr::coalesce(all_reg_data_yy$`NationalPop15-49yM+F`,0),na.rm=T)
        MedianIncMF <- median(all_reg_data_yy$'MedianIncM+F',na.rm = T)

        vprevMF <- ((all_reg_data_yy$`IncUB_97.5%M+F`-all_reg_data_yy$`IncLB_2.5%M+F`)/1.96/2)^2
        sprevMF <- sqrt(sum(vprevMF*all_reg_data_yy$`NationalPop15-49yM+F`^2,na.rm=T)/sum(all_reg_data_yy$`NationalPop15-49yM+F`,na.rm=T)^2)

        IncLB_MF <- max(EstimateIncMF-1.96*sprevMF,0)
        IncUB_MF <- min(EstimateIncMF+1.96*sprevMF,1)

        #
        CaseInc_EstF <- sum(all_reg_data_yy$CaseInc_EstF, na.rm=T)
        CaseInc_MedF <- sum(all_reg_data_yy$CaseInc_MedF, na.rm=T)

        vprevF <- ((all_reg_data_yy$`CaseInc_UB_97.5%F`-all_reg_data_yy$`CaseInc_LB_2.5%F`)/1.96/2)^2
        sprevF <- sqrt(sum(vprevF,na.rm=T))

        CaseInc_LB_F <- max(CaseInc_EstF - 1.96*sprevF,0)
        CaseInc_UB_F <- CaseInc_EstF + 1.96*sprevF

        #
        CaseInc_EstM <- sum(all_reg_data_yy$CaseInc_EstM, na.rm=T)
        CaseInc_MedM <- sum(all_reg_data_yy$CaseInc_MedM, na.rm=T)

        vprevM <- ((all_reg_data_yy$`CaseInc_UB_97.5%M`-all_reg_data_yy$`CaseInc_LB_2.5%M`)/1.96/2)^2
        sprevM <- sqrt(sum(vprevM,na.rm=T))

        CaseInc_LB_M <- max(CaseInc_EstM - 1.96*sprevM,0)
        CaseInc_UB_M <- CaseInc_EstM + 1.96*sprevM

        #
        CaseInc_EstMF <- sum(all_reg_data_yy$'CaseInc_EstM+F', na.rm=T)
        CaseInc_MedMF <- sum(all_reg_data_yy$'CaseInc_MedM+F', na.rm=T)

        vprevMF <- ((all_reg_data_yy$`CaseInc_UB_97.5%M+F`-all_reg_data_yy$`CaseInc_LB_2.5%M+F`)/1.96/2)^2
        sprevMF <- sqrt(sum(vprevMF,na.rm=T))

        CaseInc_LB_MF <- max(CaseInc_EstMF - 1.96*sprevMF,0)
        CaseInc_UB_MF <- CaseInc_EstMF + 1.96*sprevMF

        #
        NationalPop1549yF <- sum(all_reg_data_yy$'NationalPop15-49yF', na.rm=T)
        NationalPop1549yM <- sum(all_reg_data_yy$'NationalPop15-49yM', na.rm=T)
        NationalPop1549yMF <- sum(all_reg_data_yy$'NationalPop15-49yM+F', na.rm=T)

        temp <- data.matrix(matrix(c(EstimatePrevF, MedianPrevF, PrevLB_F, PrevUB_F, EstimatePrevM, MedianPrevM, PrevLB_M, PrevUB_M, EstimatePrevMF, MedianPrevMF,
                                     PrevLB_MF, PrevUB_MF, CasePrev_EstF, CasePrev_MedF, CasePrevLB_F, CasePrevUB_F, CasePrev_EstM, CasePrev_MedM, CasePrevLB_M,
                                     CasePrevUB_M, CasePrev_EstMF, CasePrev_MedMF, CasePrevLB_MF, CasePrevUB_MF, EstimateIncF, MedianIncF, IncLB_F, IncUB_F, EstimateIncM,
                                     MedianIncM, IncLB_M, IncUB_M, EstimateIncMF, MedianIncMF, IncLB_MF, IncUB_MF, CaseInc_EstF, CaseInc_MedF, CaseInc_LB_F, CaseInc_UB_F,
                                     CaseInc_EstM, CaseInc_MedM, CaseInc_LB_M, CaseInc_UB_M, CaseInc_EstMF, CaseInc_MedMF, CaseInc_LB_MF, CaseInc_UB_MF, NationalPop1549yF,
                                     NationalPop1549yM, NationalPop1549yMF), nrow=1))

        colnames(temp) <- c("EstimatePrevF", "MedianPrevF",
                            "PrevLB_2.5%F", "PrevUB_97.5%F", "EstimatePrevM",
                            "MedianPrevM", "PrevLB_2.5%M", "PrevUB_97.5%M", "EstimatePrevM+F","MedianPrevM+F",
                            "PrevLB_2.5%M+F", "PrevUB_97.5%M+F", "CasePrev_EstF", "CasePrev_MedF", "CasePrev_LB_2.5%F",
                            "CasePrev_UB_97.5%F","CasePrev_EstM", "CasePrev_MedM", "CasePrev_LB_2.5%M", "CasePrev_UB_97.5%M",
                            "CasePrev_EstM+F", "CasePrev_MedM+F", "CasePrev_LB_2.5%M+F", "CasePrev_UB_97.5%M+F", "EstimateIncF",
                            "MedianIncF", "IncLB_2.5%F", "IncUB_97.5%F", "EstimateIncM", "MedianIncM", "IncLB_2.5%M",
                            "IncUB_97.5%M","EstimateIncM+F", "MedianIncM+F", "IncLB_2.5%M+F", "IncUB_97.5%M+F",
                            "CaseInc_EstF", "CaseInc_MedF", "CaseInc_LB_2.5%F", "CaseInc_UB_97.5%F", "CaseInc_EstM",
                            "CaseInc_MedM", "CaseInc_LB_2.5%M", "CaseInc_UB_97.5%M", "CaseInc_EstM+F", "CaseInc_MedM+F",
                            "CaseInc_LB_2.5%M+F", "CaseInc_UB_97.5%M+F","NationalPop15-49yF", "NationalPop15-49yM",
                            "NationalPop15-49yM+F")

        temp <- as.data.frame(temp)
        alltemp <- data.frame(Country=NA,SDG_Region=reg,'ISO code'=NA, Year=yy, check.names = F)

        alltemp <- cbind(alltemp,temp)
        RegDataPrevInc <- rbind(RegDataPrevInc,alltemp)
      }
    }
  }#End if(!is.null(list_countries))

  ################################################################################
  ################################################################################
  var_names_for_long <- c("Prevalence (%)", "Prevalence (%)","Prevalence (%)",
                          "Prevalence cases (#)", "Prevalence cases (#)", "Prevalence cases (#)",
                          "Incidence rate", "Incidence rate","Incidence rate",
                          "Incidence cases", "Incidence cases", "Incidence cases")

  var_names_sexes <- c("Males", "Females","Both sexes",
                       "Males", "Females","Both sexes",
                       "Males", "Females","Both sexes",
                       "Males", "Females","Both sexes")

  idxtoplot <- list(c("EstimatePrevF", "PrevLB_2.5%F", "PrevUB_97.5%F"),
                    c("EstimatePrevM", "PrevLB_2.5%M", "PrevUB_97.5%M"),
                    c("EstimatePrevM+F", "PrevLB_2.5%M+F", "PrevUB_97.5%M+F"),
                    c("CasePrev_EstF", "CasePrev_LB_2.5%F",  "CasePrev_UB_97.5%F"),
                    c("CasePrev_EstM", "CasePrev_LB_2.5%M",  "CasePrev_UB_97.5%M"),
                    c("CasePrev_EstM+F", "CasePrev_LB_2.5%M+F",  "CasePrev_UB_97.5%M+F"),
                    c("EstimateIncF", "IncLB_2.5%F", "IncUB_97.5%F"),
                    c("EstimateIncM", "IncLB_2.5%M", "IncUB_97.5%M"),
                    c("EstimateIncM+F", "IncLB_2.5%M+F", "IncUB_97.5%M+F"),
                    c("CaseInc_EstF", "CaseInc_LB_2.5%F", "CaseInc_UB_97.5%F"),
                    c("CaseInc_EstM", "CaseInc_LB_2.5%M", "CaseInc_UB_97.5%M"),
                    c("CaseInc_EstM+F", "CaseInc_LB_2.5%M+F", "CaseInc_UB_97.5%M+F")

  )

  LongRegDataPrevIncForPlots <- data.frame()
  if(is.null(list_countries))
  {
    for(ii in seq_len(length(idxtoplot)))
    {
      v_name <- var_names_for_long[ii];
      idxb <- which(names(RegDataPrevInc)==idxtoplot[[ii]][1])
      idxlb <- which(names(RegDataPrevInc)==idxtoplot[[ii]][2])
      idxub <- which(names(RegDataPrevInc)==idxtoplot[[ii]][3])

      temp_data <- RegDataPrevInc[,which(is.element(names(RegDataPrevInc),c("Country", "SDG_Region","ISO code","Year")))]
      temp_data$indicator <- v_name
      temp_data$sex <- var_names_sexes[ii]
      temp_data$value <- RegDataPrevInc[,idxb]
      temp_data$lower <- RegDataPrevInc[,idxlb]
      temp_data$upper <- RegDataPrevInc[,idxub]
      LongRegDataPrevIncForPlots <- rbind(LongRegDataPrevIncForPlots,temp_data)
    }
  }#End if(!is.null(list_countries))

  ##############################################################################
  ##############################################################################
  results$CongenDataOut <- CongenDataOut
  results$LongCongenDataOutForPlots <- LongCongenDataOutForPlots
  results$RegCSABO <- RegCSABO
  results$LongRegCSABO <- LongRegCSABO
  results$RegDataPrevInc <- RegDataPrevInc
  results$LongRegDataPrevIncForPlots <- LongRegDataPrevIncForPlots
  results$all_iso = all_countries_iso
  results$CongenDataIn <- CongenDataIn
  results$CongenDataRaw <- CongenDataRaw

  wb <- openxlsx::createWorkbook()
  nnames <- c("CS Estimates", "CS Estimates long format", "Regional CS ABO", "Regional CS ABO long",
              "Reg. Prev. and Inc. Est.", "Reg. Prev. and Inc Est. long")
  for(ii in 1:6)
  {
    sheet <- openxlsx::addWorksheet(wb,sheetName=nnames[ii])
    temp_dat = results[[ii]]
    openxlsx::writeData(wb,nnames[ii],data.matrix(temp_dat))
  }

  #Simplimfied file for Jane
  CongenDataOutJane <- CongenDataOut[,names(CongenDataOut)%in%c("Country","ISO code", "Year", "Live Births", "Still births", "Pregnancies", "Women with >= 1 ANC visit (%)",
                                                                "Syphilis-tested (any ANC visit,%)", "EstimatePrevPregWom", "CS cases","CS / 100,000 live births",
                                                                "ABO = CS, ANC women treated", "CS cases, ANC-screened women not treated", "ABO cases",
                                                                "ABO, ANC-screened women not treated", "ABO, not seen in ANC", "ABO, ANC women not screened",
                                                                "CS cases, not seen in ANC","CS cases, ANC women not screened", "CS cases, ANC-screened women not treated")]

  wbsimple <- openxlsx::createWorkbook()
  nnames <- c("CS Estimates", "CS Estimates long format", "Regional CS ABO", "Regional CS ABO long",
              "Reg. Prev. and Inc. Est.", "Reg. Prev. and Inc Est. long")
  for(ii in 1:6)
  {
    sheet <- openxlsx::addWorksheet(wbsimple,sheetName=nnames[ii])
    temp_dat = results[[ii]]
    if(ii ==1) temp_dat = CongenDataOutJane
    openxlsx::writeData(wbsimple,nnames[ii],data.matrix(temp_dat))
  }

  results$wb <- wb
  results$wbsimple <- wbsimple
  class(results) <- "CSProj"
  return(results)
}

#' Congenital Syphilis estimation and projection using results from the procedure RunFitSyphilis
#'
#' @param syphfitobj an object of class "Syph-fit", output from RunFitSyphilis.
#' @param list_countries Vector containing the alpha-numeric ISO codes of countries Congenital Syphilis estimates are to be estimated.
#' @param proj_years integer vector of years for which estimates and projections are to be made. By default, this is the sequence 1990:2025
#' @param min_year non-zero real number indicating the earliest year for data to be used. This is 2011 by default.
#' @param rSyph_preg Ratio of Syphilis prevalence among all women to prevalence among pregnant women. This is 1 by default.
#' @param CSinputfiles list of four file names for (1) prevalence: syphilis prevalence used for prevalence estimates, (2) screening: data for Syphilis screen among pregnant women, (3) csdb: congenital syphilis prevalence, and (4) demographics: all demographic data, including all births by country. By default, this is set to NULL, which means data 2023 default data bases are used.
#' @return A list of class "CSProj", of dataframes CongenDataOut, LongCongenDataOutForPlots, RegCSABO, LongRegCSABO, RegDataPrevInc, LongRegDataPrevIncForPlots containing national and regional congenital Syphilis estimates.
#' @examples Not available
CalcCS <- function(syphfitobj, list_countries=NULL, proj_years=1990:2025, min_year=2011, rSyph_preg=1, CSinputfiles=NULL, f_DiagnosticTest)
{
  atmfn <- c(LETTERS[sample(1:15)],letters[sample(1:15)])
  tmfn <- atmfn[1]
  for(xx in atmfn[-1]) tmfn <- paste(tmfn,xx, sep="")
  tmfn <- paste(tmfn,".xlsx", sep="")
  saveSyphfit(syphfitobj,tmfn)
  result <- CalcCS_p(tmfn, list_countries, proj_years,min_year, rSyph_preg, CSinputfiles, f_DiagnosticTest)
  if(file.exists(tmfn)) file.remove(tmfn)
  return(result)
}


#' Save CS estimates to xlsx file
#'
#' @param xCSProj object of class "CSProj". This is a list containing Congenital Syphilis estimates.
#' @param fname name of the xlsx file to which the results are to be saved.
#' @return boolean, TRUE upon completion of the task
#' @examples Not available
saveCS <- function(xCSProj, fname=NULL)
{
  if(class(xCSProj)!="CSProj") stop("class(xCSProj) should be xCSProj")
  fnameout <- fname
  if(is.null(fname))
  {
    chardate <- Sys.Date();
    chardate <- gsub("-","", chardate)
    fnameout <- paste(paste("CSProj", chardate, sep="_"),".xlsx",sep="")
  }

  openxlsx::saveWorkbook(xCSProj$wb,file=fnameout, overwrite = T)
  return("File saved")
}

#' Plot Syphilis prevalence trends estimates by country
#'
#' @param list of class "Syph-fit" (containing a file name and the data used to run the fits)
#' @param ctrname Country's alpha-numeric ISO code.
#' @param sex A character indicating the sex for which estimates are to be drawn. Accepted values are "both" (both sexes, default), "males" and "females".
#' @param years A numeric vector of positions for years on the y-axis. This set to 2010:2021 by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_ctr_SyphPrev <- function(syphfits, ctr_iso3, sex="both", years= 2010:2021, fn_population="All")
{
  if(is.null(syphfits)) return(NULL)

  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  library(cowplot)

  if(length(ctr_iso3)!=1)
  {
    stop("length(ctr_iso3) must be 1")
  }

  llev_y <- as.factor(years)
  num_y <- as.numeric(as.character(llev_y))

  mout_file <- syphfits$wb;
  SyphData <- syphfits$SyphData
  #df=xCSProj$LongCongenDataOutForPlots
  ctr <- ctr_iso3

  all_res <- openxlsx::read.xlsx(mout_file,sheet="SYPH_RBootstrap_All")
  if(is.null(all_res)) return(NULL)
  all_res$datatype <- "Model"
  all_res$weight <- 1


  rgp <- syphfits$ListRiskGroupsComp

  #Prevalence rate
  names(all_res)[names(all_res)=="MedianPrevF"] = "PrevMed_Females"
  names(all_res)[names(all_res)=="EstimatePrevF"] = "PrevEst_Females"
  names(all_res)[names(all_res)=="PrevLB_2.5%F"] = "PrevLB_Females"
  names(all_res)[names(all_res)=="PrevUB_97.5%F"] = "PrevUB_Females"

  names(all_res)[names(all_res)=="MedianPrevM"] = "PrevMed_Males"
  names(all_res)[names(all_res)=="EstimatePrevM"] = "PrevEst_Males"
  names(all_res)[names(all_res)=="PrevLB_2.5%M"] = "PrevLB_Males"
  names(all_res)[names(all_res)=="PrevUB_97.5%M"] = "PrevUB_Males"

  names(all_res)[names(all_res)=="MedianPrevM+F"] = "PrevMed_BothSexes"
  names(all_res)[names(all_res)=="EstimatePrevM+F"] = "PrevEst_BothSexes"
  names(all_res)[names(all_res)=="PrevLB_2.5%M+F"] = "PrevLB_BothSexes"
  names(all_res)[names(all_res)=="PrevUB_97.5%M+F"] = "PrevUB_BothSexes"

  #Case prevalence
  names(all_res)[names(all_res)=="CasePrev_EstF"] = "CasePrevEst_Females"
  names(all_res)[names(all_res)=="CasePrev_MedF"] = "CasePrevMed_Females"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%F"] = "CasePrevLB_Females"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%F"] = "CasePrevUB_Females"

  names(all_res)[names(all_res)=="CasePrev_EstM"] = "CasePrevEst_Males"
  names(all_res)[names(all_res)=="CasePrev_MedM"] = "CasePrevMed_Males"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%M"] = "CasePrevLB_Males"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%M"] = "CasePrevUB_Males"

  names(all_res)[names(all_res)=="CasePrev_EstM+F"] = "CasePrevEst_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_MedM+F"] = "CasePrevMed_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%M+F"] = "CasePrevLB_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%M+F"] = "CasePrevUB_BothSexes"

  #*Incidence rates
  names(all_res)[names(all_res)=="EstimateIncF"] = "InciMed_Females"
  names(all_res)[names(all_res)=="MedianIncF"] = "InciEst_Females"
  names(all_res)[names(all_res)=="IncLB_2.5%F"] = "InciLB_Females"
  names(all_res)[names(all_res)=="IncUB_97.5%F"] = "InciUB_Females"

  names(all_res)[names(all_res)=="EstimateIncM"] = "InciMed_Males"
  names(all_res)[names(all_res)=="MedianIncM"] = "InciEst_Males"
  names(all_res)[names(all_res)=="IncLB_2.5%M"] = "InciLB_Males"
  names(all_res)[names(all_res)=="IncUB_97.5%M"] = "InciUB_Males"

  names(all_res)[names(all_res)=="EstimateIncM+F"] = "InciMed_BothSexes"
  names(all_res)[names(all_res)=="MedianIncM+F"] = "InciEst_BothSexes"
  names(all_res)[names(all_res)=="IncLB_2.5%M+F"] = "InciLB_BothSexes"
  names(all_res)[names(all_res)=="IncUB_97.5%M+F"] = "InciUB_BothSexes"

  #*Incidence cases
  names(all_res)[names(all_res)=="CaseInc_EstF"] = "CaseInciMed_Females"
  names(all_res)[names(all_res)=="CaseInc_MedF"] = "CaseInciEst_Females"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%F"] = "CaseInciLB_Females"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%F"] = "CaseInciUB_Females"

  names(all_res)[names(all_res)=="CaseInc_EstM"] = "CaseInciMed_Males"
  names(all_res)[names(all_res)=="CaseInc_MedM"] = "CaseInciEst_Males"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%M"] = "CaseInciLB_Males"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%M"] = "CaseInciUB_Males"

  names(all_res)[names(all_res)=="CaseInc_EstM+F"] = "CaseInciMed_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_MedM+F"] = "CaseInciEst_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%M+F"] = "CaseInciLB_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%M+F"] = "CaseInciUB_BothSexes"

  #KP
  names(all_res)[names(all_res)=="EstimatePrevFSW"] = "PrevEst_FSW"
  names(all_res)[names(all_res)=="MedianPrevFSW"] = "PrevMed_FSW"
  names(all_res)[names(all_res)=="PrevLB_2.5%FSW"] = "PrevLB_FSW"
  names(all_res)[names(all_res)=="PrevUB_97.5%FSW"] = "PrevUB_FSW"

  names(all_res)[names(all_res)=="EstimatePrevMSM"] = "PrevEst_MSM"
  names(all_res)[names(all_res)=="MedianPrevMSM"] = "PrevMed_MSM"
  names(all_res)[names(all_res)=="PrevLB_2.5%MSM"] = "PrevLB_MSM"
  names(all_res)[names(all_res)=="PrevUB_97.5%MSM"] = "PrevUB_MSM"


  names(all_res)[names(all_res)=="CasePrev_EstFSW"] = "CasePrevEst_FSW"
  names(all_res)[names(all_res)=="CasePrev_MedFSW"] = "CasePrevMed_FSW"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%FSW"] = "CasePrevLB_FSW"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%FSW"] = "CasePrevUB_FSW"


  names(all_res)[names(all_res)=="CasePrev_EstMSM"] = "CasePrevEst_MSM"
  names(all_res)[names(all_res)=="CasePrev_MedMSM"] = "CasePrevMed_MSM"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%MSM"] = "CasePrevLB_MSM"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%MSM"] = "CasePrevUB_MSM"

  names(all_res)[names(all_res)=="EstimatePrevPregWom"] = "PrevEst_PregWom"
  names(all_res)[names(all_res)=="MedianPrevPregWom"] = "PrevMed_PregWom"
  names(all_res)[names(all_res)=="PrevLB_2.5%PregWom"] = "PrevLB_PregWom"
  names(all_res)[names(all_res)=="PrevUB_97.5%PregWom"] = "PrevUB_PregWom"

  names(all_res)[names(all_res)=="EstimatePrevGeneWom"] = "PrevEst_GeneWom"
  names(all_res)[names(all_res)=="MedianPrevGeneWom"] = "PrevMed_GeneWom"
  names(all_res)[names(all_res)=="PrevLB_2.5%GeneWom"] = "PrevLB_GeneWom"
  names(all_res)[names(all_res)=="PrevUB_97.5%GeneWom"] = "PrevUB_GeneWom"

  names(all_res)[names(all_res)=="EstimatePrevGeneMen"] = "PrevEst_GeneMen"
  names(all_res)[names(all_res)=="MedianPrevGeneMen"] = "PrevMed_GeneMen"
  names(all_res)[names(all_res)=="PrevLB_2.5%GeneMen"] = "PrevLB_GeneMen"
  names(all_res)[names(all_res)=="PrevUB_97.5%GeneMen"] = "PrevUB_GeneMen"

  names(all_res)[names(all_res)=="EstimatePrevBloodDo"] = "PrevEst_BloodDo"
  names(all_res)[names(all_res)=="MedianPrevBloodDo"] = "PrevMed_BloodDo"
  names(all_res)[names(all_res)=="PrevLB_2.5%BloodDo"] = "PrevLB_BloodDo"
  names(all_res)[names(all_res)=="PrevUB_97.5%BloodDo"] = "PrevUB_BloodDo"

  if(!(sex%in%c("both", "males", "females"))) return(NULL)
  temp_long_ctr <- ctr_df <- subset(all_res, ISO3==ctr_iso3 & Year%in%years)


  if(nrow(ctr_df)==0) return(NULL)
  ctr <- ctr_df$Country[1]

  long_ctr_men <- data.frame(Country = ctr,
                             sex = "Males",
                             population = "All",
                             datatype = "Model",
                             weight = 1,
                             Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                             indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                             Median = c(temp_long_ctr$CaseInciMed_Males,temp_long_ctr$InciMed_Males,temp_long_ctr$CasePrevMed_Males,temp_long_ctr$PrevMed_Males),
                             BestFit = c(temp_long_ctr$CaseInciEst_Males,temp_long_ctr$InciEst_Males,temp_long_ctr$CasePrevEst_Males,temp_long_ctr$PrevEst_Males),
                             Lower = c(temp_long_ctr$CaseInciLB_Males,temp_long_ctr$InciLB_Males,temp_long_ctr$CasePrevLB_Males,temp_long_ctr$PrevLB_Males),
                             Upper = c(temp_long_ctr$CaseInciUB_Males,temp_long_ctr$InciUB_Males,temp_long_ctr$CasePrevUB_Males,temp_long_ctr$PrevUB_Males)
  )

  long_ctr_msm <- data.frame(Country = ctr,
                             sex = "Males",
                             population = "MSM",
                             datatype = "Model",
                             weight = 1,
                             Year =c(temp_long_ctr$Year,temp_long_ctr$Year),
                             indicator = rep(c("CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),2)),
                             Median = c(temp_long_ctr$CasePrevMed_MSM,temp_long_ctr$PrevMed_MSM),
                             BestFit = c(temp_long_ctr$CasePrevEst_MSM,temp_long_ctr$PrevEst_MSM),
                             Lower = c(temp_long_ctr$CasePrevLB_MSM,temp_long_ctr$PrevLB_MSM),
                             Upper = c(temp_long_ctr$CasePrevUB_MSM,temp_long_ctr$PrevUB_MSM)
  )

  long_ctr_genemen <- data.frame(Country = ctr,
                             sex = "Males",
                             population = "General-men",
                             datatype = "Model",
                             weight = 1,
                             Year =c(temp_long_ctr$Year),
                             indicator = rep(c("PrevalenceRate"), rep(nrow(temp_long_ctr),1)),
                             Median = c(temp_long_ctr$PrevMed_GeneMen),
                             BestFit = c(temp_long_ctr$PrevEst_GeneMen),
                             Lower = c(temp_long_ctr$PrevLB_GeneMen),
                             Upper = c(temp_long_ctr$PrevUB_GeneMen)
  )

  long_ctr_women <- data.frame(Country = ctr,
                               sex = "Females",
                               population = "All",
                               datatype = "Model",
                               weight = 1,
                               Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                               indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                               Median = c(temp_long_ctr$CaseInciMed_Females,temp_long_ctr$InciMed_Females,temp_long_ctr$CasePrevMed_Females,temp_long_ctr$PrevMed_Females),
                               BestFit = c(temp_long_ctr$CaseInciEst_Females,temp_long_ctr$InciEst_Females,temp_long_ctr$CasePrevEst_Females,temp_long_ctr$PrevEst_Females),
                               Lower = c(temp_long_ctr$CaseInciLB_Females,temp_long_ctr$InciLB_Females,temp_long_ctr$CasePrevLB_Females,temp_long_ctr$PrevLB_Females),
                               Upper = c(temp_long_ctr$CaseInciUB_Females,temp_long_ctr$InciUB_Females,temp_long_ctr$CasePrevUB_Females,temp_long_ctr$PrevUB_Females)
  )

  long_ctr_fsw <- data.frame(Country = ctr,
                               sex = "Females",
                               population = "FSW",
                               datatype = "Model",
                               weight = 1,
                               Year =c(temp_long_ctr$Year,temp_long_ctr$Year),
                               indicator = rep(c("CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),2)),
                               Median = c(temp_long_ctr$CasePrevMed_FSW,temp_long_ctr$PrevMed_FSW),
                               BestFit = c(temp_long_ctr$CasePrevEst_FSW,temp_long_ctr$PrevEst_FSW),
                               Lower = c(temp_long_ctr$CasePrevLB_FSW,temp_long_ctr$PrevLB_FSW),
                               Upper = c(temp_long_ctr$CasePrevUB_FSW,temp_long_ctr$PrevUB_FSW)
  )

  long_ctr_pregwom <- data.frame(Country = ctr,
                             sex = "Females",
                             population = "Pregnant women",
                             datatype = "Model",
                             weight = 1,
                             Year =c(temp_long_ctr$Year),
                             indicator = rep(c("PrevalenceRate"), rep(nrow(temp_long_ctr),1)),
                             Median = c(temp_long_ctr$PrevMed_PregWom),
                             BestFit = c(temp_long_ctr$PrevEst_PregWom),
                             Lower = c(temp_long_ctr$PrevLB_PregWom),
                             Upper = c(temp_long_ctr$PrevUB_PregWom)
  )


  long_ctr_genewom <- data.frame(Country = ctr,
                                 sex = "Females",
                                 population = "General-women",
                                 datatype = "Model",
                                 weight = 1,
                                 Year =c(temp_long_ctr$Year),
                                 indicator = rep(c("PrevalenceRate"), rep(nrow(temp_long_ctr),1)),
                                 Median = c(temp_long_ctr$PrevMed_GeneWom),
                                 BestFit = c(temp_long_ctr$PrevEst_GeneWom),
                                 Lower = c(temp_long_ctr$PrevLB_GeneWom),
                                 Upper = c(temp_long_ctr$PrevUB_GeneWom)
  )

  long_ctr_both <- data.frame(Country = ctr,
                              sex = "BothSexes",
                              population = "All",
                              datatype = "Model",
                              weight = 1,
                              Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                              indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                              Median = c(temp_long_ctr$CaseInciMed_BothSexes,temp_long_ctr$InciMed_BothSexes,temp_long_ctr$CasePrevMed_BothSexes,temp_long_ctr$PrevMed_BothSexes),
                              BestFit = c(temp_long_ctr$CaseInciEst_BothSexes,temp_long_ctr$InciEst_BothSexes,temp_long_ctr$CasePrevEst_BothSexes,temp_long_ctr$PrevEst_BothSexes),
                              Lower = c(temp_long_ctr$CaseInciLB_BothSexes,temp_long_ctr$InciLB_BothSexes,temp_long_ctr$CasePrevLB_BothSexes,temp_long_ctr$PrevLB_BothSexes),
                              Upper = c(temp_long_ctr$CaseInciUB_BothSexes,temp_long_ctr$InciUB_BothSexes,temp_long_ctr$CasePrevUB_BothSexes,temp_long_ctr$PrevUB_BothSexes)
  )

  long_ctr_blooddon <- data.frame(Country = ctr,
                                  sex = "Both",
                                  population = "BloodDonor",
                                  datatype = "Model",
                                  weight = 1,
                                  Year =c(temp_long_ctr$Year),
                                  indicator = rep(c("PrevalenceRate"), rep(nrow(temp_long_ctr),1)),
                                  Median = c(temp_long_ctr$PrevMed_BloodDo),
                                  BestFit = c(temp_long_ctr$PrevEst_BloodDo),
                                  Lower = c(temp_long_ctr$PrevLB_BloodDo),
                                  Upper = c(temp_long_ctr$PrevUB_BloodDo)
  )

  long_ctr <- rbind(long_ctr_men,long_ctr_women,long_ctr_both, long_ctr_blooddon)

  min_year <- min(long_ctr$Year)
  #Adding data points
  temp_ctr <- subset(SyphData,Country==ctr & !is.na(Weight_for_Spectrum_fitting))
  ppui <- temp_ctr$Prevalence
  ppui[ppui<=0] = 1/100;

  temp_all_res <- data.frame(Country = ctr,
                             sex = NA,
                             population = temp_ctr$Data_type,
                             datatype = "Reported",#datatype = temp_ctr$Data_type,
                             weight = temp_ctr$Weight_for_Spectrum_fitting,
                             Year=temp_ctr$Year,
                             indicator="PrevalenceRate",
                             Median = temp_ctr$Prevalence/100,
                             BestFit = temp_ctr$Prevalence/100,
                             Lower = NA,
                             Upper = NA
  )

  sdall <- sqrt(ppui/100*(1-ppui/100)/temp_ctr$N_tested)
  sdall[is.na(sdall)] <- 0
  temp_all_res$Lower <- pmax(temp_all_res$Median-1.96*sdall,0)
  temp_all_res$Upper <- pmin(temp_all_res$Median+1.96*sdall,1)

  temp_all_res$sex[temp_all_res$population=="ANC Routine screening"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="ANC Survey"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="BloodDonor Screening Men"] <- "Males"
  temp_all_res$sex[temp_all_res$population=="FSW"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="MSM"] <- "Males"
  temp_all_res$sex[temp_all_res$population=="Survey LowRisk Men"] <- "Males"
  temp_all_res$sex[temp_all_res$population=="BloodDonor Screening Men + Women"] <- "BothSexes"
  temp_all_res$sex[temp_all_res$population=="Male Sex Workers"] <- "Males"
  temp_all_res$sex[temp_all_res$population=="MSM + MSW combined"] <- "Males"
  temp_all_res$sex[temp_all_res$population=="PWID-Female"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="PWID-Male"] <- "Males"
  temp_all_res$sex[temp_all_res$population=="Survey LowRisk Women"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="Trans-Genders"] <- "BothSexes"
  temp_all_res$sex[temp_all_res$population=="BloodDonor Screening Women"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="Prisoners, Men"] <- "Males"
  temp_all_res$sex[temp_all_res$population=="Prisoners, Women"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="Wives of PWID"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="Survey LowRisk Men+Women"] <- "BothSexes"
  temp_all_res$sex[temp_all_res$population=="MSM"] <- "Males"
  temp_all_res$sex[temp_all_res$population=="FSW"] <- "Females"
  temp_all_res$sex[temp_all_res$population=="Other"] <- "BothSexes"

  long_ctr <- data.frame();
  mtitle <- "Syphilis prevalence trend among adults (15-49 y)"

  if(sex=="both")
  {
    long_ctr <- rbind(long_ctr_both, temp_all_res, long_ctr_blooddon)
  } else if(sex=="males")
  {
    temp_all_res$sex[temp_all_res$population%in%c(rgp$GeneMen, "All", "MSM")] <- "Males"
    long_ctr <- rbind(long_ctr_men, long_ctr_msm, long_ctr_genemen, temp_all_res)
    long_ctr <- subset(long_ctr, sex=="Males")
    mtitle <- "Syphilis prevalence trend among males (15-49 y)"
  } else if (sex=="females")
  {
    temp_all_res$sex[temp_all_res$population%in%c(rgp$GeneWom, rgp$PregWom, "All", "Pregnant women", "FSW")] <- "Females"
    long_ctr <- rbind(long_ctr_women,long_ctr_fsw, long_ctr_pregwom, long_ctr_genewom, temp_all_res)
    long_ctr <- subset(long_ctr, sex=="Females")
    mtitle <- "Syphilis prevalence trend among females (15-49 y)"
  }

  long_ctr <- subset(long_ctr,Year>=min_year)
  long_ctr$population <- factor(long_ctr$population, levels=c("All","ANC Routine screening","ANC Survey","BloodDonor Screening Men",
                                                          "FSW","MSM","Survey LowRisk Men","BloodDonor Screening Men + Women",
                                                          "Male Sex Workers","MSM + MSW combined","PWID-Female","PWID-Male",
                                                          "Survey LowRisk Women","Trans-Genders","BloodDonor Screening Women",
                                                          "Prisoners, Men","Prisoners, Women","Wives of PWID","Survey LowRisk Men+Women",
                                                          "Other", "BloodDonor",
                                                          "Pregnant women/General-women","Blood donors/General-men", "Blood donors",
                                                          "Blood donors/General-women","Blood donors/General-women/General-men",
                                                          "General-men","General-women/General-men","General-women", "Pregnant women"))

  show_chart <- TRUE
  if(fn_population!="All")
  {
    temp_pop = NULL
    if(fn_population=="General-men" )
    {
      temp_pop = c(rgp$GeneMen,"General-men")
      mtitle <- "Syphilis prevalence trend among General-men (15-49 y)"
    } else if(fn_population=="General-women" )
    {
      temp_pop = c(rgp$GeneWom, "General-women")
      mtitle <- "Syphilis prevalence trend among General-women (15-49 y)"
    } else if(fn_population=="Pregnant women" )
    {
      temp_pop = c(rgp$PregWom,"Pregnant women")
      mtitle <- "Syphilis prevalence trend among pregnant women (15-49 y)"
      if(sum(temp_all_res$population%in%c("ANC Routine screening","ANC Routine screening"))==0)
      {
        show_chart <- FALSE
      }
    } else if(fn_population=="FSW" )
    {
      temp_pop = "FSW"
      mtitle <- "Syphilis prevalence trend among Females Sex Workers (FSW)"
    } else if(fn_population=="MSM" )
    {
      temp_pop = "MSM"
      mtitle <- "Syphilis prevalence trend among Men who have Sex with Men (MSM)"
    } else if(fn_population=="BloodDonor" )
    {
      temp_pop = c("BloodDonor Screening Men","BloodDonor Screening Men + Women",
                   "BloodDonor Screening Women", "BloodDonor",
                   "Blood donors", "Blood donors/General-women","Blood donors/General-women/General-men")
      mtitle <- "Syphilis prevalence trend among Blood Donors"

      if(sum(temp_all_res$population%in%temp_pop)==0)
      {
        show_chart <- FALSE
      }
    }

    if(!is.null(temp_pop))
    {
      long_ctr <- subset(long_ctr,population%in%temp_pop)
    }
  } else
  {
    long_ctr <- subset(long_ctr,!population%in%c("FSW","Pregnant women","MSM", "General-women","General-men"))
  }

  if(fn_population=="None")
  {
    long_ctr <- subset(long_ctr,datatype%in%"Reported")
  }

  if(nrow(long_ctr)==0) return(NULL)
  if(!show_chart) return(NULL)

  long_ctr$weight <- ifelse(long_ctr$weight==1,16,1)
  long_ctr$weight <- factor(long_ctr$weight, levels=c("16","1"))

  all_plots_fig_S1 <- long_ctr%>%filter(Country==ctr & indicator=="PrevalenceRate")%>%ggplot(aes(x=Year, y=100*Median, ymin = 100*Lower, ymax= 100*Upper,color=population,fill=population))+
    geom_line(data =. %>% filter(Country==ctr& datatype=="Model"), aes(Year, 100*BestFit)) +
    geom_ribbon(data =. %>% filter(Country==ctr& datatype=="Model"),alpha=0.2,linetype=0)+
    geom_pointrange(data =. %>% filter(Country==ctr& datatype!="Model"), position=position_jitter(w = 0.05, h = 0),#position=position_jitter(width=0.5),
                    linetype='solid', aes(shape=weight)) + guides(shape = "none") +
    expand_limits(y = 0) +
    labs(title = mtitle, x = "Year", y="Test-adjusted prevalence, active syphilis (%)") +
    theme_minimal() +
    theme(legend.position = "bottom")

  all_plots_fig_S1$labels$colour="Source"
  all_plots_fig_S1$labels$fill="Source"

  title <- ggdraw() +
    draw_label(
      ctr,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(2, 2, 2, 7)
    )

  all_plots_fig_S2 <- plot_grid(
    title, all_plots_fig_S1,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
    #rel_heights = c(1,2)
  )

  all_plots_fig_S2
}


#' Plot Syphilis incident cases trends estimates by country
#'
#' @param list of class "Syph-fit" (containing a file name and the data used to run the fits)
#' @param ctrname Country's alpha-numeric ISO code.
#' @param sex A character indicating the sex for which estimates are to be drawn. Accepted values are "both" (both sexes, default), "males" and "females".
#' @param years A numeric vector of positions for years on the y-axis. This set to 2010:2021 by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_ctr_SyphIncCases <- function(syphfits, ctr_iso3, sex="both", years= 2010:2021)
{
  if(is.null(syphfits)) return(NULL)
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  library(cowplot)

  if(length(ctr_iso3)!=1)
  {
    stop("length(ctr_iso3) must be 1")
  }

  llev_y <- as.factor(years)
  num_y <- as.numeric(as.character(llev_y))

  mout_file <- syphfits$wb;
  SyphData <- syphfits$SyphData
  #df=xCSProj$LongCongenDataOutForPlots
  ctr <- ctr_iso3

  all_res <- openxlsx::read.xlsx(mout_file,sheet="SYPH_RBootstrap_All")
  if(is.null(all_res)) return(NULL)
  all_res$datatype <- "Model"
  all_res$weight <- 1

  #Prevalence rate
  names(all_res)[names(all_res)=="MedianPrevF"] = "PrevMed_Females"
  names(all_res)[names(all_res)=="EstimatePrevF"] = "PrevEst_Females"
  names(all_res)[names(all_res)=="PrevLB_2.5%F"] = "PrevLB_Females"
  names(all_res)[names(all_res)=="PrevUB_97.5%F"] = "PrevUB_Females"

  names(all_res)[names(all_res)=="MedianPrevM"] = "PrevMed_Males"
  names(all_res)[names(all_res)=="EstimatePrevM"] = "PrevEst_Males"
  names(all_res)[names(all_res)=="PrevLB_2.5%M"] = "PrevLB_Males"
  names(all_res)[names(all_res)=="PrevUB_97.5%M"] = "PrevUB_Males"

  names(all_res)[names(all_res)=="MedianPrevM+F"] = "PrevMed_BothSexes"
  names(all_res)[names(all_res)=="EstimatePrevM+F"] = "PrevEst_BothSexes"
  names(all_res)[names(all_res)=="PrevLB_2.5%M+F"] = "PrevLB_BothSexes"
  names(all_res)[names(all_res)=="PrevUB_97.5%M+F"] = "PrevUB_BothSexes"

  #Case prevalence
  names(all_res)[names(all_res)=="CasePrev_EstF"] = "CasePrevEst_Females"
  names(all_res)[names(all_res)=="CasePrev_MedF"] = "CasePrevMed_Females"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%F"] = "CasePrevLB_Females"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%F"] = "CasePrevUB_Females"

  names(all_res)[names(all_res)=="CasePrev_EstM"] = "CasePrevEst_Males"
  names(all_res)[names(all_res)=="CasePrev_MedM"] = "CasePrevMed_Males"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%M"] = "CasePrevLB_Males"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%M"] = "CasePrevUB_Males"

  names(all_res)[names(all_res)=="CasePrev_EstM+F"] = "CasePrevEst_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_MedM+F"] = "CasePrevMed_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%M+F"] = "CasePrevLB_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%M+F"] = "CasePrevUB_BothSexes"

  #Incidence rates
  names(all_res)[names(all_res)=="EstimateIncF"] = "InciMed_Females"
  names(all_res)[names(all_res)=="MedianIncF"] = "InciEst_Females"
  names(all_res)[names(all_res)=="IncLB_2.5%F"] = "InciLB_Females"
  names(all_res)[names(all_res)=="IncUB_97.5%F"] = "InciUB_Females"

  names(all_res)[names(all_res)=="EstimateIncM"] = "InciMed_Males"
  names(all_res)[names(all_res)=="MedianIncM"] = "InciEst_Males"
  names(all_res)[names(all_res)=="IncLB_2.5%M"] = "InciLB_Males"
  names(all_res)[names(all_res)=="IncUB_97.5%M"] = "InciUB_Males"

  names(all_res)[names(all_res)=="EstimateIncM+F"] = "InciMed_BothSexes"
  names(all_res)[names(all_res)=="MedianIncM+F"] = "InciEst_BothSexes"
  names(all_res)[names(all_res)=="IncLB_2.5%M+F"] = "InciLB_BothSexes"
  names(all_res)[names(all_res)=="IncUB_97.5%M+F"] = "InciUB_BothSexes"

  #Incidence cases
  names(all_res)[names(all_res)=="CaseInc_EstF"] = "CaseInciMed_Females"
  names(all_res)[names(all_res)=="CaseInc_MedF"] = "CaseInciEst_Females"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%F"] = "CaseInciLB_Females"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%F"] = "CaseInciUB_Females"

  names(all_res)[names(all_res)=="CaseInc_EstM"] = "CaseInciMed_Males"
  names(all_res)[names(all_res)=="CaseInc_MedM"] = "CaseInciEst_Males"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%M"] = "CaseInciLB_Males"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%M"] = "CaseInciUB_Males"

  names(all_res)[names(all_res)=="CaseInc_EstM+F"] = "CaseInciMed_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_MedM+F"] = "CaseInciEst_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%M+F"] = "CaseInciLB_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%M+F"] = "CaseInciUB_BothSexes"

  names(all_res)[names(all_res)=="EstimatePrevFSW"] = "PrevEst_FSW"
  names(all_res)[names(all_res)=="MedianPrevFSW"] = "PrevMed_FSW"
  names(all_res)[names(all_res)=="PrevLB_2.5%FSW"] = "PrevLB_FSW"
  names(all_res)[names(all_res)=="PrevUB_97.5%FSW"] = "PrevUB_FSW"

  names(all_res)[names(all_res)=="EstimatePrevMSM"] = "PrevEst_MSM"
  names(all_res)[names(all_res)=="MedianPrevMSM"] = "PrevMed_MSM"
  names(all_res)[names(all_res)=="PrevLB_2.5%MSM"] = "PrevLB_MSM"
  names(all_res)[names(all_res)=="PrevUB_97.5%MSM"] = "PrevUB_MSM"

  names(all_res)[names(all_res)=="CasePrev_EstFSW"] = "CasePrevEst_FSW"
  names(all_res)[names(all_res)=="CasePrev_MedFSW"] = "CasePrevMed_FSW"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%FSW"] = "CasePrevLB_FSW"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%FSW"] = "CasePrevUB_FSW"

  names(all_res)[names(all_res)=="CasePrev_EstMSM"] = "CasePrevEst_MSM"
  names(all_res)[names(all_res)=="CasePrev_MedMSM"] = "CasePrevMed_MSM"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%MSM"] = "CasePrevLB_MSM"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%MSM"] = "CasePrevUB_MSM"

  if(!(sex%in%c("both", "males", "females"))) return(NULL)

  #long_all_res <- data.frame()
  temp_long_ctr <- ctr_df <- subset(all_res, ISO3==ctr_iso3 & Year%in%years)

  if(nrow(ctr_df)==0) return(NULL)
  ctr <- ctr_df$Country[1]

  long_ctr_men <- data.frame(Country = ctr,
                             sex = "Males",
                             population = "All",
                             datatype = "Model",
                             weight = 1,
                             Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                             indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                             Median = c(temp_long_ctr$CaseInciMed_Males,temp_long_ctr$InciMed_Males,temp_long_ctr$CasePrevMed_Males,temp_long_ctr$PrevMed_Males),
                             BestFit = c(temp_long_ctr$CaseInciEst_Males,temp_long_ctr$InciEst_Males,temp_long_ctr$CasePrevEst_Males,temp_long_ctr$PrevEst_Males),
                             Lower = c(temp_long_ctr$CaseInciLB_Males,temp_long_ctr$InciLB_Males,temp_long_ctr$CasePrevLB_Males,temp_long_ctr$PrevLB_Males),
                             Upper = c(temp_long_ctr$CaseInciUB_Males,temp_long_ctr$InciUB_Males,temp_long_ctr$CasePrevUB_Males,temp_long_ctr$PrevUB_Males)
  )

  long_ctr_women <- data.frame(Country = ctr,
                               sex = "Females",
                               population = "All",
                               datatype = "Model",
                               weight = 1,
                               Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                               indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                               Median = c(temp_long_ctr$CaseInciMed_Females,temp_long_ctr$InciMed_Females,temp_long_ctr$CasePrevMed_Females,temp_long_ctr$PrevMed_Females),
                               BestFit = c(temp_long_ctr$CaseInciEst_Females,temp_long_ctr$InciEst_Females,temp_long_ctr$CasePrevEst_Females,temp_long_ctr$PrevEst_Females),
                               Lower = c(temp_long_ctr$CaseInciLB_Females,temp_long_ctr$InciLB_Females,temp_long_ctr$CasePrevLB_Females,temp_long_ctr$PrevLB_Females),
                               Upper = c(temp_long_ctr$CaseInciUB_Females,temp_long_ctr$InciUB_Females,temp_long_ctr$CasePrevUB_Females,temp_long_ctr$PrevUB_Females)
  )


  long_ctr_both <- data.frame(Country = ctr,
                              sex = "BothSexes",
                              population = "All",
                              datatype = "Model",
                              weight = 1,
                              Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                              indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                              Median = c(temp_long_ctr$CaseInciMed_BothSexes,temp_long_ctr$InciMed_BothSexes,temp_long_ctr$CasePrevMed_BothSexes,temp_long_ctr$PrevMed_BothSexes),
                              BestFit = c(temp_long_ctr$CaseInciEst_BothSexes,temp_long_ctr$InciEst_BothSexes,temp_long_ctr$CasePrevEst_BothSexes,temp_long_ctr$PrevEst_BothSexes),
                              Lower = c(temp_long_ctr$CaseInciLB_BothSexes,temp_long_ctr$InciLB_BothSexes,temp_long_ctr$CasePrevLB_BothSexes,temp_long_ctr$PrevLB_BothSexes),
                              Upper = c(temp_long_ctr$CaseInciUB_BothSexes,temp_long_ctr$InciUB_BothSexes,temp_long_ctr$CasePrevUB_BothSexes,temp_long_ctr$PrevUB_BothSexes)
  )

  long_ctr <- data.frame();
  mtitle <- "Syphilis incident cases trend among adults (15-49 y)"

  if(sex=="both")
  {
    long_ctr <- rbind(long_ctr_both)
  } else if(sex=="males")
  {
    long_ctr <- rbind(long_ctr_men)
    long_ctr <- subset(long_ctr, sex=="Males")
    #mtitle <- paste(ctr," Syphilis prevalence trend among males (15-49 y)", sep=",")
    mtitle <- "Syphilis incident cases trend among males (15-49 y)"
  } else if (sex=="females")
  {
    long_ctr <- rbind(long_ctr_women)
    long_ctr <- subset(long_ctr, sex=="Females")
    #mtitle <- paste(ctr," Syphilis prevalence trend among females (15-49 y)", sep=",")
    mtitle <- "Syphilis incident cases trend among females (15-49 y)"
  }

  #long_ctr <- subset(long_ctr,Year>=min_year)
  long_ctr$population <- factor(long_ctr$population, levels=c("All","ANC Routine screening","ANC Survey","BloodDonor Screening Men",
                                                              "FSW","MSM","Survey LowRisk Men","BloodDonor Screening Men + Women",
                                                              "Male Sex Workers","MSM + MSW combined","PWID-Female","PWID-Male",
                                                              "Survey LowRisk Women","Trans-Genders","BloodDonor Screening Women",
                                                              "Prisoners, Men","Prisoners, Women","Wives of PWID","Survey LowRisk Men+Women",
                                                              "Other"))

  all_plots_fig_S1 <- long_ctr%>%filter(Country==ctr & indicator=="CaseIncidence")%>%ggplot(aes(x=Year, y=Median, ymin = Lower, ymax= Upper,color=population,fill=population))+
    geom_line(data =. %>% filter(Country==ctr& datatype=="Model"), aes(Year, BestFit)) +
    geom_ribbon(data =. %>% filter(Country==ctr& datatype=="Model"),alpha=0.2,linetype=0)+
    geom_pointrange(data =. %>% filter(Country==ctr& datatype!="Model"), position=position_jitter(w = 0.05, h = 0),#position=position_jitter(width=0.5),
                    linetype='solid', aes(size=weight),size=0.75)+
    expand_limits(y = 0) +
    labs(title = mtitle, x = "Year", y="New syphilis cases (/1000)") +
    theme_minimal() +
    theme(legend.position = "bottom")

  all_plots_fig_S1$labels$colour="Source"
  all_plots_fig_S1$labels$fill="Source"

  title <- ggdraw() +
    draw_label(
      ctr,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(2, 2, 2, 7)
    )

  all_plots_fig_S2 <- plot_grid(
    title, all_plots_fig_S1,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
    #rel_heights = c(1,2)
  )

  all_plots_fig_S2
}

#' Plot EMTCT overview by country
#'
#' @param xCSProj object of class "CSProj". This is a list containing COngenital Syphilis estimates.
#' @param ctr_iso3 Country's alpha-numeric ISO code.
#' @param years A numeric vector of positions for years on the y-axis. This set to 2010:2021 by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_ctr_EMTCT <- function(xCSProj, ctr_iso3, years= 2015:2021)
{
  if(is.null(xCSProj)) return(NULL)

  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  library(cowplot)

  if(length(ctr_iso3)!=1)
  {
    stop("length(ctr_iso3) must be 1")
  }

  llev_y <- as.factor(years)
  num_y <- as.numeric(as.character(llev_y))
  df=xCSProj$LongCongenDataOutForPlots
  ctr <- ctr_iso3

  dff_a <- subset(df,df$`ISO code`==ctr & indicator%in%c("Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)",
                                                         "Syphilis-tested (1st ANC, %)","Women with >= 1 ANC visit (%)","Treated (%)") & Year%in%num_y & !is.na(value))
  ctr_name <- df$Country[df$`ISO code`==ctr][1]

  dff_a$value[dff_a$indicator=="Treated (%)"] <- dff_a$value[dff_a$indicator=="Treated (%)"]/100
  dff_a$lower[dff_a$indicator=="Treated (%)"] <- dff_a$lower[dff_a$indicator=="Treated (%)"]/100
  dff_a$upper[dff_a$indicator=="Treated (%)"] <- dff_a$upper[dff_a$indicator=="Treated (%)"]/100

  dff_a$value[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$value[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100
  dff_a$lower[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$lower[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100
  dff_a$upper[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$upper[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100

  dff_a$value[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$value[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100
  dff_a$lower[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$lower[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100
  dff_a$upper[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$upper[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100

  dff_a$Year <- factor(dff_a$Year, levels=llev_y)
  if(nrow(dff_a)<=0) return(NULL);
  for(yy in llev_y)
  {
    if((!yy%in% as.character(unique(dff_a$Year))))
    {
      t_a <-  rbind(dff_a[1,],dff_a[1,],dff_a[1,],dff_a[1,])
      t_a$Year <- yy
      t_a$value <- 0
      t_a$lower <- NA
      t_a$upper <- NA
      t_a$indicator=c("Treated (%)","Maternal syphilis prevalence","Syphilis-tested (1st ANC, %)","Women with >= 1 ANC visit (%)")
      dff_a <- rbind(dff_a,t_a)
    }
  }#End for(yy in llev_y)

  dff_b <- subset(df,df$`ISO code`==ctr & indicator%in%c("CS / 100,000 live births") & Year%in%num_y & !is.na(value))
  dff_b$indicator[dff_b$indicator=="CS / 100,000 live births"] <- "CS case rate\n per 100,000 live births"
  dff_b$Year <- factor(dff_b$Year, levels=llev_y)

  nff_b <- nrow(dff_b)

  for(yy in llev_y)
  {
    if((!yy%in% as.character(unique(dff_b$Year))))
    {
      t_b <-  rbind(dff_b[1,])
      t_b$Year <- yy
      t_b$value <- 0
      t_b$lower <- NA
      t_b$upper <- NA
      t_b$indicator=c("CS case rate\n per 100,000 live births")
      dff_b <- rbind(dff_b,t_b)
    }
  }#End for(yy in llev_y)

  dff_a$indicator[dff_a$indicator=="Treated (%)"] <- "ANC-based syphilis\n treatment"
  dff_a$indicator[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- "ANC-based screening"
  dff_a$indicator[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- "ANC-1"

  dff_a$mcolor <- as.factor(dff_a$Year)
  dff_a$Year <- as.numeric(as.character(dff_a$Year))
  dff_a$indicator[dff_a$indicator=="Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)"] = "Maternal syphilis prevalence"
  dff_a$indicatortype <- ifelse(dff_a$indicator!="Maternal syphilis prevalence", "Treatment Cascade",
                                "Maternal syphilis prevalence")

  dff_a$mgroup <- paste(dff_a$indicator, dff_a$datatype, sep=", ")
  dff_a$mgroup[dff_a$mgroup=="Maternal syphilis prevalence, ANC Routine screening"] = "ANC Routine screening"
  dff_a$mgroup[dff_a$mgroup=="Maternal syphilis prevalence, ANC Survey"] = "ANC Survey"
  dff_a$mgroup[dff_a$mgroup=="Maternal syphilis prevalence, Projected"] = "Prevalence estimates"
  dff_a$mgroup[dff_a$mgroup=="Syphilis-tested (1st ANC, %), Projected"] = "% tested (1st ANC), Projected"
  dff_a$mgroup[dff_a$mgroup=="Syphilis-tested (1st ANC, %), Reported"] = "% tested (1st ANC), Reported"
  dff_a$mgroup[dff_a$mgroup=="Women with >= 1 ANC visit (%), Projected"] = "% with >= 1 ANC visit, Projected"
  dff_a$mgroup[dff_a$mgroup=="Women with >= 1 ANC visit (%), Reported"] = "% with >= 1 ANC visit, Reported"

  v_levels <- c("Prevalence estimates",
                "Treated (%), Projected","Treated (%), Reported",
                "% tested (1st ANC), Projected", "% with >= 1 ANC visit, Projected",
                "ANC Routine screening", "ANC Survey", "% tested (1st ANC), Reported",
                "% with >= 1 ANC visit, Reported")

  dff_a$mgroup <- factor(dff_a$mgroup, levels=v_levels)

  dff_a$mshape <- "+";
  dff_a$mshape[dff_a$mgroup%in%c("Treated (%), Reported","Syphilis-tested (1st ANC, %), Reported", "Women with >= 1 ANC visit (%), Reported",
                                 "ANC Routine screening","ANC Survey")] = "*"

  npregwomdata <- sum(dff_a$mgroup%in%c("ANC Routine screening","ANC Survey"))
  if(npregwomdata==0) return(NULL)

  dff_ap <- subset(dff_a, indicator%in%c("ANC-based syphilis\n treatment","ANC-based screening","ANC-1"))
  dff_ap$indicator <- factor(dff_ap$indicator,levels=c("ANC-1","ANC-based screening","ANC-based syphilis\n treatment"))

  dff_a <- data.frame()
  for(xx in unique(dff_ap$indicator))
  {
    temp_xx <- subset(dff_ap, indicator%in%xx)
    yy_xx <- temp_xx$Year[temp_xx$datatype=="Reported"]

    if(length(yy_xx)>=1)
    {
      idxrm <- sapply(yy_xx, function(iy){
        which(temp_xx$datatype=="Projected" & temp_xx$Year==iy)
      })

      idxrm <- unlist(idxrm)
      if(length(idxrm)>=1)
      {
        temp_xx <- temp_xx[-idxrm,]
      }#End if(length(idxrm)>=1)
    }#End if(length(yy_xx)>=1)
    dff_a <- rbind(dff_a,temp_xx)
  }

  dff_ar <- subset(dff_a,datatype=="Reported")
  dff_ap <- subset(dff_a,datatype!="Reported")


  dff_arb <- data.frame()
  for(mgp in unique(dff_ar$indicator))
  {
    tmp <- subset(dff_ar, indicator%in%mgp)
    if(nrow(tmp)>=1)
    {
      for(yy in unique(tmp$Year))
      {
        tmp_yy <- subset(tmp,Year==yy)
        if(nrow(tmp_yy)>=2)
        {
          tmp_yy$value <- mean(tmp_yy$value, na.rm=T)
          tmp_yy <- tmp_yy[1,]
        }
        dff_arb <- rbind(dff_arb,tmp_yy)
      }
    }
  }

  dff_a <- rbind(dff_arb,dff_ap)

  p_xa <- NULL
  #ANC Coverage
  p_xa <- dff_a%>%ggplot2::ggplot(aes(group=datatype, fill=datatype, color=datatype)) +
    geom_bar(aes(x=Year, y=100*value),stat="identity")+
    facet_grid(cols=vars(indicator))+
    labs(title = "", y = "Coverage (%)",x="")+
    ylim(0,NA)+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      #legend.position="none",
      legend.position="bottom",
      #legend.position="top",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=12),
      legend.text=element_text(size=11)
    )

  p_xa$labels$fill <- p_xa$labels$colour <- "Source"

  grobs <- ggplotGrob(p_xa)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

  p_xa <- p_xa +theme(legend.position = "none" )

  ##############################################################################
  ##############################################################################
  dff_bp <- dff_b
  dff_b <- data.frame()
  for(xx in unique(dff_bp$indicator))
  {
    temp_xx <- subset(dff_bp, indicator%in%xx)
    yy_xx <- temp_xx$Year[temp_xx$datatype=="Reported"]

    if(length(yy_xx)>=1)
    {
      idxrm <- sapply(yy_xx, function(iy){
        which(temp_xx$datatype=="Projected" & temp_xx$Year==iy)
      })

      idxrm <- unlist(idxrm)
      if(length(idxrm)>=1)
      {
        temp_xx <- temp_xx[-idxrm,]
      }#End if(length(idxrm)>=1)
    }#End if(length(yy_xx)>=1)
    dff_b <- rbind(dff_b,temp_xx)
  }

  f_bks <- sort(floor(seq(max(as.numeric(as.character(dff_b$Year))),min(as.numeric(as.character(dff_b$Year))), length=4)))
  ###Testing Coverage
  p_xb <- dff_b%>%ggplot2::ggplot(aes(group=datatype, fill=datatype, color=datatype)) +
    geom_bar(aes(x=Year, y=value),stat="identity")+
    facet_grid(cols=vars(indicator))+
    labs(title = "", y = "CS case rate\n per 100,000 live births",x="")+
    scale_x_discrete(breaks = f_bks)+
    scale_y_continuous(position = "right")+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      legend.position="none",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=12),
      legend.text=element_text(size=11)
    )

  p_xb$labels$fill <- p_xb$labels$colour <- "Source"

  all_plots_fig_S1 <- p_xa
  if(nff_b!=0){
    all_plots_fig_S1 <- plot_grid(p_xa, p_xb,ncol = 2, labels = NULL, rel_widths = c(3,1))
  }

  all_plots_fig_S1 <- plot_grid(all_plots_fig_S1, legend,ncol = 1, labels = NULL, rel_heights = c(1,.1))

  ttle <- paste(ctr_name, ", EMTCT", sep="")

  title <- ggdraw() +
    draw_label(
      ttle,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(2, 2, 2, 7)
    )

  all_plots_fig_S2 <- plot_grid(
    title, all_plots_fig_S1,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
    #rel_heights = c(1,2)
  )

  all_plots_fig_S2
}#End plot_ctr_EMTCT

#' Plot CS reported by country
#'
#' @param xCSProj object of class "CSProj". This is a list containing COngenital Syphilis estimates.
#' @param ctr_iso3 Country's alpha-numeric ISO code.
#' @param years A numeric vector of positions for years on the y-axis. This set to 2010:2021 by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_ctr_csreport <- function(xCSProj, ctr_iso3, years= 2010:2021)
{
  if(is.null(xCSProj)) return(NULL)
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  library(cowplot)

  if(length(ctr_iso3)!=1)
  {
    stop("length(ctr_iso3) must be 1")
  }

  llev_y <- as.factor(years)
  num_y <- as.numeric(as.character(llev_y))
  #df=xCSProj$CongenDataOut
  df=xCSProj$CongenDataIn
  ctr <- ctr_iso3

  short_dfc <- subset(df,df$`ISO code`==ctr & df$Year%in%years)
  nobs <- nrow(short_dfc)
  if(nobs==0) return(NULL)

  ctr_name <- df$Country[df$`ISO code`==ctr][1]

  long_dfc <- data.frame(
    Country = ctr_name,
    indicator = c(rep("Clinical CS", nobs),rep("Non-Clinical CS", nobs),rep("Stillbirths", nobs),rep("Premature/low birth weight", nobs),rep("Neonatal deaths", nobs)),
    Year = c(short_dfc$Year,short_dfc$Year,short_dfc$Year,short_dfc$Year,short_dfc$Year),
    value = c(short_dfc$`CS cases`-short_dfc$`Asymptomatic CS`,short_dfc$`Asymptomatic CS`, short_dfc$`Still births`, short_dfc$`Prematurity or LBW due to CS`, short_dfc$`Neonatal death due to CS`),#?
    datatype = c(short_dfc$`N any ANC visits`)
  )

  long_dfc <- subset(long_dfc, !is.na(value))
  if(nrow(long_dfc)==0) return(NULL)

  #CS
  p_xa <- long_dfc%>%ggplot2::ggplot(aes(group=indicator, fill=indicator, color=indicator)) +
    geom_bar(aes(x=Year, y=value),stat="identity")+
    labs(title = "CS-cases and outcomes", y = "Cases numbers",x="")+
    ylim(0,NA)+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      #legend.position="none",
      legend.position="bottom",
      #legend.position="top",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=12),
      legend.text=element_text(size=11)
    )

  p_xa$labels$fill <- p_xa$labels$colour <- "Source"

  #Combined
  #all_plots_fig_S1 <- plot_grid(p_a, p_b,ncol = 1, labels = "AUTO", rel_heights = c(2,3))
  all_plots_fig_S1 <- p_xa

  title <- ggdraw() +
    draw_label(
      ctr_name,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(2, 2, 2, 7)
    )

  all_plots_fig_S2 <- plot_grid(
    title, all_plots_fig_S1,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
    #rel_heights = c(1,2)
  )

  all_plots_fig_S2
}#End plot_ctr_csreport

#' Plot CS estimates only by country
#'
#' @param xCSProj object of class "CSProj". This is a list containing COngenital Syphilis estimates.
#' @param ctr_iso3 Country's alpha-numeric ISO code.
#' @param years A numeric vector of positions for years on the y-axis. This set to 2010:2021 by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_ctr_CSProj <- function(xCSProj, ctr_iso3, years= 2010:2021)
{
  if(is.null(xCSProj)) return(NULL)
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  library(cowplot)

  if(length(ctr_iso3)!=1)
  {
    stop("length(ctr_iso3) must be 1")
  }

  llev_y <- as.factor(years)
  num_y <- as.numeric(as.character(llev_y))
  df=xCSProj$CongenDataOut
  ctr <- ctr_iso3

  #
  short_dfc <- subset(df,df$`ISO code`==ctr & df$Year%in%years)
  nobs <- nrow(short_dfc)
  if(nobs==0) return(NULL)

  ctr_name <- df$Country[df$`ISO code`==ctr][1]

  long_dfc <- data.frame(
    Country = ctr_name,
    indicator = c(rep("Cases not seen in ANC", nobs),rep("Screened but  not treated", nobs),rep("Attended ANC but not screened", nobs),rep("Not in ANC", nobs)),
    Year = c(short_dfc$Year,short_dfc$Year,short_dfc$Year,short_dfc$Year),
    value = c(short_dfc$`CS cases, not seen in ANC`,short_dfc$`CS cases, ANC-screened women not treated`, short_dfc$`CS cases, ANC women not screened`, short_dfc$`ABO, not seen in ANC`),#?
    datatype = "Estimated"
  )

  long_dfc <- subset(long_dfc, !is.na(value))
  if(nrow(long_dfc)==0) return(NULL)

  #CS and outcome
  p_xa <- long_dfc%>%ggplot2::ggplot(aes(group=indicator, fill=indicator, color=indicator)) +
    geom_bar(aes(x=Year, y=value),stat="identity")+
    labs(title = "CS-Cases and missed opportunities", y = "Cases numbers",x="")+
    ylim(0,NA)+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      #legend.position="none",
      legend.position="bottom",
      #legend.position="top",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=12),
      legend.text=element_text(size=11)
    )

  p_xa$labels$fill <- p_xa$labels$colour <- "Source"

  #Combined
  #all_plots_fig_S1 <- plot_grid(p_a, p_b,ncol = 1, labels = "AUTO", rel_heights = c(2,3))
  all_plots_fig_S1 <- p_xa

  title <- ggdraw() +
    draw_label(
      ctr_name,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(2, 2, 2, 7)
    )

  all_plots_fig_S2 <- plot_grid(
    title, all_plots_fig_S1,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
    #rel_heights = c(1,2)
  )

  all_plots_fig_S2
} #plot_ctr_CSProj


#' Plot CS estimates by country
#'
#' @param xCSProj object of class "CSProj". This is a list containing COngenital Syphilis estimates.
#' @param ctr_iso3 Country's alpha-numeric ISO code.
#' @param fbreaks A numeric vector of positions for breaks on the y-axis. This set to c(2010,2012,2016,2020) by default.
#' @param years A numeric vector of positions for years on the y-axis. This set to 2010:2021 by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_all_CSProj <- function(xCSProj, ctr_iso3, fbreaks=c(2010,2012,2016,2020), years= 2010:2021)
{
  if(is.null(xCSProj)) return(NULL)
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  library(cowplot)

  if(length(ctr_iso3)!=1)
  {
    stop("length(ctr_iso3) must be 1")
  }

  llev_y <- as.factor(years)
  num_y <- as.numeric(as.character(llev_y))
  df=xCSProj$LongCongenDataOutForPlots
  ctr <- ctr_iso3

  if(ctr_iso3%in%c("AFR","AMR","EMR","WPR","EUR","SEAR"))
  {
    dff_a <- subset(df,df$Country==ctr & indicator%in%c("Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)",
                                                           "Syphilis-tested (1st ANC, %)","Women with >= 1 ANC visit (%)","Treated (%)") & Year%in%num_y & !is.na(value))

  } else
  {
    dff_a <- subset(df,df$`ISO code`==ctr & indicator%in%c("Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)",
                                                           "Syphilis-tested (1st ANC, %)","Women with >= 1 ANC visit (%)","Treated (%)") & Year%in%num_y & !is.na(value))
    ctr_name <- df$Country[df$`ISO code`==ctr][1]
  }

  dff_a$value[dff_a$indicator=="Treated (%)"] <- dff_a$value[dff_a$indicator=="Treated (%)"]/100
  dff_a$lower[dff_a$indicator=="Treated (%)"] <- dff_a$lower[dff_a$indicator=="Treated (%)"]/100
  dff_a$upper[dff_a$indicator=="Treated (%)"] <- dff_a$upper[dff_a$indicator=="Treated (%)"]/100

  dff_a$value[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$value[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100
  dff_a$lower[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$lower[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100
  dff_a$upper[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$upper[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100

  dff_a$value[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$value[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100
  dff_a$lower[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$lower[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100
  dff_a$upper[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$upper[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100

  dff_a$Year <- factor(dff_a$Year, levels=llev_y)
  if(nrow(dff_a)<=0) return(NULL);
  for(yy in llev_y)
  {
    if((!yy%in% as.character(unique(dff_a$Year))))
    {
      t_a <-  rbind(dff_a[1,],dff_a[1,],dff_a[1,],dff_a[1,])
      t_a$Year <- yy
      t_a$value <- 0
      t_a$lower <- NA
      t_a$upper <- NA
      t_a$indicator=c("Treated (%)","Maternal syphilis prevalence","Syphilis-tested (1st ANC, %)","Women with >= 1 ANC visit (%)")
      dff_a <- rbind(dff_a,t_a)
    }
  }#End for(yy in llev_y)
  dff_a$mcolor <- as.factor(dff_a$Year)
  dff_a$Year <- as.numeric(as.character(dff_a$Year))
  dff_a$indicator[dff_a$indicator=="Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)"] = "Maternal syphilis prevalence"
  dff_a$indicatortype <- ifelse(dff_a$indicator!="Maternal syphilis prevalence", "Treatment Cascade",
                                "Maternal syphilis prevalence")

  dff_a$mgroup <- paste(dff_a$indicator, dff_a$datatype, sep=", ")
  dff_a$mgroup[dff_a$mgroup=="Maternal syphilis prevalence, ANC Routine screening"] = "ANC Routine screening"
  dff_a$mgroup[dff_a$mgroup=="Maternal syphilis prevalence, ANC Survey"] = "ANC Survey"
  dff_a$mgroup[dff_a$mgroup=="Maternal syphilis prevalence, Projected"] = "Prevalence estimates"
  dff_a$mgroup[dff_a$mgroup=="Syphilis-tested (1st ANC, %), Projected"] = "% tested (1st ANC), Projected"
  dff_a$mgroup[dff_a$mgroup=="Syphilis-tested (1st ANC, %), Reported"] = "% tested (1st ANC), Reported"
  dff_a$mgroup[dff_a$mgroup=="Women with >= 1 ANC visit (%), Projected"] = "% with >= 1 ANC visit, Projected"
  dff_a$mgroup[dff_a$mgroup=="Women with >= 1 ANC visit (%), Reported"] = "% with >= 1 ANC visit, Reported"
  v_levels <- c("Prevalence estimates",
                "Treated (%), Projected","Treated (%), Reported",
                "% tested (1st ANC), Projected", "% with >= 1 ANC visit, Projected",
                "ANC Routine screening", "ANC Survey", "% tested (1st ANC), Reported",
                "% with >= 1 ANC visit, Reported")

  dff_a$mgroup <- factor(dff_a$mgroup, levels=v_levels)

  dff_a$mshape <- "+";
  dff_a$mshape[dff_a$mgroup%in%c("Treated (%), Reported","Syphilis-tested (1st ANC, %), Reported", "Women with >= 1 ANC visit (%), Reported",
                                 "ANC Routine screening","ANC Survey")] = "*"

  p_a <- ggplot2::ggplot(dff_a,aes(group=mgroup, fill=mgroup, color=mgroup)) +
    geom_line(data= .%>%subset(datatype=="Projected"), aes(x=Year, y=100*value)) +#"#619CFF"
    #scale_fill_manual(values=c("#CD9600", "#7CAE00", "#00BE67", "#00A9FF", "#000000",
    #                           "#F8766D", "#C77CFF","#00BFC4", "#FF61CC")[9:1])+
    geom_ribbon(data=.%>%subset(datatype=="Projected"),
                aes(x=Year, y=100*value, ymin=100*lower, ymax=100*upper), alpha=0.25, color=NA) +
    geom_errorbar(data=.%>%subset(datatype!="Projected"),
                  aes(x=Year, y=100*value, ymin=100*lower, ymax=100*upper), width=0.05, alpha=0.5, size=0.9) +
    geom_point(data=.%>%subset(datatype!="Projected"), aes(x=Year, y=100*value)) +
    labs(title = "", y = "Value (%)",x="")+
    ylim(0,NA)+
    scale_x_continuous(breaks=fbreaks)+
    facet_wrap(.~indicatortype,scales="free_y")+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      #legend.position="none",
      legend.position="bottom",
      #legend.position="top",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=12),
      legend.text=element_text(size=11)
    )

  p_a$labels$fill <- p_a$labels$colour <- "Source"
  dff_b <- subset(df,df$`ISO code`==ctr & indicator%in%c("ABO cases","ABO/100,000 live births",
                                                   "CS cases","CS / 100,000 live births") & Year%in%num_y & !is.na(value))
  if(nrow(dff_b)<=0) return(NULL);
  for(yy in llev_y)
  {
    if((!yy%in%as.character(unique(dff_b$Year))))
    {
      t_a <-  rbind(dff_b[1,],dff_b[1,],dff_b[1,],dff_b[1,])
      t_a$Year <- yy
      t_a$value <- 0
      t_a$lower <- NA
      t_a$upper <- NA
      t_a$indicator=c("ABO cases","ABO/100,000 live births","CS cases","CS / 100,000 live births")
      dff_b <- rbind(dff_b,t_a)
    }
  }#End for(yy in llev_y)
  dff_b$mcolor <- as.factor(dff_b$Year)
  dff_b$cases <- "Cases"
  dff_b$cases[dff_b$indicator%in%c("ABO/100,000 live births","CS / 100,000 live births")] <- "Rates (/100,000 live births)"
  dff_b$dtype <- "ABO"
  dff_b$dtype[dff_b$indicator%in%c("CS cases","CS / 100,000 live births")] <- "CS"
  dff_b$Year = as.numeric(as.character(dff_b$Year))
  p_b <- ggplot2::ggplot(dff_b) +
    geom_line( aes(x=Year, y=value), color = "#619CFF") +
    geom_ribbon( aes(x=Year, y=value, ymin=lower, ymax=upper), fill = "#619CFF", alpha=0.5, col=NA) +
    labs(title = "", y = "Value",x="Time (years)")+
    facet_grid(cases~dtype,scales="free_y")+
    ylim(0, NA)+
    scale_x_continuous(breaks=fbreaks)+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      legend.position="none",
      #legend.position="bottom",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=11)
    )

  #Combined
  all_plots_fig_S1 <- plot_grid(p_a, p_b,ncol = 1, labels = "AUTO", rel_heights = c(2,3))

  if(ctr_name=="AFR")
  {
    ctr_name = "African Region (AFR)"
  } else if(ctr_name =="AMR")
  {
    ctr_name = "American region"
  } else if(ctr_name=="EMR")
  {
    ctr_name = "Easter Mediterranean Region (EMR)"
  } else if(ctr_name=="WPR")
  {
    ctr_name = "Westen Pacific region"
  } else if(ctr_name=="EUR")
  {
    ctr_name= "European Region (EUR)"
  } else if(ctr_name=="SEAR")
  {
    ctr_name = "South East Asian region (SEAR)"
  }

  title <- ggdraw() +
    draw_label(
      ctr_name,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(2, 2, 2, 7)
    )

  all_plots_fig_S2 <- plot_grid(
    title, all_plots_fig_S1,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
    #rel_heights = c(1,2)
  )

  all_plots_fig_S2
}

#' Plot ANC, treatment or Syphilis testing coverage by country
#'
#' @param xCSProj object of class "CSProj". This is a list containing COngenital Syphilis estimates.
#' @param ctr_iso3 Country's alpha-numeric ISO code.
#' @param years A numeric vector of positions for years on the y-axis. This set to 2010:2021 by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_ctr_Coverage <- function(xCSProj, ctr_iso3, indicator="ANC Coverage",years= 2010:2021)
{
  if(is.null(xCSProj)) return(NULL)
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  library(cowplot)

  if(length(ctr_iso3)!=1)
  {
    stop("length(ctr_iso3) must be 1")
  }

  llev_y <- as.factor(years)
  num_y <- as.numeric(as.character(llev_y))
  df=xCSProj$LongCongenDataOutForPlots
  ctr <- ctr_iso3

  dff_a <- subset(df,df$`ISO code`==ctr & indicator%in%c("Syphilis-tested (1st ANC, %)","Women with >= 1 ANC visit (%)","Treated (%)") & Year%in%num_y & !is.na(value))
  ctr_name <- df$Country[df$`ISO code`==ctr][1]

  dff_a$value[dff_a$indicator=="Treated (%)"] <- dff_a$value[dff_a$indicator=="Treated (%)"]/100
  dff_a$lower[dff_a$indicator=="Treated (%)"] <- dff_a$lower[dff_a$indicator=="Treated (%)"]/100
  dff_a$upper[dff_a$indicator=="Treated (%)"] <- dff_a$upper[dff_a$indicator=="Treated (%)"]/100

  dff_a$value[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$value[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100
  dff_a$lower[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$lower[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100
  dff_a$upper[dff_a$indicator=="Syphilis-tested (1st ANC, %)"] <- dff_a$upper[dff_a$indicator=="Syphilis-tested (1st ANC, %)"]/100

  dff_a$value[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$value[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100
  dff_a$lower[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$lower[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100
  dff_a$upper[dff_a$indicator=="Women with >= 1 ANC visit (%)"] <- dff_a$upper[dff_a$indicator=="Women with >= 1 ANC visit (%)"]/100

  dff_a$Year <- factor(dff_a$Year, levels=llev_y)
  if(nrow(dff_a)<=0) return(NULL);
  for(yy in llev_y)
  {
    if((!yy%in% as.character(unique(dff_a$Year))))
    {
      t_a <-  rbind(dff_a[1,],dff_a[1,],dff_a[1,],dff_a[1,])
      t_a$Year <- yy
      t_a$value <- -1
      t_a$lower <- NA
      t_a$upper <- NA
      t_a$indicator=c("Treated (%)","Maternal syphilis prevalence","Syphilis-tested (1st ANC, %)","Women with >= 1 ANC visit (%)")
      dff_a <- rbind(dff_a,t_a)
    }
  }#End for(yy in llev_y)

  dff_a$mcolor <- as.factor(dff_a$Year)
  dff_a$Year <- as.numeric(as.character(dff_a$Year))

  dff_a$mgroup <- paste(dff_a$indicator, dff_a$datatype, sep=", ")
  dff_a$mgroup[dff_a$mgroup=="Syphilis-tested (1st ANC, %), Projected"] = "% tested (1st ANC), Projected"
  dff_a$mgroup[dff_a$mgroup=="Syphilis-tested (1st ANC, %), Reported"] = "% tested (1st ANC), Reported"
  dff_a$mgroup[dff_a$mgroup=="Women with >= 1 ANC visit (%), Projected"] = "% with >= 1 ANC visit, Projected"
  dff_a$mgroup[dff_a$mgroup=="Women with >= 1 ANC visit (%), Reported"] = "% with >= 1 ANC visit, Reported"
  v_levels <- c("Treated (%), Projected","Treated (%), Reported",
                "% tested (1st ANC), Projected", "% with >= 1 ANC visit, Projected",
                "ANC Routine screening", "ANC Survey", "% tested (1st ANC), Reported",
                "% with >= 1 ANC visit, Reported")

  dff_a$mgroup <- factor(dff_a$mgroup, levels=v_levels)

  dff_ap <- data.frame()
  for(mgp in unique(dff_a$mgroup))
  {
    tmp <- subset(dff_a, mgroup%in%mgp)
    if(nrow(tmp)>=1)
    {
      for(yy in unique(tmp$Year))
      {
        tmp_yy <- subset(tmp,Year==yy)
        if(nrow(tmp_yy)>=2)
        {
          tmp_yy$value <- mean(tmp_yy$value, na.rm=T)
          tmp_yy <- tmp_yy[1,]
        }
        dff_ap <- rbind(dff_ap,tmp_yy)
      }
    }
  }

  #ANC Coverage
  if(nrow(dff_ap)==0)
  {
    return(NULL)
  }

  p_xa <- dff_ap%>%ggplot2::ggplot(aes(group=mgroup, fill=mgroup, color=mgroup)) +
    geom_bar(data= .%>%subset(datatype=="Reported" & indicator=="Women with >= 1 ANC visit (%)"), aes(x=Year, y=100*value),stat="identity")+
    geom_line(data= .%>%subset(datatype=="Projected" & indicator=="Women with >= 1 ANC visit (%)"), aes(x=Year, y=100*value))+
    geom_ribbon(data=.%>%subset(datatype=="Projected" & indicator=="Women with >= 1 ANC visit (%)"),
                aes(x=Year, y=100*value, ymin=100*lower, ymax=100*upper), alpha=0.25, color=NA)+
    geom_point(data=.%>%subset(datatype!="Projected" & indicator=="Women with >= 1 ANC visit (%)"), aes(x=Year, y=100*value)) +
    labs(title = "", y = "Value (%)",x="")+
    ylim(0,NA)+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      #legend.position="none",
      legend.position="bottom",
      #legend.position="top",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=12),
      legend.text=element_text(size=11)
    )

  p_xa$labels$fill <- p_xa$labels$colour <- "Source"

  ###Testing Coverage
  p_xb <- dff_ap%>%ggplot2::ggplot(aes(group=mgroup, fill=mgroup, color=mgroup)) +
    geom_bar(data= .%>%subset(datatype=="Reported" & indicator=="Syphilis-tested (1st ANC, %)"), aes(x=Year, y=100*value),stat="identity")+
    geom_line(data= .%>%subset(datatype=="Projected" & indicator=="Syphilis-tested (1st ANC, %)"), aes(x=Year, y=100*value))+
    geom_ribbon(data=.%>%subset(datatype=="Projected" & indicator=="Syphilis-tested (1st ANC, %)"),
                aes(x=Year, y=100*value, ymin=100*lower, ymax=100*upper), alpha=0.25, color=NA)+
    geom_point(data=.%>%subset(datatype!="Projected" & indicator=="Syphilis-tested (1st ANC, %)"), aes(x=Year, y=100*value)) +
    labs(title = "", y = "Value (%)",x="")+
    ylim(0,NA)+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      #legend.position="none",
      legend.position="bottom",
      #legend.position="top",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=12),
      legend.text=element_text(size=11)
    )

  p_xb$labels$fill <- p_xb$labels$colour <- "Source"

  #Treatment Coverage
  p_xc <- dff_ap%>%ggplot2::ggplot(aes(group=mgroup, fill=mgroup, color=mgroup)) +
    geom_bar(data= .%>%subset(datatype=="Reported" & indicator=="Treated (%)"), aes(x=Year, y=100*value),stat="identity")+
    geom_line(data= .%>%subset(datatype=="Projected" & indicator=="Treated (%)"), aes(x=Year, y=100*value))+
    geom_ribbon(data=.%>%subset(datatype=="Projected" & indicator=="Treated (%)"),
                aes(x=Year, y=100*value, ymin=100*lower, ymax=100*upper), alpha=0.25, color=NA)+
    geom_point(data=.%>%subset(datatype!="Projected" & indicator=="Treated (%)"), aes(x=Year, y=100*value)) +
    labs(title = "", y = "Value (%)",x="")+
    ylim(0,NA)+
    theme(
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      #legend.position="none",
      legend.position="bottom",
      #legend.position="top",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=12),
      legend.text=element_text(size=11)
    )

  p_xc$labels$fill <- p_xc$labels$colour <- "Source"


  all_plots_fig_S1 <- NULL
  ttle <- ctr_name
  ##Combined
  #all_plots_fig_S1 <- plot_grid(p_a, p_b,ncol = 1, labels = "AUTO", rel_heights = c(2,3))
  if(indicator=="ANC Coverage")
  {
    all_plots_fig_S1 <- p_xa
    ttle <- paste(ttle, ", ANC coverage", sep="")
  } else if (indicator=="Testing Coverage")
  {
    all_plots_fig_S1 <- p_xb
    ttle <- paste(ttle, ", testing coverage", sep="")
  } else if (indicator=="Treatement Coverage")
  {
    all_plots_fig_S1 <- p_xc
    ttle <- paste(ttle, ", treatment coverage", sep="")
  }

  title <- ggdraw() +
    draw_label(
      ttle,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(2, 2, 2, 7)
    )

  all_plots_fig_S2 <- plot_grid(
    title, all_plots_fig_S1,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.05, 1)
    #rel_heights = c(1,2)
  )

  all_plots_fig_S2
}



#' Plot Global and Regional CS-ABO estimates
#'
#' @param xCSProj object of class "CSProj". This is a list containing COngenital Syphilis estimates.
#' @param proj_years Years for which estimates are displayed in charts.
#' @return plot, object from ggplot
#' @examples Not available
plot_GlobSyphPerBirth <- function(xCSProj, proj_years=c(2012,2016,2020))
{
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  require(cowplot)

  if(nrow(xCSProj$LongRegCSABO)==0) return(NULL)

  dcsabo <- subset(xCSProj$LongRegCSABO, indicator%in%c("Liveborn with clinical CS",
                                                "Prematurity or LBW due to CS",
                                                "Neonatal death due to CS",
                                                "Stillbirth due to CS",
                                                "Asymptomatic CS") & Year%in%proj_years)

  dcsabo$indicator[dcsabo$indicator=="Asymptomatic CS"] = "Non-clinical CS"
  dcsabo$indicator <- factor(dcsabo$indicator,levels=c("Non-clinical CS", "Stillbirth due to CS","Neonatal death due to CS",
                                                       "Prematurity or LBW due to CS", "Liveborn with clinical CS"))

  names(dcsabo)[names(dcsabo)=="SDG_Region"] <- "region"
  dcsabo <- subset(dcsabo,!is.na(dcsabo$region))

  dcsabo$region[dcsabo$region=="Oceania (excluding Australia & New Zealand)"] <- "Oceania"

  dcsabo$region <- factor(dcsabo$region, levels=c("Global", sort(unique(dcsabo$region)[unique(dcsabo$region)!="Global"])))

  levels(dcsabo$region)[levels(dcsabo$region)=="Australia & New Zealand"] <- "Australia &\n New Zealand"
  levels(dcsabo$region)[levels(dcsabo$region)=="Sub-Saharan Africa"] <- "Sub-Saharan\n Africa"
  levels(dcsabo$region)[levels(dcsabo$region)=="Northern Africa"] <- "Northern\n Africa"
  levels(dcsabo$region)[levels(dcsabo$region)=="South-eastern Asia"] <- "South-eastern\n Asia"
  levels(dcsabo$region)[levels(dcsabo$region)=="Central Asia"] <- "Central\n Asia"
  levels(dcsabo$region)[levels(dcsabo$region)=="Southern Asia"] <- "Southern\n Asia"
  levels(dcsabo$region)[levels(dcsabo$region)=="Eastern Asia"] <- "Eastern\n Asia"
  levels(dcsabo$region)[levels(dcsabo$region)=="Western Asia"] <- "Western\n Asia"
  levels(dcsabo$region)[levels(dcsabo$region)=="Northern America"] <- "Northern\n America"
  levels(dcsabo$region)[levels(dcsabo$region)=="Central America"] <- "Central\n America"
  levels(dcsabo$region)[levels(dcsabo$region)=="South America"] <- "South\n America"

  p_cs <- ggplot2::ggplot(dcsabo, aes(x = factor(Year), y=value , fill = indicator))+#ggplot(dcsabo, aes(x = factor(Year), y=value, fill = indicator))+
    geom_bar(stat = "identity")+
    labs(title = "Congenital Syphilis cases per 100,000 live births", y = "",x="")+
    facet_grid(col=vars(region),scales="free_y")+
    scale_fill_discrete(drop=FALSE) +
    theme(
      axis.title.y = element_text(size = rel(0.8)),
      axis.title.x = element_text(size = rel(0.8)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.05, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      legend.position="bottom",
      axis.line = element_line(linewidth = 0.5, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.01, linetype = 2,
                                      colour = "grey")
    )

  p_cs$labels$colour=""
  p_cs$labels$fill=""
  p_cs
}


#' Plot Global and Regional ABO estimates
#'
#' @param xCSProj object of class "CSProj". This is a list containing Congenital Syphilis estimates.
#' @param proj_years Years for which estimates are displayed in charts.
#' @return plot, object from ggplot
#' @examples Not available
plot_GlobABO <- function(xCSProj, proj_years=c(2012,2016,2020))
{
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  require(cowplot)

  if(nrow(xCSProj$LongRegCSABO)==0) return(NULL)

  dabo <- subset(xCSProj$LongRegCSABO, indicator%in%c("ABO, not seen in ANC", "ABO, ANC women not screened",
                                              "ABO, ANC-screened women not treated",
                                              "ABO, ANC women treated") & Year%in%proj_years)

  dabo$indicator[dabo$indicator=="ABO, not seen in ANC"] = "Not in ANC"
  dabo$indicator[dabo$indicator=="ABO, ANC women not screened"] = "In ANC, not tested"
  dabo$indicator[dabo$indicator=="ABO, ANC-screened women not treated"] = "ANC-tested but not treated"
  dabo$indicator[dabo$indicator=="ABO, ANC women treated"] = "ANC-tested and treated"
  dabo$indicator <- factor(dabo$indicator,levels=c("ANC-tested and treated","ANC-tested but not treated","In ANC, not tested", "Not in ANC"))

  names(dabo)[names(dabo)=="SDG_Region"] <- "region"

  dabo <- subset(dabo,!is.na(dabo$region))
  dabo$region[dabo$region=="Oceania (excluding Australia & New Zealand)"] <- "Oceania"
  dabo$region <- factor(dabo$region, levels=c("Global", sort(unique(dabo$region)[unique(dabo$region)!="Global"])))

  dabo <- subset(dabo, region!="Global")

  sdg_regions_order = c("Sub-Saharan Africa", "Northern Africa",  "Central Asia", "Southern Asia", "Eastern Asia",
                        "South-eastern Asia",
                        "Western Asia", "Caribbean", "Central America", "South America", "Australia & New Zealand",
                        "Oceania", "Northern America", "Europe", "Global")
  #levels(dabo$region) <- sdg_regions_order[sdg_regions_order%in%levels(dabo$region)]
  dabo$region <- factor(dabo$region, levels=sdg_regions_order[sdg_regions_order%in%levels(dabo$region)])
  levels(dabo$region)[levels(dabo$region)=="Australia & New Zealand"] <- "Australia &\n New Zealand"

  levels(dabo$region)[levels(dabo$region)=="Australia & New Zealand"] <- "Australia &\n New Zealand"
  levels(dabo$region)[levels(dabo$region)=="Sub-Saharan Africa"] <- "Sub-Saharan\n Africa"
  levels(dabo$region)[levels(dabo$region)=="Northern Africa"] <- "Northern\n Africa"
  levels(dabo$region)[levels(dabo$region)=="South-eastern Asia"] <- "South-eastern\n Asia"
  levels(dabo$region)[levels(dabo$region)=="Central Asia"] <- "Central\n Asia"
  levels(dabo$region)[levels(dabo$region)=="Southern Asia"] <- "Southern\n Asia"
  levels(dabo$region)[levels(dabo$region)=="Eastern Asia"] <- "Eastern\n Asia"
  levels(dabo$region)[levels(dabo$region)=="Western Asia"] <- "Western\n Asia"
  levels(dabo$region)[levels(dabo$region)=="Northern America"] <- "Northern\n America"
  levels(dabo$region)[levels(dabo$region)=="Central America"] <- "Central\n America"
  levels(dabo$region)[levels(dabo$region)=="South America"] <- "South\n America"

  dabo$region <- droplevels(dabo$region)

  p_abo <- ggplot2::ggplot(dabo, aes(x = factor(Year), y=value , fill = indicator))+#ggplot(dcsabo, aes(x = factor(Year), y=value, fill = indicator))+
    geom_bar(stat = "identity")+
    labs(title = "Adverse Birth Outcomes", y = "",x="")+
    facet_grid(col=vars(region),scales="free_y")+
    scale_fill_discrete(drop=FALSE) +
    theme(
      axis.title.y = element_text(size = rel(0.8)),
      axis.title.x = element_text(size = rel(0.8)),
      axis.text.x = element_text(angle = 45),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.05, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      legend.position="bottom",
      axis.line = element_line(linewidth = 0.5, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.01, linetype = 2,
                                      colour = "grey"),
      strip.text=element_text(size=7),
      legend.text=element_text(size=10)
    )

  p_abo$labels$colour=""
  p_abo$labels$fill=""
  p_abo
}

#' Plot Global and Regional Maternal Syphilis estimates
#'
#' @param xCSProj object of class "CSProj". This is a list containing COngenital Syphilis estimates.
#' @param proj_years A numeric vercot of ears for which estimates are displayed in charts. This is 2010:2021 by default.
#' @param fbreaks A numeric vector of positions for breaks on the y-axis. This set to c(2010,2016,2021) by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_GlobMatSypPrev <- function(xCSProj, proj_years=c(2010:2021), fbreaks=c(2010,2016,2021))
{
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  require(cowplot)

  if(nrow(xCSProj$LongRegCSABO)==0) return(NULL)

  temp_maternal_syph <- subset(xCSProj$LongCongenDataOutForPlots,Country%in%c("Global", unique(xCSProj$LongCongenDataOutForPlots$SDG_Region))
                               &indicator=="Maternal syphilis prevalence (= F adult, from Spectrum minus 10%)"
                               &Year%in%proj_years)
  temp_maternal_syph$`WHO Region` <- factor(temp_maternal_syph$`WHO Region`,
                                            levels=c("Global", "AFR", "AMR", "EMR",  "EUR",  "SEAR",  "WPR"))

  temp_maternal_syph$region <- temp_maternal_syph$SDG_Region
  temp_maternal_syph$region[temp_maternal_syph$region=="Oceania (excluding Australia & New Zealand)"] <- "Oceania"

  temp_maternal_syph$region <- factor(temp_maternal_syph$region, levels=c("Global", sort(unique(temp_maternal_syph$region)[unique(temp_maternal_syph$region)!="Global"])))

  sdg_regions_order = c("Sub-Saharan Africa", "Northern Africa",  "Central Asia", "Southern Asia", "Eastern Asia",
                        "South-eastern Asia",
                        "Western Asia", "Caribbean", "Central America", "South America", "Australia & New Zealand",
                        "Oceania", "Northern America", "Europe", "Global")

  temp_maternal_syph$region <- factor(temp_maternal_syph$region, levels=sdg_regions_order[sdg_regions_order%in%levels(temp_maternal_syph$region)])

  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Australia & New Zealand"] <- "Australia &\n New Zealand"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Sub-Saharan Africa"] <- "Sub-Saharan\n Africa"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Northern Africa"] <- "Northern\n Africa"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="South-eastern Asia"] <- "South-eastern\n Asia"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Central Asia"] <- "Central\n Asia"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Southern Asia"] <- "Southern\n Asia"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Eastern Asia"] <- "Eastern\n Asia"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Western Asia"] <- "Western\n Asia"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Northern America"] <- "Northern\n America"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="Central America"] <- "Central\n America"
  levels(temp_maternal_syph$region)[levels(temp_maternal_syph$region)=="South America"] <- "South\n America"

  p_maternal <- ggplot2::ggplot(temp_maternal_syph, aes(x = Year, y=100*value))+#ggplot(dcsabo, aes(x = factor(Year), y=value, fill = indicator))+
    geom_line( aes(x=Year, y=100*value), color = "#619CFF") +
    geom_ribbon( aes(x=Year, y=100*value, ymin=100*lower, ymax=100*upper), fill = "#619CFF", alpha=0.5, col=NA)+
    labs(title = "", y = "Prevalence (%)",x="Year")+
    facet_grid(col=vars(region),scales="free_y")+
    scale_x_continuous(breaks = fbreaks)+
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=12),
      axis.text.y = element_text(size=12),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8),
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      legend.position="none",
      #legend.position="bottom",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      #strip.text=element_text(size=5),
      legend.text=element_text(size=10)
    )
  p_maternal$labels$colour=""
  p_maternal$labels$fill=""
  p_maternal
}

#' Plot Global and Regional Syphilis Incidence and Prevalence Rates estimates
#'
#' @param xCSProj object of class "CSProj". This is a list containing Congenital Syphilis estimates.
#' @param proj_years Years for which estimates are displayed in charts.
#' @param fbreaks A numeric vector of positions for breaks on the y-axis. This set to c(2010,2016,2021) by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_GlobPrevIncR <- function(xCSProj, proj_years=c(2010:2021), fbreaks=c(2010,2016,2021))
{
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  require(cowplot)

  if(nrow(xCSProj$LongRegDataPrevIncForPlots)==0) return(NULL)

  df_a <- subset(xCSProj$LongRegDataPrevIncForPlots,Year%in%proj_years & indicator%in%c("Prevalence (%)","Incidence rate") & sex%in%"Both sexes")
  df_a$indicator[df_a$indicator=="Incidence rate"] <- "Incidence rate (%)"
  names(df_a)[names(df_a)=="SDG_Region"] <- "region"

  df_a$region[df_a$region=="Oceania (excluding Australia & New Zealand)"] <- "Oceania"
  df_a$region[df_a$region=="Global"] <- "All regions"
  df_a$region <- factor(df_a$region, levels=c("All regions", sort(unique(df_a$region)[unique(df_a$region)!="All regions"])))

  sdg_regions_order = c("Sub-Saharan Africa", "Northern Africa",  "Central Asia", "Southern Asia", "Eastern Asia",
                        "South-eastern Asia",
                        "Western Asia", "Caribbean", "Central America", "South America", "Australia & New Zealand",
                        "Oceania", "Northern America", "Europe", "All regions")

  df_a$region <- factor(df_a$region, levels=sdg_regions_order[sdg_regions_order%in%levels(df_a$region)])
  levels(df_a$region)[levels(df_a$region)=="Australia & New Zealand"] <- "Australia &\n New Zealand"
  levels(df_a$region)[levels(df_a$region)=="Sub-Saharan Africa"] <- "Sub-Saharan\n Africa"
  levels(df_a$region)[levels(df_a$region)=="Northern Africa"] <- "Northern\n Africa"
  levels(df_a$region)[levels(df_a$region)=="Central Asia"] <- "Central\n Asia"
  levels(df_a$region)[levels(df_a$region)=="Eastern Asia"] <- "Eastern\n Asia"
  levels(df_a$region)[levels(df_a$region)=="South-eastern Asia"] <- "South-eastern\n Asia"
  levels(df_a$region)[levels(df_a$region)=="Southern Asia"] <- "Southern\n Asia"
  levels(df_a$region)[levels(df_a$region)=="Western Asia"] <- "Western\n Asia"
  levels(df_a$region)[levels(df_a$region)=="Northern America"] <- "Northern\n America"
  levels(df_a$region)[levels(df_a$region)=="Central America"] <- "Central\n America"
  levels(df_a$region)[levels(df_a$region)=="South America"] <- "South\n America"

  df_a <- subset(df_a, !is.na(region))
  if(nrow(df_a)==0) return(NULL)

  p_a <- ggplot2::ggplot(df_a) +
    geom_line(aes(x=Year, y=100*value), color = "#619CFF", alpha=0.95) +
    geom_ribbon(aes(x=Year, y=100*value, ymin=100*lower, ymax=100*upper, color="#619CFF"), fill="#619CFF", alpha=0.25, colour = NA)+
    labs(title = "", y = "",x="Time (Years)")+
    facet_grid(row=vars(indicator),col=vars(region),scales="free_y")+
    scale_fill_discrete(drop=FALSE) +
    scale_x_continuous(breaks = fbreaks)+
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=8),
      axis.text.y = element_text(size=10),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8),
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      legend.position="none",
      #legend.position="bottom",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      legend.text=element_text(size=10)
    )
  p_a
}

#' Plot Global and Regional Syphilis Incidence and Prevalence Cases estimates
#'
#' @param xCSProj object of class "CSProj". This is a list containing Congenital Syphilis estimates.
#' @param proj_years Years for which estimates are displayed in charts.
#' @param fbreaks A numeric vector of positions for breaks on the y-axis. This set to c(2010,2016,2021) by default.
#' @return plot, object from ggplot
#' @examples Not available
plot_GlobPrevIncCases <- function(xCSProj, proj_years=c(2010:2021), fbreaks=c(2010,2016,2021))
{
  require(ggplot2)
  library(tidyverse)
  library(scales)
  library(gridExtra)
  require(cowplot)

  if(nrow(xCSProj$LongRegDataPrevIncForPlots)==0) return(NULL)

  df_b <- subset(xCSProj$LongRegDataPrevIncForPlots,Year%in%proj_years & indicator%in%c("Incidence cases","Prevalence cases (#)") & sex%in%"Both sexes")
  df_b$indicator[df_b$indicator=="Incidence cases"] <- "Incident cases (x1000)"
  df_b$indicator[df_b$indicator=="Prevalence cases (#)"] <- "Prevalent cases (x1000)"

  names(df_b)[names(df_b)=="SDG_Region"] <- "region"
  df_b$region[df_b$region=="Oceania (excluding Australia & New Zealand)"] <- "Oceania"

  df_b$region[df_b$region=="Global"] <- "All regions"
  df_b$region <- factor(df_b$region, levels=c("All regions", sort(unique(df_b$region)[unique(df_b$region)!="All regions"])))

  sdg_regions_order = c("Sub-Saharan Africa", "Northern Africa",  "Central Asia", "Southern Asia", "Eastern Asia",
                        "South-eastern Asia",
                        "Western Asia", "Caribbean", "Central America", "South America", "Australia & New Zealand",
                        "Oceania", "Northern America", "Europe", "All regions")

  df_b$region <- factor(df_b$region, levels=sdg_regions_order[sdg_regions_order%in%levels(df_b$region)])
  levels(df_b$region)[levels(df_b$region)=="Australia & New Zealand"] <- "Australia &\n New Zealand"
  levels(df_b$region)[levels(df_b$region)=="Sub-Saharan Africa"] <- "Sub-Saharan\n Africa"
  levels(df_b$region)[levels(df_b$region)=="Northern Africa"] <- "Northern\n Africa"
  levels(df_b$region)[levels(df_b$region)=="Central Asia"] <- "Central\n Asia"
  levels(df_b$region)[levels(df_b$region)=="Eastern Asia"] <- "Eastern\n Asia"
  levels(df_b$region)[levels(df_b$region)=="South-eastern Asia"] <- "South-eastern\n Asia"
  levels(df_b$region)[levels(df_b$region)=="Southern Asia"] <- "Southern\n Asia"
  levels(df_b$region)[levels(df_b$region)=="Western Asia"] <- "Western\n Asia"
  levels(df_b$region)[levels(df_b$region)=="Northern America"] <- "Northern\n America"
  levels(df_b$region)[levels(df_b$region)=="Central America"] <- "Central\n America"
  levels(df_b$region)[levels(df_b$region)=="South America"] <- "South\n America"

  df_b <- subset(df_b, !is.na(region))
  if(nrow(df_b)==0) return(NULL)

  #Combined
  p_b <- ggplot2::ggplot(df_b) +
    geom_line(aes(x=Year, y=value/1000), color = "#619CFF", alpha=0.95) +
    geom_ribbon(aes(x=Year, y=value/1000, ymin=lower/1000, ymax=upper/1000), fill = "#619CFF", alpha=0.25, colour = NA)+
    labs(title = "Syphilis prevalence and incidence estimates among adults", y = "",x="Time (Years)")+
    facet_grid(row=vars(indicator),col=vars(region),scales="free_y")+
    scale_fill_discrete(drop=FALSE) +
    scale_x_continuous(breaks = fbreaks)+
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=8),
      axis.text.y = element_text(size=10),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(size = 8),
      axis.title.y = element_text(size = rel(1.2)),
      axis.title.x = element_text(size = rel(1.2)),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(linewidth = 0.5, colour = "grey",linetype = 2),
      panel.ontop = TRUE,
      legend.position="none",
      axis.line = element_line(linewidth = 1, colour = "grey"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 2,
                                      colour = "grey"),
      legend.text=element_text(size=10)
    )
  p_b
}
