#include <Rcpp.h>
#include <iostream>
#include <omp.h>

using namespace Rcpp;
using namespace std;

double const pi=3.14159265358979323846;

double expit(double x)
{
    double result;
    if (x <= 0)
    {
        result = std::exp(x);
        result = result / (1 + result);
    }
    else
    {
        result = exp(-x);
        result = 1 / (1 + result);
    }
    return(result);
}

//#{Calculate the integral from 0 to x>0, of P, where P is a second order polynomial }
double IntegrePol2(const double&x, const double&a0, const double&b0, const double&c0) //#:double;// this function works only for x>0
{
  double result=a0*x+b0*x*x/2.0+c0*x*x*x/3;
  return(result);
}

//#{Calculate the integral from 0 to x>0, of max(P,0), where P is a second order spline }
// [[Rcpp::export]]
double IntegrePol2CPos(const double&x, const double&a0, const double&b0, const double&c0) //#:double;// this function works only for x>0
{
    double res = 0;
    double x1, x2, Delta, gepsilon;
    if (c0 == 0)
    {
        if ((a0 > 0) & (b0 > 0))
        {
            res = a0 * x + 0.5 * b0 * x * x;
        }
        else if ((a0 > 0) & (b0 < 0))
        {
            x1 = -a0 / b0;
            if (x < x1)
            {
                res = x * (a0 + 0.5 * b0 * x);
            }
            else
            {
                res = x1 * (a0 + 0.5 * b0 * x1);
            }
        }
        else if ((a0 > 0) & (b0 == 0))
        {
            res = a0 * x;
        }
        else if ((a0 < 0) & (b0 < 0))
        {
            res = 0;
        }
        else if ((a0 < 0) & (b0 > 0))
        {
            x1 = -a0 / b0;
            if (x < x1)
            {
                res = 0;
            }
            else
            {
                res = (x - x1) * (a0 + 0.5 * b0 * (x + x1));
            }
        }
        else if ((a0 <= 0) & (b0 == 0))
        {
            res = 0;
        }
        else if ((a0 == 0) & (b0 > 0))
        {
            res = x * (a0 + 0.5 * b0 * x);
        }
    }
    else if (c0 != 0)
    {
        Delta = b0 * b0 - 4.0 * a0 * c0;
        if (ISNAN(Delta))
        {
            res = 100000;
        };
        if (Delta < 0)
        {
            if (c0 <= 0)
            {
                res = 0;
            }
            else
            {
                res = x * (a0 + x * (0.5 * b0 + c0 * x / 3));
            }

        }
        else if (Delta >= 0)//#Delta>=0
        {
            if (Delta == 0)
            {
                if (c0 > 0)
                {
                    res = x * (a0 + x * (0.5 * b0 + c0 * x / 3));
                }
                else //#c0<=0
                {
                    res = 0;
                };
            }
            else
            {
                gepsilon = 1.0;
                if (b0 < 0)
                {
                    gepsilon = -1.0;
                }
                x1 = (-b0 - gepsilon * sqrt(Delta)) / 2.0 / c0;
                x2 = 2.0 * a0 / (-b0 - gepsilon * sqrt(Delta));

                if (x1 > x2)
                {
                    x1 = 2.0 * a0 / (-b0 - gepsilon * sqrt(Delta));
                    x2 = (-b0 - gepsilon * sqrt(Delta)) / 2.0 / c0;
                }

                if (c0 > 0)
                {
                    if (x2 < 0)
                    {
                        res = x * (a0 + x * (0.5 * b0 + c0 * x / 3));
                    }
                    else
                    {
                        if (x <= x1)
                        {
                            res = x * (a0 + x * (0.5 * b0 + c0 * x / 3));
                        }
                        else if ((x > x1) & (x <= x2))
                        {
                            if (x1 > 0)
                            {
                                res = x1 * (a0 + x1 * (0.5 * b0 + c0 * x1 / 3));
                            }
                            else
                            {
                                res = 0;
                            }
                        }
                        else if (x > x2)
                        {
                            res = (x - x2) * (a0 + 0.5 * b0 * (x + x2) + c0 / 3.0 * (x * x + x * x2 + x2 * x2));
                            if (x1 > 0)
                            {
                                res = res + x1 * (a0 + x1 * (0.5 * b0 + c0 * x1 / 3.0));
                            }
                        }
                    }

                }
                else //#c0<=0
                {
                    if (x2 < 0)
                    {
                        res = 0;
                    }
                    else
                    {
                        if (x <= x1)
                        {
                            res = 0;
                        }
                        else if ((x > x1) & (x <= x2))
                        {
                            if (x1 > 0)
                            {
                                res = (x - x1) * (a0 + 0.5 * b0 * (x + x1) + c0 / 3.0 * (x * x + x * x1 + x1 * x1));
                            }
                            else
                            {
                                res = x * (a0 + x * (0.5 * b0 + c0 * x / 3.0));
                            };

                        }
                        else if (x > x2)
                        {
                            if (x1 > 0)
                            {
                                res = (x2 - x1) * (a0 + 0.5 * b0 * (x2 + x1) + c0 / 3.0 * (x2 * x2 + x2 * x1 + x1 * x1));
                            }
                            else
                            {
                                res = x2 * (a0 + x2 * (0.5 * b0 + c0 * x2 / 3.0));
                            }
                        }
                    }
                }
            }
        }
    }
    return(res);
}

NumericVector VectLogit(const NumericVector&x)
{
  int n =x.size();
  NumericVector result(n);
  double sumx=0;
  for(int i=0;i<n;i++)
  {
    sumx+=x[i];
  }

  for(int i=0;i<n;i++)
  {
    result[i]=log(x[i]/(1-sumx));
  }
  return(result);
}

NumericVector Vectexpit(const NumericVector&vx)
{
  NumericVector result(vx.size());
  double xmax = max(vx);
  double denom=1;

  if(vx.size()==1)
  {
    result[0] = expit(vx[0]);
  } else
  {
    if(xmax>0)
    {
      denom = exp(-xmax);
      for(int i=0; i<vx.size(); i++)
      {
        result[i]=vx[i]-xmax;
      }
    }else
    {
      for(int i=0; i<vx.size(); i++)
      {
        result[i]=vx[i];
      }
    }
    for(int i=0; i<vx.size(); i++)
    {
      result[i]=exp(result[i]);
      denom+=result[i];
    }
    for(int i=0; i<vx.size(); i++)
    {
      result[i]=result[i]/denom;
    }
  }
  return(result);
}

NumericVector InciSegPol2PID(const NumericVector&vectx, const NumericVector&xpars, const double&fresknots)
{
  int n0, nknots;
  double ix;
  n0 = xpars.size()-3;
  n0 = n0/2;
  nknots = n0+1;

  NumericVector result(vectx.size());

  NumericVector valpha(nknots);
  NumericVector vbeta(nknots);
  NumericVector vgamma(nknots);

  valpha[0] = expit(xpars[0]);
  vbeta[0] = xpars[1];
  vgamma[0] = xpars[2];

  NumericVector vh(n0);
  NumericVector posknots(nknots+1);
  posknots[0]=0;
  posknots[nknots]=1;
  if(n0>=1)
  {
    for(int i=0; i<n0; i++)
    {
      vh[i]=xpars[3+2*i];
      vgamma[1+i]= xpars[3+2*i+1];
    }
    vh=Vectexpit(vh);
    for(int i=0; i<vh.size(); i++)
    {
      vh[i]=fresknots*vh[i];
    }
    for(int i=1; i<nknots; i++)
    {
      vbeta[i]=vbeta[i-1]+2.0*vgamma[i-1]*vh[i-1];
      valpha[i]=valpha[i-1]+vbeta[i-1]*vh[i-1]+vgamma[i-1]*vh[i-1]*vh[i-1];
    }

    for(int i=0; i<vh.size(); i++)
    {
      posknots[1+i]=posknots[i]+vh[i];
    }
  }

  int kx;
  double dhx;
  for(int i=0; i<vectx.size(); i++)
  {
    kx=-1;
    ix = vectx[i];
    if (ix>1) ix=1;
    while(TRUE)
    {
      kx++;
      if(((ix>=posknots[kx]) & (ix<posknots[kx+1])) | (kx>=(nknots-1)))
      {
        break;
      }
    }
    dhx = ix-posknots[kx];
    result[i]=valpha[kx]+vbeta[kx]*dhx+vgamma[kx]*dhx*dhx;
    if(result[i]<0) result[i]=0;
  }
  return(result);
}

// [[Rcpp::export]]
NumericVector VectIntegreSegPol2(const NumericVector&xpars, const NumericVector&vectx,const double&frestrknots)
{
  int n0, nknots;
  double ix;
  n0 = xpars.size()-3;
  n0 = n0/2;
  nknots = n0+1;

  NumericVector result(vectx.size());

  NumericVector valpha(nknots);
  NumericVector vbeta(nknots);
  NumericVector vgamma(nknots);

  valpha[0] = expit(xpars[0]);
  vbeta[0] = xpars[1];
  vgamma[0] = xpars[2];

  NumericVector vh(n0);
  NumericVector posknots(nknots+1);
  posknots[0]=0;
  posknots[nknots]=1;
  if(n0>=1)
  {
    for(int i=0; i<n0; i++)
    {
      vh[i]=xpars[3+2*i];
      vgamma[1+i]= xpars[3+2*i+1];
    }
    vh=Vectexpit(vh);
    for(int i=0; i<vh.size(); i++)
    {
      vh[i]=frestrknots*vh[i];
    }
    for(int i=1; i<nknots; i++)
    {
      vbeta[i]=vbeta[i-1]+2.0*vgamma[i-1]*vh[i-1];
      valpha[i]=valpha[i-1]+vbeta[i-1]*vh[i-1]+vgamma[i-1]*vh[i-1]*vh[i-1];
    }

    for(int i=0; i<vh.size(); i++)
    {
      posknots[1+i]=posknots[i]+vh[i];
    }
  }

  NumericVector Integralknots(nknots);
  Integralknots[0]=0;

  if(n0>=1)
  {
    for(int i=0; i<n0; i++)
    {
      Integralknots[1+i]= Integralknots[i]+valpha[i]*vh[i]+ vbeta[i]*vh[i]*vh[i]/2+vgamma[i]*vh[i]*vh[i]*vh[i]/3;
    }
  }

  int kx;
  double dhx;
  for(int i=0; i<vectx.size(); i++)
  {
    kx=-1;
    ix = vectx[i];
    if (ix>1) ix=1;
    while(TRUE)
    {
      kx++;
      if(((ix>=posknots[kx]) & (ix<posknots[kx+1])) | (kx>=(nknots-1)))
      {
        break;
      }
    }
    dhx = ix-posknots[kx];
    result[i]=Integralknots[kx]+valpha[kx]*dhx+vbeta[kx]*dhx*dhx/2+vgamma[kx]*dhx*dhx*dhx/3;
    if(result[i]<0) result[i]=0;
  }
  return(result);
}

// [[Rcpp::export]]
NumericVector VectIntegrePol2CPos(const NumericVector&xpars, const NumericVector&vectx, const double&frestrknots)
{
  int n0, nknots;
  int kx;
  double dhx, ix;
  n0 = xpars.size()-3;
  n0 = n0/2;
  nknots = n0+1;

  NumericVector result(vectx.size());

  NumericVector valpha(nknots);
  NumericVector vbeta(nknots);
  NumericVector vgamma(nknots);

  valpha[0] = expit(xpars[0]);
  vbeta[0] = xpars[1];
  vgamma[0] = xpars[2];

  NumericVector vh(n0);
  NumericVector posknots(nknots+1);
  posknots[0]=0;
  posknots[nknots]=1;
  if(n0>=1)
  {
    for(int i=0; i<n0; i++)
    {
      vh[i]=xpars[3+2*i];
      vgamma[1+i]= xpars[3+2*i+1];
    }
    vh=Vectexpit(vh);
    for(int i=0; i<vh.size(); i++)
    {
      vh[i]=frestrknots*vh[i];
    }

    for(int i=1; i<nknots; i++)
    {
      vbeta[i]=vbeta[i-1]+2.0*vgamma[i-1]*vh[i-1];
      valpha[i]=valpha[i-1]+vbeta[i-1]*vh[i-1]+vgamma[i-1]*vh[i-1]*vh[i-1];
    }

    for(int i=0; i<vh.size(); i++)
    {
      posknots[1+i]=posknots[i]+vh[i];
    }
  }

  NumericVector Integralknots(nknots);
  Integralknots[0]=0;

  if(n0>=1)
  {
    for(int i=0; i<n0; i++)
    {
      Integralknots[1+i]= Integralknots[i]+IntegrePol2CPos(vh[i],valpha[i],vbeta[i],vgamma[i]);//valpha[i]*vh[i]+ vbeta[i]*vh[i]*vh[i]/2+vgamma[i]*vh[i]*vh[i]*vh[i]/3;
    }
  }

//#pragma omp parallel for num_threads(2)
  for(int i=0; i<vectx.size(); i++)
  {
    kx=-1;
    ix = vectx[i];
    if (ix>1) ix=1;
    while(TRUE)
    {
      kx++;
      if(((ix>=posknots[kx]) & (ix<posknots[kx+1])) | (kx>=(nknots-1)))
      {
        break;
      }
    }
    dhx = ix-posknots[kx];
    result[i]=Integralknots[kx]+IntegrePol2CPos(dhx,valpha[kx],vbeta[kx],vgamma[kx]);
    if(result[i]<0) result[i]=0.0;
  }
  return(result);
}

extern "C" {
SEXP ProjSegPol2PID(SEXP s_xpars, SEXP s_vectx, SEXP s_grecov, SEXP s_uscmax, SEXP s_frestrknots, SEXP s_scfrac)
{
  Rcpp::NumericVector xpars(s_xpars);
  Rcpp::NumericVector vectx(s_vectx);
  Rcpp::NumericVector v_grecov(s_grecov);
  double grecov=v_grecov[0];
  Rcpp::NumericVector v_uscmax(s_uscmax);
  double uscmax=v_uscmax[0];
  Rcpp::NumericVector v_frestrknots(s_frestrknots);
  double frestrknots = v_frestrknots[0];
  Rcpp::NumericVector v_scfrac(s_scfrac);
  double scfrac = v_scfrac[0];
  double xpar, dl, maxy;
  int npoints = vectx.size();

  xpar = uscmax;
  if(xpar<0) xpar = -uscmax;
  if (xpar<1)
  {
    xpar = 1;
  }
  dl = 0.25/xpar;
  maxy = 1.0+max(vectx)/dl;

  npoints = int(maxy);
  int inc_npoints = 2*npoints+1;

  NumericVector all_x(inc_npoints);
  all_x[0]=0.0;
  for(int i=1; i<inc_npoints; i++)
  {
    all_x[i]=all_x[i-1]+dl/2;
  }

  NumericVector vectInci = InciSegPol2PID(all_x,xpars, frestrknots);
  NumericVector vectIntegralInci = VectIntegreSegPol2(xpars,all_x,frestrknots);

  for( int i=0; i< all_x.size(); i++)
  {
    vectIntegralInci[i]=vectIntegralInci[i]+grecov*all_x[i];
  }

  NumericVector all_res(npoints);
  all_res[0]= vectInci[0]/(vectInci[0]+grecov);

  double at0, at1, at2;

  for(int i=1; i<npoints; i++)
  {
    all_res[i]= all_res[i-1]*exp(-scfrac*(vectIntegralInci[2*i]-vectIntegralInci[2*(i-1)]));

    at0 = vectInci[2*i];
    at1 = vectInci[2*i-1]*exp(-scfrac*(vectIntegralInci[2*i]-vectIntegralInci[2*i-1]));
    at2 = vectInci[2*i-2]*exp(-scfrac*(vectIntegralInci[2*i]-vectIntegralInci[2*i-2]));

    all_res[i] = all_res[i]+scfrac*dl/6.0*(at0+4.0*at1+at2);
  }

  NumericVector Result(vectx.size());
  int indx;

  for(int i=0; i<vectx.size(); i++)
  {
    indx = int(vectx[i]/dl);
    Result[i] = all_res[indx];
  }
  return(Result);
}

SEXP LBCheckCalc(SEXP s_XPars, SEXP s_frestrknots)
{
  Rcpp::NumericVector XPars(s_XPars);
  Rcpp::NumericVector v_frestrknots(s_frestrknots);
  double frestrknots = v_frestrknots[0];
  int n0, nknots;
  n0 = XPars.size()-3;
  n0 = n0/2;
  nknots = n0+1;

  NumericVector result(1);

  NumericVector valpha(nknots);
  NumericVector vbeta(nknots);
  NumericVector vgamma(nknots);

  valpha[0] = expit(XPars[0]);
  vbeta[0] = XPars[1];
  vgamma[0] = XPars[2];

  NumericVector vh(n0);
  NumericVector posknots(nknots+1);
  posknots[0]=0;
  posknots[nknots]=1;
  if(n0>=1)
  {
    for(int i=0; i<n0; i++)
    {
      vh[i]=XPars[3+2*i];
      vgamma[1+i]= XPars[3+2*i+1];
    }
    vh=Vectexpit(vh);
    for(int i=0; i<vh.size(); i++)
    {
      vh[i]=frestrknots*vh[i];
    }
    for(int i=1; i<nknots; i++)
    {
      vbeta[i]=vbeta[i-1]+2.0*vgamma[i-1]*vh[i-1];
      valpha[i]=valpha[i-1]+vbeta[i-1]*vh[i-1]+vgamma[i-1]*vh[i-1]*vh[i-1];
    }

    for(int i=0; i<vh.size(); i++)
    {
      posknots[1+i]=posknots[i]+vh[i];
    }
  }

  result[0]=0;
  double dhx, p_res;
  for(int i=1; i<posknots.size(); i++)
  {
    dhx = posknots[i]-posknots[i-1];
    p_res = IntegrePol2CPos(dhx,-valpha[i-1],-vbeta[i-1],-vgamma[i-1]);//

    result[0]+= p_res*p_res;
  }
  return(result);
}

SEXP UBCheckCalc(SEXP s_XPars, SEXP s_frestrknots, SEXP s_upb)
{
  Rcpp::NumericVector XPars(s_XPars);
  Rcpp::NumericVector v_frestrknots(s_frestrknots);
  Rcpp::NumericVector v_upb(s_upb);
  double upb=v_upb[0];
  double frestrknots=v_frestrknots[0];
  int n0, nknots;
  n0 = XPars.size()-3;
  n0 = n0/2;
  nknots = n0+1;

  NumericVector result(1);

  NumericVector valpha(nknots);
  NumericVector vbeta(nknots);
  NumericVector vgamma(nknots);

  valpha[0] = expit(XPars[0])-upb;
  vbeta[0] = XPars[1];
  vgamma[0] = XPars[2];

  NumericVector vh(n0);
  NumericVector posknots(nknots+1);
  posknots[0]=0;
  posknots[nknots]=1;
  if(n0>=1)
  {
    for(int i=0; i<n0; i++)
    {
      vh[i]=XPars[3+2*i];
      vgamma[1+i]= XPars[3+2*i+1];
    }
    vh=Vectexpit(vh);
    for(int i=0; i<vh.size(); i++)
    {
      vh[i]=frestrknots*vh[i];
    }
    for(int i=1; i<nknots; i++)
    {
      vbeta[i]=vbeta[i-1]+2.0*vgamma[i-1]*vh[i-1];
      valpha[i]=valpha[i-1]+vbeta[i-1]*vh[i-1]+vgamma[i-1]*vh[i-1]*vh[i-1];
    }

    for(int i=0; i<vh.size(); i++)
    {
      posknots[1+i]=posknots[i]+vh[i];
    }
  }

  result[0]=0;
  double dhx, p_res;
  for(int i=1; i<posknots.size(); i++)
  {
    dhx = posknots[i]-posknots[i-1];
    p_res = IntegrePol2CPos(dhx,valpha[i-1],vbeta[i-1],vgamma[i-1]);// valpha[kx]+vbeta[kx]*dhx+vgamma[kx]*dhx*dhx;
    result[0]+= p_res*p_res;
  }
  return(result);
}

SEXP InciSegPol2PID_Rv(SEXP s_vectx, SEXP s_xpars, SEXP s_fresknots)
{
  Rcpp::NumericVector vectx(s_vectx);
  Rcpp::NumericVector xpars(s_xpars);
  Rcpp::NumericVector v_fresknots(s_fresknots);
  double fresknots=v_fresknots[0];
  Rcpp::NumericVector result = InciSegPol2PID(vectx, xpars, fresknots);
  return(result);
}

//Derivatives Random effect model on the prevalence
List PartialPrevRandEff(SEXP s_prev, SEXP s_sigma)
{
  NumericVector prev(s_prev);
  NumericVector sigma(s_sigma);
  const int len_prev = prev.size();

  List Result=List::create(Named("Const_dp")=NumericVector(len_prev),Named("dsigma")=NumericVector(len_prev));

  NumericVector Const_dp(len_prev);
  NumericVector dsigma(len_prev);

  const double dx=0.01;
  NumericVector x_disc(101);

  for(int i=0; i<101;i++)
  {
    x_disc[i]=dx*i;
  }

  NumericVector t1x_disc(101);
  NumericVector t2x_disc(101);
  for(int jj=1;jj<100;jj++)
  {
    t1x_disc[jj] = x_disc[jj]/(1-x_disc[jj]);
    t2x_disc[jj] = 1.0/pow(1-x_disc[jj],2.0);
  }

  NumericVector internal_sigma(len_prev);
  if(len_prev==sigma.size())
  {
    for(int i=0;i<len_prev; i++)
    {
      internal_sigma[i]=sigma[i];
    }
  } else
  {
    for(int i=0;i<len_prev; i++)
    {
      internal_sigma[i]=sigma[0];
    }
    if(sigma.size()>1)
    {
      Rcpp::Rcout << "prev and sigma do not have the same length; only the first element of sigma was used" << endl;
    }
  }

  double const scs = sqrt(2.0*pi);
  double eu1, eu2;

#pragma omp parallel for num_threads(2)
  for(int ii=0;ii<len_prev; ii++)
  {
    NumericVector fx_disc(101);
    fx_disc[0] = 0;
    fx_disc[100] = 0;
    for(int jj=1; jj<100; jj++)
    {
      fx_disc[jj] = exp(-0.5*pow(t1x_disc[jj]/internal_sigma[ii],2.0))/scs/internal_sigma[ii]*t2x_disc[jj];
    }

    Const_dp[ii]=0;
    dsigma[ii]=0;
    for(int ll=0; ll<50; ll++)
    {
      eu1= exp(-t1x_disc[2*ll+1]);
      dsigma[ii]+= 4.0*(-1.0/internal_sigma[ii]*(1.0/(eu1-prev[ii]*eu1+prev[ii])+eu1/(1-prev[ii]+prev[ii]*eu1)))*fx_disc[2*ll+1];
      dsigma[ii]+= 4.0*(1/pow(internal_sigma[ii],3.0)*(1.0/(eu1-prev[ii]*eu1+prev[ii])+eu1/(1-prev[ii]+prev[ii]*eu1)))*pow(t1x_disc[2*ll+1],2.0)*fx_disc[2*ll+1];

      eu2= exp(-t1x_disc[2*ll+2]);
      dsigma[ii]+= 2.0*(-1.0/internal_sigma[ii]*(1.0/(eu2-prev[ii]*eu2+prev[ii])+eu2/(1-prev[ii]+prev[ii]*eu2)))*fx_disc[2*ll+2];
      dsigma[ii]+= 2.0*(1/pow(internal_sigma[ii],3.0)*(1.0/(eu2-prev[ii]*eu2+prev[ii])+eu2/(1-prev[ii]+prev[ii]*eu2)))*pow(t1x_disc[2*ll+2],2.0)*fx_disc[2*ll+2];

      Const_dp[ii]+=4.0*(eu1/pow(eu1-prev[ii]*eu1+prev[ii],2.0)+eu1/pow(1-prev[ii]+prev[ii]*eu1,2.0))*fx_disc[2*ll+1];
      Const_dp[ii]+= 2.0*(eu2/pow(eu2-prev[ii]*eu2+prev[ii],2.0)+eu2/pow(1-prev[ii]+prev[ii]*eu2,2.0))*fx_disc[2*ll+2];
    }

    eu1= exp(-t1x_disc[99]);
    dsigma[ii]+= 4.0*(-1.0/internal_sigma[ii]*(1.0/(eu1-prev[ii]*eu1+prev[ii])+eu1/(1-prev[ii]+prev[ii]*eu1)))*fx_disc[99];
    dsigma[ii]+= 4.0*(1/pow(internal_sigma[ii],3.0)*(1.0/(eu1-prev[ii]*eu1+prev[ii])+eu1/(1-prev[ii]+prev[ii]*eu1)))*pow(t1x_disc[99],2.0)*fx_disc[99];

    dsigma[ii] = dsigma[ii]*prev[ii]*dx/3.0;

    Const_dp[ii]+= 4.0*(eu1/pow(eu1-prev[ii]*eu1+prev[ii],2.0)+eu1/pow(1-prev[ii]+prev[ii]*eu1,2.0))*fx_disc[99];
    Const_dp[ii] = Const_dp[ii]*dx/3.0;
  }
  Result["Const_dp"] = Const_dp;
  Result["dsigma"] = dsigma;
  return(Result);
}

//Random effect model on the prevalence
SEXP PrevRandEff(SEXP s_prev, SEXP s_sigma)
{
  Rcpp::NumericVector prev(s_prev);
  Rcpp::NumericVector sigma(s_sigma);
  const double dx=0.01;
  NumericVector x_disc(101);

  NumericVector t1x_disc(101);
  NumericVector t2x_disc(101);
  for(int jj=1;jj<100;jj++)
  {
    t1x_disc[jj] = x_disc[jj]/(1-x_disc[jj]);
    t2x_disc[jj] = 1.0/pow(1-x_disc[jj],2.0);
  }

  const int len_prev = prev.size();
  NumericVector internal_sigma(len_prev);
  if(len_prev==sigma.size())
  {
    for(int i=0;i<len_prev; i++)
    {
      internal_sigma[i]=sigma[i];
    }
  } else
  {
    for(int i=0;i<len_prev; i++)
    {
      internal_sigma[i]=sigma[0];
    }
    if(sigma.size()>1)
    {
      Rcpp::Rcout << "prev and sigma do not have the same length; only the first element of sigma was used" << endl;
    }
  }

  double const scs = sqrt(2.0*pi);
  NumericVector result(len_prev);

  for (int i = 0; i < 101; i++)
  {
      x_disc[i] = dx * i;
  }

#pragma omp parallel for num_threads(2)
  for(int ii=0;ii<len_prev; ii++)
  {
    //NumericVector fx_disc(101);
    std::vector<double> fx_disc(101);
    fx_disc[0] = 0;
    fx_disc[100] = 0;

    for(int jj=1; jj<100; jj++)
    {
      fx_disc[jj] = exp(-0.5*pow(t1x_disc[jj]/internal_sigma[ii],2.0))/scs/internal_sigma[ii]*t2x_disc[jj];
    }

    result[ii]=0;
    for(int ll=0; ll<50; ll++)
    {
      result[ii]+= 4.0*(1.0/(exp(-t1x_disc[2*ll+1])-prev[ii]*exp(-t1x_disc[2*ll+1])+prev[ii])+exp(-t1x_disc[2*ll+1])/(1-prev[ii]+prev[ii]*exp(-t1x_disc[2*ll+1])))*fx_disc[2*ll+1];
      result[ii]+= 2.0*(1.0/(exp(-t1x_disc[2*ll+2])-prev[ii]*exp(-t1x_disc[2*ll+2])+prev[ii])+exp(-t1x_disc[2*ll+2])/(1-prev[ii]+prev[ii]*exp(-t1x_disc[2*ll+2])))*fx_disc[2*ll+2];

    }
    result[ii] += 4.0*(1.0/(exp(-t1x_disc[99])-prev[ii]*exp(-t1x_disc[99])+prev[ii])+exp(-t1x_disc[99])/(1-prev[ii]+prev[ii]*exp(-t1x_disc[99])))*fx_disc[99];

    result[ii] = result[ii]*prev[ii]*dx/3;
  }
  return(result);
}

//Random effect model on the prevalence
SEXP PrevRandEff_2(SEXP s_prev, SEXP s_sigma)
{
  NumericVector prev(s_prev);
  NumericVector sigma(s_sigma);
  const double dx =1.0/15;
  const double x0=-2.0;
  const double x1=-1.8;
  const double x2=-1.6;
  const double x3=-1.4;
  const double x4=-1.2;
  const double x5=-1.0;
  const double x6=-1.0+dx;
  const double x7=-1.0+2*dx;
  const double x8=-1.0+3*dx;
  const double x9=-1.0+4*dx;
  const double x10=-1.0+5*dx;
  const double x11=-1.0+6*dx;
  const double x12=-1.0+7*dx;
  const double x13=-1.0+8*dx;
  const double x14=-1.0+9*dx;
  const double x15=-1.0+10*dx;
  const double x16=-1.0+11*dx;
  const double x17=-1.0+12*dx;
  const double x18=-1.0+13*dx;
  const double x19=-1.0+14*dx;
  const double x20= 0.0;
  const double x21=0.0+dx;
  const double x22=0.0+2*dx;
  const double x23=0.0+3*dx;
  const double x24=0.0+4*dx;
  const double x25=0.0+5*dx;
  const double x26=0.0+6*dx;
  const double x27=0.0+7*dx;
  const double x28=0.0+8*dx;
  const double x29=0.0+9*dx;
  const double x30=0.0+10*dx;
  const double x31=0.0+11*dx;
  const double x32=0.0+12*dx;
  const double x33=0.0+13*dx;
  const double x34=0.0+14*dx;
  const double x35=0.0+15*dx;
  const double x36=1.2;
  const double x37=1.4;
  const double x38=1.6;
  const double x39=1.8;
  const double x40=2.0;

  const int len_prev = prev.size();
  NumericVector internal_sigma(len_prev);
  if(len_prev==sigma.size())
  {
    for(int i=0;i<len_prev; i++)
    {
      internal_sigma[i]=sigma[i];
    }
  } else
  {
    for(int i=0;i<len_prev; i++)
    {
      internal_sigma[i]=sigma[0];
    }
    if(sigma.size()>1)
    {
      Rcpp::Rcout << "prev and sigma do not have the same length; only the first element of sigma was used" << endl;
    }
  }

  double const scs = sqrt(2.0*pi);

  double const h0=(2.0-1.0)/5;
  double const h1=(1-(-1.0))/30;

  NumericVector result(len_prev);

#pragma omp parallel for num_threads(2)
  for(int ii=0;ii<len_prev; ii++)
  {
    //NumericVector fxi(41);
    std::vector<double> fxi(41);
    fxi[0]= exp(-0.5*pow(x0/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[1]= exp(-0.5*pow(x1/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[2]= exp(-0.5*pow(x2/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[3]= exp(-0.5*pow(x3/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[4]= exp(-0.5*pow(x4/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[5]= exp(-0.5*pow(x5/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[6]= exp(-0.5*pow(x6/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[7]= exp(-0.5*pow(x7/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[8]= exp(-0.5*pow(x8/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[10]= exp(-0.5*pow(x9/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[11]= exp(-0.5*pow(x11/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[12]= exp(-0.5*pow(x12/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[13]= exp(-0.5*pow(x13/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[14]= exp(-0.5*pow(x14/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[15]= exp(-0.5*pow(x15/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[16]= exp(-0.5*pow(x16/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[17]= exp(-0.5*pow(x17/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[18]= exp(-0.5*pow(x18/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[19]= exp(-0.5*pow(x19/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[20]= exp(-0.5*pow(x20/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[21]= exp(-0.5*pow(x21/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[22]= exp(-0.5*pow(x22/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[23]= exp(-0.5*pow(x23/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[24]= exp(-0.5*pow(x24/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[25]= exp(-0.5*pow(x25/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[26]= exp(-0.5*pow(x26/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[27]= exp(-0.5*pow(x27/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[28]= exp(-0.5*pow(x28/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[29]= exp(-0.5*pow(x29/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[30]= exp(-0.5*pow(x30/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[31]= exp(-0.5*pow(x31/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[32]= exp(-0.5*pow(x32/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[33]= exp(-0.5*pow(x33/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[34]= exp(-0.5*pow(x34/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[35]= exp(-0.5*pow(x35/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[36]= exp(-0.5*pow(x36/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[37]= exp(-0.5*pow(x37/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[38]= exp(-0.5*pow(x38/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[39]= exp(-0.5*pow(x39/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    fxi[40]= exp(-0.5*pow(x40/internal_sigma[ii],2.0))/scs/internal_sigma[ii];

    double fx4_5 = exp(-0.5*pow(-1.1/internal_sigma[ii],2.0))/scs/internal_sigma[ii];
    double fx35_5 = exp(-0.5*pow(1.1/internal_sigma[ii],2.0))/scs/internal_sigma[ii];

    double approx_integral1;
    double approx_integral2;
    double approx_integral3;

    approx_integral1 = exp(x0)/(1-prev[ii]+prev[ii]*exp(x0))*fxi[0]+4*exp(x1)/(1-prev[ii]+prev[ii]*exp(x1))*fxi[1]+2*exp(x2)/(1-prev[ii]+prev[ii]*exp(x2))*fxi[2];
    approx_integral1 += 4*exp(x3)/(1-prev[ii]+prev[ii]*exp(x3))*fxi[3]+exp(x4)/(1-prev[ii]+prev[ii]*exp(x4))*fxi[4];
    approx_integral1 = approx_integral1*h0/3;

    approx_integral1 = approx_integral1+(x5-x4)/6*(exp(x4)/(1-prev[ii]+prev[ii]*exp(x4))*fxi[4]+4*exp(-1.1)/(1-prev[ii]+prev[ii]*exp(-1.1))*fx4_5+exp(x5)/(1-prev[ii]+prev[ii]*exp(x5))*fxi[5]);

    approx_integral2 = exp(x5)/(1-prev[ii]+prev[ii]*exp(x5))*fxi[5]+4*exp(x6)/(1-prev[ii]+prev[ii]*exp(x6))*fxi[6];
    approx_integral2 += 2*exp(x7)/(1-prev[ii]+prev[ii]*exp(x7))*fxi[7]+4*exp(x8)/(1-prev[ii]+prev[ii]*exp(x8))*fxi[8];
    approx_integral2 += 2*exp(x9)/(1-prev[ii]+prev[ii]*exp(x9))*fxi[9]+4*exp(x10)/(1-prev[ii]+prev[ii]*exp(x10))*fxi[10];
    approx_integral2 += 2*exp(x11)/(1-prev[ii]+prev[ii]*exp(x11))*fxi[11]+4*exp(x12)/(1-prev[ii]+prev[ii]*exp(x12))*fxi[12];
    approx_integral2 += 2*exp(x13)/(1-prev[ii]+prev[ii]*exp(x13))*fxi[13]+4*exp(x14)/(1-prev[ii]+prev[ii]*exp(x14))*fxi[14];
    approx_integral2 += 2*exp(x15)/(1-prev[ii]+prev[ii]*exp(x15))*fxi[15]+4*exp(x16)/(1-prev[ii]+prev[ii]*exp(x16))*fxi[16];
    approx_integral2 += 2*exp(x17)/(1-prev[ii]+prev[ii]*exp(x17))*fxi[17]+4*exp(x18)/(1-prev[ii]+prev[ii]*exp(x18))*fxi[18];
    approx_integral2 += 2*exp(x19)/(1-prev[ii]+prev[ii]*exp(x19))*fxi[19]+4*exp(x20)/(1-prev[ii]+prev[ii]*exp(x20))*fxi[20];
    approx_integral2 += 2*exp(x21)/(1-prev[ii]+prev[ii]*exp(x21))*fxi[21]+4*exp(x22)/(1-prev[ii]+prev[ii]*exp(x22))*fxi[22];
    approx_integral2 += 2*exp(x23)/(1-prev[ii]+prev[ii]*exp(x23))*fxi[23]+4*exp(x24)/(1-prev[ii]+prev[ii]*exp(x24))*fxi[24];
    approx_integral2 += 2*exp(x25)/(1-prev[ii]+prev[ii]*exp(x25))*fxi[25]+4*exp(x26)/(1-prev[ii]+prev[ii]*exp(x26))*fxi[26];
    approx_integral2 += 2*exp(x27)/(1-prev[ii]+prev[ii]*exp(x27))*fxi[27]+4*exp(x28)/(1-prev[ii]+prev[ii]*exp(x28))*fxi[28];
    approx_integral2 += 2*exp(x29)/(1-prev[ii]+prev[ii]*exp(x29))*fxi[29]+4*exp(x30)/(1-prev[ii]+prev[ii]*exp(x30))*fxi[30];
    approx_integral2 += 2*exp(x31)/(1-prev[ii]+prev[ii]*exp(x31))*fxi[31]+4*exp(x32)/(1-prev[ii]+prev[ii]*exp(x32))*fxi[32];
    approx_integral2 += 2*exp(x33)/(1-prev[ii]+prev[ii]*exp(x33))*fxi[33]+4*exp(x34)/(1-prev[ii]+prev[ii]*exp(x34))*fxi[34]+exp(x35)/(1-prev[ii]+prev[ii]*exp(x35))*fxi[35];
    approx_integral2 = approx_integral2*h1/3;

    approx_integral3=(x36-x35)/6*(exp(x35)/(1-prev[ii]+prev[ii]*exp(x35))*fxi[35]+4*exp(1.1)/(1-prev[ii]+prev[ii]*exp(1.1))*fx35_5+exp(x36)/(1-prev[ii]+prev[ii]*exp(x36))*fxi[36]);


    result[ii]= exp(x36)/(1-prev[ii]+prev[ii]*exp(x36))*fxi[36]+4*exp(x37)/(1-prev[ii]+prev[ii]*exp(x37))*fxi[37];
    result[ii] += 2*exp(x38)/(1-prev[ii]+prev[ii]*exp(x38))*fxi[38]+4*exp(x39)/(1-prev[ii]+prev[ii]*exp(x39))*fxi[39];
    result[ii] += exp(x40)/(1-prev[ii]+prev[ii]*exp(x40))*fxi[40];
    result[ii] = result[ii]*h0/3;

    result[ii]+= approx_integral3+approx_integral2+approx_integral1;
    result[ii]= prev[ii]*result[ii];
  }
  return(result);
}

SEXP logit(SEXP s_x)
{
    Rcpp::NumericVector x(s_x);
    int n = x.size();
    NumericVector result(n);
    for (int i = 0; i < n; i++)
    {
        result[i] = log(x[i] / (1 - x[i]));
    }
    return(result);
}

}






