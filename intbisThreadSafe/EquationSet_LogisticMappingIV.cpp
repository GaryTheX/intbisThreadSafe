#include "EquationSet_LogisticMappingIV.h"


EquationSet_LogisticMappingIV::EquationSet_LogisticMappingIV(void)
  {
  pX0 = 0;
  nx0 = 0;
  }


EquationSet_LogisticMappingIV::~EquationSet_LogisticMappingIV(void)
  {
  delete[] pX0;
  }

void EquationSet_LogisticMappingIV::setInitialValue(long nx, double *x0)
  {
  if(nx0<nx)
    {
    delete[] pX0;
    pX0 = new double[nx];
    }
  nx0 = nx;
  for(int i=0; i<nx0; ++i) pX0[i] = x0[i];
  }

void EquationSet_LogisticMappingIV::fun(intBox *pBox, long idx)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;
// Feigenbaum example
  long i, n = pBox->szX;
  if(idx<0)
    {
    f[0] = -param[0]*pX0[0]*pX0[0] + param[0]*pX0[0]   - x[0];
    for(i=1; i<n; ++i)
      f[i] = -param[0]*power(x[i-1]  ,2) + param[0]*x[i-1]   - x[i];
    }
  else
    {
    i = idx;
    if(i>0)
      f[i] = -param[0]*power(x[i-1]  ,2) + param[0]*x[i-1]   - x[i];
    else
      f[i] = -param[0]*pX0[0]*pX0[0] + param[0]*pX0[0]   - x[0];
    }
  }

void EquationSet_LogisticMappingIV::jac(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Interval *df= pBox->df;
  // Feigenbaum example
  long i, j, n = pBox->szX;
  if(n==1)
    {
    df[0] = - 1.0;
    return;
    }

  if(idxEq<0)
    {
    for(j=0; j<n; ++j) df[j] = 0.0;
    df[0] = -1.0;
    for(i=1; i<n; ++i)
      {
      //f[i] = -r_Fgb*power(x[i]  ,2) + r_Fgb*x[i]   - x[i+1];
      for(j=0; j<n; ++j) df[i*n+j] = 0.0;
      df[i*n+i-1] = -param[0]*2.0*x[i-1] + param[0];
      df[i*n+i]   = -1.0;
      }
    }
  else
    {
    if(idxEq>0)
      {
      if(idxVr==idxEq-1)
        df[idxEq*n+idxVr] = -param[0]*2.0*x[idxVr] + param[0];
      else if(idxVr==idxEq)
        df[idxEq*n+idxVr] = -1.0;
      else
        df[idxEq*n+idxVr] = 0.0;
      }
    else
      {
      if(idxVr==0)
        df[0] = -1.0;
      else
        df[idxVr] = 0.0;
      }
    }
  }

void EquationSet_LogisticMappingIV::tens(intBox *pBox, long idxEq, long idxVr)
  {
  Taylor *ts = pBox->ts;

// Feigenbaum example
  if(idxVr==idxEq-1)
    ts->Taylor_fatParam = -2.0*param[0];
  else
    ts->Taylor_fatParam = 0.0;
  ts->Remain_fatParam = 0.0;
  ts->isRemainConst = true;
  }