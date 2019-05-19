#include "EquationSet_Economic.h"


EquationSet_Economic::EquationSet_Economic(void)
  {
  pc = 0;
  }


EquationSet_Economic::~EquationSet_Economic(void)
  {
  delete[] pc;
  }

void EquationSet_Economic::setBenchMarkID(int id)
  {
  icase = id;
  switch(icase)
    {
    case 1:
      {
      neq = 4;
      break;
      }
    case 2:
      {
      neq = 5;
      break;
      }
    case 3:
      {
      neq = 6;
      break;
      }
    case 4:
      {
      neq = 7;
      break;
      }
    default:
      break;
    }
  }

bool EquationSet_Economic::obt_buildInBound(intBox *pBox)
  {
  pBox->initialize(neq);
  int i;
  
  for(i=0; i<neq; ++i)
    //pBox->X[i] = 1.0;
    pBox->X[i].setInterval(-1.0e2, 1.0e2);

  pc = new double[neq-1];
  for(i=0; i<neq-1; ++i) pc[i] = -1.0;

  return true;
  }

void EquationSet_Economic::fun(intBox *pBox, long idx)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;

  long i, n = pBox->szX, is=0, ie = n, neq_1 = n - 1;

  if(idx>=0)
    {
    is = idx; ie = is+1;
    }

  for(i=is; i<ie; ++i)
    {
    if(i==neq_1)
      {
      f[i] = 1.0;
      for(long l=0; l<neq_1-1; ++l)
        f[i] += x[l];
      }
    else
      {
      f[i] = x[i];
      for(long l=0; l<neq_1-i; ++l)
        f[i] += x[l]*x[l+i];
      f[i] *= x[neq_1];
      f[i] -= pc[i];
      }
    }
  }

void EquationSet_Economic::jacDigit(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;
  Interval *df= pBox->df;

  long i, l, n = pBox->szX, is=0, ie = n, js=0, je=n, neq_1 = n - 1;

  if(idxEq>=0)
    {
    is = idxEq; ie = is+1;
    }
  if(idxVr>=0)
    {
    js = idxVr; je = js+1;
    }

  for(i=is; i<ie; ++i)
    {
    if(i==neq_1)
      {
      for(l=0; l<neq_1-1; ++l)
        df[i*n+l] = 1.0;
      df[i*n+neq_1] = 0.0;
      }
    else
      {
      for(l=0; l<neq_1; ++l)
        {
        if(i==l)
          df[i*n+l] = 1.0;
        else
          df[i*n+l] = 0.0;
        for(long ii=0; ii<neq_1-i; ++ii)
          {
          if(l==ii)
            df[i*n+l] += x[i+ii];
          if(l==ii+i)
            df[i*n+l] += x[ii];
          }
        df[i*n+l] *= x[neq_1];
        }
      df[i*n+neq_1] = x[i];
      for(l=0; l<neq_1-i; ++l)
        df[i*n+neq_1] += x[l]*x[l+i];
      }
    }
  }

void EquationSet_Economic::jac(intBox *pBox, long idxEq, long idxVr)
  {
  return jacDigit(pBox, idxEq, idxVr);

  Interval *X = pBox->X;
  Interval *F = pBox->f;
  Interval *df= pBox->df;

  IntervalD *x = pBox->vx;
  IntervalD *f = pBox->vf;

  long i, n = pBox->szX, is=0, ie = n, neq_1 = n-1;
  long j, js = 0, je = n;

  for(i=0; i<n; ++i)
    {
    x[i].setValue(X[i]);
    x[i].setIdx_Var(i);
    }

  if(idxEq>=0)
    {
    is = idxEq; ie = is+1;
    }
  if(idxVr>=0)
    {
    js = idxVr; je = js+1;
    }
  for(i=is; i<ie; ++i)
    {
    if(i==neq_1)
      {
      f[i] = x[0];
      for(long l=1; l<neq_1-1; ++l)
        f[i] += x[l];
      }
    else
      {
      f[i] = x[i];
      for(long l=0; l<neq_1-i; ++l)
        f[i] += x[l]*x[l+i];
      f[i] *= x[neq_1];
      }
    for(j=js; j<je; ++j)
      df[i*n+j] = f[i].getDerivative(j);
    }
  }