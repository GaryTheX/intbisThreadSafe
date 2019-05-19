#include "EquationSet_Broyden.h"


EquationSet_Broyden::EquationSet_Broyden(void)
  {
  }


EquationSet_Broyden::~EquationSet_Broyden(void)
  {
  }

char* EquationSet_Broyden::getName()
  {
  if(icase==1)
    return "Broyden_1";
  else if(icase==2)
    return "Broyden_2";  
  else if (icase == 3)
    return "Broyden_3";
  else
    return "";
  }

void EquationSet_Broyden::setBenchMarkID(int id)
  {
  icase = id;
  switch(icase)
    {
    case 1:
      {
      neq = 10;
      break;
      }
    case 2:
      {
      neq = 20;
      break;
      }
    case 3:
    case 4:
      {
      neq = 64;
      break;
      }
    default:
      break;
    }
  }
bool EquationSet_Broyden::obt_buildInBound(intBox *pBox)
  {
  pBox->initialize(neq);
  int i;
  if(icase==1 || icase==2 || icase==3)
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(-1.0, 1.0);
    }
  else if(icase==4)
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(-1.0e8, 1.0e8);
    }

  return true;
  }

void EquationSet_Broyden::fun(intBox *pBox, long idx)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;

  long i, n = pBox->szX, is=0, ie = n, j, jmin, jmax;

  if(idx>=0)
    {
    is = idx; ie = is+1;
    }

  for(i=is; i<ie; ++i)
    {
    f[i] = x[i]*(2.0+5.0*power(x[i],2)) + 1.0;
    jmin = (1>i-5)?1:i-5;
    jmax = (n-1<i+1)?n-1:i+1;
    for(j=jmin; j<jmax; ++j)
      {
      if(j==i) continue;
      f[i] += x[j]*(1.0+x[j]);
      }
    }
  }

void EquationSet_Broyden::jac(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *X = pBox->X;
  Interval *df= pBox->df;

  IntervalD *x = pBox->vx;
  IntervalD *f = pBox->vf;

  long n = pBox->szX, is=0, ie = n, k, kmin, kmax;
  
  long i, j, js = 0, je = n;

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
    f[i] = x[i]*(2.0+5.0*x[i]*x[i]) + 1.0;
    kmin = (1>i-5)?1:i-5;
    kmax = (n-1<i+1)?n-1:i+1;
    for(k=kmin; k<kmax; ++k)
      {
      if(k==i) continue;
      f[i] += x[k]*(1.0+x[k]);
      }
    for(j=js; j<je; ++j)
      df[i*n+j] = f[i].getDerivative(j);
    }
  }