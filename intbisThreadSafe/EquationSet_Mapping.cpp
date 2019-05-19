#include "EquationSet_Mapping.h"

EquationSet_Mapping::EquationSet_Mapping(void)
  {
  nparam = 0;
  param  = 0;
  }

EquationSet_Mapping::~EquationSet_Mapping(void)
  {
  flushMemory();
  }

void EquationSet_Mapping::flushMemory()
  {
  delete[] param; param = 0;
  nparam = 0;
  }
char* EquationSet_Mapping::getName()
  {
  if(icase==1)
    return "Feigenbaum_1";
  else if(icase==2)
    return "Feigenbaum_2";
  else if(icase==3)
    return "Feigenbaum_3";
  else
    return "";
  }

void EquationSet_Mapping::setParameters(long np, double *p)
  {
  if(nparam < np)
    {
    flushMemory();
    nparam = np;
    param = new double[np];
    }
  for(int i=0; i<np; ++i) param[i] = p[i];
  }

void EquationSet_Mapping::setBenchMarkID(int id)
  {
  icase = id;
  switch(icase)
    {
    case 1:
      {
      neq = 3;
      nparam = 1;
      if(!param) param = new double[nparam];
      param[0] = 3.84;
      break;
      }
    case 2:
      {
      neq = 5;
      nparam = 1;
      if(!param) param = new double[nparam];
      param[0] = 3.84;
      break;
      }
    case 3:
      {
      neq = 13;
      nparam = 1;
      if(!param) param = new double[nparam];
      param[0] = 3.84;
      break;
      }
    case 11:
      {
      neq = 18;
      nparam = 2;
      if(!param) param = new double[nparam];
      param[0] = 1.4; param[1] = 0.3;
      break;
      }
    default:
      break;
    }
  }

bool EquationSet_Mapping::obt_buildInBound(intBox *pBox)
  {
  pBox->initialize(neq);
  int i;
  if(icase==1 || icase==2 || icase==3)
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(0.0, 1.0e2);
    }
  else if(icase==11)
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(-1.6, 1.6);
    }

  return true;
  }

void EquationSet_Mapping::fun(intBox *pBox, long idx)
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
      f[neq_1] = -param[0]*power(x[neq_1],2) + param[0]*x[neq_1] - x[0];
      }
    else
      {
      f[i] = -param[0]*power(x[i]  ,2) + param[0]*x[i]   - x[i+1];
      }
    }
  }

void EquationSet_Mapping::jac(intBox *pBox, long idxEq, long idxVr)
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
      f[neq_1] = x[neq_1]*param[0] - x[0] - param[0]*power(x[neq_1],2);
      }
    else
      {
      f[i] = -param[0]*power(x[i]  ,2) + param[0]*x[i]   - x[i+1];
      }
    for(j=js; j<je; ++j)
      df[i*n+j] = f[i].getDerivative(j);
    }
  }

void EquationSet_Mapping::jacDigit(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Interval *df= pBox->df; 
  long i, n = pBox->szX, is=0, ie = n, neq_1 = n-1;
  long j, js = 0, je = n;

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
      for(j=js; j<je; ++j)
        df[i*n+j] = 0.0;
      df[i*n] = -1.0;
      df[i*n+i] = param[0] * (1.0-2.0*x[i]);
      //f[neq_1] = x[neq_1]*param[0] - x[0] - param[0]*power(x[neq_1],2);
      }
    else
      {
      for(j=js; j<je; ++j) df[i*n+j] = 0.0;
      df[i*n+i]   = param[0]*(1.0-2.0*x[i]);
      df[i*n+i+1] = -1.0;
      //f[i] = -param[0]*power(x[i]  ,2) + param[0]*x[i]   - x[i+1];
      }
    }
  }


void EquationSet_Mapping::tens(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Taylor *ts = pBox->ts;
  Interval *f = pBox->f;

// Feigenbaum example
  if(idxVr==idxEq)
    ts->Taylor_fatParam = -2.0*param[0];
  else
    ts->Taylor_fatParam = 0.0;
  ts->Remain_fatParam = 0.0;
  ts->isRemainConst = true;
  }

