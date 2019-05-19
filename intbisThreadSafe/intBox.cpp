#include "intBox.h"


intBox::intBox(void)
  {
  nullpointer();

  }


intBox::~intBox(void)
  {
  flushMemory();
  }

// *****************************************************************************

void intBox::flushMemory()
  {
  if(X)
    {
    delete[] X; X=0;
    }
  if(Xpt)
    {
    delete[] Xpt; Xpt = 0;
    }

  if(f) delete[] f;
  if(df) delete[] df;
  if(fpt) delete[] fpt;
  if(dfpt) delete[] dfpt;

  delete ts;
  delete[] vx;
  delete[] vf;
  nullpointer();
  }

void intBox::nullpointer()
  {
  X   = 0;
  Xpt = 0;

  f   = 0;
  df  = 0;
  fpt = 0;
  dfpt= 0;
  
  szX = 0;
  
  fobjpt = 1.0e+100;

  ts = 0;

  vx = 0;
  vf = 0;
  }

void intBox::initialize(long n)
  {
  X = new Interval[n];
  Xpt = new double[n];
  f = new Interval[n];
  df = new Interval[n*n];
  fpt = new double[n];
  dfpt = new double[n*n];

  ts = new Taylor();
  vx = new IntervalD[n];
  vf = new IntervalD[n];
  for(int i=0; i<n; ++i) 
    {
    vx[i].setDependency(INDEPENDENT);
    vf[i].setDependency(DEPENDENT);
    vx[i].resize(n);
    vf[i].resize(n);
    }

  szX = n;
  }

// *****************************************************************************

long intBox::operator&(const intBox& copyfrom) const
  {
  long nx, ret;
  
  nx = copyfrom.szX;
  if(nx!=szX)
    {
    ret = 0;
    }
  else
    {
    long i;
    ret = 1;
    for(i=0; i<nx; ++i)
      {
      if(!(X[i] & copyfrom.X[i]))
        {
        ret = 0; break;
        }
      }    
    }
  return ret;
  }

// *****************************************************************************

void intBox::setValue(long n, Interval *x, double *xpt)
  {
  long i;  
  for(i=0; i<n; ++i) 
    X[i] = x[i];
  if(xpt)
    {
    for(i=0; i<n; ++i) Xpt[i] = xpt[i];
    }
  else
    {
    for(i=0; i<n; ++i) Xpt[i] = mid(X[i]);
    }
  }

// *****************************************************************************

void intBox::getValue(long n, Interval *x, double *xpt)
  {
  for(long i=0;i<n;++i)
    {
    x[i] = X[i];
    xpt[i] = Xpt[i];
    }
  }

void intBox::setValue(long i, double Lb, double Ub)
  {
  X[i].setInterval(Lb, Ub);
  Xpt[i] = 0.5*(Lb+Ub);
  }

// *****************************************************************************

intBox& intBox::operator=(const intBox& copyfrom)
  {
  if(szX!=copyfrom.szX)
    {
    flushMemory();
    szX = copyfrom.szX;
    initialize(szX);
    }
  setValue(szX, copyfrom.X, copyfrom.Xpt);
  return *this;
  }

