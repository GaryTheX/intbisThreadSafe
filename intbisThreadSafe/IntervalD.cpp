#include "IntervalD.h"
#include <math.h>
#include <ostream>

extern "C"
  {
  extern void ROUNDDOWN(double& left);
  extern void ROUNDUP(double& right);
  extern double gettiny();
  }
const double huge_val = 1.0e20;
const double negInfy  = -230.0;
const double ln10 = 2.30258509299405;

#define doDyn

IntervalD::IntervalD()
  {
  Ndim = 0;
  DvDx = 0;
  appended = false;
  idx_Var = 0;
  idx_Der = -1;
  adType = INDEPENDENT;
#ifdef doDyn
  aPool = intervalDPool::getInstance();
#endif
  }

IntervalD::IntervalD(int n)
{
  Ndim = n;
  DvDx = new Interval[n];
  for (int i = 0; i < n; ++i) 
    DvDx[i].setInterval(0.0, 0.0);
  appended = false;
  idx_Var = 0;
  idx_Der = -1;
  adType = INDEPENDENT;
#ifdef doDyn
  aPool = intervalDPool::getInstance();
#endif
}
// *********************************************

IntervalD::IntervalD(double left, double right)
  {
  VF.setInterval(left, right);
  appended = false;
  Ndim = 0;
  DvDx = 0;
  adType = INDEPENDENT;
#ifdef doDyn
  aPool = intervalDPool::getInstance();
#endif
  }

IntervalD::~IntervalD(void)
  {
  }

void IntervalD::killPointers()
  {
  delete[] DvDx; DvDx=0;
  intervalDPool::deleteInstance();
  }

// *********************************************

void IntervalD::resize(long N)
  {
  if(Ndim<N)
    {
    delete[] DvDx;
    Ndim = N;
    DvDx = new Interval[Ndim];       
    }
  }

// *********************************************

void IntervalD::reappend(intervalDPool *outsidePool)
  {
  //aPool = outsidePool;
  }

void IntervalD::setValue(Interval x)
  {
  VF = x;
  }

// *********************************************

void IntervalD::setIdx_Var(long idx)
  {
  idx_Var = idx;
  DvDx[idx_Var] = 1.0;
  }

// *********************************************

IntervalD& IntervalD::operator = (const IntervalD &x)
{
  VF = x.VF;
  
  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  resize(n);
  for(int i=0;i<n;i++)
    DvDx[i]=x.DvDx[i];
#ifdef doDyn
  (const_cast<IntervalD*>(&x))->getIntervalDPool()->setFree();
#endif
  return *this;
}

IntervalD &IntervalD::operator = (double x)
{
  VF.setInterval(x,x);
  return *this;
}

// *********************************************

IntervalD &IntervalD::operator = (long x)
{
  VF.setInterval(x,x);
  return *this;
}
// *********************************************
#ifndef doDyn
IntervalD operator + (const IntervalD& x, const IntervalD& y)
  {
  IntervalD tmp;
  tmp.VF = x.VF + y.VF;
  long n  = (const_cast<IntervalD*>(&x))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=x.DvDx[i]+y.DvDx[i];
  return tmp;
  }
 
// *********************************************

IntervalD operator + (const IntervalD& x, double y)
  {
  IntervalD tmp;
  tmp.VF = x.VF + y;
  
  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=x.DvDx[i];
  return tmp;
  }
// *********************************************

IntervalD operator + (double x, const IntervalD& y)
  {
  IntervalD tmp;
  tmp.VF = x + y.VF;

  long n = (const_cast<IntervalD*>(&y))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=y.DvDx[i];
  return tmp;
  }

// *********************************************

IntervalD operator - (const IntervalD& x, const IntervalD& y)
  {
  IntervalD tmp;
  tmp.VF = x.VF - y.VF;

  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=x.DvDx[i]-y.DvDx[i];
  return tmp;
  }

// *********************************************

IntervalD operator - (double x, const IntervalD& y)
  {
  IntervalD tmp;
  tmp.VF = x - y.VF;

  long n = (const_cast<IntervalD*>(&y))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=0.0-y.DvDx[i];
  return tmp;
  }
// *********************************************

IntervalD operator - (const IntervalD& x, double y)
  {
  IntervalD tmp;

  tmp.VF = x.VF - y;
  
  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=x.DvDx[i];
  return tmp;
  }

// *********************************************

IntervalD operator * (const IntervalD& x, const IntervalD& y)
  {
  IntervalD tmp;

  tmp.VF = x.VF * y.VF;

  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=x.DvDx[i]*y.VF+x.VF*y.DvDx[i];
  return tmp;
  }
 
// *********************************************

IntervalD operator * (double x, const IntervalD& y)
  {
  IntervalD tmp;
 
  tmp.VF = x * y.VF;

  long n = (const_cast<IntervalD*>(&y))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=x*y.DvDx[i];
  return tmp;
  }

// *********************************************

IntervalD operator * (const IntervalD& x, double y)
  {
  IntervalD tmp;
  
  tmp.VF = x.VF * y;

  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=y * x.DvDx[i];
  return tmp;
  }

// *********************************************

IntervalD operator / (const IntervalD& x, const IntervalD& y)
  {
  IntervalD tmp;
 
  tmp.VF = x.VF / y.VF;

  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=(x.DvDx[i] - tmp.VF * y.DvDx[i]/y.VF)/y.VF;
  return tmp;
  }

// *********************************************

IntervalD operator / (double x, const IntervalD& y)
  {
  IntervalD tmp;
  
  tmp.VF = x / y.VF;

  long n = (const_cast<IntervalD*>(&y))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=-x/(y.VF*y.VF) *y.DvDx[i];
  return tmp;
  }
// *********************************************

IntervalD operator / (const IntervalD& x, double y)
  {
  IntervalD tmp;
  
  tmp.VF = x.VF / y;

  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  tmp.resize(n);
  for(int i=0;i<n;i++)
    tmp.DvDx[i]=x.DvDx[i]/y;
  return tmp;
  }

void IntervalD::operator += (const IntervalD &x)
{
  VF+=x.VF;
long n = (const_cast<IntervalD*>(&x))->getDimension();
for (int i = 0; i < n; i++)
  DvDx[i] += x.DvDx[i];
}

void IntervalD::operator -= (const IntervalD &x)
{
  VF -= x.VF;
  long n = (const_cast<IntervalD*>(&x))->getDimension();
  for (int i = 0; i < n; i++)
    DvDx[i] -= x.DvDx[i];
}

void IntervalD::operator *= (const IntervalD &y)
{
  VF = VF * y.VF;
  long n = (const_cast<IntervalD*>(&y))->getDimension();
  for (int i = 0; i < n; i++)
    DvDx[i] = VF * y.DvDx[i] + DvDx[i] * y.VF;
}

//------------------------------------------
IntervalD ilog(const IntervalD &x)
{
  IntervalD tmp;

  tmp.VF = ilog(x.VF);
  long n = (const_cast<IntervalD*>(&x))->getDimension();
  tmp.resize(n);
  for (int i = 0; i < n; i++)
    tmp.DvDx[i] = x.DvDx[i] / x.VF;
  return tmp;
}

//------------------------------------------
IntervalD ilog10(const IntervalD &x)
{
  IntervalD tmp;

  tmp.VF = ilog10(x.VF);
  long n = (const_cast<IntervalD*>(&x))->getDimension();
  tmp.resize(n);

  for (int i = 0; i < n; i++)
    tmp.DvDx[i] = x.DvDx[i] / (x.VF*ln10);
  return tmp;
}


IntervalD iexp(const IntervalD &x)
{
  IntervalD tmp;

  tmp.VF = iexp(x.VF);

  long n = (const_cast<IntervalD*>(&x))->getDimension();
  tmp.resize(n);

  for (int i = 0; i < n; i++)
    tmp.DvDx[i] = tmp.VF * x.DvDx[i];
  return tmp;
}

IntervalD pow(const IntervalD &x, double p)
{
  IntervalD tmp;

  tmp.VF = iipow(x.VF, p);
  Interval tmp2 = p * iipow(x.VF, p - 1.0);

  long n = (const_cast<IntervalD*>(&x))->getDimension();
  tmp.resize(n);

  for (int i = 0; i < n; i++)
    tmp.DvDx[i] = tmp2 * x.DvDx[i];
  return tmp;
}

IntervalD power(const IntervalD &x, long p)
{
  IntervalD tmp;

  tmp.VF = power(x.VF, p);
  Interval tmp2 = p * power(x.VF, p - 1);

  long n = (const_cast<IntervalD*>(&x))->getDimension();
  tmp.resize(n);

  for (int i = 0; i < n; i++)
    tmp.DvDx[i] = tmp2 * x.DvDx[i];
  return tmp;
}
#else
// *********************************************

IntervalD operator + (const IntervalD& x, const IntervalD& y)
{
  IntervalD *tmp = 0;
  int idx;
  if((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);
    tmp->VF = x.VF + y.VF;
    long n  = (const_cast<IntervalD*>(&x))->getDimension(); 
    //tmp->resize(n);
    for(int i=0;i<n;i++)
      tmp->DvDx[i]=x.DvDx[i]+y.DvDx[i];
  }

  return *tmp;
  }
// *********************************************

IntervalD operator + (const IntervalD& x, double y)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->VF = x.VF + y;
    tmp->reappend(x.aPool);
    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = x.DvDx[i];
  }
  return *tmp;
  }
// *********************************************

IntervalD operator + (double x, const IntervalD& y)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&y))->getIntervalDPool()->obtainItem(y, tmp, idx))
  {
    tmp->reappend(y.aPool);
    tmp->VF = x + y.VF;

    long n = (const_cast<IntervalD*>(&y))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = y.DvDx[i];
  }
  return *tmp;
  }

// *********************************************

IntervalD operator - (const IntervalD& x, const IntervalD& y)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);
    tmp->VF = x.VF - y.VF;

    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = x.DvDx[i] - y.DvDx[i];
  }
  return *tmp;
  }

// *********************************************

IntervalD operator - (double x, const IntervalD& y)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&y))->getIntervalDPool()->obtainItem(y, tmp, idx))
  {
    tmp->reappend(y.aPool);
    tmp->VF = x - y.VF;

    long n = (const_cast<IntervalD*>(&y))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = 0.0 - y.DvDx[i];
  }
  return *tmp;
  }
// *********************************************

IntervalD operator - (const IntervalD& x, double y)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->VF = x.VF - y;
    tmp->reappend(x.aPool);
    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = x.DvDx[i];
  }
  return *tmp;
  }

// *********************************************

IntervalD operator * (const IntervalD& x, const IntervalD& y)
  {
  IntervalD *tmp=0;
  int idx;
  if( (const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
    {
      tmp->reappend(x.aPool);
      tmp->VF = x.VF * y.VF;

      long n = (const_cast<IntervalD*>(&x))->getDimension();
      //tmp->resize(n);
      for (int i = 0; i < n; i++)
        tmp->DvDx[i] = x.DvDx[i] * y.VF + x.VF*y.DvDx[i];
    
    }  
  return *tmp;
  }

// *********************************************

IntervalD operator * (double x, const IntervalD& y)
  {
  IntervalD *tmp=0;
  int idx;
  bool success = (const_cast<IntervalD*>(&y))->getIntervalDPool()->obtainItem(y, tmp, idx);
  if (success) 
  {
    tmp->reappend(y.aPool);
    tmp->VF = x * y.VF;

    long n = (const_cast<IntervalD*>(&y))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = x * y.DvDx[i];
  }
  return *tmp;
  }

// *********************************************

IntervalD operator * (const IntervalD& x, double y)
  {
  IntervalD *tmp=0;
  int idx;
  bool success = (const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx);
  if (success) 
  {
    tmp->reappend(x.aPool);
    tmp->VF = x.VF * y;

    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = y * x.DvDx[i];
  }
  return *tmp;
  }

// *********************************************

IntervalD operator / (const IntervalD& x, const IntervalD& y)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);
    tmp->VF = x.VF / y.VF;

    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = (x.DvDx[i] - tmp->VF * y.DvDx[i] / y.VF) / y.VF;
  }
  return *tmp;
  }

// *********************************************

IntervalD operator / (double x, const IntervalD& y)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&y))->getIntervalDPool()->obtainItem(y, tmp, idx))
  {
    tmp->reappend(y.aPool);
    tmp->VF = x / y.VF;

    long n = (const_cast<IntervalD*>(&y))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = -x / (y.VF*y.VF) *y.DvDx[i];
  }
  return *tmp;
  }
// *********************************************

IntervalD operator / (const IntervalD& x, double y)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);
    tmp->VF = x.VF / y;

    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = x.DvDx[i] / y;
  }
  return *tmp;
  }

void IntervalD::operator += (const IntervalD &x)
{
  VF+=x.VF;
  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  for(int i=0;i<n;i++)
    DvDx[i]+=x.DvDx[i]; 
}

void IntervalD::operator -= (const IntervalD &x)
{
  VF-=x.VF;
  long n = (const_cast<IntervalD*>(&x))->getDimension(); 
  for(int i=0;i<n;i++)
    DvDx[i]-=x.DvDx[i]; 
}

void IntervalD::operator *= (const IntervalD &y)
{
  VF=VF*y.VF;  
  long n = (const_cast<IntervalD*>(&y))->getDimension(); 
  for(int i=0;i<n;i++)
    DvDx[i]=VF*y.DvDx[i]+DvDx[i]*y.VF;
}

//------------------------------------------
IntervalD ilog(const IntervalD &x)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);
    tmp->VF = ilog(x.VF);
    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);
    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = x.DvDx[i] / x.VF;
  }
  return *tmp;
  }

//------------------------------------------
IntervalD ilog10(const IntervalD &x)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);
    tmp->VF = ilog10(x.VF);
    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);

    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = x.DvDx[i] / (x.VF*ln10);
  }
  return *tmp;
  }


IntervalD iexp(const IntervalD &x)
  {
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);

    tmp->VF = iexp(x.VF);

    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);

    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = tmp->VF * x.DvDx[i];
  }
  return *tmp;
  }

IntervalD pow(const IntervalD &x,double p)
  { 
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);
    tmp->VF = iipow(x.VF, p);
    Interval tmp2 = p * iipow(x.VF, p - 1.0);

    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);

    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = tmp2 * x.DvDx[i];
  }
  return *tmp;
  }

IntervalD power(const IntervalD &x,long p)
  { 
  IntervalD *tmp=0;
  int idx;
  if ((const_cast<IntervalD*>(&x))->getIntervalDPool()->obtainItem(x, tmp, idx))
  {
    tmp->reappend(x.aPool);

    tmp->VF = power(x.VF, p);
    Interval tmp2 = p * power(x.VF, p - 1);

    long n = (const_cast<IntervalD*>(&x))->getDimension();
    //tmp->resize(n);

    for (int i = 0; i < n; i++)
      tmp->DvDx[i] = tmp2 * x.DvDx[i];
  }
  return *tmp;
  }
#endif