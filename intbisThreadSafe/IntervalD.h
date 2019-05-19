#pragma once

class intervalDPool;

#include "interval.h"
#include "intervalDPool.h"
enum DType{INDEPENDENT, DEPENDENT, INTERMIDIATE};

class IntervalD :  public Interval
  {
  private:
    DType adType;
    bool appended;
    long Ndim, idx_Var, idx_Der;
    Interval *DvDx, VF;

    intervalDPool *aPool;

  public:
    IntervalD();
    IntervalD(int n);
    IntervalD(double left, double right);
    virtual ~IntervalD(void);

    void killPointers();
    void setDependency(DType dType){adType = dType;}

    void resize(long N);
    void reappend(intervalDPool *outsidePool);
    intervalDPool* getIntervalDPool(){return aPool;};

    void setIdx_Var(long idx);
    void setValue(Interval x);

    IntervalD &operator = (const IntervalD&);
    virtual IntervalD &operator = (double);
    virtual IntervalD &operator = (long);

    friend IntervalD operator + (const IntervalD&,const IntervalD&);
    friend IntervalD operator + (double,const IntervalD&);
    friend IntervalD operator + (const IntervalD&,double);

    friend IntervalD operator - (const IntervalD&,const IntervalD&);
    friend IntervalD operator - (double,const IntervalD&);    
    friend IntervalD operator - (const IntervalD&,double);
    
    friend IntervalD operator * (const IntervalD&,const IntervalD&);
    friend IntervalD operator * (double,const IntervalD&);    
    friend IntervalD operator * (const IntervalD&,double);

    friend IntervalD operator / (const IntervalD&,const IntervalD&);
    friend IntervalD operator / (double,const IntervalD&);    
    friend IntervalD operator / (const IntervalD&,double);

    void  operator +=(const IntervalD&);
    void  operator -=(const IntervalD&);
    void  operator *=(const IntervalD&);

    friend IntervalD ilog(const IntervalD&);
    friend IntervalD ilog10(const IntervalD&);
    friend IntervalD iexp(const IntervalD&);
    friend IntervalD pow(const IntervalD&,double);
    friend IntervalD power(const IntervalD& x, long n);

    long getDimension(){return Ndim;};
    Interval getValue(){return VF;};
    Interval getDerivative(int i){return DvDx[i];};
    Interval* getDerivative(){return DvDx;};
  };

