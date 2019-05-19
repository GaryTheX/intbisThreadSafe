#pragma once
#include "equationset.h"

class EquationSet_iBench :  public EquationSet
  {
  private:
    int icase;
    long neq;
    double *pa, *pb;
    void jacDigit(intBox *pBox, long idxEq, long idxVr);

  public:
    EquationSet_iBench(void);
    virtual ~EquationSet_iBench(void);
    
    virtual void setBenchMarkID(int bmID);
    virtual void setParameters(long np, double *p);
    virtual void fun(intBox *pBox, long idx=-1);
    virtual void jac(intBox *pBox, long idxEq=-1, long idxVr=-1);

    virtual bool obt_buildInBound(intBox *pBox);

    char* getName();
  };

