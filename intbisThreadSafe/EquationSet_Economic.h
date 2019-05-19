#pragma once
#include "equationset.h"
class EquationSet_Economic :  public EquationSet
  {
  private:
    int icase;
    long neq;
    double *pc;

    void jacDigit(intBox *pBox, long idxEq, long idxVr);

  public:
    EquationSet_Economic(void);
    virtual ~EquationSet_Economic(void);

    virtual void setBenchMarkID(int bmID);

    virtual void fun(intBox *pBox, long idx=-1);
    virtual void jac(intBox *pBox, long idxEq=-1, long idxVr=-1);

    virtual bool obt_buildInBound(intBox *pBox);
    char* getName(){return "Economic";};
  };

