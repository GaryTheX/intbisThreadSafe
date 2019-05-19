#pragma once
#include "equationset.h"
class EquationSet_Neurophysiology :  public EquationSet
  {
  private:
    int icase;
    long neq;
    double *pc;

    void jacDigit(intBox *pBox, long idxEq, long idxVr);

  public:
    EquationSet_Neurophysiology(void);
    virtual ~EquationSet_Neurophysiology(void);

    virtual void fun(intBox *pBox, long idx=-1);
    virtual void jac(intBox *pBox, long idxEq=-1, long idxVr=-1);
    virtual bool obt_buildInBound(intBox *pBox);

    virtual void tens(intBox *pBox, long idxEq, long idxVr);
    //virtual bool isTensProvided(){return true;};
    char* getName();
  };

