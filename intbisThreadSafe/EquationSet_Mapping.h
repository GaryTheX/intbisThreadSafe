#pragma once
#include "equationset.h"
class EquationSet_Mapping :  public EquationSet
  {
  protected:
    long nparam, icase, neq;
    double *param;

    void flushMemory();

    void jacDigit(intBox *pBox, long idxEq, long idxVr);

  public:
    EquationSet_Mapping(void);
    virtual ~EquationSet_Mapping(void);

    virtual void setParameters(long np, double *p);
    virtual void setBenchMarkID(int bmID);

    virtual void fun(intBox *pBox, long idx=-1);
    virtual void jac(intBox *pBox, long idxEq=-1, long idxVr=-1);
    virtual void tens(intBox *pBox, long idxEq, long idxVr);
    //virtual bool isTensProvided(){return true;};

    virtual bool obt_buildInBound(intBox *pBox);

    char* getName();
  };

