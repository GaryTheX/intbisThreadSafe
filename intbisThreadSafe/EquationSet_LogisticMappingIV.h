#pragma once
#include "equationset_mapping.h"
class EquationSet_LogisticMappingIV :  public EquationSet_Mapping
  {
  private:
    long nx0;
    double *pX0;
  public:
    EquationSet_LogisticMappingIV(void);
    virtual ~EquationSet_LogisticMappingIV(void);

    void setInitialValue(long nx, double *x0);

    void fun(intBox *pBox, long idx=-1);
    void jac(intBox *pBox, long idxEq=-1, long idxVr=-1);
    virtual void tens(intBox *pBox, long idxEq, long idxVr);
    //virtual bool isTensProvided(){return true;};
  };

