#pragma once
#include "equationset.h"
class EquationSet_ChemicalEquilib :  public EquationSet
  {
  private:
    int icase;
    long neq;
    double R, R5, R6, R7, R8, R9, R10, sqrt40;

    void jacDigit(intBox *pBox, long idxEq, long idxVr);

  public:
    EquationSet_ChemicalEquilib(void);
    virtual ~EquationSet_ChemicalEquilib(void);

    virtual void setBenchMarkID(int bmID);

    virtual void fun(intBox *pBox, long idx=-1);
    virtual void jac(intBox *pBox, long idxEq=-1, long idxVr=-1);
    virtual bool obt_buildInBound(intBox *pBox);
    char* getName();
  };

