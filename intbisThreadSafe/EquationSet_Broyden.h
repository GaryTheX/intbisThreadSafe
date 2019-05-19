#pragma once
#include "equationset.h"
class EquationSet_Broyden :
  public EquationSet
  {
  private:
    int icase;
    long neq;

  public:
    EquationSet_Broyden(void);
    virtual ~EquationSet_Broyden(void);
    
    virtual void setBenchMarkID(int bmID);
    
    virtual void fun(intBox *pBox, long idx=-1);
    virtual void jac(intBox *pBox, long idxEq=-1, long idxVr=-1);

    virtual bool obt_buildInBound(intBox *pBox);
		char* getName();
  };

