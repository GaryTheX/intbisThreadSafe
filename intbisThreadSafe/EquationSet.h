#pragma once

#include "intBox.h"

class EquationSet
  {  
  public:
    EquationSet(void);
    virtual ~EquationSet(void);

    virtual void setBenchMarkID(int bmID){;};
    virtual void setParameters(long np, double *p){;}
    virtual void setInitialValue(long nx, double *x0){;}

    virtual void fun(intBox *pBox, long idx=-1){;};
    virtual void jac(intBox *pBox, long idxEq=-1, long idxVr=-1){;};
    virtual void tens(intBox *pBox, long idxEq, long idxVr){;};
    virtual bool isTensProvided()
      {
      return false;
      };

    virtual bool obt_buildInBound(intBox *pBox){return false;};

    virtual char* getName(){return "";};

  };

