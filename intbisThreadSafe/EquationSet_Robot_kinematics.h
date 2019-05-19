#pragma once
#include "EquationSet.h"

class EquationSet_Robot_kinematics : public EquationSet
  {
  private:
    int icase;

    long neq;
    double *pa, *pBL, *pBU;
    void jacDigit(intBox *pBox, long idxEq, long idxVr);

  public:
    EquationSet_Robot_kinematics(void);
    virtual ~EquationSet_Robot_kinematics(void);

    virtual void fun(intBox *pBox, long idx=-1);
    virtual void jac(intBox *pBox, long idxEq=-1, long idxVr=-1);

    virtual bool obt_buildInBound(intBox *pBox);
    char* getName(){return "RobotKinematic";};
  };

