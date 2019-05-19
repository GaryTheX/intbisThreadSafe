#pragma once
#include "Interval.h"
#include "intBox.h"
#include "EquationSet.h"
#include "EquationSet_Mapping.h"
#include "EquationSet_LogisticMappingIV.h"
#include "EquationSet_Robot_kinematics.h"
#include "EquationSet_Neurophysiology.h"
#include "EquationSet_iBench.h"
#include "EquationSet_ChemicalEquilib.h"
#include "EquationSet_Broyden.h"
#include "EquationSet_Economic.h"
#include "intbis.h"
#include "intbisThread.h"

class intSolEngine
  {
  private:
    long neq, nThreads;
    bool doMultiThread;
    EquationSet *eqSetI;
    intBox *pX;
    intbisThread *pSolThread;
    intbis *pSol;

    void nullPointers();
    void killPointers();

  public:
    intSolEngine(void);
    virtual ~intSolEngine(void);

    void reinitialize();

    void setMultiThreadingConfig(int nThread=1);
    void setEquationSet(int eqId, int caseId=-1);
    bool isReady();

    void solve();
    long solutionsFound(); 
    void getSolution(double *pl, double *pu);
    long getDimension();

    char* getName(){return eqSetI->getName();};
    int getNumThreads(){return nThreads;};
  };

