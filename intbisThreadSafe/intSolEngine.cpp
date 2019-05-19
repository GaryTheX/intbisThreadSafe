#include <iostream>
#include "intSolEngine.h"

intSolEngine::intSolEngine(void)
  {
  doMultiThread = false;
  nullPointers();
  }


intSolEngine::~intSolEngine(void)
  {
  killPointers();
  }

void intSolEngine::nullPointers()
  {
  nThreads = 0;
  neq = 0;
  eqSetI = 0;
  pX = 0;
  pSolThread = 0;
  pSol = 0;
  }

void intSolEngine::killPointers()
  {
  delete eqSetI;
  delete pSolThread;
  delete pSol;
  // pX deleted as part of pSol or pSolThread
  nullPointers();
  }

void intSolEngine::reinitialize()
  {
  killPointers();
  }

void intSolEngine::setMultiThreadingConfig(int nThread)
  {
  nThreads = nThread;
  if(nThread>1)
    {
    doMultiThread = true;
    pSolThread = new intbisThread();
    pSolThread->multiThreading(nThread);      
    }
  else
    {
    doMultiThread = false;
    pSol = new intbis();
    }
  }

void intSolEngine::setEquationSet(int eqId, int caseId)
  {
  switch(eqId)
    {
    case 1:
      {
      eqSetI = new EquationSet_iBench();  
      eqSetI->setBenchMarkID(caseId);
      break;
      }
    case 2:
      {
      eqSetI = new EquationSet_Economic(); 
      eqSetI->setBenchMarkID(caseId);  // ID=4 will exhaust memory
      break;
      }
    case 3:
      {
      eqSetI = new EquationSet_Mapping();  
      eqSetI->setBenchMarkID(caseId);
      break;
      }
    case 4:
      {
      eqSetI = new EquationSet_ChemicalEquilib();
      eqSetI->setBenchMarkID(caseId);
      break;
      }
    case 5:
      {
      eqSetI = new EquationSet_Robot_kinematics();
      break;
      }
		case 6:
			{
			eqSetI = new EquationSet_Neurophysiology();
			break;
			}
		case 7:
			{
			eqSetI = new EquationSet_Broyden(); 
			eqSetI->setBenchMarkID(caseId);
			break;
			}
    default:
      break;
    }
  // not so well
  //eqSetI = new EquationSet_Neurophysiology();
  //eqSetI = new EquationSet_Broyden(); eqSetI->setBenchMarkID(3);
  }

long intSolEngine::getDimension()
  {
  return pX->szX;
  }

bool intSolEngine::isReady()
  {
  bool ready = false;
  if(eqSetI && ((doMultiThread && pSolThread) || (!doMultiThread && pSol)))
    {
    ready = true;
    pX = new intBox();
    eqSetI->obt_buildInBound(pX);
    if(doMultiThread) pSolThread->appendItems(pX, eqSetI);
    else pSol->appendItems(pX, eqSetI);
    }
  return ready;
  }

void intSolEngine::solve()
  {
  try {
    if (doMultiThread) pSolThread->solve();
    else pSol->solve();
  }
  catch (...) {
    printf("Exception occurred");
  }
  }

long intSolEngine::solutionsFound()
  {
  if(doMultiThread) return pSolThread->SolutionFound();
  else return  pSol->SolutionFound();
  }

void intSolEngine::getSolution(double *pl, double *pu)
  {
  if(doMultiThread) pSolThread->getSolution(pl, pu);
  else pSol->getSolution(pl,pu);
  }