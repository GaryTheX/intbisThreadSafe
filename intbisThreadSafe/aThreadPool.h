#pragma once
#include <windows.h>
#include <process.h>
#include "intbisThread.h"

class intbisThread;

class aThreadPool
  {
  private:
    int poolSz, initSz, activeID;    
    bool isFull; 
    HANDLE  hRunMutex;
  public:
    intbisThread **athread;   
    
  public:
    aThreadPool();
    virtual ~aThreadPool();

    bool isPoolFull(){return isFull;};
    void initialize(int n)
      {
      poolSz = n;
      athread = new intbisThread*[n];
      for(int i = 0; i<n; ++i) athread[i] = 0;
      }
    int poolSize(){return poolSz;};

    intbisThread* createObj(nativeStack *stack_wrk, EquationSet *eqSet, bool full=false);
  };

