#pragma once

#include "intbis.h"
#include "aThreadPool.h"
enum ThreadStatus{UNKNOWN=0, STARTING, EXECUTING, WAITING, TERMINATE};

class aThreadPool;

class intbisThread :  public intbis
  {
  private:
    ThreadStatus intbisStatus;
    aThreadPool *intbisPool;
    bool doMultiThreading;
    long ID, iterCnt;
    HANDLE  hRunMutex;

  public:
    static intbisThread* current;

    aThreadPool* getPool(){return intbisPool;};
    
  public:
    intbisThread(void);   
    virtual ~intbisThread(void);

    void assignPoolAddress(aThreadPool *intbisPooladd){intbisPool = intbisPooladd;};

    intbisThread* createAnObj(aThreadPool *intbisPool, bool full = false);
    void setID(long id){ID = id;};
    long getID(){return ID;};
    //virtual void appendItems(intBox *pBox, EquationSet *eqSet){pBoxExtr = pBox; eqSetExtr = eqSet;};
    ThreadStatus getStatus(){return intbisStatus;};
    void setStatus(ThreadStatus askedStatus);

    void multiThreading(int numOfThr);

    virtual void appendItems(intBox *pBox, EquationSet *eqSet);

    virtual void appendItems(nativeStack *astack, EquationSet *eqSet, bool shareHalf=false);

    void ftest(intBox *pBox, long& UNKNOWN, long& SIGRT, EquationSet *eqSet);    
    bool solveOneIterationWithConclusion();
    void solve();
    void solveByThis();

    void setMutex()
      {
      doMultiThreading = true;
      hRunMutex = CreateMutex( NULL, FALSE, NULL );
      };
  };



