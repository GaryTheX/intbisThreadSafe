#include<stdio.h>

#include "intbisThread.h"


intbisThread::intbisThread(void)
  {
  intbisStatus = UNKNOWN;
  intbisPool = 0;
  ID = 0;  iterCnt = 0;
  doMultiThreading = false;
  }

intbisThread::~intbisThread(void)
  {
  if(ID==0)//intbisPool) 
    {
    delete intbisPool;
    intbisPool = 0;
    }

  if(doMultiThreading && ID!=0)
    {
    WaitForSingleObject( hRunMutex, INFINITE );
    CloseHandle( hRunMutex );
    }
  }

void intbisThread::multiThreading(int numOfThr)
  {
  if(numOfThr>=2)
    {
    setMutex();
    intbisPool = new aThreadPool();
    intbisPool->initialize(numOfThr);
    }
  }

void intbisThread::setStatus(ThreadStatus askedStatus)
  {
  bool setInitial = false;
  while(!setInitial)
    {
    DWORD mutexStatus = WaitForSingleObject(hRunMutex, 75L);
    switch(mutexStatus)
      {
      case WAIT_OBJECT_0:
        {    
        intbisStatus = askedStatus;
        setInitial = true;
        ReleaseMutex(hRunMutex); 
        break;
        }
      default:
        break;
      }
    }
  }

void intbisThread::appendItems(intBox *pBox, EquationSet *eqSet)
  {
  if(!doMultiThreading)
    {
    pBoxExtr=pBox; 
    eqSetExtr=eqSet;
    initialize(pBoxExtr->szX);
    return;
    }
  bool setInitial = false;
  while(!setInitial)
    {
    DWORD mutexStatus = WaitForSingleObject(hRunMutex, 75L);
    switch(mutexStatus)
      {
      case WAIT_OBJECT_0:
        {    
        pBoxExtr=pBox; 
        eqSetExtr=eqSet;
        initialize(pBoxExtr->szX);
        ReleaseMutex(hRunMutex);  
        setInitial = true;
        break;
        }
      default:
        break;
      }
    }
  }

void intbisThread::appendItems(nativeStack *astack, EquationSet *eqSet, bool shareHalf)
  {
  if(!doMultiThreading)
    {
    pBoxExtr = astack->pop();
    eqSetExtr= eqSet;
    initialize(pBoxExtr->szX);
    while(astack->is_empty()==false)
      {
      stack_wrk->push(astack->pop());
      }
    return;
    }
  bool setInitial = false;
  while(!setInitial)
    {
    DWORD mutexStatus = WaitForSingleObject(hRunMutex, 75L);
    switch(mutexStatus)
      {
      case WAIT_OBJECT_0:
        {   
        if(shareHalf)
          {
          pBoxExtr = astack->pop();
          eqSetExtr= eqSet;
          initialize(pBoxExtr->szX);

          long itemCnt = astack->numberOfItems();
          long itemHalf= (itemCnt+1)/2, i=0;
          while(i<itemHalf)
            {
            stack_wrk->push(astack->pop());
            ++i;
            }          
          }
        else
          {
          pBoxExtr = astack->pop();
          eqSetExtr= eqSet;
          initialize(pBoxExtr->szX);
          while(astack->is_empty()==false)
            {
            stack_wrk->push(astack->pop());
            }
          }
        ReleaseMutex(hRunMutex);
        setInitial = true;
        break;
        }
      default:
        break;
      }
    }
  }

// *****************************

void intbisThread::ftest(intBox *pBox, long& UNKNOWN, long &SIGRT, EquationSet *eqSet)
  {
  long retcon;
  
  retcon=hsing(pBox,eqSet);
  //printf("\n *** ID == %2d is executing with status %2d",ID, getStatus());

  if(retcon>=iLUNKNW && retcon<=iUUNKNW) // roots inclusion
    {
    UNKNOWN = iTRUE;
    
    long ip = pvslct(pBox);
    if(ip<0)
      {
      UNKNOWN = iFALSE; 
      SIGRT   = iTRUE;
      }
    else
      {      
      double  diamReal = diamcp(pBox->szX, pBox->X);        
      if(diamReal < EPSBOX)
        {
        UNKNOWN = iFALSE;
        SIGRT   = iFALSE;
        }        
      }
    }
  else if(retcon>=iLFALSE && retcon<=iUFALSE) // no root included
    {
    UNKNOWN = iFALSE;
    SIGRT   = iFALSE;
    }
  else if(retcon>=iLTRUE && retcon<=iUTRUE) // unique root
    {
    UNKNOWN = iFALSE;
    SIGRT   = iTRUE;
    }
  else
    {
    UNKNOWN = iTRUE;
    //printf("\n retcon == %4d",retcon);
    }
  }

// *****************************

bool intbisThread::solveOneIterationWithConclusion()
  {
  bool concluded = false;
  //long N = pBoxExtr->szX;
  long nleaves = 0, indentCount = 1;
  bool branchDone = 0, branchSolve = 0, initSplit = 0;
  double r = 0.5; // r is a parameter used in the expansion/deletion process
  long UNKNOWN, SIGRT;

  double error = EPSBOX * 2.5;

  long testFlg; // = 0, no more boxes in the stack; 1, boxes remain in the stack.
  long testNum=0; // total number of ftest call

  ftest(pBoxExtr,UNKNOWN,SIGRT,eqSetExtr);
  testNum++;

  if(!UNKNOWN)
    {
    concluded = true;
    testFlg = 0;
    if(SIGRT)
      {
      long totalSol = stack_sln->numberOfItems();
      bool dupNode = false;
      if(totalSol>0) 
        {
        expand(pBoxExtr,r);
        stack_sln->dellst(pBoxExtr,error,r, dupNode);
        }
      if(!dupNode) 
        stack_sln->push(pBoxExtr);
      }
    else
      stack_rst->push(pBoxExtr);
    }
    
  testFlg = 0;

  return concluded;
}

// ************************

void intbisThread::solve()
  {
  iterCnt = 0;
  if(doMultiThreading)
    {
    bool concluded = solveOneIterationWithConclusion();
    if(!concluded)
      {
      stack_wrk->push(pBoxExtr);
      }

    if(stack_wrk->is_empty()==false)
      {
      if(!concluded)
        {
        pBoxExtr = stack_wrk->pop();
        Box2 = generateBox(pBoxExtr);
        bisect(pBoxExtr, Box2);
        stack_wrk->push(Box2);
        stack_wrk->push(pBoxExtr);
        }
      intbisPool->createObj(stack_wrk, eqSetExtr, true);
      }
    else
      {
      long retcon = newtonRaphson(pBoxExtr, eqSetExtr);
      }

    long i;
    bool allDone = false;

    while(!allDone)
      {
      long waitingCount = 0, threadCount = 0;

      for(i=0; i<intbisPool->poolSize(); ++i)
        {
        if(intbisPool->athread[i]!=0) // these athreads are created be the call above "intbisPool->createObj"
          {
          threadCount++;
          if(intbisPool->athread[i]->getStatus()==WAITING)
            waitingCount++;
          else
            break;
          }
        }

      if(threadCount==waitingCount)
        allDone = true;

      if(allDone)
        {
        for(i=0; i<intbisPool->poolSize(); ++i)
          {
          if(intbisPool->athread[i]!=0)
            {
            intbisPool->athread[i]->setStatus(TERMINATE);
            intbisPool->athread[i]->getherSolution(stack_sln);
            }
          }
        WaitForSingleObject( hRunMutex, INFINITE );      
        CloseHandle( hRunMutex );

        break;
        }
      }
    }
  else
    {
    while(true)
      {
      while(!solveOneIterationWithConclusion())
        {
        bisect(pBoxExtr);
        }
      iterCnt++;
      if(stack_wrk->is_empty())
        {
        long retcon = newtonRaphson(pBoxExtr, eqSetExtr);
        break;
        }
      pBoxExtr = stack_wrk->pop();
      }
    }
  }

// *****************************

void intbisThread::solveByThis()
  {
  iterCnt = 0;
  while(true)
    {    
    switch(getStatus())
      {
      case EXECUTING:
        {
        bool isConcluded = solveOneIterationWithConclusion();  iterCnt++;
        if(!isConcluded)
          {
          Box2 = generateBox(pBoxExtr);
          bisect(pBoxExtr, Box2);
          stack_wrk->push( Box2);
          intbisPool->createObj(stack_wrk, eqSetExtr, false);
          }
        else
          {
          if(stack_wrk->is_empty())
            {
            //long retcon = newtonRaphson(pBoxExtr, eqSetExtr);
            setStatus(WAITING);
            }
          else
            {
            intbisPool->createObj(stack_wrk, eqSetExtr, false);
            if(stack_wrk->is_empty())
              setStatus(WAITING);
            else
              pBoxExtr = stack_wrk->pop();
            }
          }

        break;
        }
      case WAITING:
        break;
      case TERMINATE:
        break;
      default:
        break;
      }
    if(getStatus()==TERMINATE)
      break;
    }
  }
