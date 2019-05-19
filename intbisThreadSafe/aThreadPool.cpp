#include<stdio.h>
#include "aThreadPool.h"

intbisThread* intbisThread::current = 0;
void CallBackFn( void * MyID )
{
  long i = (long) MyID;
  intbisThread::current->getPool()->athread[i]->solveByThis();
}

aThreadPool::aThreadPool(void)
  {
  athread = 0;
  poolSz = initSz = 0;
  isFull = false;
  hRunMutex = CreateMutex( NULL, FALSE, NULL );
  }

aThreadPool::~aThreadPool(void)
  {
  int i;
  if(athread!=0)
    {
    for(i=0; i<poolSz; ++i)
      {
      if(athread[i]) delete athread[i];
      }
    }
  delete athread;
      
  WaitForSingleObject( hRunMutex, INFINITE );
  CloseHandle( hRunMutex );
  }

intbisThread* aThreadPool::createObj(nativeStack *stack_wrk, EquationSet *eqSet, bool full)
  {
  intbisThread *aThread = 0;

  long i;
  DWORD mutexStatus = WaitForSingleObject(hRunMutex, 75L);
  switch(mutexStatus)
    {
    case WAIT_OBJECT_0:
      {
      bool  created = false;       
        
      if(!isPoolFull())
        {
        for(i=0; i<poolSz; ++i)
          {
          if(athread[i]==0)
            {
            athread[i] = new intbisThread();
            athread[i]->assignPoolAddress(this);
            aThread = athread[i];
            aThread->setMutex();
            aThread->setID(i+1);

            if(stack_wrk->is_empty())
              {
              aThread->setStatus(WAITING);
              }
            else
              {
              aThread->setStatus(EXECUTING);
              if(full)
                {
                aThread->appendItems(stack_wrk, eqSet);
                }
              else
                {
                aThread->appendItems(stack_wrk->pop(), eqSet);
                }
              }
            intbisThread::current = aThread;
            _beginthread(CallBackFn, 0, (void*) i);           
            initSz++;
            if(initSz>=poolSz) 
              isFull = true;
            created = true;
            break;
            }
          }
        }
      if(!created)
        {
        for(i=0; i<poolSz; ++i)
          {
          if(athread[i]->getStatus()==WAITING)
            {
            aThread = athread[i];
            aThread->appendItems(stack_wrk, eqSet, true);
            aThread->setStatus(EXECUTING);
            break;
            }
          }
        }

      break;
      }
    case WAIT_ABANDONED:
    case WAIT_TIMEOUT:
    case WAIT_FAILED:
    default:
      break;
    }
  ReleaseMutex(hRunMutex);
  return aThread;
  };
