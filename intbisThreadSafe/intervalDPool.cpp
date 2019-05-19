#include "intervalDPool.h"

intervalDPool* intervalDPool::Me=0;
int intervalDPool::refs=0;

intervalDPool* intervalDPool::getInstance()
  {
  if(Me==0) Me = new intervalDPool();
  refs++;
  return Me;
  }

bool intervalDPool::releaseInstance()
  {
  bool deleted = false;
  if(refs>0)
    {
    if(--refs==0)
      {
      delete Me; Me = 0; deleted = true;
      }
    }
  return deleted;
  }
void intervalDPool::deleteInstance()
  {
  delete Me; Me = 0; refs = 0;
  }

intervalDPool::intervalDPool()
  {
  MaxSize = 1000;
  currSize= currFreeIdx = 0;
  anItem = new IntervalD*[MaxSize];
  isFree = new bool[MaxSize];
  for(int i=0; i<MaxSize; ++i)
    {
    isFree[i] = true;
    anItem[i] = 0;
    }
  }

intervalDPool::~intervalDPool(void)
  {
//  /*
  for(int i=0; i<currSize; ++i)
    {
    delete anItem[i]; anItem[i] = 0;
    }
//    */
  delete[] anItem; anItem = 0;
  delete[] isFree; isFree = 0;
  }

bool intervalDPool::obtainItem(const IntervalD& x, IntervalD* &obt, int &idx)
  {
  bool successful = false;
  if (currSize >= MaxSize)
  {
    IntervalD tmp;
    tmp.resize((const_cast<IntervalD*>(&x))->getDimension());
    obt = &tmp;
    return true;
  }

  while(!successful && currSize<MaxSize)
    {
    if(isFree[currFreeIdx])
      {
      if(!anItem[currFreeIdx]) 
        {
        anItem[currFreeIdx] = new IntervalD((const_cast<IntervalD*>(&x))->getDimension());
        anItem[currFreeIdx]->setDependency(INTERMIDIATE);
        currSize++;
        }

      obt = anItem[currFreeIdx];
      isFree[currFreeIdx] = false;
      successful = true;
      break;
      }
    else
      currFreeIdx++;
    }
  if (successful) idx = currFreeIdx;

  return successful;
  }

void intervalDPool::setFree()
  {
  for(int i=currFreeIdx; i>=0; --i) 
    isFree[i] = true;
  currFreeIdx = 0;
  }
