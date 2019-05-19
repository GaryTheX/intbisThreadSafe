#pragma once

class IntervalD;

#include "IntervalD.h"

class intervalDPool
  {
  private:
    static intervalDPool *Me;
    static int refs;

    IntervalD **anItem;
    bool *isFree;

    int MaxSize, currSize, currFreeIdx;

  protected:
    intervalDPool();
    virtual ~intervalDPool(void);

  public:
    static intervalDPool* getInstance();
    static bool releaseInstance();
    static void deleteInstance();
    
    bool obtainItem(const IntervalD& x, IntervalD* &obt, int &idx);
    void setFree();
  };

