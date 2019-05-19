#pragma once

#include "Interval.h"
#include "IntervalD.h"

struct Taylor
  {
  Interval Taylor_fatParam, Remain_fatParam;
  bool isRemainConst;
  };

class intBox
  {
  private:
    void nullpointer();
    void flushMemory();

  public:
    intBox(void);
    intBox(const intBox& copyfrom);
    virtual ~intBox(void);

    long szX;
    Interval *X, *f, *df, fobj;
    double   *Xpt, *fpt, *dfpt, fobjpt;

    IntervalD *vx, *vf;

    Taylor *ts;
    virtual intBox &operator = (const intBox& copyfrom);
    void initialize(long n);
    virtual long operator & (const intBox& copyfrom) const;        
    void setValue(long n, Interval *x, double *xpt=0);
    void getValue(long n, Interval *x, double *xpt);

    void setValue(long i, double Lb, double Ub);

    Interval* getIntervalAddress(){return X;};

    void setToKillADPool(){
      if(vx) vx[0].killPointers();
      }
  };

