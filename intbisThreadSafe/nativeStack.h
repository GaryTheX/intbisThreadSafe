#pragma once
#include "intBox.h"

class nativeNode
  {
  public:
    intBox *abox;
    nativeNode *prev;
  };

class nativeStack
  {
  private:
    long listlen;
    nativeNode *top;

  public:
    nativeStack(void);
    virtual ~nativeStack(void);

    void push(intBox *aBox);
    intBox* pop();

    bool is_empty();
    void dellst(intBox* aItem, double error,double r,bool& dupNode);
    void remove(bool killADPool = false);
    bool itemCopy(long ith, intBox &aItem);
    long numberOfItems(){return listlen;}
  };

