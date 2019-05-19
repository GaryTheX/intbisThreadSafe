#pragma once
#include "intBox.h"

class anItem{
  friend class myStack;
private:

  intBox *val;
protected:
  //anItem *next;
  anItem *prev;
  
public:
  anItem(){;};
  ~anItem(){;};
};

class myStack
  {
  public:
    myStack(void);
    ~myStack(void);
    void remove();
    void dellst(intBox* aItem, double error,double r,bool& dupNode);
    virtual void push(intBox *aItem);
    virtual intBox* pop();
    virtual void pop(intBox **aItem);
    //virtual intBox* queueOut();
    bool itemCopy(long ith, intBox &aItem);
    bool is_empty();
    long numberOfItems(){return listlen;}
    anItem *at_end;

private:
  anItem *mylist;
  anItem ** PoolItem;
  long PoolLength, PoolTop;
protected:
  anItem *top;
  long listlen;
  };
