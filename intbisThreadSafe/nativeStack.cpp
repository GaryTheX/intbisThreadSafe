#include "nativeStack.h"

nativeStack::nativeStack(void)
  {
  listlen = 0;
  top = 0;
  }

nativeStack::~nativeStack(void)
  {
  remove();
  }

void nativeStack::push(intBox *aBox)
  {
  nativeNode *newTop = new nativeNode;
  newTop->abox = aBox;  // pointer
  newTop->prev = top;

  top = newTop;
  listlen++;
  }

intBox* nativeStack::pop()
  {
  intBox* aBox = 0;
  if(top==0)
    return aBox;
  else
    {
    nativeNode* old = top;
    top = top->prev;
    listlen--;
    aBox = old->abox;
    delete old;
    }
  return aBox;
  }

//test if the link list is empty
bool nativeStack::is_empty()
{
  //return mylist==0 ? true:false;
  if(listlen==0) return true;
  else if(top==0) return true;
  else return false;
}

//delete the redounracy item from the link list
void nativeStack::dellst(intBox* aItem,double error,double r, bool &dupNode)
{
  nativeNode *pt = top;
  long intsecFlg;
  
  if(pt)
    {//check the whole
    intsecFlg=1;
    if(!((*aItem) & (*(pt->abox))))
      intsecFlg = 0;
    if(intsecFlg)
      {
      if(aItem == pt->abox)
        {
        dupNode = true;
        return;
        }
      nativeNode *tmp=pt->prev;
      delete pt->abox; 
      delete pt;      
      pt=tmp;
      listlen--; 
      top = pt;
      }
    }
  
  if(top==0){//it at the end, exit
    return;
  }

  nativeNode *prv=pt;pt=pt->prev;
 
  while(pt){ //check other items
    intsecFlg=1;
    if(!((*aItem) & (*(pt->abox))) && aItem!=pt->abox)
      intsecFlg = 0;
    if (intsecFlg)
      {
      if(aItem == pt->abox)
        {
        dupNode = true;
        return;
        }
      prv->prev=pt->prev;

      delete pt->abox;
      delete pt;
      pt=prv->prev;
      listlen--;
      }
    else
      {
      prv=pt;
      pt=pt->prev;
      }
    }
  }

//remove the whole link list
void nativeStack::remove(bool killADPool)
{
  //Item *pt=mylist;
  nativeNode *pt = top;
  while(pt){
    nativeNode *tmp=pt;
    //pt=pt->next;
    pt = pt->prev;
    if(killADPool)
      tmp->abox->setToKillADPool();
    delete tmp->abox;
    delete tmp;
  }
  listlen=0; top = 0;
}
//copy a item from the link list
bool nativeStack::itemCopy(long ith, intBox &aItem)
{
  bool succeed=0;
  nativeNode *mynode;
  if(listlen!=0 && ith<=listlen)
    {
    long k=0;
    mynode = top;
    for(nativeNode *pt=mynode;pt;pt=pt->prev)
      {
      k++;
      if(k==ith)
        {
        aItem=(*(pt->abox));
        succeed=1;
        break;
        }
      }
    } 

  return succeed;
}
