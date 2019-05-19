#include "myStack.h"

myStack::myStack(void)
  {
  mylist=0;
  top = mylist;
  listlen = 0;
  }

myStack::~myStack(void)
  {
  remove();
  }

//test if the link list is empty
bool myStack::is_empty()
{
  //return mylist==0 ? true:false;
  if(listlen==0) return true;
  else return false;
}

//intert a item into the link list

void myStack::push(intBox *aItem)
{
  if(aItem!=0)
    {
    mylist = top;
    anItem *pt = new anItem();

    pt->val = aItem;

    top = pt;
    if(listlen==0)
      {
      top->prev=0;
      }
    else
      {
      top->prev=mylist;
      }
    listlen++;
    }
}

intBox* myStack::pop()
{
  intBox *pt;
  if(listlen!=0)
    {
    listlen--;    
    if(!top) pt = 0;
    else
      {
      pt = top->val;
      top = top->prev;
      }
    }
  else
    {
    pt = 0;
    }
  return pt;
}

void myStack::pop(intBox **item)
{
  *item = pop();
}

//copy a item from the link list
bool myStack::itemCopy(long ith, intBox &aItem)
{
  bool succeed=0;
  if(listlen!=0 && ith<=listlen)
    {
    long k=0;
    mylist = top;
    for(anItem *pt=mylist;pt;pt=pt->prev)
      {
      k++;
      if(k==ith)
        {
        aItem=(*(pt->val));
        succeed=1;
        break;
        }
      }
    } 

  return succeed;
}

//remove the whole link list
void myStack::remove()
{
  //Item *pt=mylist;
  anItem *pt = top;
  while(pt){
    anItem *tmp=pt;
    //pt=pt->next;
    pt = pt->prev;
    delete tmp;
  }
  listlen=0; top = 0;
}


//delete the redounracy item from the link list
void myStack:: dellst(intBox* aItem,double error,double r, bool &dupNode)
{
  anItem *pt = top;
  long intsecFlg;
  
  double err = error;//*r/2.0;
  if(pt)
    {//check the whole
    intsecFlg=1;
    if(!((*aItem) & (*(pt->val))))
      intsecFlg = 0;
    if(intsecFlg)
      {
      if(aItem == pt->val)
        {
        dupNode = true;
        return;
        }
      anItem *tmp=pt->prev;
      delete pt->val; 
      //delete pt;      
      pt=tmp;
      listlen--; 
      if(listlen>0) top = pt;
      }
    }
  
  if((mylist=pt)==0){//it at the end, exit
    at_end=0;
    return;
  }

  anItem *prv=pt;pt=pt->prev;
 
  while(pt){ //check other items
    intsecFlg=1;
    if(!((*aItem) & (*(pt->val))) && aItem!=pt->val)
      intsecFlg = 0;
    if (intsecFlg)
      {
      if(aItem == pt->val)
        {
        dupNode = true;
        return;
        }
      prv->prev=pt->prev;
      if(at_end==pt)
	      at_end=prv;
      delete pt->val;
      //delete pt;
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
