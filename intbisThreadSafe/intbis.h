#pragma once

#include "intBox.h"
#include "EquationSet.h"
#include "myStack.h"
#include "nativeStack.h"

class intbis
  {  
  protected:
    long iLUNKNW,iUUNKNW,iLTRUE,iUTRUE,iLFALSE,iUFALSE,iTRUE,iFALSE;

    long do_mono, probsize;
    nativeStack *stack_wrk, *stack_rst, *stack_sln; 

    double EPSBOX, EPSFUN, *pwork, *ymat;
    Interval *tmpm1, *wiv1, *wiv2;
    Interval *FuncL, *JacbL;
    intBox *Box1, *Box2, *pBox1, *pBox2;
    long *ipvt;
    long JACFLG, PIVFLG, DERIVFLG,RCUFLG, BTreeFLG;

    void Deallocate();
    intBox* generateBox(intBox *refBox);
    long dvsbin(double xpt,Interval &x,Interval &nome,Interval &deno,
           Interval &imag1,Interval &imag2); 

    void expand(intBox *varBox, double r);

    virtual long componentWiseGS(intBox *pBox, EquationSet *eqSet, Interval *func=0, Interval *jacb=0);
    virtual long componentWiseFi(intBox *pBox, EquationSet *eqSet, Interval *jaci, long i, long k);
    virtual long parametric_poly(intBox *pBox, EquationSet *eqSet, long idxF, long idxx,
                                 Interval *Fi, Interval xi, Interval f0, Interval Df0, Interval D2f0, Interval D3f0);
    long preconditionC(intBox *pBox, Interval *matr=0);
    //virtual long LuDcmp(double *a, long n, long *idx, double *vv);
    long singleStepGS(long idmen, intBox *pBox, intBox *pBoxThin, Interval &imag1, Interval &imag2);
    long newtonRaphson(intBox *pBox, EquationSet *eqSet, double tol=1.0e-9);

  public:
    intBox *pBoxExtr;
    EquationSet *eqSetExtr;   

  public:
    intbis(void);
    virtual ~intbis(void);

    virtual void initialize(long N);
    void allocate(long N);
    long pvslct(intBox *pBox);
    void bisect(intBox *varBox, long ip=-1);
    void bisect(intBox *varBox, intBox *varBox2);
    long hsing(intBox *pBox, EquationSet *eqSet);
    virtual void ftest(intBox *pBox, long& UNKNOWN, long& SIGRT, EquationSet *eqSet);
    virtual void appendItems(intBox *pBox, EquationSet *eqSet)
      {
      pBoxExtr=pBox; 
      eqSetExtr=eqSet;
      initialize(pBoxExtr->szX);
      };
    virtual void appendItems(nativeStack *astack, EquationSet *eqSet)
      {
      pBoxExtr = astack->pop();
      eqSetExtr= eqSet;
      initialize(pBoxExtr->szX);
      while(astack->is_empty()==false)
        {
        stack_wrk->push(astack->pop());
        }
      }
    virtual void solve();
    virtual void multiThreading(int numOfThr){;};

    long SolutionFound(){return stack_sln->numberOfItems();};
    bool getSolution(double *lb, double *ub)
      {
      bool gotIt = false;
      if(!stack_sln->is_empty())
        {
        intBox *aBox = stack_sln->pop();
        for(long i=0; i<aBox->szX; ++i)
          {
          lb[i] = aBox->X[i].getleft();
          ub[i] = aBox->X[i].getright();
          }
        stack_rst->push(aBox);
        gotIt = true;
        }
      return gotIt;
      }

    void getherSolution(nativeStack *to)
      {
      while(!stack_sln->is_empty())
        {
        to->push(stack_sln->pop());
        }
      }
  };

