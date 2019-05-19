#include <math.h>
#include <iostream>
#include <ostream>

#include "intbis.h"

extern "C"
  {
  extern long ludcmp(double *a, long n, long *indx, double *d, double *vv);
	extern long LuDcmp(double *a, long n, long *indx, double *vv);
  extern void lubksb(long n,double *a,double *b,long *ipvt);
  extern void inverse(long n,double *a,double *b,double *wrk,long *ipvt);
  extern void tranverse(long n,double *a);
  extern void USEINTLIB();
  extern double gettiny();
  extern void ROUNDDOWN(double& left);
  extern void ROUNDUP(double& right);
  }

intbis::intbis(void)
  {
  probsize = 0;
  EPSBOX = 1.0e-6;
  EPSFUN = 1.0e-6;
  USEINTLIB();
  tmpm1 = 0;
  ipvt  = 0;
  pwork = 0;
  FuncL = JacbL = 0;

  do_mono = 1;

  iLUNKNW=0;
  iUUNKNW=2;
  iLTRUE=10;
  iUTRUE=12;
  iLFALSE=20;
  iUFALSE=26;
  iTRUE=1;
  iFALSE=0;

  stack_wrk = new nativeStack();
  stack_rst = new nativeStack();
  stack_sln = new nativeStack();

  //stack_Sln = new std::stack<intBox>();
  }


intbis::~intbis(void)
  {
  Deallocate();
  stack_wrk->remove(true);
  stack_rst->remove(true);
  stack_sln->remove(true);
  delete stack_wrk;
  delete stack_rst;
  delete stack_sln;
  }

void intbis::initialize(long n)
  {
  allocate(n);
  }

void intbis::allocate(long n)
  {
  if(probsize<n)
    {
    Deallocate();
    probsize = n;
    tmpm1 = new Interval[n*n + n*3];
    ipvt = new long[n];
    pwork = new double[2*n+n*n];  
    }
  
  wiv1  = &tmpm1[n*n];
  wiv2  = &wiv1[n];
  ymat = &pwork[2*n];
  }

void intbis::Deallocate()
  {
  delete[] ipvt;
  delete[] tmpm1;
  delete[] pwork;
  }

// ***************

intBox* intbis::generateBox(intBox *refBox)
  {
  intBox *aBox;
  while(true)
    {
    if(stack_rst->is_empty())
      {
      aBox = new intBox();
      break;
      }
    else
      {
      aBox = stack_rst->pop();
      if(aBox!=0 && aBox != refBox)
        break;
      }
    }
  *aBox = *refBox;
  return aBox;
  }

// *****************************

long intbis::dvsbin(double xpt, Interval &x, Interval &nomi, Interval &deno,
                    Interval &imag1, Interval &imag2)
  {
  long xcase,icase;
  Interval imag1_xpt,imag2_xpt;
  // Do the extended-value interval division.
  xdiv(xcase,nomi,deno,imag1_xpt,imag2_xpt);

  //If the result is a single interval, do the subtraction and
  // the intersection, then return.
  if(xcase!=5)
		{
    xsclsb(xcase,xpt,imag1_xpt,imag1);//If the result is a single interval, do the subtraction and
    xint(xcase,x,imag1);              // the intersection, then return.
    icase=xcase;                      // place the result in imag1
		}
  else
		{
    xcase=3;                          //for the two interval result, do the
    xsclsb(xcase,xpt,imag1_xpt,imag1);//subtraction and intersection with the first one 
    xint(xcase,x,imag1);              //place the result in imag1
    icase=xcase;                      // if the result is empty,  
    xcase=2;
    if(icase==0)
			{                     
      xsclsb(xcase,xpt,imag2_xpt,imag1);//do the subtraction and intersection with the second interval
      xint(xcase,x,imag1);              //place the result using the second interval in imag1
      icase=xcase;
		  }
    else
			{                                 //otherwise,
      xsclsb(xcase,xpt,imag2_xpt,imag2);//do the subtraction and intersection with the second interval
      xint(xcase,x,imag2);              //place the result using the second interval in imag2
      icase=xcase+icase; 
			}
		}
  return icase;
  }

void intbis::bisect(intBox *pBox, long ipi)
  {
  long ip = ipi;

  if(ipi<0) ip = pvslct(pBox);

  Box2 = generateBox(pBox);
    
  Interval *X1,*X2;
      
  X1 = pBox->getIntervalAddress();
  X2 = Box2->getIntervalAddress();

  double tmidpt = mid(X1[ip]);
	X1[ip].setright(tmidpt);
	X2[ip].setleft(tmidpt);

  stack_wrk->push(Box2);
  }

void intbis::bisect(intBox *pBox, intBox *pBox2)
  {
  long ip = pvslct(pBox);
    
  Interval *X1,*X2;
      
  X1 = pBox->getIntervalAddress();
  X2 = pBox2->getIntervalAddress();

  double tmidpt = mid(X1[ip]);
	X1[ip].setright(tmidpt);
	X2[ip].setleft(tmidpt);
  }


long intbis::pvslct(intBox *pBox)
  {
  long ip = -1, i;
  double wid = 1.0e-10, tmp;
  Interval *X = pBox->X;
  for(i=0;i<pBox->szX;i++)
    {
    tmp=width(X[i]);
    if(tmp>wid)
      {
	    ip =i; wid=tmp;
      }
    }

  return ip;
  }

// *****************************

void intbis::expand(intBox *varBox,double r)
  {
  double ci, lb, rb;
  Interval *x = varBox->X;
  for(long i=0;i<varBox->szX;i++)
    {
    ci=mid(x[i]);
    lb=x[i].getleft() -EPSBOX/2;
    rb=x[i].getright()+EPSBOX/2;
    x[i]=ci+2.0/r*(x[i]-ci);
    if(x[i].getleft()<lb) x[i].setleft(lb);
    if(x[i].getright()>rb) x[i].setright(rb);
    }
  }

// *****************************

void intbis::ftest(intBox *pBox, long& UNKNOWN, long &SIGRT, EquationSet *eqSet)
  {
  long retcon;
  
  retcon=hsing(pBox,eqSet);

  if(retcon>=iLUNKNW && retcon<=iUUNKNW) // roots inclusion
    {
    UNKNOWN = iTRUE;
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
  }
// *****************************

void intbis::solve()
  {
  long N = pBoxExtr->szX, ip;
  long nleaves = 0, indentCount = 1;
  bool branchDone = 0, branchSolve = 0, initSplit = 0;
  double r = 0.5, diamReal; // r is a parameter used in the expansion/deletion process
  long UNKNOWN, SIGRT;
  
  allocate(N);

  double error = EPSBOX * 2.5;

  long testFlg; // = 0, no more boxes in the stack; 1, boxes remain in the stack.
  long testNum=0; // total number of ftest call

  do{
    do{
      ftest(pBoxExtr,UNKNOWN,SIGRT,eqSetExtr);
      testNum++;

      diamReal = diamcp(N, pBoxExtr->X);

      if(UNKNOWN)
        {
        ip = pvslct(pBoxExtr);
        if(ip<0)
          {
          UNKNOWN = 0; SIGRT = 1;
          }
        else
          {
          bisect(pBoxExtr, ip);
          }

        if(diamReal >= EPSBOX) 
          testFlg = 1;
        else
          {
          testFlg = 0;
          stack_rst->push(pBoxExtr);
          }
        }

      if(!UNKNOWN)
        {
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
          if(!dupNode) stack_sln->push(pBoxExtr);
          }
        else
          {
          stack_rst->push(pBoxExtr);
          }
        }
      }while(testFlg);
    
    if(stack_wrk->is_empty()==false)
      {
#ifdef _DEBUG
      nleaves++;
      long nItems = stack_sln->numberOfItems();
#endif
      pBoxExtr = stack_wrk->pop();
      testFlg=1;
      }
    else
      {
      long noteGen = stack_rst->numberOfItems()+stack_sln->numberOfItems();
      long NRTest  = testNum;
      printf("\n\n*** Number of generated nodes \t\t %ld ***\n\n",noteGen);
      printf("\n\n*** Number of root inclusion test \t %ld ***\n",NRTest);
      testFlg = 0;
      newtonRaphson(pBoxExtr, eqSetExtr);
      }
    }while(testFlg);
  }

// ***************
long diagFlag=0;

long intbis::newtonRaphson(intBox *pBox, EquationSet *eqSet, double tol)
  {
  long i, n, retcon = 1, iter, iterMax=30;
  double *fpt, *vv, *jac = pBox->dfpt;
  try
    {
    n = pBox->szX;

    fpt  = &pwork[0];
    vv    = &fpt[n];

    while(!stack_sln->is_empty())
      {
      pBox = stack_sln->pop();
      for(i=0; i<n; ++i)
        {
        wiv2[i] = pBox->X[i];
        pBox->Xpt[i] = mid(pBox->X[i]);
        wiv1[i] = pBox->Xpt[i];
        }
      pBox->setValue(n, wiv1);
      iter = 0;
      while(iter<iterMax)
        {
        iter++;
        eqSet->fun(pBox);
        for(i=0; i<n; ++i)
          {
          if(fabs(pBox->f[i].getleft())>tol || fabs(pBox->f[i].getright())>tol)
            {
            retcon = 0; // not yet
            break;
            }
          }
        if(retcon!=0)
          {// tight enough
          for(i=0; i<n; ++i) pBox->Xpt[i] = mid(wiv1[i]);
          pBox->setValue(n, wiv2);
          stack_wrk->push(pBox);
          break;
          }
        else
          {
          long info, j;
          double dxMax;
          eqSet->jac(pBox);
          for(i=0; i<n; ++i) vv[i] = mid(pBox->f[i]);
          for(i=0; i<n*n; ++i) jac[i] = mid(pBox->df[i]);
          info = LuDcmp(jac,n,ipvt,vv);
          inverse(n, jac, ymat, vv, ipvt);
  #ifdef _DEBUG
          if(diagFlag)
            {
            long m;
            printf("\nCheck J*J^-1 =[\n");
            for(i=0; i<n; ++i)
              {
              for(j=0; j<n; ++j)
                {
                dxMax = 0.0;
                for(m=0; m<n; ++m)
                  dxMax += mid(pBox->df[i*n+m]) * ymat[m*n+j];
                printf("\t%-16.8g",dxMax);
                }
              printf("\n");
              }
            }
  #endif
          dxMax = 0.0;
          for(i=0; i<n; ++i)
            {
            vv[i] = 0.0;
            for(j=0; j<n; ++j) vv[i] += ymat[i*n+j] * mid(pBox->f[j]);
            pBox->Xpt[i] -= vv[i];
            if(pBox->Xpt[i]<wiv2[i].getleft())
              pBox->Xpt[i] = wiv2[i].getleft();
            if(pBox->Xpt[i]>wiv2[i].getright())
              pBox->Xpt[i] = wiv2[i].getright();
            wiv1[i] = pBox->Xpt[i];
            if(dxMax<fabs(vv[i]))
              dxMax = fabs(vv[i]);
            }
          if(dxMax<tol)
            {
            retcon = 1;
            for(i=0; i<n; ++i) pBox->Xpt[i] = mid(wiv1[i]);
            pBox->setValue(n, wiv2);
            stack_wrk->push(pBox);
            break;
            }
          else
            pBox->setValue(n, wiv1);
          }
        }
      if(retcon==0)
        {// newton did not converge, but solution does exist
        retcon = 1;
        for(i=0; i<n; ++i) pBox->Xpt[i] = mid(wiv1[i]);
        pBox->setValue(n, wiv2);
        stack_wrk->push(pBox);
        }
      }

    n = stack_wrk->numberOfItems();
    while(!stack_wrk->is_empty())
      {
      pBox = stack_wrk->pop();
      stack_sln->push(pBox);
      }
    }
  catch(...)
    {
    printf("\nwrong in newton\n");
    }
  return retcon;
  }

// **********************

long intbis::preconditionC(intBox *pBox, Interval *matr)
  {
  long info=0, n=pBox->szX;
  long dummy=0, i;

  double *cwrk, *vv, *jac = pBox->dfpt;

  cwrk  = &pwork[0];
  vv    = &cwrk[n];

  for(i=0; i<n*n; ++i)
    jac[i] = mid(pBox->df[i]);

  info = LuDcmp(jac,n,ipvt,vv);

  if(info!=0)
    return info;

  inverse(n, jac, ymat, cwrk, ipvt);
  
  return info;
  }
// *****************************

long intbis::singleStepGS(long i, intBox *pBox, intBox *pBoxthin, Interval &imag1, Interval &imag2)
  {
  long info,j,iN, N=pBox->szX;
  Interval *fint = pBoxthin->f, *x=pBox->X;
  double *xpt=pBoxthin->Xpt;
  //                   y_if + Sum_j(j!=i) y_iA_j(X_j- xpt_j)
  //  wiv1_i = xpt_i - -----------------------------------------
  //                             y_iA_i
  //
  //
  //  while :     y_iA_j = tmpm1[i+j*n]
  //
  ///////////////////////////////////////
  iN=i*N;

  //
  // single step of interval Gauss-Seidel (i.e. Hansen-Sengupta) starts here:
  //
  wiv1[i]=0.0;
  for(j=0;j<N;j++)
    {
    if((fint[j].getleft()==0.0 || fint[j].getright()==0.0) && ymat[j+iN]==0.0) 
      continue;
    wiv1[i] += ymat[j+iN] * fint[j];
    }
  //
  // this is the first term on the nomenator (y_if)
  //
  ///////////////////////////////////////

  //
  //one step of interval Newton method
  for(j=0;j<N;j++)
    {
    if(j!=i) wiv1[i] += tmpm1[j+iN]*(x[j]-xpt[j]);
    }
  //
  // this is the second term on the nomenator.
  //
  //////////////////////////////////////

  //
  info=dvsbin(xpt[i],x[i],wiv1[i],tmpm1[iN+i],imag1,imag2); 
  //
  // this is to perform xpt_i - Q_i/D_i, where 
  //
  //                             Q_i is the nomenator;
  //                             D_i is the denomenator.
  //
  // single step of interval Gauss-Seidel (i.e. Hansen-Sengupta) ends here:
  //
  //////////////////////////////////////
  return info;
  }

// ****************************

long intbis::hsing(intBox *pBox, EquationSet *eqSetI)
  {
  long info, retcon, i, j, k, N = pBox->szX;
  double *xpt = pBox->Xpt;
  Interval *x = pBox->X, *fint = pBox->f, *jacint = pBox->df;
  intBox *pBoxL;

  bool locIterDone = 0, monotonic = 1;
  long maxlocIter = 3, locIter = 0;
  double decreaseRate = 1.0, funcWidth;

  do
    {
    locIterDone = 1;
    locIter++;

    retcon = 1;

    eqSetI->fun(pBox);

    funcWidth = 0;
    for(i=0; i<N; ++i)
      {
      if(!(fint[i].getleft()<=0.0 && fint[i].getright()>=0.0))
        {
        retcon = 20;
        break;
        }
      else
        funcWidth += width(fint[i]);
      }

    if(retcon != 20)
      {
      if(funcWidth<EPSFUN)
        {
        retcon = 10;
        break;
        }
      else if(funcWidth>1e+10)
        {
        retcon = 1;
        return retcon;
        }
      }
    if(retcon != 20)
      {
      eqSetI->jac(pBox);

      if(do_mono)
        {
        retcon = componentWiseGS(pBox, eqSetI);
        if(eqSetI->isTensProvided() && retcon!=20 && diamcp(N, x)>0.5)
          {// this helps when the box is still wide
          k=0;
          double widXk = width(x[k]);
          for(i=1; i<N; ++i)
            {
            if(widXk>width(x[i]))
              {
              k = i;
              widXk = width(x[i]);
              }
            }
          for(i=0; i<N; ++i)
            {
            retcon = componentWiseFi(pBox, eqSetI, &jacint[i*N+k], i, k);
            if(retcon==20)
              break;
            }
          }
        }
      }
    }while(!locIterDone);

  if(retcon==20 || retcon==10) 
    return retcon;

  // obtain point value function evaluation

  pBoxL = generateBox(pBox);

  for(i=0;i<N;++i) 
    {
    xpt[i] = mid(x[i]);
    wiv2[i] = xpt[i];
    }

  pBoxL->setValue(N, wiv2,xpt);
  eqSetI->fun(pBoxL);
  eqSetI->jac(pBoxL);
 
  bool cntn=true;   

  Interval imag1, imag2;

  if(N>1)
    {
    info = preconditionC(pBoxL);
    }
  else
    {
    info = 0;
    }

  if(info!=0 && info!=23)
    {
    retcon = 1;
    return retcon;
    }

  Interval tmp;
  for(i=0;i<N;++i)
    {
    for(j=0;j<N;j++)
      {
      tmpm1[j+i*N]=0.0;
      for(k=0;k<N;++k)
        {
        if((jacint[j+k*N].getleft()==0.0 || jacint[j+k*N].getright()==0.0) && ymat[k+i*N]==0.0)
          continue;
        if(jacint[j+k*N].getleft()==0 && jacint[j+k*N].getright()==0.0)
          continue;
        tmp = ymat[k+i*N]*jacint[j+k*N];
        if(tmp.getleft()==0.0 && tmp.getright()==0.0)
          continue;
        tmpm1[j+i*N] += tmp;
        }
      }
    }

  long sameWidth = 0;
  for(i=0;i<N;i++)
    {
    info = singleStepGS(i,pBox, pBoxL,imag1,imag2);

    if(info==0)
      {// no overlap
      retcon = 20;
      return retcon;
      }
    else if(info==2)
      {   //jacobian include zero, then there are 2 imagine
      
      x[i]=imag2;// save one part

      Box2 = generateBox(pBox);
      stack_wrk->push(Box2); 
      x[i]=imag1;
      
      //retcon = 1;
      //return retcon;
      }

    if((imag1.getleft()<=x[i].getleft()) ||(imag1.getright()>=x[i].getright()))
      {
      cntn=false;  //no division 
      }

    if(imag1.getleft()==x[i].getleft() && imag1.getright()==x[i].getright())
      sameWidth++;
  
    //xboxWidthi = width(x[i]);

    x[i]=imag1;    //update the box on current dimension
    }

  if(cntn)        //unique root
    retcon=11;
  else            //unknow unique root 
    {
    retcon=1;
    double Xwidth = diamcp(N, x);
    if(Xwidth<EPSBOX)
      retcon = 10;
    else if(Xwidth<EPSBOX*10*N && sameWidth==N)
      retcon = 10;
    }

  stack_rst->push(pBoxL);
  
  return retcon;
  }

// *****************************
long intbis::componentWiseGS(intBox *pBox, EquationSet *eqSet, Interval *func, Interval *jacb)
  {
  long k, i, j, n=pBox->szX, retcon = 1;
  func = pBox->f;
  jacb = pBox->df;
  intBox *pBoxL;

  double widXk, xmid;

  Interval *x = pBox->X;
  for(i=0; i<n; ++i) wiv2[i] = x[i];

  Interval dfdx, xn, funci;

  pBoxL = generateBox(pBox);

  for(k=0; k<n; ++k)
    {// variable index
    widXk = width(x[k]);
    //if(widXk<EPSBOX) continue;

    for(i=0; i<n; ++i)
      {// function index
      dfdx = jacb[i*n+k];
      if(!(dfdx.getleft()<=0.0 && dfdx.getright()>=0.0))
        {
        for(j=0; j<n; ++j) wiv1[j] = x[j];

        xmid = mid(x[k]);
        wiv1[k] = xmid;
        pBoxL->setValue(n, wiv1, 0);

        eqSet->fun(pBoxL, i);
        funci = pBoxL->f[i];

        xn = xmid - funci/dfdx;

        if(!(xn.getleft()>x[k].getright() || xn.getright()<x[k].getleft()))
          {
          if(xn.getleft()>x[k].getleft() || xn.getright()<x[k].getright())
            {
            double left = xn.getleft();
            if(left<x[k].getleft())
              left = x[k].getleft();
            double right = xn.getright();
            if(right>x[k].getright())
              right = x[k].getright();
            x[k].setInterval(left,right);
            }
          }
        else
          {
          retcon = 20; i = k = n;
          }
        }
      else
        {
        if(zeroContain(dfdx) && width(dfdx)>1e-100 && k<pBox->szX)
          {
          Interval imag1, imag2;
          xmid = mid(x[k]);

          for(j=0; j<n; ++j) wiv1[j] = x[j];

          xmid = mid(x[k]);
          wiv1[k] = xmid;
          pBoxL->setValue(n, wiv1, 0);

          eqSet->fun(pBoxL, i);
          funci = pBoxL->f[i];

          long info = dvsbin(xmid,x[k],funci,dfdx,imag1,imag2); 
          if(info==2)
            {   //jacobian include zero, then there are 2 imagine
            
            x[k]=imag2;// save one part
            Box2 = generateBox(pBox);
            stack_wrk->push(Box2); 
            x[k]=imag1;          
            //retcon = 1;
            }
          else if(info==0)
            {
            retcon = 20; i = k = n;
            }
          }
        }
      }
    }

  if(retcon!=20)
    retcon = 1;

  stack_rst->push(pBoxL);

  return retcon;
  }

// *****************************

long intbis::componentWiseFi(intBox *pBox, EquationSet *eqSet, Interval *Jaci, long i, long k)
  {
  long retcon = 1, n=pBox->szX, j;
  Interval *x=pBox->X, funci, jaci;

  double xmid;
  intBox *pBoxL;
  Interval JikTaylor, xTaylor, One(0,1), xDiff, FiTaylor;

  pBoxL = generateBox(pBox);

  xmid = mid(x[k]);
  xDiff = x[k]-xmid;
  xTaylor = xmid + One*xDiff;
  for(j=0; j<n; ++j) wiv1[j] = x[j];
  wiv1[k] = xmid;

  pBoxL->setValue(n, wiv1,0);
  eqSet->fun(pBoxL, i);
  eqSet->jac(pBoxL, i, k);

  funci = pBoxL->f[i];
  jaci  = pBoxL->df[i*n+k];

  wiv1[k] = xTaylor;
  pBoxL->setValue(n, wiv1,0);
  eqSet->tens(pBoxL, i, k);

  Interval D2fi = pBoxL->ts->Taylor_fatParam;
  Interval D3fi = pBoxL->ts->Remain_fatParam;
  JikTaylor = jaci + xDiff * (D2fi + xDiff/2.0*D3fi);
  if(pBoxL->ts->isRemainConst)
    {
    parametric_poly(pBoxL, eqSet, i, k, &FiTaylor, xDiff, funci, jaci, D2fi, D3fi);
    }
  else
    {
    FiTaylor = funci + xDiff * JikTaylor;
    }

  if(FiTaylor.getleft()>0.0 || FiTaylor.getright()<0.0)
    retcon = 20;
  else
    {
    if(FiTaylor.getleft()>pBox->f[i].getleft())
      pBox->f[i].setleft(FiTaylor.getleft());

    if(FiTaylor.getright()<pBox->f[i].getright())
      pBox->f[i].setright(FiTaylor.getright());
    }

  if(JikTaylor.getleft() >Jaci->getleft())
    Jaci->setleft (JikTaylor.getleft() );

  if(JikTaylor.getright()<Jaci->getright())
    Jaci->setright(JikTaylor.getright());

  stack_rst->push(pBoxL);

  return retcon;
  }

// *****************************

long intbis::parametric_poly(intBox *pBox, EquationSet *eqSet, long idxF, long idxx,
                                   Interval *Fi, Interval xi, 
                                   Interval f0, Interval Df0, Interval D2f0, Interval D3f0)
  {
  intBox *pBoxL;

  pBoxL = generateBox(pBox);

  Interval DFDx = Df0 + xi*(D2f0 + xi/2.0*D3f0);
  if(DFDx.getleft()>0.0 || DFDx.getright()<0.0)
    {
    Interval tU, tL;
    long n = pBox->szX, i;
    for(i=0; i<n; ++i) wiv1[i] = pBox->X[i];
    tU = pBox->X[idxx].getright();
    wiv1[idxx] = tU;
    pBoxL->setValue(n, wiv1,0);
    eqSet->fun(pBoxL, idxF);
    tU = pBoxL->f[idxF];
    tL = pBox->X[idxx].getleft();
    wiv1[idxx] = tL;
    pBoxL->setValue(n, wiv1,0);
    eqSet->fun(pBoxL, idxF);
    tL = pBoxL->f[idxF];
    if(DFDx.getleft()>0.0)
      Fi->setInterval(tL.getleft(), tU.getright());
    else
      Fi->setInterval(tU.getleft(), tL.getright());
    }
  else
    {
    if(D3f0.getleft()==0.0 && D3f0.getright()==0.0 && (D2f0.getleft()>0.0 || D2f0.getright()<0.0))
      {
      double left, right;
      Interval x1 = 0.0 - Df0/D2f0, xL = xi.getleft(), xR=xi.getright(), F1, F2;
      F1 = f0 + xL*(Df0 + xL/2.0*D2f0);
      F2 = f0 + xR*(Df0 + xR/2.0*D2f0);
      if(x1.getleft()>=xL.getright() && x1.getright()<=xR.getleft())
        {
        Interval Fm = f0 + x1*(Df0 + x1/2.0*D2f0);
        left = Fm.getleft();
        if(left>F1.getleft()) left = F1.getleft();
        if(left>F2.getleft()) left = F2.getleft();
        right = Fm.getright();
        if(right<F1.getright()) right = F1.getright();
        if(right<F2.getright()) right = F2.getright();
        }
      else
        {
        left = F1.getleft();
        if(left>F1.getleft()) left = F2.getleft();
        right = F2.getright();
        if(right<F1.getright()) right = F1.getright();
        }
      Fi->setInterval(left, right);
      }
    else
      *Fi = f0 + xi*(Df0 + xi/2.0*(D2f0 + xi/3.0*D3f0));
    }

  stack_rst->push(pBoxL);
  return 0;
  }
