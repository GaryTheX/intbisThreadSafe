#include <iostream>
#include <ostream>
using namespace std;

#include "miscmach.h"
const double dmach1 = 2.225073858507201e-308;
const double dtest1 = 2.220446049250313e-016;

extern "C"
  {
  extern void USEINTLIB();
  extern void ROUNDDOWN(double& left);
  extern void ROUNDUP(double& right);
  double thetiny, theMXULP;
  extern long ludcmp(double *a, long n, long *indx, double *d, double *vv);  
	extern long LuDcmp(double *a, long n, long *indx, double *vv);
  extern void lubksb(long n,double *a,double *b,long *ipvt);
  extern void inverse(long n,double *a,double *b,double *wrk,long *ipvt);
  extern void tranverse(long n,double *a);
  extern double gettiny();

  }

void USEINTLIB()
  {
  theMXULP = dtest1;
  thetiny = dmach1/(1.0-3.0*theMXULP);
//#ifdef _DEBUG
  /*
  char mychar1[50]="MACHINE CONSTANTS FOR IEEE ARITHMETIC";
  char mychar2[50]="\t\t\t\t\tMACHINES AND 8087-BASED MICROS";
  char mychar3[50]="\t\t\t\t\tSUCH AS THE IBM PC AND AT&T 6300";
  char mychar4[50]="\t\t\t\t\t\tTHE LEAST SIGNIFICANT BYTE IS";
  cout << mychar1 << "\n" << mychar2 << "\n" << mychar3 << "\n" << mychar4 << "\n\n\t\t\t\t\t\t" << thetiny << "\n";
*/
  //#endif
  }

void ROUNDDOWN(double& x)
  {
  if(x>=thetiny)
    x = (1.0e0 - theMXULP) * x;
  else if(x<= -thetiny)
    x = (1.0e0 + theMXULP) * x;
  else if(x==0.0e0)
    x = -thetiny;
  else
    x = 0.0e0;
  }

void ROUNDUP(double& x)
  {
  if(x>=thetiny)
    x = (1.0e0 + theMXULP) * x;
  else if(x<=-thetiny)
    x = (1.0e0 - theMXULP) * x;
  else if(x>=0.0)
    x = thetiny;
  else
    x = 0.0e0;
  }

double gettiny()
  {
  return thetiny;
  }

#define TINY 1.0e-20;

long ludcmp(double *a, long n, long *indx, double *d, double *vv)
  {
  long i,imax=-1,j,k;
  double big,dum,sum=0.0,temp;
  *d=1.0;
  for (i=0;i<n;i++)
    {
    big=0.0;
    for (j=0;j<n;j++)
      {
      if ((temp=fabs(a[i*n+j])) > big) big=temp;
      }
    if (big == 0.0)
      return 1;//nrerror("Singular matrix in routine ludcmp")
    vv[i]=1.0/big;
    }
  for (j=0;j<n;j++)
    {
    for (i=0;i<j;i++)
      {
      sum=a[i*n+j];
      for (k=0;k<i;k++)
        sum -= a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
      }
    big=0.0;
    for (i=j;i<n;i++)
      {
      sum=a[i*n+j];
      for (k=0;k<j;k++)
        sum -= a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big)
        {
        big=dum;
        imax=i;
        }
      }
    if (j != imax)
      {
      for (k=0;k<n;k++)
        {
        dum=a[imax*n+k];
        a[imax*n+k]=a[j*n+k];
        a[j*n+k]=dum;
        }
      *d = -(*d);
      vv[imax]=vv[j];
      }
    indx[j]=imax;
    if (a[j*n+j] == 0.0) 
      a[j*n+j]=TINY;
    if (j != n-1)
      {
      dum=1.0/(a[j*n+j]);
      for (i=j+1;i<n;i++) a[i*n+j] *= dum;
      }
    }
  return 0;
  }
#undef TINY
// *****************************
//
long LuDcmp(double *a, long n, long *indx, double *vv)
  {
  long i,imax,j,k;
  double big,dum,sum,temp, sd;
  sd=1.0;

  for (i=0;i<n;i++)
    {
    big=0.0;
    for (j=0;j<n;j++)
      {
      if ((temp=fabs(a[i*n+j])) > big) 
        {
        big=temp; break;
        }
      }
    if (big == 0.0)
      return 1+i;//nrerror("Singular matrix in routine ludcmp")
    vv[i]=1.0/big;
    }

  for (j=0;j<n;j++)
    {
    for (i=0;i<j;i++)
      {
      sum=a[i*n+j];
      for (k=0;k<i;k++)
        sum -= a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
      }
    big=0.0;
    for (i=j;i<n;i++)
      {
      sum=a[i*n+j];
      for (k=0;k<j;k++)
        sum -= a[i*n+k]*a[k*n+j];
      a[i*n+j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big)
        {
        big=dum;
        imax=i;
        }
      }
    if (j != imax)
      {
      for (k=0;k<n;k++)
        {
        dum=a[imax*n+k];
        a[imax*n+k]=a[j*n+k];
        a[j*n+k]=dum;
        }
      sd = -(sd);
      vv[imax]=vv[j];
      }
    indx[j]=imax;
    //if (j != n-1)
      {
      if (a[j*n+j] == 0.0) 
        {
        return 1; // singular matrix, has to restart from different initial values
        }
      dum=1.0/(a[j*n+j]);
      for (i=j+1;i<n;i++) a[i*n+j] *= dum;
      }
    }
  return 0;
  }

//---------------------------------------
/*  this routine is the back sunstitution 
 *  solve the equation 
 *    a x= b
 *    place the result into b
 *  Argument
 *    n     the dimension of the vector
 *                 (input)
 *    a     the coefficent matrix (n*n)
 *                 (input)
 *    b     the right side vector (n*1)
 *                 (input/output)
 *    ipvt  pivot working vector (n*1)
 *                 (output)
 */
void lubksb(long n,double *a,double *b,long *ipvt)
{
  long i,k;
  double temp;

  for(i=0;i<n;i++)//get the exchanged vector
    if(ipvt[i]!=i){
      temp=b[i];
      b[i]=b[ipvt[i]];
      b[ipvt[i]]=temp;
    }

  for(i=0;i<n;i++)  //slove " L y=B "
    for(k=0;k<i;k++)
      b[i]-=a[i*n+k]*b[k];
  
  for(i=n-1;i>=0;i--){ //slove "U x=y "
    for(k=i+1;k<n;k++)
      b[i]-=a[i*n+k]*b[k];
    b[i]=b[i]/a[i*n+i];
  }
}

// *****************************
//---------------------------------------
/*  this routine is to get a tranverse matrix
 *  
 *    place the result back into a
 *  Argument
 *    n     the dimension of the vector
 *                 (input)
 *    a     the square matrix (n*n)
 *                 (input/output)
 */
void tranverse(long n,double *a)
{
  long i,j;
  double t1;
  for(i=0;i<n-1;i++)
    for(j=i+1;j<n;j++){
	t1=a[i*n+j];
	a[i*n+j]=a[j*n+i];
	a[j*n+i]=t1;
    }
}
void inverse(long n,double *a,double *b,double *c,long *ipvt)
  {
  long i,j;
  for(i=0;i<n;i++)
    {
    for(j=0;j<n;j++) //get one column vector of unit matrix
      {
      if(j==i)
        c[j]=1;
      else
        c[j]=0;
      }
    lubksb(n,a,c,ipvt); //slove a x= c, then c=x
    for(j=0;j<n;++j)
      {
      b[j*n+i]=c[j];
      }
    }
  }
/*
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
*/