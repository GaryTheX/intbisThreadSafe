#include <stdio.h>
#include <math.h>
#include "Interval.h"
//#include "numeric.h"

extern "C"
  {
  extern void ROUNDDOWN(double& left);
  extern void ROUNDUP(double& right);
  extern double gettiny();
  }
const double huge_val = 1.0e20;
const double negInfy  = -230.0;

Interval::Interval(double left, double right)
  {
  this->left  = left ;
  this->right = right;
  }

// *********************************************

Interval::~Interval(void)
  {
  }

// *********************************************

void Interval::printInt()
  {
  printf("\n\t[%-16.8g, %-16.8g]\n",left, right);
  }

// *********************************************

long sign(const Interval& x)
  {
  long retval = 0;
  if(x.left > 0.0)
    {
    retval = 1;
    }
  else if(x.right < 0.0)
    {
    retval = -1;
    }
  return retval;
  }

// *********************************************

double mid(const Interval &x)
  {
  double retval;
  if(x.left>=0.0 && x.right<=gettiny())
    return 0.0;
  else    
    return retval = 0.5*(x.left+x.right);
  }

// *********************************************

double mid(const double &x)
  {
  return x;
  }

bool zeroContain(const Interval &x)
  {
  bool within = false;
  double xL = x.left;
  double xR = x.right;
  //if((-xL<=gettiny() && xL<=0.0) &&
  //   (xR<=gettiny()  && xR>=0.0))
  if(xL<0.0) ROUNDUP(xL);
  if(xR>0.0) ROUNDDOWN(xR);
  if(-xL>gettiny() && xR>gettiny())
    within = true;
  return within;
  }

bool NeqZero(const Interval &x)
  {
  double thetiny = gettiny();
  double xL = x.left;
  double xR = x.right;
  if(xL<-thetiny || xR>thetiny)
    return true;
  else
    return false;
  }

// ***************************************************

double diamcp(long N,Interval *x)
{
  double diam,w,big,temp;
  diam=0;

  for(long i=0;i<N;i++)
    {
    w=width(x[i]);
    big = fabs(x[i].left);
    temp= fabs(x[i].right);
    if(big<temp) big = temp;
    if(big<1.0) big = 1.0;
    temp = w/big;
    if(diam<temp) diam = temp;
    }
  return diam;
}
// *********************************************

void Interval::operator *=(const Interval &x)
  {
  Interval t1(left, right), t2;
  t2 = t1 * x;
  left = t2.left;
  right= t2.right;
  }
// *********************************************

void Interval::operator *=( double x)
  {
  Interval t1(left, right), t2;
  t2 = t1 * x;
  left = t2.left;
  right= t2.right;
  }

// *********************************************

void Interval::operator +=(const Interval &x)
  {
  left += x.left;
  right+= x.right;
  ROUNDDOWN(left);
  ROUNDUP(right);
  }

// *********************************************

void Interval::operator -=(const Interval &x)
  {
  left -= x.right;
  right-= x.left;
  ROUNDDOWN(left);
  ROUNDUP(right);
  }
// *********************************************

void Interval::operator +=(double x)
  {
  left += x;
  right+= x;
  ROUNDDOWN(left);
  ROUNDUP(right);
  }

// *********************************************

void Interval::operator -=(double x)
  {
  left -= x;
  right-= x;
  ROUNDDOWN(left);
  ROUNDUP(right);
  }
// *********************************************

void Interval::operator +=(long x)
  {
  left += x;
  right+= x;
  ROUNDDOWN(left);
  ROUNDUP(right);
  }

// *********************************************

void Interval::operator -=(long x)
  {
  left -= x;
  right-= x;
  ROUNDDOWN(left);
  ROUNDUP(right);
  }

// *********************************************
Interval& Interval::operator = (const Interval &x)
{
  left=x.left;right=x.right;
  return *this;
}

Interval &Interval::operator = (double x)
{
  left=right=x;
  return *this;
}

// *********************************************

Interval &Interval::operator = (long x)
{
  left=right=x; 
  return *this;
}

int operator > (const Interval &x,const Interval &y)
{
  //return x.right>y.right;
  return x.left > y.right;
}

int operator > (const Interval &x, double y)
{
  return x.left>y;
}
int operator > (double y, const Interval &x)
{
  return y>x.right;
}
int operator < (const Interval &x,const Interval &y)
{
  //return x.right>y.right;
  return x.right < y.left;
}

int operator < (const Interval &x, double y)
{
  return x.right<y;
}
int operator < (double y, const Interval &x)
{
  return y<x.left;
}

int operator <= (const Interval &x,const Interval &y)
{
  //return (x.left>=y.left && x.right<=y.right);
  return (x.right <= y.left);
}

int operator <= (double x,const Interval &y)
{
  //return x>=y.left && x<=y.right;
  return x<=y.left;
}

int operator <= (const Interval &x, double y)
{
  //return (x.right<=y || x.left<=y);
  return x.right<=y;
}

// *********************************************

Interval operator + (const Interval& x, const Interval& y)
  {
  double left = x.left  + y.left;
  if(y.left) ROUNDDOWN(left);
  double right= x.right + y.right;
  if(y.right) ROUNDUP(right);
  return Interval(left,right);
  }

// *********************************************

Interval operator + (const Interval& x, double y)
  {
  double left = x.left  + y;
  if(y) ROUNDDOWN(left);
  double right= x.right + y;
  if(y) ROUNDUP(right);
  return Interval(left,right);
  }

// *********************************************

Interval operator + (double y, const Interval& x)
  {
  double left = x.left  + y;
  if(y) ROUNDDOWN(left);
  double right= x.right + y;
  if(y) ROUNDUP(right);
  return Interval(left,right);
  }
// *********************************************

Interval operator + (const Interval& x, long y)
  {
  double left = x.left  + y;
  if(y) ROUNDDOWN(left);
  double right= x.right + y;
  if(y) ROUNDUP(right);
  return Interval(left,right);
  }

// *********************************************

Interval operator + (long y, const Interval& x)
  {
  double left = x.left  + y;
  if(y) ROUNDDOWN(left);
  double right= x.right + y;
  if(y) ROUNDUP(right);
  return Interval(left,right);
  }

// *********************************************

Interval operator - (const Interval &x, double y)
  {
  double left = x.left - y;
  if(y) ROUNDDOWN(left);
  double right= x.right- y;
  if(y) ROUNDUP(right);
  return Interval(left,right);
  }

// *********************************************

Interval operator - (double y, const Interval &x)
  {
  double left = y - x.right; //x.left - y;
  if(x.right) ROUNDDOWN(left);
  double right= y - x.left; //x.right- y;
  if(x.left) ROUNDUP(right);
  return Interval(left,right);
  }
// *********************************************

Interval operator - (const Interval &x, long y)
  {
  double left = x.left - y;
  if(y) ROUNDDOWN(left);
  double right= x.right- y;
  if(y) ROUNDUP(right);
  return Interval(left,right);
  }

// *********************************************

Interval operator - (long y, const Interval &x)
  {
  double left = y - x.right; //x.left - y;
  if(x.right) ROUNDDOWN(left);
  double right= y - x.left; //x.right- y;
  if(x.left) ROUNDUP(right);
  return Interval(left,right);
  }

// *********************************************

Interval operator - (const Interval &x, const Interval &y)
  {
  double left = x.left - y.right;
  if(y.right) ROUNDDOWN(left);
  double right= x.right- y.left;
  if(y.left) ROUNDUP(right);
  return Interval(left,right);
  }

// ***************************************************

Interval operator * (const Interval &x, double y)
  {
  double left, right;
  if(y==0.0 || (x.left == 0.0 && x.right == 0.0))
    {
    left = right = 0.0;
    }
  else
    {
    if(x.left==0.0)
      {
      if(y<0.0)
        {
        left = y * x.right;
        ROUNDDOWN(left);
        right = 0.0;
        }
      else
        {
        left = 0.0;
        right = y * x.right;
        ROUNDUP(right);
        }
      }
    else if(x.right==0.0)
      {
      if(y<0.0)
        {
        left = 0.0;
        right= y * x.left;
        ROUNDUP(right);
        }
      else
        {
        left = y * x.left;
        ROUNDDOWN(left);
        right = 0.0;
        }
      }
    else if(y>0.0)
      {
      left = y * x.left;
      right= y * x.right;
      ROUNDDOWN(left);
      ROUNDUP(right);
      }
    else
      {
      left = y * x.right;
      right= y * x.left;
      ROUNDDOWN(left);
      ROUNDUP(right);
      }
    }
  return Interval(left,right);
  }

// ***************************************************

Interval operator * (double y, const Interval &x)
  {
  Interval temp;
  temp = x * y;
  return temp;
  }

// ***************************************************

Interval operator * (const Interval &x,const Interval &y)
{
  double left,right,temp;
  if(sign(x)==1 && sign(y)==1)
    { // case 1. x,y both positive
    left = x.left * y.left;
    right= x.right* y.right;
    }
  else if(sign(x)==0 && sign(y)==1)
    { // case 2. y positive, x contain zero
    left = x.left * y.right;
    right= x.right* y.right;
    }
  else if(sign(x)==-1 && sign(y)==1)
    { // case 3. x negative, y positive
    left = x.left * y.right;
    right= x.right* y.left;
    }
  else if(sign(x)==1 && sign(y)==0)
    { // case 4. x positive, y contain zero
    left = x.right * y.left;
    right= x.right * y.right;
    }
  else if(sign(x)==-1 && sign(y)==0)
    { // case 5. x negative, y contain zero
    left = x.left * y.right;
    right= x.left * y.left;
    }
  else if(sign(x)==1 && sign(y)==-1)
    { // case 6. x positive, y negative
    left = x.right * y.left;
    right= x.left  * y.right;
    }
  else if(sign(x)==0 && sign(y)==-1)
    { // case 7, x contain zero, y negative
    left = x.right * y.left;
    right= x.left  * y.left;
    }
  else if(sign(x)==-1 && sign(y)==-1)
    { // case 8, both x and y negative
    left = x.right * y.right;
    right= x.left  * y.left;
    }
  else
    {
    left = x.left * y.right;
    temp = x.right* y.left;
    if(temp < left) left = temp;
    right = x.left * y.left;
    temp = x.right * y.right;
    if(temp > right) right = temp;
    }
  ROUNDDOWN(left);
  ROUNDUP(right);
  return Interval(left,right);
  }

// ***************************************************

Interval operator * (long y, const Interval &x)
  {
  Interval temp;
  double yy = y;
  temp = x * yy;
  return temp;
  }
// ***************************************************

Interval operator * (const Interval &x, long y)
  {
  Interval temp;
  double yy = y;
  temp = x * yy;
  return temp;
  }

// ***************************************************

Interval operator / (const Interval &x, long y)
  {
  Interval temp;
  double y_inv = 1.0/y;
  temp = x * y_inv;
  return temp;
  }

// ***************************************************

Interval operator / (long y, const Interval &x)
  {
  Interval temp;
  double yy = y;
  temp = yy / x;
  return temp;
  }

// ***************************************************

Interval operator / (const Interval &x, double y)
  {
  Interval temp;
  double y_inv = 1.0/y;
  temp = x * y_inv;
  return temp;
  }

// ***************************************************

Interval operator / (double y, const Interval &x)
  {
  double left, right, temp;
  if(sign(x)!=0)
    {
    if(y==0.0)
      {
      left = right = 0.0;
      }
    else
      {
      left = y/x.right;
      right= y/x.left;
      if(left>right)
        {
        temp = right;
        right= left;
        left = temp;
        }
      ROUNDDOWN(left);
      ROUNDUP(right);
      }
    }
  else
    {
    printf("divided by interval contains zero!\n");
    }
  return Interval(left,right);   
  }

// ***************************************************

Interval operator / (const Interval &x, const Interval &y)
  {
  Interval temp, retval;
  if(sign(y)!=0)
    {
    temp = 1.0/y;
    retval = x * temp;
    }
  else
    {
    retval.left = -1.0e100;
    retval.right= 1.0e100;
    }
  return retval;
  }

// ***************************************************

bool operator == (const Interval &x, long y)
  {
  bool retval;
  if(x.left==y && x.right==y) retval = 1;
  else retval = 0;
  return retval;
  }

// ***************************************************

bool operator == (long y, const Interval &x)
  {
  bool retval;
  if(x.left==y && x.right==y) retval = 1;
  else retval = 0;
  return retval;
  }
bool operator == (const Interval &x, double y)
  {
  bool retval;
  if(x.left==y && x.right==y) retval = 1;
  else retval = 0;
  return retval;
  }

// ***************************************************

bool operator == (double y, const Interval &x)
  {
  bool retval;
  if(x.left==y && x.right==y) retval = 1;
  else retval = 0;
  return retval;
  }

// ***************************************************

bool operator == (const Interval &x, const Interval &y)
  {
  bool retval;
  if(x.left==y.left && x.right==y.right) retval = 1;
  else retval = 0;
  return retval;
  }
// ***************************************************

long operator & (const Interval &x, const Interval &y)
  {  // interval intersection operation
  double xrp = x.right;
  double yrp = y.right;
  ROUNDUP(xrp);
  ROUNDUP(yrp);
  return !(xrp<y.left||x.left>yrp);
  }
/*
 * This routine performs the extended interval arithmetic division
 * x / y, and places the result(s) in r1 (and r2)
 * Arguments:
 *   x     is the numerator of the quotient. 
 *                (input)
 *   y     is the denominator of the quotient.
 *                (input)
 *   r1    is the interval quotient
 *   r2    is the second element of the interval quotient
 */
void xdiv(long &xcase,const Interval &x,const Interval &y,
	  Interval &r1,Interval &r2)
{
  Interval tmp;
  double tiny = gettiny();
  if(y.left>tiny || y.right<-tiny)
    {
    r1 = x/y; xcase=1;
    }
  else if(y.right>tiny && (y.left>=-tiny && y.left<=tiny))
    {                      // this means y.left==0)
    if(x.left>=0)
      {
      r1.left=x.left/y.right;// nominator is positive
      ROUNDDOWN(r1.left);
      r1.right=huge_val;
      xcase=2;
      }
    else
      {
      r1.right=x.right/y.right;  // nominator is not positive
      ROUNDUP(r1.right);
      r1.left=-huge_val;
      xcase=3;
      }
    }
  else if(y.left<-tiny && (y.right>=-tiny && y.right<=tiny))
    {
    if(x.right<=0)
      {
      r1.left=x.right/y.left;// nominator is negiitive
      ROUNDDOWN(r1.left);
      r1.right=huge_val;
      xcase=2;
      }
    else
      {
      r1.right=x.left/y.left;// nominator is not posiitive
      ROUNDUP(r1.right);
      r1.left=-huge_val;
      xcase=3;
      }
    }     
  else if(y.left<-tiny && y.right>tiny)
    {//the denominator include zero
    if(x.left>0)
      {  
      //if(x.left!=0)rounddown(); //nominator is positive
      r1.right=x.left/y.left;
      if(x.left!=0) ROUNDDOWN(r1.right);
      //if(x.left!=0)roundup();
      r2.left=x.left/y.right;
      if(x.left!=0) ROUNDUP(r2.left);
      r1.left=-huge_val;
      r2.right=huge_val;
      if(r2.left>r1.right)
        xcase=5;
      else 
        {
        r1.setleft(-huge_val);
        r1.setright(huge_val);
        xcase=4;
        }
      }
    else if(x.left<0)
      {
      //if(x.right!=0)rounddown();//nominator is not positive
      r1.right=x.right/y.right;
      if(x.right!=0) ROUNDDOWN(r1.right);
      //if(x.right!=0)roundup();
      r2.left=x.right/y.left;
      if(x.right!=0) ROUNDUP(r2.left);
      r1.left=-huge_val;
      r2.right=huge_val;
      if(r2.left>r1.right)
        xcase=5;
      else 
        {
        r1.setleft(-huge_val);
        r1.setright(huge_val);
        xcase=4;
        }
      }
    else
      {
      r1.setleft(-huge_val);
      r1.setright(huge_val);
      xcase=4;
      }
    }
  else
    {
    r1.setleft(-huge_val);
    r1.setright(huge_val);
    xcase=4;
    }
  }

//------------------------------------------
/*
 *  This routine subtracts an extended-value interval image_xpt, whose type
 *  is indicated by XCASE, from a point xpt.  The result is an extended-
 *  value interval image, whose type is stored in XCASE.
 */
void xsclsb(long& xcase,double xpt,const Interval &image_xpt,Interval &image)
{
  
  if(xcase==1)//the quantity is a single finite interval
    image=xpt-image_xpt;
  else if(xcase==2){//the quantity is a single semi-infinite interval
    //roundup();      // of the form [finite, +infinity]
    image.right=xpt-image_xpt.left;
    ROUNDUP(image.right);
    image.left=-huge_val;
    xcase=3;
  }
  else if(xcase==3){//the quantity is a single semi-infinite interval
    //rounddown();    //of the form [ - infinity , finite ]
    image.left=xpt-image_xpt.right;
    ROUNDDOWN(image.left);
    image.right=huge_val;
    xcase=2;
  }  
}

//------------------------------------------
/*  This routine intersects a finite interval A with an extended-value
*  interval B whose type is given by CASE.  It places the resulting
*  interval in B, if it is not null;  the routine resets the
*  variable CASE to indicate the type of B.
*/
double fmax(double a, double b)
  {
  if(a>b) return a;
  else return b;
  }

double fmin(double a, double b)
  {
  if(a>b) return b;
  else return a;
  }

void xint(long& xcase,const Interval &a,Interval &b)
{

  if(xcase==1)//the quantity is a single finite interval;
    if(a&b){
      b.left=fmax(a.left,b.left);
      b.right=fmin(a.right,b.right);
    }
    else xcase=0;
  else if(xcase==2)//if the quantity is a single semi-infinite interval
    if(a.right<b.left)xcase=0;// of the form [finite, +infinity];
    else {
      b.left=fmax(a.left,b.left);
      b.right=a.right;
      xcase=1;
    }
  else if(xcase==3)//if the quantity is a single semi-infinite interval
    if(a.left>b.right)xcase=0;//of the form [ - infinity , finite ];
    else {
      b.left=a.left;
      b.right=fmin(a.right,b.right);
      xcase=1;
    }
  else if(xcase==4){//if the quantity is the real line.
    b=a;
    xcase=1;
  }
}

//------------------------------------------
Interval ilog(const Interval &x)
  {
  double left,right;
  if(x.left<=0.0)
    left = negInfy;
  else 
    left =log(x.left);

  if(x.right<=0.0) 
    right = negInfy;
  else 
    right=log(x.right);

  ROUNDDOWN(left);
  ROUNDUP(right);
  return Interval(left,right);
  }

// ***************************************************

double ilog(const double &x)
  {
  double retval;
  if(x<=0.0) 
    retval = negInfy;
  else 
    retval = log(x);

  return retval;
  }
// ***************************************************

Interval ilog10(const Interval &x)
  {
  double left,right;

  if(x.left<=0.0)
    left = negInfy;
  else
    left = log10(x.left);

  if(x.right<=0.0)
    right = negInfy;
  else
    right = log10(x.right);

  ROUNDDOWN(left);
  ROUNDUP(right);

  return Interval(left,right);
  }

// ***************************************************

double ilog10(const double &x)
  {
  double retval;
  if(x<=0.0)
    retval = negInfy;
  else
    retval = log10(x);
  return retval;
  }
// ***************************************************
/*
Interval iexp(const Interval &x)
  {
  double left, right;
  
  double *ODF = Interval::ODF;
  Interval pval, G, tmp, xx, E14, A3;
  long ipass, np;
  ipass = 0;
  E14.setInterval(1.28402541668774,1.28402541668774);
  A3.setInterval(1.13314845306683, 1.13314845306683);

  while(ipass<2)
    {
    if(ipass==0) xx.setInterval(x.left, x.left);
    else xx.setInterval(x.right, x.right);
    if(xx < 1.0e-12) xx = 1.0e-12;
    if(xx > 708.0)   xx = 708.0;
    tmp = (4.0*xx);
    np = long(tmp.getleft());
    tmp = tmp/4.0;
    G   = x - tmp;
    
    pval = ODF[0] + G*(ODF[0] +
                    G*(ODF[1] +
                    G*(ODF[2] +
                    G*(ODF[3] +
                    G*(ODF[4] + 
                    G*(ODF[5] +
                    G*(ODF[6] +
                    G*(ODF[7] +
                    G*(ODF[8] +
                    G*(ODF[9] +
                    G*(ODF[10])))))))))));

    tmp = power(G, 12);
    tmp = ODF[11] * A3 * tmp;
    pval += tmp;
    pval = pval * power(E14,np);

    if(ipass==0) left = pval.getleft();
    else right = pval.getright();
    ipass++;
    }
  

  left = iexp(x.left);
  right = iexp(x.right);

  ROUNDDOWN(left);
  ROUNDUP(right);
  return Interval(left,right);
  }
  */

// ***************************************************

double iexp(const double &x)
{
    double y, ret;
    if(x>230.0) y = 230.0;
    else y = x;
    ret = exp(y);
    ret = (ret>0.0)? ret : 1.0e-100;
    return ret;
}
// ***************************************************
Interval iexp(const Interval &x)
  {
  double left, right, middle, exp_m;
  Interval xx, Left, Right;
  long ipass;
  ipass = 0;

  while(ipass<2)
    {
    if(ipass==0) xx.setInterval(x.left, x.left);
    else xx.setInterval(x.right, x.right);

    middle = mid(xx);
    exp_m = iexp(middle);
    if(ipass==0)
      Left  = exp_m*(1.0 + (xx-middle)*(1.0 + (xx-middle)*(0.5 + (xx-middle)*(1.0/3.0)))) + power(xx-middle,4)/24.0*iexp(xx.left);
    else
      Right = exp_m*(1.0 + (xx-middle)*(1.0 + (xx-middle)*(0.5 + (xx-middle)*(1.0/3.0)))) + power(xx-middle,4)/24.0*iexp(xx.right);
    ipass++;
    }

  left = Left.getleft();
  right = Right.getright();

  ROUNDDOWN(left);
  ROUNDUP(right);
  return Interval(left,right);
  }

// ***************************************************

Interval power(const Interval& x, long n)
  {
  Interval t1, t2, Fm, FmI;
  double x0, f0, fm;
  long np = n, m, mt, nt, i;
  bool neg = false;

  if(n==0)
    {
    Fm = 0.0;
    return Fm;
    }
  else if(n<0)
    {
    neg = true; np = -n;
    }
  
  FmI = x;
  for(i=1; i<np; ++i) FmI *= x;
  if(neg) FmI = 1.0/FmI;

  if((x.left<=0 && x.right>=0) || fabs(x.right-x.left)>0.1)
    return FmI;

  x0 = x.left;
  f0 = pow(x0,(double)np);

  Fm = f0;
  t1 = x - x0;
  for(m=1; m<=np; ++m)
    {
    mt = 1;
    nt = np;
    t2 = t1;
    for(i=1; i<=m; ++i) 
      {
      mt *= i;
      if(np!=m)
        nt *= m;
      t2 *= t2;
      }

    fm = nt * pow(x0, (double)(np-m));

    Fm += fm/mt * t2;
    }
  if(neg)
    {
    Fm = 1.0/Fm;
    ROUNDUP(Fm.right);
    ROUNDDOWN(Fm.left);
    }

  if(Fm.getleft()>=FmI.getleft() && Fm.getright()<=FmI.getright())
    return Fm;
  else
    return FmI;
  }

// ***************************************************

Interval iipow(const Interval& x, Interval y)
  {
  Interval tmpX, tmpY, ret, tmp;
  tmpX.setInterval(x.left,x.right);
  tmpY.setInterval(y.left,y.right);

  if(tmpX.left > tmpX.right || tmpX.left<0.0)//tmpX.left<=MTHCNS.ZERO[1])
    {
    printf("bad interval in iipow()!\n");
    //exit(1);
    }
  
  ret = ilog(tmpX);
  tmp = ret * tmpY;
  ret = iexp(tmp);
  return ret;
  }

// ***************************************************

Interval iipow(const Interval& x, double y)
  {
  Interval tmpY, ret;
  double err = 1.0e-100;
  if(fabs(y-(long)(y)) > err)
    {
    tmpY.setInterval(y,y);

    ret = iipow(x,tmpY);
    }
  else
    {
    ret = power(x, (long)(y));
    }
  return ret;
  }

// ***************************************************

Interval isqrt(const Interval &x)
  {
  Interval retval;
  double half = 0.5;
  retval = iipow(x,half);
  return retval;
  }

double isqrt(const double &x)
  {
  double retval = sqrt(fabs(x));
  return retval;
  }
// ***************************************************
double ixlnx(const double &x, long np)
{
    double retval = x*log(x);
    return retval;
}

Interval ixlnx(const Interval &x, long np)
  {
  Interval retval, Fpt, Dx, DX, myOne(0,1), Rf;
  double xpt = mid(x);
  double fpt = ixlnx(xpt);

  double dfpt, d2fpt, d3fpt;
  dfpt = log(xpt) + 1.0;
  d2fpt = 1.0/xpt;
  d3fpt = -d2fpt*d2fpt;
  
  Dx = x - xpt;
  Fpt = fpt + Dx*(dfpt + Dx/2.0*(d2fpt + Dx/3.0*d3fpt));

  DX = xpt + Dx * myOne;
  Rf = Dx*Dx*Dx*Dx/24.0 * 2.0/(DX*DX*DX);

  retval = Fpt + Rf;

  Fpt = x * ilog(x);

  if(retval.getleft()>Fpt.getleft() && retval.getright()<Fpt.getright())
    return retval;
  else
    return Fpt;
  }

// ***************************************************

Interval ixexpx(const Interval &x, long np)
  {
  Interval retval, Fpt, Dx, DX, myOne(0,1), Rf;
  double xpt = mid(x);
  double expx = iexp(xpt), fpt = xpt*expx;

  double dfpt, d2fpt, d3fpt;
  dfpt  = expx*(1.0+xpt);
  d2fpt = expx*(2.0+xpt);
  d3fpt = expx*(3.0+xpt);
  
  Dx = x - xpt;
  Fpt = fpt + Dx*(dfpt + Dx/2.0*(d2fpt + Dx/3.0*d3fpt));

  DX = xpt + Dx * myOne;
  Rf = power(DX,4)/24.0 * (4.0+x)*iexp(x);

  retval = Fpt + Rf;

  Fpt = x * iexp(x);

  if(retval.getleft()>Fpt.getleft() && retval.getright()<Fpt.getright())
    return retval;
  else
    return Fpt;
  }

double ixexpx(const double &x, long np)
  {
  double retval = x*iexp(x);
  return retval;
  }
