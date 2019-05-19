
// this routine is the declaration of the first derivative class
#pragma once
#include <math.h>


template <class Type>
class Value_Derivatives{
public:
  int idx_Der, idx_Var;
  Type f;
  Type *df;
  int size;
  void init(Type,int);

public:
  Value_Derivatives<Type>(Type vf,int sz);
  //Value_Derivatives<Type>(int sz);
  Value_Derivatives<Type>(int N=1);
  Value_Derivatives<Type>(const Value_Derivatives<Type>&);
  ~Value_Derivatives<Type>(){delete [] df;}
  int getsize(){return size;}
  void resize(int N){init(0,N);}
  void setIdx_Var(int ix){idx_Var = ix;}
  void setIdx_Der(int id){idx_Der = id;}
  
  Value_Derivatives<Type>& operator = (const Value_Derivatives<Type>&);  
  Value_Derivatives<Type>& operator = (double); 

friend Value_Derivatives<Type> operator + (const Value_Derivatives<Type>&,double);
friend Value_Derivatives<Type> operator + (double,const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> operator + (const Value_Derivatives<Type>&,const Value_Derivatives<Type>&);

friend Value_Derivatives<Type> operator - (const Value_Derivatives<Type>&,double);
friend Value_Derivatives<Type> operator - (double,const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> operator - (const Value_Derivatives<Type>&,const Value_Derivatives<Type>&);

friend Value_Derivatives<Type> operator * (const Value_Derivatives<Type>&,double);         
friend Value_Derivatives<Type> operator * (double,const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> operator * (const Value_Derivatives<Type>&,const Value_Derivatives<Type>&);

friend Value_Derivatives<Type> operator / (const Value_Derivatives<Type>&,double);
friend Value_Derivatives<Type> operator / (double,const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> operator / (const Value_Derivatives<Type>&,const Value_Derivatives<Type>&);
 

  void  operator +=(const Value_Derivatives<Type>&);
  void  operator -=(const Value_Derivatives<Type>&);
  void  operator *=(const Value_Derivatives<Type>&);
  void  operator /=(const Value_Derivatives<Type>&);

friend Value_Derivatives<Type> sin(const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> cos(const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> arctan(const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> exp(const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> log(const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> log10(const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> pow(const Value_Derivatives<Type>&,double);
friend Value_Derivatives<Type> sqrt(const Value_Derivatives<Type>&);
friend Value_Derivatives<Type> sqr(const Value_Derivatives<Type>&);
//friend ostream& operator << (ostream&,const Value_Derivatives<Type>&); 

};

template <class Type>
inline  Value_Derivatives<Type>:: Value_Derivatives(int N)
{
  idx_Der = -1;
  idx_Var = 0;
  init(0,N);
}

/*
template <class Type>
inline  Value_Derivatives<Type>:: Value_Derivatives(int sz)
{
  init(0,0,sz);
  }*/

template <class Type>
inline  Value_Derivatives<Type>:: Value_Derivatives(Type vf,int sz)
{
  idx_Der = -1;
  idx_Var = 0;
  init(vf,sz);
}


template <class Type>
inline void  Value_Derivatives<Type>:: init(Type vf,int sz)
{
  df=new Type[size=sz];
  f=vf;
  for(int ix=0;ix<sz;++ix)
    {
    if(ix==idx_Var)
      df[ix] = 1.0;
    else
      df[ix] = 0.0;
    }
}

template <class Type>
inline  Value_Derivatives<Type>:: Value_Derivatives(const Value_Derivatives<Type> &vfvdf)
{
  idx_Der = vdvdf.idx_Der;
  idx_Var = vdvdf.idx_Var;
  init(vfvdf.f,vfvdf.size);
  for(int ix =0;ix<vdvdv.size; ++ix)
    df[ix] = vfvdf.df[ix];
}

template <class Type>
inline Value_Derivatives<Type>& Value_Derivatives<Type>::operator = (const Value_Derivatives<Type> &x)
{
  if(this==&x)return *this;
  idx_Der = x.idx_Der;
  idx_Var = x.idx_Var;
  f=x.f;
  for(int i=0;i<size;++i)
    df[i]=x.df[i]; 
  return *this;
}

template <class Type>
inline Value_Derivatives<Type>& Value_Derivatives<Type>::operator = (double x)
{
  f=x;
  for(int i=0;i<size;++i)
    df[i]=0;   
  return *this;
}



template <class Type>
inline Value_Derivatives<Type> operator + (const Value_Derivatives<Type> &x,const Value_Derivatives<Type> &y)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp(sz);
  tmp.f=x.f+y.f;
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.df[i]+y.df[i];
  return tmp;
}

template <class Type>
inline Value_Derivatives<Type> operator + (const Value_Derivatives<Type> &x,double y)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x.f+y;
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.df[i];
  return tmp;
}


template <class Type>
inline Value_Derivatives<Type> operator + (double x,const Value_Derivatives<Type> &y)
{
  int sz=y.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x+y.f;
  for(int i=0;i<sz;i++)
    tmp.df[i]=y.df[i];
  return tmp;

}


template <class Type>
inline Value_Derivatives<Type> operator - (const Value_Derivatives<Type> &x,const Value_Derivatives<Type> &y)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp; 
  tmp.f=x.f-y.f;
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.df[i]-y.df[i];
  return tmp;
}


template <class Type>
inline Value_Derivatives<Type> operator - (const Value_Derivatives<Type> &x,double y)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x.f-y;
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.df[i];
  return tmp;
}


template <class Type>
inline Value_Derivatives<Type> operator - (double x,const Value_Derivatives<Type> &y)
{
  int sz=y.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x-y.f;
  for(int i=0;i<sz;i++)
    tmp.df[i]=-y.df[i];
  return tmp;
}


template <class Type>
inline void Value_Derivatives<Type>::operator += (const Value_Derivatives<Type> &x)
{
  f+=x.f;
  for(int i=0;i<size;i++)
    df[i]+=x.df[i]; 
}


template <class Type>
inline void Value_Derivatives<Type>::operator -= (const Value_Derivatives<Type> &x)
{
  f-=x.f;
  for(int i=0;i<size;i++)
    df[i]-=x.df[i]; 
}

template <class Type>
inline Value_Derivatives<Type> sqrt(const Value_Derivatives<Type> &x)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=sqrt(x.f);
  for(int i=0;i<sz;i++)
    tmp.df[i]=0.5*x.df[i]/tmp.f;
  return tmp;
}


template <class Type>
inline Value_Derivatives<Type> exp(const Value_Derivatives<Type> &x)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=exp(x.f);
  for(int i=0;i<sz;i++)
    tmp.df[i]=tmp.f*x.df[i];
  return tmp;  
}

template <class Type>
inline Value_Derivatives<Type> log(const Value_Derivatives<Type> &x)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=log(x.f);
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.df[i]/tmp.f;
  return tmp;
}


template <class Type>
inline Value_Derivatives<Type> log10(const Value_Derivatives<Type> &x)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=log10(x.f);
  for(int i=0;i<sz;i++)
   tmp.df[i]=2.336*x.df[i]/tmp.f;
  return tmp;
}


template <class Type>
inline void Value_Derivatives<Type>::operator *= (const Value_Derivatives<Type> &y)
{
  for(int i=0;i<size;i++)
   df[i]=f*y.df[i]+df[i]*y.f;
  f=f*y.f;
}

template <class Type>
inline Value_Derivatives<Type> operator * (const Value_Derivatives<Type> &x,const Value_Derivatives<Type> &y)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x.f*y.f;
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.f*y.df[i]+x.df[i]*y.f;
  return tmp;
}

template <class Type>
inline Value_Derivatives<Type> operator * (const Value_Derivatives<Type> &x,double y)
{

  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x.f*y;
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.df[i]*y;

  return tmp;
}


template <class Type>
inline Value_Derivatives<Type> operator * (double x,const Value_Derivatives<Type> &y)
{
  int sz=y.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x*y.f;
  for(int i=0;i<sz;i++)
    tmp.df[i]=x*y.df[i];

  return tmp; 
}


template <class Type>
inline void Value_Derivatives<Type>::operator /= (const Value_Derivatives<Type> &y)
{ 
  f=f/y.f;
  for(int i=0;i<size;i++)
    df[i]=(df[i]-f*y.df[i])/y.f;
  
}


template <class Type>
inline Value_Derivatives<Type> operator / (const Value_Derivatives<Type> &x,const Value_Derivatives<Type> &y)
{ 
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x.f/y.f;
  for(int i=0;i<sz;i++)
    tmp.df[i]=(x.df[i]-tmp.f*y.df[i])/y.f;
  return tmp; 
}


template <class Type>
inline Value_Derivatives<Type> operator / (const Value_Derivatives<Type> &x,double y)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=x.f/y;
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.df[i]/y;

  return tmp;
}


template <class Type>
inline Value_Derivatives<Type> operator / (double x,const Value_Derivatives<Type> &y)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  Type y1=1/y.f;
  Type y2=pow(y1,2);
  tmp.f=x*y1;
  for(int i=0;i<sz;i++)
    tmp.df[i]=-x*ydf[i]*y2;
  return tmp;

}

template <class Type>
inline Value_Derivatives<Type> sin(const Value_Derivatives<Type> &x)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=sin(x.f);
  for(int i=0;i<sz;i++)
    tmp.df[i]=x.df[i]*cos(x.f);
  return tmp;
}


template <class Type>
inline Value_Derivatives<Type> cos(const Value_Derivatives<Type> &x)
{
  int sz=x.size;
  Value_Derivatives<Type> tmp;
  tmp.f=cos(x.f);
  for(int i=0;i<sz;i++)
    tmp.df[i]=-x.df[i]*sin(x.f);
  return tmp;
}

template <class Type>
inline Value_Derivatives<Type> pow(const Value_Derivatives<Type> &x,double p)
{ 

  int sz=x.size;
  Value_Derivatives<Type> tmp;
  
  tmp.f=pow(x.f,p);
  Type tmp2=p*pow(x.f,p-1);
  for(int i=0;i<sz;i++)
    tmp.df[i]=tmp2*x.df[i];
  return tmp;
}














































































































































