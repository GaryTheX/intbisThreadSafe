#pragma once
class Interval
  {
  private:
    double left, right;

  public:
    Interval(void){left=right=0;};
    Interval(double d){left=right=d;};
    Interval(double left, double right);
    virtual ~Interval(void);

    virtual void setInterval(double l, double r)
      {
      this->left = l; this->right = r;
      }
    virtual void setleft(double x){left = x;};
    virtual void setright(double x){right = x;};
    virtual double getleft(){return left;};
    virtual double getright(){return right;};

    void printInt();

    virtual Interval &operator = (const Interval&);
    virtual Interval &operator = (double);
    virtual Interval &operator = (long);

    friend long sign(const Interval&);
    friend double mid(const Interval&);
    friend double mid(const double&);
    friend bool zeroContain(const Interval&);
    friend bool NeqZero(const Interval&);

    virtual void  operator +=(const Interval&);
    virtual void  operator -=(const Interval&);
    virtual void  operator +=(double);
    virtual void  operator -=(double);
    virtual void  operator +=(long);
    virtual void  operator -=(long);
    virtual void  operator *=(const Interval&);
    virtual void  operator *=(double);

    friend Interval operator + (const Interval&,const Interval&);
    friend Interval operator + (const Interval&,double);
    friend Interval operator + (double, const Interval&);
    friend Interval operator + (const Interval&,long);
    friend Interval operator + (long, const Interval&);

    friend Interval operator - (const Interval&,const Interval&);
    friend Interval operator - (const Interval&,double);
    friend Interval operator - (double,const Interval&);
    friend Interval operator - (const Interval&,long);
    friend Interval operator - (long,const Interval&);

    friend Interval operator * (const Interval&,double);         
    friend Interval operator * (double, const Interval&);
    friend Interval operator * (const Interval&,const Interval&);
    friend Interval operator * (const Interval&,long);         
    friend Interval operator * (long, const Interval&);

    friend Interval operator / (const Interval&,const Interval&);
    friend Interval operator / (const Interval&,double);
    friend Interval operator / (double,const Interval&);
    friend Interval operator / (const Interval&,long);
    friend Interval operator / (long,const Interval&);

    friend double width(const Interval& x){return x.right - x.left;}
    friend double width(const double&){return 0.0;}
    friend double diamcp(long N, Interval *x);

    friend bool operator ==(const Interval&, const Interval&);
    friend bool operator ==(const Interval&, double);
    friend bool operator ==(double, const Interval&);
    friend bool operator ==(const Interval&, long);
    friend bool operator ==(long, const Interval&);

// --------------------------
    friend long operator & (const Interval&,const Interval&);
    friend int operator > (const Interval&,const Interval&);
    friend int operator > (const Interval&,double);
    friend int operator > (double, const Interval&);
    friend int operator < (const Interval&,const Interval&);
    friend int operator < (const Interval&,double);
    friend int operator < (double,const Interval&);
    friend int operator <=(const Interval&,const Interval&);
    friend int operator <=(double,const Interval&);
    friend int operator <=(const Interval&, double);

    friend void xdiv(long &xcase,const Interval &x,const Interval &y,
	            Interval &r1,Interval &r2);

    friend void xsclsb(long& xcase,double xpt,const Interval &image_xpt,
	            Interval &image);
    friend void xint(long& xcase,const Interval &a,Interval &b);
// ---------------------------------

    friend Interval ilog(const Interval&);
    friend double   ilog(const double&);
    friend Interval ilog10(const Interval&);
    friend double   ilog10(const double&);
    friend Interval iexp(const Interval&);
    friend double   iexp(const double&);
    friend Interval power(const Interval& x, long n);
    friend Interval iipow(const Interval& x, Interval y);
    friend Interval iipow(const Interval& x, double y);
    friend Interval isqrt(const Interval&);
    friend double   isqrt(const double&);

    friend Interval ixlnx(const Interval&, long np=3);
    friend double   ixlnx(const double&, long np=3);
    friend Interval ixexpx(const Interval&, long np=3);
    friend double   ixexpx(const double&, long np=3);

  };

