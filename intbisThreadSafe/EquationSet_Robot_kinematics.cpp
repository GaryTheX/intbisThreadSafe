#include "EquationSet_Robot_kinematics.h"


EquationSet_Robot_kinematics::EquationSet_Robot_kinematics(void)
  {
  pa = new double[19];
  pa[0] = 4.731e-3;
  pa[1] = -0.3578;
  pa[2] = -0.1238;

  pa[3] = -1.637e-3;
  pa[4] = -0.9338;
  pa[5] = 1.0;

  pa[6] = -0.3571;
  pa[7] = 0.2238;
  pa[8] = 0.7623;

  pa[9] = 0.2638;
  pa[10]= -0.7745e-1;
  pa[11]= -0.6734;

  pa[12] = -0.6022;
  pa[13] = 1.0;
  pa[14] = 0.3578;

  pa[15] = 4.731e-3;
  pa[16]= -0.7623;
  pa[17]= 0.2238;

  pa[18] = 0.3461;

  neq = 8;  
  icase = 1;

  }


EquationSet_Robot_kinematics::~EquationSet_Robot_kinematics(void)
  {
  delete[] pa;
  }

bool EquationSet_Robot_kinematics::obt_buildInBound(intBox *pBox)
  {
  pBox->initialize(neq);
  int i;
  for(i=0; i<neq; ++i)
    pBox->X[i].setInterval(-1.0, 1.0);

  return true;
  }

void EquationSet_Robot_kinematics::fun(intBox *pBox, long idx)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;

  long i, n = pBox->szX, is=0, ie = n;

  if(idx>=0)
    {
    is = idx; ie = is+1;
    }

  for(i=is; i<ie; ++i)
    {
    switch(i)
      {
      case 0:
        f[i] = pa[0]*x[0]*x[2] +pa[1]*x[1]*x[2] + pa[2]*x[0] + pa[3]*x[1] + pa[4]*x[3] + pa[5]*x[6] + pa[6];
        break;
      case 1:
        f[i] = pa[7]*x[0]*x[2]+pa[8]*x[1]*x[2] +pa[9]*x[0] + pa[10]*x[1] + pa[11]*x[3] + pa[12];
        break;
      case 2:
        f[i] = pa[13]*x[5]*x[7] + pa[14]*x[0] + pa[15]*x[1];
        break;
      case 3:
        f[i] = pa[16]*x[0] + pa[17]*x[1] + pa[18];
        break;
      case 4:
        f[i] = x[0]*x[0] + x[1]*x[1] - 1.0;
        break;
      case 5:
        f[i] = x[2]*x[2] + x[3]*x[3] - 1.0;
        break;
      case 6:
        f[i] = x[4]*x[4]+x[5]*x[5] - 1.0;
        break;
      case 7:
        f[i] = x[6]*x[6]+x[7]*x[7] - 1.0;
        break;
      default:
        break;
      }
    }
  }

void EquationSet_Robot_kinematics::jacDigit(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;
  Interval *df= pBox->df;
  long i, j, n = pBox->szX, is=0, ie = n, js=0, je=n;

  if(idxEq>=0)
    {
    is = idxEq; ie = is+1;
    }
  if(idxVr>=0)
    {
    js = idxVr; je = js+1;
    }

  for(i=is; i<ie; ++i)
    {
    switch(i)
      {
      case 0:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+0] = pa[0]*x[2]+pa[2];
        df[i*n+1] = pa[1]*x[2]+pa[3];
        df[i*n+2] = pa[0]*x[0] + pa[1]*x[1];
        df[i*n+3] = pa[4];
        df[i*n+6] = pa[5];
        //f[i] = pa[0]*x[0]*x[2] +pa[1]*x[1]*x[2] + pa[2]*x[0] + pa[3]*x[1] + pa[4]*x[3] + pa[5]*x[6] + pa[6];
        break;
        }
      case 1:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+0] = pa[7]*x[2]+pa[9];
        df[i*n+1] = pa[8]*x[2]+pa[10];
        df[i*n+2] = pa[7]*x[0]+pa[8]*x[1];
        df[i*n+3] = pa[11];        
        //f[i] = pa[7]*x[0]*x[2]+pa[8]*x[1]*x[2] +pa[9]*x[0] + pa[10]*x[1] + pa[11]*x[3] + pa[12];
        break;
        }
      case 2:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+0] = pa[14];
        df[i*n+1] = pa[15];
        df[i*n+5] = pa[13]*x[7];
        df[i*n+7] = pa[13]*x[5];        
        //f[i] = pa[13]*x[5]*x[7] + pa[14]*x[0] + pa[15]*x[1];
        break;
        }
      case 3:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+0] = pa[16];
        df[i*n+1] = pa[17]; 
        //f[i] = pa[16]*x[0] + pa[17]*x[1] + pa[18];
        break;
        }
      case 4:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+0] = 2.0*x[0];
        df[i*n+1] = 2.0*x[1]; 
        //f[i] = x[0]*x[0] + x[1]*x[1] - 1.0;
        break;
        }
      case 5:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+2] = 2.0*x[2];
        df[i*n+3] = 2.0*x[3]; 
        //f[i] = x[2]*x[2] + x[3]*x[3] - 1.0;
        break;
        }
      case 6:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+4] = 2.0*x[4];
        df[i*n+5] = 2.0*x[5]; 
        //f[i] = x[4]*x[4]+x[5]*x[5] - 1.0;
        break;
        }
      case 7:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+6] = 2.0*x[6];
        df[i*n+7] = 2.0*x[7]; 
        //f[i] = x[6]*x[6]+x[7]*x[7] - 1.0;
        break;
        }
      default:
        break;
      }
    }
  }

void EquationSet_Robot_kinematics::jac(intBox *pBox, long idxEq, long idxVr)
  {
  return jacDigit(pBox, idxEq, idxVr);

  Interval *X = pBox->X;
  Interval *df= pBox->df;

  IntervalD *x = pBox->vx;
  IntervalD *f = pBox->vf;

  long i, n = pBox->szX, is=0, ie = n;
  long j, js = 0, je = n;

  for(i=0; i<n; ++i)
    {
    x[i].setValue(X[i]);
    x[i].setIdx_Var(i);
    }

  if(idxEq>=0)
    {
    is = idxEq; ie = is+1;
    }
  if(idxVr>=0)
    {
    js = idxVr; je = js+1;
    }

  for(i=is; i<ie; ++i)
    {
    switch(i)
      {
      case 0:
        {
        f[i] = pa[0]*x[0]*x[2] +pa[1]*x[1]*x[2] + pa[2]*x[0] + pa[3]*x[1] + pa[4]*x[3] + pa[5]*x[6] + pa[6];
        for(j=js; j<je; ++j)
          df[i*n+j] = f[i].getDerivative(j);
        break;
        }
      case 1:  
        f[i] = pa[7]*x[0]*x[2]+pa[8]*x[1]*x[2] +pa[9]*x[0] + pa[10]*x[1] + pa[11]*x[3] + pa[12];
        for(j=js; j<je; ++j)
          df[i*n+j] = f[i].getDerivative(j);
        break;
      case 2:
        f[i] = pa[13]*x[5]*x[7] + pa[14]*x[0] + pa[15]*x[1];
        for(j=js; j<je; ++j)
          df[i*n+j] = f[i].getDerivative(j);
        break;
      case 3:        
        f[i] = pa[16]*x[0] + pa[17]*x[1] + pa[18];
        for(j=js; j<je; ++j)
          df[i*n+j] = f[i].getDerivative(j);
        break;
      case 4:
        f[i] = x[0]*x[0] + x[1]*x[1] - 1.0;
        for(j=js; j<je; ++j)
          df[i*n+j] = f[i].getDerivative(j);
        break;
      case 5:
        f[i] = x[2]*x[2] + x[3]*x[3] - 1.0;
        for(j=js; j<je; ++j)
          df[i*n+j] = f[i].getDerivative(j);
        break;
      case 6:
        f[i] = x[4]*x[4]+x[5]*x[5] - 1.0;
        for(j=js; j<je; ++j)
          df[i*n+j] = f[i].getDerivative(j);
        break;
      case 7:
        f[i] = x[6]*x[6]+x[7]*x[7] - 1.0;
        for(j=js; j<je; ++j)
          df[i*n+j] = f[i].getDerivative(j);
        break;
      default:
        break;
      }
    }
  }
