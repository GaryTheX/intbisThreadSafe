#include "EquationSet_Neurophysiology.h"


EquationSet_Neurophysiology::EquationSet_Neurophysiology(void)
  {
  neq = 6;
  pc = new double[4];
  pc[0] = pc[2] = pc[1] = pc[3] = 0.0;
  icase = 1;
  }


EquationSet_Neurophysiology::~EquationSet_Neurophysiology(void)
  {
  delete[] pc;
  }

char* EquationSet_Neurophysiology::getName()
{
  if (icase == 1)
    return "Neurophysiology_1";
  else if (icase == 2)
    return "Neurophysiology_2";
  else if (icase == 3)
    return "Neurophysiology_3";
  else
    return "";
}

bool EquationSet_Neurophysiology::obt_buildInBound(intBox *pBox)
  {
  pBox->initialize(neq);
  int i;
  for(i=0; i<neq; ++i)
    pBox->X[i].setInterval(-1.0e0, 1.0e0);
  pBox->X[4].setInterval(-1.0,1.0);
  pBox->X[5].setInterval(-1.0,1.0);

  return true;
  }

void EquationSet_Neurophysiology::fun(intBox *pBox, long idx)
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
        f[i] = power(x[0],2) + power(x[2],2) - 1.0;
        break;
      case 1:
        f[i] = power(x[1],2)+ power(x[3],2) - 1.0;
        break;
      case 2:
        f[i] = x[4]*power(x[2],3) + x[5]*power(x[3],3) - pc[0];
        break;
      case 3:
        f[i] = x[4]*power(x[0],3) + x[5]*power(x[1],3) - pc[1];
        break;
      case 4:
        f[i] = x[4]*x[0]*power(x[2],2) + x[5]*power(x[3],2)*x[1] - pc[2];
        break;
      case 5:
        f[i] = x[4]*power(x[0],2)*x[2] + x[5]*power(x[1],2)*x[3] - pc[3];
        break;
      default:
        break;
      }
    }
  }

void EquationSet_Neurophysiology::jac(intBox *pBox, long idxEq, long idxVr)
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
        f[i] = power(x[0],2) + power(x[2],2)- 1.0;
        break;
      case 1:
        f[i] = power(x[1],2) + power(x[3],2) - 1.0;
        break;
      case 2:
        f[i] = x[4]*power(x[2],3) + x[5]*power(x[3],3) - pc[0];
        break;
      case 3:
        f[i] = x[4]*power(x[0],3) + x[5]*power(x[1],3) - pc[1];
        break;
      case 4:
        f[i] = x[4]*x[0]*power(x[2],2) + x[5]*power(x[3],2)*x[1] - pc[2];
        break;
      case 5:
        f[i] = x[4]*power(x[0],2)*x[2] + x[5]*power(x[1],2)*x[3] - pc[3];
        break;
      default:
        break;
      }
    for(j=js; j<je; ++j)
      df[i*n+j] = f[i].getDerivative(j);
    }
  }

void EquationSet_Neurophysiology::jacDigit(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Interval *df= pBox->df;

  long i, n = pBox->szX, is=0, ie = n;
  long j, js = 0, je = n;

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
        df[i*n]   = 2.0*x[0];
        df[i*n+2] = 2.0*x[2];
        //f[i] = power(x[0],2) + power(x[2],2) - 1.0;
        break;
        }
      case 1:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+1] = 2.0*x[1];
        df[i*n+3] = 2.0*x[3];
        //f[i] = power(x[1],2)+ power(x[3],2) - 1.0;
        break;
        }
      case 2:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n+2] = 3.0*power(x[2],2)*x[4];
        df[i*n+3] = 3.0*power(x[3],2)*x[5];
        df[i*n+4] = power(x[2],3);
        df[i*n+5] = power(x[3],3);
        //f[i] = x[4]*power(x[2],3) + x[5]*power(x[3],3) - pc[0];
        break;
        }
      case 3:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n] = x[4]*3.0*power(x[0],2);
        df[i*n+1] = x[5]*3.0*power(x[1],2);
        df[i*n+4] = power(x[0],3);
        df[i*n+5] = power(x[1],3);
        //f[i] = x[4]*power(x[0],3) + x[5]*power(x[1],3) - pc[1];
        break;
        }
      case 4:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n]   = x[4]*power(x[2],2);
        df[i*n+1] = x[5]*power(x[3],2);
        df[i*n+2] = x[4]*x[0]*2.0*x[2];
        df[i*n+3] = x[5]*2.0*x[3]*x[1];
        df[i*n+4] = x[0]*power(x[2],2);
        df[i*n+5] = power(x[3],2)*x[1];
        //f[i] = x[4]*x[0]*power(x[2],2) + x[5]*power(x[3],2)*x[1] - pc[2];
        break;
        }
      case 5:
        {
        for(j=js; j<je; ++j) df[i*n+j] = 0.0;
        df[i*n]   = x[4]*2.0*x[0]*x[2];
        df[i*n+1] = x[5]*2.0*x[1]*x[3];
        df[i*n+2] = x[4]*power(x[0],2);
        df[i*n+3] = x[5]*power(x[1],2);
        df[i*n+4] = power(x[0],2)*x[2];
        df[i*n+5] = power(x[1],2)*x[3];

        //f[i] = x[4]*power(x[0],2)*x[2] + x[5]*power(x[1],2)*x[3] - pc[3];
        break;
        }
      default:
        break;
      }
    }
  }


void EquationSet_Neurophysiology::tens(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Taylor *ts = pBox->ts;

    switch(idxEq)
      {
      case 0:
        {        
        if(idxVr==0)
          {
          ts->Taylor_fatParam = 2.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }
        else if(idxVr==1)
          {
          ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }
        else if(idxVr==2)
          {
          ts->Taylor_fatParam = 2.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }
        else 
          {
          ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }

        break;
        }
      case 1:        
        {
        if(idxVr==0)
          {
          ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }
        else if(idxVr==1)
          {
          ts->Taylor_fatParam = 2.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }
        else if(idxVr==2)
          {
          ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }
        else if(idxVr==3)
          {
          ts->Taylor_fatParam = 2.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }
        else 
          {
          ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
          }
        break;
        }     
      case 2:
          {
          if(idxVr==0)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==1)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==2)
            {
            ts->Taylor_fatParam = x[4]*6.0*x[2]; ts->Remain_fatParam = x[4]*6.0;  ts->isRemainConst = false;
            }
          else if(idxVr==3)
            {
            ts->Taylor_fatParam = x[5]*6.0*x[3]; ts->Remain_fatParam = x[5]*6.0;  ts->isRemainConst = false;
            }
          else if(idxVr==4)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==5)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          break;   
          }
 
      case 3:        
          {
          if(idxVr==0)
            {
            ts->Taylor_fatParam = x[4]*6.0*x[0]; ts->Remain_fatParam = x[4]*6.0;  ts->isRemainConst = false;
            }
          else if(idxVr==1)
            {
            ts->Taylor_fatParam = x[5]*6.0*x[1]; ts->Remain_fatParam = x[5]*6.0;  ts->isRemainConst = false;
            }
          else if(idxVr==2)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==3)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==4)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==5)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          break;        
          }
      case 4:
          {
          if(idxVr==0)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==1)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==2)
            {
            ts->Taylor_fatParam = x[4]*x[0]*2.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==3)
            {
            ts->Taylor_fatParam = x[5]*x[1]*2.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else 
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          break;
          }
      case 5:
          {
          if(idxVr==0)
            {
            ts->Taylor_fatParam = x[4]*2.0*x[2]; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==1)
            {
            ts->Taylor_fatParam = x[5]*2.0*x[3]; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==2)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==3)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==4)
            {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          else if(idxVr==5)
             {
            ts->Taylor_fatParam = 0.0; ts->Remain_fatParam = 0.0;  ts->isRemainConst = true;
            }
          break;
          }

      default:
        break;
      }
  }
