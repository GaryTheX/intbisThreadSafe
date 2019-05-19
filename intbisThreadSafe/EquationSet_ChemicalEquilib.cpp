#include <math.h>
#include "EquationSet_ChemicalEquilib.h"


EquationSet_ChemicalEquilib::EquationSet_ChemicalEquilib(void)
  {
  }


EquationSet_ChemicalEquilib::~EquationSet_ChemicalEquilib(void)
  {
  }

char* EquationSet_ChemicalEquilib::getName()
  {
  if(icase==1)
    return "Chemical Eqlib";
  else if(icase==2)
    return "Combustion";  
  else
    return "";
  }

void EquationSet_ChemicalEquilib::setBenchMarkID(int id)
  {
  icase = id;
  switch(icase)
    {
    case 1:
      {
      sqrt40 = sqrt(40.0);
      R = 10.0;
      R5 = 0.193;
      R6 = 0.002597/sqrt40;
      R7 = 0.003448/sqrt40;
      R8 = 0.00001799/40.0;
      R9 = 0.0002155/sqrt40;
      R10= 0.00003846/40.0;
      
      neq = 5;
      break;
      }
    case 2:    
      {
      neq = 10;
      break;
      }
    default:
      break;
    }
  }
bool EquationSet_ChemicalEquilib::obt_buildInBound(intBox *pBox)
  {
  pBox->initialize(neq);
  int i; 
  if(icase==1)
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(0.0, 1.0e+8);
    }
  else
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(-1.0, 1.0);
    }

  return true;
  }

void EquationSet_ChemicalEquilib::fun(intBox *pBox, long idx)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;

  long i, n = pBox->szX, is=0, ie = n;

  if(idx>=0)
    {
    is = idx; ie = is+1;
    }
  switch(icase)
      {
      case 1:
        {
        for(i=is; i<ie; ++i)
          {
          switch(i)
            {
            case 0:
              f[i] = x[0]*x[1] + x[0] - 3.0*x[4];
              break;
            case 1:
              f[i] = 2.0*x[0]*x[1]+x[0]+x[1]*power(x[2],2)+R8*x[1]-R*x[4]+2.0*R10*x[1]*x[1]+R7*x[1]*x[2]+R9*x[1]*x[3];
              break;
            case 2:
              f[i] = 2.0*x[1]*power(x[2],2)+2.0*R5*x[2]*x[2] - 8.0*x[4] + R6*x[2] + R7*x[1]*x[2];
              break;
            case 3:
              f[i] = R9*x[1]*x[3] + 2.0*power(x[3],2) - 4.0*R*x[4];
              break;
            case 4:
              f[i] = x[0]*(x[1]+1.0)+R10*x[1]*x[1]+x[1]*x[2]*x[2]+R8*x[1]+R5*x[2]*x[2]+x[3]*x[3]-1.0+R6*x[2]+R7*x[1]*x[2]+R9*x[1]*x[3];
              break;
            default:
              break;
            }
          }
        break;
        }
      case 2:
        {
        f[0] = x[1]+2.0*x[5]+x[8]+2.0*x[9]-1e-5;
        f[1] = x[2]+x[7]-3e-5;
        f[2] = x[0]+x[2]+2.0*x[4]+2.0*x[7]+x[8]+x[9]-5e-5;
        f[3] = x[3]+2.0*x[6]-1e-5;
        f[4] = 0.5140437*1e-7*x[4]-power(x[0],2);
        f[5] = 0.1006932*1e-6*x[5]-2.0*power(x[1],2);
        f[6] = 0.7816278*1e-15*x[6]-power(x[3],2);
        f[7] = 0.1496236*1e-6*x[7]-x[0]*x[2];
        f[8] = 0.6194411*1e-7*x[8]-x[0]*x[1];
        f[9] = 0.2089296*1e-14*x[9]-x[0]*power(x[1],2);
        break;
        }
      default:
        break;
    }
  }

void EquationSet_ChemicalEquilib::jac(intBox *pBox, long idxEq, long idxVr)
  {
  return jacDigit(pBox, idxEq, idxVr);

  Interval *X = pBox->X;
  Interval *F = pBox->f;
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
        f[i] = x[0]*x[1] + x[0] - 3.0*x[4];
        break;
      case 1:
        f[i] = 2.0*x[0]*x[1]+x[0]+x[1]*power(x[2],2)+R8*x[1]-R*x[4]+2.0*R10*x[1]*x[1]+R7*x[1]*x[2]+R9*x[1]*x[3];
        break;
      case 2:
        f[i] = 2.0*x[1]*power(x[2],2)+2.0*R5*x[2]*x[2] - 8.0*x[4] + R6*x[2] + R7*x[1]*x[2];
        break;
      case 3:
        f[i] = R9*x[1]*x[3] + 2.0*power(x[3],2) - 4.0*R*x[4];
        break;
      case 4:
        f[i] = x[0]*(x[1]+1.0)+R10*x[1]*x[1]+x[1]*x[2]*x[2]+R8*x[1]+R5*x[2]*x[2]+x[3]*x[3]-1.0+R6*x[2]+R7*x[1]*x[2]+R9*x[1]*x[3];
        break;
      default:
        break;
      }
    for(j=js; j<je; ++j)
      df[i*n+j] = f[i].getDerivative(j);
    }
  }

void EquationSet_ChemicalEquilib::jacDigit(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Interval *df= pBox->df; 
  long i, n = pBox->szX, is=0, ie = n, neq_1 = n-1;
  long j, js = 0, je = n;

  if(idxEq>=0)
    {
    is = idxEq; ie = is+1;
    }
  if(idxVr>=0)
    {
    js = idxVr; je = js+1;
    }
  switch(icase)
    {
    case 1:
      {
      for(i=is; i<ie; ++i)
        {
        switch(i)
          {
          case 0:
            {
            for(j=js; j<je; ++j) df[i*n+j] = 0.0;
            df[i*n] = x[1]+1.0;
            df[i*n+1] = x[0];
            df[i*n+4] = -3.0;
            //f[i] = x[0]*x[1] + x[0] - 3.0*x[4];
            break;
            }
          case 1:
            {
            for(j=js; j<je; ++j) df[i*n+j] = 0.0;
            df[i*n] = 2.0*x[1] + 1.0;
            df[i*n+1] = 2.0*x[0] + power(x[2],2) + R8 + 4.0*R10*x[1] + R7*x[2]+R9*x[3];
            df[i*n+2] = 2.0*x[1]*x[2] + R7*x[1];
            df[i*n+3] = R9*x[1];
            df[i*n+4] = -R;
            //f[i] = 2.0*x[0]*x[1]+x[0]+x[1]*power(x[2],2)+R8*x[1]-R*x[4]+2.0*R10*x[1]*x[1]+R7*x[1]*x[2]+R9*x[1]*x[3];
            break;
            }
          case 2:
            {
            for(j=js; j<je; ++j) df[i*n+j] = 0.0;
            df[i*n+1] = 2.0*power(x[2],2) + R7*x[2];
            df[i*n+2] = 4.0*x[1]*x[2] + 4.0*R5*x[2] + R6 + R7*x[1];        
            df[i*n+4] = -8.0;
            //f[i] = 2.0*x[1]*power(x[2],2)+2.0*R5*x[2]*x[2] - 8.0*x[4] + R6*x[2] + R7*x[1]*x[2];
            break;
            }
          case 3:
            {
            for(j=js; j<je; ++j) df[i*n+j] = 0.0;
            df[i*n+1] = R9*x[3];
            df[i*n+3] = R9*x[1] + 4.0*x[3];        
            df[i*n+4] = -4.0*R;
            //f[i] = R9*x[1]*x[3] + 2.0*power(x[3],2) - 4.0*R*x[4];
            break;
            }
          case 4:
            {
            for(j=js; j<je; ++j) df[i*n+j] = 0.0;
            df[i*n]   = x[1]+1.0;
            df[i*n+1] = x[0] +2.0*R10*x[1] + power(x[2],2) + R8 + R7*x[2] + R9*x[3];
            df[i*n+2] = 2.0*x[1]*x[2]+R5*2.0*x[2]+R6+R7*x[1];
            df[i*n+3] = 2.0*x[3] + R9*x[1];               
            //f[i] = x[0]*(x[1]+1.0)+R10*x[1]*x[1]+x[1]*x[2]*x[2]+R8*x[1]+R5*x[2]*x[2]+x[3]*x[3]-1.0+R6*x[2]+R7*x[1]*x[2]+R9*x[1]*x[3];
            break;
            }
          default:
            break;
          }
        }
      break;
      }
    case 2:
      {
      for(j=js; j<je; ++j) df[j] = 0.0;
      df[1] = 1.0;
      df[5] = 2.0;
      df[8] = 1.0;
      df[9] = 2.0;

      for(j=js; j<je; ++j) df[n+j] = 0.0;
      df[n+2] = 1.0;
      df[n+7] = 1.0;

      for(j=js; j<je; ++j) df[2*n+j] = 0.0;
      df[2*n+0] = 1.0;
      df[2*n+2] = 1.0;
      df[2*n+4] = 2.0;
      df[2*n+7] = 2.0;
      df[2*n+8] = 1.0;
      df[2*n+9] = 1.0;

      for(j=js; j<je; ++j) df[3*n+j] = 0.0;
      df[3*n+3] = 1.0;
      df[3*n+6] = 2.0;

      for(j=js; j<je; ++j) df[4*n+j] = 0.0;
      df[4*n+0] = -2.0*x[0];
      df[4*n+4] = 0.5140437*1e-7;

      for(j=js; j<je; ++j) df[5*n+j] = 0.0;
      df[5*n+1] = -4.0*x[1];
      df[5*n+5] = 0.1006932*1e-6;

      for(j=js; j<je; ++j) df[6*n+j] = 0.0;
      df[6*n+3] = -2.0*x[3];
      df[6*n+6] = 0.7816278*1e-15;

      for(j=js; j<je; ++j) df[7*n+j] = 0.0;
      df[7*n+0] = -1.0*x[2];
      df[7*n+2] = -1.0*x[0];
      df[7*n+7] = 0.1496236*1e-6;

      for(j=js; j<je; ++j) df[8*n+j] = 0.0;
      df[8*n+0] = -1.0*x[1];
      df[8*n+1] = -1.0*x[0];
      df[8*n+8] = 0.6194411*1e-7;

      for(j=js; j<je; ++j) df[9*n+j] = 0.0;
      df[9*n+0] = -1.0*power(x[1],2);
      df[9*n+1] = -2.0*x[0]*x[1];
      df[9*n+8] = 0.6194411*1e-7;

      break;
      }
    default:
      break;
    }
  }