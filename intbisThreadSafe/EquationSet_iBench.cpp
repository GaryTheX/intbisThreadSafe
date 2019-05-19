#include "EquationSet_iBench.h"


EquationSet_iBench::EquationSet_iBench(void)
  {
  pa = 0;
  pb = 0;
  }


EquationSet_iBench::~EquationSet_iBench(void)
  {
  delete[] pa;
  delete[] pb;
  }

char* EquationSet_iBench::getName()
  {
  if(icase==1)
    return "Benchmark_1";
  else if(icase==2)
    return "Benchmark_2";
  else if(icase==3)
    return "Benchmark_3";
  else if(icase==4)
    return "Benchmark_4";
  else
    return "";
  }

void EquationSet_iBench::setBenchMarkID(int id)
  {
  icase = id;
  switch(icase)
    {
    case 1:
    case 4:
      {
      neq = 10;
      pa = new double[neq];
      pb = new double[neq];
      pa[0] = -0.25428722;    pb[0] = -0.18324757;
      pa[1] = -0.37842197;    pb[1] = -0.16275449;
      pa[2] = -0.27162577;    pb[2] = -0.16955071;
      pa[3] = -0.19807914;    pb[3] = -0.15585316;
      pa[4] = -0.44166728;    pb[4] = -0.19950920;
      pa[5] = -0.14654113;    pb[5] = -0.18922793;
      pa[6] = -0.42937161;    pb[6] = -0.21180486;
      pa[7] = -0.07056438;    pb[7] = -0.17081208;
      pa[8] = -0.34504906;    pb[8] = -0.19612740;
      pa[9] = -0.42651102;    pb[9] = -0.21466544;
      break;
      }
    case 2:
    case 3:
      {
      neq = 20;
      pa = new double[neq];
      pb = new double[neq];

      pa[0] = -0.24863995;    pb[0] = -0.19594124;
      pa[1] = -0.87528587;    pb[1] = -0.05612619;
      pa[2] = -0.23939835;    pb[2] = -0.20177810;
      pa[3] = -0.47620128;    pb[3] = -0.16497518;
      pa[4] = -0.24711044;    pb[4] = -0.20198178;
      pa[5] = -0.33565227;    pb[5] = -0.15724045;
      pa[6] = -0.13128974;    pb[6] = -0.12384342;
      pa[7] = -0.45937304;    pb[7] = -0.18180253;
      pa[8] = -0.46896600;    pb[8] = -0.21241045;
      pa[9] = -0.57596835;    pb[9] = -0.16522613;

      pa[10] = -0.56896263;    pb[10] = -0.17221383;
      pa[11] = -0.70561396;    pb[11] = -0.23556251;
      pa[12] = -0.59642512;    pb[12] = -0.24475135;
      pa[13] = -0.46588640;    pb[13] = -0.21790395;
      pa[14] = -0.10607114;    pb[14] = -0.20920602;
      pa[15] = -0.26516898;    pb[15] = -0.21037773;
      pa[16] = -0.20436664;    pb[16] = -0.19838792;
      pa[17] = -0.56003141;    pb[17] = -0.18114505;
      pa[18] = -0.92894617;    pb[18] = -0.04417537;
      pa[19] = -0.57001682;    pb[19] = -0.17949149;
      break;
      }
    default:
      break;
    }
  }

void EquationSet_iBench::setParameters(long np, double *p)
  {
  ;
  }

bool EquationSet_iBench::obt_buildInBound(intBox *pBox)
  {
  pBox->initialize(neq);
  int i;
  if(icase==1 || icase==3)
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(-2.0, 2.0);
    }
  else if(icase==2)
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(-1.0, 2.0);
    }
  else if(icase==4)
    {
    for(i=0; i<neq; ++i)
      pBox->X[i].setInterval(-1.0, 1.0);
    }

  return true;
  }

void EquationSet_iBench::fun(intBox *pBox, long idx)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;

  long i, n = pBox->szX, is=0, ie = n;

  if(idx>=0)
    {
    is = idx; ie = is+1;
    }

  if(icase==1)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          f[i] = x[0] + pa[0] + pb[0]*x[3]*x[2]*x[8];
          break;
        case 1:
          f[i] = x[1] + pa[1] + pb[1]*x[0]*x[9]*x[5];
          break;
        case 2:
          f[i] = x[2] + pa[2] + pb[2]*x[0]*x[1]*x[9];
          break;
        case 3:
          f[i] = x[3] + pa[3] + pb[3]*x[6]*x[0]*x[5];
          break;
        case 4:
          f[i] = x[4] + pa[4] + pb[4]*x[6]*x[5]*x[2];
          break;
        case 5:
          f[i] = x[5] + pa[5] + pb[5]*x[7]*x[4]*x[9];
          break;
        case 6:
          f[i] = x[6] + pa[6] + pb[6]*x[1]*x[4]*x[7];
          break;
        case 7:
          f[i] = x[7] + pa[7] + pb[7]*x[0]*x[6]*x[5];
          break;
        case 8:
          f[i] = x[8] + pa[8] + pb[8]*x[9]*x[5]*x[7];
          break;
        case 9:
          f[i] = x[9] + pa[9] + pb[9]*x[3]*x[7]*x[0];
          break;
        default:
          break;
        }
      }
    }
  else if(icase==2 || icase==3)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          f[i] = x[0] + pa[0] + pb[0]*x[6]*x[9]*x[15];
          break;
        case 1:
          f[i] = x[1] + pa[1] + pb[1]*x[17]*x[7]*x[10];
          break;
        case 2:
          f[i] = x[2] + pa[2] + pb[2]*x[9]*x[6]*x[10];
          break;
        case 3:
          f[i] = x[3] + pa[3] + pb[3]*x[11]*x[14]*x[0];
          break;
        case 4:
          f[i] = x[4] + pa[4] + pb[4]*x[7]*x[8]*x[15];
          break;
        case 5:
          f[i] = x[5] + pa[5] + pb[5]*x[15]*x[17]*x[10];
          break;
        case 6:
          f[i] = x[6] + pa[6] + pb[6]*x[11]*x[12]*x[14];
          break;
        case 7:
          f[i] = x[7] + pa[7] + pb[7]*x[18]*x[14]*x[17];
          break;
        case 8:
          f[i] = x[8] + pa[8] + pb[8]*x[12]*x[1]*x[16];
          break;
        case 9:
          f[i] = x[9] + pa[9] + pb[9]*x[11]*x[8]*x[12];
          break;
        case 10:
          f[i] = x[10] + pa[9] + pb[9]*x[15]*x[16]*x[7];
          break;
        case 11:
          f[i] = x[11] + pa[9] + pb[9]*x[13]*x[10]*x[3];
          break;
        case 12:
          f[i] = x[12] + pa[9] + pb[9]*x[6]*x[15]*x[19];
          break;
        case 13:
          f[i] = x[13] + pa[9] + pb[9]*x[12]*x[2]*x[9];
          break;
        case 14:
          f[i] = x[14] + pa[9] + pb[9]*x[0]*x[8]*x[9];
          break;
        case 15:
          f[i] = x[15] + pa[9] + pb[9]*x[3]*x[18]*x[8];
          break;
        case 16:
          f[i] = x[16] + pa[9] + pb[9]*x[19]*x[9]*x[13];
          break;
        case 17:
          f[i] = x[17] + pa[9] + pb[9]*x[5]*x[12]*x[7];
          break;
        case 18:
          f[i] = x[18] + pa[9] + pb[9]*x[6]*x[12]*x[15];
          break;
        case 19:
          f[i] = x[19] + pa[9] + pb[9]*x[0]*x[2]*x[10];
          break;

        default:
          break;
        }
      }
    }  
  else if(icase==4)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          f[i] = x[0]*x[0] + pa[0] + pb[0]*power(x[3]*x[2]*x[8],2);
          break;
        case 1:
          f[i] = x[1]*x[1] + pa[1] + pb[1]*power(x[0]*x[9]*x[5],2);
          break;
        case 2:
          f[i] = x[2]*x[2] + pa[2] + pb[2]*power(x[0]*x[1]*x[9],2);
          break;
        case 3:
          f[i] = x[3]*x[3] + pa[3] + pb[3]*power(x[6]*x[0]*x[5],2);
          break;
        case 4:
          f[i] = x[4]*x[4] + pa[4] + pb[4]*power(x[6]*x[5]*x[2],2);
          break;
        case 5:
          f[i] = x[5]*x[5] + pa[5] + pb[5]*power(x[7]*x[4]*x[9],2);
          break;
        case 6:
          f[i] = x[6]*x[6] + pa[6] + pb[6]*power(x[1]*x[4]*x[7],2);
          break;
        case 7:
          f[i] = x[7]*x[7] + pa[7] + pb[7]*power(x[0]*x[6]*x[5],2);
          break;
        case 8:
          f[i] = x[8]*x[8] + pa[8] + pb[8]*power(x[9]*x[5]*x[7],2);
          break;
        case 9:
          f[i] = x[9]*x[9] + pa[9] + pb[9]*power(x[3]*x[7]*x[0],2);
          break;
        default:
          break;
        }
      }
    }
  }

void EquationSet_iBench::jacDigit(intBox *pBox, long idxEq, long idxVr)
  {
  Interval *x = pBox->X;
  Interval *f = pBox->f;
  Interval *df= pBox->df;

  long i, j, n = pBox->szX, is=0, ie = n, js=0, je = n;

  if(idxEq>=0)
    {
    is = idxEq; ie = is+1;
    }
  if(idxVr>=0)
    {
    js = idxVr; je = js+1;
    }

  if(icase==1)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = 1.0;
          df[i*n+2] = pb[0]*x[3]*x[8];
          df[i*n+3] = pb[0]*x[2]*x[8];
          df[i*n+8] = pb[0]*x[3]*x[2];
          //f[i] = x[0] + pa[0] + pb[0]*x[3]*x[2]*x[8];
          
          break;
          }
        case 1:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[1]*x[9]*x[5];
          df[i*n+1] = 1.0;
          df[i*n+5] = pb[1]*x[0]*x[9];
          df[i*n+9] = pb[1]*x[0]*x[5];
          //f[i] = x[1] + pa[1] + pb[1]*x[0]*x[9]*x[5];
          break;
          }
        case 2:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[2]*x[1]*x[9];
          df[i*n+1] = pb[2]*x[0]*x[9];
          df[i*n+2] = 1.0;
          df[i*n+9] = pb[2]*x[0]*x[1];
          //f[i] = x[2] + pa[2] + pb[2]*x[0]*x[1]*x[9];
          break;
          }
        case 3:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[3]*x[6]*x[5];
          df[i*n+3] = 1.0;
          df[i*n+5] = pb[3]*x[6]*x[0];
          df[i*n+6] = pb[3]*x[0]*x[5];
          //f[i] = x[3] + pa[3] + pb[3]*x[6]*x[0]*x[5];
          break;
          }
        case 4:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+2] = pb[4]*x[6]*x[5];
          df[i*n+4] = 1.0;
          df[i*n+5] = pb[4]*x[6]*x[2];
          df[i*n+6] = pb[4]*x[5]*x[2];
          //f[i] = x[4] + pa[4] + pb[4]*x[6]*x[5]*x[2];
          break;
          }
        case 5:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+4] = pb[5]*x[7]*x[9];
          df[i*n+5] = 1.0;
          df[i*n+7] = pb[5]*x[4]*x[9];
          df[i*n+9] = pb[5]*x[7]*x[4];
          //f[i] = x[5] + pa[5] + pb[5]*x[7]*x[4]*x[9];
          break;
          }
        case 6:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+1] = pb[6]*x[4]*x[7];
          df[i*n+4] = pb[6]*x[1]*x[7];
          df[i*n+7] = pb[6]*x[1]*x[4];
          df[i*n+6] = 1.0;
          break;
          }
          //f[i] = x[6] + pa[6] + pb[6]*x[1]*x[4]*x[7];
        case 7:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[7]*x[6]*x[5];
          df[i*n+5] = pb[7]*x[0]*x[6];
          df[i*n+6] = pb[7]*x[0]*x[5];
          df[i*n+7] = 1.0;
          //f[i] = x[7] + pa[7] + pb[7]*x[0]*x[6]*x[5];
          break;
          }
        case 8:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+5] = pb[8]*x[9]*x[7];
          df[i*n+7] = pb[8]*x[9]*x[5];
          df[i*n+8] = 1.0;
          df[i*n+9] = pb[8]*x[5]*x[7];
          //f[i] = x[8] + pa[8] + pb[8]*x[9]*x[5]*x[7];
          break;
          }
        case 9:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[9]*x[3]*x[7];
          df[i*n+3] = pb[9]*x[7]*x[0];
          df[i*n+7] = pb[9]*x[3]*x[0];
          df[i*n+9] = 1.0;
          //f[i] = x[9] + pa[9] + pb[9]*x[3]*x[7]*x[0];
          break;
          }
        default:
          break;
        }
      }
    }
  else if(icase==2 || icase==3)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = 1.0;
          df[i*n+6] = pb[0]*x[9]*x[15];
          df[i*n+9] = pb[0]*x[6]*x[15];
          df[i*n+15] = pb[0]*x[6]*x[9];
          //f[i] = x[0] + pa[0] + pb[0]*x[6]*x[9]*x[15];
          break;
          }
        case 1:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+1] = 1.0;
          df[i*n+7] = pb[1]*x[17]*x[10];
          df[i*n+10] = pb[1]*x[17]*x[7];
          df[i*n+17] = pb[1]*x[7]*x[10];
          //f[i] = x[1] + pa[1] + pb[1]*x[17]*x[7]*x[10];
          break;
          }
        case 2:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+2] = 1.0;
          df[i*n+6] = pb[2]*x[9]*x[10];
          df[i*n+9] = pb[2]*x[6]*x[10];
          df[i*n+10] = pb[2]*x[9]*x[6];
          //f[i] = x[2] + pa[2] + pb[2]*x[9]*x[6]*x[10];
          break;
          }
        case 3:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[3]*x[11]*x[14];
          df[i*n+3] = 1.0;
          df[i*n+11] = pb[3]*x[14]*x[0];
          df[i*n+14] = pb[3]*x[11]*x[0];
          //f[i] = x[3] + pa[3] + pb[3]*x[11]*x[14]*x[0];
          break;
          }
        case 4:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+4] = 1.0;
          df[i*n+7] = pb[4]*x[8]*x[15];
          df[i*n+8] = pb[4]*x[7]*x[15];
          df[i*n+15] = pb[4]*x[7]*x[8];
          //f[i] = x[4] + pa[4] + pb[4]*x[7]*x[8]*x[15];
          break;
          }
        case 5:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+5] = 1.0;
          df[i*n+10] = pb[5]*x[15]*x[17];
          df[i*n+15] = pb[5]*x[17]*x[10];
          df[i*n+17] = pb[5]*x[15]*x[10];
          //f[i] = x[5] + pa[5] + pb[5]*x[15]*x[17]*x[10];
          break;
          }
        case 6:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+6] = 1.0;
          df[i*n+11] = pb[6]*x[12]*x[14];
          df[i*n+12] = pb[6]*x[11]*x[14];
          df[i*n+14] = pb[6]*x[11]*x[12];
          //f[i] = x[6] + pa[6] + pb[6]*x[11]*x[12]*x[14];
          break;
          }
        case 7:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+7] = 1.0;
          df[i*n+14] = pb[7]*x[18]*x[17];
          df[i*n+17] = pb[7]*x[18]*x[14];
          df[i*n+18] = pb[7]*x[14]*x[17];
          //f[i] = x[7] + pa[7] + pb[7]*x[18]*x[14]*x[17];
          break;
          }
        case 8:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+1] = pb[8]*x[12]*x[16];
          df[i*n+8] = 1.0;
          df[i*n+12] = pb[8]*x[1]*x[16];
          df[i*n+16] = pb[8]*x[12]*x[1];
          //f[i] = x[8] + pa[8] + pb[8]*x[12]*x[1]*x[16];
          break;
          }
        case 9:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+8] = pb[9]*x[11]*x[12];
          df[i*n+9] = 1.0;
          df[i*n+11] = pb[9]*x[8]*x[12];
          df[i*n+12] = pb[9]*x[11]*x[8];
          //f[i] = x[9] + pa[9] + pb[9]*x[11]*x[8]*x[12];
          break;
          }
        case 10:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+10] = 1.0;
          df[i*n+15] = pb[9]*x[16]*x[7];
          df[i*n+16] = pb[9]*x[15]*x[7];
          df[i*n+7] = pb[9]*x[15]*x[16];
          //f[i] = x[10] + pa[9] + pb[9]*x[15]*x[16]*x[7];
          break;
          }
        case 11:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+3] = pb[9]*x[13]*x[10];
          df[i*n+10] = pb[9]*x[13]*x[3];
          df[i*n+11] = 1.0;
          df[i*n+13] = pb[9]*x[10]*x[3];
          //f[i] = x[11] + pa[9] + pb[9]*x[13]*x[10]*x[3];
          break;
          }
        case 12:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+6] = pb[9]*x[15]*x[19];
          df[i*n+12] = 1.0;
          df[i*n+15] = pb[9]*x[6]*x[19];
          df[i*n+19] = pb[9]*x[6]*x[15];
          //f[i] = x[12] + pa[9] + pb[9]*x[6]*x[15]*x[19];
          break;
          }
        case 13:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+2] = pb[9]*x[12]*x[9];
          df[i*n+9] = pb[9]*x[12]*x[2];
          df[i*n+12] = pb[9]*x[2]*x[9];
          df[i*n+13] = 1.0;
          //f[i] = x[13] + pa[9] + pb[9]*x[12]*x[2]*x[9];
          break;
          }
        case 14:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[9]*x[8]*x[9];
          df[i*n+8] = pb[9]*x[0]*x[9];
          df[i*n+9] = pb[9]*x[0]*x[8];
          df[i*n+14] = 1.0;
          //f[i] = x[14] + pa[9] + pb[9]*x[0]*x[8]*x[9];
          break;
          }
        case 15:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+3] = pb[9]*x[18]*x[8];
          df[i*n+8] = pb[9]*x[3]*x[18];
          df[i*n+15] = 1.0;
          df[i*n+18] = pb[9]*x[3]*x[8];
          //f[i] = x[15] + pa[9] + pb[9]*x[3]*x[18]*x[8];
          break;
          }
        case 16:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+9] = pb[9]*x[19]*x[13];
          df[i*n+13] = pb[9]*x[19]*x[9];
          df[i*n+16] = 1.0;
          df[i*n+19] = pb[9]*x[9]*x[13];
          //f[i] = x[16] + pa[9] + pb[9]*x[19]*x[9]*x[13];
          break;
          }
        case 17:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+5] = pb[9]*x[12]*x[7];
          df[i*n+7] = pb[9]*x[5]*x[12];
          df[i*n+12] = pb[9]*x[5]*x[7];
          df[i*n+17] = 1.0;
          //f[i] = x[17] + pa[9] + pb[9]*x[5]*x[12]*x[7];
          break;
          }
        case 18:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+6] = pb[9]*x[12]*x[15];
          df[i*n+12] = pb[9]*x[6]*x[15];
          df[i*n+15] = pb[9]*x[6]*x[12];
          df[i*n+18] = 1.0;
          //f[i] = x[18] + pa[9] + pb[9]*x[6]*x[12]*x[15];
          break;
          }
        case 19:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] =  pb[9]*x[2]*x[10];
          df[i*n+2] =  pb[9]*x[0]*x[10];
          df[i*n+10] =  pb[9]*x[0]*x[2];
          df[i*n+19] = 1.0;
          //f[i] = x[19] + pa[9] + pb[9]*x[0]*x[2]*x[10];
          break;
          }
        default:
          break;
        }
      }
    }  
  else if(icase==4)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = 2.0*x[0];
          df[i*n+2] =  pb[0]*power(x[3]*x[8],2)* 2.0*x[2];
          df[i*n+3] =  pb[0]*power(x[2]*x[8],2)* 2.0*x[3];
          df[i*n+8] =  pb[0]*power(x[3]*x[2],2)* 2.0*x[8];
          //f[i] = x[0]*x[0] + pa[0] + pb[0]*power(x[3]*x[2]*x[8],2);
          break;
          }
        case 1:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[1]*power(x[9]*x[5],2)* 2.0*x[0];
          df[i*n+1] = 2.0*x[1];
          df[i*n+5] = pb[1]*power(x[0]*x[9],2)* 2.0*x[5];
          df[i*n+9] = pb[1]*power(x[0]*x[5],2)* 2.0*x[9];
          //f[i] = x[1]*x[1] + pa[1] + pb[1]*power(x[0]*x[9]*x[5],2);
          break;
          }
        case 2:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[2]*power(x[1]*x[9],2)* 2.0*x[0];
          df[i*n+1] = pb[2]*power(x[0]*x[9],2)* 2.0*x[1];
          df[i*n+2] = 2.0*x[2];
          df[i*n+9] = pb[2]*power(x[0]*x[1],2)* 2.0*x[9];
          //f[i] = x[2]*x[2] + pa[2] + pb[2]*power(x[0]*x[1]*x[9],2);
          break;
          }
        case 3:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[3]*power(x[6]*x[5],2)* 2.0*x[0];
          df[i*n+3] = 2.0*x[3];
          df[i*n+5] = pb[3]*power(x[6]*x[0],2)* 2.0*x[5];
          df[i*n+6] = pb[3]*power(x[0]*x[5],2)* 2.0*x[6];
          //f[i] = x[3]*x[3] + pa[3] + pb[3]*power(x[6]*x[0]*x[5],2);
          break;
          }
        case 4:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+2] = pb[4]*power(x[6]*x[5],2)*2.0*x[2];
          df[i*n+4] = 2.0*x[4];
          df[i*n+5] = pb[4]*power(x[6]*x[2],2)*2.0*x[5];
          df[i*n+6] = pb[4]*power(x[5]*x[2],2)*2.0*x[6];
          //f[i] = x[4]*x[4] + pa[4] + pb[4]*power(x[6]*x[5]*x[2],2);
          break;
          }
        case 5:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+4] = pb[5]*power(x[7]*x[9],2)*2.0*x[4];
          df[i*n+5] = 2.0*x[5];
          df[i*n+7] = pb[5]*power(x[4]*x[9],2)*2.0*x[7];
          df[i*n+9] = pb[5]*power(x[7]*x[4],2)*2.0*x[9];
          //f[i] = x[5]*x[5] + pa[5] + pb[5]*power(x[7]*x[4]*x[9],2);
          break;
          }
        case 6:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+1] = pb[6]*power(x[4]*x[7],2)*2.0*x[1];
          df[i*n+4] = pb[6]*power(x[1]*x[7],2)*2.0*x[4];
          df[i*n+6] = 2.0*x[6];
          df[i*n+7] = pb[6]*power(x[1]*x[4],2)*2.0*x[7];
          //f[i] = x[6]*x[6] + pa[6] + pb[6]*power(x[1]*x[4]*x[7],2);
          break;
          }
        case 7:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[7]*power(x[6]*x[5],2)*2.0*x[0];
          df[i*n+5] = pb[7]*power(x[0]*x[6],2)*2.0*x[5];
          df[i*n+6] = pb[7]*power(x[0]*x[5],2)*2.0*x[6];
          df[i*n+7] = 2.0*x[7];
          //f[i] = x[7]*x[7] + pa[7] + pb[7]*power(x[0]*x[6]*x[5],2);
          break;
          }
        case 8:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+5] = pb[8]*power(x[9]*x[7],2)*2.0*x[5];
          df[i*n+7] = pb[8]*power(x[9]*x[5],2)*2.0*x[7];
          df[i*n+8] = 2.0*x[8];
          df[i*n+9] = pb[8]*power(x[5]*x[7],2)*2.0*x[9];
          //f[i] = x[8]*x[8] + pa[8] + pb[8]*power(x[9]*x[5]*x[7],2);
          break;
          }
        case 9:
          {
          for(j=js; j<je; ++j) df[i*n+j] = 0.0;
          df[i*n+0] = pb[9]*power(x[3]*x[7],2)*2.0*x[0];
          df[i*n+3] = pb[9]*power(x[7]*x[0],2)*2.0*x[3];
          df[i*n+7] = pb[9]*power(x[3]*x[0],2)*2.0*x[7];
          df[i*n+9] = 2.0*x[9];
          //f[i] = x[9]*x[9] + pa[9] + pb[9]*power(x[3]*x[7]*x[0],2);
          break;
          }
        default:
          break;
        }
      }
    }
  }

void EquationSet_iBench::jac(intBox *pBox, long idxEq, long idxVr)
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

  if(icase==1)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          {
          f[i] = x[0] + pa[0] + pb[0]*x[3]*x[2]*x[8];
          break;
          }
        case 1:  
          f[i] = x[1] + pa[1] + pb[1]*x[0]*x[9]*x[5];
          break;
        case 2:
          f[i] = x[2] + pa[2] + pb[2]*x[0]*x[1]*x[9];
          break;
        case 3:        
          f[i] = x[3] + pa[3] + pb[3]*x[6]*x[0]*x[5];
          break;
        case 4:
          f[i] = x[4] + pa[4] + pb[4]*x[6]*x[5]*x[2];
          break;
        case 5:
          f[i] = x[5] + pa[5] + pb[5]*x[7]*x[4]*x[9];
          break;
        case 6:
          f[i] = x[6] + pa[6] + pb[6]*x[1]*x[4]*x[7];
          break;
        case 7:
          f[i] = x[7] + pa[7] + pb[7]*x[0]*x[6]*x[5];
          break;
        case 8:
          f[i] = x[8] + pa[8] + pb[8]*x[9]*x[5]*x[7];
          break;
        case 9:
          f[i] = x[9] + pa[9] + pb[9]*x[3]*x[7]*x[0];
          break;
        default:
          break;
        }
      for(j=js; j<je; ++j)
        df[i*n+j] = f[i].getDerivative(j);
      }
    }
  else if(icase==2 || icase==3)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          f[i] = x[0] + pa[0] + pb[0]*x[6]*x[9]*x[15];
          break;
        case 1:
          f[i] = x[1] + pa[1] + pb[1]*x[17]*x[7]*x[10];
          break;
        case 2:
          f[i] = x[2] + pa[2] + pb[2]*x[9]*x[6]*x[10];
          break;
        case 3:
          f[i] = x[3] + pa[3] + pb[3]*x[11]*x[14]*x[0];
          break;
        case 4:
          f[i] = x[4] + pa[4] + pb[4]*x[7]*x[8]*x[15];
          break;
        case 5:
          f[i] = x[5] + pa[5] + pb[5]*x[15]*x[17]*x[10];
          break;
        case 6:
          f[i] = x[6] + pa[6] + pb[6]*x[11]*x[12]*x[14];
          break;
        case 7:
          f[i] = x[7] + pa[7] + pb[7]*x[18]*x[14]*x[17];
          break;
        case 8:
          f[i] = x[8] + pa[8] + pb[8]*x[12]*x[1]*x[16];
          break;
        case 9:
          f[i] = x[9] + pa[9] + pb[9]*x[11]*x[8]*x[12];
          break;
        case 10:
          f[i] = x[10] + pa[9] + pb[9]*x[15]*x[16]*x[7];
          break;
        case 11:
          f[i] = x[11] + pa[9] + pb[9]*x[13]*x[10]*x[3];
          break;
        case 12:
          f[i] = x[12] + pa[9] + pb[9]*x[6]*x[15]*x[19];
          break;
        case 13:
          f[i] = x[13] + pa[9] + pb[9]*x[12]*x[2]*x[9];
          break;
        case 14:
          f[i] = x[14] + pa[9] + pb[9]*x[0]*x[8]*x[9];
          break;
        case 15:
          f[i] = x[15] + pa[9] + pb[9]*x[3]*x[18]*x[8];
          break;
        case 16:
          f[i] = x[16] + pa[9] + pb[9]*x[19]*x[9]*x[13];
          break;
        case 17:
          f[i] = x[17] + pa[9] + pb[9]*x[5]*x[12]*x[7];
          break;
        case 18:
          f[i] = x[18] + pa[9] + pb[9]*x[6]*x[12]*x[15];
          break;
        case 19:
          f[i] = x[19] + pa[9] + pb[9]*x[0]*x[2]*x[10];
          break;

        default:
          break;
        }
      for(j=js; j<je; ++j)
        df[i*n+j] = f[i].getDerivative(j);
      }
    }
  else if(icase==4)
    {
    for(i=is; i<ie; ++i)
      {
      switch(i)
        {
        case 0:
          f[i] = x[0]*x[0] + pa[0] + pb[0]*power(x[3]*x[2]*x[8],2);
          break;
        case 1:
          f[i] = x[1]*x[1] + pa[1] + pb[1]*power(x[0]*x[9]*x[5],2);
          break;
        case 2:
          f[i] = x[2]*x[2] + pa[2] + pb[2]*power(x[0]*x[1]*x[9],2);
          break;
        case 3:
          f[i] = x[3]*x[3] + pa[3] + pb[3]*power(x[6]*x[0]*x[5],2);
          break;
        case 4:
          f[i] = x[4]*x[4] + pa[4] + pb[4]*power(x[6]*x[5]*x[2],2);
          break;
        case 5:
          f[i] = x[5]*x[5] + pa[5] + pb[5]*power(x[7]*x[4]*x[9],2);
          break;
        case 6:
          f[i] = x[6]*x[6] + pa[6] + pb[6]*power(x[1]*x[4]*x[7],2);
          break;
        case 7:
          f[i] = x[7]*x[7] + pa[7] + pb[7]*power(x[0]*x[6]*x[5],2);
          break;
        case 8:
          f[i] = x[8]*x[8] + pa[8] + pb[8]*power(x[9]*x[5]*x[7],2);
          break;
        case 9:
          f[i] = x[9]*x[9] + pa[9] + pb[9]*power(x[3]*x[7]*x[0],2);
          break;
        default:
          break;
        }

      for(j=js; j<je; ++j)
        df[i*n+j] = f[i].getDerivative(j);
      }
    }
  }
