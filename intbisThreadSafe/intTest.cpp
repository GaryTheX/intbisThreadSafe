#ifdef _DEBUG
  #include <crtdbg.h>
#endif // _DEBUG

#include <windows.h>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <time.h>

#include "intSolEngine.h"

clock_t start, finish;
double duration;

//#define DOMULTI
//#define INITIALVALUE
long main(long nargs, char *arg1[])
  {
#ifdef _DEBUG
  // Specify memory debugging options
  _CrtSetDbgFlag(
                 _CRTDBG_ALLOC_MEM_DF   // | // Enable debug heap allocations and use of memory block type identifiers, such as _CLIENT_BLOCK.
//               _CRTDBG_CHECK_ALWAYS_DF   | // Call _CrtCheckMemory at every allocation and deallocation request.
//               _CRTDBG_CHECK_CRT_DF      | // Include _CRT_BLOCK types in leak detection and memory state difference operations.
//               _CRTDBG_DELAY_FREE_MEM_DF | // Keep freed memory blocks in the heap's linked list, assign them the _FREE_BLOCK type, and fill them with the byte value 0xDD.
//               _CRTDBG_LEAK_CHECK_DF       // Perform automatic leak checking at program exit.
                 );
#endif

  FILE *freport, *freportSol;
  fopen_s(&freport, "multiReport.txt","w");
  fprintf(freport, "\nProblemName\tprobDimension\tnumThreads\tnumSolution\tCPU(seconds)");
  fopen_s(&freportSol, "multiReportSol.txt", "w");

  intSolEngine *intSolEng = new intSolEngine();

  long iEq, icase, icaseEnd, ipass;

  for(iEq=1; iEq<=6; iEq++)
    {		
    if (iEq == 1)
      icaseEnd = 4;
    else if (iEq == 2)
      icaseEnd = 3;
    else if (iEq == 3)
      icaseEnd = 3;
    else if (iEq == 4)
      icaseEnd = 2;
    else if (iEq == 5)
      icaseEnd = 1;
    else if (iEq == 6)
      continue;
    else if(iEq==7)
      icaseEnd = 2;

    for(icase=1; icase<=icaseEnd; icase++)
      {    
      ipass = 0;
      while(ipass<=2)  // <= 2 means using up to 4 cores
        {
        intSolEng->reinitialize();
        if(ipass==0) intSolEng->setMultiThreadingConfig(1);
        else if(ipass==1) intSolEng->setMultiThreadingConfig(2);
        else if(ipass==2) intSolEng->setMultiThreadingConfig(4);
        else if(ipass==3) intSolEng->setMultiThreadingConfig(8);
        else if(ipass==4) intSolEng->setMultiThreadingConfig(16);

        intSolEng->setEquationSet(iEq, icase);
        if(intSolEng->isReady())
          {
					printf("\nStarted solving...\n");
          start = clock();
          intSolEng->solve();
          finish = clock();
          printf("\nFinished solving...\n");
          duration = (double)(finish - start);
          duration /= CLOCKS_PER_SEC;

          long nsol = intSolEng->solutionsFound();

          if(nsol>0)
            {
        #ifdef INITIALVALUE
            printf("\nInitial value at: %-16.8g\n",initialVal);
        #endif
            fprintf(freportSol,"\n******\t%s\n******\n",intSolEng->getName());

            long neq = intSolEng->getDimension();
            double *pl = new double[neq];
            double *pu = new double[neq];
            for(int i=0; i<nsol; ++i)
              {
              intSolEng->getSolution(pl,pu);
              printf("Solution # %4d:\n",i+1);
              fprintf(freportSol,"Solution # %4d:\n",i+1);
              for(int j=0; j<neq; ++j)
                {
                printf("\t[\t%-16.8g,\t%-16.8g\t]\n",pl[j],pu[j]);
                fprintf(freportSol,"\t[\t%-16.8g,\t%-16.8g\t]\n",pl[j],pu[j]);
                }
              }
            delete[] pl;
            delete[] pu;
            }
    
          printf("\n\nFinished at clock ticks %f\n\n",duration);
          printf("Interval Newton spends %f seconds to solve this problem.\n\t\t------\n",duration);
          printf("\n%s\t%4d\t%2d\t%6d\t%16.8f",intSolEng->getName()
                                                        ,intSolEng->getDimension()
                                                        ,intSolEng->getNumThreads()
                                                        ,nsol, duration);

          fprintf(freportSol,"\n\nFinished at clock ticks %f\n\n",duration);
          fprintf(freportSol,"Interval Newton spends %f seconds to solve this problem.\n\t\t------\n",duration);
          fprintf(freportSol,"\n%s\t%4d\t%2d\t%6d\t%16.8f",intSolEng->getName()
                                                        ,intSolEng->getDimension()
                                                        ,intSolEng->getNumThreads()
                                                        ,nsol, duration);

          fprintf(freport, "\n%s\t%4d\t\t%2d\t\t%6d\t\t%16.8f",intSolEng->getName()
                                                        ,intSolEng->getDimension()
                                                        ,intSolEng->getNumThreads()
                                                        ,nsol, duration);
          }
        printf("\n");
        ipass++;
        }
      fprintf(freport, "\n");
			}
    }

  delete intSolEng;
  fclose(freportSol);
  fclose(freport);

#ifdef _DEBUG
  _CrtDumpMemoryLeaks();
#endif
  return 0;
  }