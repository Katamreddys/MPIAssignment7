#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <mpi.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

float f1(float x, int intensity);
float f2(float x, int intensity);
float f3(float x, int intensity);
float f4(float x, int intensity);

#ifdef __cplusplus
}
#endif

  
int main (int argc, char* argv[]) {
  
  if (argc < 6) {
    std::cerr<<"usage: "<<argv[0]<<" <functionid> <a> <b> <n> <intensity>"<<std::endl;
    return -1;
  }
  
  MPI_Init(&argc,&argv);
  int fid = atoi(argv[1]);
  float a = atof(argv[2]);
  float b = atof(argv[3]);
  int n = atoi(argv[4]);
  int intensity = atoi(argv[5]);
  int rank,p;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);
  float (*func)(float,int);
  float limit = (b-a)/n;
  int itstart,itend,partition;
  double start = MPI_Wtime(); 
  partition = n/p;
  itstart = rank*partition; itend = (rank+1)*partition;
  if(rank == p-1)
	  itend = n;
  switch(fid)
  {
     case 1 : func = &f1;
                
                break;

     case 2 : func = &f2;
                break;

     case 3 : func = &f3;
                break;

     case 4 : func = &f4;
                break;

     default :  cerr<<"Functionid does not exists, Enter 1, 2, 3 or 4" <<endl; 
                return -1;                

   };
   float par_integralval= 0.0,integralval = 0.0;
   for(int i = itstart ; i<itend ;i++)
        {
           float x = a + ( (i + 0.5) *limit );
	       par_integralval += (*func)(x,intensity); 
        }
   MPI_Reduce(&par_integralval,&integralval,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
   if(rank == 0)
   {
	   integralval = integralval * limit;
	   cout<<integralval<<endl;
           double end = MPI_Wtime(); 
           cerr<<end-start<<endl;
   }
   
  
   
   MPI_Finalize();
   
  return 0;
}