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
  int function_id = atoi(argv[1]);
  float a = atof(argv[2]);
  float b = atof(argv[3]);
  int n = atoi(argv[4]);
  int intensity = atoi(argv[5]);
  int rank_mpi,size_mpi;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_mpi);
  MPI_Comm_size(MPI_COMM_WORLD,&size_mpi);
  float (*function)(float,int);
  float preCalMul = (b-a)/n;
  int startLoop,tendLoop,pt;
  double startTime = MPI_Wtime(); 
  pt = n/size_mpi;
  startLoop = rank_mpi*pt; tendLoop = (rank_mpi+1)*pt;
  if(rank_mpi == size_mpi-1)
	  tendLoop = n;
  switch(function_id)
  {
     case 1 : function = &f1;
                
                break;

     case 2 : function = &f2;
                break;

     case 3 : function = &f3;
                break;

     case 4 : function = &f4;
                break;

     default :  cerr<<"Please Enter 1, 2, 3 or 4" <<endl; 
                return -1;                

   };
   float partialInt= 0.0,finalInt = 0.0;
   for(int i = startLoop ; i<tendLoop ;i++)
        {
           float x = a + ( (i + 0.5) *preCalMul );
	       partialInt += (*function)(x,intensity); 
        }
   MPI_Reduce(&partialInt,&finalInt,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
   if(rank_mpi == 0)
   {
	   finalInt = finalInt * preCalMul;
	   cout<<finalInt<<endl;
           double endTime = MPI_Wtime(); 
           cerr<<endTime-startTime<<endl;
   }
   
  
   
   MPI_Finalize();
   
  return 0;
}