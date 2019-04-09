#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <mpi.h>
#include <math.h>

using namespace std;

float genA (int row, int col) {
  if (row > col)
    return 1.0;
  else
    return 0.0;
}

float genx0 (int i) {
  return 1.0;
}


void checkx (int iter, long i, float xval) {
  if (iter == 1) {
    float shouldbe = i;
    if (fabs(xval/shouldbe) > 1.01 || fabs(xval/shouldbe) < .99 )
      cout<<"incorrect : x["<<i<<"] at iteration "<<iter<<" should be "<<shouldbe<<" not "<<xval<<endl;
  }

  if (iter == 2) {
    float shouldbe =(i-1)*i/2;
    if (fabs(xval/shouldbe) > 1.01 || fabs(xval/shouldbe) < .99)
      cout<<"incorrect : x["<<i<<"] at iteration "<<iter<<" should be "<<shouldbe<<" not "<<xval<<endl;
  }
}

//perform dense y=Ax on an n \times n matrix
void matmul(float*A, float*x, float*y, long n) {
  
  for (long row = 0; row<n; ++row) {
    float sum = 0;
    
    for (long col = 0; col<n ; ++col) {
      //sum += x[col] *A[row][col]
       sum += x[col] * A[row*n+col];
    }

       y[row] = sum;
  }
}

int main (int argc, char* argv[]) {

  if (argc < 3) {
    std::cout<<"usage: "<<argv[0]<<" <n> <iteration>"<<std::endl;
  }
  //initialize data
  
   MPI_Init(&argc,&argv);

   long n = atol(argv[1]);

   int iteration = atoi(argv[2]);
   bool check = true;

   int rankmpi,sizenpi;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&rankmpi);
   MPI_Comm_size(MPI_COMM_WORLD,&sizenpi);

   int p = sqrt(sizenpi);
   long divInd = n/p; 
   
   int divRow = rankmpi/p,divCol = rankmpi%p;
   
   MPI_Comm commRow;
   int rankRow;
   MPI_Comm_split(MPI_COMM_WORLD, divRow, rankmpi, &commRow);
   MPI_Comm_rank(commRow,&rankRow);

   MPI_Comm commCol;
   int rankCol;
   MPI_Comm_split(MPI_COMM_WORLD, divCol, rankmpi, &commCol);
   MPI_Comm_rank(commCol,&rankCol);

   float* mainArray = new float[divInd*divInd];
   long beginRow = (divRow*divInd),startCol = (divCol*divInd);
   long endRow = beginRow+divInd,colend = startCol+divInd;

  for (long row = beginRow,rset=0; row<endRow; row++,rset++) {
    for (long col= startCol,cset=0; col<colend; col++,cset++) {
      mainArray[(rset*divInd)+cset] = genA(row, col);
    }
  }

  float* x = new float[divInd];

  for (long i=0; i<divInd; i++)
    {
     x[i] = genx0(i);
    }
  float* y = new float[divInd];
  for (long i=0; i<divInd; i++)
    y[i] = 0.0;
 
   double beginTime = MPI_Wtime(); 
   for (int it = 0; it<iteration; it++) 
   {
      matmul(mainArray,x,y,divInd);
      MPI_Reduce(y,x,divInd,MPI_FLOAT,MPI_SUM,divRow,commRow);
      MPI_Bcast(x,divInd,MPI_FLOAT,divCol, commCol);
  if (check){
    for (long i = startCol,p=0; i<colend; ++i,++p){
        checkx (it+1, i, x[p]);
             }
          }
  }
  
  if(rankmpi == 0){
        double endTime = MPI_Wtime(); 
        cerr<<endTime-beginTime<<endl;
   }
 
  MPI_Comm_free(&commRow);
  MPI_Comm_free(&commCol);

  delete[] mainArray;
  delete[] x;
  delete[] y;

  MPI_Finalize();
 
  return 0;
}