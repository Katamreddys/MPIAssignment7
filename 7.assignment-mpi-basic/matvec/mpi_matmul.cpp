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

   int iter = atoi(argv[2]);
   bool check = true;

   int worldrank,np;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
   MPI_Comm_size(MPI_COMM_WORLD,&np);

   int p = sqrt(np);
   long each_div = n/p; 
   
   int row_div = worldrank/p,col_div = worldrank%p;
   
   MPI_Comm row_comm;
   int rowrank;
   MPI_Comm_split(MPI_COMM_WORLD, row_div, worldrank, &row_comm);
   MPI_Comm_rank(row_comm,&rowrank);

   MPI_Comm col_comm;
   int colrank;
   MPI_Comm_split(MPI_COMM_WORLD, col_div, worldrank, &col_comm);
   MPI_Comm_rank(col_comm,&colrank);

   float* A = new float[each_div*each_div];
   long rowstart = (row_div*each_div),colstart = (col_div*each_div);
   long rowend = rowstart+each_div,colend = colstart+each_div;

  for (long row = rowstart,rset=0; row<rowend; row++,rset++) {
    for (long col= colstart,cset=0; col<colend; col++,cset++) {
      A[(rset*each_div)+cset] = genA(row, col);
    }
  }

  float* x = new float[each_div];

  for (long i=0; i<each_div; i++)
    {
     x[i] = genx0(i);
    }
  float* y = new float[each_div];
  for (long i=0; i<each_div; i++)
    y[i] = 0.0;
 
   double start = MPI_Wtime(); 
   // cout<<"y---------------"<<worldrank<<endl;
   for (int it = 0; it<iter; it++) 
   {
      matmul(A,x,y,each_div);
      /* cout<<"iteration"<<it<<endl;
       for (long row = 0; row<each_div; row=row+1) {
         float sum = 0.0;
    
          for (long col = 0; col<each_div ; col=col+1) {
          sum = sum + (x[col] * A[(row*each_div)+col]);
          }
             y[row] = sum;
       }*/
      MPI_Reduce(y,x,each_div,MPI_FLOAT,MPI_SUM,row_div,row_comm);
      MPI_Bcast(x,each_div,MPI_FLOAT,col_div, col_comm);
  if (check){
    for (long i = colstart,p=0; i<colend; ++i,++p){
        checkx (it+1, i, x[p]);
             }
          }
  }
  
  if(worldrank == 0){
        double end = MPI_Wtime(); 
        cerr<<end-start<<endl;
   }
 
  MPI_Comm_free(&row_comm);
  MPI_Comm_free(&col_comm);

  delete[] A;
  delete[] x;
  delete[] y;

  MPI_Finalize();
 
  return 0;
}