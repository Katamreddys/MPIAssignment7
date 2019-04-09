#include <mpi.h>
#include <unistd.h>
#include <iostream>

int main(int argc, char*argv[]) {
int MPI_Rank , MPI_Size;
char hostName[1024];
   hostName[1023] = '\0';
   gethostname(hostName, 1023);
MPI_Init(&argc,& argv);
MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);
MPI_Comm_size(MPI_COMM_WORLD, &MPI_Size);
printf(" i am process %d out of %d. i am running on %s. ",MPI_Rank, MPI_Size,hostName);
MPI_Finalize();

}