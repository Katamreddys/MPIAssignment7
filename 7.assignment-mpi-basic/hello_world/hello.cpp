#include <mpi.h>
#include <unistd.h>
#include <iostream>

int main(int argc, char*argv[]) {
int MPI_Rank , MPI_Size;
char Host_name[1024];
   Host_name[1023] = '\0';
   gethostname(Host_name, 1023);
MPI_Init(&argc,& argv);
MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);
MPI_Comm_size(MPI_COMM_WORLD, &MPI_Size);
printf(" i am process %d out of %d. i am running on %s. ",MPI_Rank, MPI_Size,Host_name);
MPI_Finalize();

}