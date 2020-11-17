#include <iostream>
#include <mpich/mpi.h>
int main(int argc, char **argv)
{
    int * p = new int;
	 MPI_Init(&argc, &argv);
	 MPI_Comm_size(MPI_COMM_WORLD, p);
  {
    std::cout << "Hello World!\n"<<*p;
  }
	MPI_Finalize();
}
