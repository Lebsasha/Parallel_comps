#include <iostream>
#include <mpich/mpi.h>
int main(int argc, char **argv)
{
	 MPI_Init(&argc, &argv);
  {
    std::cout << "Hello World!\n";
  }
	MPI_Finalize();
}
