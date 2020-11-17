#include <iostream>
#ifdef CBB
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#include </usr/include/x86_64-linux-gnu/mpich/mpi.h>
#pragma GCC diagnostic pop
#else // CBB
#include <mpich/mpi.h>
#endif // CBB

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
