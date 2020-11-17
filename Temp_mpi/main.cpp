#include <iostream>
#include <cstring>
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
    std::string message;
    MPI_Status status = MPI_Status();
	 MPI_Init(&argc, &argv);
	 MPI_Comm_size(MPI_COMM_WORLD, p);
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank );
    if (myrank == 0)
    MPI_Send("Hello", strlen("Hello"), MPI_CHAR, 1, 0, MPI_COMM_WORLD);
    else
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        std::cout<<"Hello"<<std::endl;
    }
	MPI_Finalize();
}
