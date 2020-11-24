#include <iostream>
#include <cstring>
#include <cmath>
#include <cassert>
#ifdef CBB
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#include </usr/include/x86_64-linux-gnu/mpich/mpi.h>
#pragma GCC diagnostic pop
#else // CBB

#include <mpich/mpi.h>

#endif // CBB

int main(int argc, char** argv)
{
    std::string message="Info";
    MPI_Init(&argc, &argv);
    MPI_Status status = MPI_Status();
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int n = static_cast<int>(sqrt(p));
    assert(n*n==p);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(myrank%n!=0)
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-1, 0, MPI_COMM_WORLD, &status);
    if(myrank>n-1)
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-n, 0, MPI_COMM_WORLD, &status);
    std::cout<<message<<myrank<<std::endl;
    if((myrank+1)%n!=0)
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+1, 0, MPI_COMM_WORLD);
    if(myrank<n*(n-1))
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+n, 0, MPI_COMM_WORLD);
    MPI_Finalize();
}
