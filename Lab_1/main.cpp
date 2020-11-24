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

//1 2 3
//4 5 6
//7 8 9
int main(int argc, char** argv)
{
    std::string message="Info";
    MPI_Status status = MPI_Status();
    MPI_Init(&argc, &argv);
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int n = sqrt(p);
    assert(n*n==p);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank==0)//1
    {
        std::cout<<message<<myrank<<std::endl;
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+1, 0, MPI_COMM_WORLD);
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+n, 0, MPI_COMM_WORLD);
    }
    else if(myrank+1==n)//3
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-1, 0, MPI_COMM_WORLD, &status);
        std::cout<<message<<myrank<<std::endl;
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+n, 0, MPI_COMM_WORLD);
    }
    else if(myrank==n*(n-1))//7
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-n, 0, MPI_COMM_WORLD, &status);
        std::cout<<message<<myrank<<std::endl;
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+1, 0, MPI_COMM_WORLD);
    }
    else if(myrank==n*n-1)//9
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-n, 0, MPI_COMM_WORLD, &status);
        std::cout<<message<<myrank<<std::endl;
    }
    else if(myrank > 0 && myrank <n-1)//2
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-1, 0, MPI_COMM_WORLD, &status);
        std::cout<<message<<myrank<<std::endl;
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+1, 0, MPI_COMM_WORLD);
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+n, 0, MPI_COMM_WORLD);
    }
    else if((myrank)%n == 0)//4
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-n, 0, MPI_COMM_WORLD, &status);
        std::cout<<message<<myrank<<std::endl;
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+1, 0, MPI_COMM_WORLD);
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+n, 0, MPI_COMM_WORLD);
    }
    else if ((myrank+1)%n == 0)//6
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-n, 0, MPI_COMM_WORLD, &status);
        std::cout<<message<<myrank<<std::endl;
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+n, 0, MPI_COMM_WORLD);
    }
    else if (myrank - n*(n-1)>=0&&myrank - n*(n-1)<n-1)//8
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-n, 0, MPI_COMM_WORLD, &status);
        std::cout<<message<<myrank<<std::endl;
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+1, 0, MPI_COMM_WORLD);
    }
    else //5
    {
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, myrank-n, 0, MPI_COMM_WORLD, &status);
        std::cout<<message<<myrank<<std::endl;
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+1, 0, MPI_COMM_WORLD);
        MPI_Send(message.c_str(), message.length(), MPI_CHAR, myrank+n, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}
