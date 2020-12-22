#include <iostream>
#include <mpich/mpi.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Status status = MPI_Status();
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    char buffer[4] = {0};
    if (my_rank == 0)
    {
        MPI_Send("1", 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
        system("sleep 2");
        MPI_Send("2", 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
        system("sleep 2");
        MPI_Send("3", 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
        std::cout<<"first"<<std::endl;
    }
    else
    {
        system("sleep 7");
        MPI_Recv(buffer, 30, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        for(int i=0; i < 4; ++i)
        std::cout << buffer[i] << std::endl;
        MPI_Recv(buffer, 30, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        for(int i=0; i < 4; ++i)
        std::cout << buffer[i] << std::endl;
        MPI_Recv(buffer, 30, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        for(int i=0; i < 4; ++i)
        std::cout << buffer[i] << std::endl;
    }
    MPI_Finalize();
}
