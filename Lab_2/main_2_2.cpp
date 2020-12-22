#include <iostream>
#include <cstring>
#include <cassert>
#include <iomanip>
#include <mpich/mpi.h>

using namespace std;

/**
 * @result Prints pMatrix in matrix n_row × n_col form numbers.
 */
template<typename T>
void View(T* pMatrix, size_t n_row, size_t n_col);

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Status status = MPI_Status();
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    const int n = 4;
    const int m = 5;
    const int k = 2;
    const int q = 2;
    assert(k < n);
    assert(p <= m);
    assert(q < p);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int* A = new int[n * m];
    int* B;
    int* recvbuf = new int[m];
    int* counts = new int[p];
    int* displs = new int[p];
    for (int* ptr = A; ptr < A + n * m; ++ptr)
        *ptr = my_rank;
    for (int i = 0; i < p; ++i)
    {
        counts[i] = m - i;
        displs[i] = i * m;
    }
    for (int i = 0; i < m; ++i)
        recvbuf[i] = 0;
    if (my_rank == q)
    {
        B = new int[p * m];
        for (int i = 0; i < p; ++i)
            for (int j = 0; j < m; ++j)
                B[i * m + j] = i;
        for (int i = 0; i < p; ++i)
        {
            if (i != q)
                MPI_Send(B + displs[i], counts[i], MPI_INT, i, 0, MPI_COMM_WORLD);
            else
                MPI_Sendrecv(B + displs[i], counts[i], MPI_INT, i, 0, recvbuf, counts[i], MPI_INT, i, 0, MPI_COMM_WORLD, &status);
        }
    }
//    View(A, n, m);
    //MPI_Scatterv(B, counts, displs, MPI_INT, recvbuf, counts[my_rank], MPI_INT, q, MPI_COMM_WORLD);
    if (my_rank != q)
        MPI_Recv(recvbuf, counts[my_rank], MPI_INT, q, 0, MPI_COMM_WORLD, &status);
    for (int* ptr = A; ptr < A + n * m; ++ptr)
        *ptr = 0;
    for (int* ptr = A + k * m; ptr < A + (k + 1) * m; ++ptr, ++recvbuf)
        *ptr = *recvbuf;
    char* temp = new char[2];
    if (my_rank != 0) // Следующий процесс будет выводить на экран только после прерыдущего
        MPI_Recv((void*) temp, 1, MPI_CHAR, my_rank - 1, 0, MPI_COMM_WORLD, &status);
    View(A, n, m);
    if (my_rank != p - 1)
        MPI_Send("1", 1, MPI_CHAR, my_rank + 1, 0, MPI_COMM_WORLD);
    MPI_Finalize();
}

template<typename T>
void View(T* pMatrix, size_t n_row, const size_t n_col)
{
    while (n_row--)
    {
        for (size_t j = 0; j < n_col - 1; ++j)
        {
            std::cout << std::setw(3) << (*pMatrix) << " ";
            pMatrix++;
        }
        std::cout << std::setw(3) << *pMatrix++;
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
