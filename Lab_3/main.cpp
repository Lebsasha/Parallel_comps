#include <iostream>
#include <cstring>
#include <cassert>
#include <iomanip>

#ifdef CBB
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#include </usr/include/x86_64-linux-gnu/mpich/mpi.h>
#pragma GCC diagnostic pop
#else // CBB

#include <mpich/mpi.h>

#endif // CBB
using namespace std;

/**
 * @result Prints pMatrix in matrix n_row × n_col form numbers.
 */
template<typename T>
void View(T* pMatrix, size_t n_row, size_t n_col) noexcept;
MPI_Status status;
///root process
int q;
void send(int n, const void* data, int process, MPI_Datatype type=MPI_INT);

void receive(int n, void* data, MPI_Datatype type=MPI_INT, int root=q);

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    status = MPI_Status();
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    const int n = 4;
    q = 2;
    assert(q < p);
    int* A;
    int* B;
    int* C;
    int* D;
    int* B_col = new int[n];
    int* C_col = new int[n];
    int* A_row = new int[n];
    int d = 0;
    int* buffer = new int[p];
    /// If i have work
    bool signal_j = true;
    bool signal_i = true;
    if (my_rank == q)
    {
        A = new int[n * n];
        B = new int[n * n];
        C = new int[n * n];
        D = new int[n * n];
        for (int i = 0; i < n * n; ++i)
        {
            A[i] = i;
            B[i] = i;
            C[i] = i;
            D[i] = 0;
        }
    }

    int general_offset_i = 0;
    {
        int i = 0;
        int j = 0;
        while (j < n)
        {
            if (my_rank == q)
                for (int k = 0; k < p; ++k)
                {
                    for (int m = 0; m < n; ++m)
                    {
                        B_col[m] = B[m * n + j];
                    }
                    if (k != q)
                        MPI_Send(B_col, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                    else//TODO Remove as unused?
                        MPI_Sendrecv(B_col, n, MPI_INT, k, 0, B_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
                }
            MPI_Recv(B_col, n, MPI_INT, q, 0, MPI_COMM_WORLD, &status);
            if (my_rank == q)
                for (int k = 0; k < p; ++k)
                {
                    for (int m = 0; m < n; ++m)
                    {
                        C_col[m] = C[m * n + j];
                    }
                    if (k != q)
                        MPI_Send(C_col, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                    else//TODO Remove as unused?
                        MPI_Sendrecv(C_col, n, MPI_INT, k, 0, C_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
                }
            MPI_Recv(C_col, n, MPI_INT, q, 0, MPI_COMM_WORLD, &status);

            while (i < n)
            {
                if (my_rank == q)
                    for (int k = 0; k < p; ++k)
                    {
                        if (k != q)
                            MPI_Send(A + k * n, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                        else
                            MPI_Sendrecv(A + k * n, n, MPI_INT, k, 0, A_row, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
                    }
                MPI_Recv(A_row, n, MPI_INT, q, 0, MPI_COMM_WORLD, &status);
                for (int k = 0; k < n; ++k)
                {
                    d += A_row[k] * (B_col[k] + C_col[k]);
                }
                MPI_Gather(&d, 1, MPI_INT, D + general_offset_i, p, 1, q, MPI_COMM_WORLD);
                if (my_rank == q)
                {
                    general_offset_i += p;
                    ++i;
                }
//            MPI_Scatter(B, n, MPI_INT, B_col, n, MPI_INT, q, MPI_COMM_WORLD);
//            MPI_Scatter(B, n, MPI_INT, C_col, n, MPI_INT, q, MPI_COMM_WORLD);
//            MPI_Scatter(B, n, MPI_INT, A_row, n, MPI_INT, q, MPI_COMM_WORLD);
            }
            ++j;
        }
        if (my_rank == q)
            View(D, n, n);
        MPI_Finalize();
    }

    if (my_rank != q)
    {
        MPI_Recv(&signal_j, 1, MPI_CHAR, q, 0, MPI_COMM_WORLD, &status);
        while (signal_j)//if (have work)
        {
            receive(n, B_col);
            MPI_Recv(C_col, n, MPI_INT, q, 0, MPI_COMM_WORLD, &status);
            while(signal_i)
            {
                /// C2
                MPI_Recv(A_row, n, MPI_INT, q, 0, MPI_COMM_WORLD, &status);
                for (int k = 0; k < n; ++k)
                {
                    d += A_row[k] * (B_col[k] + C_col[k]);
                }
                MPI_Gather(&d, 1, MPI_INT, buffer, p, 1, q, MPI_COMM_WORLD);
                /// EC2
                receive(1, &signal_i, MPI_INT);
            }
            receive(1, &signal_j, MPI_INT);
        }
        //else
        {
            goto computatuions_off;
        }
    }
    else // I am q
    {
        //TODO signal_j
        //int i = 0;
        //int j = 0;
        for (int k = 0; k < p; ++k)/// B_Col
        {
            for (int m = 0; m < n; ++m)
                B_col[m] = B[m * n + j];
            if (k != q)
                send(n, B_col, k);
            else//TODO Remove as unused?
                MPI_Sendrecv(B_col, n, MPI_INT, k, 0, B_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
        }
        for (int k = 0; k < p; ++k)/// C_Col
        {
            for (int m = 0; m < n; ++m)
                C_col[m] = C[m * n + j];
            if (k != q)
                MPI_Send(C_col, n, MPI_INT, k, 0, MPI_COMM_WORLD);
            else//TODO Remove as unused?
                MPI_Sendrecv(C_col, n, MPI_INT, k, 0, C_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
        }
        for (int i = 0; i < n; ++i)
        {
            /// C2
            for (int k = 0; k < p; ++k)/// A_row
            {
                if (k != q)
                    MPI_Send(A + k * n, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                else
                    MPI_Sendrecv(A + k * n, n, MPI_INT, k, 0, A_row, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
            }
            for (int k = 0; k < n; ++k)/// Compute
            {
                d += A_row[k] * (B_col[k] + C_col[k]);
            }
            MPI_Gather(&d, 1, MPI_INT, buffer, p, 1, q, MPI_COMM_WORLD);
            for (int k = 0; k < p; ++k)
            {
                D[general_offset_i + k] = buffer[k];
            }
            general_offset_i += p;
            assert(general_offset_i < n);
            /// EC2
        }
    }

    computatuions_off:
    MPI_Finalize();

//    const int m = 5;
//    const int k = 2;
//    const int q = 2;
//    assert(k < n);
//    assert(p <= m);
//
//    int my_rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//    int* A = new int[n * m];
//    int* B;
//    int* recvbuf = new int[m];
//    int* counts = new int[p];
//    int* displs = new int[p];
//    for (int* ptr = A; ptr < A + n * m; ++ptr)
//    {
//        *ptr = my_rank;
//    }
//    for (int i = 0; i < p; ++i)
//    {
//        counts[i] = m - i;
//    }
//    for (int i = 0; i < p; ++i)
//    {
//        displs[i] = i * m;
//    }
//    for (int i = 0; i < m; ++i)
//    {
//        recvbuf[i] = 0;
//    }
//    if (my_rank == q)
//    {
//        B = new int[p * m];
//        for (int i = 0; i < p; ++i)
//            for (int j = 0; j < m; ++j)
//            {
//                B[i * m + j] = i;
//            }
//        for (int i = 0; i < p; ++i)
//        {
//            if (i != q)
//                MPI_Send(B + displs[i], counts[i], MPI_INT, i, 0, MPI_COMM_WORLD);
//            else
//                MPI_Sendrecv(B + displs[i], counts[i], MPI_INT, i, 0, recvbuf, counts[i], MPI_INT, i, 0, MPI_COMM_WORLD, &status);
//        }
//    }
////    View(A, n, m);
//    //MPI_Scatterv(B, counts, displs, MPI_INT, recvbuf, counts[my_rank], MPI_INT, q, MPI_COMM_WORLD);
//    if (my_rank != q)
//        MPI_Recv(recvbuf, counts[my_rank], MPI_INT, q, 0, MPI_COMM_WORLD, &status);
//    for (int* ptr = A; ptr < A + n * m; ++ptr)
//        *ptr = 0;
//    for (int* ptr = A + k * m; ptr < A + (k + 1) * m; ++ptr, ++recvbuf)
//        *ptr = *recvbuf;
//    char* temp = new char[2];
//    if (my_rank != 0) // Следующий процесс будет выводить на экран только после прерыдущего
//        MPI_Recv((void*) temp, 1, MPI_CHAR, my_rank - 1, 0, MPI_COMM_WORLD, &status);
//    View(A, n, m);
//    if (my_rank != p - 1)
//        MPI_Send("1", 1, MPI_CHAR, my_rank + 1, 0, MPI_COMM_WORLD);
}

void receive(const int n, void* data, MPI_Datatype type, const int root)
{ MPI_Recv(data, n, type, root, 0, MPI_COMM_WORLD, &status); }

void send(const int n, const void* data, int process, MPI_Datatype type)
{ MPI_Send(data, n, type, process, 0, MPI_COMM_WORLD); }

template<typename T>
void View(T* pMatrix, size_t n_row, const size_t n_col) noexcept
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
