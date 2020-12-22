#include <iostream>
#include <cstring>
#include <cassert>
#include <iomanip>

#include </usr/include/x86_64-linux-gnu/mpich/mpi.h>
using namespace std;

/**
 * @result Prints pMatrix in matrix n_row × n_col form numbers.
 */
template<typename T>
void View(T* pMatrix, size_t n_row, size_t n_col);

MPI_Status status;

///root process
int q;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    status = MPI_Status();
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    const size_t n = 2000;
    q = 0;
    assert(q < p);
    int* A, *B, *C, *D, *buffer;
    int* A_row = new int[n];
    int* B_col = new int[n];
    int* C_col = new int[n];
    int result = 0;
    bool signal_i = true;
    bool signal_j = true;
    int last_active_process;
    size_t general_offset_i;
    size_t general_offset_j;
    double computation_time;

    if (my_rank != q)
    {
        while (signal_j)//if (have work)
        {
            MPI_Recv(B_col, n, MPI_INT, q, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(C_col, n, MPI_INT, q, 0, MPI_COMM_WORLD, &status);
            while (signal_i)
            {
                MPI_Recv(A_row, n, MPI_INT, q, 0, MPI_COMM_WORLD, &status);
                result = 0;
                for (size_t k = 0; k < n; ++k)
                {
                    result += A_row[k] * (B_col[k] + C_col[k]);
                }
                MPI_Send(&result, 1, MPI_INT, q, 0, MPI_COMM_WORLD);
                MPI_Recv(&signal_i, 1, MPI_CXX_BOOL, q, 0, MPI_COMM_WORLD, &status);
            }
            signal_i=true;
            MPI_Recv(&signal_j, 1, MPI_CXX_BOOL, q, 0, MPI_COMM_WORLD, &status);
        }

        /// Work is done
    }
    else // I am q
    {

        A = new int[n * n];
        B = new int[n * n];
        C = new int[n * n];
        D = new int[n * n];
        for (size_t i = 0; i < n * n; ++i)
        {
            A[i] = i + 1;
            B[i] = i + 1;
            C[i] = i + 1;
            D[i] = 0;
        }
        buffer = new int[n];
        last_active_process=p-1;
        general_offset_i = 0;
        general_offset_j = 0;
        computation_time=MPI_Wtime();

        while (true)
        {
            if(general_offset_j%100 == 0)
                cout<<general_offset_j<<endl;
            for (int k = 0; k <= last_active_process; ++k)/// B_Col
            {
                for (size_t m = 0; m < n; ++m)
                    buffer[m] = B[m * n + general_offset_j + k];
                if (k != q)
                    MPI_Send(buffer, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                else
                    MPI_Sendrecv(buffer, n, MPI_INT, k, 0, B_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
            }
            for (int k = 0; k <= last_active_process; ++k)/// C_Col
            {
                for (size_t m = 0; m < n; ++m)
                    buffer[m] = C[m * n + general_offset_j + k];
                if (k != q)
                    MPI_Send(buffer, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                else
                    MPI_Sendrecv(buffer, n, MPI_INT, k, 0, C_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
            }
            while (true)
            {
                for (int k = 0; k <= last_active_process; ++k)/// A_row
                {
                    if (k != q)
                        MPI_Send(A + (general_offset_i) * n, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                    else
                        MPI_Sendrecv(A + (general_offset_i) * n, n, MPI_INT,k, 0,
                                        A_row, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
                }
                if (q <= last_active_process)
                {
                    result = 0;
                    for (size_t k = 0; k < n; ++k)/// Compute
                    {
                        result += A_row[k] * (B_col[k] + C_col[k]);
                    }
                }
                for (int k = 0; k <= last_active_process; ++k)
                {
                    if (k != q)
                    {
                        MPI_Recv(D + (general_offset_i) * n + general_offset_j + k, 1, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
                    }
                    else
                    {
                        MPI_Sendrecv(&result, 1, MPI_INT, k, 0, D + (general_offset_i) * n + general_offset_j+k, 1, MPI_INT, k, 0,
                                     MPI_COMM_WORLD,
                                     &status);
                    }

                }
                ++general_offset_i;
                if (general_offset_i >= n)/// If rows end
                {
                    signal_i = false;
                    for (int k = 0; k <= last_active_process; ++k)/// If last signal is needed
                    {
                        if (k != q)
                        {
                            MPI_Send(&signal_i, 1, MPI_CXX_BOOL, k, 0, MPI_COMM_WORLD);
                        }
                    }
                    break;
                }
                signal_i=true;
                for (int k = 0; k <= last_active_process; ++k)/// Continue
                {
                    if (k != q)
                    {
                        MPI_Send(&signal_i, 1, MPI_CXX_BOOL, k, 0, MPI_COMM_WORLD);
                    }
                }
            }
            general_offset_i = 0;
            general_offset_j += p;
            if (general_offset_j >= n) /// If work done
            {
                signal_j=false;
                for (int k = 0; k <= last_active_process; ++k)/// If last signal is needed
                {
                    if (k != q)
                    {
                        MPI_Send(&signal_j, 1, MPI_CXX_BOOL, k, 0, MPI_COMM_WORLD);
                    }
                }
                break;
            }
            last_active_process = n - general_offset_j > p ? p-1 : n - general_offset_j-1;
            for (int k = 0; k < p; ++k)/// If process is needed
            {
                signal_j= k<=last_active_process;
                if (k != q)
                {
                    MPI_Send(&signal_j, 1, MPI_CXX_BOOL, k, 0, MPI_COMM_WORLD);
                }
            }
        }

        /// Work is done in root process
        computation_time=MPI_Wtime()-computation_time;
        // View(D, n, n);
        cout<<"Work is done in time: "<<computation_time<<endl;
        delete[] buffer;
        delete[] A;
        delete[] B;
        delete[] C;
        delete[] D;
    }

    delete[] A_row;
    delete[] B_col;
    delete[] C_col;
    MPI_Finalize();
    return 0;
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
