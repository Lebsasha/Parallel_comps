#include <iostream>
#include <cstring>
#include <cassert>
#include <iomanip>

#include </usr/include/x86_64-linux-gnu/mpich/mpi.h> //TODO
using namespace std;

/**
 * @result Prints pMatrix in matrix n_row × n_col form numbers.
 */
template<typename T>
void View(T* pMatrix, size_t n_row, size_t n_col);

MPI_Status status;

///root process
int q;

void send(int n, const void* data, int process, MPI_Datatype type = MPI_INT);

void receive(int n, void* data, MPI_Datatype type = MPI_INT, int root = q);

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    status = MPI_Status();
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    const int n = 1257;
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
    int general_offset_i;
    int general_offset_j;
    double computation_time;
    if (my_rank == q)
    {
        A = new int[n * n];
        B = new int[n * n];
        C = new int[n * n];
        D = new int[n * n];
        for (int i = 0; i < n * n; ++i)
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
    }

    if (my_rank != q)
    {
        while (signal_j)//if (have work)
        {
            receive(n, B_col);
            receive(n, C_col);
            while (signal_i)
            {
                receive(n, A_row);
                result = 0;
                for (int k = 0; k < n; ++k)
                {
                    result += A_row[k] * (B_col[k] + C_col[k]);
                }
                send(1, &result, q);
                receive(1, &signal_i, MPI_CXX_BOOL);
            }
            signal_i=true;
            receive(1, &signal_j, MPI_CXX_BOOL);
        }

        /// Work is done
    }
    else // I am q
    {
        while (true)
        {
            if(general_offset_j%100 == 0)
            cout<<general_offset_j<<endl;
            for (int k = 0; k <= last_active_process; ++k)/// B_Col
            {
                for (int m = 0; m < n; ++m)
                    buffer[m] = B[m * n + general_offset_j + k];
                if (k != q)
                    MPI_Send(buffer, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                else
                    MPI_Sendrecv(buffer, n, MPI_INT, k, 0, B_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
            }
            for (int k = 0; k <= last_active_process; ++k)/// C_Col
            {
                for (int m = 0; m < n; ++m)
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
                        MPI_Sendrecv(A + (general_offset_i) * n, n, MPI_INT, k, 0, A_row, n, MPI_INT, k, 0, MPI_COMM_WORLD,
                                     &status);
                }
                if (q <= last_active_process)
                {
                    result = 0;
                    for (int k = 0; k < n; ++k)/// Compute
                    {
                        result += A_row[k] * (B_col[k] + C_col[k]);
                    }
                }
                for (int k = 0; k <= last_active_process; ++k)
                {
                    if (k != q)
                        receive(1, D + (general_offset_i) * n + general_offset_j+k, MPI_INT, k);
                    else
                    {
                        MPI_Sendrecv(&result, 1, MPI_INT, k, 0, D + (general_offset_i) * n + general_offset_j+k, 1, MPI_INT, k, 0,
                                     MPI_COMM_WORLD,
                                     &status);
                    }

                }
                general_offset_i += 1;
                if (general_offset_i >= n)/// If process is ended by i
                {
                    signal_i = false;
                    for (int k = 0; k <= last_active_process; ++k)/// If last signal is needed
                    {
                        if (k != q)
                            send(1, &signal_i, k, MPI_CXX_BOOL);
                    }
                    break;
                }
                signal_i=true;
                for (int k = 0; k <= last_active_process; ++k)/// Continue
                {
                    if (k != q)
                        send(1, &signal_i, k, MPI_CXX_BOOL);
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
                        send(1, &signal_j, k, MPI_CXX_BOOL);
                }
                break;
            }
            last_active_process = n - general_offset_j > p ? p-1 : n - general_offset_j-1;
            for (int k = 0; k < p; ++k)/// If j is needed
            {
                signal_j= k<=last_active_process;
                if (k != q)
                    send(1, &signal_j, k, MPI_CXX_BOOL);
            }
        }

        /// Work is done in root process
        computation_time=MPI_Wtime()-computation_time;
        // View(D, n, n);
        cout<<"Work is off in time: "<<computation_time<<endl;
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

void receive(const int n, void* data, const MPI_Datatype type, const int root)
{
    MPI_Recv(data, n, type, root, 0, MPI_COMM_WORLD, &status);
}

void send(const int n, const void* data, int process, const MPI_Datatype type)
{
    MPI_Send(data, n, type, process, 0, MPI_COMM_WORLD);
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
