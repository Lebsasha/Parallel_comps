#include <iostream>
#include <cstring>
#include <cassert>
#include <iomanip>

#include <omp.h>

using namespace std;

/**
 * @result Prints pMatrix in matrix n_row × n_col form numbers.
 */
template<typename T>
void View(T* pMatrix, size_t n_row, size_t n_col);


///root process
int q;

int main(int argc, char** argv)
{
    omp_set_num_threads(4);
    int p = 1;
    const size_t n = 4;
    q = 0;
    assert(q < p);
    int* A, * B, * C, * D, * buffer;
    int* A_row = new int[n];
    int* B_col = new int[n];
    int* C_col = new int[n];
    bool signal_i = true;
    bool signal_j = true;
    int last_active_process;
    double computation_time;

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
    last_active_process = p - 1;
    computation_time = omp_get_wtime();
// TODO parallel by j

#pragma omp parallel sections default(none) shared(A, B, C, D, cout)
    {
#pragma omp section
        {
            for (size_t i = 0; i < n/2; ++i)
            {
                for (size_t j = 0; j < n; ++j)
                    for (size_t k = 0; k < n; ++k)
                    {
                        D[i * n + j] += A[i * n + k] * (B[k * n + j] + C[k * n + j]);
                    }
            }
        }
#pragma omp section
        {
            for (size_t i = n/2; i < n; ++i)
            {
                for (size_t j = 0; j < n; ++j)
                    for (size_t k = 0; k < n; ++k)
                    {
                        D[i * n + j] += A[i * n + k] * (B[k * n + j] + C[k * n + j]);
                    }
            }
        }
    }

    computation_time = omp_get_wtime() - computation_time;
    View(D, n, n);
    cout << "Work is done in time: " << computation_time << endl;
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;

//    180 200 220 240
//    404 456 508 560
//    628 712 796 880
//    852 968 1084 1200


    delete[] buffer;
    delete[] A_row;
    delete[] B_col;
    delete[] C_col;

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
