#include <iostream>
#include <iomanip>
#include <omp.h>

using namespace std;

/**
 * @result Prints pMatrix in matrix n_row × n_col form numbers.
 */
template<typename T>
void View(T* pMatrix, size_t n_row, size_t n_col);


int main(int argc, char** argv)
{
    omp_set_num_threads(4);
    const size_t n = 1000;
    int* A, * B, * C, * D;
    double computation_time;
    A = new int[n * n];
    B = new int[n * n];
    C = new int[n * n];
    D = new int[n * n];
    for (int i = 0; i < (int) (n * n); ++i)
    {
        A[i] = i + 1;
        B[i] = i + 1;
        C[i] = i + 1;
        D[i] = 0;
    }
    computation_time = omp_get_wtime();

#pragma omp parallel default(none) shared(A, B, C, D)
    {
#pragma omp for schedule(static, 3)
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                for (size_t k = 0; k < n; ++k)
                    D[i * n + j] += A[i * n + k] * (B[k * n + j] + C[k * n + j]);
    }

    computation_time = omp_get_wtime() - computation_time;
    View(D, n, n);
    cerr << "Work is done in time: " << computation_time << endl;
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;
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
