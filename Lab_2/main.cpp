#include <iostream>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iomanip>
#include <mutex>
#ifdef CBB
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#include </usr/include/x86_64-linux-gnu/mpich/mpi.h>
#pragma GCC diagnostic pop

#else // CBB
#include <mpich/mpi.h>

#endif // CBB
#define SINGL_INIT(number) bool i ## number=true
#define SINGL_B(number) if(i ## number){i##number=false
#define SINGL_E(number) }i ## number
using namespace std;
template<typename T>
void View (T* pArray, size_t Nstr, const size_t Nstb) noexcept;

int main(int argc, char** argv)
{
    SINGL_INIT(66);
    SINGL_B(66);
    SINGL_E(66);
    MPI_Init(&argc, &argv);
    MPI_Status status = MPI_Status();
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    const int n = 4;
    const int m = 4;
    const int k = 2;
    const int q = 1;
    assert(p>1);//TODO Temp
    assert(k < n);
    assert(p<=m);
    assert(q<p);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int* A=new int [n*m];
    int* B;
    int* temp=new int [m];
    int* counts=new int[p];
    int* smeschenie=new int[p];
    if(my_rank==q)
    {
        B = new int [p*m];
        for(int i=0; i <p; ++i)
        for(int j=0; j<m; ++j)
        {
            B[i*m+j]=i;
        }
    }
    for (int* ptr = A + n * m - 1; ptr >= A; --ptr)
    {
        *ptr=my_rank;
    }
    for (int i = 0; i < p; ++i)
    {
        counts[i]=m-i;
    };;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    for (int i = 0; i < m; ++i)
    {
        temp[i]=0;
    }
    for (int i = 0; i < p; ++i)
    {
        smeschenie[i]=i*m;
    }
//    View(A, n, m);
    MPI_Scatterv(B, counts, smeschenie, MPI_INT, temp, counts[my_rank-1], MPI_INT, q, MPI_COMM_WORLD);

    for (int* ptr = A + n * m - 1; ptr >= A; --ptr)
    {
        *ptr=0;
    }
    for(int* ptr=A+k*m; ptr < A+(k+1)*m; ++ptr, ++temp)
        *ptr=*temp;
    system(("sleep "+to_string(my_rank)).c_str());
    View(A, n, m);
//    if(my_rank%n!=0)
//        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, my_rank-1, 0, MPI_COMM_WORLD, &status);
//    if(my_rank>n-1)
//        MPI_Recv((void*) message.c_str(), 20, MPI_CHAR, my_rank-n, 0, MPI_COMM_WORLD, &status);
//    std::cout<<message<<my_rank<<std::endl;
//    if((my_rank+1)%n!=0)
//        MPI_Send(message.c_str(), message.length(), MPI_CHAR, my_rank+1, 0, MPI_COMM_WORLD);
//    if(my_rank<n*(n-1))
//        MPI_Send(message.c_str(), message.length(), MPI_CHAR, my_rank+n, 0, MPI_COMM_WORLD);
    MPI_Finalize();
}
template<typename T>
void View (T* pArray, size_t Nstr, const size_t Nstb) noexcept
{
    //s.lock();
    while (Nstr--)
    {
        for (size_t j = 0; j < Nstb-1; ++j)
        {
            std::cout<<std::setw(3)<<(*pArray)<<" ";
            pArray++;
        }
        std::cout<<std::setw(3)<<*pArray++;
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    //s.unlock();
}
