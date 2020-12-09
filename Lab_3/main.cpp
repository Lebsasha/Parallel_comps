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

#include </usr/include/x86_64-linux-gnu/mpich/mpi.h>
#include <cmath>
///usr/include/x86_64-linux-gnu/mpich/mpi.h
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

void send(int n, const void* data, int process, MPI_Datatype type = MPI_INT);

void receive(int n, void* data, MPI_Datatype type = MPI_INT, int root = q);


/// 180	200	 220	 240
/// 404	456	 508	 560
/// 628	712	 796	 880
/// 852	968	1084	1200


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    status = MPI_Status();
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    const int n = 4;
    q = 0;
    assert(q < p);
    int* A;
    int* B;
    int* C;
    int* D;
    int* A_row = new int[n];
    int* B_col = new int[n];
    int* C_col = new int[n];
    int d = 0;
    int* buffer;
    bool signal_i = true;
    bool signal_j = true;
    bool* process_on_i_active;
    bool* process_on_j_active;
    int general_offset_i;
    int general_offset_j;
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
            D[i] = 666;
        }
        buffer = new int[n];
        process_on_i_active = new bool[p];
        process_on_j_active = new bool[p];
        for (int k = 0; k < p; ++k)
        {
            process_on_i_active[k] = true;
            process_on_j_active[k] = true;
        }
        general_offset_i = 0;
        general_offset_j = 0;
    }

    begg:
    if (my_rank != q)
    {
        //MPI_Recv(&signal_j, 1, MPI_CHAR, q, 0, MPI_COMM_WORLD, &status);
        while (signal_j)//if (have work)
        {
            /// C1
            receive(n, B_col);
            receive(n, C_col);
            View(B_col, 1, n);
            View(C_col, 1, n);
            /// EC1
            while (signal_i)
            {
                /// C2
                receive(n, A_row);
                d = 0;
                for (int k = 0; k < n; ++k)
                {
                    d += A_row[k] * (B_col[k] + C_col[k]);
                }
                cout << my_rank << " = " << d << endl;
                send(1, &d, q);
                //MPI_Gather(&d, 1, MPI_INT, buffer, p, 1, q, MPI_COMM_WORLD);
                /// EC2
                /// C3
                receive(1, &signal_i, MPI_CHAR);
                /// EC3
            }
            ///C4
            receive(1, &signal_j, MPI_CHAR);
            /// EC4
        }
        //else
        // goto computatuions_off;
        // goto begg;
    }
    else // I am q
    {
        while (true)
        {
//            int i = 0;
//            int j = 0;
            /// C1
            cout << "C1 in q" << endl;
            for (int k = 0; process_on_j_active[k] && k < p; ++k)/// B_Col
            {
                for (int m = 0; m < n; ++m)
                    buffer[m] = B[m * n + general_offset_j + k];
                if (k != q)
                    send(n, buffer, k);
                else
                    MPI_Sendrecv(buffer, n, MPI_INT, k, 0, B_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
            }
            for (int k = 0; process_on_j_active[k] && k < p; ++k)/// C_Col
            {
                for (int m = 0; m < n; ++m)
                    buffer[m] = C[m * n + general_offset_j + k];
                if (k != q)
                    MPI_Send(buffer, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                else
                    MPI_Sendrecv(buffer, n, MPI_INT, k, 0, C_col, n, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
            }
            /// EC1
            while (true)// TODO
            {
                /// C2
                cout << "C2 in q" << endl;
                for (int k = 0; process_on_i_active[k] && k < p; ++k)/// A_row
                {
                    if (k != q)
                        MPI_Send(A + (general_offset_i + k) * n, n, MPI_INT, k, 0, MPI_COMM_WORLD);
                    else
                        MPI_Sendrecv(A + (general_offset_i + k) * n, n, MPI_INT, k, 0, A_row, n, MPI_INT, k, 0, MPI_COMM_WORLD,
                                     &status);
                }
                if (process_on_i_active[q])
                {
                    d = 0;
                    for (int k = 0; k < n; ++k)/// Compute
                    {
                        d += A_row[k] * (B_col[k] + C_col[k]);
                    }
                }
                for (int k = 0; process_on_i_active[k] && k < p; ++k)
                {
                    if (k != q)
                        receive(1, D + (general_offset_i + k) * n + general_offset_j, MPI_INT, k);
                    else
                    {
                        MPI_Sendrecv(&d, 1, MPI_INT, k, 0, D + (general_offset_i + k) * n + general_offset_j, 1, MPI_INT, k, 0,
                                     MPI_COMM_WORLD,
                                     &status);
                    }

                }
                //MPI_Gather(&d, 1, MPI_INT, buffer, p, 1, q, MPI_COMM_WORLD);
//            for (int k = 0; k < p; ++k)// Gather
//            {
//                D[general_offset_i + k] = buffer[k];
//            }
                /// EC2
                /// C3
                general_offset_i += p;
//            assert(general_offset_i < n);
                if (general_offset_i >= n)
                {
                    for (int k = 0; process_on_i_active[k] && k < p; ++k)/// If last signal is needed
                    {
                        process_on_i_active[k] = false;
                        if (k != q)
                            send(1, &process_on_i_active[k], k, MPI_CHAR); //MPI_CXX_BOOL
                    }
                    break;
                }
                for (int k = 0; k < p; ++k)/// If i is needed
                {
                    if (process_on_j_active[k] && n - general_offset_i - k > 0)
                        process_on_i_active[k] = true;
                    else
                        process_on_i_active[k] = false;

                    if (k != q)
                        send(1, &process_on_i_active[k], k, MPI_CHAR);
                }
                /// EC3
            }
            /// C4
            general_offset_i = 0;
            general_offset_j += p;
//            assert(general_offset_j < n);
            if (general_offset_j >= n)
            {
                for (int k = 0; process_on_j_active[k] && k < p; ++k)/// If last signal is needed
                {
                    process_on_j_active[k] = false;
                    if (k != q)
                        send(1, &process_on_j_active[k], k, MPI_CHAR); //MPI_CXX_BOOL
                }
                break;
            }
            for (int k = 0; k < p; ++k)/// If j is needed
            {
                if (n - general_offset_j - k > 0)
                    process_on_j_active[k] = true;
                else
                    process_on_j_active[k] = false;

                if (k != q)
                    send(1, &process_on_j_active[k], k, MPI_CHAR);
            }
            /// EC4
            copy(process_on_j_active, process_on_j_active + p, process_on_i_active);
        }
    }

    if (my_rank == q)
    {
        View(D, n, n);
        int B_D[] = {180, 200, 220, 240, 404, 456, 508, 560, 628, 712, 796, 880, 852, 968, 1084, 1200};
        assert(std::equal(D, D + n * n, B_D));
    }

    computatuions_off:
    MPI_Finalize();
    signal_i = true;
    if (signal_i)
        return 0;


    general_offset_i = 0;
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

}

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
// //    View(A, n, m);
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

string to_str(const MPI_Datatype& type)
{
    switch (type)
    {
        case MPI_INT:
            return "int";
        case MPI_CHAR:
            return "char";
        default:
            return "None";
    }
}

void receive(const int n, void* data, const MPI_Datatype type, const int root)
{
    cout << "trying receive " << n << ' ' << to_str((MPI_Datatype) type);
    if (root != q)
        cout << " from " << root;
    cout << endl;
    MPI_Recv(data, n, type, root, 0, MPI_COMM_WORLD, &status);
}

void send(const int n, const void* data, int process, const MPI_Datatype type)
{
    cout << "trying send " << n << " to process " << process << " with type " << to_str(type) << endl;
    MPI_Send(data, n, type, process, 0, MPI_COMM_WORLD);
}

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
