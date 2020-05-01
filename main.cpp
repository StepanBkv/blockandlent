#include <iostream>
#include "math.h"
#include <omp.h>

using namespace std;
int N;

long double Sum(int** a, int size) {
    long double sum = 0;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            sum += a[i][j];
        }
    }
    return sum;
}

void multipl(int** C, int** A, int** B)
{
#pragma omp parallel shared(C,B)
    {
        int numj = omp_get_num_threads();
        int kol = N / numj;
        int tekj = omp_get_thread_num();
        if (tekj != numj)
        {

            for (int i = kol * tekj; i < kol * tekj + kol; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int c = 0; c < N; c++)
                    {
                        C[i][j] += A[i][c] * B[c][j];
                    }

                }
            }
        }
        else
        {

            for (int i = kol * tekj; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int c = 0; c < N; c++)
                    {
                        C[i][j] += A[i][c] * B[c][j];
                    }
                }
            }
        }
    }
}

void multipl1(int** C, int** A, int** B)
{

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int c = 0; c < N; c++)
            {
                C[i][j] += A[i][c] * B[c][j];
            }
        }
    }
}

void Checkerboard_Block(int** matr1, int** matr2, int** matrRez)
{
    int i, j, k, thread;
    int number_threads = 1;
    omp_set_num_threads(number_threads);
    int Set = int(sqrt((double)number_threads));
    int Block = N / Set;

#pragma omp parallel shared(matrRez) private(i,j,k,thread)
    {
        thread = omp_get_thread_num();
        int Col = thread % Set;
        int Row = thread / Set;
        for (int s = 0; s < Set; s++) {
            for (i = Row * Block; i < (Row + 1) * Block; i++)
                for (j = Col * Block; j < (Col + 1) * Block; j++)
                    for (k = s * Block; k < (s + 1) * Block; k++)
                        matrRez[i][j] += matr1[i][k] * matr2[k][j];
        }
    }

}

void obnul_element(int **mas, int N){
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mas[i][j] = 0;
        }
    }
}
int main()
{
    setlocale(LC_ALL, "RUS");
    cout << "Введите размер матрицы: ";
    cin >> N;
    int** mas_A = new int* [N];
    int** mas_B = new int* [N];
    int** mas_Res = new int* [N];
    for (int i = 0; i < N; i++)
    {
        mas_A[i] = new int[N];
        mas_B[i] = new int[N];
        mas_Res[i] = new int[N];
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mas_A[i][j] = 1;
            mas_B[i][j] = 1;
            mas_Res[i][j] = 0;
        }
    }
    double t0, t1;
    long double sum = 0;

    t0 = omp_get_wtime();
    multipl1(mas_Res, mas_A, mas_B);
    t1 = omp_get_wtime();
    cout << "Время выполнения последовательного алгоритма: " << t1 - t0 << endl;
    sum = Sum(mas_Res, N);
    cout << "Сумма элементов = " << sum << "\n\n";
    obnul_element(mas_Res, N);
    t0 = omp_get_wtime();
    multipl(mas_Res, mas_A, mas_B);
    t1 = omp_get_wtime();
    cout << "Время выполнения ленточного разбиения: " << t1 - t0 << endl;
    sum = Sum(mas_Res, N);
    cout << "Сумма элементов = " << sum << "\n\n";
    obnul_element(mas_Res, N);
    t0 = omp_get_wtime();
    Checkerboard_Block(mas_A, mas_B, mas_Res);
    t1 = omp_get_wtime();
    cout << "Время выполнения блочного разбиения: " << t1 - t0 << endl;
    sum = Sum(mas_Res, N);
    cout << "Сумма элементов = " << sum << "\n\n";
    return 0;
}