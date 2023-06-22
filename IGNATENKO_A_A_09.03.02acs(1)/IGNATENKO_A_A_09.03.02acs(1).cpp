#include <iostream>
#include <ctime>
#include <cstdlib>
#include <omp.h>

const int MATRIX_SIZE = 1000; // Размер матрицы

void matrixMultiplicationSequential(float** A, float** B, float** C) {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            C[i][j] = 0.0f;
            for (int k = 0; k < MATRIX_SIZE; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void matrixMultiplicationParallel(float** A, float** B, float** C) {
#pragma omp parallel for
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            C[i][j] = 0.0f;
            for (int k = 0; k < MATRIX_SIZE; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main() {
    // Создание и заполнение матрицы A
    float** A = new float* [MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        A[i] = new float[MATRIX_SIZE];
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            A[i][j] = static_cast<float>(rand()) / RAND_MAX; // Заполнение случайным числом в диапазоне [0, 1]
        }
    }

    // Создание и заполнение матрицы B для одинарной точности
    float** Bsgl = new float* [MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        Bsgl[i] = new float[MATRIX_SIZE];
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            Bsgl[i][j] = static_cast<float>(rand()) / RAND_MAX; // Заполнение случайным числом в диапазоне [0, 1]
        }
    }

    // Создание и заполнение матрицы B для двойной точности
    double** Bdbl = new double* [MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        Bdbl[i] = new double[MATRIX_SIZE];
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            Bdbl[i][j] = static_cast<double>(rand()) / RAND_MAX; // Заполнение случайным числом в диапазоне [0, 1]
        }
    }

    // Создание матрицы C для одинарной точности
    float** Csgl = new float* [MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        Csgl[i] = new float[MATRIX_SIZE];
    }

    // Создание матрицы C для двойной точности
    double** Cdbl = new double* [MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        Cdbl[i] = new double[MATRIX_SIZE];
    }

    // Массивы для хранения времени выполнения для разных значений количества потоков
    double elapsedTimesSeq[5];
    double elapsedTimesPar[5];

    // Измерение времени умножения матриц для разных значений количества потоков
    int numThreads[] = { 2, 4, 8, 16, 32 };
    for (int i = 0; i < 5; ++i) {
        int numThread = numThreads[i];
        omp_set_num_threads(numThread); // Установка количества потоков

        // Умножение матриц для одинарной точности (последовательно)
        clock_t start = clock();
        matrixMultiplicationSequential(A, Bsgl, Csgl);
        clock_t end = clock();
        double elapsedTimeSeq = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        // Умножение матриц для одинарной точности (параллельно)
        start = clock();
        matrixMultiplicationParallel(A, Bsgl, Csgl);
        end = clock();
        double elapsedTimePar = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        // Запись времени выполнения в массивы
        elapsedTimesSeq[i] = elapsedTimeSeq;
        elapsedTimesPar[i] = elapsedTimePar;
    }

    // Вывод результатов в таблицу
    std::cout << "Num Threads\tElapsed Time (Sequential)\tElapsed Time (Parallel)" << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::cout << numThreads[i] << "\t\t" << elapsedTimesSeq[i] << " seconds\t\t" << elapsedTimesPar[i] << " seconds" << std::endl;
    }

    // Освобождение памяти
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        delete[] A[i];
        delete[] Bsgl[i];
        delete[] Bdbl[i];
        delete[] Csgl[i];
        delete[] Cdbl[i];
    }
    delete[] A;
    delete[] Bsgl;
    delete[] Bdbl;
    delete[] Csgl;
    delete[] Cdbl;

    return 0;
}