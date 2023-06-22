#include <iostream>
#include <ctime>
#include <cstdlib>
#include <omp.h>

const int ARRAY_SIZE = 10e6; // Размер массива

int main() {
    // Создание и заполнение массива одинарной точности
    float* Asgl = new float[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; ++i) {
        Asgl[i] = static_cast<float>(rand()) / RAND_MAX; // Заполнение случайным числом в диапазоне [0, 1]
    }

    // Создание и заполнение массива двойной точности
    double* Adbl = new double[ARRAY_SIZE];
    for (int i = 0; i < ARRAY_SIZE; ++i) {
        Adbl[i] = static_cast<double>(rand()) / RAND_MAX; // Заполнение случайным числом в диапазоне [0, 1]
    }

    // Массивы для хранения времени выполнения для разных значений количества потоков
    double elapsedTimesAsgl[5];
    double elapsedTimesAdbl[5];

    // Измерение времени суммирования массивов для разных значений количества потоков
    int numThreads[] = { 2, 4, 8, 16, 32 };
    for (int i = 0; i < 5; ++i) {
        int numThread = numThreads[i];
        omp_set_num_threads(numThread); // Установка количества потоков

        // Измерение времени суммирования массива одинарной точности
        clock_t start = clock();
        float sumAsgl = 0.0f;
#pragma omp parallel for reduction(+:sumAsgl)
        for (int j = 0; j < ARRAY_SIZE; ++j) {
            sumAsgl += Asgl[j];
        }
        clock_t end = clock();
        double elapsedTimeAsgl = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        // Измерение времени суммирования массива двойной точности
        start = clock();
        double sumAdbl = 0.0;
#pragma omp parallel for reduction(+:sumAdbl)
        for (int j = 0; j < ARRAY_SIZE; ++j) {
            sumAdbl += Adbl[j];
        }
        end = clock();
        double elapsedTimeAdbl = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        // Запись времени выполнения в массивы
        elapsedTimesAsgl[i] = elapsedTimeAsgl;
        elapsedTimesAdbl[i] = elapsedTimeAdbl;
    }

    // Вывод результатов в таблицу
    std::cout << "Num Threads\tElapsed Time (single precision)\tElapsed Time (double precision)" << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::cout << numThreads[i] << "\t\t" << elapsedTimesAsgl[i] << " seconds\t\t" << elapsedTimesAdbl[i] << " seconds" << std::endl;
    }

    // Освобождение памяти
    delete[] Asgl;
    delete[] Adbl;

    return 0;
}
