#include <iostream>
#include <cmath>

using namespace std;

void PrintVector(double* x, int N) {
    for (int i = 0; i < N; ++i) {
        cout << x[i] << " ";
    }
    cout << endl;
    cout << endl;
}

class eigenvector {
private:
    double lamda; // собственное число
    double counter; // кол-во итераций потребовавшихся на вычисление
    double* vector;
    int N;

public:
    eigenvector(double l, double cou, double* m_X, int N) : lamda(l), counter(cou), vector(m_X), N(N) {};

    double GetLamda() {
        return lamda;
    }
    double GetCounter() {
        return counter;
    }
    double* GetVector() {
        return vector;
    }

    void Print() {
        cout << "steps:" << counter << endl;
        cout << "lambda:" << lamda << endl;
        cout << "vector: ";
        PrintVector(vector, N);
    }
};



void Printmatrix(double** m_Arr, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << m_Arr[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// Функция для создания единичной матрицы
double** createIdentityMatrix(int size) {
    double** identityMatrix = new double* [size];
    for (int i = 0; i < size; ++i) {
        identityMatrix[i] = new double[size];
        for (int j = 0; j < size; ++j) {
            identityMatrix[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    return identityMatrix;
}
double** vector_to_matrix(double* m_x, int N) {
    double** m_Arr = new double* [N];
    for (int i = 0; i < N; ++i) {
        m_Arr[i] = new double[N];
        for (int j = 0; j < N; ++j) {
            m_Arr[i][j] = m_x[i] * m_x[j];
        }
    }
    return m_Arr;
}

// Функция для освобождения памяти, выделенной под матрицу
void deleteMatrix(double** matrix, int size) {
    for (int i = 0; i < size; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}



double Scalar_Product(double* m_X, double* m_Y, int N) {
    double result = 0;
    for (int i = 0; i < N; ++i) {
        result += m_X[i] * m_Y[i];
    }
    return result;
}
double Norma(double* m_X, double* m_Y, int N) {
    return sqrt(Scalar_Product(m_X, m_Y, N));
}
double** MultiplyMatrixScalar(double** m_Arr, double m_X, int N) {
    double** result = new double* [N];
    for (int i = 0; i < N; ++i) {
        result[i] = new double[N];
        for (int j = 0; j < N; ++j) {
            result[i][j] = m_Arr[i][j] * m_X;
        }
    }
    return result;
}

double* MultiplyMatrixVector(double** m_Arr, double* m_X, int N) {
    double* result = new double[N];
    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for (int j = 0; j < N; ++j) {
            result[i] += m_Arr[i][j] * m_X[j];
        }
    }
    return result;
}
double* MultiplyScalarVector(double a, double* m_X, int N) {
    double* result = new double[N];
    for (int i = 0; i < N; ++i) {
        result[i] = m_X[i] * a;
    }
    return result;
}
double* Discrepancy(double** m_Arr, double* uk1, double lamda, int N) {
    double* a = MultiplyMatrixVector(m_Arr, uk1, N);
    double* b = MultiplyScalarVector(lamda, uk1, N);
    double* result = new double[N];
    for (int i = 0; i < N; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

double Scalar_Product(double* m_X, int N) {
    double result = 0;
    for (int i = 0; i < N; ++i) {
        result += m_X[i] * m_X[i];
    }
    return result;
}
double Norma(double* m_X, int N) {
    return sqrt(Scalar_Product(m_X, N));
}
double** Matrix_Difference(double** m_Arr, double** m_Brr, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            m_Arr[i][j] = m_Arr[i][j] - m_Brr[i][j];
        }
    }
    return m_Arr;
}

double Power_Method(double** m_Arr, int N, double* m_Y, eigenvector** vectors, int id, double eps) {
    int counter = 0;
    double* uk = MultiplyScalarVector(1 / Norma(m_Y, N), m_Y, N);
    double* yk1 = MultiplyMatrixVector(m_Arr, uk, N);
    double lamda = Scalar_Product(yk1, uk, N);
    double* uk1 = MultiplyScalarVector(1 / Norma(yk1, N), yk1, N);
    while (Norma(Discrepancy(m_Arr, uk1, lamda, N), N) > eps) {
        if (counter > 1000000000) {
            cerr << "ne shod";
            exit(0);
        }
        else {
            counter++;
            uk = MultiplyScalarVector(1 / Norma(yk1, N), yk1, N);
            yk1 = MultiplyMatrixVector(m_Arr, uk, N);
            lamda = Scalar_Product(yk1, uk, N);
            uk1 = MultiplyScalarVector(1 / Norma(yk1, N), yk1, N);
        }
    }
    vectors[id] = new eigenvector(lamda, counter, uk1, N);

    return lamda;
}
double Power_Method_for_minim(double** m_Arr, int N, double* m_Y, double a, eigenvector** vectors, double eps) {
    double** E = createIdentityMatrix(N);
    E = MultiplyMatrixScalar(E, a, N);
    E = Matrix_Difference(E, m_Arr, N);

    int counter = 0;
    double* uk = MultiplyScalarVector(1 / Norma(m_Y, N), m_Y, N);
    double* yk1 = MultiplyMatrixVector(E, uk, N);
    double lamda = Scalar_Product(yk1, uk, N);
    double* uk1 = MultiplyScalarVector(1 / Norma(yk1, N), yk1, N);
    while (Norma(Discrepancy(E, uk1, lamda, N), N) > eps) {
        if (counter > 1000) {
            cerr << "Метод не сходится";
            exit(0);
        }
        else {
            counter++;
            uk = MultiplyScalarVector(1 / Norma(yk1, N), yk1, N);
            yk1 = MultiplyMatrixVector(E, uk, N);
            lamda = Scalar_Product(yk1, uk, N);
            uk1 = MultiplyScalarVector(1 / Norma(yk1, N), yk1, N);
        }
    }

    vectors[1] = new eigenvector(-(lamda - a), counter, uk1, N);
}

double** MultiplyMatrixMatrix(double** m_arr, double** m_barr, int N) {
    // Создаем новую матрицу для результата
    double** result = new double* [N];
    for (int i = 0; i < N; ++i) {
        result[i] = new double[N];
    }

    // Выполняем умножение матриц
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i][j] = 0.0;
            for (int k = 0; k < N; ++k) {
                result[i][j] += m_arr[i][k] * m_barr[k][j];
            }
        }
    }

    return result;
}
void Power_method_for_all_numbers(double** m_Arr, int N, double* Y, eigenvector** vectors, double eps) {
    for (int i = 0; i < N; ++i) {
        Power_Method(m_Arr, N, Y, vectors, i,eps);
        vectors[i]->Print();
        double* b = vectors[i]->GetVector();
        double** P = Matrix_Difference(createIdentityMatrix(N), vector_to_matrix(b, N), N);

        double** A = MultiplyMatrixMatrix(MultiplyMatrixMatrix(P, m_Arr, N), P, N);
        m_Arr = A;
    }
}

void maxim_minim(double** m_Arr, int N, double* Y, eigenvector** vectors,double eps) {
    cout << "max lambda:" << endl;
    Power_Method(m_Arr, N, Y, vectors, 0, eps);
    vectors[0]->Print();
    double a = vectors[0]->GetLamda();
    Power_Method_for_minim(m_Arr, N, Y, a, vectors,eps);
    cout << "min lambda" << endl;
    vectors[1]->Print();
}

int main() {
   
    system("chcp 65001");
    int N;
    double eps;
    cin >> N>>eps;
    double** m_Arr = new double* [N];
    for (int i = 0; i < N; ++i) {
        m_Arr[i] = new double[N];
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cin>> m_Arr[i][j];
        }
    }

    double* m_Y = new double[N];
    for (int i = 0; i < N; ++i) {
        cin>> m_Y[i];
    }
    eigenvector** Sov = new eigenvector * [N];
    //Power_method_for_all_numbers(m_Arr, N, m_Y, Sov,eps);
    maxim_minim(m_Arr, N, m_Y, Sov, eps);


    return 0;
}