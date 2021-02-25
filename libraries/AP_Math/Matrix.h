#pragma GCC optimize("O2")

#define ALLOW_DOUBLE_MATH_FUNCTIONS
#include <stdint.h>
#include "AP_Math.h"
#include <AP_HAL/AP_HAL.h>
#include <AP_HAL/AP_HAL_Boards.h>

template <typename T>
class Matrix
{
    private:
    T* elements;
    size_t row;
    size_t col;
    size_t op_counter;
    
    public:
    Matrix<T>(size_t r, size_t c);
    Matrix<T>(size_t size);
    Matrix<T>(const Matrix<T>& matrix);
    ~Matrix<T>();

    Matrix<T>& operator= (const Matrix<T>& matrix);
    Matrix<T> operator+ (const Matrix<T>& matrix);
    Matrix<T> operator- (const Matrix<T>& matrix);
    Matrix<T> operator* (const Matrix<T>& matrix);

    Matrix<T> transpose();

    void operator*= (const Matrix<T>& matrix);
    void operator*= (const T& value);
    void operator/= (const T& value);
    Matrix<T> operator* (const T& value);
    Matrix<T> operator/ (const T& value);
    friend Matrix<T> operator* (const T& value, const Matrix<T>& matrix)
    {
        Matrix result(matrix.row, matrix.col);
        for(size_t i = 0; i < matrix.row * matrix.col; i++)
        {
        result.elements[i] = matrix.elements[i] * value;
        }

        return result;   
    }

    float norm();
    void set(size_t i, size_t j, T value);
    T get(size_t i, size_t j);
    Matrix<T>& operator<< (T x);
    Matrix<T>& operator, (T x);
    Matrix<T>& identity();
    Matrix<T>& zero();
    Matrix<T> minor_matrix(size_t r, size_t c);
    T minor(size_t r, size_t c);
    T algcompl(size_t r, size_t c);
    T det();
    Matrix<T> inv();
};

