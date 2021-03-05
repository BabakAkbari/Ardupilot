#pragma once


#include <cmath>
#include <AP_Common/AP_Common.h>
#include <stdio.h>
#include <AP_Math/AP_Math.h>



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

    Matrix<T> transpose() const;

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

    float norm() const;
    void set(size_t i, size_t j, T value);
    T get(size_t i, size_t j) const;
    Matrix<T>& operator<< (T x);
    Matrix<T>& operator, (T x);
    Matrix<T>& identity();
    Matrix<T>& zero();
    Matrix<T> minor_matrix(size_t r, size_t c) const;
    T minor_a(size_t r, size_t c) const;
    T algcompl(size_t r, size_t c) const;
    T det() const;
    Matrix<T> inv() const;
    T maxCoeff() const;
    bool hasNaN() const;
    void printmat() const;
};

