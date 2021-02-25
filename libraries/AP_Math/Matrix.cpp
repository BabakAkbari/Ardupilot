#include "Matrix.h"

template <typename T>
Matrix<T>::Matrix(size_t r, size_t c) : elements{new T[r * c]}, row{r}, col{c}
{}

template <typename T>
Matrix<T>::Matrix(size_t size) : elements{new T[size * size]}, row{size}, col{size}
{}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& matrix)
{
    row = matrix.row;
    col = matrix.col;

    elements = new T[row * col];

    for (size_t i=0; i < row * col; i++)
    {
        elements[i] = matrix.elements[i];
    }
}

template <typename T>
Matrix<T>::~Matrix()
{
    delete[] elements;
}

template <typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& matrix)
{
    this->~Matrix();

    row = matrix.row;
    col = matrix.col;

    elements = new T[row * col];

    for (size_t i=0; i < row * col; i++)
    {
        elements[i] = matrix.elements[i];
    }   

    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& matrix)
{
    if(row == matrix.row && col == matrix.col)
    {
        Matrix result(row, col);
        
        for (size_t i=0; i < row * col; i++)
        {
            result.elements[i] = elements[i] + matrix.elements[i];
        }

        return result;
    }
}

template <typename T>
Matrix<T> Matrix<T>::operator- (const Matrix<T>& matrix)
{
    if(row == matrix.row && col == matrix.col)
    {
        Matrix result(row, col);
        
        for (size_t i=0; i < row * col; i++)
        {
            result.elements[i] = elements[i] - matrix.elements[i];
        }

        return result;
    }
}

template <typename T>
Matrix<T> Matrix<T>::operator* (const Matrix<T>& matrix)
{
    if(col == matrix.row)
    {
        Matrix result(row, matrix.col);

        for(size_t i = 0; i < row; i++)
        {
            for(size_t j = 0; j < matrix.col; j++)
            {
                result.elements[i * result.col + j] = 0;
                for(size_t k = 0; k < col; k++)
                {
                    result.elements[i * result.col + j] += elements[i * col + k] * matrix.elements[k * matrix.col + j];
                }
            }
        }
        return result;
    }
}

template <typename T>
Matrix<T> Matrix<T>::transpose()
{
    Matrix result(col, row);
    
    for(size_t i=0; i < row; i++)
    {
        for(size_t j=0; j < col; j++)
        {
            result.elements[j * row + i] = elements[i * col + j]; 
        }  
    }

    return result;
}

template <typename T>
void Matrix<T>::operator*= (const Matrix<T>& matrix)
{
    *this = *this * matrix;
}

template <typename T>
void Matrix<T>::operator*= (const T& value)
{
    for(size_t i = 0; i < row * col; i++)
    {
        elements[i] *= value;
    }
}

template <typename T>
void Matrix<T>::operator/= (const T& value)
{
    for(size_t i = 0; i < row * col; i++)
    {
        elements[i] /= value;
    }
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T& value)
{
    Matrix result(row, col);
    for(size_t i = 0; i < row * col; i++)
    {
        result.elements[i] = elements[i] * value;
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(const T& value)
{
    Matrix result(row, col);
    for(size_t i = 0; i < row * col; i++)
    {
        result.elements[i] = elements[i] / value;
    }
    
    return result;
}

template <typename T>
float Matrix<T>::norm()
{
    if(col == 1 || row == 1)
    {
        T sum;
        for(size_t i = 0; i < row * col; i++)
        {
            sum += elements[i] ^ 2;
        }

        return sqrt(sum);
    }
}

template <typename T>
void Matrix<T>::set(size_t i, size_t j, T value)
{
    elements[i * col + j] = value;
}

template <typename T>
T Matrix<T>::get(size_t i, size_t j)
{
    return elements[i * col + j];
}

template <typename T>
Matrix<T>& Matrix<T>::operator<<(T x)
{
    op_counter = 0;
    if(row > 0 && col >0)
    {
        elements[0] = x;
        op_counter ++;
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::operator,(T x)
{
    if(row * col > op_counter)
    {
        elements[op_counter] = x;
        op_counter ++;
    }
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::identity()
{
    if (row == col)
    {
        for(size_t i=0; i < row; i++)
            set(i, i, 1);
    }
    
    return *this;
}

template <typename T>
Matrix<T>& Matrix<T>::zero()
{
    for(size_t i = 0; i < row * col; i++)
        elements[i] = 0;
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::minor_matrix(size_t r, size_t c)
{
    if (r <= row && c <= col)
    {
        Matrix result(row - 1, col - 1);
        int curi = 0;
        int curj = 0;
        for(size_t i = 0; i < row; i++)
        {
            if(i == r - 1)
                curi --;
            for(size_t j = 0; j < col; j++)
            {
                if (j == c - 1)
                    curj --;
                if(i != r - 1 && j != c -1)
                    result.elements[curi * (col - 1) + curj] = elements[i * col + j];

                curj ++;
            }

            curj = 0;
            curi ++;
        }

        return result;
    }
}

template <typename T>
T Matrix<T>::minor(size_t r, size_t c)
{
    if(row == col && r <= row && c <= col)
    {
        return this->minor_matrix(r, c).det();
    }
}

template <typename T>
T Matrix<T>::algcompl(size_t r, size_t c)
{
    if(row == col && r <= row && c <= col)
    {
        return (((r + c) % 2 == 0) ? 1 : -1) * this->minor(r, c);
    }
}

template <typename T>
T Matrix<T>::det()
{
    if(row == col)
    {
        if(row == 1)
        {
            return elements[0];
        }
        else if(row == 2)
        {
            return elements[0] * elements[3] - elements[1] * elements[2];
        }
        else
        {
            T det = elements[0] * algcompl(1, 1);
            for(size_t i = 2; i <= row; i++)
            {
                det += elements[(i - 1) * col] * algcompl(i, 1);
            }

            return det;
        }
    }
}

template <typename T>
Matrix<T> Matrix<T>::inv()
{
    if(row == col)
    {
        T det = this->det();
        if(det != 0)
        {
            Matrix result(row);

            for(size_t i=0; i < row; i++)
            {
                for(size_t j=0; j < col; j++)
                {
                    result.elements[i * col + j] = this->algcompl(j + 1, i + 1);
                }
            }
            result /= det;

            return result;
        }
    }
}