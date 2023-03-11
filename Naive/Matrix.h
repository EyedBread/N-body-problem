//
// DD1388 - Lab 2: The matrix
//

#ifndef MATRIX_H
#define MATRIX_H

#include <initializer_list>
#include <iostream>

#include <stdio.h>
#include <tgmath.h>

template <typename T>
class Matrix {

    static_assert( std::is_move_constructible_v<T>, "T must be move-constructible");
    static_assert( std::is_move_assignable_v<T>, "T must be move-assignable");

public:
    // constructors and assignment operators
    /* TODO: Make the appropriate constructor(s) explicit */
    Matrix();
    explicit Matrix(size_t dim);
    Matrix(size_t rows, size_t cols);
    Matrix(const std::initializer_list<T> & list);
    Matrix(const Matrix<T> & other);
    Matrix(Matrix<T> && other) noexcept;

    Matrix<T> & operator=(const Matrix<T> & other);
    Matrix<T> & operator=(Matrix<T> && other) noexcept;

    ~Matrix();

    // accessors
    size_t rows() const;
    size_t cols() const;

    T & operator()(size_t row, size_t col);
    const T & operator()(size_t row, size_t col) const;

    // operators
    Matrix<T> operator*(const Matrix<T> & other) const;
    Matrix<T> operator*(const size_t & other) const;
    Matrix<T> operator+(const Matrix<T> & other) const;
    Matrix<T> operator-(const Matrix<T> & other) const;
    bool operator==(const Matrix<T> & other) const;

    void operator*=(const Matrix<T> & other);
    void operator*=(const size_t & other);
    void operator+=(const Matrix<T> & other);
    void operator-=(const Matrix<T> & other);

    // methods
    void reset();

    void insert_row(size_t row);
    void append_row(size_t row);
    void remove_row(size_t row);
    void insert_column(size_t col);
    void append_column(size_t col);
    void remove_column(size_t col);

    // iterators
    typedef T* iterator;

    iterator begin();
    iterator end();

private:
    size_t m_rows;
    size_t m_cols;
    size_t m_capacity;
    T * m_vec;
};

// input/output operators
template<typename T>
std::istream & operator>>(std::istream & is, Matrix<T> & m);

template<typename T>
std::ostream & operator<<(std::ostream & os, const Matrix<T> & m);

// functions
template<typename T>
Matrix<T> identity(size_t dim);


//
// Implementations
//

template<typename T>
Matrix<T>::Matrix()
//Default constructor
 : m_rows(0), m_cols(0), m_capacity(0), m_vec(new T[0]) {
    static T def;
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            m_vec[i*m_cols + j] = def;
        }
    }
 }

template<typename T>
Matrix<T>::Matrix(size_t dim)
 : m_rows(dim), m_cols(dim), m_capacity(dim*dim), m_vec(new T[dim*dim]) {
    static T def;
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            m_vec[i*m_cols + j] = def;
        }
    }
 }

template<typename T>
Matrix<T>::Matrix(size_t rows, size_t cols)
 : m_rows(rows), m_cols(cols), m_capacity(rows*cols), m_vec(new T[rows*cols]) {
    static T def;
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            m_vec[i*m_cols + j] = def;
        }
    }
 }

template<typename T>
Matrix<T>::Matrix(const std::initializer_list<T> & list)
/*:  Initialize members here */ {
    int size = list.size();
    int sr = std::sqrt(size);
    if (sr*sr == size) {
        m_rows = sr;
        m_cols = sr;
        m_capacity = sr*sr;
        m_vec = new T[sr*sr];
        typename std::initializer_list<T>::iterator it;
        int i = 0;
        for (it = list.begin(); it!=list.end(); ++it) {
            m_vec[i] = *it;
            i++;
        }        
    }
    else {
        throw std::out_of_range("List doesn't have an even square root!"); 
    }
        
}

//Copy constructor
template<typename T>
Matrix<T>::Matrix(const Matrix<T> & other)
 : Matrix(other.m_rows, other.m_cols) {
        std::copy(other.m_vec, other.m_vec + other.m_capacity, m_vec);
    }

//Move constructor
template<typename T>
Matrix<T>::Matrix(Matrix<T> && other) noexcept
 : m_rows(other.m_rows), m_cols(other.m_cols), m_capacity(other.m_capacity), m_vec(other.m_vec) {
    other.m_rows = 0;
    other.m_cols = 0;
    other.m_capacity = 0;
    other.m_vec = nullptr;
}

template<typename T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> & other) {
    if (this == &other)
            return *this;


    m_cols = other.m_cols;
    m_rows = other.m_rows;
    m_capacity = other.m_capacity;
    T* new_m_vec = new T[m_capacity];            // allocate
    memcpy(new_m_vec, other.m_vec, m_capacity*sizeof(T)); // populate
    delete[] m_vec;                           // deallocate

    m_vec = new_m_vec;
    return *this;
}
// Move assignment
template<typename T>
Matrix<T> & Matrix<T>::operator=(Matrix<T> && other) noexcept {
    if (this == &other)
            return *this;

    delete[] m_vec; // deallocate
    m_cols = other.m_cols;
    m_rows = other.m_rows;
    m_capacity = other.m_capacity;
    m_vec = other.m_vec;

    other.m_cols = 0; //Make other hollow
    other.m_rows = 0;
    other.m_capacity = 0;
    other.m_vec = nullptr;

    return *this;
}

template<typename T>
Matrix<T>::~Matrix() {
    delete[] m_vec;
}

template<typename T>
size_t Matrix<T>::rows() const {
    return m_rows;
}

template<typename T>
size_t Matrix<T>::cols() const {
    return m_cols;
}

template<typename T>
T & Matrix<T>::operator()(size_t row, size_t col) {
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
        throw std::out_of_range("Index out of range!");
    return m_vec[row * m_cols + col];
}

template<typename T>
const T & Matrix<T>::operator()(size_t row, size_t col) const {
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
        throw std::out_of_range("Index out of range!");
    return m_vec[row * m_cols + col];
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> & other) const {
    if (m_cols != other.m_rows)
        throw std::out_of_range("Dimensions don't match!");
    Matrix<T> result(m_rows, other.m_cols);
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < other.m_cols; j++) {
            for (int k = 0; k < m_cols; k++) {
                result(i,j) += (*this)(i,k) * other(k,j);
            }
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const size_t& other) const {
    Matrix<T> result(m_rows, m_cols);
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            result(i,j) = (*this)(i,j) * other;
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> & other) const {
    if (m_rows != other.m_rows || m_cols != other.m_cols)
        throw std::out_of_range("Dimensions don't match!");
    Matrix<T> result(m_rows, m_cols);
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            result(i, j) = other(i, j) + (*this)(i, j);
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> & other) const {
    // Implementation goes here
    if (m_rows != other.m_rows || m_cols != other.m_cols)
        throw std::out_of_range("Dimensions don't match!");
    Matrix<T> result(m_rows, m_cols);
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            result(i, j) = (*this)(i, j) - other(i, j);
        }
    }
    return result;
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T> & other) const {
    if (m_rows != other.m_rows || m_cols != other.m_cols)
        return false;
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            if ((*this)(i,j) != other(i,j)) return false;
        }
    }
    return true;
}

template<typename T>
void Matrix<T>::operator*=(const Matrix<T> & other) {
    *this = *this * other;
}

template<typename T>
void Matrix<T>::operator*=(const size_t & other) {
    *this = *this * other;
}

template<typename T>
void Matrix<T>::operator+=(const Matrix<T> & other) {
    *this = *this + other;
}

template<typename T>
void Matrix<T>::operator-=(const Matrix<T> & other) {
    *this = *this - other;
}

template<typename T>
void Matrix<T>::reset() {
    // Set all elements to the default value of the type
    static T def;
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            (*this)(i,j) = def;
        }
    }
}

template<typename T>
void Matrix<T>::insert_row(size_t row) {
    static T def;
    // Implementation goes here
    if (row > m_rows - 1 || row < 0)
        throw std::out_of_range("Invalid row!");

    Matrix<T> tmp(m_rows+1,m_cols);
    for (int i = 0; i < m_cols; i++) {

        for (int j = 0; j < row; j++) {
            tmp(j,i) = (*this)(j,i);
        }  
        tmp(row, i) = def;
        for (int j = row+1; j < m_rows+1; j++) {
            tmp(j,i) = (*this)(j-1,i);
        }

    }

    std::swap(*this,tmp);
    // delete *this.begin();
}

template<typename T>
void Matrix<T>::append_row(size_t row) {
    static T def;
    if (row > m_rows - 1 || row < 0)
        throw std::out_of_range("Invalid row!");

    Matrix<T> tmp(m_rows+1,m_cols);
    for (int i = 0; i < m_cols; i++) {

        for (int j = 0; j < row+1; j++) {
            tmp(j,i) = (*this)(j,i);
        }  
        tmp(row+1, i) = def;
        for (int j = row+2; j < m_rows+1; j++) {
            tmp(j,i) = (*this)(j-1,i);
        }
    }
    std::swap(*this,tmp);
    //tmp is automatically deallocated since declared locally
}

template<typename T>
void Matrix<T>::remove_row(size_t row) {
    if (row > m_rows - 1 || row < 0)
        throw std::out_of_range("Invalid row!");

    Matrix<T> tmp(m_rows-1,m_cols);
    for (int i = 0; i < m_cols; i++) {

        for (int j = 0; j < row; j++) {
            tmp(j,i) = (*this)(j,i);
        }  
        for (int j = row; j < m_rows-1; j++) {
            tmp(j,i) = (*this)(j+1,i);
        }

    }
    std::swap(*this,tmp);
    //tmp is automatically deallocated since declared locally
}

template<typename T>
void Matrix<T>::insert_column(size_t col) {
    static T def;
    if (col > m_cols - 1 || col < 0)
        throw std::out_of_range("Invalid column!");

    Matrix<T> tmp(m_rows,m_cols+1);
    for (int i = 0; i < m_rows; i++) {

        for (int j = 0; j < col; j++) {
            tmp(i,j) = (*this)(i,j);
        }  
        tmp(i, col) = def;
        for (int j = col+1; j < m_cols+1; j++) {
            tmp(i,j) = (*this)(i,j-1);
        }

    }

    std::swap(*this,tmp);
}

template<typename T>
void Matrix<T>::append_column(size_t col) {
    static T def;
    if (col > m_cols - 1 || col < 0)
        throw std::out_of_range("Invalid column!");

    Matrix<T> tmp(m_rows,m_cols+1);
    for (int i = 0; i < m_rows; i++) {

        for (int j = 0; j < col+1; j++) {
            tmp(i,j) = (*this)(i,j);
        }  
        tmp(i, col+1) = def;
        for (int j = col+2; j < m_cols+1; j++) {
            tmp(i,j) = (*this)(i,j-1);
        }

    }

    std::swap(*this,tmp);
}

template<typename T>
void Matrix<T>::remove_column(size_t col) {
    if (col > m_cols - 1 || col < 0)
        throw std::out_of_range("Invalid column!");

    Matrix<T> tmp(m_rows,m_cols-1);
    for (int i = 0; i < m_rows; i++) {

        for (int j = 0; j < col; j++) {
            tmp(i,j) = (*this)(i,j);
        }  
        for (int j = col; j < m_cols-1; j++) {
            tmp(i,j) = (*this)(i,j+1);
        }
    }
    std::swap(*this,tmp);
    //tmp is automatically deallocated since declared locally
}

template<typename T>
typename Matrix<T>::iterator Matrix<T>::begin() {
    return m_vec;
}

template<typename T>
typename Matrix<T>::iterator Matrix<T>::end() {
    return m_vec + m_capacity - 1;
}

template<typename T>
std::istream & operator>>(std::istream & is, Matrix<T> & m) {

    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            is >> m(i,j);
        }
    }

    return is;
}

template<typename T>
std::ostream & operator<<(std::ostream & os, const Matrix<T> & m) {

    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            
            os << m(i,j);
            if (j != m.cols()-1)
                os << " ";
        }
        os << std::endl;
    }

    return os;
    
}

template<typename T>
Matrix<T> identity(size_t dim) {
    Matrix<T> m(dim);
    for (int i = 0; i < dim; i++)
        m(i,i) = 1;
    return m;
}



#endif //MATRIX_H