#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <matrix.h>

// matrix object constructor
matrix::matrix(size_t nRows, size_t nCols)
{
    m_nCols = nCols;
    m_nRows = nRows;
    m_M.resize(m_nRows, std::vector<double>(m_nCols, 0.0));
}

// addition operator for our matrices
matrix matrix::operator+(matrix& rhs)
{
    // sanity check for coherent matrix dimensions, throws exception otherwise
    if ((this->m_nRows != rhs.m_nRows) || (this->m_nCols != rhs.m_nCols)) {
        throw std::runtime_error("Matrix dimensions must match for addition.");
    }

    // initializing the resulting matrix of the addition (homogeneity dimension-wise)
    matrix res(this->m_nRows, this->m_nCols);

    for (int i = 0; i < this->m_nRows; i++) {
        for (int j = 0; j < this->m_nCols; j++) {
            res.m_M[i][j] = this->m_M[i][j] + rhs.m_M[i][j];
        }
    }
    return res;
}

// term by term multiplication (dot product) operator for our matrices
matrix matrix::operator*(matrix& rhs)
{
    // sanity check for coherent matrix dimensions, throws exception otherwise
    if ((this->m_nCols != rhs.m_nRows)) {
        throw std::runtime_error("Columns number of first matrix must match rows number of second matrix.");
    }

    // initializing the resulting matrix of the dot product (adapting dimensions)
    matrix res(this->m_nRows, rhs.m_nCols);

    for (int i = 0; i < this->m_nRows; i++) {
        for (int j = 0; j < rhs.m_nCols; j++) {
            double temp = 0;
            for (int h = 0; h < this->m_nCols; h++) {
                temp += this->m_M[i][h] * rhs.m_M[h][j];
            }
            res.m_M[i][j] = temp;
        }
    }
    return res;
}

// matrix determinant builder
double matrix::determinant(const std::vector<std::vector<double>>& matrix)
{
    int det = 0;
    if (this->m_nCols == 2) {return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];}

    for (int x = 0; x < this->m_nCols; ++x) {det += ((x % 2 == 0 ? 1 : -1) * matrix[0][x] * determinant(submatrix(matrix, x, 0, this->m_nCols-1)));}
    return det;
}

// submatrix builder
std::vector<std::vector<double>> submatrix(const std::vector<std::vector<double>>& matrix, int x, int y, int n)
{
    std::vector<std::vector<double>> submatrix;
    int subm_i = 0;

    for (int i = 0; i < n; i++) {
        int subm_j = 0;
        if (i == y) {continue;}

        for (int j = 0; j < n; j++) {
            if (j == x) {continue;}
            submatrix[subm_i][subm_j] = matrix[i][j];
            subm_j++;
        }
        subm_i++;
    }
    return submatrix;
}

// inversion of a given squared matrix
matrix matrix::inversion()
{
    // augmented matrix initialization
    matrix augmented(this->m_nRows, 2*this->m_nCols);

    // squared sanity check
    if (this->m_nRows != this->m_nCols) {
        throw std::runtime_error("Inversion requires a squared matrix.");
    }

    // invertibility sanity check
    if (determinant(this->m_M) != 0) {
        throw std::runtime_error("Determinant is not zero");
    }

    // creation of the augmented matrix
    for (int i = 0; i < this->m_nCols; i++) {
        for (int j = 0; j < this->m_nCols; j++) {
            augmented.m_M[i][j] = this->m_M[i][j];
        }
        augmented.m_M[i][i + this->m_nCols] = 1;
    }
}