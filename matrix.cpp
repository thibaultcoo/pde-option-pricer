#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "matrix.h"

// matrix object constructor
matrix::matrix(size_t nRows, size_t nCols)
{
    m_nCols = nCols;
    m_nRows = nRows;
    m_M.resize(m_nRows, std::vector<double>(m_nCols, 0.0));
}

// addition operator for our matrices
matrix matrix::operator+(const matrix& rhs)
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
matrix matrix::operator*(const matrix& rhs)
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

// multiplies the specific line of a matrix by a scalar
matrix matrix::lineMultiplier(matrix& augmented, int lineIdx, double scalar)
{
    for (int j = 0; j < augmented.m_nCols; j++) {
        augmented.m_M[lineIdx][j] *= scalar;
    }
    return augmented;
}

// swaps two lines of a matrix
matrix matrix::lineSwapper(matrix& augmented, int upperIdx, int lowerIdx) 
{
    matrix upperLine(0, augmented.m_nCols);

    for (int j = 0; j < augmented.m_nCols; j++) {
        upperLine.m_M[upperIdx][j] = augmented.m_M[upperIdx][j];
        augmented.m_M[upperIdx][j] = augmented.m_M[lowerIdx][j];
        augmented.m_M[lowerIdx][j] = upperLine.m_M[upperIdx][j];
    }
    return augmented;
}

// matrix determinant builder
double matrix::determinant(const matrix& mat)
{
    int det = 0;
    if (this->m_nCols == 2) {return mat.m_M[0][0] * mat.m_M[1][1] - mat.m_M[1][0] * mat.m_M[0][1];}

    for (int x = 0; x < this->m_nCols; ++x) {det += ((x % 2 == 0 ? 1 : -1) * mat.m_M[0][x] * determinant(submatrix(mat, x, 0, this->m_nCols-1)));}
    return det;
}

// submatrix builder
matrix matrix::submatrix(const matrix& mat, int x, int y, int n)
{
    matrix submatrix(n, n);
    int subm_i = 0;

    for (int i = 0; i < n; i++) {
        int subm_j = 0;
        if (i == y) {continue;}

        for (int j = 0; j < n; j++) {
            if (j == x) {continue;}
            submatrix.m_M[subm_i][subm_j] = mat.m_M[i][j];
            subm_j++;
        }
        subm_i++;
    }
    return submatrix;
}

matrix matrix::lineIsolator(const matrix& augmented, int lineIdx)
{
    matrix isolatedLine(augmented.m_nRows, augmented.m_nCols);
    for (int j = 0; j < augmented.m_nCols; j++) {
        isolatedLine.m_M[lineIdx][j] = augmented.m_M[lineIdx][j];
    }

    return isolatedLine;
}

// inversion of a given squared matrix
matrix matrix::inversion()
{
    // augmented matrix initialization
    matrix augmented(this->m_nRows, 2*this->m_nCols);

    // inversion resulting matrix initialization
    matrix inverted(this->m_nRows, this->m_nCols);

    // zero-filled matrix with the exception of a unique line
    matrix isolatedLine(augmented.m_nRows, augmented.m_nCols);

    // squared sanity check
    if (this->m_nRows != this->m_nCols) {
        throw std::runtime_error("Inversion requires a squared matrix.");
    }

    // invertibility sanity check
    if (determinant(*this) != 0) {
        throw std::runtime_error("Determinant is not zero");
    }

    // creation of the augmented matrix
    for (int i = 0; i < this->m_nRows; i++) {
        for (int j = 0; j < this->m_nCols; j++) {
            augmented.m_M[i][j] = this->m_M[i][j];
        }
        augmented.m_M[i][i + this->m_nCols] = 1;
    }

    // threshold definition for gauss elimination pivot handling
    double threshold = 1e-4;

    // forward elimination (to construct an upper triangular matrix on the left side)
    for (int j = 0; j < augmented.m_nRows - 1; j++) {

        // need to check if the pivot element is significantly non-zero
        double pivot = augmented.m_M[j][j];
        int nextLineIdx = j + 1;

        while (pivot < threshold) {
            augmented = lineSwapper(augmented, j, nextLineIdx);
            pivot  = augmented.m_M[j][j];
            nextLineIdx++;
        }

        // now that the pivot can be handled, we want all its below elements to be zero
        for (int i = j + 1; i < augmented.m_nRows; i++) {
            isolatedLine = lineIsolator(augmented, i);

            isolatedLine = lineMultiplier(isolatedLine, i, -1 * pivot / isolatedLine.m_M[i][j]);
            augmented = augmented + isolatedLine;
        }
    }

    // check for the significance of the final diagonal element
    while (augmented.m_M[augmented.m_nRows][augmented.m_nRows] < threshold) {
        for (int j = 0; j < augmented.m_nCols; j++)
            augmented.m_M[augmented.m_nRows][j] *= 10;
    }

    // backward elimination (make the lower triangle zero to end up with a diagonal matrix)
    for (int j = augmented.m_nRows; j > 1; j--) {

        // defining the pivot once again
        double pivot = augmented.m_M[j][j];

        for (int i = augmented.m_nRows - 1; i > 0; i--) {
            
            // we already treated the significance issue for pivot above so no need for redundance here
            isolatedLine = lineIsolator(augmented, i);

            isolatedLine = lineMultiplier(isolatedLine, i, -1 * pivot / isolatedLine.m_M[i][j]);
            augmented = augmented + isolatedLine;
        }
    }

    // normalization of diagonal elements (divide all elements by their scalar)
    for (int i = 0; i < augmented.m_nRows; i++) {
        double diag = augmented.m_M[i][i];

        for (int j = 0; j < augmented.m_nCols; j++) {
            augmented.m_M[i][j] /= diag;
        }
    }

    // isolation of the inversion resulting matrix
    for (int i = 0; i < inverted.m_nRows; i++) {
        for (int j = 0; j < inverted.m_nCols; j++) {
            inverted.m_M[i][j] = augmented.m_M[i][j + this->m_nCols];
        }
    }
    return inverted;
}