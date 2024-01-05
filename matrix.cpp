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

// converting vector into matrix object
matrix matrix::toMatrix(const std::vector<double> vec)
{
    // initialization of the future resulting matrix
    matrix res(vec.size(), 1);

    // filling the newly created matrix with the vector elements
    for (int i = 0; i < vec.size(); i++) {
        res.m_M[i][0] = vec[i];
    }
    return res;
}

// checks if the matrix is entirely empty
bool matrix::isEmpty()
{
    for (int i = 0; i < this->m_nRows; i++) {
        for (int j = 0; j < this->m_nCols; j++) {
            if (this->m_M[i][j] != 0.0) {return false;};
        }
    }
    return true;
}

// addition operator for our matrices
matrix matrix::operator+(const matrix& mat)
{
    // sanity check for coherent matrix dimensions, throws exception otherwise
    if ((this->m_nRows != mat.m_nRows) || (this->m_nCols != mat.m_nCols)) {
        throw std::runtime_error("Matrix dimensions must match for addition.");
    }

    // initializing the resulting matrix of the addition (homogeneity dimension-wise)
    matrix res(this->m_nRows, this->m_nCols);

    for (int i = 0; i < this->m_nRows; i++) {
        for (int j = 0; j < this->m_nCols; j++) {
            res.m_M[i][j] = this->m_M[i][j] + mat.m_M[i][j];
        }
    }
    return res;
}

// term by term multiplication (dot product) operator for our matrices
matrix matrix::operator*(const matrix& mat)
{
    // sanity check for coherent matrix dimensions, throws exception otherwise
    if ((this->m_nCols != mat.m_nRows)) {
        throw std::runtime_error("Columns number of first matrix must match rows number of second matrix.");
    }

    // initializing the resulting matrix of the dot product (adapting dimensions)
    matrix res(this->m_nRows, mat.m_nCols);

    // performing the dot product
    for (int i = 0; i < this->m_nRows; i++) {
        for (int j = 0; j < mat.m_nCols; j++) {
            double temp = 0;
            for (int h = 0; h < this->m_nCols; h++) {
                temp += this->m_M[i][h] * mat.m_M[h][j];
            }
            // visual simplification used to round to zero when arguably very close to zero
            res.m_M[i][j] = zeroRounder(temp);
        }
    }
    return res;
}

// scalar product operator for our matrices
matrix matrix::operator*(double lambda)
{
    // initializing the resulting matrix of the scalar product
    matrix res(this->m_nRows, this->m_nCols);

    for (int i = 0; i < this->m_nRows; i++) {
        for (int j = 0; j < this->m_nCols; j++) {
            res.m_M[i][j] = this->m_M[i][j] * lambda;
        }
    }
    return res;
}

// makes our matrix easier to read by simplifying close to zero terms
double matrix::zeroRounder(double value)
{
    double threshold = 1e-10;
    if ((value < threshold) && (-value < threshold)) {
        return 0.0;
    } else {
        return value;
    }
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
    matrix upperLine(1, augmented.m_nCols);

    for (int j = 0; j < augmented.m_nCols; j++) {
        upperLine.m_M[0][j] = augmented.m_M[upperIdx][j];
        augmented.m_M[upperIdx][j] = augmented.m_M[lowerIdx][j];
        augmented.m_M[lowerIdx][j] = upperLine.m_M[0][j];
    }
    return augmented;
}

// matrix determinant builder
double matrix::determinant(const matrix& mat)
{
    double det = 0;
    // straightforward determinant formula for 2 by 2 matrices
    if (mat.m_nCols == 2) {return mat.m_M[0][0] * mat.m_M[1][1] - mat.m_M[1][0] * mat.m_M[0][1];}

    for (int x = 0; x < mat.m_nCols; x++) {
        // iteratively breaks down sub-matrices into smaller matrices until size is 2 by 2: determinant is then immediate
        det += ((x % 2 == 0 ? 1 : -1) * mat.m_M[0][x] * determinant(submatrix(mat, x, 0, mat.m_nCols-1)));
    }
    return det;
}

// submatrix builder
matrix matrix::submatrix(const matrix& mat, int x, int y, int n)
{
    matrix submatrix(n, n);
    int subm_i = 0;

    // simply skips both the line and the column when the co-factor is located
    for (int i = 0; i <= n; i++) {
        int subm_j = 0;
        if (i == y) {continue;}

        for (int j = 0; j <= n; j++) {
            if (j == x) {continue;}
            submatrix.m_M[subm_i][subm_j] = mat.m_M[i][j];
            subm_j++;
        }
        subm_i++;
    }
    return submatrix;
}

// isolates a specific line from a matrix (wraps it around zeros)
matrix matrix::lineIsolator(const matrix& augmented, int lineIdx, int pivotIdx)
{
    matrix isolatedLine(augmented.m_nRows, augmented.m_nCols);
    for (int j = 0; j < augmented.m_nCols; j++) {
        isolatedLine.m_M[lineIdx][j] = augmented.m_M[pivotIdx][j];
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
        throw std::runtime_error("Inversion requires a squared matrix");
    }

    // invertibility assumed for very large matrices (computationally too expensive)
    if (this->m_nCols < 11) {
        // invertibility sanity check
        if (determinant(*this) == 0) {
            throw std::runtime_error("Determinant is zero");
        }
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

        while ((pivot < threshold) && (-pivot < threshold) && (nextLineIdx < augmented.m_nRows)) {
            augmented = lineSwapper(augmented, j, nextLineIdx);
            pivot  = augmented.m_M[j][j];
            nextLineIdx++;
        }

        // now that the pivot can be handled, we want all its below elements to be zero
        for (int i = j + 1; i < augmented.m_nRows; i++) {
            isolatedLine = lineIsolator(augmented, i, j);

            isolatedLine = lineMultiplier(isolatedLine, i, -augmented.m_M[i][j] / pivot);
            augmented = augmented + isolatedLine;
        }
    }

    // backward elimination (make the lower triangle zero to end up with a diagonal matrix)
    for (int j = augmented.m_nRows - 1; j > 0; j--) {

        // defining the pivot once again
        double pivot = augmented.m_M[j][j];
        int nextLineIdx = j - 1;

        for (int i = nextLineIdx; i >= 0; i--) {
            
            // we already treated the significance issue for pivot above so no need for redundance here
            isolatedLine = lineIsolator(augmented, i, j);

            isolatedLine = lineMultiplier(isolatedLine, i, -augmented.m_M[i][j] / pivot);
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