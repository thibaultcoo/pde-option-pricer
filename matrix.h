#ifndef MATRIX_H
#define MATRIX_H

#include <cstddef>
#include <vector>

class matrix
{
public:
    matrix(size_t nRows = 0, size_t nCols = 0);
    matrix operator+(const matrix& mat);
    matrix operator*(const matrix& mat);
    matrix operator*(double lambda);
    matrix lineIsolator(const matrix& augmented, int lineIdx, int pivotIdx);
    matrix lineMultiplier(matrix& augmented, int lineIdx, double scalar);
    matrix lineSwapper(matrix& augmented, int upperIdx, int lowerIdx);
    matrix inversion();
    matrix submatrix(const matrix& mat, int x, int y, int n);
    double zeroRounder(double value);
    double determinant(const matrix& mat);

    std::vector<std::vector<double>> m_M;

private:
    size_t m_nCols;
    size_t m_nRows;
};

#endif