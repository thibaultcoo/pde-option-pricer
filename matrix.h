#include <cstddef>
#include <vector>

class matrix
{
public:
    matrix(size_t nRows, size_t nCols);
    matrix operator+(matrix& rhs);
    matrix operator*(matrix& rhs);
    matrix lineScalar(matrix& augmented, int lineIdx, double scalar);
    matrix lineIsolator(matrix& augmented, int lineIdx);
    matrix lineSwapper(matrix& augmented, int upperIdx, int lowerIdx);
    matrix inversion();
    matrix submatrix(const matrix& mat, int x, int y, int n);
    double determinant(const matrix& mat);

private:
    size_t m_nCols;
    size_t m_nRows;
    std::vector<std::vector<double>> m_M;
};