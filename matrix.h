#include <cstddef>
#include <vector>

class matrix
{
public:
    matrix(size_t nRows, size_t nCols);
    matrix operator+(matrix& rhs);
    matrix operator*(matrix& rhs);
    matrix lineScalar(matrix& augmented, int lineIdx, double scalar);
    matrix lineSwapper(matrix& augmented, int upperIdx, int lowerIdx);
    matrix inversion();
    double determinant(const std::vector<std::vector<double>>& matrix);
    std::vector<std::vector<double>> submatrix(const std::vector<std::vector<double>>& matrix, int x, int y, int n);

private:
    size_t m_nCols;
    size_t m_nRows;
    std::vector<std::vector<double>> m_M;
};