#include <cstddef>
#include <vector>

class matrix
{
public:
    matrix(size_t nRows, size_t nCols);
    matrix operator+(matrix& rhs);
    matrix operator*(matrix& rhs);
    int determinant(int **matrix, int n);
    int **submatrix(int **matrix, int n, int x, int y);
    matrix inversion(int **matrix, int n);

private:
    size_t m_nCols;
    size_t m_nRows;
    std::vector<std::vector<double>> m_M;
};