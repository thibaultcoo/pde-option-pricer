# Pricing vanilla options through PDEs
**<u> Group1 </u>:** *Garriga Maxime, Soughati Kenza, Saulis Lukas, Collin Thibault*, within the scope of the M203 C++ class

The goal of this project is to build a pricer of vanilla option from the PDE of the derivatives payoff. Key implementations include matrix inversion algorithm and boundaries condition verification. Comments are written along the way, and a detailed documentation will follow soon.

### Creation of a matrix class
*Matrix construction*: builds a matrix element given its dimensions. Provides a representation.

```cpp
matrix::matrix(size_t nRows, size_t nCols)
```

*Matrix inversion*: uses the **Gauss-Jordan** pivot algorithm to find the inverse of a matrix through its augmented form. Checks for legality of inversion *(matrix is squared and determinant is zero)*.

```cpp
matrix::inversion()
```

*Matrix determinant calculator*: uses the **Laplace expansion** algorithm to break down the n-squared matrices iteratively into 2 by 2 sub-matrices, whose determinant is easy to compute.

```cpp
matrix::determinant(const matrix& mat)
```

*Submatrix calculator*: breaks down an initial n-squared matrix into a reduced (n-1)-squared form, removing the line and column where a given element is located.

```cpp
matrix::submatrix(const matrix& mat, int x, int y, int n)
```

*Matrix dot multiplier*: replacing the built-in multiplier operator to make it adapted for matrix dot product.

```cpp
matrix::operator*(const matrix& rhs)
```

*Matrix addition*: replacing the built-in addition operator to make it adapted for matrix dot product.

```cpp
matrix::operator+(const matrix& rhs)
```
