# Pricing vanilla options through PDEs
**<u> Group1 </u>:** *Garriga Maxime, Soughati Kenza, Saulis Lukas, Collin Thibault*, within the scope of the ```M203 C++ class```

The aim of this project is to build a pricer of vanilla option from the PDE of the derivatives payoff. Key implementations include matrix inversion algorithm and Crank-Nicholson scheme resolution. Comments are written along the way, and a detailed documentation can be found inside the repository. *(last update: Jan 5th)*


*Matrix inversion*: uses the **Gauss-Jordan** pivot algorithm to find the inverse of a matrix through its augmented form. Checks for legality of inversion *(matrix is squared and determinant is zero)*.

```cpp
matrix::inversion()
```

*Matrix determinant calculator*: uses the **Laplace expansion** algorithm to break down the n-squared matrices iteratively into 2 by 2 sub-matrices, whose determinant is easy to compute. Used to check if a matrix can be inverted.

```cpp
matrix::determinant(const matrix& mat)
```

