# Pricing vanilla options through PDEs
**<u> Group1 </u>:** *Garriga Maxime, Soughati Kenza, Saulis Lukas, Collin Thibault*, within the scope of the ```M203 C++ class```

The aim of this project is to build a pricer of vanilla option from the PDE of the derivatives payoff. Key implementations include matrix inversion algorithm and Crank-Nicholson scheme resolution. Comments are written along the way, and a detailed documentation can be found inside the repository. *(last update: Jan 7th)*


*Matrix inversion*: uses the **Gauss-Jordan** pivot algorithm to find the inverse of a matrix through its augmented form. Checks for legality of inversion *(matrix is squared and determinant is zero)*.

```cpp
matrix::inversion()
```

*Matrix determinant calculator*: uses the **Laplace expansion** algorithm to break down the n-squared matrices iteratively into 2 by 2 sub-matrices, whose determinant is easy to compute. Used to check if a matrix can be inverted.

```cpp
matrix::determinant(const matrix& mat)
```

*Theoretical call option pricing*: uses **Black-Scholes** formula with dividends and repo to compute the price of a vanilla call option. Used to study the convergence of a price obtained through finite difference.

```cpp
blackScholes::callOptionPrice()
```

*Finite difference grid initialization*: prepares the time and spot grids for our finite difference resolution. Also prepares the boundary and terminal conditions, default values are used if none are provided by the user.

```cpp
pricerPDE::setupGrid()
```

*Finite difference call option pricing*: uses **Crank-Nicholson** scheme to build the whole final prices grid backwards in time, and an extraction function to isolate the price corresponding to our input level of spot.

```cpp
pricerPDE::callOptionPrice()
```

The scope of changeable lines is labelled with ```Prof.``` and encompasses the resolution parameters, option parameters and boundary/terminal conditions. All code is ours except when stated otherwise. Thanks for this semester, we all wish you a happy New Year !
