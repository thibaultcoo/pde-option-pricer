#include <vector>
#include <functional>
#include "matrix.h"

class pricerPDE
{
public:
    using coeffFunction = std::function<double(double, double)>;

    pricerPDE(double spot, double strike, double matu, double vol, 
              double rate, double divs, double repo, double multiplier, 
              const matrix& terminalCondition, 
              const matrix& lowerBoundaries, const matrix& upperBoundaries,
              coeffFunction a, coeffFunction b, coeffFunction c, coeffFunction d);

    double callOptionPrice();

private:
    void setupGrid();
    void applyCrankNicholson();
    double extractPrice();
    int findClosestIdx(const matrix& grid, double value);
    double interpo(int xIdx);

    matrix p_timeGrid, p_spotGrid, p_priceGrid;
    matrix p_terminalCondition ,p_lowerBoundaries, p_upperBoundaries;
    coeffFunction p_a, p_b, p_c, p_d;

    double p_multiplier;
    double p_spot;
    double p_strike;
    double p_matu;
    double p_vol;
    double p_rate;
    double p_divs;
    double p_repo;
    double p_dt;
    double p_dS;
    int p_m;
    int p_n;
};
