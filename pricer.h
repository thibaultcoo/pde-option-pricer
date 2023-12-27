#include <vector>
#include <functional>
#include "matrix.h"

class pricerPDE
{
public:
    using coeffFunction = std::function<double(double, double)>;

    pricerPDE(double spot, double strike, double matu, double vol, 
              double rate, double divs, double repo, double multiplier, 
              const matrix& terminalCondition, const matrix& boundaryConditions,
              coeffFunction a, coeffFunction b, coeffFunction c, coeffFunction d);

    double callOptionPrice();

private:
    void setupGrid();
    void setTerminalCondition(const matrix& terminalCondition);
    void setBoundaryConditions(const matrix& boundaryConditions);
    void applyCrankNicholson();
    double extractPrice();
    int findClosestIdx(const matrix& grid, double value);
    double interpo(int tIdx, int xIdx);

    matrix p_priceGrid, p_timeGrid, p_spotGrid;
    coeffFunction a, b, c, d;

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
