#include <vector>
#include <functional>
#include "matrix.h"

class pricerPDE
{
public:
    using coeffFunction = std::function<double(double, double)>;

    pricerPDE(double strike, double matu, double vol, 
              double rate, double divs, double repo,
              double multiplier, const matrix& terminalCondition,
              const matrix& boundaryConditions,
              coeffFunction a, coeffFunction b, coeffFunction c);

    double callOptionPrice();

private:
    void setupGrid();
    void setTerminalCondition(const matrix& terminalCondition);
    void setBoundaryConditions(const matrix& boundaryConditions);
    void applyCrankNicholson();

    matrix p_priceGrid;
    matrix p_timeGrid;
    matrix p_spotGrid;

    coeffFunction a, b, c;

    double p_multiplier;
    double p_strike;
    double p_matu;
    double p_vol;
    double p_rate;
    double p_divs;
    double p_repo;
    int p_m;
    int p_n;
};
