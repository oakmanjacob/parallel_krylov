#include <vector>
#include <functional>

#include <parallel_krylov/matrix_csr.h>

/**
 * @brief This is a class which is used to generate
 * a sparse matrix and right hand side from the
 * 2-dimensional convection-diffusion equation
 * 
 * (p*u_x)_x - (q*u_y)_y + r*u_x + s*u_y + t*u = f
 * 
 * using 2nd order central finite differences
 * discretization. Here p and q are function of x,y,
 * and u and therefore the convection diffusion
 * equation above in nonlinear in diffusion.
 * 
 * This was translated from matlab code written by
 * Arielle Carr
 * 
 * @author Jacob Oakman
 */
class GeneratorCD2D {
    /**
     * @brief Number of grid points in x-direction
     * not counting boundary points for Dir bc
     * counting boundary points for Neumann/Robin bc
     */
    size_t nx;

    /**
     * @brief Number of grid points in y direction
     * not counting boundary points for Dir bc
     * counting boundary points for Neumann/Robin bc
     */
    size_t ny;

    /**
     * @brief Fixed mesh width in x-direction; depends on bc
     */
    size_t dx;

    /**
     * @brief idem
     */
    size_t dy;

    size_t x0;
    size_t y0;

    /**
     * @brief Solution at the previous Newton iteration,
     * stored as a vector
     * 
     */
    std::vector<double> u;

    /**
     * @brief diffusion coefficient (positive) in x direction
     */
    std::function<double(size_t x, size_t y, std::vector<double> &u)> p;
    
    /**
     * @brief idem y direction
     */
    std::function<double(size_t x, size_t y, std::vector<double> &u)> q;
    
    /**
     * @brief derivative of p
     */
    std::function<double(size_t x, size_t y, std::vector<double> &u)> pdx;
    
    /**
     * @brief derivative of q
     */
    std::function<double(size_t x, size_t y, std::vector<double> &u)> pqx;


    /**
     * @brief convection in x direction (flow velocity)
     */
    std::function<double(size_t x, size_t y)> r;

    /**
     * @brief idem y direction
     */
    std::function<double(size_t x, size_t y)> s;

    /**
     * @brief reaction or absorption rate
     */
    std::function<double(size_t x, size_t y)> t;

    /**
     * @brief forcing function (rhs); source or sink
     */
    std::function<double(size_t x, size_t y)> f;

    // Output fields

    std::vector<double> F;

    MatrixCSR<double> J;

    MatrixCSR<double> Fmat;

    std::vector<double> rhsF;


    GeneratorCD2D();

    ~GeneratorCD2D() {}
};