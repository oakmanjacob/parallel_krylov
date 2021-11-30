#pragma once
#include <parallel_krylov/gmres.h>

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * Allow complex values
 * 
 * @param in
 * @param out
 */
void gmres(GMRES_In<complex<double>> &in, GMRES_Out<complex<double>> &out);

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in 
 * @param out 
 */
void ogmres(GMRES_In<complex<double>> &in, GMRES_Out<complex<double>> &out);