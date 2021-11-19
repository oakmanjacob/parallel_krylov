#pragma once

#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <complex>

#include <spdlog/spdlog.h>
#include <parallel_krylov/matrix_csr.h>
#include <parallel_krylov/utility.h>

using namespace std;

template <typename V> struct GMRES_In {
    MatrixCSR<V> A;
    vector<V> b;
    vector<V> x0;
    double tol;
    size_t max_it;
    size_t restart;
};

template <typename V> struct GMRES_Out {
    vector<V> x;
    vector<V> r;
    vector<double> r_nrm;
    size_t iter;
    bool converged;
};

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning or optimizations
 * 
 * @param in 
 * @param out 
 */
void sgmres_old(GMRES_In<complex<double>> &in, GMRES_Out<complex<double>> &out);

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning or optimizations
 * 
 * @param in 
 * @param out 
 */
void sgmres_old(GMRES_In<double> &in, GMRES_Out<double> &out);

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in 
 * @param out 
 */
void sgmres_new(GMRES_In<complex<double>> &in, GMRES_Out<complex<double>> &out);

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in 
 * @param out 
 */
void sgmres_new(GMRES_In<double> &in, GMRES_Out<double> &out);