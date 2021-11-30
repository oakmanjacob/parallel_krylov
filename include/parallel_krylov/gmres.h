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

/**
 * @brief Struct for holding the imput values to GMRES implementations
 * 
 * @tparam V the datatype of the values in A, b, and x0
 * @param A input sparse matrix A
 * @param b input vector b
 * @param x0 input initial guess
 * @param tol input tolerance
 * @param max_it input max number of iterations
 * @param restart input max number of iterations before restarting
 */
template <typename V> struct GMRES_In {
    MatrixCSR<V> A;
    vector<V> b;
    vector<V> x0;
    double tol;
    size_t max_it;
    size_t restart;
};

/**
 * @brief Struct for holding the output values to GMRES implementations
 * 
 * @tparam V the datatype of the values in x and r
 * @param x output guess for x
 * @param r output estimated residual r
 * @param r_nrm output vector containing the estimated 2-norm of the residual at each step
 * @param iter output containing the number of iterations performed
 * @param converged output true if converged
 */
template <typename V> struct GMRES_Out {
    vector<V> x;
    vector<V> r;
    vector<double> r_nrm;
    size_t iter = 0;
    bool converged = false;
};

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning or optimizations
 * 
 * @param in 
 * @param out 
 */
void gmres(GMRES_In<double> &in, GMRES_Out<double> &out);

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in 
 * @param out 
 */
void ogmres(GMRES_In<double> &in, GMRES_Out<double> &out);

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in 
 * @param out 
 */
void ogmres_simd(GMRES_In<double> &in, GMRES_Out<double> &out);

/**
 * @brief Perform restarted gmres in parallel to solve the equation Ax = b without preconditioning
 * 
 * @param in 
 * @param out 
 */
void pgmres_sync(GMRES_In<double> &in, GMRES_Out<double> &out, size_t thread_count);

/**
 * @brief Perform restarted gmres in parallel to solve the equation Ax = b without preconditioning
 * 
 * @param in 
 * @param out 
 */
// void pgmres_async(GMRES_In<double> &in, GMRES_Out<double> &out, size_t threads);