#pragma once

#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <complex>

#include <spdlog/spdlog.h>
#include <parallel_krylov/matrix_csr.h>

using namespace std;

struct GMRES_In {
    MatrixCSR<complex<double>> A;
    vector<complex<double>> b;
    vector<complex<double>> x0;
    double tol;
    size_t max_it;
    size_t restart;
};

struct GMRES_Out {
    vector<complex<double>> x;
    vector<complex<double>> r;
    vector<double> r_nrm;
    size_t iter;
    bool converged;
};

/**
 * @brief Solve the system Ax = b where we consider A to be the n x n square at the top left of the provided A. 
 * This is to prevent copying.
 * 
 * @param A 
 * @param n 
 * @param b 
 * @param x 
 */
void backsub(const vector<vector<complex<double>>> &A, const size_t n, const vector<complex<double>> &b, vector<complex<double>> &x);

/**
 * @brief 
 * 
 * @param mat 
 * @param vec 
 * @param out 
 */
void matvec(MatrixCSR<complex<double>> &mat, vector<complex<double>> &vec, vector<complex<double>> &out);

/**
 * @brief Used to calculate out += Ax with an input matrix which is transposed
 * 
 * @param mat 
 * @param m 
 * @param n 
 * @param vec 
 * @param out 
 */
void matvecaddT_emplace(const vector<vector<complex<double>>> &mat, size_t m, size_t n, vector<complex<double>> &vec, vector<complex<double>> &out);

/**
 * @brief out = first - second
 * 
 * @param first 
 * @param second 
 * @param out 
 */
void vecsub(const vector<complex<double>> &first, const vector<complex<double>> &second, vector<complex<double>> &out);

/**
 * @brief Compute the operation such that
 * first[i] = first[i] + second[i] * scalar for i in [0,first.size())
 * 
 * Note: this is likely redundant with vecaddmult_emplace because we could just make
 * scalar negative.
 * 
 * @param first 
 * @param second 
 * @param scalar 
 */
void vecsubmult_emplace(vector<complex<double>> &first, const vector<complex<double>> &second, const complex<double> scalar);

/**
 * @brief Compute the operation such that
 * first[i] = first[i] + second[i] * scalar for i in [0,first.size())
 * 
 * @param first 
 * @param second 
 * @param scalar 
 */
void vecaddmult_emplace(vector<complex<double>> &first, const vector<complex<double>> &second, const complex<double> scalar);

/**
 * @brief Element-wise vector division such that
 * second[i] = first[i] * second[i] for i in [0,second.size())
 * 
 * @param first 
 * @param second 
 */
void vecdiv_emplace(const complex<double> first, vector<complex<double>> &second);

/**
 * @brief Element-wise vector multiplication such that
 * second[i] = first[i] * second[i] for i in [0,second.size())
 * 
 * @param first 
 * @param second 
 */
void vecmult_emplace(const complex<double> first, vector<complex<double>> &second);

/**
 * @brief Compute the Givens rotation matrix parameters for a and b
 * 
 * @param a 
 * @param b 
 * @param c 
 * @param s 
 */
void rotmat(const complex<double> a, const complex<double> b, complex<double> &c, complex<double> &s);

/**
 * @brief Compute the dot product of two vectors
 * 
 * @param first the first vector
 * @param second the second vector
 * @return complex<double> first dot second
 */
complex<double> dot(const vector<complex<double>> &first, const vector<complex<double>> &second);

/**
 * @brief Compute the 2-Norm of a vector
 * 
 * @param vec the vector to compute
 * @return double the resulting norm
 */
double norm(const vector<complex<double>> &vec);

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param A input sparse matrix A
 * @param b input vector b
 * @param x0 input initial guess
 * @param tol input tolerance
 * @param max_it input max number of iterations
 * @param restart input max number of iterations before restarting
 * @param x output guess for x
 * @param r output estimated residual r
 * @param r_nrm output vector containing the estimated 2-norm of the residual at each step
 * @param iter output containing the number of iterations performed
 * @param converged output true if converged
 */
void gmres(MatrixCSR<complex<double>> &A, vector<complex<double>> &b, vector<complex<double>> &x0,
           double tol, size_t max_it, size_t restart, vector<complex<double>> &x, vector<complex<double>> &r,
           vector<double> &r_nrm, size_t &iter, bool &converged);

/**
 * @brief Overloaded GMRES function simplifying interface by taking structs instead
 * of a ton of parameters
 * 
 * @param in 
 * @param out 
 */
void gmres(GMRES_In &in, GMRES_Out &out);