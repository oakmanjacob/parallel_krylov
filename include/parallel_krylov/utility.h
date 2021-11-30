#pragma once

#include <vector>
#include <string>

#include <parallel_krylov/matrix_csr.h>

using namespace std;

/**
 * @brief Solve the system Ax = b where we consider A to be the n x n square at the top left of the provided A. 
 * This is to prevent copying.
 * 
 * @param A 
 * @param n 
 * @param b 
 * @param x 
 */
void backsub(const vector<vector<double>> &A, const size_t n, const vector<double> &b, vector<double> &x);

/**
 * @brief Compute a sparse matvec operation using a CSR matrix and a dense vector
 * 
 * out = mat * vec
 * 
 * @param mat sparse matrix
 * @param vec dense vector
 * @param out output vector
 */
void matvec(MatrixCSR<double> const &mat, const vector<double> &vec, vector<double> &out);

/**
 * @brief Used to calculate out += Ax with an input matrix which is transposed
 * 
 * @param mat 
 * @param m 
 * @param n 
 * @param vec 
 * @param out 
 */
void matvecadd_emplace(const vector<vector<double>> &mat, size_t m, size_t n, vector<double> &vec, vector<double> &out);

/**
 * @brief Used to calculate out += Ax with an input matrix which is transposed
 * 
 * @param mat 
 * @param m 
 * @param n 
 * @param vec 
 * @param out 
 */
void matvecaddT_emplace(const vector<vector<double>> &mat, size_t m, size_t n, vector<double> &vec, vector<double> &out);

/**
 * @brief out = first - second
 * 
 * @param first 
 * @param second 
 * @param out 
 */
void vecsub(const vector<double> &first, const vector<double> &second, vector<double> &out);

/**
 * @brief Compute the operation such that
 * first[i] = first[i] + second[i] * scalar for i in [0,first.size())
 * 
 * @param first 
 * @param second 
 * @param scalar 
 */
void vecaddmult_emplace(vector<double> &first, const vector<double> &second, const double scalar);

/**
 * @brief Compute the operation such that
 * first[i] = first[i] + second[i] * scalar for i in [0,first.size())
 * 
 * @param first 
 * @param second 
 * @param scalar 
 */
void vecaddmult_emplace_simd(vector<double> &first, const vector<double> &second, const double scalar);

/**
 * @brief Element-wise vector division such that
 * second[i] = first[i] * second[i] for i in [0,second.size())
 * 
 * @param first 
 * @param second 
 */
void vecdiv_emplace(const double first, vector<double> &second);

/**
 * @brief Element-wise vector multiplication such that
 * second[i] = first[i] * second[i] for i in [0,second.size())
 * 
 * @param first 
 * @param second 
 */
void vecmult_emplace(const double first, vector<double> &second);

/**
 * @brief Compute the Givens rotation matrix parameters for a and b
 * 
 * @param a 
 * @param b 
 * @param c 
 * @param s 
 */
void rotmat(const double a, const double b, double &c, double &s);

/**
 * @brief Compute the dot product of two vectors
 * 
 * @param first the first vector
 * @param second the second vector
 * @return double first dot second
 */
double dot(const vector<double> &first, const vector<double> &second);

/**
 * @brief Compute the dot product of two vectors taking advantage of AVX2 SIMD instructions
 * 
 * @param first the first vector
 * @param second the second vector
 * @return double first dot second
 */
double dot_simd(const vector<double> &first, const vector<double> &second);

/**
 * @brief Compute the 2-Norm of a vector
 * 
 * @param vec the vector to compute
 * @return double the resulting norm
 */
double norm(const vector<double> &vec);