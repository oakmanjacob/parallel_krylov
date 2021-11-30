#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <complex>

#include <spdlog/spdlog.h>
#include <parallel_krylov/matrix_csr.h>
#include <parallel_krylov/utility.h>
#include <parallel_krylov/utility_complex.h>

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
void backsub(const vector<vector<complex<double>>> &A, const size_t n, const vector<complex<double>> &b, vector<complex<double>> &x) {
    assert(n > 0);
    assert(A.size() >= n);
    assert(A[0].size() >= n);
    assert(b.size() >= n);

    x.resize(n);
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];

    for (int64_t i = n - 2; i >= 0; i--) {
        complex<double> sum = 0;
        for (size_t j = i + 1; j < n; j++) {
            sum += A[i][j]*x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
}

/**
 * @brief Compute a sparse matvec operation using a CSR matrix and a dense vector
 * 
 * out = mat * vec
 * 
 * @param mat sparse matrix
 * @param vec dense vector
 * @param out output vector
 */
void matvec(MatrixCSR<complex<double>> &mat, vector<complex<double>> &vec, vector<complex<double>> &out) {
    assert(vec.size() == mat.get_col_count());

    out.clear();
    out.resize(mat.get_row_count());
    size_t row_index = 0;
    for (size_t i = 0; i < mat._values.size(); i++) {
        while (row_index + 1 < mat._row_pointers.size() && mat._row_pointers[row_index + 1] <= i) {
            row_index++;
        }

        out[row_index] += mat._values[i] * vec[mat._col_indexes[i]];
    }
}

/**
 * @brief Used to calculate out += Ax with an input matrix which is transposed
 * 
 * @param mat 
 * @param m 
 * @param n 
 * @param vec 
 * @param out 
 */
void matvecadd_emplace(const vector<vector<complex<double>>> &mat, size_t m, size_t n, vector<complex<double>> &vec, vector<complex<double>> &out) {
    assert(m > 0 && n > 0);
    assert(mat.size() >= m);
    assert(vec.size() >= n);
    assert(mat[0].size() >= n);

    out.resize(m);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            out[i] += mat[i][j] * vec[j];
        }
    }
}

/**
 * @brief Used to calculate out += Ax with an input matrix which is transposed
 * 
 * @param mat 
 * @param m 
 * @param n 
 * @param vec 
 * @param out 
 */
void matvecaddT_emplace(const vector<vector<complex<double>>> &mat, size_t m, size_t n, vector<complex<double>> &vec, vector<complex<double>> &out) {
    assert(m > 0 && n > 0);
    assert(mat.size() >= n);
    assert(vec.size() >= n);

    size_t new_n = n;
    while (new_n > 0 && mat[new_n].size() < m) {
        new_n--;
    }

    if (new_n != n) {
        spdlog::warn("matvecaddT: Averted seg fault by reducing n {} m {}, new n {}", n, m, new_n);
        n = new_n;
    }

    out.resize(m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            out[j] += mat[i][j] * vec[i];
        }
    }
}

/**
 * @brief out = first - second
 * 
 * @param first 
 * @param second 
 * @param out 
 */
void vecsub(const vector<complex<double>> &first, const vector<complex<double>> &second, vector<complex<double>> &out) {
    assert(first.size() == second.size());
    out = first;

    for (size_t i = 0; i < out.size(); i++) {
        out[i] -= second[i];
    }
}

/**
 * @brief Compute the operation such that
 * first[i] = first[i] + second[i] * scalar for i in [0,first.size())
 * 
 * @param first 
 * @param second 
 * @param scalar 
 */
void vecaddmult_emplace(vector<complex<double>> &first, const vector<complex<double>> &second, const complex<double> scalar) {
    assert(first.size() == second.size());

    for (size_t i = 0; i < first.size(); i++) {
        first[i] += second[i] * scalar;
    }
}

/**
 * @brief Element-wise vector division such that
 * second[i] = first[i] * second[i] for i in [0,second.size())
 * 
 * @param first 
 * @param second 
 */
void vecdiv_emplace(const complex<double> first, vector<complex<double>> &second) {
    assert(first != 0.0);

    for (size_t i = 0; i < second.size(); i++) {
        second[i] /= first;
    }
}

/**
 * @brief Element-wise vector multiplication such that
 * second[i] = first[i] * second[i] for i in [0,second.size())
 * 
 * @param first 
 * @param second 
 */
void vecmult_emplace(const complex<double> first, vector<complex<double>> &second) {
    for (size_t i = 0; i < second.size(); i++) {
        second[i] *= first;
    }
}

/**
 * @brief Compute the Givens rotation matrix parameters for a and b
 * 
 * @param a 
 * @param b 
 * @param c 
 * @param s 
 */
void rotmat(const complex<double> a, const complex<double> b, complex<double> &c, complex<double> &s) {
    assert(b == 0.0 || abs(a) != 0);

    if (b == 0.0) {
        c = 1.0;
        s = 0.0;
    }
    else if ( abs(b) > abs(a)) {
        auto temp = abs(a) / abs(b);
        auto tempsq = temp * temp;
        c = sqrt(tempsq / (tempsq + 1));
        s = conj(a/(temp*b)) / sqrt(1 + tempsq);
    }
    else {
        auto temp = abs(b) / abs(a);
        auto tempsq = temp * temp;
        c = 1.0 / sqrt(1 + tempsq);
        s = sqrt(tempsq/(tempsq+1)) * (b/(a*temp));
    }
}

/**
 * @brief Compute the dot product of two vectors
 * 
 * @param first the first vector
 * @param second the second vector
 * @return complex<double> first dot second
 */
complex<double> dot(const vector<complex<double>> &first, const vector<complex<double>> &second) {
    assert(first.size() == second.size());
    complex<double> sum = 0;

    for (size_t i = 0; i < first.size(); i++) {
        sum += first[i] * second[i];
    }

    return sum;
}

/**
 * @brief Compute the 2-Norm of a vector
 * 
 * @param vec the vector to compute
 * @return double the resulting norm
 */
double norm(const vector<complex<double>> &vec) {
    return sqrt(real(dot(vec, vec)));
}