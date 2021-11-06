#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <complex>

#include <spdlog/spdlog.h>
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
 * @brief 
 * 
 * @param mat 
 * @param vec 
 * @param out 
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
 * Note: this is likely redundant with vecaddmult_emplace because we could just make
 * scalar negative.
 * 
 * @param first 
 * @param second 
 * @param scalar 
 */
void vecsubmult_emplace(vector<complex<double>> &first, const vector<complex<double>> &second, const complex<double> scalar) {
    assert(first.size() == second.size());

    for (size_t i = 0; i < first.size(); i++) {
        first[i] -= second[i] * scalar;
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
           vector<double> &r_nrm, size_t &iter, bool &converged) {
    assert(A.get_row_count() > 0 && A.get_col_count() > 0);
    assert(A.get_col_count() == x0.size());
    assert(A.get_row_count() == b.size());

    x = x0;
    r_nrm.clear();
    r_nrm.resize(max_it + 1);
    converged = false;

    // Compute the initial value for Ax
    vector<complex<double>> Ax;
    matvec(A, x0, Ax);
    
    // Compute the initial residual
    vecsub(b, Ax, r);
    r_nrm[0] = norm(r);

    if (r_nrm[0] <= tol) {
        return;
    }

    // Normalize Tolerance
    tol = norm(b)*tol;

    vector<vector<complex<double>>> V(restart + 1);
    vector<vector<complex<double>>> H(restart + 1);
    for (size_t i = 0; i < restart + 1; i++) {
        H[i].resize(restart);
    }

    vector<complex<double>> cs(restart);
    vector<complex<double>> sn(restart);

    vector<complex<double>> e1(restart + 1);
    e1[0] = 1.0;

    iter = 0;
    bool notconv = true;
    while (notconv && iter < max_it) {
        // Arnoldi (Construct orthonormal basis)
        int64_t i = 0;
        V[0] = r;
        vecdiv_emplace(r_nrm[iter], V[0]);

        vector<complex<double>> s = e1;
        vecmult_emplace(r_nrm[iter], s);
        while (notconv && iter < max_it && i < (int64_t)restart) {
            iter++;
            matvec(A, V[i], V[i + 1]);
            for (int64_t k = 0; k <= i; k++) {
                H[k][i] = dot(V[k], V[i + 1]);

                // w = w - H(k,i)*V(:,k);
                vecsubmult_emplace(V[i + 1], V[k], H[k][i]);
            }

            H[i+1][i] = norm(V[i + 1]);
            if (real(H[i+1][i]) > 0) {
                vecdiv_emplace(H[i+1][i], V[i + 1]);
                
                // Apply Givens rotation
                for (int64_t k = 0; k < i; k++) {
                    auto temp = cs[k]*H[k][i] + conj(sn[k])*H[k + 1][i];
                    H[k + 1][i] = -sn[k]*H[k][i] + conj(cs[k])*H[k + 1][i];
                    H[k][i] = temp;
                }
                // Form the i-th Givens rotation matrix
                rotmat(H[i][i], H[i+1][i], cs[i], sn[i]);

                // Approximate residual norm
                auto temp = cs[i]*s[i];
                s[i + 1] = -sn[i]*s[i];
                s[i] = temp;

                H[i][i] = cs[i]*H[i][i] + conj(sn[i])*H[i + 1][i];
                H[i + 1][i] = 0.0;
                r_nrm[iter] = abs(s[i + 1]);

                if (r_nrm[iter] <= tol) {
                    notconv = false;
                }
                else {
                    i++;
                }
            }
            else {
                notconv = false;
            }
        }

        // Form the (approximate) solution
        if (!notconv) {
            vector<complex<double>> y;
            backsub(H, i+1, s, y);
            // x = x + V(:,1:i)*y;
            matvecaddT_emplace(V, x0.size(), i+1, y, x);
        }
        else {
            vector<complex<double>> y;
            backsub(H, restart, s, y);
            matvecaddT_emplace(V, x0.size(), restart, y, x);
            i--;
        }

        // Compute new Ax
        matvec(A, x, Ax);

        // Compute residual
        vecsub(b, Ax, r);

        // Compute norm of residual
        r_nrm[iter] = norm(r);

        if (r_nrm[iter] <= tol) {
            notconv = false;
        }
        else {
            if (!notconv) {
                spdlog::warn("Encountered false convergence at iteration {}", iter);
            }
            notconv = true;
        }
    }

    // Did we converge?
    if (r_nrm[iter] <= tol) {
        converged = true;
    }

    iter++;

    // Eleminate excess size in residual
    r_nrm.resize(iter);
}

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

void gmres(GMRES_In &in, GMRES_Out &out) {
    gmres(in.A, in.b, in.x0, in.tol, in.max_it, in.restart, out.x, out.r, out.r_nrm, out.iter, out.converged);
}