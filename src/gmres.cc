#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <complex>

#include <spdlog/spdlog.h>
#include <parallel_krylov/matrix_csr.h>
#include <parallel_krylov/gmres.h>

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning or optimizations
 * Allow complex values
 * 
 * @param in &GMRES_In<complex<double>> input struct
 * @param out &GMRES_Out<complex<double>> output struct 
 */
void sgmres_old(GMRES_In<complex<double>> &in, GMRES_Out<complex<double>> &out) {
    assert(in.A.get_row_count() > 0 && in.A.get_col_count() > 0);
    assert(in.A.get_col_count() == in.x0.size());
    assert(in.A.get_row_count() == in.b.size());

    out.x = in.x0;
    out.r_nrm.clear();
    out.r_nrm.resize(in.max_it + 1);
    out.converged = false;

    // Compute the initial value for Ax
    vector<complex<double>> Ax;
    matvec(in.A, in.x0, Ax);
    
    // Compute the initial residual
    vecsub(in.b, Ax, out.r);
    out.r_nrm[0] = norm(out.r);

    if (out.r_nrm[0] <= in.tol) {
        return;
    }

    // Normalize Tolerance
    double rel_tol = norm(in.b)*in.tol;

    vector<vector<complex<double>>> V(out.x.size());
    for (size_t i = 0; i < V.size(); i++) {
        V[i].resize(in.restart + 1);
    }

    vector<complex<double>> column(out.x.size());

    vector<vector<complex<double>>> H(in.restart + 1);
    for (size_t i = 0; i < H.size(); i++) {
        H[i].resize(in.restart);
    }

    vector<complex<double>> cs(in.restart);
    vector<complex<double>> sn(in.restart);

    vector<complex<double>> e1(in.restart + 1);
    e1[0] = 1.0;

    out.iter = 0;
    bool notconv = true;
    while (notconv && out.iter < in.max_it) {
        // Arnoldi (Construct orthonormal basis)
        int64_t i = 0;

        for (size_t j = 0; j < out.r.size(); j++) {
            V[j][0] = out.r[j] / out.r_nrm[out.iter];
        }

        vector<complex<double>> s = e1;
        vecmult_emplace(out.r_nrm[out.iter], s);
        while (notconv && out.iter < in.max_it && i < (int64_t)in.restart) {
            out.iter++;
            vector<complex<double>> w(in.restart + 1);
            for (size_t j = 0; j < V.size(); j++) {
                column[j] = V[j][i];
            }
            matvec(in.A, column, w);
            for (int64_t k = 0; k <= i; k++) {
                for (size_t j = 0; j < w.size(); j++) {
                    H[k][i] += V[j][k] * w[j];
                }

                // w = w - H(k,i)*V(:,k);
                for (size_t j = 0; j < w.size(); j++) {
                    w[j] -= H[k][i] * V[j][k];
                }
            }

            H[i+1][i] = norm(w);
            if (real(H[i+1][i]) > 0) {
                for (size_t j = 0; j < w.size(); j++) {
                    V[j][i+1] = w[j] / H[i+1][i];
                }
                
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
                out.r_nrm[out.iter] = abs(s[i + 1]);

                if (out.r_nrm[out.iter] <= rel_tol) {
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

            matvecadd_emplace(V, in.x0.size(), i+1, y, out.x);
        }
        else {
            vector<complex<double>> y;
            backsub(H, in.restart, s, y);

            matvecadd_emplace(V, in.x0.size(), in.restart, y, out.x);
            i--;
        }

        // Compute new Ax
        matvec(in.A, out.x, Ax);

        // Compute residual
        vecsub(in.b, Ax, out.r);

        // Compute norm of residual
        out.r_nrm[out.iter] = norm(out.r);

        if (out.r_nrm[out.iter] <= rel_tol) {
            notconv = false;
        }
        else {
            if (!notconv) {
                spdlog::warn("Encountered false convergence at iteration {}", out.iter);
            }
            notconv = true;
        }
    }

    // Did we converge?
    if (out.r_nrm[out.iter] <= rel_tol) {
        out.converged = true;
    }

    out.iter++;

    // Eleminate excess size in residual
    out.r_nrm.resize(out.iter);
}

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning or optimizations
 * 
 * @param in &GMRES_In<double> input struct
 * @param out &GMRES_Out<double> output struct 
 */
void sgmres_old(GMRES_In<double> &in, GMRES_Out<double> &out) {
    assert(in.A.get_row_count() > 0 && in.A.get_col_count() > 0);
    assert(in.A.get_col_count() == in.x0.size());
    assert(in.A.get_row_count() == in.b.size());

    out.x = in.x0;
    out.r_nrm.clear();
    out.r_nrm.resize(in.max_it + 1);
    out.converged = false;

    // Compute the initial value for Ax
    vector<double> Ax;
    matvec(in.A, in.x0, Ax);
    
    // Compute the initial residual
    vecsub(in.b, Ax, out.r);
    out.r_nrm[0] = norm(out.r);

    if (out.r_nrm[0] <= in.tol) {
        return;
    }

    // Normalize Tolerance
    double rel_tol = norm(in.b)*in.tol;

    vector<vector<double>> V(out.x.size());
    for (size_t i = 0; i < V.size(); i++) {
        V[i].resize(in.restart + 1);
    }

    vector<double> column(out.x.size());

    vector<vector<double>> H(in.restart + 1);
    for (size_t i = 0; i < H.size(); i++) {
        H[i].resize(in.restart);
    }

    vector<double> cs(in.restart);
    vector<double> sn(in.restart);

    vector<double> e1(in.restart + 1);
    e1[0] = 1.0;

    out.iter = 0;
    bool notconv = true;
    while (notconv && out.iter < in.max_it) {
        // Arnoldi (Construct orthonormal basis)
        int64_t i = 0;

        for (size_t j = 0; j < out.r.size(); j++) {
            V[j][0] = out.r[j] / out.r_nrm[out.iter];
        }

        vector<double> s = e1;
        vecmult_emplace(out.r_nrm[out.iter], s);
        while (notconv && out.iter < in.max_it && i < (int64_t)in.restart) {
            out.iter++;
            vector<double> w(in.restart + 1);
            for (size_t j = 0; j < V.size(); j++) {
                column[j] = V[j][i];
            }
            matvec(in.A, column, w);
            for (int64_t k = 0; k <= i; k++) {
                for (size_t j = 0; j < w.size(); j++) {
                    H[k][i] += V[j][k] * w[j];
                }

                // w = w - H(k,i)*V(:,k);
                for (size_t j = 0; j < w.size(); j++) {
                    w[j] -= H[k][i] * V[j][k];
                }
            }

            H[i+1][i] = norm(w);
            if (real(H[i+1][i]) > 0) {
                for (size_t j = 0; j < w.size(); j++) {
                    V[j][i+1] = w[j] / H[i+1][i];
                }
                
                // Apply Givens rotation
                for (int64_t k = 0; k < i; k++) {
                    auto temp = cs[k]*H[k][i] + sn[k]*H[k + 1][i];
                    H[k + 1][i] = -sn[k]*H[k][i] + cs[k]*H[k + 1][i];
                    H[k][i] = temp;
                }
                // Form the i-th Givens rotation matrix
                rotmat(H[i][i], H[i+1][i], cs[i], sn[i]);

                // Approximate residual norm
                auto temp = cs[i]*s[i];
                s[i + 1] = -sn[i]*s[i];
                s[i] = temp;

                H[i][i] = cs[i]*H[i][i] + sn[i]*H[i + 1][i];
                H[i + 1][i] = 0.0;
                out.r_nrm[out.iter] = abs(s[i + 1]);

                if (out.r_nrm[out.iter] <= rel_tol) {
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
            vector<double> y;
            backsub(H, i+1, s, y);
            // x = x + V(:,1:i)*y;

            matvecadd_emplace(V, in.x0.size(), i+1, y, out.x);
        }
        else {
            vector<double> y;
            backsub(H, in.restart, s, y);

            matvecadd_emplace(V, in.x0.size(), in.restart, y, out.x);
            i--;
        }

        // Compute new Ax
        matvec(in.A, out.x, Ax);

        // Compute residual
        vecsub(in.b, Ax, out.r);

        // Compute norm of residual
        out.r_nrm[out.iter] = norm(out.r);

        if (out.r_nrm[out.iter] <= rel_tol) {
            notconv = false;
        }
        else {
            if (!notconv) {
                spdlog::warn("Encountered false convergence at iteration {}", out.iter);
            }
            notconv = true;
        }
    }

    // Did we converge?
    if (out.r_nrm[out.iter] <= rel_tol) {
        out.converged = true;
    }

    out.iter++;

    // Eleminate excess size in residual
    out.r_nrm.resize(out.iter);
}

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in &GMRES_In<complex<double>> input struct
 * @param out &GMRES_Out<complex<double>> output struct 
 */
void sgmres_new(GMRES_In<complex<double>> &in, GMRES_Out<complex<double>> &out) {
    assert(in.A.get_row_count() > 0 && in.A.get_col_count() > 0);
    assert(in.A.get_col_count() == in.x0.size());
    assert(in.A.get_row_count() == in.b.size());

    out.x = in.x0;
    out.r_nrm.clear();
    out.r_nrm.resize(in.max_it + 1);
    out.converged = false;

    // Compute the initial value for Ax
    vector<complex<double>> Ax;
    matvec(in.A, in.x0, Ax);
    
    // Compute the initial residual
    vecsub(in.b, Ax, out.r);
    out.r_nrm[0] = norm(out.r);

    if (out.r_nrm[0] <= in.tol) {
        return;
    }

    // Normalize Tolerance
    double rel_tol = norm(in.b)*in.tol;

    vector<vector<complex<double>>> V(in.restart + 1);
    vector<vector<complex<double>>> H(in.restart + 1);
    for (size_t i = 0; i < in.restart + 1; i++) {
        H[i].resize(in.restart);
    }

    vector<complex<double>> cs(in.restart);
    vector<complex<double>> sn(in.restart);

    vector<complex<double>> e1(in.restart + 1);
    e1[0] = 1.0;

    out.iter = 0;
    bool notconv = true;
    while (notconv && out.iter < in.max_it) {
        // Arnoldi (Construct orthonormal basis)
        int64_t i = 0;
        V[0] = out.r;
        vecdiv_emplace(out.r_nrm[out.iter], V[0]);

        vector<complex<double>> s = e1;
        vecmult_emplace(out.r_nrm[out.iter], s);
        while (notconv && out.iter < in.max_it && i < (int64_t)in.restart) {
            out.iter++;
            matvec(in.A, V[i], V[i + 1]);
            for (int64_t k = 0; k <= i; k++) {
                H[k][i] = dot(V[k], V[i + 1]);

                // w = w - H(k,i)*V(:,k);
                vecaddmult_emplace(V[i + 1], V[k], -H[k][i]);
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
                out.r_nrm[out.iter] = abs(s[i + 1]);

                if (out.r_nrm[out.iter] <= rel_tol) {
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
            matvecaddT_emplace(V, in.x0.size(), i+1, y, out.x);
        }
        else {
            vector<complex<double>> y;
            backsub(H, in.restart, s, y);
            matvecaddT_emplace(V, in.x0.size(), in.restart, y, out.x);
            i--;
        }

        // Compute new Ax
        matvec(in.A, out.x, Ax);

        // Compute residual
        vecsub(in.b, Ax, out.r);

        // Compute norm of residual
        out.r_nrm[out.iter] = norm(out.r);

        if (out.r_nrm[out.iter] <= rel_tol) {
            notconv = false;
        }
        else {
            if (!notconv) {
                spdlog::warn("Encountered false convergence at iteration {}", out.iter);
            }
            notconv = true;
        }
    }

    // Did we converge?
    if (out.r_nrm[out.iter] <= rel_tol) {
        out.converged = true;
    }

    out.iter++;

    // Eleminate excess size in residual
    out.r_nrm.resize(out.iter);
}

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in &GMRES_In<double> input struct
 * @param out &GMRES_Out<double> output struct 
 */
void sgmres_new(GMRES_In<double> &in, GMRES_Out<double> &out) {
    assert(in.A.get_row_count() > 0 && in.A.get_col_count() > 0);
    assert(in.A.get_col_count() == in.x0.size());
    assert(in.A.get_row_count() == in.b.size());

    out.x = in.x0;
    out.r_nrm.clear();
    out.r_nrm.resize(in.max_it + 1);
    out.converged = false;

    // Compute the initial value for Ax
    vector<double> Ax;
    matvec(in.A, in.x0, Ax);
    
    // Compute the initial residual
    vecsub(in.b, Ax, out.r);
    out.r_nrm[0] = norm(out.r);

    if (out.r_nrm[0] <= in.tol) {
        return;
    }

    // Normalize Tolerance
    double rel_tol = norm(in.b)*in.tol;

    vector<vector<double>> V(in.restart + 1);
    vector<vector<double>> H(in.restart + 1);
    for (size_t i = 0; i < in.restart + 1; i++) {
        H[i].resize(in.restart);
    }

    vector<double> cs(in.restart);
    vector<double> sn(in.restart);

    vector<double> e1(in.restart + 1);
    e1[0] = 1.0;

    out.iter = 0;
    bool notconv = true;
    while (notconv && out.iter < in.max_it) {
        // Arnoldi (Construct orthonormal basis)
        int64_t i = 0;
        V[0] = out.r;
        vecdiv_emplace(out.r_nrm[out.iter], V[0]);

        vector<double> s = e1;
        vecmult_emplace(out.r_nrm[out.iter], s);
        while (notconv && out.iter < in.max_it && i < (int64_t)in.restart) {
            out.iter++;
            matvec(in.A, V[i], V[i + 1]);
            for (int64_t k = 0; k <= i; k++) {
                H[k][i] = dot_simd(V[k], V[i + 1]);

                // w = w - H(k,i)*V(:,k);
                vecaddmult_emplace(V[i + 1], V[k], -H[k][i]);
            }

            H[i+1][i] = norm(V[i + 1]);
            if (real(H[i+1][i]) > 0) {
                vecdiv_emplace(H[i+1][i], V[i + 1]);
                
                // Apply Givens rotation
                for (int64_t k = 0; k < i; k++) {
                    auto temp = cs[k]*H[k][i] + sn[k]*H[k + 1][i];
                    H[k + 1][i] = -sn[k]*H[k][i] + cs[k]*H[k + 1][i];
                    H[k][i] = temp;
                }
                // Form the i-th Givens rotation matrix
                rotmat(H[i][i], H[i+1][i], cs[i], sn[i]);

                // Approximate residual norm
                auto temp = cs[i]*s[i];
                s[i + 1] = -sn[i]*s[i];
                s[i] = temp;

                H[i][i] = cs[i]*H[i][i] + sn[i]*H[i + 1][i];
                H[i + 1][i] = 0.0;
                out.r_nrm[out.iter] = abs(s[i + 1]);

                if (out.r_nrm[out.iter] <= rel_tol) {
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
            vector<double> y;
            backsub(H, i+1, s, y);
            // x = x + V(:,1:i)*y;
            matvecaddT_emplace(V, in.x0.size(), i+1, y, out.x);
        }
        else {
            vector<double> y;
            backsub(H, in.restart, s, y);
            matvecaddT_emplace(V, in.x0.size(), in.restart, y, out.x);
            i--;
        }

        // Compute new Ax
        matvec(in.A, out.x, Ax);

        // Compute residual
        vecsub(in.b, Ax, out.r);

        // Compute norm of residual
        out.r_nrm[out.iter] = norm(out.r);

        if (out.r_nrm[out.iter] <= rel_tol) {
            notconv = false;
        }
        else {
            if (!notconv) {
                spdlog::warn("Encountered false convergence at iteration {}", out.iter);
            }
            notconv = true;
        }
    }

    // Did we converge?
    if (out.r_nrm[out.iter] <= rel_tol) {
        out.converged = true;
    }

    out.iter++;

    // Eleminate excess size in residual
    out.r_nrm.resize(out.iter);
}