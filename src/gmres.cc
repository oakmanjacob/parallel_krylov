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
 * 
 * @param in &GMRES_In<double> input struct
 * @param out &GMRES_Out<double> output struct 
 */
void gmres(const GMRES_In<double> &in, GMRES_Out<double> &out) {
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

void ogmres_helper(const GMRES_In<double> &in, GMRES_Out<double> &out,
    const function<double(const vector<double>&,const vector<double>&)> dotf) {
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
                H[k][i] = dotf(V[k], V[i + 1]);

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

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in &GMRES_In<double> input struct
 * @param out &GMRES_Out<double> output struct 
 */
void ogmres(const GMRES_In<double> &in, GMRES_Out<double> &out) {
    ogmres_helper(in, out, dot);
}

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in &GMRES_In<double> input struct
 * @param out &GMRES_Out<double> output struct 
 */
void ogmres_simd(const GMRES_In<double> &in, GMRES_Out<double> &out) {
    ogmres_helper(in, out, dot_simd);
}

struct PGMRES_TConfig {
    vector<vector<double>> H;
    vector<vector<double>> V;
    vector<double> s;
    vector<double> cs;
    vector<double> sn;
    vector<double> Ax;
    bool notconv;
    double rel_tol;
    int64_t i = 0;
};

void pgmres_sync_helper(PGMRES_TConfig &tconf, const GMRES_In<double> &in, GMRES_Out<double> &out) {
    tconf.V[0] = out.r;

    vecdiv_emplace(out.r_nrm[out.iter], tconf.V[0]);
    vecmult_emplace(out.r_nrm[out.iter], tconf.s);

    while (tconf.notconv && out.iter < in.max_it && tconf.i < (int64_t)in.restart) {
        out.iter++;
        matvec(in.A, tconf.V[tconf.i], tconf.V[tconf.i + 1]);
        for (int64_t k = 0; k <= tconf.i; k++) {
            tconf.H[k][tconf.i] = dot_simd(tconf.V[k], tconf.V[tconf.i + 1]);

            // w = w - H(k,i)*V(:,k);
            vecaddmult_emplace(tconf.V[tconf.i + 1], tconf.V[k], -tconf.H[k][tconf.i]);
        }

        tconf.H[tconf.i+1][tconf.i] = norm(tconf.V[tconf.i + 1]);
        if (real(tconf.H[tconf.i+1][tconf.i]) > 0) {
            vecdiv_emplace(tconf.H[tconf.i+1][tconf.i], tconf.V[tconf.i + 1]);
            
            // Apply Givens rotation
            for (int64_t k = 0; k < tconf.i; k++) {
                auto temp = tconf.cs[k]*tconf.H[k][tconf.i] + tconf.sn[k]*tconf.H[k + 1][tconf.i];
                tconf.H[k + 1][tconf.i] = -tconf.sn[k]*tconf.H[k][tconf.i] + tconf.cs[k]*tconf.H[k + 1][tconf.i];
                tconf.H[k][tconf.i] = temp;
            }
            // Form the i-th Givens rotation matrix
            rotmat(tconf.H[tconf.i][tconf.i], tconf.H[tconf.i+1][tconf.i], tconf.cs[tconf.i], tconf.sn[tconf.i]);

            // Approximate residual norm
            auto temp = tconf.cs[tconf.i]*tconf.s[tconf.i];
            tconf.s[tconf.i + 1] = -tconf.sn[tconf.i]*tconf.s[tconf.i];
            tconf.s[tconf.i] = temp;

            tconf.H[tconf.i][tconf.i] = tconf.cs[tconf.i]*tconf.H[tconf.i][tconf.i] + tconf.sn[tconf.i]*tconf.H[tconf.i + 1][tconf.i];
            tconf.H[tconf.i + 1][tconf.i] = 0.0;
            out.r_nrm[out.iter] = abs(tconf.s[tconf.i + 1]);

            if (out.r_nrm[out.iter] <= tconf.rel_tol) {
                tconf.notconv = false;
            }
            else {
                tconf.i++;
            }
        }
        else {
            tconf.notconv = false;
        }
    }

    // Form the (approximate) solution
    if (!tconf.notconv) {
        vector<double> y;
        backsub(tconf.H, tconf.i+1, tconf.s, y);
        // x = x + V(:,1:i)*y;
        matvecaddT_emplace(tconf.V, in.x0.size(), tconf.i+1, y, out.x);
    }
    else {
        vector<double> y;
        backsub(tconf.H, in.restart, tconf.s, y);
        matvecaddT_emplace(tconf.V, in.x0.size(), in.restart, y, out.x);
        tconf.i--;
    }

    // Compute new Ax
    matvec(in.A, out.x, tconf.Ax);

    // Compute residual
    vecsub(in.b, tconf.Ax, out.r);

    // Compute norm of residual
    out.r_nrm[out.iter] = norm(out.r);

    if (out.r_nrm[out.iter] <= tconf.rel_tol) {
        tconf.notconv = false;
    }
    else {
        if (!tconf.notconv) {
            spdlog::warn("Encountered false convergence at iteration {}", out.iter);
        }
        tconf.notconv = true;
    }
}

/**
 * @brief Perform GMRES to solve the equation Ax = b without preconditioning
 * 
 * @param in &GMRES_In<double> input struct
 * @param out &GMRES_Out<double> output struct 
 */
void pgmres_sync(const GMRES_In<double> &in, GMRES_Out<double> &out, size_t thread_count, double priority) {
    assert(in.A.get_row_count() > 0 && in.A.get_col_count() > 0);
    assert(in.A.get_col_count() == in.x0.size());
    assert(in.A.get_row_count() == in.b.size());

    
    out.r_nrm.clear();
    out.r_nrm.resize(in.max_it + 1);
    out.converged = false;

    vector<thread> threads;

    bool notconv = true;
    vector<PGMRES_TConfig> tconfigs(thread_count);

    vector<GMRES_Out<double>> outs(thread_count);

    vector<double> e1(in.restart + 1);
    e1[0] = 1.0;

    for (size_t i = 0; i < thread_count; i++) {
        outs[i].x.clear();
        outs[i].x.resize(in.x0.size());
        outs[i].x[i] = 1;
        outs[i].r_nrm.resize(in.max_it + 1);

        // Compute the initial value for Ax
        matvec(in.A, in.x0, tconfigs[i].Ax);
        
        // Compute the initial residual
        vecsub(in.b, tconfigs[i].Ax, outs[i].r);
        outs[i].r_nrm[0] = norm(outs[i].r);

        if (outs[i].r_nrm[0] <= in.tol) {
            out = outs[i];
            return;
        }

        // Normalize Tolerance
        tconfigs[i].rel_tol = norm(in.b)*in.tol;

        tconfigs[i].V.resize(in.restart + 1);
        tconfigs[i].H.resize(in.restart + 1);
        for (size_t j = 0; j < in.restart + 1; j++) {
            tconfigs[i].H[j].resize(in.restart);
        }

        tconfigs[i].cs.resize(in.restart);
        tconfigs[i].sn.resize(in.restart);

        tconfigs[i].notconv = true;
    }

    while (notconv && out.iter < in.max_it) {
        for (size_t i = 0; i < thread_count; i++) {
            tconfigs[i].s = e1;
            threads.emplace_back(pgmres_sync_helper, ref(tconfigs[i]), ref(in), ref(outs[i]));
        }

        for (size_t i = 0; i < thread_count; i++) {
            threads.at(i).join();
        }
        threads.clear();

        for (size_t i = 0; i < thread_count; i++) {
            if (tconfigs[i].notconv && outs[i].iter < in.max_it) {
                tconfigs[i].s = e1;

                // Do some weird ass shit
                for (size_t j = 0; j < thread_count; j++) {
                    if (i != j) {
                        vecaddmult_emplace(outs[i].x, outs[j].x, priority);
                    }
                }
            }
            else {
                outs[i].iter++;
                out = outs[i];

                // Eleminate excess size in residual
                out.r_nrm.resize(out.iter);

                // Did we converge?
                if (out.r_nrm[out.iter] <= tconfigs[i].rel_tol) {
                    out.converged = true;
                }

                return;

            }
        }
    }
}