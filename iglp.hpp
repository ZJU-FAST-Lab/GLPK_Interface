#ifndef IGLP_HPP
#define IGLP_HPP

#include <glpk.h>
#include <Eigen/Eigen>

#include <cstdio>
#include <cstdlib>
#include <string>

namespace iglp
{

inline double linprog(const Eigen::VectorXd &c,
                      const Eigen::MatrixXd &A,
                      const Eigen::VectorXd &b,
                      Eigen::VectorXd &x,
                      bool ipm = false,
                      bool verbose = false)
// linprog:
//         min cTx s.t. Ax<=b
// input:
//         c: d*1 objective coeffs
//         A: m*d constraint matrix
//         b: m*1 constraint bound
//         ipm: use interior point method
//              or simplex method
//         verbose: show details
// output:
//         x: d*1 decision variables
// return:
//         inf: No feasible solution or fail
//        -inf: Unbounded problem
//         real: minimum objective function
{
    int d = c.size();
    int m = b.size();
    int dm = d * m;
    x = Eigen::VectorXd::Zero(d);

    glp_prob *lp;
    int *ia = new int[dm + 1];
    int *ja = new int[dm + 1];
    double *ar = new double[dm + 1];
    int s;
    double z;

    lp = glp_create_prob();
    glp_set_prob_name(lp, "lp");
    glp_set_obj_dir(lp, GLP_MIN);

    glp_add_rows(lp, m);
    for (int i = 1; i <= m; i++)
    {
        glp_set_row_name(lp, i, (std::to_string(i) + "y").c_str());
        glp_set_row_bnds(lp, i, GLP_UP, 0.0, b(i - 1));
    }

    glp_add_cols(lp, d);
    for (int i = 1; i <= d; i++)
    {
        glp_set_col_name(lp, i, (std::to_string(i) + "x").c_str());
        glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
        glp_set_obj_coef(lp, i, c(i - 1));
    }

    int k = 1;
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= d; j++)
        {
            ia[k] = i;
            ja[k] = j;
            ar[k] = A(i - 1, j - 1);
            k++;
        }
    }
    glp_load_matrix(lp, dm, ia, ja, ar);

    if (!ipm)
    {
        glp_smcp param;
        glp_init_smcp(&param);
        param.msg_lev = verbose ? GLP_MSG_ALL : GLP_MSG_OFF;
        glp_simplex(lp, &param);
        s = glp_get_status(lp);
        z = INFINITY;
        if (s == GLP_OPT || s == GLP_UNBND)
        {
            z = (s == GLP_UNBND) ? -INFINITY : glp_get_obj_val(lp);
            for (int i = 1; i <= d; i++)
            {
                x(i - 1) = glp_get_col_prim(lp, i);
            }
        }
    }
    else
    {
        glp_iptcp param;
        glp_init_iptcp(&param);
        param.msg_lev = verbose ? GLP_MSG_ALL : GLP_MSG_OFF;
        glp_interior(lp, &param);
        s = glp_ipt_status(lp);
        z = INFINITY;
        if (s == GLP_OPT)
        {
            z = glp_ipt_obj_val(lp);
            for (int i = 1; i <= d; i++)
            {
                x(i - 1) = glp_ipt_col_prim(lp, i);
            }
        }
    }

    glp_delete_prob(lp);
    delete[] ia;
    delete[] ja;
    delete[] ar;

    return z;
}

} // namespace iglp

#endif