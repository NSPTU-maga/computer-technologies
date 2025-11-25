#include "CGM_LLT.h"
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <cstdio>
#include <cmath>

void CGM_LLT::init(size_t * gi_s, size_t * gj_s, double * di_s,
                          double * gg_s, size_t n_s)
{
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gg = gg_s;
    n = n_s;

    size_t m = gi[n];
    if(!r) r = new double [n];
    if(!z) z = new double [n];
    if(!p) p = new double [n];
    if(!s) s = new double [n];

    if(!L_di) L_di = new double [n];
    if(!L_gg) L_gg = new double [m];

    for(size_t i = 0; i < n; i++)
    {
        L_di[i] = di[i];
        //x0[i] = 0.0;
    }

    for(size_t i = 0 ; i < m ; i++)
        L_gg[i] = gg[i];
}

void CGM_LLT::make_LLT_decomposition()
{
    double sum_d, sum_l;

    for(size_t k = 0; k < n ; k++)
    {
        sum_d = 0;
        size_t i_s = gi[k], i_e = gi[k+1];

        for(size_t i = i_s; i < i_e ; i++)
        {
            sum_l = 0;
            size_t j_s = gi[gj[i]], j_e = gi[gj[i]+1];

            for(size_t m = i_s; m < i; m++)
            {
                for(size_t j = j_s; j < j_e; j++)
                {
                    if(gj[m] == gj[j])
                    {
                        sum_l += L_gg[m] * L_gg[j];
                        j_s++;
                    }
                }
            }
            L_gg[i] = (L_gg[i] -  sum_l) / L_di[gj[i]];

            sum_d += L_gg[i] * L_gg[i];
        }
        L_di[k] = sqrt(L_di[k] - sum_d);
    }
}

double CGM_LLT::dot_prod(const double * a, const double * b) const
{
    double d_p = 0.0;
    for(size_t i = 0; i < n; i++)
        d_p += a[i] * b[i];
    return d_p;
}

void CGM_LLT::mul_matrix(const double * f, double * x) const
{
    for(size_t i = 0; i < n; i++)
    {
        double v_el = f[i];
        x[i] = di[i] * v_el;
        for(size_t k = gi[i], k1 = gi[i+1]; k < k1; k++)
        {
            size_t j = gj[k];
            x[i] += gg[k] * f[j];
            x[j] += gg[k] * v_el;
        }
    }
}

void CGM_LLT::solve_L(const double * f, double * x) const
{
    for(size_t k = 1, k1 = 0; k <= n; k++, k1++)
    {
        double sum = 0.0;

        for(size_t i = gi[k1]; i < gi[k]; i++)
            sum += L_gg[i] * x[gj[i]];

        x[k1] = (f[k1] - sum) / L_di[k1];
    }
}

void CGM_LLT::solve_LT(const double * f, double * x) const
{
    for(size_t k = n, k1 = n-1; k > 0; k--, k1--)
    {
        x[k1] = f[k1] / L_di[k1];
        double v_el = x[k1];

        for(size_t i = gi[k1]; i < gi[k]; i++)
            x[gj[i]] -= L_gg[i] * v_el;
    }
}

void CGM_LLT::solve_LLT(const double * f, double * x) const
{
    solve_L(f, x);
    solve_LT(x, x);
}

bool CGM_LLT::is_fpu_error(double x) const
{
    double y = x - x;
    return x != x || y != y;
}

void CGM_LLT::solve(double * solution, double * rp_s, double eps)
{
    // Параметры решателя
    size_t max_iter = 15000;

    rp = rp_s;

    x0 = new double [n];
    for(size_t i = 0; i < n; i++)
        x0[i] = solution[i];

    mul_matrix(x0, r);
    make_LLT_decomposition();

    for(size_t i = 0; i < n ; i++)
        r[i] = rp[i] - r[i];

    solve_LLT(r, z);
    for(size_t i = 0; i < n; i++)
        p[i] = z[i];

    double alpha, beta, prod_1, prod_2;
    double discr, rp_norm;

    rp_norm = sqrt(dot_prod(rp, rp));
    if(is_fpu_error(rp_norm))
    {
        fprintf(stderr, "Error: FPU error detected in right part!\n");
        for(size_t i = 0; i < n; i++)
            solution[i] = x0[i];
        delete [] x0;
        return;
    }

    prod_1 = dot_prod(p, r);

    bool finished = false;

    size_t iter;
    for(iter = 1; iter <= max_iter && !finished; iter++)
    {
        discr = sqrt(dot_prod(r, r));
        if(is_fpu_error(discr))
        {
            fprintf(stderr, "Error: FPU error detected in (r, r)!\n");
            for(size_t i = 0; i < n; i++)
                solution[i] = x0[i];
            delete [] x0;
            return;
        }

        if(iter%10 == 0)
        {
            printf("CGM_LLT Residual:\t%5lu\t%.3e\r", (unsigned long)iter, discr / rp_norm);
            fflush(stdout);
        }

        if(eps < discr / rp_norm)
        {
            mul_matrix(z, s);

            alpha = prod_1 / dot_prod(s, z);

            for(size_t i = 0; i < n ; i++)
            {
                x0[i] += alpha * z[i];
                r[i] -= alpha * s[i];
            }

            solve_LLT(r, p);
            prod_2 = dot_prod(p, r);

            beta = prod_2 / prod_1;

            prod_1 = prod_2;

            for(size_t i = 0; i < n; i++)
                z[i] = p[i] + beta * z[i];
        }
        else
            finished = true;
    }

    mul_matrix(x0, r);
    for(size_t i = 0; i < n; i++)
        r[i] = rp[i] - r[i];
    discr = sqrt(dot_prod(r, r));
    printf("CGM_LLT Residual:\t%5lu\t%.3e\n", (unsigned long)iter, discr / rp_norm);

    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");

    for(size_t i = 0; i < n; i++)
        solution[i] = x0[i];
    delete [] x0;
}

CGM_LLT::CGM_LLT()
{
    r = x0 = z = p = s = L_di = L_gg = NULL;
}

CGM_LLT::~CGM_LLT()
{
    if(r) delete [] r;
    if(z) delete [] z;
    if(p) delete [] p;
    if(s) delete [] s;
    if(L_di) delete [] L_di;
    if(L_gg) delete [] L_gg;
}
