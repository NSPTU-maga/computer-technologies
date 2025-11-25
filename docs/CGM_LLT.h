#ifndef CGM_LLT_H_INCLUDED
#define CGM_LLT_H_INCLUDED

#include <cstdlib>

class CGM_LLT
{
public:
    void init(size_t * gi_s, size_t * gj_s, double * di_s,
              double * gg_s, size_t n_s);
    void solve(double * solution, double * rp_s, double eps);

    CGM_LLT();
    ~CGM_LLT();
private:
    void make_LLT_decomposition();
    void mul_matrix(const double * f, double * x) const;
    void solve_L(const double * f, double * x) const;
    void solve_LT(const double * f, double * x) const;
    void solve_LLT(const double * f, double * x) const;
    double dot_prod(const double * a, const double * b) const;
    bool is_fpu_error(double x) const;

    size_t n;
    size_t * gi, * gj;
    double * di, * gg, * rp, * r, * x0, * z, * p, * s;
    double * L_di, * L_gg;
};

#endif // CGM_LLT_H_INCLUDED
