#ifndef FEM_H
#define FEM_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "solvers/CGM_LLT.h"

using namespace std;

// Класс точка
class point
{
public:
    double x;
    double y;
    size_t num;
    point();
    point(double x, double y);
    point(double x, double y, size_t num);
    point(const point & p);
    double & operator [] (size_t i);
    double operator [] (size_t i) const;
    bool operator < (const point & t) const;
    bool operator == (const point & t) const;
    point & operator = (const point & other);
    friend ostream & operator << (ostream & os, point a);
};

// Класс физическая область
class phys_area
{
public:
    double lambda;
    double gamma;
    size_t num;
    size_t gmsh_phys_num;
};

// Класс четырехугольник
class quadrilateral
{
public:
    friend class glwidget;

    point * nodes[4];
    phys_area * ph;
    void init();
    double get_local_G(size_t i, size_t j) const;
    double get_local_M(size_t i, size_t j) const;
    double get_local_rp(size_t i) const;
    bool inside(const point & p) const;
    double bfunc_2d(size_t num, const point & p) const;
    point grad_bfunc_2d(size_t num, const point & p) const;

    size_t degree_of_freedom[9];
private:
    point to_local(const point & p) const;
    point to_global(const point & p) const;
    point gauss_points[9];
    double gauss_weights[9];
    double det_J_local(const point & p) const;

    double alpha0, alpha1, alpha2;
    double beta1, beta2, beta3, beta4, beta5, beta6;
    double S;

    void two2one(size_t two, size_t & one1, size_t & one2) const;
    double bfunc_1d(size_t func_n, double ksi) const;
    double dbfunc_1d(size_t func_n, double ksi) const;
    double bfunc_2d_local(size_t num, const point & p) const;
    point grad_bfunc_2d_local(size_t num, const point & p) const;
};

// Класс СЛАУ
class SLAE
{
public:
    SLAE();
    ~SLAE();
    void solve(double eps);
    void alloc_all(size_t n_size, size_t gg_size);
    void add(size_t i, size_t j, double elem);
    double * gg, * di, * rp, * x;
    size_t * ig, * jg;
    size_t n;
private:
    CGM_LLT solver;
};

// Класс ребро
class edge
{
public:
    const point * nodes[2];
    size_t gmsh_phys_num;
    size_t degree_of_freedom[3];
    bool operator < (const edge & t) const;
    bool operator == (const edge & t) const;
    edge();
    edge(const point * begin, const point * end);
    point get_freedom_position(size_t num) const;
    double get_matrix_M(size_t i, size_t j);
    double bfunc_1d(size_t func_n, double x) const;

    size_t num;
    const quadrilateral * fes[2];
};

// Класс МКЭ
class FEM
{
public:
    friend class glwidget;

    FEM();
    ~FEM();
    SLAE slae;

    void input();
    void make_portrait();
    void assembling_global();
    void applying_bounds();

    double get_solution(const point & p) const;
    double get_solution(const point & p, const quadrilateral * fe) const;
    point get_grad(const point & p) const;
    point get_grad(const point & p, const quadrilateral * fe) const;

    const quadrilateral * get_fe(const point & p) const;

    // lab02
    vector<edge>::const_iterator get_edge_by_nodes(const point & begin, const point & end) const;
    vector<edge>::const_iterator get_edge_by_fe(const quadrilateral & fe, size_t edge_num) const;
    const point * get_nodes_by_edge(const edge & e, size_t node_num) const;
    const quadrilateral * get_fe_by_edge(const edge & e, size_t fe_num) const;

private:
    point * nodes;
    size_t nodes_num;
    quadrilateral * qls;
    size_t qls_num;
    phys_area * phs;
    size_t phs_num;
    edge * bounds;
    size_t bounds_num;

    size_t degree_of_freedom_num;
    vector<edge> edges_freedom;
    void add_edge_freedom(const point * begin, const point * end, size_t ql_num);
};

double func_rp(const point & p, const quadrilateral * r);
size_t get_type_of_bound(size_t gmsh_num);
double get_bound_value(const point & p, const edge & e);

#endif // FEM_H
