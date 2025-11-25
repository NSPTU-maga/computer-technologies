#include "fem.h"

//#define TEST1
#define TEST2
//#define TEST3

double func_rp(const point & p, const quadrilateral * r)
{
    (void)p;
    (void)r;
#if defined(TEST1)
    return p.x;
#endif
#if defined(TEST2)
    return -4.0 + p.x * p.x + p.y * p.y - p.x * p.y;//-2.0*exp(p.x+p.y);
#endif
#if defined(TEST3)
    return -0.02 * exp(0.1 * (p.x + p.y)) + exp(0.1 * (p.x + p.y));
#endif
}

size_t get_type_of_bound(size_t gmsh_num)
{
    (void)gmsh_num;
#if defined(TEST1)
    switch(gmsh_num)
    {
    case 600:
        return 2;
    case 601:
        return 1;
    case 602:
        return 2;
    case 603:
        return 1;
    case 604:
        return 2;
    case 605:
        return 1;
    case 606:
        return 2;
    case 607:
        return 1;
    }
    return 1;
#endif
#if defined(TEST2) || defined(TEST3)
    return 1;
#endif
}

double get_bound_value(const point & p, const edge & e)
{
    (void)p;
    (void)e;
#if defined(TEST1)
    switch(e.gmsh_phys_num)
    {
    case 600:
        return 0.0;
    case 601:
        return p.x;
    case 602:
        return 0.0;
    case 603:
        return p.x;
    case 604:
        return 0.0;
    case 605:
        return p.x;
    case 606:
    {
        double nx = - e.nodes[1]->y + e.nodes[0]->y;
        double ny = e.nodes[1]->x - e.nodes[0]->x;
        double norm_len = sqrt(nx * nx + ny * ny);
        double cosa = nx / norm_len;
        return cosa;
    }
    case 607:
        return p.x;
    }
    return p.x;
#endif
#if defined(TEST2)
    return p.x * p.x + p.y * p.y - p.x * p.y;//exp(p.x + p.y);
#endif
#if defined(TEST3)
    return exp(0.1 * (p.x + p.y));
#endif
}

#if !defined(HAVE_QT)
int main()
{
    FEM fem;
    fem.input();
    fem.make_portrait();
    fem.assembling_global();
    fem.applying_bounds();
    fem.slae.solve(1e-16);

    cout.precision(17);
    for(double i = 0.0; i <= 2.0; i += 0.25)
	{
        for(double j = 0.0; j <= 2.0; j += 0.25)
		{
            cout << i << " \t" << j << " \t" << fem.get_solution(point(i, j))
                 << " \t" << get_bound_value(point(i, j), edge()) << endl;
		}
	}

    ofstream result("result.txt");
    result.precision(17);
    for(double i = 0.0; i <= 12.0; i += 0.25)
    {
        for(double j = 0.0; j <= 4.0; j += 0.25)
        {
            result << i << " \t" << j << " \t" << fem.get_solution(point(i, j))
                 << " \t" << get_bound_value(point(i, j), edge()) << endl;
        }
    }
    result.flush();
    result.close();

#if defined(_WIN32)
	system("pause");
#endif

    return 0;
}
#endif
