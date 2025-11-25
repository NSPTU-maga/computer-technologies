#include "fem.h"

// ============================================================================

SLAE::SLAE()
{
    gg = di = rp = x = NULL;
    ig = jg = NULL;
    n = 0;
}

SLAE::~SLAE()
{
    if(gg) delete [] gg;
    if(di) delete [] di;
    if(rp) delete [] rp;
    if(x)  delete [] x;
    if(ig) delete [] ig;
    if(jg) delete [] jg;
}

void SLAE::solve(double eps)
{
    cout << "Solving SLAE ..." << endl;
    solver.init(ig, jg, di, gg, n);
    solver.solve(x, rp, eps);
}

void SLAE::alloc_all(size_t n_size, size_t gg_size)
{
    n = n_size;
    ig = new size_t [n + 1];
    jg = new size_t [gg_size];
    di = new double [n];
    gg = new double [gg_size];
    rp = new double [n];
    x  = new double [n];

    memset(di, 0, sizeof(double) * n);
    memset(rp, 0, sizeof(double) * n);
    memset(x,  0, sizeof(double) * n);
    memset(gg, 0, sizeof(double) * gg_size);
}

void SLAE::add(size_t i, size_t j, double elem)
{
    if(i == j)
    {
        di[j] += elem;
        return;
    }

    if(j > i)
        swap(i, j);
    size_t ind = 0;
    bool flag = false;
    for(size_t k = ig[i]; k < ig[i + 1] && !flag; k++)
    {
        if(jg[k] == j)
        {
            ind = k;
            flag = true;
        }
    }
    gg[ind] += elem;
}

// ============================================================================

FEM::FEM()
{
    nodes = NULL;
    qls = NULL;
    phs   = NULL;
    bounds = NULL;

    nodes_num = 0;
    qls_num = 0;
    phs_num   = 0;
    bounds_num = 0;
}

FEM::~FEM()
{
    if(nodes) delete [] nodes;
    if(qls) delete [] qls;
    if(phs)   delete [] phs;
    if(bounds) delete [] bounds;
}

void FEM::add_edge_freedom(const point * begin, const point * end, size_t ql_num)
{
    edge ed(begin, end);
    edge ei(end, begin);
    vector<edge>::iterator id = find(edges_freedom.begin(), edges_freedom.end(), ed);
    vector<edge>::iterator ii = find(edges_freedom.begin(), edges_freedom.end(), ei);
    if(id == edges_freedom.end() && ii == edges_freedom.end())
    {
        edges_freedom.push_back(ed);
        edges_freedom[edges_freedom.size()-1].fes[0] = qls + ql_num;
    }
    else
    {
        vector<edge>::iterator it;
        edge et;
        if(ii == edges_freedom.end())
        {
            it = id;
            et = ed;
        }
        else
        {
            it = ii;
            et = ei;
        }
        if(!it->fes[0])
            it->fes[0] = qls + ql_num;
        else if(!it->fes[1])
            it->fes[1] = qls + ql_num;
        else
            assert(!it->fes[0] || !it->fes[1]);
    }
}

void FEM::input()
{
    ifstream ifs;
    cout << "Reading data..." << endl;

    ifs.open("data/phys.txt", ios::in);
    assert(ifs.good());
    ifs >> phs_num;
    cout << " > detected " << phs_num << " physical areas" << endl;
    phs = new phys_area [phs_num];
    for(size_t i = 0; i < phs_num; i++)
    {
        ifs >> phs[i].gmsh_phys_num;
        ifs >> phs[i].lambda;
        ifs >> phs[i].gamma;
        phs[i].num = i;
    }
    ifs.close();

    string line;
    double garbage;
    ifstream gmsh_file;
    gmsh_file.open("data/mesh.msh", ios::in);
    assert(gmsh_file.good());
    do
        getline(gmsh_file, line);
    while(line.find("$Nodes") == string::npos && gmsh_file.good());

    gmsh_file >> nodes_num;
    cout << " > detected " << nodes_num << " nodes" << endl;
    nodes = new point[nodes_num];
    for(size_t i = 0; i < nodes_num; i++)
    {
        gmsh_file >> garbage >> nodes[i].x >> nodes[i].y >> garbage;
        nodes[i].num = i;
    }

    do
        getline(gmsh_file, line);
    while(line.find("$Elements") == string::npos && gmsh_file.good());

    size_t num_elem;
    size_t type_elem;
    vector<edge> tmp_bounds;
    vector<quadrilateral> tmp_rects;
    gmsh_file >> num_elem;
    for(size_t i = 0; i < num_elem; i++)
    {
        gmsh_file >> garbage >> type_elem;

        if(type_elem == 1)
        {
            edge tmp_edge;
            gmsh_file >> garbage >> tmp_edge.gmsh_phys_num >> garbage;
            for(int j = 0; j < 2; j++)
            {
                size_t tmp_node;
                gmsh_file >> tmp_node;
                tmp_edge.nodes[j] = nodes + tmp_node - 1;
            }
            if(tmp_edge.nodes[0]->num > tmp_edge.nodes[1]->num)
                swap(tmp_edge.nodes[0], tmp_edge.nodes[1]);
            tmp_bounds.push_back(tmp_edge);
        }
        else if(type_elem == 3)
        {
            quadrilateral tmp_rect;

            size_t gmsh_phys_num;
            gmsh_file >> garbage >> gmsh_phys_num >> garbage;
            tmp_rect.ph = NULL;
            for(size_t j = 0; j < phs_num; j++)
                if(phs[j].gmsh_phys_num == gmsh_phys_num)
                    tmp_rect.ph = phs + j;

            for(int j = 0; j < 4; j++)
            {
                size_t tmp_node;
                gmsh_file >> tmp_node;
                tmp_rect.nodes[j] = nodes + tmp_node - 1;
            }
            swap(tmp_rect.nodes[2], tmp_rect.nodes[3]);

            tmp_rects.push_back(tmp_rect);
        }
    }
    gmsh_file.close();

    qls_num = tmp_rects.size();
    cout << " > detected " << qls_num << " rectangles" << endl;
    qls = new quadrilateral [qls_num];
    for(size_t i = 0; i < qls_num; i++)
        qls[i] = tmp_rects[i];
    tmp_rects.clear();

    bounds_num = tmp_bounds.size();
    cout << " > detected " << bounds_num << " bounds" << endl;
    bounds = new edge [bounds_num];
    for(size_t i = 0; i < bounds_num; i++)
        bounds[i] = tmp_bounds[i];
    tmp_bounds.clear();

    // Вот тут будем делать степени свободы
    for(size_t i = 0; i < qls_num; i++)
    {
        add_edge_freedom(qls[i].nodes[0], qls[i].nodes[1], i);
        add_edge_freedom(qls[i].nodes[0], qls[i].nodes[2], i);
        add_edge_freedom(qls[i].nodes[1], qls[i].nodes[3], i);
        add_edge_freedom(qls[i].nodes[2], qls[i].nodes[3], i);
    }
    sort(edges_freedom.begin(), edges_freedom.end());
    vector<edge>::iterator it;
    it = unique(edges_freedom.begin(), edges_freedom.end());
    edges_freedom.resize(distance(edges_freedom.begin(), it));

    // Занумеруем ребра
    for(size_t i = 0; i < edges_freedom.size(); i++)
        edges_freedom[i].num = i;

    // У узлов пусть будут номера узлов
    for(size_t i = 0; i < edges_freedom.size(); i++)
    {
        edges_freedom[i].degree_of_freedom[0] = edges_freedom[i].nodes[0]->num;
        edges_freedom[i].degree_of_freedom[2] = edges_freedom[i].nodes[1]->num;
    }
    // Переменная, содержащая свободную степень свободы
    size_t curr_freedom = nodes_num;
    // Теперь заполним степени свободы на ребрах
    for(size_t i = 0; i < edges_freedom.size(); i++)
    {
        edges_freedom[i].degree_of_freedom[1] = curr_freedom;
        curr_freedom++;
    }
    // Дальше заполним центры и узлы, их сразу
    for(size_t i = 0; i < qls_num; i++)
    {
        qls[i].degree_of_freedom[4] = curr_freedom;
        curr_freedom++;
        qls[i].degree_of_freedom[0] = qls[i].nodes[0]->num;
        qls[i].degree_of_freedom[2] = qls[i].nodes[1]->num;
        qls[i].degree_of_freedom[6] = qls[i].nodes[2]->num;
        qls[i].degree_of_freedom[8] = qls[i].nodes[3]->num;
    }
    // Теперь и то, что на ребрах
    for (size_t i = 0; i < qls_num; i++)
    {
        vector<edge>::iterator it;
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[0], qls[i].nodes[1]));
        assert(it != edges_freedom.end());
        qls[i].degree_of_freedom[1] = it->degree_of_freedom[1];
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[0], qls[i].nodes[2]));
        assert(it != edges_freedom.end());
        qls[i].degree_of_freedom[3] = it->degree_of_freedom[1];
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[1], qls[i].nodes[3]));
        assert(it != edges_freedom.end());
        qls[i].degree_of_freedom[5] = it->degree_of_freedom[1];
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[2], qls[i].nodes[3]));
        assert(it != edges_freedom.end());
        qls[i].degree_of_freedom[7] = it->degree_of_freedom[1];
    }

    // Заполним также инфу о степенях свободы в ребрах для краевых
    for(size_t i = 0; i < bounds_num; i++)
    {
        vector<edge>::iterator it;
        it = find(edges_freedom.begin(), edges_freedom.end(), bounds[i]);
        if(it == edges_freedom.end()) // Здравствуй жопа новый год
        {
            swap(bounds[i].nodes[0], bounds[i].nodes[1]);
            it = find(edges_freedom.begin(), edges_freedom.end(), bounds[i]);
        }
        assert(it != edges_freedom.end());
        for(size_t j = 0; j < 3; j++)
            bounds[i].degree_of_freedom[j] = it->degree_of_freedom[j];
        // За одно и инфу о КЭ
        for(size_t j = 0; j < 2; j++)
            bounds[i].fes[j] = it->fes[j];
    }

    degree_of_freedom_num = curr_freedom;
    cout << " > detected " << degree_of_freedom_num << " degree of freedom" << endl;

    // Перенумерация степеней свободы
    map<size_t, size_t> convert;
    size_t curr_deg = 0;
    for(size_t i = 0; i < qls_num; i++)
    {
        for(size_t j = 0; j < 9; j++)
        {
            if(convert.find(qls[i].degree_of_freedom[j]) == convert.end())
            {
                convert[qls[i].degree_of_freedom[j]] = curr_deg;
                curr_deg++;
            }
        }
    }

    // Применение перенумерованных значений
    // В четырехугольниках
    for(size_t i = 0; i < qls_num; i++)
        for(size_t j = 0; j < 9; j++)
            qls[i].degree_of_freedom[j] = convert[qls[i].degree_of_freedom[j]];
    // В ребрах
    for(size_t i = 0; i < edges_freedom.size(); i++)
        for(size_t j = 0; j < 3; j++)
            edges_freedom[i].degree_of_freedom[j] = convert[edges_freedom[i].degree_of_freedom[j]];
    // В краевых
    for(size_t i = 0; i < bounds_num; i++)
        for(size_t j = 0; j < 3; j++)
            bounds[i].degree_of_freedom[j] = convert[bounds[i].degree_of_freedom[j]];

    for(size_t i = 0; i < qls_num; i++)
        qls[i].init();
}

void FEM::make_portrait()
{
    // Формирование профиля (портрета)
    cout << "Generating profile..." << endl;
    size_t gg_size = 0;
    // Создаем массив списков для хранения связей
    set<size_t> * profile = new set<size_t> [degree_of_freedom_num];

    // Связь есть, если узлы принадлежат одному КЭ
    // Поэтому обходим конечные элементы и добавляем в список общие вершины
    for(size_t i = 0; i < qls_num; i++)
        for(size_t j = 0; j < 9; j++)
            for(size_t k = 0; k < j; k++)
            {
                size_t ik = qls[i].degree_of_freedom[k];
                size_t ij = qls[i].degree_of_freedom[j];
                if(ik > ij) swap(ik, ij);
                profile[ik].insert(ij);
            }

    // Удаляем повторяющиеся записи в спиках и сортируем их, попутно считая размер матрицы
    for(size_t i = 0; i < degree_of_freedom_num; i++)
        gg_size += profile[i].size();
    slae.alloc_all(degree_of_freedom_num, gg_size);

    cout << " > slae.n_size = " << degree_of_freedom_num << endl;
    cout << " > slae.gg_size = " << gg_size << endl;

    // Заполнение профиля (портрета)
    slae.ig[0] = 0;
    slae.ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < slae.n; i++)
    {
        size_t k = 0;
        for(size_t j = 0; j <= i; j++)
        {
            // Если есть связь между i и j, значит в этом месте матрицы будет ненулевой элемент
            // занесем информацию об этом в jg
            if(profile[j].find(i) != profile[j].end())
            {
                slae.jg[tmp] = j;
                tmp++;
                k++;
            }
        }
        // а в ig занесем информацию о количестве ненулевых элементов в строке
        slae.ig[i + 1] = slae.ig[i] + k;
    }

    // Очистка списков
    for(size_t i = 0; i < degree_of_freedom_num; i++)
        profile[i].clear();
    delete [] profile;
}

void FEM::assembling_global()
{
    cout << "Assembling global matrix..." << endl;
    for(size_t k = 0; k < qls_num; k++)
    {
        for(size_t i = 0; i < 9; i++)
        {
            for(size_t j = 0; j <= i; j++)
            {
                slae.add(qls[k].degree_of_freedom[i], qls[k].degree_of_freedom[j],
                         qls[k].get_local_G(i, j) * qls[k].ph->lambda);
                slae.add(qls[k].degree_of_freedom[i], qls[k].degree_of_freedom[j],
                         qls[k].get_local_M(i, j) * qls[k].ph->gamma);
            }

            slae.rp[qls[k].degree_of_freedom[i]] += qls[k].get_local_rp(i);
        }
    }
/*
    size_t n = slae.n;
    double ** A = new double * [n];
    for(size_t i = 0; i < n; i++)
    {
        A[i] = new double [n + 1];
        memset(A[i], 0, sizeof(double) * (n + 1));
    }
    for(size_t i = 0; i < n; i++)
    {
        for(size_t k = slae.ig[i]; k < slae.ig[i + 1]; k++)
        {
            A[i][slae.jg[k]] = slae.gg[k];
            A[slae.jg[k]][i] = slae.gg[k];
        }
        A[i][i] = slae.di[i];
        A[i][n] = slae.rp[i];
    }

    ofstream ofs("matrix.txt");
    for(size_t i = 0; i < 9; i++)
    {
        for(size_t j = 0; j < 9; j++)
            ofs << A[i][j] << " \t";
        ofs << endl;
    }
    ofs.flush();
    ofs.close();
*/
}

void FEM::applying_bounds()
{
    cout << "Applying bounds..." << endl;
    for(size_t bi = 0; bi < bounds_num; bi++)
    {
        if(get_type_of_bound(bounds[bi].gmsh_phys_num) == 2)
        {
            double theta[3] =
            {
                get_bound_value(bounds[bi].get_freedom_position(0), bounds[bi]),
                get_bound_value(bounds[bi].get_freedom_position(1), bounds[bi]),
                get_bound_value(bounds[bi].get_freedom_position(2), bounds[bi])
            };

            double lambda = bounds[bi].fes[0]->ph->lambda;
            for(size_t i = 0; i < 3; i++)
            {
                double vfrr = 0.0;
                for(size_t j = 0; j < 3; j++)
                    vfrr += bounds[bi].get_matrix_M(i, j) * theta[j] / lambda;
                slae.rp[bounds[bi].degree_of_freedom[i]] += vfrr;
            }
        }
    }

    for(size_t bi = 0; bi < bounds_num; bi++)
    {
        if(get_type_of_bound(bounds[bi].gmsh_phys_num) == 1)
        {
            for(size_t j = 0; j < 3; j++)
            {
                size_t freedom = bounds[bi].degree_of_freedom[j];
                double val = get_bound_value(bounds[bi].get_freedom_position(j), bounds[bi]);
                // Учет первых краевых "по-хорошему"
                // В диагональ пишем 1
                slae.di[freedom] = 1.0;
                // В правую часть пишем значение краевого
                slae.rp[freedom] = val;
                // А вот тут все веселье
                // Нам надо занулить строку, а у нас симметричная матрица
                // Поэтому будем бегать по матрице, занулять стоки
                // А то, что было в столбцах - выкидывать в правую часть
                size_t i_s = slae.ig[freedom], i_e = slae.ig[freedom + 1];
                for(size_t i = i_s; i < i_e; i++)
                {
                    slae.rp[slae.jg[i]] -= slae.gg[i] * val;
                    slae.gg[i] = 0.0;
                }
                for(size_t p = freedom + 1; p < degree_of_freedom_num; p++)
                {
                    size_t i_s = slae.ig[p], i_e = slae.ig[p + 1];
                    for(size_t i = i_s; i < i_e; i++)
                    {
                        if(slae.jg[i] == freedom)
                        {
                            slae.rp[p] -= slae.gg[i] * val;
                            slae.gg[i] = 0.0;
                        }
                    }
                }
            }
        }
    }
}

const quadrilateral * FEM::get_fe(const point & p) const
{
    for(size_t i = 0; i < qls_num; i++)
    {
        if(qls[i].inside(p))
            return qls + i;
    }
    cerr << "Warning: Target point " << p << " is outside of area" << endl;
    return NULL;
}

double FEM::get_solution(const point & p) const
{
    return get_solution(p, get_fe(p));
}

double FEM::get_solution(const point & p, const quadrilateral * fe) const
{
    if(fe)
    {
        double solution = 0.0;
        for(size_t i = 0; i < 9; i++)
            solution += fe->bfunc_2d(i, p) * slae.x[fe->degree_of_freedom[i]];
        return solution;
    }
    return 0.0;
}

point FEM::get_grad(const point & p) const
{
    return get_grad(p, get_fe(p));
}

point FEM::get_grad(const point & p, const quadrilateral * fe) const
{
    if(fe)
    {
        point solution(0.0, 0.0);
        for(size_t i = 0; i < 9; i++)
        {
            double q = slae.x[fe->degree_of_freedom[i]];
            point gbf = fe->grad_bfunc_2d(i, p);
            solution.x += q * gbf.x;
            solution.y += q * gbf.y;
        }
        return solution;
    }
    return point(0.0, 0.0);
}

// ============================================================================

vector<edge>::const_iterator FEM::get_edge_by_nodes(const point & begin, const point & end) const
{
    vector<edge>::const_iterator it;
    it = find(edges_freedom.begin(), edges_freedom.end(), edge(&begin, &end));
    if(it == edges_freedom.end())
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(&end, &begin));
    return it;
}

vector<edge>::const_iterator FEM::get_edge_by_fe(const quadrilateral & fe, size_t edge_num) const
{
    if(edge_num == 0) return get_edge_by_nodes(* fe.nodes[0], * fe.nodes[1]);
    if(edge_num == 1) return get_edge_by_nodes(* fe.nodes[0], * fe.nodes[2]);
    if(edge_num == 2) return get_edge_by_nodes(* fe.nodes[1], * fe.nodes[3]);
    if(edge_num == 3) return get_edge_by_nodes(* fe.nodes[2], * fe.nodes[3]);
    assert(edge_num < 4);
    return edges_freedom.end();
}

const point * FEM::get_nodes_by_edge(const edge & e, size_t node_num) const
{
    assert(node_num < 2);
    return e.nodes[node_num];
}

const quadrilateral * FEM::get_fe_by_edge(const edge & e, size_t fe_num) const
{
    assert(fe_num < 2);
    return e.fes[fe_num];
}

