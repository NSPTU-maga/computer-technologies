#!/usr/bin/env python3
"""
ЛР3. Объединённый решатель для тестов в четверти круга на сетке radial_mesh.txt.

Поддерживаемые тесты (выбор через аргумент командной строки):

    1) lin_y  : u = y,   -Δu + u = y,
                дуга r=1  — Дирихле, оси — Нейман (g=-1 на Ox, g=0 на Oy).

    2) poly   : u = x^2 + y^2,
                -Δu + u = -4 + (x^2 + y^2),
                Дирихле на всей границе (дуга + оси).

    3) exp    : u = exp(0.1(x+y)),
                -Δu + u = (-0.02+1) * exp(0.1(x+y)),
                Дирихле на всей границе (дуга + оси).

Запуск примеров:
    python build_portrait_all.py lin_y radial_mesh.txt
    python build_portrait_all.py poly  radial_mesh.txt
    python build_portrait_all.py exp   radial_mesh.txt

Результаты:
    portrait_*.txt   — портрет матрицы,
    solution_*.txt   — таблица значений (u_num, u_exact, ошибки).
"""

import sys
import math
from collections import defaultdict
from typing import List, Sequence, Tuple, Callable, Dict

import numpy as np

# ---------- чтение сетки ----------

def read_mesh(filename: str):
    with open(filename, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    if lines[i] != "NODES":
        raise ValueError("Ожидалась секция NODES")
    i += 1

    n_nodes = int(lines[i])
    i += 1

    nodes: List[Tuple[float, float]] = []
    for _ in range(n_nodes):
        parts = lines[i].split()
        x = float(parts[1])
        y = float(parts[2])
        nodes.append((x, y))
        i += 1

    if lines[i] != "ELEMENTS":
        raise ValueError("Ожидалась секция ELEMENTS")
    i += 1

    n_elems = int(lines[i])
    i += 1

    elements: List[Tuple[List[int], int]] = []
    for _ in range(n_elems):
        parts = lines[i].split()
        n1 = int(parts[1])
        n2 = int(parts[2])
        n3 = int(parts[3])
        n4 = int(parts[4])
        mat = int(parts[5])
        elements.append(([n1, n2, n3, n4], mat))
        i += 1

    return n_nodes, nodes, elements

def build_portrait(n_nodes: int, elements):
    row_to_cols = defaultdict(set)

    for elem_nodes, _mat in elements:
        k = len(elem_nodes)
        for a in range(k):
            ia = elem_nodes[a]
            for b in range(a + 1, k):
                ib = elem_nodes[b]
                if ia == ib:
                    continue
                i_min = min(ia, ib)
                i_max = max(ia, ib)
                row_to_cols[i_min].add(i_max)

    ig = [0] * (n_nodes + 1)
    jg: List[int] = []

    nnz = 0
    for i in range(n_nodes):
        cols = row_to_cols.get(i, set())
        cols_sorted = sorted(cols)
        nnz += len(cols_sorted)
        ig[i + 1] = nnz
        jg.extend(cols_sorted)

    return ig, jg


def init_global_matrix(ig, jg):
    n = len(ig) - 1
    nnz = len(jg)

    di = [0.0] * n
    ggl = [0.0] * nnz
    ggu = [0.0] * nnz
    b  = [0.0] * n

    return di, ggl, ggu, b


def save_portrait(filename, ig, jg):
    with open(filename, "w", encoding="utf-8") as f:
        f.write(f"N = {len(ig) - 1}\n")
        f.write("ig:\n")
        f.write(" ".join(str(x) for x in ig) + "\n")
        f.write("jg:\n")
        f.write(" ".join(str(x) for x in jg) + "\n")


def add_local_matrix(
    di: List[float],
    ggl: List[float],
    ggu: List[float],
    ig: Sequence[int],
    jg: Sequence[int],
    L: Sequence[int],
    A_loc: Sequence[Sequence[float]],
):
    k = len(L)

    for a in range(k):
        i_glob = L[a]
        di[i_glob] += A_loc[a][a]

    for a in range(k):
        i_glob = L[a]
        row_start = ig[i_glob]
        row_end   = ig[i_glob + 1]
        for b in range(a + 1, k):
            j_glob = L[b]
            pos = -1
            for p in range(row_start, row_end):
                if jg[p] == j_glob:
                    pos = p
                    break
            if pos == -1:
                print(f"WARNING: пара ({i_glob}, {j_glob}) не найдена в портрете")
                continue
            ggu[pos] += A_loc[a][b]
            ggl[pos] += A_loc[b][a]

def _element_geometry(elem_nodes: Sequence[int],
                      nodes_coords: Sequence[Tuple[float, float]]):
    pts = [(nodes_coords[nid][0], nodes_coords[nid][1], nid) for nid in elem_nodes]

    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    def closest(pt_list, x0, y0):
        return min(pt_list, key=lambda p: (abs(p[0] - x0) + abs(p[1] - y0)))

    bl = closest(pts, min_x, min_y)
    br = closest(pts, max_x, min_y)
    tr = closest(pts, max_x, max_y)
    tl = closest(pts, min_x, max_y)

    pts_std = [bl, br, tr, tl]
    x_std = [p[0] for p in pts_std]
    y_std = [p[1] for p in pts_std]
    nid_std = [p[2] for p in pts_std]
    return x_std, y_std, nid_std


def compute_local_stiffness_quad(elem_nodes: Sequence[int],
                                 nodes_coords: Sequence[Tuple[float, float]],
                                 lam: float = 1.0) -> List[List[float]]:
    assert len(elem_nodes) == 4
    x_std, y_std, nid_std = _element_geometry(elem_nodes, nodes_coords)

    a = 1.0 / math.sqrt(3.0)
    gauss_pts = [(-a, -a), ( a, -a), ( a,  a), (-a,  a)]
    gauss_w   = [1.0, 1.0, 1.0, 1.0]

    K_std = [[0.0] * 4 for _ in range(4)]

    for (xi, eta), w in zip(gauss_pts, gauss_w):
        dN1_dxi  = -0.25 * (1 - eta)
        dN1_deta = -0.25 * (1 - xi)
        dN2_dxi  =  0.25 * (1 - eta)
        dN2_deta = -0.25 * (1 + xi)
        dN3_dxi  =  0.25 * (1 + eta)
        dN3_deta =  0.25 * (1 + xi)
        dN4_dxi  = -0.25 * (1 + eta)
        dN4_deta =  0.25 * (1 - xi)

        dN_dxi  = [dN1_dxi,  dN2_dxi,  dN3_dxi,  dN4_dxi]
        dN_deta = [dN1_deta, dN2_deta, dN3_deta, dN4_deta]

        dx_dxi = dx_deta = dy_dxi = dy_deta = 0.0
        for iN in range(4):
            dx_dxi  += dN_dxi[iN]  * x_std[iN]
            dx_deta += dN_deta[iN] * x_std[iN]
            dy_dxi  += dN_dxi[iN]  * y_std[iN]
            dy_deta += dN_deta[iN] * y_std[iN]

        J11 = dx_dxi
        J12 = dx_deta
        J21 = dy_dxi
        J22 = dy_deta
        detJ = J11 * J22 - J12 * J21

        invJ11 =  J22 / detJ
        invJ12 = -J12 / detJ
        invJ21 = -J21 / detJ
        invJ22 =  J11 / detJ

        dN_dx = [0.0] * 4
        dN_dy = [0.0] * 4
        for iN in range(4):
            dxi  = dN_dxi[iN]
            deta = dN_deta[iN]
            dN_dx[iN] = invJ11 * dxi + invJ12 * deta
            dN_dy[iN] = invJ21 * dxi + invJ22 * deta

        for iN in range(4):
            for jN in range(4):
                grad_dot = dN_dx[iN] * dN_dx[jN] + dN_dy[iN] * dN_dy[jN]
                K_std[iN][jN] += lam * grad_dot * detJ * w

    pos_in_std = {nid: i for i, nid in enumerate(nid_std)}
    K_elem = [[0.0] * 4 for _ in range(4)]
    for i_e, nid_i in enumerate(elem_nodes):
        i_s = pos_in_std[nid_i]
        for j_e, nid_j in enumerate(elem_nodes):
            j_s = pos_in_std[nid_j]
            K_elem[i_e][j_e] = K_std[i_s][j_s]
    return K_elem


def compute_local_mass_quad(elem_nodes: Sequence[int],
                            nodes_coords: Sequence[Tuple[float, float]]) -> List[List[float]]:
    assert len(elem_nodes) == 4
    x_std, y_std, nid_std = _element_geometry(elem_nodes, nodes_coords)

    a = 1.0 / math.sqrt(3.0)
    gauss_pts = [(-a, -a), ( a, -a), ( a,  a), (-a,  a)]
    gauss_w   = [1.0, 1.0, 1.0, 1.0]

    M_std = [[0.0] * 4 for _ in range(4)]

    for (xi, eta), w in zip(gauss_pts, gauss_w):
        N1 = 0.25 * (1 - xi) * (1 - eta)
        N2 = 0.25 * (1 + xi) * (1 - eta)
        N3 = 0.25 * (1 + xi) * (1 + eta)
        N4 = 0.25 * (1 - xi) * (1 + eta)
        N = [N1, N2, N3, N4]

        dN1_dxi  = -0.25 * (1 - eta)
        dN1_deta = -0.25 * (1 - xi)
        dN2_dxi  =  0.25 * (1 - eta)
        dN2_deta = -0.25 * (1 + xi)
        dN3_dxi  =  0.25 * (1 + eta)
        dN3_deta =  0.25 * (1 + xi)
        dN4_dxi  = -0.25 * (1 + eta)
        dN4_deta =  0.25 * (1 - xi)

        dN_dxi  = [dN1_dxi,  dN2_dxi,  dN3_dxi,  dN4_dxi]
        dN_deta = [dN1_deta, dN2_deta, dN3_deta, dN4_deta]

        dx_dxi = dx_deta = dy_dxi = dy_deta = 0.0
        for iN in range(4):
            dx_dxi  += dN_dxi[iN]  * x_std[iN]
            dx_deta += dN_deta[iN] * x_std[iN]
            dy_dxi  += dN_dxi[iN]  * y_std[iN]
            dy_deta += dN_deta[iN] * y_std[iN]

        J11 = dx_dxi
        J12 = dx_deta
        J21 = dy_dxi
        J22 = dy_deta
        detJ = J11 * J22 - J12 * J21

        for iN in range(4):
            for jN in range(4):
                M_std[iN][jN] += N[iN] * N[jN] * detJ * w

    pos_in_std = {nid: i for i, nid in enumerate(nid_std)}
    M_elem = [[0.0] * 4 for _ in range(4)]
    for i_e, nid_i in enumerate(elem_nodes):
        i_s = pos_in_std[nid_i]
        for j_e, nid_j in enumerate(elem_nodes):
            j_s = pos_in_std[nid_j]
            M_elem[i_e][j_e] = M_std[i_s][j_s]
    return M_elem


def compute_local_rhs_quad(elem_nodes: Sequence[int],
                           nodes_coords: Sequence[Tuple[float, float]],
                           f_rhs: Callable[[float, float], float]) -> List[float]:
    assert len(elem_nodes) == 4
    x_std, y_std, nid_std = _element_geometry(elem_nodes, nodes_coords)

    a = 1.0 / math.sqrt(3.0)
    gauss_pts = [(-a, -a), ( a, -a), ( a,  a), (-a,  a)]
    gauss_w   = [1.0, 1.0, 1.0, 1.0]

    b_std = [0.0] * 4

    for (xi, eta), w in zip(gauss_pts, gauss_w):
        N1 = 0.25 * (1 - xi) * (1 - eta)
        N2 = 0.25 * (1 + xi) * (1 - eta)
        N3 = 0.25 * (1 + xi) * (1 + eta)
        N4 = 0.25 * (1 - xi) * (1 + eta)
        N = [N1, N2, N3, N4]

        x = sum(N[i] * x_std[i] for i in range(4))
        y = sum(N[i] * y_std[i] for i in range(4))

        dN1_dxi  = -0.25 * (1 - eta)
        dN1_deta = -0.25 * (1 - xi)
        dN2_dxi  =  0.25 * (1 - eta)
        dN2_deta = -0.25 * (1 + xi)
        dN3_dxi  =  0.25 * (1 + eta)
        dN3_deta =  0.25 * (1 + xi)
        dN4_dxi  = -0.25 * (1 + eta)
        dN4_deta =  0.25 * (1 - xi)

        dN_dxi  = [dN1_dxi,  dN2_dxi,  dN3_dxi,  dN4_dxi]
        dN_deta = [dN1_deta, dN2_deta, dN3_deta, dN4_deta]

        dx_dxi = dx_deta = dy_dxi = dy_deta = 0.0
        for iN in range(4):
            dx_dxi  += dN_dxi[iN]  * x_std[iN]
            dx_deta += dN_deta[iN] * x_std[iN]
            dy_dxi  += dN_dxi[iN]  * y_std[iN]
            dy_deta += dN_deta[iN] * y_std[iN]

        J11 = dx_dxi
        J12 = dx_deta
        J21 = dy_dxi
        J22 = dy_deta
        detJ = J11 * J22 - J12 * J21

        f_val = f_rhs(x, y)

        for iN in range(4):
            b_std[iN] += f_val * N[iN] * detJ * w

    pts_std = list(zip(x_std, y_std, nid_std, range(4)))
    nid_to_std_index = {p[2]: p[3] for p in pts_std}

    b_elem = [0.0] * 4
    for i_e, nid in enumerate(elem_nodes):
        s_idx = nid_to_std_index[nid]
        b_elem[i_e] = b_std[s_idx]
    return b_elem

def build_masks_lin_y(nodes_coords, radius=1.0,
                      eps_axis=1e-8, eps_r=1e-6):
    n = len(nodes_coords)
    mask_D = [False] * n
    mask_N = [False] * n
    for i, (x, y) in enumerate(nodes_coords):
        r = math.hypot(x, y)
        on_arc    = abs(r - radius) < eps_r
        on_x_axis = abs(y) < eps_axis
        on_y_axis = abs(x) < eps_axis
        if on_arc:
            mask_D[i] = True
        elif on_x_axis or on_y_axis:
            mask_N[i] = True
    return mask_D, mask_N


def build_masks_dirichlet_all(nodes_coords, radius=1.0,
                              eps_axis=1e-8, eps_r=1e-6):
    n = len(nodes_coords)
    mask_D = [False] * n
    for i, (x, y) in enumerate(nodes_coords):
        r = math.hypot(x, y)
        on_arc = abs(r - radius) < eps_r
        on_x   = abs(y) < eps_axis
        on_y   = abs(x) < eps_axis
        if on_arc or on_x or on_y:
            mask_D[i] = True
    return mask_D, [False] * n

def add_neumann_edges(b: List[float],
                      elements,
                      nodes_coords,
                      mask_neumann,
                      g_on_edge: Callable[[float, float, float, float], float],
                      eps_axis=1e-8):
    for elem_nodes, _mat in elements:
        edges = [
            (elem_nodes[0], elem_nodes[1]),
            (elem_nodes[1], elem_nodes[2]),
            (elem_nodes[2], elem_nodes[3]),
            (elem_nodes[3], elem_nodes[0]),
        ]
        for ni, nj in edges:
            if not (mask_neumann[ni] and mask_neumann[nj]):
                continue
            x1, y1 = nodes_coords[ni]
            x2, y2 = nodes_coords[nj]

            g_val = g_on_edge(x1, y1, x2, y2)
            if g_val is None:
                continue

            length = math.hypot(x2 - x1, y2 - y1)
            if length == 0.0:
                continue

            a = 1.0 / math.sqrt(3.0)
            gauss_pts = [-a, a]
            gauss_w   = [1.0, 1.0]

            b1_loc = 0.0
            b2_loc = 0.0
            for ksi, w in zip(gauss_pts, gauss_w):
                phi1 = 0.5 * (1.0 - ksi)
                phi2 = 0.5 * (1.0 + ksi)
                b1_loc += g_val * phi1 * w
                b2_loc += g_val * phi2 * w

            J = length / 2.0
            b1_loc *= J
            b2_loc *= J

            b[ni] += b1_loc
            b[nj] += b2_loc


def apply_dirichlet(di, ggl, ggu, ig, jg, b,
                    nodes_coords,
                    mask_dirichlet,
                    u_exact: Callable[[float, float], float]):
    n = len(di)
    for i in range(n):
        if not mask_dirichlet[i]:
            continue
        x, y = nodes_coords[i]
        u_val = u_exact(x, y)

        for p in range(ig[i], ig[i + 1]):
            ggu[p] = 0.0
        for row in range(n):
            if row == i:
                continue
            for p in range(ig[row], ig[row + 1]):
                if jg[p] == i:
                    ggl[p] = 0.0

        di[i] = 1.0
        b[i]  = u_val

class TestDef:
    def __init__(self,
                 name: str,
                 u_exact: Callable[[float, float], float],
                 f_rhs: Callable[[float, float], float],
                 build_masks: Callable,
                 use_mass: bool,
                 use_neumann: bool,
                 g_on_edge: Callable[[float, float, float, float], float] | None):
        self.name = name
        self.u_exact = u_exact
        self.f_rhs = f_rhs
        self.build_masks = build_masks
        self.use_mass = use_mass
        self.use_neumann = use_neumann
        self.g_on_edge = g_on_edge


def make_tests() -> Dict[str, TestDef]:
    def u_lin_y(x, y): return y
    def f_lin_y(x, y): return y

    def g_lin_y(x1, y1, x2, y2):
        eps = 1e-8
        on_x = abs(y1) < eps and abs(y2) < eps
        on_y = abs(x1) < eps and abs(x2) < eps
        if on_x:
            return -1.0
        if on_y:
            return 0.0
        return None
    
    def u_poly(x, y): return x*x + y*y
    def f_poly(x, y): return -4.0 + (x*x + y*y)

    def u_exp(x, y): return math.exp(0.1 * (x + y))
    def f_exp(x, y): return (-0.02 + 1.0) * math.exp(0.1 * (x + y))

    return {
        "lin_y": TestDef("lin_y", u_lin_y, f_lin_y,
                         build_masks_lin_y, use_mass=True,
                         use_neumann=True, g_on_edge=g_lin_y),
        "poly":  TestDef("poly",  u_poly,  f_poly,
                         build_masks_dirichlet_all, use_mass=True,
                         use_neumann=False, g_on_edge=None),
        "exp":   TestDef("exp",   u_exp,   f_exp,
                         build_masks_dirichlet_all, use_mass=True,
                         use_neumann=False, g_on_edge=None),
    }


def main():

    test_name = sys.argv[1]
    mesh_file = sys.argv[2]

    tests = make_tests()

    test = tests[test_name]

    n_nodes, nodes, elements = read_mesh(mesh_file)

    ig, jg = build_portrait(n_nodes, elements)
    save_portrait(f"portrait_{test_name}.txt", ig, jg)
    di, ggl, ggu, b = init_global_matrix(ig, jg)

    for elem_nodes, _mat in elements:
        K = compute_local_stiffness_quad(elem_nodes, nodes, lam=1.0)
        A_loc = K
        if test.use_mass:
            M = compute_local_mass_quad(elem_nodes, nodes)
            A_loc = [[K[i][j] + M[i][j] for j in range(4)] for i in range(4)]

        b_loc = compute_local_rhs_quad(elem_nodes, nodes, test.f_rhs)

        L = sorted(elem_nodes)
        index_in_L = {node: idx for idx, node in enumerate(L)}
        k = len(elem_nodes)

        A_reordered = [[0.0] * k for _ in range(k)]
        for i_loc, node_i in enumerate(elem_nodes):
            i_new = index_in_L[node_i]
            for j_loc, node_j in enumerate(elem_nodes):
                j_new = index_in_L[node_j]
                A_reordered[i_new][j_new] = A_loc[i_loc][j_loc]

        b_reordered = [0.0] * k
        for i_loc, node_i in enumerate(elem_nodes):
            i_new = index_in_L[node_i]
            b_reordered[i_new] = b_loc[i_loc]

        add_local_matrix(di, ggl, ggu, ig, jg, L, A_reordered)
        for a in range(k):
            i_glob = L[a]
            b[i_glob] += b_reordered[a]

    mask_D, mask_N = test.build_masks(nodes)
    if test.use_neumann and test.g_on_edge is not None:
        add_neumann_edges(b, elements, nodes, mask_N, test.g_on_edge)
    apply_dirichlet(di, ggl, ggu, ig, jg, b, nodes, mask_D, test.u_exact)

    n = len(di)
    A = np.zeros((n, n))
    for i in range(n):
        A[i, i] = di[i]
        for p in range(ig[i], ig[i + 1]):
            j = jg[p]
            A[i, j] = ggu[p]
            A[j, i] = ggl[p]

    u = np.linalg.solve(A, np.array(b))
    max_err = 0.0

    out_name = f"solution_{test_name}.txt"
    with open(out_name, "w", encoding="utf-8") as fsol:
        fsol.write("i  x  y  u_num  u_exact  err  |err|\n")
        for i, (x, y) in enumerate(nodes):
            ue = test.u_exact(x, y)
            err = u[i] - ue
            aerr = abs(err)
            if aerr > max_err:
                max_err = aerr
            fsol.write(
                f"{i:4d}  {x: .6e}  {y: .6e}  {u[i]: .6e}  {ue: .6e}  {err: .3e}  {aerr: .3e}\n"
            )

    print(f" Максимальная погрешность = {max_err:.3e}")


if __name__ == "__main__":
    main()