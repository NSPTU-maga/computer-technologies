import sys
from collections import defaultdict
from typing import List, Sequence, Tuple
import numpy as np
import math

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
        row_end = ig[i_glob + 1]

        for b in range(a + 1, k):
            j_glob = L[b]

            pos = -1
            for p in range(row_start, row_end):
                if jg[p] == j_glob:
                    pos = p
                    break

            ggu[pos] += A_loc[a][b]
            ggl[pos] += A_loc[b][a]


def compute_local_stiffness_quad(elem_nodes: Sequence[int],
                                 nodes_coords: Sequence[Tuple[float, float]],
                                 lam: float = 1.0) -> List[List[float]]:

    assert len(elem_nodes) == 4

    pts = []
    for nid in elem_nodes:
        x, y = nodes_coords[nid]
        pts.append((x, y, nid))

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

    a = 1.0 / math.sqrt(3.0)
    gauss_pts = [(-a, -a), (a, -a), (a, a), (-a, a)]
    gauss_w = [1.0, 1.0, 1.0, 1.0]

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

def build_boundary_masks(nodes_coords,
                         radius=1.0,
                         eps_axis=1e-8,
                         eps_r=1e-6):

    n = len(nodes_coords)
    mask_dirichlet = [False] * n
    mask_neumann   = [False] * n

    for i, (x, y) in enumerate(nodes_coords):
        r = (x * x + y * y) ** 0.5
        on_arc = abs(r - radius) < eps_r

        on_x_axis = abs(y) < eps_axis   
        on_y_axis = abs(x) < eps_axis  

        if on_arc:
            mask_dirichlet[i] = True
        elif on_x_axis or on_y_axis:
            mask_neumann[i] = True

    return mask_dirichlet, mask_neumann

def add_neumann_on_edge(b: List[float],
                        node_i: int,
                        node_j: int,
                        nodes_coords: Sequence[Tuple[float, float]],
                        g_value: float):

    x1, y1 = nodes_coords[node_i]
    x2, y2 = nodes_coords[node_j]

    length = math.hypot(x2 - x1, y2 - y1)
    if length == 0.0:
        return

    a = 1.0 / math.sqrt(3.0)
    gauss_pts = [-a, a]
    gauss_w   = [1.0, 1.0]

    b1_loc = 0.0
    b2_loc = 0.0
    for ksi, w in zip(gauss_pts, gauss_w):
        phi1 = 0.5 * (1.0 - ksi)
        phi2 = 0.5 * (1.0 + ksi)
        b1_loc += g_value * phi1 * w
        b2_loc += g_value * phi2 * w

    J = length / 2.0
    b1_loc *= J
    b2_loc *= J

    b[node_i] += b1_loc
    b[node_j] += b2_loc


def add_neumann_boundary_contributions(b: List[float],
                                       elements,
                                       nodes_coords,
                                       mask_neumann,
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

            on_x_axis = abs(y1) < eps_axis and abs(y2) < eps_axis
            on_y_axis = abs(x1) < eps_axis and abs(x2) < eps_axis

            if not (on_x_axis or on_y_axis):
                continue

            g_value = -1.0
            add_neumann_on_edge(b, ni, nj, nodes_coords, g_value)


def apply_dirichlet_conditions(di, ggl, ggu, ig, jg, b,
                               nodes_coords,
                               mask_dirichlet):
    n = len(di)
    for i in range(n):
        if not mask_dirichlet[i]:
            continue

        x, y = nodes_coords[i]
        u_val = x + y

        for p in range(ig[i], ig[i + 1]):
            ggu[p] = 0.0

        for row in range(n):
            if row == i:
                continue
            for p in range(ig[row], ig[row + 1]):
                if jg[p] == i:
                    ggl[p] = 0.0

        di[i] = 1.0
        b[i] = u_val

def main():
    mesh_file = sys.argv[1]

    n_nodes, nodes, elements = read_mesh(mesh_file)

    ig, jg = build_portrait(n_nodes, elements)
    save_portrait("portrait_lr3.txt", ig, jg)
    di, ggl, ggu, b = init_global_matrix(ig, jg)


    for elem_nodes, _mat in elements:
        K_geom = compute_local_stiffness_quad(elem_nodes, nodes, lam=1.0)

        L = sorted(elem_nodes)
        index_in_L = {node: idx for idx, node in enumerate(L)}
        k = len(elem_nodes)
        K_reordered = [[0.0] * k for _ in range(k)]
        for i_loc, node_i in enumerate(elem_nodes):
            i_new = index_in_L[node_i]
            for j_loc, node_j in enumerate(elem_nodes):
                j_new = index_in_L[node_j]
                K_reordered[i_new][j_new] = K_geom[i_loc][j_loc]

        add_local_matrix(di, ggl, ggu, ig, jg, L, K_reordered)

    mask_dirichlet, mask_neumann = build_boundary_masks(nodes, radius=1.0)
    add_neumann_boundary_contributions(b, elements, nodes, mask_neumann)
    apply_dirichlet_conditions(di, ggl, ggu, ig, jg, b, nodes, mask_dirichlet)

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
    for i, (x, y) in enumerate(nodes):
        u_exact = x + y
        err = abs(u[i] - u_exact)
        if err > max_err:
            max_err = err
    print(f"Максимальная ошибка по тесту u = x + y: {max_err:.3e}")

    with open("solution_lr3_test_linear.txt", "w", encoding="utf-8") as fsol:
        fsol.write("i  x  y  u_num  u_exact  |err|\n")
        for i in range (n):
            x, y = nodes[i]
            u_exact = x + y
            err = abs(u[i] - u_exact)
            fsol.write(f"{i:3d}  {x:.6f}  {y:.6f}  {u[i]:.6e}  {u_exact:.6e}  {err:.3e}\n")

if __name__ == "__main__":
    main()