import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from mesh_io import read_mesh, read_solution_table
from triangulate_quad_mesh import quads_to_tris


def build_u(nodes, sol, solution_type: str) -> np.ndarray:
    n = len(nodes)
    u = np.zeros(n)

    if solution_type == "from_file":
        missing = 0
        for i, (x, y) in enumerate(nodes):
            if i in sol:
                u_num, u_exact = sol[i]
                u[i] = u_num
            else:
                u[i] = x + y
                missing += 1
        if missing > 0:
            print(f"Предупреждение: для {missing} узлов нет численного решения, "
                  f"использовано u_exact = x + y.")
        return u

    # чисто аналитические функции
    for i, (x, y) in enumerate(nodes):
        if solution_type == "xy":
            u[i] = x + y
        elif solution_type == "r2":
            u[i] = x * x + y * y
        elif solution_type == "xy_prod":
            u[i] = x * y
        else:
            raise ValueError(f"Неизвестный тип решения: {solution_type}")

    return u


def main():
    mesh_file = "radial_mesh.txt"
    solution_file = "solution_lr3_test_r2.txt"

    if len(sys.argv) >= 2:
        solution_type = sys.argv[1]
    else:
        solution_type = "from_file" 

    nodes, quads = read_mesh(mesh_file)
    print(f"Узлов: {len(nodes)}, четырёхугольников: {len(quads)}")

    if solution_type == "from_file":
        sol = read_solution_table(solution_file)
        print("Размер словаря sol =", len(sol))
    else:
        sol = {}

    x = np.array([p[0] for p in nodes])
    y = np.array([p[1] for p in nodes])

    u = build_u(nodes, sol, solution_type)

    tris = quads_to_tris(quads)
    triangles = np.array(tris, dtype=int)
    triang = mtri.Triangulation(x, y, triangles)

    fig, ax = plt.subplots(figsize=(6, 5))

    tpc = ax.tripcolor(triang, u, shading="gouraud", cmap="viridis")
    fig.colorbar(tpc, ax=ax, label="u")

    n_levels = 10
    levels = np.linspace(u.min(), u.max(), n_levels)
    cs = ax.tricontour(triang, u, levels=levels, colors="k", linewidths=0.5)

    ax.triplot(triang, color="black", linewidth=0.2, alpha=0.5)

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()