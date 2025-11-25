from typing import List, Tuple


def read_mesh(filename: str):
    """
    Чтение файла вида radial_mesh.txt:
      NODES
      <n_nodes>
      id x y
      ...
      ELEMENTS
      <n_elems>
      id n1 n2 n3 n4 mat

    Возвращает:
      nodes: List[Tuple[float,float]]  # координаты (x, y) по глобальному номеру
      elements: List[List[int]]        # список четырёхугольников по номерам узлов
    """
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
        # parts: id, x, y
        x = float(parts[1])
        y = float(parts[2])
        nodes.append((x, y))
        i += 1

    if lines[i] != "ELEMENTS":
        raise ValueError("Ожидалась секция ELEMENTS")
    i += 1

    n_elems = int(lines[i])
    i += 1

    elements: List[List[int]] = []
    for _ in range(n_elems):
        parts = lines[i].split()
        # parts: elem_id, n1, n2, n3, n4, mat
        n1 = int(parts[1])
        n2 = int(parts[2])
        n3 = int(parts[3])
        n4 = int(parts[4])
        elements.append([n1, n2, n3, n4])
        i += 1

    if len(nodes) != n_nodes:
        raise ValueError("Несовпадение числа узлов при чтении")

    if len(elements) != n_elems:
        raise ValueError("Несовпадение числа элементов при чтении")

    return nodes, elements


from typing import Dict, Tuple, List


def read_solution_table(filename: str) -> Dict[int, Tuple[float, float]]:
    """
    Читает файл solution_lr3_test_r2.txt в формате:

      i  x  y  u_num  u_exact  |err|

    Возвращает:
      sol[i] = (u_num, u_exact)
    """
    sol: Dict[int, Tuple[float, float]] = {}

    with open(filename, "r", encoding="utf-8") as f:
        header = f.readline()  # "i  x  y  u_num  u_exact  |err|"
        if not header:
            raise ValueError("Пустой файл решения")

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            # parts: i, x, y, u_num, u_exact, err
            i_node = int(parts[0])
            u_n = float(parts[3])
            u_e = float(parts[4])
            sol[i_node] = (u_n, u_e)

    print("Прочитано строк решения:", len(sol))
    return sol