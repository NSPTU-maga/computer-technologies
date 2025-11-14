import matplotlib.pyplot as plt

def read_mesh(filename):
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    nodes = {}
    elements = []

    i = 0
    if lines[i] != "NODES":
        raise ValueError("Ожидалась секция NODES")
    i += 1

    n_nodes = int(lines[i])
    i += 1

    for _ in range(n_nodes):
        parts = lines[i].split()
        nid = int(parts[0])
        x = float(parts[1])
        y = float(parts[2])
        nodes[nid] = (x, y)
        i += 1

    if lines[i] != "ELEMENTS":
        raise ValueError("Ожидалась секция ELEMENTS")
    i += 1

    n_elems = int(lines[i])
    i += 1

    for _ in range(n_elems):
        parts = lines[i].split()
        elem_id = int(parts[0])
        n1 = int(parts[1])
        n2 = int(parts[2])
        n3 = int(parts[3])
        n4 = int(parts[4])
        mat = int(parts[5])
        elements.append([n1, n2, n3, n4, mat])
        i += 1

    return nodes, elements

def build_edges(elements):
    edge_dict = {}         
    edge_list = []          
    elem_edges = []        
    edge_to_elems = {}    

    for elem_id, elem in enumerate(elements):
        n1, n2, n3, n4, _ = elem

        local_edges = [
            (n1, n2),
            (n2, n3),
            (n3, n4),
            (n4, n1)
        ]

        local_ids = []

        for a, b in local_edges:
            key = tuple(sorted((a, b)))   

            if key not in edge_dict:
                edge_dict[key] = len(edge_list)
                edge_list.append(key)
                edge_to_elems[edge_dict[key]] = []

            eid = edge_dict[key]
            local_ids.append(eid)
            edge_to_elems[eid].append(elem_id)

        elem_edges.append(local_ids)

    return edge_list, elem_edges, edge_dict, edge_to_elems

def find_edge_id(edge_dict, a, b):
    key = tuple(sorted((a, b)))
    return edge_dict.get(key, None)


def get_element_edges(elem_edges_table, elem_id):
    return elem_edges_table[elem_id]


def get_edge_info(edge_list, edge_to_elems, edge_id):
    nodes = edge_list[edge_id]
    elems = edge_to_elems[edge_id]
    return nodes, elems

def plot_edges(nodes, elements, edges):
    plt.figure(figsize=(10, 10))

    # элементы
    for elem in elements:
        n1, n2, n3, n4, _ = elem
        xs = [nodes[n][0] for n in (n1, n2, n3, n4, n1)]
        ys = [nodes[n][1] for n in (n1, n2, n3, n4, n1)]
        plt.plot(xs, ys, "b-", linewidth=1)

    # номера ребер
    for eid, (a, b) in enumerate(edges):
        x = (nodes[a][0] + nodes[b][0]) / 2
        y = (nodes[a][1] + nodes[b][1]) / 2
        plt.text(x, y, str(eid), fontsize=7, color="red")

    # узлы
    for nid, (x, y) in nodes.items():
        plt.plot(x, y, "ko", markersize=3)

    plt.gca().set_aspect("equal", adjustable="box")
    plt.grid(True, alpha=0.3)
    plt.show()

def menu(nodes, elements, edges, elem_edges_table, edge_dict, edge_to_elems):
    while True:
        print("\n=============")
        print("1 – Номер ребра по двум узлам")
        print("2 – Номера рёбер по номеру элемента")
        print("3 – Информация о ребре")
        print("4 – Визуализировать сетку с рёбрами")
        print("0 – Выход")

        cmd = input("Введите номер команды: ").strip()

        if cmd == "1":
            a = int(input("Узел A: "))
            b = int(input("Узел B: "))
            eid = find_edge_id(edge_dict, a, b)
            if eid is None:
                print("Такого ребра нет.")
            else:
                print(f"Ребро между узлами {a} и {b} имеет номер {eid}")

        elif cmd == "2":
            elem_id = int(input("Введите номер элемента: "))
            print("Рёбра элемента:", get_element_edges(elem_edges_table, elem_id))

        elif cmd == "3":
            e = int(input("Введите номер ребра: "))
            (a, b), elems = get_edge_info(edges, edge_to_elems, e)
            print(f"Ребро {e}: узлы ({a}, {b})")
            print("Принадлежит элементам:", elems)

        elif cmd == "4":
            plot_edges(nodes, elements, edges)

        elif cmd == "0":
            print("Выход.")
            break

        else:
            print("Неверная команда!")

if __name__ == "__main__":
    nodes, elements = read_mesh("radial_mesh.txt")
    edges, elem_edges_table, edge_dict, edge_to_elems = build_edges(elements)
    menu(nodes, elements, edges, elem_edges_table, edge_dict, edge_to_elems)
