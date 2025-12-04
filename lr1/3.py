import matplotlib.pyplot as plt
import math

class MeshData:
    def __init__(self):
        self.radius = 1.0
        self.square_size = 0.7
        self.Kx = 0
        self.Ky = 0
        self.coord_lines = []  
        self.subregions = []
        self.nodes = []        
        self.elements = []    
        self.n_radial_div = 5 
        self.spacing = 1.1     

class Subregion:
    def __init__(self):
        self.material = 0
        self.nxb = self.nxe = 0
        self.nyb = self.nye = 0

def read_area_file(filename):
    mesh = MeshData()
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]

    cur = 0
    if 'radius' in lines[cur]:
        params = lines[cur].split()
        for i in range(0, len(params), 2):
            if params[i] == 'radius':
                mesh.radius = float(params[i+1])
            elif params[i] == 'square_size':
                mesh.square_size = float(params[i+1])
        cur += 1

    mesh.Kx, mesh.Ky = map(int, lines[cur].split())
    cur += 1

    mesh.coord_lines = []
    for i in range(mesh.Ky):
        coords = list(map(float, lines[cur].split()))
        line_points = []
        for j in range(0, len(coords), 2):
            line_points.append([coords[j], coords[j+1]])
        mesh.coord_lines.append(line_points)  
        cur += 1

    num_sub = int(lines[cur]); cur += 1
    mesh.subregions = []
    for _ in range(num_sub):
        data = list(map(int, lines[cur].split()))
        sr = Subregion()
        sr.material = data[0]
        sr.nxb, sr.nxe, sr.nyb, sr.nye = data[1], data[2], data[3], data[4]
        mesh.subregions.append(sr)
        cur += 1
    return mesh

def merge_duplicate_nodes(mesh, eps=1e-12):
    unique = {}
    new_nodes = []
    map_idx = {}
    for old_id, (x, y) in enumerate(mesh.nodes):
        key = (round(x/eps), round(y/eps))
        if key not in unique:
            unique[key] = len(new_nodes)
            new_nodes.append([x, y])
        map_idx[old_id] = unique[key]

    new_elements = []
    for n1,n2,n3,n4,mat in mesh.elements:
        new_elements.append([map_idx[n1], map_idx[n2], map_idx[n3], map_idx[n4], mat])

    mesh.nodes = new_nodes
    mesh.elements = new_elements

def generate_radial_mesh(mesh):
    s = mesh.square_size
    R = mesh.radius
    Kx, Ky = mesh.Kx, mesh.Ky
    L = mesh.n_radial_div      
    assert L >= 2, "n_radial_div должно быть >=2"

    nodes = []
    for i in range(Ky):
        row = mesh.coord_lines[i]
        row_sorted = sorted(row, key=lambda p: p[0])
        nodes.extend(row_sorted)
    base_count = len(nodes)  

    right_sq = [mesh.coord_lines[i][-1] for i in range(Ky)]        
    top_sq   = mesh.coord_lines[-1]                                

    arc_right = [] 
    for i in range(Ky):
        y = right_sq[i][1]
        ratio = y / s                      
        ang = (math.pi/4) * ratio           
        arc_right.append([R*math.cos(ang), R*math.sin(ang)])

    arc_top = []  
    for j in range(Kx):
        x = top_sq[j][0]
        ratio = 1 - x/s                 
        ang = (math.pi/4) + (math.pi/4)*ratio 
        arc_top.append([R*math.cos(ang), R*math.sin(ang)])

    def interp_t(layer_idx):
        if mesh.spacing == 1.0:
            return layer_idx/(L-1)
        total = (mesh.spacing**(L-1) - 1)/(mesh.spacing - 1)
        l = (mesh.spacing**layer_idx - 1)/(mesh.spacing - 1)
        return l/total

    right_layers = [] 
    top_layers   = [] 

    for layer_idx in range(1, L-1):
        t = interp_t(layer_idx)

        rlayer = []
        for i in range(Ky):
            x0,y0 = right_sq[i]
            xr,yr = arc_right[i]
            rlayer.append([x0 + (xr-x0)*t, y0 + (yr-y0)*t])
        right_layers.append(rlayer) 

        tlayer = []
        for j in range(Kx):
            x0,y0 = top_sq[j]
            xr,yr = arc_top[j]
            tlayer.append([x0 + (xr-x0)*t, y0 + (yr-y0)*t])
        top_layers.append(tlayer)   

    idx_right_layer = []  
    idx_top_layer   = []  

    for rlayer in right_layers:
        start = len(nodes)
        nodes.extend(rlayer)
        idx_right_layer.append(list(range(start, start+Ky)))

    for tlayer in top_layers:
        start = len(nodes)
        nodes.extend(tlayer)
        idx_top_layer.append(list(range(start, start+Kx)))

    start_arc_right = len(nodes)
    nodes.extend(arc_right)                   
    idx_arc_right = list(range(start_arc_right, start_arc_right+Ky))

    start_arc_top = len(nodes)
    nodes.extend(arc_top)                      
    idx_arc_top = list(range(start_arc_top, start_arc_top+Kx))

    mesh.nodes = nodes

    elements = []

    for i in range(Ky-1):
        for j in range(Kx-1):
            n1 = i*Kx + j
            n2 = i*Kx + j + 1
            n3 = (i+1)*Kx + j + 1
            n4 = (i+1)*Kx + j
            elements.append([n1,n2,n3,n4, 1])

    col_sq_right = [i*Kx + (Kx-1) for i in range(Ky)] 

    def right_col(layer_step):  
        if layer_step == 0:
            return col_sq_right
        if layer_step == (L-1):
            return idx_arc_right
        return idx_right_layer[layer_step-1]

    for s in range(L-1):
        lower = right_col(s)
        upper = right_col(s+1)
        for i in range(Ky-1):
            n1 = lower[i]
            n2 = lower[i+1]
            n3 = upper[i+1]
            n4 = upper[i]
            elements.append([n1,n2,n3,n4, 2])

    row_sq_top = [(Ky-1)*Kx + j for j in range(Kx)]  

    def top_row(layer_step):   # 0..(L-1)
        if layer_step == 0:
            return row_sq_top
        if layer_step == (L-1):
            return idx_arc_top
        return idx_top_layer[layer_step-1]

    for s in range(L-1):
        lower = top_row(s)
        upper = top_row(s+1)
        for j in range(Kx-1):
            n1 = lower[j]
            n2 = lower[j+1]
            n3 = upper[j+1]
            n4 = upper[j]
            elements.append([n1,n2,n3,n4, 2])

    mesh.elements = elements

    merge_duplicate_nodes(mesh)

    return mesh

def visualize_radial_mesh(mesh, show_node_ids=True, show_elem_ids=False):
    plt.figure(figsize=(8, 8))

    for eid, (n1,n2,n3,n4,_) in enumerate(mesh.elements):
        x = [mesh.nodes[n][0] for n in (n1,n2,n3,n4,n1)]
        y = [mesh.nodes[n][1] for n in (n1,n2,n3,n4,n1)]
        plt.plot(x, y, 'b-', linewidth=1)

        if show_elem_ids:
            cx = sum(x[:-1])/4
            cy = sum(y[:-1])/4
            plt.text(cx, cy, f"E{eid}", color='blue', fontsize=7, ha='center', va='center')

    X = [p[0] for p in mesh.nodes]
    Y = [p[1] for p in mesh.nodes]
    plt.plot(X, Y, 'ro', ms=3)

    if show_node_ids:
        for i,(x,y) in enumerate(mesh.nodes):
            plt.text(x, y, str(i), fontsize=7, color='black', ha='left', va='bottom')

    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend([])
    plt.tight_layout()
    plt.show()

def save_mesh_to_file(mesh, filename):
    with open(filename, 'w') as f:
        f.write("NODES\n")
        f.write(f"{len(mesh.nodes)}\n")
        for i,(x,y) in enumerate(mesh.nodes):
            f.write(f"{i} {x:.6f} {y:.6f}\n")
        f.write("\nELEMENTS\n")
        f.write(f"{len(mesh.elements)}\n")
        for i,(n1,n2,n3,n4,mat) in enumerate(mesh.elements):
            f.write(f"{i} {n1} {n2} {n3} {n4} {mat}\n")


if __name__ == "__main__":
    mesh = read_area_file("computational_domain3.txt")
    mesh = generate_radial_mesh(mesh)
    save_mesh_to_file(mesh, "radial_mesh.txt")
    visualize_radial_mesh(mesh, show_node_ids=True, show_elem_ids=False)
