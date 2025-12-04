import math 
def print_square_coords(square_size: float, nx: int) -> None:
    h = square_size / (nx - 1)

    for iy in range(nx):
        y = iy * h
        row_parts = []
        for ix in range(nx):
            x = ix * h
            row_parts.append(f"{x:.7f} {y:.7f}")
        print("  ".join(row_parts))


if __name__ == "__main__":
    square_size =   math.sqrt(2)*1e-6
    nx = 10            

    print_square_coords(square_size, nx)