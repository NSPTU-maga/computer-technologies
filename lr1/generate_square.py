def print_square_coords(square_size: float, nx: int) -> None:
    h = square_size / (nx - 1)

    for iy in range(nx):
        y = iy * h
        row_parts = []
        for ix in range(nx):
            x = ix * h
            row_parts.append(f"{x:.2f} {y:.2f}")
        print("  ".join(row_parts))


if __name__ == "__main__":
    square_size = 0.4  
    nx = 5            

    print_square_coords(square_size, nx)