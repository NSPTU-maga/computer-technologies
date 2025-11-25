from typing import List, Tuple


def quads_to_tris(elements: List[List[int]]) -> List[Tuple[int, int, int]]:
    """
    Преобразует список четырёхугольных элементов в список треугольников.
    Для каждого [n1, n2, n3, n4] создаём два треугольника:
      (n1, n2, n3) и (n1, n3, n4).

    Возвращает:
      tris: List[Tuple[int,int,int]]
    """
    tris: List[Tuple[int, int, int]] = []
    for quad in elements:
        if len(quad) != 4:
            raise ValueError(f"Ожидался четырёхугольник, получили {quad}")
        n1, n2, n3, n4 = quad
        tris.append((n1, n2, n3))
        tris.append((n1, n3, n4))
    return tris