# Лабораторная работа №2  
Алгоритмы нумерации базисных функций в МКЭ  
(нумерация рёбер)

## Цель работы

Реализовать алгоритмы нумерации глобальных базисных функций в конечноэлементной сетке.

### Чтение сетки

```python
nodes, elements = read_mesh("radial_mesh.txt")
```

- `read_mesh(filename)`:
  - читает секции `NODES` и `ELEMENTS`;
  - возвращает:
    - `nodes: dict[int, (float, float)]` — словарь `id -> (x, y)`;
    - `elements: list[list[int]]` — каждая запись вида `[n1, n2, n3, n4, mat]`.

### Построение списка рёбер и таблиц связей

```python
edges, elem_edges_table, edge_dict, edge_to_elems = build_edges(elements)
```

Функция `build_edges(elements)`:

- перебирает все элементы;
- для каждого элемента берёт его 4 ребра:
  - `(n1, n2)`, `(n2, n3)`, `(n3, n4)`, `(n4, n1)`;
- каждое ребро хранится как **упорядоченная пара** узлов `key = (min(a,b), max(a,b))`,
  чтобы исключить дублирование;
- если ребро ещё не встречалось, оно добавляется в `edge_list` и получает новый номер `eid`;
- для каждого элемента запоминаются номера его рёбер.

Возвращаемые структуры:

- `edges: list[tuple[int, int]]` — глобальный список рёбер, `edges[eid] = (a, b)`;
- `elem_edges_table: list[list[int]]` — для каждого элемента `elem_id`:
  ```python
  elem_edges_table[elem_id] = [e1, e2, e3, e4]
  ```
- `edge_dict: dict[(int,int), int]` — словарь для быстрого поиска:
  ```python
  edge_dict[(min(a,b), max(a,b))] = edge_id
  ```
- `edge_to_elems: dict[int, list[int]]` — для каждого ребра `eid` хранит список элементов, которым оно принадлежит:
  ```python
  edge_to_elems[eid] = [elem_id1, elem_id2, ...]
  ```

### Вспомогательные функции

- `find_edge_id(edge_dict, a, b)`  
  По двум узлам `a` и `b` возвращает номер ребра или `None`, если ребро не найдено.

- `get_element_edges(elem_edges_table, elem_id)`  
  Возвращает список номеров рёбер элемента `elem_id`.

- `get_edge_info(edge_list, edge_to_elems, edge_id)`  
  Возвращает:
  ```python
  nodes, elems = get_edge_info(...)
  # nodes = (a, b) — узлы ребра
  # elems = [ ... ] — номера элементов, содержащих это ребро
  ```

### Визуализация сетки с номерами рёбер

```python
plot_edges(nodes, elements, edges)
```

### Вывод результатов в файл

```python
write_results_to_file(edges, elem_edges_table, edge_to_elems,
                      filename="lr2_output.txt")
```

Формат `lr2_output.txt`:

- сначала список всех рёбер и элементов, которым они принадлежат:
  ```text
  eid : a b : [elem_ids...]
  ```
- затем список рёбер по элементам:
  ```text
  elem_id: [e1, e2, e3, e4]
  ```
  
### Интерактивное меню

```python
menu(nodes, elements, edges, elem_edges_table, edge_dict, edge_to_elems)
```

1. **Номер ребра по двум узлам**  
   Вводятся узлы `A` и `B` → выводится номер ребра `eid` или сообщение, что такого ребра нет.  
   Используется `find_edge_id`.

2. **Номера рёбер по номеру элемента**  
   Вводится `elem_id` → выводится список `elem_edges_table[elem_id]`.

3. **Информация о ребре**  
   Вводится `eid` → выводятся:
   - его узлы `(a, b)`;
   - список элементов `edge_to_elems[eid]`.

4. **Визуализировать сетку с рёбрами**  
   Запускает `plot_edges`.

0. **Выход** — завершение программы.
