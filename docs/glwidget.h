#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QString>
#include <QtOpenGL>
#include <QMutex>
#include "fem.h"

// Класс треугольник
class triangle
{
public:
    point nodes[3];
    double color[3];
    double solution[3];
    triangle() {}
    triangle(const point & node1, const point & node2, const point & node3)
    {
        nodes[0] = node1;
        nodes[1] = node2;
        nodes[2] = node3;
    }
};

// Класс OpenGL виджет
class glwidget : public QGLWidget
{
    Q_OBJECT
protected:
    // Событие таймера
    void timerEvent(QTimerEvent *);
public:
    // Инициализация сцены
    void initializeGL();
    // Действие при изменении размеров виджета
    void resizeGL(int nWidth, int nHeight);
    // Отрисовка сцены
    void paintGL();
    // Конструктор
    glwidget(QWidget* parent = 0);

    // Флаг отрисовки конечноэлементной сетки
    bool draw_mesh;
    // Флаг отрисовки изолиний
    bool draw_isolines;
    // Флаг закраски цветом
    bool draw_color;
    // Пересчет значений изолиний
    void set_isolines_num(size_t isolines_num);
    // Изменение количества сегментов, на которые разбивается каждый КЭ
    void set_div_num(size_t num);
private:
    // Мьютекс чтобы не случилось странностей при изменении данных извне
    QMutex mtx;
    // Класс МКЭ
    FEM fem;

    // Минимальные и максимальные значения геометрии + размер
    double min_x, max_x, size_x;
    double min_y, max_y, size_y;
    // Количество шагов координатной сетки
    size_t num_ticks_x, num_ticks_y;
    // Подгонка осей под реальность и вычисление шагов координатной сетки
    void adjustAxis(double & min, double & max, size_t & numTicks);

    // Минимальное и максимальное значения решения
    double min_u, max_u;
    // Значения изолиний
    set<double> isolines;
    // Вспомогательные шаги по цвету для закраски
    double step_u_big, step_u_small;

    // Треугольники, которые будем рисовать
    vector<triangle> triangles;
};

#endif // GLWIDGET_H
