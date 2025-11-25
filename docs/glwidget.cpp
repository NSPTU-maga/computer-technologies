#if defined(HAVE_QT5)
#include <QtWidgets>
#else
#include <QtGui>
#endif
#include "glwidget.h"
#include <QFont>
#include <QString>

#if !defined(GL_MULTISAMPLE)
#define GL_MULTISAMPLE  0x809D
#endif

// Конструктор
glwidget::glwidget(QWidget* parent) : QGLWidget(parent)
{
    // Пусть будем рисовать все
    draw_mesh = true;
    draw_isolines = true;
    draw_color = true;

    // Решаем МКЭ задачу
    fem.input();
    fem.make_portrait();
    fem.assembling_global();
    fem.applying_bounds();
    fem.slae.solve(1e-16);

    // Ищем минимальные и максимальные значения координат
    max_x = min_x = fem.nodes[0].x;
    max_y = min_y = fem.nodes[0].y;
    for(size_t i = 1; i < fem.nodes_num; i++)
    {
        if(fem.nodes[i].x > max_x) max_x = fem.nodes[i].x;
        if(fem.nodes[i].y > max_y) max_y = fem.nodes[i].y;
        if(fem.nodes[i].x < min_x) min_x = fem.nodes[i].x;
        if(fem.nodes[i].y < min_y) min_y = fem.nodes[i].y;
    }
    size_x = max_x - min_x;
    size_y = max_y - min_y;
    min_x -= size_x * 0.01;
    max_x += size_x * 0.01;
    min_y -= size_y * 0.01;
    max_y += size_y * 0.01;

    // Поправляем значения мин / макс чтобы влазило в сетку
    adjustAxis(min_x, max_x, num_ticks_x);
    adjustAxis(min_y, max_y, num_ticks_y);
    size_x = max_x - min_x;
    size_y = max_y - min_y;

    // Ищем максимальные и минимальные значения, чтобы нормировать цвет
    max_u = min_u = fem.slae.x[0];
    for(size_t i = 1; i < fem.slae.n; i++)
    {
        if(fem.slae.x[i] > max_u) max_u = fem.slae.x[i];
        if(fem.slae.x[i] < min_u) min_u = fem.slae.x[i];
    }

    // Рассчитываем изолинии
    const size_t isolines_num = 25;
    set_isolines_num(isolines_num);

    // Шаги для разбиения по цветовым областям
    step_u_big = (max_u - min_u) / 4.0;
    step_u_small = step_u_big / 256.0;

    // Число разбиений КЭ на сегменты
    set_div_num(1);

    // Запускаем таймер отрисовки
    startTimer(500);
}

// Пересчет значений изолиний
void glwidget::set_isolines_num(size_t isolines_num)
{
    mtx.lock();
    isolines.clear();
    double isolines_step = (max_u - min_u) / (double)(isolines_num + 1);
    for(size_t i = 0; i < isolines_num; i++)
        isolines.insert(min_u + isolines_step * (double)(i + 1));
    mtx.unlock();
}

// Изменение количества сегментов, на которые разбивается каждый КЭ
void glwidget::set_div_num(size_t num)
{
    mtx.lock();

    vector<triangle> tmp1;
    vector<triangle> tmp2;
    // Посчитаем в локальных координатах
    tmp1.push_back(triangle(point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)));
    tmp1.push_back(triangle(point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)));
    for(size_t i = 0; i < num; i++)
    {
        for(size_t j = 0; j < tmp1.size(); j++)
        {
            point middles[3] =
            {
                point((tmp1[j].nodes[0].x + tmp1[j].nodes[1].x) / 2.0, (tmp1[j].nodes[0].y + tmp1[j].nodes[1].y) / 2.0),
                point((tmp1[j].nodes[0].x + tmp1[j].nodes[2].x) / 2.0, (tmp1[j].nodes[0].y + tmp1[j].nodes[2].y) / 2.0),
                point((tmp1[j].nodes[1].x + tmp1[j].nodes[2].x) / 2.0, (tmp1[j].nodes[1].y + tmp1[j].nodes[2].y) / 2.0)
            };
            tmp2.push_back(triangle(tmp1[j].nodes[0], middles[0], middles[1]));
            tmp2.push_back(triangle(middles[0], tmp1[j].nodes[1], middles[2]));
            tmp2.push_back(triangle(middles[1], middles[2], tmp1[j].nodes[2]));
            tmp2.push_back(triangle(middles[0], middles[2], middles[1]));
        }
        tmp2.swap(tmp1);
        tmp2.clear();
    }

    // Заполняем вектор из треугольников переводя координаты в глобальные и считая цвет
    triangles.clear();
    for(size_t i = 0; i < fem.qls_num; i++)
    {
        for(size_t j = 0; j < tmp1.size(); j++)
        {
            triangle tmp_tr;
            // Переводим координаты в глобальные
            for(size_t k = 0; k < 3; k++)
                tmp_tr.nodes[k] = fem.qls[i].to_global(tmp1[j].nodes[k]);
            // Занесем значение решения в узлах
            for(size_t k = 0; k < 3; k++)
                tmp_tr.solution[k] = fem.get_solution(tmp_tr.nodes[k], fem.qls + i);

            // Барицентр треугольника
            double cx = 0.0, cy = 0.0;
            for(size_t k = 0; k < 3; k++)
            {
                cx += tmp_tr.nodes[k].x;
                cy += tmp_tr.nodes[k].y;
            }
            point center(cx / 3.0, cy / 3.0);
            // Решение в барицентре
            double center_u = fem.get_solution(center, fem.qls + i);

            // Ищем цвет решения по алгоритму заливки радугой (Rainbow colormap)
            unsigned short r_color = 0, g_color = 0, b_color = 0;
            if(center_u > min_u + step_u_big * 3.0)
            {
                r_color = 255;
                g_color = 255 - (unsigned short)((center_u - (min_u + step_u_big * 3.0)) / step_u_small);
                b_color = 0;
            }
            else if(center_u > min_u + step_u_big * 2.0)
            {
                r_color = (unsigned short)((center_u - (min_u + step_u_big * 2.0)) / step_u_small);
                g_color = 255;
                b_color = 0;
            }
            else if(center_u > min_u + step_u_big)
            {
                unsigned short tmp = (unsigned short)((center_u - (min_u + step_u_big)) / step_u_small);
                r_color = 0;
                g_color = tmp;
                b_color = 255 - tmp;
            }
            else
            {
                unsigned short tmp = 76 - (unsigned short)((center_u - min_u) / (step_u_small * (255.0 / 76.0)));
                r_color = tmp;
                g_color = 0;
                b_color = 255 - tmp;
            }

            // Приглушаем кислотные цвета
            r_color = r_color * 3 / 4 + 64;
            g_color = g_color * 3 / 4 + 64;
            b_color = b_color * 3 / 4 + 64;

            // Задаем посчитанный цвет
            tmp_tr.color[0] = (double)r_color / 255.0;
            tmp_tr.color[1] = (double)g_color / 255.0;
            tmp_tr.color[2] = (double)b_color / 255.0;

            // И заносим в вектор
            triangles.push_back(tmp_tr);
        }
    }

    mtx.unlock();
}

// Подгонка осей под реальность и вычисление шагов координатной сетки
void glwidget::adjustAxis(double & min, double & max, size_t & numTicks)
{
    static const double axis_epsilon = 1.0 / 10000.0;
    if(max - min < axis_epsilon)
    {
        min -= 2.0 * axis_epsilon;
        max += 2.0 * axis_epsilon;
    }

    static const size_t MinTicks = 10;
    double grossStep = (max - min) / MinTicks;
    double step = pow(10, floor(log10(grossStep)));

    if (5 * step < grossStep)
        step *= 5;
    else if (2 * step < grossStep)
        step *= 2;

    numTicks = (size_t)(ceil(max / step) - floor(min / step));
    min = floor(min / step) * step;
    max = ceil(max / step) * step;
}

// Событие таймера
void glwidget::timerEvent(QTimerEvent *)
{
    updateGL();
}

// Инициализация сцены
void glwidget::initializeGL()
{
    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    glEnable(GL_CULL_FACE);
    glEnable(GL_MULTISAMPLE);
}

// Действие при изменении размеров виджета
void glwidget::resizeGL(int nWidth, int nHeight)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-0.06, 1.015, -0.06, 1.02, -10.0, 1.0);
    glViewport(0, 0,(GLint)nWidth, (GLint)nHeight);
}

// Отрисовка сцены
void glwidget::paintGL()
{
    if(!mtx.tryLock()) return;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Координатные оси
    glColor3d(0.0, 0.0, 0.0);
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    glVertex3d(0.0, -0.005, -0.1);
    glVertex3d(0.0, 1.005, -0.1);
    glVertex3d(-0.005, 0.0, -0.1);
    glVertex3d(1.005, 0.0, -0.1);
    glEnd();

    // Подписи осей
    QFont fnt_mono("Courier", 8);
    QFont fnt_serif("Times", 10);
    fnt_mono.setLetterSpacing(QFont::PercentageSpacing, 75.0);
    fnt_serif.setBold(true);
    renderText(0.99f, -0.04f, 0.0f, trUtf8("x"), fnt_serif);
    renderText(-0.05f, 0.99f, 0.0f, trUtf8("y"), fnt_serif);

    // Координатная сетка
    glColor3d(0.85, 0.85, 0.85);
    glLineWidth(1.0f);
    for(size_t i = 0; i <= num_ticks_x; i++)
    {
        double x = (double)i / (double)num_ticks_x;
        glBegin(GL_LINES);
        glVertex3d(x, -0.01, -0.2);
        glVertex3d(x, 1.0, -0.2);
        glEnd();
    }
    for(size_t i = 0; i <= num_ticks_y; i++)
    {
        double y = (double)i / (double)num_ticks_y;
        glBegin(GL_LINES);
        glVertex3d(-0.01, y, -0.2);
        glVertex3d(1.0, y, -0.2);
        glEnd();
    }

    // Отрисовка шкалы
    glColor3d(0.0, 0.0, 0.0);
    glLineWidth(2.0f);
    for(size_t i = 0; i < num_ticks_x; i++)
    {
        double x = (double)i / (double)num_ticks_x;
        double x_real = (double)(floor((x * size_x + min_x) * 10000.0 + 0.5)) / 10000.0;
        QString st = QString::number(x_real);
        renderText((float)x - 0.01f, -0.04f, 0.001f, st, fnt_mono);
    }
    for(size_t i = 0; i < num_ticks_y; i++)
    {
        double y = (double)i / (double)num_ticks_y;
        double y_real = (double)(floor((y * size_y + min_y) * 10000.0 + 0.5)) / 10000.0;
        QString st = QString::number(y_real);
        renderText(-0.05f, (float)y - 0.01f, 0.001f, st, fnt_mono);
    }

    // Отрисовка конечноэлементной сетки
    if(draw_mesh)
    {
        glColor3d(0.1, 0.1, 0.1);
        glLineWidth(1.0f);
        glBegin(GL_LINES);
        for(vector<edge>::const_iterator i = fem.edges_freedom.begin(); i != fem.edges_freedom.end(); i++)
        {
            glVertex3d((i->nodes[0]->x - min_x) / size_x, (i->nodes[0]->y - min_y) / size_y, 0.1);
            glVertex3d((i->nodes[1]->x - min_x) / size_x, (i->nodes[1]->y - min_y) / size_y, 0.1);
        }
        glEnd();
    }

    // Отрисовка всех треугольников
    for(size_t i = 0; i < triangles.size(); i++)
    {
        // Раскрашивать будем если запрошено сие, иначе зальем белым цветом
        if(draw_color)
            // Задаем посчитанный цвет
            glColor3dv(triangles[i].color);
        else
            // Задаем белый цвет
            glColor3d(1.0, 1.0, 1.0);

        // Рисуем
        glBegin(GL_TRIANGLES);
        for(size_t k = 0; k < 3; k++)
            glVertex3d((triangles[i].nodes[k].x - min_x) / size_x, (triangles[i].nodes[k].y - min_y) / size_y, 0.01);
        glEnd();
    }

    for(size_t i = 0; i < triangles.size(); i++)
    {
        // Изолинии рисуем только если оно нам надо
        if(draw_isolines)
        {
            // Теперь рисуем изолинии
            // Будем искать наименьшее значение, большее или равное решению
            // Если значения на разных концах ребра будут разными, значит изолиния проходит через это ребро
            set<double>::const_iterator segment_isol[3];
            for(size_t k = 0; k < 3; k++)
                segment_isol[k] = isolines.lower_bound(triangles[i].solution[k]);

            // А теперь нарисуем, согласно вышеприведенному условию
            glColor3d(0.0, 0.0, 0.0);
            glLineWidth(1.0f);
            glBegin(GL_LINE_STRIP);
            if(segment_isol[0] != segment_isol[1])
                glVertex3d(((triangles[i].nodes[1].x + triangles[i].nodes[0].x) * 0.5 - min_x) / size_x, ((triangles[i].nodes[1].y + triangles[i].nodes[0].y) * 0.5 - min_y) / size_y, 0.02);
            if(segment_isol[0] != segment_isol[2])
                glVertex3d(((triangles[i].nodes[2].x + triangles[i].nodes[0].x) * 0.5 - min_x) / size_x, ((triangles[i].nodes[2].y + triangles[i].nodes[0].y) * 0.5 - min_y) / size_y, 0.02);
            if(segment_isol[1] != segment_isol[2])
                glVertex3d(((triangles[i].nodes[2].x + triangles[i].nodes[1].x) * 0.5 - min_x) / size_x, ((triangles[i].nodes[2].y + triangles[i].nodes[1].y) * 0.5 - min_y) / size_y, 0.02);
            glEnd();
        }
    }

    mtx.unlock();
}
