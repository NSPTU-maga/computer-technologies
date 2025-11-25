#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "glwidget.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    // Конструктор
    explicit MainWindow(QWidget *parent = 0);
    // Деструктор
    ~MainWindow();
    // Событие при изменении размера окна
    void resizeEvent(QResizeEvent *event);
private slots:
    // Событие при переключении рисования сетки
    void on_checkBox_clicked();
    // Событие при переключении рисования изолиний
    void on_checkBox_2_clicked();
    // Событие при переключении закраски цветом
    void on_checkBox_3_clicked();
    // Событие при изменении числа изолиний
    void on_spinBox_valueChanged(int arg1);
    // Событие при изменении числа внутренних сегментов каждого КЭ
    void on_spinBox_2_valueChanged(int arg1);
private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
