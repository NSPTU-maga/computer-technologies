#include "mainwindow.h"
#include "ui_mainwindow.h"

// Конструктор
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Установка минимальных размеров окна
    this->setMinimumHeight(400);
    this->setMinimumWidth(600);

    // Перемещение в центр экрана
    QPoint center = QApplication::desktop()->availableGeometry().center();
    QPoint corner = QApplication::desktop()->availableGeometry().topLeft();
    center.setX(center.x() - this->width() / 2);
    center.setY(center.y() - this->height() / 2);
    if(center.x() <= corner.x() || center.y() <= corner.y())
        this->move(corner);
    else
        this->move(center);

    // Начальные значения элементов управления
    ui->checkBox->setChecked(true);
    ui->checkBox_2->setChecked(true);
    ui->checkBox_3->setChecked(true);
    ui->spinBox->setMinimum(0);
    ui->spinBox->setMaximum(100);
    ui->spinBox->setValue(10);
    ui->spinBox_2->setMinimum(0);
    ui->spinBox_2->setMaximum(7);
    ui->spinBox_2->setValue(3);

    // Немного эстетства
    this->setWindowTitle(trUtf8("FEMA lab04"));

    // Передача начальных значений виджету
    ui->widget->draw_mesh = ui->checkBox->isChecked();
    ui->widget->draw_isolines = ui->checkBox_2->isChecked();
    ui->widget->draw_color = ui->checkBox_3->isChecked();
    ui->widget->set_isolines_num(ui->spinBox->value());
    ui->widget->set_div_num(ui->spinBox_2->value());
}

// Деструктор
MainWindow::~MainWindow()
{
    delete ui;
}

// Событие при изменении размера окна
void MainWindow::resizeEvent(QResizeEvent *event)
{
    QMainWindow::resizeEvent(event);
    // Подгонка размеров OpenGL виджета при изменении размеров окна
    QRect main = ui->centralwidget->geometry();
    QRect ogl = ui->widget->geometry();
    ui->widget->setGeometry(ogl.x(), ogl.y(), main.width() - ogl.x(), main.height() - ogl.y());
}

// Событие при переключении рисования сетки
void MainWindow::on_checkBox_clicked()
{
    ui->widget->draw_mesh = ui->checkBox->isChecked();
}

// Событие при переключении рисования изолиний
void MainWindow::on_checkBox_2_clicked()
{
    ui->widget->draw_isolines = ui->checkBox_2->isChecked();
}

// Событие при переключении закраски цветом
void MainWindow::on_checkBox_3_clicked()
{
    ui->widget->draw_color = ui->checkBox_3->isChecked();
}

// Событие при изменении числа изолиний
void MainWindow::on_spinBox_valueChanged(int arg1)
{
    if(arg1 >= ui->spinBox->minimum() && arg1 <= ui->spinBox->maximum())
        ui->widget->set_isolines_num(arg1);
}

// Событие при изменении числа внутренних сегментов каждого КЭ
void MainWindow::on_spinBox_2_valueChanged(int arg1)
{
    if(arg1 >= ui->spinBox_2->minimum() && arg1 <= ui->spinBox_2->maximum())
        ui->widget->set_div_num(arg1);
}
