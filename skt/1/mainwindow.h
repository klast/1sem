#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
QT_CHARTS_USE_NAMESPACE


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    QVector<double> tridiag(QVector<double> &a, QVector<double> &b, QVector<double> &c, QVector<double> &d);
    void first_task();
    void Graph();


private:
    Ui::MainWindow *ui;
    const int M = 18;
    const int N = pow(30 - M, log(10 + M)) / 500;
    double x0, xN, u1, u2;
    QVector<double> u;
};

#endif // MAINWINDOW_H
