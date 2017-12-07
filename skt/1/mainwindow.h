#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
#include <iostream>
#include <fstream>
QT_CHARTS_USE_NAMESPACE

using namespace std;


namespace Ui {
class MainWindow;
}

enum ConditionType { metal, source };

class Condition {
protected:
    double _x;
    ConditionType _type;
public:
    Condition() {
        _x = 0;
        _type = ConditionType::metal;
    }

    Condition(double x, ConditionType type) {
        _x = x;
        _type = type;
    }

    bool operator<(const Condition & cond) {
        return _x < cond._x;
    }

    double x() {
        return _x;
    }

    ConditionType condition_type() {
        return _type;
    }

    void print() {
        std::string type_s = "metal";
        if (_type == ConditionType::source)
            type_s = "source";
        qDebug() << type_s.c_str() << _x;
    }
};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    QVector<double> tridiag(QVector<double> &a, QVector<double> &b, QVector<double> &c, QVector<double> &d);
    void first_task();
    void Graph1();
    void Graph2();
    void second_task();


private:
    Ui::MainWindow *ui;
    const int M = 18;
    const int N = pow(30 - M, log(10 + M)) / 500;
    double x0, xN, u1, u2;
    QVector<double> u;
};

#endif // MAINWINDOW_H
