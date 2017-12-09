#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>
#include "QCustomPlot/qcustomplot.h"

namespace Ui {
class MainWindow;
}

enum ConditionType { metal, source };

class Condition
{
protected:
    double _x;
    ConditionType _type;
public:
    Condition()
    {
        _x = 0;
        _type = ConditionType::metal;
    }
    Condition(double x, ConditionType type)
    {
        _x = x;
        _type = type;
    }
    bool operator<(const Condition & cond)
    {
        return _x < cond._x;
    }
    double x()
    {
        return _x;
    }
    ConditionType condition_type()
    {
        return _type;
    }
    void print()
    {
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
    void run();
    void solve(QVector<double> &diag, QVector<double> &diag_up, QVector<double> &diag_down, QVector<double> &diag_left, QVector<double> &diag_right, QVector<double> &rhs);
    void plot();
private slots:
    void make_step();
private:
    int t;
    Ui::MainWindow *ui;
    QVector<double> u;
    QCPColorMap *colorMap;
    QCPColorScale * colorScale;
};

#endif // MAINWINDOW_H
