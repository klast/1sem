#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    x0 = 1;
    xN = 2;
    u1 = 1;
    u2 = 3;
    first_task();
}

void MainWindow::first_task()
{
    QVector<double> a(N, -2), b(N, 1), c(N, 1), d(N, 0);
    c[0] = 0;
    b[N - 1] = 0;
    a[0] = 1;
    a[N - 1] = 1;
    b[0] = 0;
    c[N - 1] = 0;
    d[0] = u1;
    d[N - 1] = u2;
    u = tridiag(a, b, c, d);
    qDebug() << u;
    Graph();
}

void MainWindow::Graph()
{
    QChart * chart = new QChart();
    QLineSeries *series = new QLineSeries();
    double dx = (xN - x0) / double(N - 1);
    for(int i = 0; i < N; i++)
    {
        series->append(x0 + i * dx, u[i]);
    }
    chart->setTitle("Метод прогонки");
    chart->addSeries(series);

    QValueAxis *axisX = new QValueAxis;
    axisX->setLabelFormat("%g");
    axisX->setTitleText("X");
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis;
    axisY->setLabelFormat("%g");
    axisY->setTitleText("Y");
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    MainWindow::setCentralWidget(chartView);
}

QVector<double> MainWindow::tridiag(QVector<double> &a, QVector<double> &b, QVector<double> &c, QVector<double> &d)
{
    int n = a.size();
    QVector<double> p(n), q(n);
    p[0] = -b[0] / a[0];
    q[0] = d[0] / a[0];
    for(int i = 1; i < n; i++)
    {
        p[i] = -b[i] / (c[i] * p[i - 1] + a[i]);
        q[i] = (d[i] - c[i] * q[i - 1]) / (c[i] * p[i - 1] + a[i]);
    }
    QVector<double> u(n);  
    u[n - 1] = (d[n - 1] - c[n - 1] * q[n - 1]) / (a[n - 1] + c[n - 1] * p[n - 1]);
    for(int i = n - 2; i >= 0; i--)
    {
        u[i] = p[i] * u[i + 1] + q[i];
    }
    return u;
}

MainWindow::~MainWindow()
{
    delete ui;
}
