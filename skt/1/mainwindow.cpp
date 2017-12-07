#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>

#define DEBUG

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    x0 = 1;
    xN = 2;
    u1 = 1;
    u2 = 3;
    //first_task();
    second_task();
}

double sq_norm(QVector<double> &a, QVector<double> &b)
{
    double norm = 0;
    for (int i = 0; i < a.size(); i++)
        norm += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(norm);
}

void calc_dx(double & dx, double & x_start, int & n, Condition left, Condition right, double pref_dx, int pref_n = -1)
{
    double delta = right.x() - left.x();
    if (delta <= 0)
    {
        dx = 0;
        n = 0;
    }

    double temp_n;

    if (pref_n > 0)
    {
        temp_n = pref_n;
        n = pref_n;
    }
    else
    {
        temp_n = int(delta / pref_dx) + 1;
        n = int(delta / pref_dx) + 1;
    }

    if (left.condition_type() == ConditionType::source)
        temp_n += 0.5;
    if (right.condition_type() == ConditionType::source)
        temp_n += 0.5;

    dx = delta / (temp_n - 1);

    x_start = left.x();
    if (left.condition_type() == ConditionType::source)
        x_start += dx / 2.0;
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
    u.resize(N);
    u = tridiag(a, b, c, d);
    qDebug() << u;
    Graph1();
}

void MainWindow::Graph1()
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

void MainWindow::Graph2()
{
    QChart * chart = new QChart();
    QLineSeries *series = new QLineSeries();
    int size = u.size();
    double dx = (xN - x0) / double(size - 1);
    for(int i = 0; i < size; i++)
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

void MainWindow::second_task()
{
    double k_1 = 237;
    double k_2 = 174;
    double L = 0.37;
    double ksi = 0.75;
    double s_l1 = 0.4, s_l2 = 0.55;
    double left_cond = 450;
    double right_cond = 1;

    double L1 = L * ksi / (1.0 + ksi);
    double L2 = L - L1;

    double s_L1 = s_l1 * L;
    double s_L2 = s_l2 * L;

    double start_u = 450;
    int min_n = 100;

    QVector<Condition> conds;
    conds.push_back(Condition(0, ConditionType::metal));
    conds.push_back(Condition(L1, ConditionType::metal));
    conds.push_back(Condition(s_L1, ConditionType::source));
    conds.push_back(Condition(s_L2, ConditionType::source));
    conds.push_back(Condition(L, ConditionType::metal));

    sort(conds.begin(), conds.end());

    int segment_ind = -1;
    double min_length = L + 1;
    for (int i = 0; i < conds.size() - 1; i++)
        if (min_length > conds[i + 1].x() - conds[i].x() && conds[i + 1].x() - conds[i].x() != 0)
        {
            min_length = conds[i + 1].x() - conds[i].x();
            segment_ind = i;
        }

    qDebug() << "Conditions:";
    for (int i = 0; i < conds.size(); i++)
    {
        conds[i].print();
    }

    qDebug() << "Shortest segment [" << conds[segment_ind].x() << "; " << conds[segment_ind + 1].x() << "]" << "\n";

    QVector<double> dx(conds.size() - 1);
    QVector<double> x_start(conds.size() - 1);
    QVector<int> nx(conds.size());
    nx[0] = 0;
    calc_dx(dx[segment_ind], x_start[segment_ind], nx[segment_ind + 1], conds[segment_ind], conds[segment_ind + 1], 0.0, min_n);
    for (int i = segment_ind + 1; i < conds.size() - 1; i++)
        calc_dx(dx[i], x_start[i], nx[i + 1], conds[i], conds[i + 1], dx[i - 1]);
    for (int i = segment_ind - 1; i > -1; i--)
        calc_dx(dx[i], x_start[i], nx[i + 1], conds[i], conds[i + 1], dx[i + 1]);

    for (int i = 1; i < conds.size(); i++)
        nx[i] += nx[i - 1];

    qDebug() << "Segments n, dx:\n";
    for (int i = 0; i < conds.size() - 1; i++)
        qDebug() << "\t[" << conds[i].x() << "; " << conds[i + 1].x() << "]  dx = " << dx[i] << " n = " << nx[i + 1] - nx[i] << " first_x = " << x_start[i] << "\n";

    int n = nx[nx.size() - 1];

    QVector<double> u_old(n);
    u.resize(n);
    for (int i = 0; i < n; i++) {
        u[i] = start_u;
        u_old[i] = u[i];
    }

    QVector<double> diag(n);
    QVector<double> diag_up(n);
    QVector<double> diag_down(n);
    QVector<double> rhs(n);

    double pos = 0;
    bool first_iteration = true;
    qDebug() << "Norm:";
    do {
        for (int i = 0; i < n; i++)
            u_old[i] = u[i];

        pos = 0;
        segment_ind = 0;

        for (int i = 1; i < n - 1; i++) {
            if (i >= nx[segment_ind + 1])
                segment_ind++;

            pos = x_start[segment_ind] + dx[segment_ind] * (i - nx[segment_ind]);

            double left_k = k_1;
            double right_k = k_1;
            if (pos > L1)
                left_k = k_2;
            if (pos >= L1)
                right_k = k_2;

            double left_dx = dx[segment_ind];
            double right_dx = dx[segment_ind];
            if (i == nx[segment_ind])
                left_dx = dx[segment_ind - 1];
            if (i == nx[segment_ind + 1])
                right_dx = dx[segment_ind + 1];

            double s_dx = 0;
            if (pos >= s_L1 && pos <= s_L2)
                s_dx = (left_dx + right_dx) / 2.0;
#ifdef DEBUG
            if (first_iteration)
                cout << i << " pos " << pos << " k " << left_k << " " << right_k << " dx " << left_dx << " " << right_dx << " s " << s_dx << endl;
#endif
            diag[i] = left_k / left_dx + right_k / right_dx - (10 + 2 * u_old[i]) * s_dx;
            diag_up[i] = -right_k / right_dx;
            diag_down[i] = -left_k / left_dx;
            rhs[i] = (5 - u_old[i] * u_old[i]) * s_dx;
        }
        diag[0] = 1;
        rhs[0] = left_cond;
        diag_up[0] = 0;

        diag[n - 1] = k_2 / dx[dx.size() - 1];
        diag_down[n - 1] = -k_2 / dx[dx.size() - 1];
        rhs[n - 1] = right_cond;
        u = tridiag(diag, diag_up, diag_down, rhs);
        qDebug() << "\t" << sq_norm(u, u_old) << "\n";
        first_iteration = false;
    } while (sq_norm(u, u_old) > 0.00001);

    QChart *chart = new QChart();
    QLineSeries *series = new QLineSeries();
    pos = 0;

    double max_u = u[0];
    double min_u = u[0];
    for (int i = 0; i < n - 1; i++)
    {
        if (i >= nx[segment_ind + 1])
            segment_ind++;

        pos = x_start[segment_ind] + dx[segment_ind] * (i - nx[segment_ind]);

        series->append(pos, u[i]);

        if (max_u < u[i])
            max_u = u[i];
        if (min_u > u[i])
            min_u = u[i];
    }
    series->append(L, u[n - 1]);
    series->setName("Температура");

    double delta = max_u - min_u;
    min_u -= 0.01 * delta;
    max_u += 0.01 * delta;

    double vals[3] = { s_L1, s_L2, L1 };
    chart->setTitle("Распределение температуры в прутке");

    chart->addSeries(series);
    QVector<QLineSeries *> tmp(3);
    for (int i = 0; i < 3; i++)
    {
        tmp[i] = new QLineSeries();
        tmp[i]->append(vals[i], min_u);
        tmp[i]->append(vals[i], max_u);
    }
    for(int i = 0; i < 3; i++)
        chart->addSeries(tmp[i]);
    tmp[0]->setName("Левая граница источника");
    tmp[1]->setName("Правая граница источника");
    tmp[2]->setName("Граница металлов");


    QValueAxis *axisX = new QValueAxis;
    axisX->setLabelFormat("%g");
    axisX->setTitleText("X");
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);
    for(int i = 0; i < 3; i++)
        tmp[i]->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis;
    axisY->setLabelFormat("%g");
    axisY->setTitleText("Y");
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);
    for(int i = 0; i < 3; i++)
        tmp[i]->attachAxis(axisY);

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
