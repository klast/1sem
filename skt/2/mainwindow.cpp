#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    run();
}

double sq_norm(QVector<double> &a, QVector<double> &b)
{
    double norm = 0;
    std::vector<double> tmp(a.size());
    std::transform(a.begin(), a.end(), b.begin(), tmp.begin(),
                   [](double &a1, double &a2){ return (a1 - a2) * (a1 - a2); });
    return sqrt(std::accumulate(tmp.begin(), tmp.end(),0));
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

void MainWindow::solve(QVector<double> &diag, QVector<double> &diag_up, QVector<double> &diag_down, QVector<double> &diag_left, QVector<double> &diag_right, QVector<double> &rhs)
{
    int n_2 = diag.size();
    int n = sqrt(n_2);
    MatrixXd a = MatrixXd::Zero(diag.size(), diag.size());
    for(int i = 0; i < diag.size(); i++)
    {
        a(i, i) = diag[i];
        if(i % (n - 1) != 0)
            a(i + 1, i) = diag_up[i];
        if(i % n != 0)
            a(i - 1, i) = diag_down[i];
    }
    for(int i = n; i < diag.size(); i++)
    {
        a(i - n, i) = diag_left[i - n];
        a(i, i - n) = diag_right[i];
    }
    qDebug() << "matrix filled";
    VectorXd rhs_solve = Map<VectorXd, Unaligned>(rhs.toStdVector().data(), rhs.toStdVector().size());

    VectorXd u_solve = a.colPivHouseholderQr().solve(rhs_solve);

    VectorXd::Map(&u[0], u_solve.size()) = u_solve;
}

void MainWindow::run()
{
    double k_1 = 110;
    double k_2 = 80;
    double L = 0.125;
    double ksi = 1.15;
    double s_l1 = 0.5, s_l2 = 0.75;
    double q_left = -3;
    double alpha_right = 0.13;
    double T_right = 412;

    double dt = 0.1;
    int n_time_steps = 10;

    double L1 = L * ksi / (1.0 + ksi);
    double L2 = L - L1;

    double s_L1 = s_l1 * L;
    double s_L2 = s_l2 * L;

    double start_u = T_right * alpha_right;
    int min_n = 1;

    QVector<Condition> conds;
    conds.push_back(Condition(0, ConditionType::metal));
    conds.push_back(Condition(s_L1, ConditionType::source));
    conds.push_back(Condition(s_L2, ConditionType::source));
    conds.push_back(Condition(L, ConditionType::metal));

    sort(conds.begin(), conds.end());

    int segment_ind = -1;
    double min_length = L + 1;
    for(int i = 0; i < conds.size() - 1; i++)
    {
        if (min_length > conds[i + 1].x() - conds[i].x() && conds[i + 1].x() - conds[i].x() != 0)
        {
            min_length = conds[i + 1].x() - conds[i].x();
            segment_ind = i;
        }
    }
    qDebug() << "Conditions:";
    for (int i = 0; i < conds.size(); i++)
    {
        conds[i].print();
    }

    qDebug() << "Shortest segment [" << conds[segment_ind].x() << "; " << conds[segment_ind + 1].x() << "]";

    QVector<double> dx(conds.size() - 1), dy(conds.size() - 1);
    QVector<double> x_start(conds.size() - 1), y_start(conds.size() - 1);
    QVector<int> nx(conds.size()), ny(conds.size());
    nx[0] = 0;
    calc_dx(dx[segment_ind], x_start[segment_ind], nx[segment_ind + 1], conds[segment_ind], conds[segment_ind + 1], 0.0, min_n);
    for (int i = segment_ind + 1; i < conds.size() - 1; i++)
        calc_dx(dx[i], x_start[i], nx[i + 1], conds[i], conds[i + 1], dx[i - 1]);
    for (int i = segment_ind - 1; i > -1; i--)
        calc_dx(dx[i], x_start[i], nx[i + 1], conds[i], conds[i + 1], dx[i + 1]);

    for (int i = 1; i < conds.size(); i++)
        nx[i] += nx[i - 1];

    qDebug() << "Segments n, dx:";
    for (int i = 0; i < conds.size() - 1; i++)
        qDebug() << "\t[" << conds[i].x() << "; " << conds[i + 1].x() << "]  dx = " << dx[i] << " n = " << nx[i + 1] - nx[i] << " first_x = " << x_start[i];

    std::copy(dx.begin(), dx.end(), dy.begin());
    std::copy(x_start.begin(), x_start.end(), y_start.begin());
    std::copy(nx.begin(), nx.end(), ny.begin());
    int num_x = nx[nx.size() - 1];
    int num_y = num_x;
    int n = num_x * num_y;
    QVector<double> u_old(n);
    u.resize(n);
    std::fill(u.begin(), u.end(), start_u);
    std::fill(u_old.begin(), u_old.end(), start_u);
    double pos_x = 0;
    double pos_y = 0;
    QVector<double> diag(n), diag_left(n), diag_right(n), diag_up(n), diag_down(n), rhs(n);
    //! считаем массив k, по ячейке определить k
    std::vector< std::vector<double> > k(num_x);
    int segment_ind_x = 0;
    int segment_ind_y = 0;
    double min_L = std::min(s_L1, s_L2);
    double max_L = std::max(s_L1, s_L2);
    double max_k = std::max(k_1, k_2);
    qDebug() << "Calc k matrix";
    for(int i = 0; i < num_x; i++)
    {
        std::vector<double> tmp(num_y, k_1);
        k[i] = tmp;
        if(i >= nx[segment_ind_x + 1])
            segment_ind_x++;
        for(int j = 0; j < num_y; j++)
        {
            if(j >= ny[segment_ind_y + 1])
                segment_ind_y++;
            pos_x = x_start[segment_ind_x] + dx[segment_ind_x] * (i - nx[segment_ind_x]);
            pos_y = y_start[segment_ind_y] + dy[segment_ind_y] * (j - ny[segment_ind_y]);
            if (pos_x <= min_L)
                k[i][j] = k_1;
            else if (pos_x >= max_L)
                k[i][j] = k_2;
            else
            {
                double x0 = s_L1;
                double y0 = 0;
                double x1 = s_L2;
                double y1 = L;
                double koeff = (y1 - y0) / (x1 - x0);
                double b = (y0 + y1 - koeff * (x0 + x1)) / 2.;
                double y_border = koeff * pos_x + b;
                if(y_border > pos_y)
                    k[i][j] = k_1;
                else
                    k[i][j] = k_2;
            }
        }
    }
    bool first_iter = true;
    qDebug() << "norm:";

    /*MatrixXd a(n, n);
    a = MatrixXd::Zero(n, n)*/;
    for(int t = 0; t < n_time_steps; t++)
    {
        double left_dx, right_dx, up_dy, down_dy;
        double left_k, right_k, up_k, down_k;
        double left, right, up, down, cur;
        int index = 0;
        std::copy(u.begin(), u.end(), u_old.begin());
        do
        {
            int i = 0;
            int j = 0;
            up_dy = dy[0];
            down_dy = dy[0];
            //! Проставим границы по [0;j]
            for(j = 1; j < num_y - 1; j++)
            {
                if(j >= ny[segment_ind_y + 1])
                    segment_ind_y++;
                if(j == ny[segment_ind_y])
                    down_dy = dy[segment_ind_y - 1];
                if(j == ny[segment_ind_y + 1])
                    up_dy = dy[segment_ind_y + 1];
                right_k = k[1][j];
                up_k = k[0][j + 1];
                down_k = k[0][j - 1];
                right_dx = dx[0];
                right = right_k / right_dx;
                up = up_k / up_dy;
                down = down_k / down_dy;
                index = i * num_y + j;
                diag[index] = right + up + down;
                diag_down[index] = -down;
                diag_up[index] = -up;
                diag_right[index] = -right;
                diag_left[index] = 0;
                rhs[index] = q_left / dx[0];
            }
            //! [0;0]
            index = 0;
            diag[0] = k[0][1] / dy[0] + k[1][0] / dx[0];
            diag_up[0] = - k[0][1] / dy[0];
            diag_right[0] = -k[1][0] / dx[0];
            diag_down[0] = 0;
            diag_left[0] = 0;
            rhs[0] = q_left / dx[0];

            //![0; num_y-1]
            index = num_y - 1;
            diag[index] = k[0][index-1] / dy[dy.size() - 1] + k[1][index] / dx[0];
            diag_down[index] = -k[0][index-1] / dy[dy.size() - 1];
            diag_right[index] = -k[1][index] / dx[0];
            diag_up[index] = 0;
            diag_left[index] = 0;
            rhs[index] = q_left / dx[0];

            i = num_x - 1;
            //! Проставим границы по [num_x-1;j]
            for(j = 1; j < num_y - 1; j++)
            {
                if(j >= ny[segment_ind_y + 1])
                    segment_ind_y++;
                if(j == ny[segment_ind_y])
                    down_dy = dy[segment_ind_y - 1];
                if(j == ny[segment_ind_y + 1])
                    up_dy = dy[segment_ind_y + 1];
                left_k = k[i - 1][j];
                up_k = k[i][j + 1];
                down_k = k[i][j - 1];
                left_dx = dx[dx.size() - 1];
                left = left_k / left_dx;
                up = up_k / up_dy;
                down = down_k / down_dy;
                index = i * num_y + j;
                diag[index] = left + up + down + alpha_right;
                diag_down[index] = -down;
                diag_up[index] = -up;
                diag_right[index] = 0;
                diag_left[index] = -left;
                rhs[index] = alpha_right * T_right;
            }
            //! [num_x-1;0]
            i = num_x - 1;
            index = i * num_y;
            diag[index] = k[i][1] / dy[0] + k[i - 1][0] / dx[dx.size() - 1];
            diag_up[index] = -k[i][1] / dy[0];
            diag_left[index] = -k[i - 1][0] / dx[dx.size() - 1];
            diag_down[index] = 0;
            diag_right[index] = 0;
            rhs[index] = alpha_right * T_right;

            //![num_x - 1; num_y-1]
            j = num_y - 1;
            index = (i) * num_y + j;
            diag[index] = k[i][j - 1] / dy[dy.size() - 1] + k[i - 1][j] / dx[dx.size() - 1];
            diag_down[index] = -k[i][j-1] / dy[dy.size() - 1];
            diag_left[index] = -k[i - 1][j] / dx[dx.size() - 1];
            diag_up[index] = 0;
            diag_right[index] = 0;
            rhs[0] = alpha_right * T_right;

            //! Теперь внутри
            left_dx = dx[0];
            right_dx = dx[0];
            for(i = 1; i < num_x - 1; i++)
            {
                if(i >= nx[segment_ind_x + 1])
                    segment_ind_x++;
                pos_x = x_start[segment_ind_x] + dx[segment_ind_x] * (i - nx[segment_ind_x]);
                if(i == nx[segment_ind_x])
                    left_dx = dx[segment_ind_x - 1];
                if(i == nx[segment_ind_x + 1])
                    right_dx = dx[segment_ind_x + 1];
                index = i * num_y + j;
                left_k = k[i-1][0];
                right_k = k[i+1][0];
                up_k = k[i][1];
                up_dy = dy[0];
                left = left_k / left_dx;
                right = right_k / right_dx;
                up = up_k / up_dy;
                //! граница по
                rhs[index] = 0;
                diag[index] = left + right + up;
                diag_up[index] = -up;
                diag_down[index] = 0;
                diag_left[index] = -left;
                diag_right[index] = -right;
                for(j = 1; j < num_y - 1; j++)
                {
                    up_k = k[i][j+1];
                    down_k = k[i][j-1];
                    if(j >= ny[segment_ind_y + 1])
                        segment_ind_y++;

                    pos_y = y_start[segment_ind_y] + dy[segment_ind_y] * (j - ny[segment_ind_y]);

                    if(j == ny[segment_ind_y])
                        down_dy = dy[segment_ind_y - 1];
                    if(j == ny[segment_ind_y + 1])
                        up_dy = dy[segment_ind_y + 1];
                    double s_dx = 0, s_dy = 0;

                    if(pos_x >= s_L1 && pos_x <= s_L2)
                        s_dx = (left_dx + right_dx) / 2.0;
                    if(pos_y >= s_L1 && pos_y <= s_L2)
                        s_dy = (down_dy + up_dy) / 2.0;

                    left = left_k / left_dx;
                    right = right_k / right_dx;
                    up = up_k / up_dy;
                    down = down_k / down_dy;
                    cur = ((right_dx + left_dx + up_dy + down_dy) / 4.0) / dt;

                    int index = i * num_y + j;
                    diag[index] = cur + left + right + up + down + s_dx * s_dy;
                    diag_left[index] = -left;
                    diag_right[index] = -right;
                    diag_up[index] = -up;
                    diag_down[index] = -down;
                    rhs[index] = cur * u_old[index] + 6 * s_dx * s_dy;
                }
                j = num_y - 1;
                index = i * num_y + j;
                //! граница по [i; num_y - 1]
                left_k = k[i-1][j];
                right_k = k[i+1][j];
                down_k = k[i][j-1];
                down_dy = dy[dy.size() - 1];
                left = left_k / left_dx;
                right = right_k / right_dx;
                down = down_k / down_dy;
                //! граница по
                rhs[index] = 0;
                diag[index] = left + right + down;
                diag_up[index] = 0;
                diag_down[index] = -down;
                diag_left[index] = -left;
                diag_right[index] = -right;
            }
            qDebug() << "Beginning to calc";
            solve(diag, diag_up, diag_down, diag_left, diag_right, rhs);
            qDebug() << "\t" << sq_norm(u, u_old);
        }while ( sq_norm(u, u_old) > 0.00001 );
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}
