import numpy as np
import numpy.linalg as ln
from numpy.linalg import norm
import matplotlib.pyplot as plt
import time
from scipy.linalg import hessenberg

np.set_printoptions(edgeitems=7, precision=3)
np.core.arrayprint._line_width = 200

def add_rectangle(actnum, min_i, min_j, max_i, max_j):
    i1 = max(1, min_i)
    i2 = min(actnum.shape[0] - 1, max_i)
    j1 = max(1, min_j)
    j2 = min(actnum.shape[1] - 1, max_j)
    actnum[i1:i2, j1:j2] = 1

def add_rhombus(actnum, center_ij, size):
    for i in range(-size, size + 1):
        for j in range(-size, size + 1):
            if abs(i) + abs(j) > size:
                continue
            x = center_ij[0] + i
            y = center_ij[1] + j
            if x < 0 or x >= actnum.shape[0] - 1 or y < 0 or y >= actnum.shape[1] - 1:
                continue
            actnum[x, y] = 1

def numerate_cells(actnum):
    m = actnum[actnum == 1].shape[0]
    actnum[actnum == 1] = np.arange(m)
    return m

def calc_a(actnum, m, dx):
    a = np.zeros((m, m))

    for i in range(1, actnum.shape[0] - 1):
        for j in range(1, actnum.shape[1] - 1):
            if actnum[i, j] == -1:
                continue
            ind = actnum[i, j]
            a[ind, ind] = -4
            if actnum[i + 1, j] != -1:
                a[ind, actnum[i + 1, j]] = 1
            if actnum[i - 1, j] != -1:
                a[ind, actnum[i - 1, j]] = 1
            if actnum[i, j - 1] != -1:
                a[ind, actnum[i, j - 1]] = 1
            if actnum[i, j + 1] != -1:
                a[ind, actnum[i, j + 1]] = 1
    return a / dx / dx

def calc_field(actnum, vec):
    m = len(vec)
    field = np.zeros(actnum.shape)
    for i in range(0, actnum.shape[0]):
        for j in range(0, actnum.shape[1]):
            if actnum[i, j] == -1:
                continue
            field[i, j] = vec[actnum[i, j]]
    return field

def calc_e_value(a, m, translate=False):
    eps = 1e-9
    x = np.ones(m)
    e = np.eye(m)

    l_curr = 0
    l_prev = 10000
    iter_num = 0
    while ln.norm(l_prev - l_curr) > eps:
        if translate:
            x = ln.solve(a + e * l_curr, x / norm(x))
        else:
            x = ln.solve(a, x / norm(x))
        l_prev = l_curr
        # отношение Рэлея
        l_curr = x.dot(a.dot(x)) / x.dot(x)
        iter_num += 1
    return l_curr, x, iter_num

def main():
    n = 10
    dx = 1.0 / (n - 1)
    actnum = -np.ones((n, n), dtype=np.int32)
    # add_rhombus(actnum, [14, 12], 10)
    add_rectangle(actnum, 1, 1, n - 1, n - 1)
    # add_rectangle(actnum, 13, 15, 19, 17)

    m = numerate_cells(actnum)
    #print(actnum)

    a = calc_a(actnum, m, dx)
    e_vals, e_vecs = ln.eigh(a)
    f1 = calc_field(actnum, e_vecs[:, 0])
    f2 = calc_field(actnum, e_vecs[:, -1])

    from mpl_toolkits.mplot3d import Axes3D

    mesh_x, mesh_y = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))

    # fig1 = plt.figure(figsize=(19, 9))
    # sub1 = fig1.add_subplot(121, projection='3d')
    # sub1.set_title("Lambda ="+  str(e_vals[0]))
    # sub1.plot_surface(mesh_x, mesh_y, f1)
    # sub2 = fig1.add_subplot(122, projection='3d')
    # sub2.set_title("Lambda =" + str(e_vals[-1]))
    # sub2.plot_surface(mesh_x, mesh_y, f2)

    # plt.show()

    # задание 5
    lambda_analytics = np.zeros((n - 1) * (n - 1))
    for i in range(1, n):
        for j in range(1, n):
            lambda_analytics[(i - 1) * (n - 1) + (j - 1)] = - np.pi**2 * (i*i + j*j)

    e_vals[::-1].sort()
    e_vecs[::-1].sort()
    lambda_analytics[::-1].sort()

    print("Задание 5: собственные значения, найденные аналитически и численно")
    print(lambda_analytics)
    print(e_vals)

    l_min, v_min, i_min = calc_e_value(a, m)
    print("Задание 6: минимальное по модулю собственное значение", l_min)

    k = -50
    a_k = a - np.eye(m) * k
    l_k, v_k, i_k = calc_e_value(a_k, m)
    l_k += k
    print("Задание 7: ближайшее по модулю к числу", k, "собственное значение равно", l_k)

    l_t, v_t, i_t = calc_e_value(a, m, True)
    print("Задание 8А: аналогично 7-му, но со смещением матрицы на очередное приближение собственного числа", l_t)

    a_tridiag = hessenberg(a)
    l_h, v_h, i_h = calc_e_value(a_tridiag, m)
    print("Задание 8Б: аналогично, но приводим к трехдиагональной", l_h)

    fig1 = plt.figure(figsize=(19, 9))
    sub1 = fig1.add_subplot(111)
    plt.plot(v_min)
    plt.plot(v_k)
    plt.plot(v_t)
    plt.plot(v_h)
    plt.legend(["Обратные итерации", "K = "+str(k), "Со смещением", "От трехдиагональной"])
    plt.show()

main()