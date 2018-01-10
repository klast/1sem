﻿import numpy as np
import scipy.optimize as sc
import matplotlib.pyplot as plt
import pandas as pd

np.set_printoptions(precision=7)

def rho(A, Za, u, ld):
    return np.linalg.norm(np.dot(A, Za) - u) ** 2 - (ld + sigma + h * np.linalg.norm(Za)) ** 2

def my_solve(A_approx, U_approx, ld, h, sigma):
    B = np.dot(A_approx.T, u_approx)
    A_k = np.dot(A_approx.T, A_approx)
    alpha_0 = 1e-10
    alpha_n = 1000
    alpha = 1
    ZAlpha = np.linalg.solve(A_k + alpha_0 * np.eye(n), B)
    left_rho = rho(A_approx, ZAlpha, u_approx, ld)
    ZAlpha = np.linalg.solve(A_k + alpha_n * np.eye(n), B)
    right_rho = rho(A_approx, ZAlpha, u_approx, ld)

    while abs(right_rho) >= 1e-8 and abs(left_rho) >= 1e-8:
        alpha = (alpha_0 + alpha_n) / 2
        ZAlpha = np.linalg.solve(A_k + alpha * np.eye(n), B)
        mid_rho = rho(A_approx, ZAlpha, u_approx, ld)
        if np.sign(mid_rho) != np.sign(right_rho):
            left_rho = mid_rho
            alpha_0 = alpha
        elif np.sign(mid_rho) != np.sign(left_rho):
            right_rho = mid_rho
            alpha_n = alpha
        else:
            print("ERROR")
    return [alpha, ZAlpha, mid_rho]

def LambdaDelta(z):
    return np.linalg.norm(np.dot(A_approx, z) - u_approx) + h * np.linalg.norm(z) + sigma

def LambdaDeltaMu(z):
    return np.linalg.norm(np.dot(A_approx, z) - u_approx)

def myplot(values, indexes, title, name):
    x = np.arange(10) + 1
    #indexses2 = indexes
    # indexses2 = [0, 1, 2, 3, 10]
    # indexses2 = [0, 4, 6, 8, 10]
    indexses2 = [0, 5, 7, 9, 10]
    for i in indexes:
        if i == (len(indexes) - 1):
            plt.plot(x, values[i], marker='o', linestyle="--")
            continue
        if i in indexses2:
            plt.plot(x, values[i])
    plt.title(title)
    plt.grid(True)
    legends = ['Ea=1e-4,Eu=1e-4', 'Ea=1e-4,Eu=1e-3', 'Ea=1e-4,Eu=1e-2', 'Ea=1e-4,Eu=1e-1', 'Ea=1e-3,Eu=1e-4',
               'Ea=1e-3,Eu =1e-3', 'Ea=1e-2,Eu =1e-4','Ea=1e-2,Eu=1e-2', 'Ea=1e-1,Eu=1e-4', 'Ea=1e-1,Eu=1e-1',
               'Точное решение']
    plt.legend([legends[i] for i in indexses2])
    plt.xlabel('i')
    plt.ylabel('Z_i')
    plt.show()
    plt.savefig(name)
    plt.clf()

m = 10
n = 10

bounds = []
for i in range(10):
    bounds.append((0, 10))

if __name__ == "__main__":
    if False:
        A = np.random.random((m, n))
        z = np.random.random((n, 1))
        while np.linalg.det(A) <= 1000:
            A += 0.1*np.eye(m)

        u = np.dot(A, z)
        np.savetxt('matrix.txt', A)
        np.savetxt('right_part.txt', u)
    else:
        A = np.loadtxt('matrix.txt', dtype=np.float, ndmin=2)
        u = np.loadtxt('right_part.txt', dtype=np.float, ndmin=1)
    z = np.linalg.solve(A, u)
    z_vec = [z for i in range(10)]
    theta_A = np.random.uniform(-1, 1, (n, m))
    theta_u = np.random.uniform(-1, 1, m)
    epsilon = []
    Zalpha = []
    number = 0
    for i in reversed(range(4)):
        for j in reversed(range(4)):
            if i == j or i == 3 or j == 3:
                epsilon.append([10 ** (-i - 1), 10 ** (-j - 1)])
    taskList = [1, 1, 1, 1]

    A0 = A.copy()
    u0 = u.copy()
    z0 = z.copy()

    # Р—Р°РґР°РЅРёРµ 1
    if taskList[0] == 1:
        print('Часть 1')
        for eps in epsilon:
            number += 1
            epsilon_A = eps[0]
            epsilon_u = eps[1]
            A_approx = np.multiply(A, 1 + theta_A * epsilon_A)
            u_approx = np.multiply(u, 1 + theta_u * epsilon_u)
            h = np.linalg.norm(A) * epsilon_A
            sigma = np.linalg.norm(u) * epsilon_u
            lambda_delta = sc.differential_evolution(LambdaDelta, bounds, tol=0.0001)
            norm_1 = (np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
            ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), h, sigma)
            alpha = ans[0]
            Zalpha.append(ans[1])
            res = np.linalg.norm(z - ans[1])
            norm_2 = np.linalg.norm(z - ans[1])
            print("в„–", number, "Ea =", epsilon_A, "Eu =", epsilon_u, "h =", h, "Sigma =", sigma,
                  "Lambda =", LambdaDelta(lambda_delta.x), "Norma1 =", norm_1, "Norma2 =", norm_2, 'alpha =', alpha)
        Zalpha.append(z)  
        z_all = [Zalpha[i] for i in range(11)] 
        myplot(z_all, np.arange(11), 'Task1', 'Task1.png')

    # Р—Р°РґР°РЅРёРµ 2
    if taskList[1] == 1:
        print('Часть 2')
        m = 12
        A = np.zeros((12, 10))
        u = np.zeros(12)
        for i in range(0, 10):
            A[i] = A0[i]
            u[i] = u0[i]
        A[10] = A[0] + A[9]
        A[11] = A[1] + A[8]
        u[10] = u[0] + u[9]
        u[11] = u[1] + u[8]
        theta_A = np.random.uniform(-1, 1, (m, n))
        theta_u = np.random.uniform(-1, 1, m)
        Zalpha = []
        number = 0
        for eps in epsilon:
            number += 1
            epsilon_A = eps[0]
            epsilon_u = eps[1]
            A_approx = np.multiply(A, 1 + theta_A * epsilon_A)
            u_approx = np.multiply(u, 1 + theta_u * epsilon_u)
            h = np.linalg.norm(A) * epsilon_A
            sigma = np.linalg.norm(u) * epsilon_u
            x0 = np.random.random(m)
            lambda_delta = sc.differential_evolution(LambdaDelta, bounds, tol=0.0001)
            norm_1 = (np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
            ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), h, sigma)
            alpha = ans[0]
            Zalpha.append(ans[1])
            res = np.linalg.norm(z - ans[1])
            norm_2 = np.linalg.norm(z - ans[1])
            mumu = lambda_delta = sc.differential_evolution(LambdaDeltaMu, bounds, tol=0.0001)
            mu = (np.linalg.norm(np.dot(A_approx, mumu.x) - u_approx))
            print(number, "Ea =", epsilon_A, "Eu =", epsilon_u, "h =", h, "Sigma =", sigma,
                  "Lambda =", LambdaDelta(lambda_delta.x), "Norma1 =", norm_1, "Norma2 =", norm_2, 'alpha =', alpha,
                  'mu =', mu)
        Zalpha.append(z)
        z_all = [Zalpha[i] for i in range(11)]
        myplot(z_all, np.arange(11), 'Task2', 'Task2.png')

    if taskList[2] == 1:
        print("Часть 3")

        m = 8
        A = np.zeros((8, 10))
        u = np.zeros(8)
        for i in range(0, 8):
            A[i] = A0[i]
            u[i] = u0[i]
        theta_A = np.random.uniform(-1, 1, (8, 10))
        theta_u = np.random.uniform(-1, 1, 8)
        Zalpha = []
        z = np.array([0.621239493461318, 0.425657429338214, 0.519183117500868, 0.300995461139938, 0.723035150122100,
                      0.207081883879460, 0.147838111664225, 1.00435100889945, 0.320579456514438, 0.422201736579000])
        number = 0
        for eps in epsilon:
            number += 1
            epsilon_A = eps[0]
            epsilon_u = eps[1]
            A_approx = np.multiply(A, 1 + theta_A * epsilon_A)
            u_approx = np.multiply(u, 1 + theta_u * epsilon_u)
            h = np.linalg.norm(A) * epsilon_A
            sigma = np.linalg.norm(u) * epsilon_u
            x0 = np.random.random(m)
            lambda_delta = sc.differential_evolution(LambdaDelta, bounds, tol=0.0001)
            norm_1 = (np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
            ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), h, sigma)
            alpha = ans[0]
            Zalpha.append(ans[1])
            res = np.linalg.norm(z - ans[1])
            norm_2 = np.linalg.norm(z - ans[1])
            print(number, "Ea =", epsilon_A, "Eu =", epsilon_u, "h =", h, "Sigma =", sigma,
                  "Lambda =", LambdaDelta(lambda_delta.x), "Norma1 =", norm_1, "Norma2 =", norm_2, 'alpha =', alpha)
        Zalpha.append(z)
        z_all = [Zalpha[i] for i in range(11)]
        myplot(z_all, np.arange(11), 'Task3', 'Task3.png')

    if taskList[3] == 1:
        print('Часть 4')
        m = 10
        A = np.zeros((10, 10))
        A = A0.copy()
        u = u0.copy()
        A[9] = A[8]
        Zalpha1 = []
        A_approx = A
        u_approx = u

        x0 = np.random.random(m)
        lambda_delta = sc.differential_evolution(LambdaDelta, bounds, tol=0.0001)
        norm_1 = (np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), 0, 0)
        alpha = ans[0]
        Zalpha1.append(ans[1])
        res = np.linalg.norm(z - ans[1])
        norm_2 = np.linalg.norm(z - ans[1])
        print("Lambda =", LambdaDelta(lambda_delta.x), "Norma1 =", norm_1, "Norma2 =", norm_2, 'alpha =', alpha)

        Zalpha1.append(z)
        z_all = [Zalpha1[i] for i in range(2)]
        x = np.arange(10)
        plt.plot(x, Zalpha1[0])
        plt.plot(x, Zalpha1[1], marker='o', linestyle="--")
        plt.title("Task4")
        plt.grid(True)
        legends = ['Arrpoximate', 'Exact']
        plt.legend(legends)
        plt.xlabel('i')
        plt.ylabel('Z_i')
        plt.show()
        plt.savefig('Task4.png')
        plt.clf()
    print('Р’С‹С‡РёСЃР»РµРЅРёСЏ РїСЂРѕС€Р»Рё СѓСЃРїРµС€РЅРѕ')