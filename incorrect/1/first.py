# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""
import numpy as np
import pandas as pd
import scipy.optimize as sc

def rho(A, Za,u,ld):
    return np.linalg.norm(np.dot(A,Za) - u)**2 - (ld + sigma + h * np.linalg.norm(Za))**2

def my_solve(A_approx, U_approx, ld, h, sigma):
    B = np.dot(A_approx.T, u_approx)
    A_k = np.dot(A_approx.T, A_approx)
    alpha_0 = 1e-10
    alpha_n = 1000
    tmp = 1
    alpha = 1
    ZAlpha = np.linalg.solve(A_k + alpha_0 * np.eye(m), B)
    left_rho = rho(A_approx, ZAlpha, u_approx, ld)
    ZAlpha = np.linalg.solve(A_k + alpha_n * np.eye(m), B)
    right_rho = rho(A_approx, ZAlpha, u_approx, ld)
    
    while abs(right_rho) >= 1e-6 and abs(left_rho) >= 1e-6:
        alpha = (alpha_0 + alpha_n) / 2
        ZAlpha = np.linalg.solve(A_k + alpha * np.eye(m), B)
        mid_rho = rho(A_approx, ZAlpha, u_approx, ld)
        if np.sign(mid_rho) != np.sign(right_rho):
            left_rho = mid_rho
            alpha_0 = alpha
        elif np.sign(mid_rho) != np.sign(left_rho):
            right_rho = mid_rho
            alpha_n = alpha
        else:
            print("ERROR")
        tmp = rho(A_approx, ZAlpha, u_approx, ld)

    return [alpha, ZAlpha, mid_rho]

def LambdaDelta(z):
    return np.linalg.norm(np.dot(A_approx, z) - u_approx) + h*np.linalg.norm(z) + sigma

n = 10
m = 10
bounds = []
for i in range(10):
    bounds.append((0, 2))

if __name__ == "__main__":
#1
    A = np.random.random((n, m))
    if np.linalg.det(A) < 1000:
        A = A + 5 * np.eye(m)
#2
    z = np.random.random(n)
#3
    u = np.dot(A, z)
#4
    theta_A = np.random.uniform(-1, 1, (n, m))
    theta_u = np.random.uniform(-1, 1, m)
    min_value = []
    real_value = []
    h_vec = []
    sigma_vec = []
    norm_1 = []
    norm_2 = []
    epsilon = []
    for i in reversed(range(4)):
        for j in reversed(range(4)):
            if i==j or i==3 or j==3:
                epsilon.append([10**(-i-1), 10**(-j-1)])
    for eps in epsilon:
        epsilon_A = eps[0]
        epsilon_u = eps[1]
        A_approx = np.multiply(A, 1 + theta_A * epsilon_A)
        u_approx = np.multiply(u, 1 + theta_u * epsilon_u)
        h = np.linalg.norm(A) * epsilon_A
        h_vec.append(h)
        sigma = np.linalg.norm(u) * epsilon_u
        sigma_vec.append(sigma)
        x0 = np.random.random(m)
#LD = lambda z: np.linalg.norm(np.dot(A_approx, z) - u_approx) + h*np.linalg.norm(z) + sigma
        lambda_delta = sc.differential_evolution(LambdaDelta, bounds)
        min_value.append(LambdaDelta(lambda_delta.x))
        print(LambdaDelta(lambda_delta.x))
        real_value.append(LambdaDelta(z))
        print(LambdaDelta(z))
#LD_val = sc.minimize(LD, x0, method = 'nelder-mead', options={'xtol': 1e-8, 'disp': True}, args=(h, sigma, A_approx, u_approx))
        print(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        norm_1.append(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), h, sigma)
        res = np.linalg.norm(z-ans[1])
        print(res)
        norm_2.append(np.linalg.norm(z-ans[1]))
    df = pd.DataFrame({
            'epsilon_A' : [epsilon[i][0] for i in range(10)],
            'epsilon_u' : [epsilon[i][1] for i in range(10)],
            'min_value' : min_value,
            'min_real_value' : real_value,
            'h' : h_vec,
            'sigma' : sigma_vec,
            'norm_1' : norm_1,
            'norm_2' : norm_2
            })
    writer = pd.ExcelWriter('output.xlsx')
    df.to_excel(writer,'sheet1')
    writer.save()
    writer.close()

