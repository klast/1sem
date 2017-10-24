# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""
import numpy as np
import pandas as pd
import scipy.optimize as sc
import logging as logger

def rho(A,Za,u,ld,h,sigma):
    return np.linalg.norm(np.dot(A,Za) - u)**2 - (ld + sigma + h * np.linalg.norm(Za))**2

def my_solve(A_approx, U_approx, ld, h, sigma):
    B = np.dot(A_approx.T, u_approx)
    A_k = np.dot(A_approx.T, A_approx)
    alpha_0 = 1e-10
    alpha_n = 1000
    tmp = 1
    alpha = 1
    ZAlpha = np.linalg.solve(A_k + alpha_0 * np.eye(m), B)
    left_rho = rho(A_approx, ZAlpha, u_approx, ld, h, sigma)
    ZAlpha = np.linalg.solve(A_k + alpha_n * np.eye(m), B)
    right_rho = rho(A_approx, ZAlpha, u_approx, ld, h, sigma)
    
    while abs(right_rho) >= 1e-6 and abs(left_rho) >= 1e-6:
        alpha = (alpha_0 + alpha_n) / 2
        ZAlpha = np.linalg.solve(A_k + alpha * np.eye(m), B)
        mid_rho = rho(A_approx, ZAlpha, u_approx, ld, h, sigma)
        if np.sign(mid_rho) != np.sign(right_rho):
            left_rho = mid_rho
            alpha_0 = alpha
        elif np.sign(mid_rho) != np.sign(left_rho):
            right_rho = mid_rho
            alpha_n = alpha
        else:
            print("ERROR")
        tmp = rho(A_approx, ZAlpha, u_approx, ld, h, sigma)

    return [alpha, ZAlpha, tmp]

def LD(z):
    return np.linalg.norm(np.dot(A_approx, z) - u_approx) + h*np.linalg.norm(z) + sigma

n = 10
m = 10
bounds = [(0,2),(0,2), (0,2), (0,2), (0,2),(0,2),(0,2), (0,2), (0,2), (0,2)]

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
theta_u = np.random.uniform(-1 ,1, m)
epsilon = 1e-4
A_approx = np.multiply(A,1+theta_A*epsilon)
u_approx = np.multiply(u,1+theta_u*epsilon)
h = np.linalg.norm(A)*epsilon
sigma = np.linalg.norm(u)*epsilon
x0 = np.random.random(m)
#LD = lambda z: np.linalg.norm(np.dot(A_approx, z) - u_approx) + h*np.linalg.norm(z) + sigma
LD_val = sc.differential_evolution(LD, bounds)
print(LD(LD_val.x))
print(LD(z))
#LD_val = sc.minimize(LD, x0, method = 'nelder-mead', options={'xtol': 1e-8, 'disp': True}, args=(h, sigma, A_approx, u_approx))
print(np.linalg.norm(np.dot(A_approx, LD_val.x) - u_approx))
ans = my_solve(A_approx, u_approx, LD(LD_val.x), h, sigma)
res = np.linalg.norm(z-ans[1])
print(res)


