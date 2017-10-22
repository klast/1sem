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
    B = np.dot(A_approx.conj().T, u_approx)
    A_k = np.dot(A_approx.conj().T, A_approx)
    alpha_0 = 1e-6
    alpha_n = 10000
    m = A_k.shape[1]
    tmp = 1
    alpha = 1
    while True:
        alpha = (alpha_0 + alpha_n) / 2
        A_k1 = A_k + np.eye(m) * alpha
        ZAlpha = A_k1 / B
        tmp = rho(A_approx, ZAlpha, u_approx, ld, h, sigma)
        if np.abs(tmp) < 0.00001:
            break
        if tmp > 0:
            alpha_n = alpha
        else:
            alpha_0 = alpha
        logger.log(logger.CRITICAL, alpha)
    return [alpha, ZAlpha, tmp]

def LD(z):
    return np.linalg.norm(np.dot(A_approx, z) - u_approx) + h*np.linalg.norm(z) + sigma

n = 10
m = 10

#1
A = np.random.random((n, m))
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
LD_val = sc.fmin(LD, x0,ftol=1e-6,xtol=1e-6, maxfun=100000)
print(np.linalg.norm(np.dot(A_approx, LD_val) - u_approx))
ans = my_solve(A_approx, u_approx, LD(LD_val), h, sigma)
res = np.linalg.norm(z-ans[1])
print(res)


