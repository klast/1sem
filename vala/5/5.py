# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 00:20:24 2017

@author: СпелеВова
"""

import scipy as sc
import numpy as np
import scipy.sparse.linalg as spla
import matplotlib.pylab as plt
import time

#входная матрица
input_matrix = sc.io.mmread('0.4solve.mtx')
#правая часть
rhs = sc.io.mmread('0.rhs')
nx = rhs.shape[0]
matrix = input_matrix.tocsc()
L = sc.sparse.tril(matrix).tocsr()
diag = matrix.diagonal()

def show_spy(input_matrix):
    plt.spy(input_matrix)
    plt.show()
#plt.spy(input_matrix)
#plt.show()
    
class Solver:
    def __init__(self, func):
        self.func = func
        
    def norm(self, x):
        return np.linalg.norm(x - x_exact) / np.linalg.norm(x)
    
    def output_results(self, name):
        print("Время решения с ",name, self.time)
        print("Норма решения ", self.norm(self.sol))
        
    def timeit(self, name, *args, **kwargs):
        start_time = time.clock()
        self.sol = self.func(*args, **kwargs)
        end_time = time.clock()
        self.time = end_time - start_time
        self.output_results(name)
        return end_time - start_time
    
    def timeit_with_info(self, name, *args, **kwargs):
        start_time = time.clock()
        self.sol, self.info = self.func(*args, **kwargs)
        end_time = time.clock()
        self.time = end_time - start_time
        self.output_results(name)
        return end_time - start_time

def S_precond(x):
    return spla.spsolve_triangular(L, x, lower = True)

def norm(x_0, x):
    return np.linalg.norm(x - x_0) / np.linalg.norm(x)

def Jacobi_precond(x):
    return x / diag



Jacobi_preconditioner = spla.LinearOperator((nx, nx), Jacobi_precond)
GaussSeidel_preconditioner = spla.LinearOperator(matrix.shape, GS_precond)
B = spla.splu(matrix)
Mz = lambda r: B.solve(r)
Cholesky_preconditioner = spla.LinearOperator(matrix.shape, Mz)   
#cond = np.linalg.cond(matrix)

start_time = time.clock()
x_exact = spla.spsolve(matrix, rhs)
time_0 = time.clock() - start_time
print("время решения с spsolve = ",time_0)

#bicg = Solver(spla.bicg)
#time_1 = bicg.timeit_with_info("bicg", matrix, rhs)
#x_1, info_1 = bicg.sol, bicg.info

#x_2, info_2, time_2 = timeit(spla.gmres, matrix, rhs)
#print("время решения gmres = ",time_2)

lgmres = Solver(spla.lgmres)
time_3 = lgmres.timeit_with_info("lgmres", matrix, rhs)
x_3, info_3 = lgmres.sol, lgmres.info

time_4 = lgmres.timeit_with_info("lgmres with jac", matrix, rhs, M=Jacobi_preconditioner)
x_4, info_4 = lgmres.sol, lgmres.info

time_5 = lgmres.timeit_with_info("lgmres with chol", matrix, rhs, M=Cholesky_preconditioner)
x_5, info_5 = lgmres.sol, lgmres.info

time_6 = lgmres.timeit_with_info("lgmres with gs ", matrix, rhs, M = GaussSeidel_preconditioner)
x_6, info_6 = lgmres.sol, lgmres.info

