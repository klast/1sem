# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 15:18:25 2017

@author: СпелеВова
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import random

eps = 1e-8

#генерирую симметрическую положительно определенную матрицу
def generate_matrix(nx, k):
    np.random.seed(153)
    B = np.random.rand(nx, nx) - 0.5
    B_trans = B.transpose()
    E = np.identity(nx)
    return B_trans @ B + k * E

class Solver:
    #конструктор
    def __init__(self, matrix):
        self.nx = matrix.shape[0]
        self.matrix = matrix
        self.init_sol = np.random.rand(self.nx)
        self.rhs = matrix @ self.init_sol
        self.eps = 1e-8
    
    #норма невязки    
    def residual_norm(self, sol):
        return np.linalg.norm(self.rhs - (self.matrix @ sol))
    
    def init_output(self):
        self.sol_norm_arr = [] 
        self.rel_sol_arr = []
        self.res_norm_arr = []
        self.rel_res_arr = []
    
    def construct_output(self, sol):
        self.sol_norm_arr.append(self.solution_norm(sol))
        self.rel_sol_arr.append(self.relative_solution(sol))
        self.res_norm_arr.append(self.residual_norm(sol))
        self.rel_res_arr.append(self.relative_residual(sol))
        
    def construct_dataframe(self):
        df = pd.DataFrame({
                '||x-x0||': self.sol_norm_arr,
                '||x-x0||/||x||': self.rel_sol_arr,
                '||b-Ax||': self.res_norm_arr,
                '||b-Ax||/||b||': self.rel_res_arr})
        return df
    
    #относительная невязка
    def relative_residual(self, sol):
        return self.residual_norm(sol) / np.linalg.norm(self.rhs)
    
    def solution_norm(self, sol):
        return np.linalg.norm(sol - self.init_sol)
    
    def relative_solution(self, sol):
        return self.solution_norm(sol) / np.linalg.norm(sol)
    
    #Петров-Галеркин
    def Petrov_Galerkin_steep(self, sol, r):
        return (self.rhs - self.matrix.dot(sol)).dot(r) / np.linalg.norm(r)
    
    def Petrov_Galerkin_minim(self, sol, r):
        return (self.rhs - self.matrix.dot(sol)).dot(self.matrix.dot(r)) / np.linalg.norm(r)
    
    #метод скорейшего спуска
    def steepest_descent(self):
        iter = 0
        sol = np.zeros(self.nx)
        r = self.rhs - (self.matrix @ sol)
        p = self.matrix @ r
        F_curr = 1e1
        self.init_output()
        residual_arr = np.array([self.relative_residual(sol)])
        while self.relative_residual(sol) > self.eps:
            iter += 1
            alpha = np.dot(r, r) / np.dot(p, r)
            sol += alpha * r
            if abs(self.Petrov_Galerkin_steep(sol, r)) > self.eps:
                print("Условие Петрова-Галеркина не выполнилось на шаге", iter)
            r -= alpha * p
            p = self.matrix @ r
            F_prev, F_curr = F_curr, np.dot(self.matrix @ sol, sol) - 2 * np.dot(self.rhs, sol)
            residual_arr = np.append(residual_arr, [self.relative_residual(sol)])
            self.construct_output(sol)
            if F_curr > F_prev:
                raise ValueError("F_curr > F_prev", F_prev)
        return sol, iter, self.construct_dataframe()
    
    #метод минимальной невязки
    def minimal_residual_iter(self):
        iter = 0
        sol = np.zeros(self.nx)
        r = self.rhs - (self.matrix @ sol)
        p = self.matrix @ r
        residual = self.residual_norm(sol)
        self.init_output()
        q_real = []
        residual_arr = np.array([self.relative_residual(sol)])
        while self.relative_residual(sol) > self.eps:
            iter += 1
            alpha = np.dot(p, r) / np.dot(p, p)
            sol += alpha * r
            if abs(self.Petrov_Galerkin_minim(sol, r)) > self.eps:
                print("Условие Петрова-Галеркина не выполнилось на шаге", iter)
            r_prev_norm = np.linalg.norm(r)
            r -= alpha * p
            q_real.append(np.linalg.norm(r) / r_prev_norm)
            p = self.matrix @ r
            residual_prev, residual = residual, self.residual_norm(sol)
            residual_arr = np.append(residual_arr, [self.relative_residual(sol)])
            self.construct_output(sol)
            if residual > residual_prev:
                raise ValueError("Норма невязки не уменьшилась на шаге", iter)
        return sol, iter, self.construct_dataframe(), q_real
    
    def conjugate_gradient(self):
        iter = 0
        sol = np.zeros(self.nx)
        r = self.rhs - (self.matrix @ sol)
        r_arr = [r]
        p = r
        p_arr = [p]
        residual = self.residual_norm(sol)
        self.init_output()
        residual_arr = np.array([self.relative_residual(sol)])
        while self.relative_residual(sol) > self.eps*1e-7:
            iter += 1
            alpha = np.dot(r, r) / np.dot(self.matrix @ p, r)
            sol += alpha * p
            r_prev , r = r, r - alpha * (self.matrix @ p)
            beta = np.dot(r, r) / np.dot(r_prev, r_prev)
            p = r + beta * p
            residual_prev, residual = residual, self.residual_norm(sol)
            self.construct_output(sol)
            residual_arr = np.append(residual_arr, [self.relative_residual(sol)])
            r_arr.append(r)
            p_arr.append(p)
            if residual > residual_prev:
                raise ValueError("Норма невязки не уменьшилась на шаге", iter)
        return sol, iter, self.construct_dataframe(), r_arr, p_arr


k = 50
nx = 500
#генериуем матрицу
mat = generate_matrix(nx, k)
print("Исходная матрица\n", mat)
#инициализируем решатель
my_solver = Solver(mat)
writer = pd.ExcelWriter('result.xlsx')
sigma = np.linalg.norm(mat)
mu = np.min(np.linalg.eigvals(mat.transpose() + mat), axis = 0)/2.
q_teor = math.sqrt(1 - mu*mu / (sigma*sigma))
#вызываем
x_1, iter_1, df_1 = my_solver.steepest_descent()
print("Метод скорейшего спуска", iter_1, "итераций")
print(df_1.tail(1))
df_1.to_excel(writer,'steepest_descent')
x_2, iter_2, df_2, q_real = my_solver.minimal_residual_iter()
print("Метод минимальной невязки", iter_2, "итераций")
print(df_2.tail(1))
df_2.to_excel(writer,'minimal_residual')
x_3, iter_3, df_3, r_arr, p_arr = my_solver.conjugate_gradient()
r_arr, p_arr = np.array(r_arr), np.array(p_arr)
print("Метод сопряженных градиентов", iter_3, "итераций")
print(df_3.tail(1))
print("Проверяем для случайных i,j")
n1 = random.randint(0, math.ceil(iter_3/2)-1)
n2 = random.randint(math.ceil(iter_3/2), iter_3-1)
print("(A*p_i,p_j) = ",(mat.dot(p_arr[n1, :])).dot(p_arr[n2, :]))
print("(A*p_i,p_i+1) = ", (mat.dot(p_arr[n1, :])).dot(p_arr[n1+1, :]))
print("(r_i,r_j) = ", r_arr[n1, :].dot(r_arr[n2, :]))
df_3.to_excel(writer,'conjugate_gradient')
writer.save()
os.system("start result.xlsx")
labels = ['Метод скорейшего спуска', 'Метод минимальной невязки', 'Метод сопряженных градиентов']
plt.semilogy(df_1['||b-Ax||/||b||'])
plt.semilogy(df_2['||b-Ax||/||b||'])
plt.semilogy(df_3['||b-Ax||/||b||'])
plt.legend(labels)
plt.show()
