# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 15:18:25 2017

@author: СпелеВова
"""
import numpy as np

eps = 1e-8

#генерирую симметрическую положительно определенную матрицу
def generate_matrix(nx, k):
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
    
    #относительная невязка
    def relative_residual(self, sol):
        return self.residual_norm(sol) / np.linalg.norm(self.rhs)
    
    #Петров-Галеркин
    def Petrov_Galerkin(self, sol, r):
        return (self.rhs - self.matrix.dot(sol)).dot(r) / np.linalg.norm(r)
    
    #метод скорейшего спуска
    def steepest_descent(self):
        iter = 0
        sol = np.zeros(self.nx)
        r = self.rhs - (self.matrix @ sol)
        p = self.matrix @ r
        F_curr = 1e6
        residual_arr = np.array([self.relative_residual(sol)])
        while self.relative_residual(sol) > self.eps:
            iter += 1
            alpha = np.dot(r, r) / np.dot(p, r)
            sol += alpha * r
            #if abs(self.Petrov_Galerkin(sol, r)) > self.eps:
                #print("Условие Петрова-Галеркина не выполнилось на шаге", iter)
            r -= alpha * p
            p = self.matrix @ r
            F_prev, F_curr = F_curr, np.dot(self.matrix @ sol, sol) - 2 * np.dot(self.rhs, sol)
            residual_arr = np.append(residual_arr, [self.relative_residual(sol)])
            if F_curr > F_prev:
                raise ValueError("F_curr > F_prev", F_prev)
        return sol, iter, residual_arr
    
    #метод минимальной невязки
    def minimal_residual_iter(self):
        iter = 0
        sol = np.zeros(self.nx)
        r = self.rhs - (self.matrix @ sol)
        p = self.matrix @ r
        residual = self.residual_norm(sol)
        while self.relative_residual(sol) > self.eps:
            iter += 1
            alpha = np.dot(p, r) / np.dot(p, p)
            sol += alpha * r
            if abs(self.Petrov_Galerkin(sol, r)) > self.eps:
                print("Условие Петрова-Галеркина не выполнилось на шаге", iter)
            r -= alpha * p
            p = self.matrix @ r
            residual_prev, residual = residual, self.residual_norm(sol)
            if residual > residual_prev:
                raise ValueError("Норма невязки не уменьшилась на шаге", iter)
        return sol, iter
    
    def conjugate_gradient(self):
        iter = 0
        sol = np.zeros(self.nx)
        r = self.rhs - (self.matrix @ sol)
        p = r
        residual = self.residual_norm(sol)
        while self.relative_residual(sol) > self.eps:
            iter += 1
            alpha = np.dot(r, r) / np.dot(self.matrix @ p, r)
            sol += alpha * p
            r_prev , r = r, r - alpha * (self.matrix @ p)
            beta = np.dot(r, r) / np.dot(r_prev, r_prev)
            p = r + beta * p
            residual_prev, residual = residual, self.residual_norm(sol)
            if residual > residual_prev:
                raise ValueError("Норма невязки не уменьшилась на шаге", iter)
        return sol, iter


k = 0.6
nx = 10
#генериуем матрицу
mat = generate_matrix(nx, k)
print("Исходная матрица\n", mat)
#инициализируем решатель
my_solver = Solver(mat)
#вызываем
x_1, iter_1, resid_arr_1 = my_solver.steepest_descent()
x_2, iter_2 = my_solver.minimal_residual_iter()
x_3, iter_3 = my_solver.conjugate_gradient()
