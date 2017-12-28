# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 17:12:42 2017

@author: СпелеВова
"""

import numpy as np
import scipy.linalg as sc
from IPython import get_ipython

np.set_printoptions(precision = 3)

# проверка матрицы на симметричность, ортогональность и собственные числа
def check_matrix(reflection_matrix):
    print("Матрица ортогональна: \n", reflection_matrix @ reflection_matrix.transpose())
    if np.array_equal(reflection_matrix, reflection_matrix.transpose()):
        print("Матрица симметрична")
    print("Собственные числа: \n", np.linalg.eig(reflection_matrix))

# оптимизированное умножение    
def optimized_mult(mat, x):
    n = x.shape[0]
    return mat - 2 * x.reshape((n,1))*(x.reshape((1,n)).dot(mat))

# построение матрицы отражений
def get_reflection_matrix(matrix_size):
    x = np.random.rand(matrix_size)
    x = (x / np.linalg.norm(x)).reshape((matrix_size, 1))
    trans_x = x.transpose()
    reflection_matrix = np.identity(matrix_size) - 2 * x@trans_x
    return x, reflection_matrix

def reflection_method(mat, x):
    
    #заполняю вспомогательные переменные
    n = mat.shape[0]
    b = mat@x.reshape(n, 1)
    mat_0 = mat.copy()
    R = mat.copy()
    for i in range(n-1):
        mat_0[:,i] = R[:,i].copy()
        mat_0[:i,i] = 0
        norm_mat_0 = np.linalg.norm(mat_0[:,i])
        y = mat_0[:,i].reshape((n,1))
        y[i,0] -= norm_mat_0
        norm_y = np.linalg.norm(y)
        x_0 = y / norm_y
        R = optimized_mult(R, x_0)        
        b = b - 2*x_0.reshape((n,1))*x_0.reshape((1,n)).dot(b)        
    
    return b, R

def reflection_method_QR(mat, x):
    
    def optimized(mat, x):
        n = x.shape[0]
        return mat - 2 * mat@(x.reshape((n,1))@x.reshape((1,n))).transpose()
    
    #заполняю вспомогательные переменные
    n = mat.shape[0]
    b = mat@x.reshape(n, 1)
    mat_0 = mat.copy()
    Q = np.identity(n)
    R = mat.copy()
    for i in range(n-1):
        mat_0[:,i] = R[:,i].copy()
        mat_0[:i,i] = 0
        norm_mat_0 = np.linalg.norm(mat_0[:,i])
        y = mat_0[:,i].reshape((n,1))
        y[i,0] -= norm_mat_0
        norm_y = np.linalg.norm(y)
        x_0 = y / norm_y
        R = optimized_mult(R, x_0)
        Q = optimized(Q, x_0)
        b = b - 2*x_0.reshape((n,1))*x_0.reshape((1,n)).dot(b)        
    
    return b, Q, R

def qr_solve(b, test_matrix):
    nx = x.shape[0]
    sol = np.zeros(nx).reshape((nx,1))
    for i in range(nx-1, -1, -1):
        s = test_matrix[i, i+1:nx]@sol[i+1:nx]
        sol[i] = (b[i] - s) / test_matrix[i,i]
    return sol 
        
    
    
# строим матрицу отражений 
x, reflection_matrix = get_reflection_matrix(3)

# проверяем ее
check_matrix(reflection_matrix)

#проверяем работу больших матрицах
nx = 500

test_matrix = np.random.rand(nx, nx)
x, reflection_matrix = get_reflection_matrix(nx)

ipython = get_ipython()

print("Время работы обычного домножения")
ipython.magic("timeit reflection_matrix.dot(test_matrix)")

print("Время работы оптимизированной процедуры")
ipython.magic("timeit optimized_mult(test_matrix, x)")

#метод отражений
nx = 2000
test_matrix = np.random.rand(nx, nx)
x = np.random.rand(nx)
x = (x / np.linalg.norm(x)).reshape((nx, 1))

print("Применяем метод отражений")
b, R = reflection_method(test_matrix, x)
print("Начальная матрица\n",test_matrix)
print("Новая матрица\n",R)

#решаем систему с использованием полученного ранее b
print("Время работы нашего решателя")
qr_x = qr_solve(b, R)
ipython.magic("timeit qr_solve(b,R)")
print("Время работы стандартного решателя")
triang_x = sc.solve_triangular(R, b)
ipython.magic("timeit sc.solve_triangular(R, b)")
print("Решение с помощью b различается на", np.linalg.norm(qr_x - triang_x, 2))

nx = 50
test_matrix = np.random.rand(nx, nx)
x = np.random.rand(nx)
x = (x / np.linalg.norm(x)).reshape((nx, 1))
b, Q, R = reflection_method_QR(test_matrix, x)

#проверяем, что QR = A
print("A - QR = \n",Q@R-test_matrix)
#проверяем, что QR == Q и R из np.linalg.qr
D = np.diag(np.sign(np.diag(R)))
Q=Q@D
R=D@R

M, N = np.linalg.qr(test_matrix)
D=np.diag(np.sign(np.diag(N)))
M=M@D
N=D@N
print("norm(M - Q) = ",np.linalg.norm(M-Q,1))
print("norm(N - R) = ",np.linalg.norm(N-R,1))
