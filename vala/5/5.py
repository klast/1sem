# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 00:20:24 2017

@author: СпелеВова
"""
import inspect
import scipy as sc
import numpy as np
import pandas as pd
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import time

#входная матрица
input_matrix = sc.io.mmread('0.4solve.mtx')
#правая часть
rhs = sc.io.mmread('0.rhs')
#количество элементов
nx = rhs.shape[0]
#перевод в csc
matrix = input_matrix.tocsc()
L = sc.sparse.tril(matrix).tocsr()
#главная диагональ
diag = matrix.diagonal()

#выводит картинку структуры разряженной матрицы
#def show_spy(input_matrix):
#    plt.spy(input_matrix, markersize=1)
#    plt.show()
#show_spy(matrix)
#raise ValueError(1)
#show_spy(L)
#raise ValueError(1)   

#для предобуславливателя Гаусса-Зейделя
def GS_precond(x):
    return spla.spsolve_triangular(L, x, lower = True)

#для вычисления нормы по решению
def norm(x_0, x):
    return np.linalg.norm(x - x_0) / np.linalg.norm(x)

#для предобуславливателя Якоби
def Jacobi_precond(x):
    return x / diag

#замер времени построения пред-ля Якоби
begin_time = time.clock()
Jacobi_preconditioner = spla.LinearOperator((nx, nx), Jacobi_precond)
end_time = time.clock()
time_Jacobi_precond = end_time - begin_time
print("Время построения предобуславливателя Якоби = ", end_time - begin_time)

#замер времени построения пред-ля Гаусса-Зейделя
begin_time = time.clock()
GaussSeidel_preconditioner = spla.LinearOperator((nx, nx), GS_precond)
end_time = time.clock()
time_GaussSeidel_precond = end_time - begin_time
print("Время построения предобуславливателя Гаусса-Зейделя = ", end_time - begin_time)

#замер времени построения пред-ля на неполном ILU разложении
begin_time = time.clock()
B = spla.spilu(matrix, drop_tol=1e-2)
Mz = lambda r: B.solve(r)
Spilu_preconditioner = spla.LinearOperator(matrix.shape, Mz)   
end_time = time.clock()
time_Spilu_precond = end_time - begin_time
print("Время построения преобуславливателя на неполном ILU разложении", end_time - begin_time)
print()

def make_precond_obj(func, tol, maxiter, time):
    return {'func' : func, 'tol' : tol, 'maxiter' : maxiter, 'time' : time}

#словарь предобуславливателей, по имени выдать [функция, tol, maxiter, время построения]
preconds = {'None'        : make_precond_obj(None, 1e-8, 1000, None),
            'Spilu'       : make_precond_obj(Spilu_preconditioner, 1e-8, 1000, time_Spilu_precond),
            'Jacobi'      : make_precond_obj(Jacobi_preconditioner, 1e-8, 1000, time_Jacobi_precond),
            'GaussSeidel' : make_precond_obj(GaussSeidel_preconditioner, 1e-8, 100, time_GaussSeidel_precond)
            }

#Находим точное решение с помощью spsolve
start_time = time.clock()
x_exact = spla.spsolve(matrix, rhs, use_umfpack = True)
time_0 = time.clock() - start_time
print("Время решения с spsolve = ",time_0)

#для отчета
report_columns = ['Число итераций', 'Время решения', 'Полное время', 'Невязка', 'Время 1 итерации', 'Отклонение решения'] 

#вспомогательный класс для автоматизации
#matrix - глобальная переменная - матрица разреяженная
#rhs - глобальная переменная - правая часть
#x_exact - глобальная переменная - решение с spsolve
class Solver:
    def __init__(self, name, func):
        self.name = name
        self.func = func
        self.res = []
        self.df = pd.DataFrame(columns=report_columns)
        
    def calc_norm(self, x):
        return np.linalg.norm(x - x_exact) / np.linalg.norm(x)
    
    def calc_resid(self, x):
        return np.linalg.norm(rhs - matrix.dot(x.reshape(rhs.shape[0],1))) / np.linalg.norm(rhs)
    
    def add_to_report(self, pr_name, *args, **kwargs):
        resid= self.calc_resid(self.sol)
        time_1_iter = self.time / self.iter
        norm = self.calc_norm(self.sol)
        if pr_name == 'None':
            full_time = self.time            
        else:
            full_time = self.time + preconds[pr_name]['time']            
        self.df.loc[pr_name] = [self.iter, self.time, full_time, resid, time_1_iter, norm]
    def callback_iter(self, x):
        global iter
        iter+=1
    
    def callback_resid(self, x):
        global iter
        if self.func == spla.gmres:
            #print(x)
            self.res.append(x)
        else:
            my_resid = np.linalg.norm(rhs - matrix.dot(x.reshape(rhs.shape[0],1))) / np.linalg.norm(rhs)
            #print(my_resid)
            self.res.append(my_resid)
        iter+=1
        
    def timeit(self, *args, **kwargs):
        start_time = time.clock()
        self.sol = self.func(*args, **kwargs)
        end_time = time.clock()
        self.time = end_time - start_time
        #self.output_results(name)
        return end_time - start_time
    
    def timeit_with_info(self, pr_name, *args, **kwargs):
        global iter
        start_time = time.clock()
        self.sol, self.info = self.func(*args, callback = self.callback_iter, **kwargs)
        end_time = time.clock()
        self.time = end_time - start_time
        self.iter = iter
        self.add_to_report(pr_name, *args, **kwargs)        
        iter = 0
        return end_time - start_time
    
    def get_residual(self, *args, **kwargs):
        global iter
        self.res = []
        self.sol, self.info = self.func(*args, callback = self.callback_resid, **kwargs)
        self.res = self.res[:iter]
        iter = 0
        return self.res

iter = 0

solvers_func = {'bicgstab' : spla.bicgstab,
                'gmres'    : spla.gmres,
                'lgmres'   : spla.lgmres,
                'cg'       : spla.cg,
                'cgs'      : spla.cgs,
                'minres'   : spla.minres}

solvers_obj = {name:Solver(name, solvers_func[name]) for name in solvers_func}

writer = pd.ExcelWriter('result.xlsx')
#основной расчет для формирования таблицы
for name in solvers_obj:
    for pr_name in preconds:
        _tol = preconds[pr_name]['tol']
        _maxiter = preconds[pr_name]['maxiter']
        _func = preconds[pr_name]['func']
        solvers_obj[name].timeit_with_info(pr_name, matrix, rhs, tol=_tol, maxiter=_maxiter, M=_func)
        print("Посчитан " + name + " с " + pr_name + " предобусл ")
    print(solvers_obj[name].df[['Число итераций', 'Время решения', 'Невязка']])
    #solvers_obj[name].df.to_excel(writer, name)
    
#сохранение в Excel    
keys=[name for name in solvers_func]
frames = [solvers_obj[name].df for name in keys]
result = pd.concat(frames, keys=keys)
result.to_excel(writer, 'result')
writer.save()
writer.close()

plt.spy(matrix, markersize=1)
plt.savefig('spy')
plt.clf()

#def save_residuals(pr_name):
    #plt.figure(figsize=(12,8))
    #if pr_name == None:
        #t

for pr_name in preconds:
    plt.figure(figsize=(12,8))
    _tol = preconds[pr_name]['tol']
    _maxiter = preconds[pr_name]['maxiter']
    _func = preconds[pr_name]['func']
    labels = []
    for name in solvers_obj:
        resid = solvers_obj[name].get_residual(matrix, rhs, tol=_tol, maxiter=_maxiter, M=_func)
        print("Найдены невязки для " + name + "_" + pr_name)
        labels.append(name+"_"+pr_name)
        plt.semilogy(resid)
    plt.legend(labels)
    plt.savefig(pr_name)
    plt.clf()
    
    

#print("bicgstab chol resid")
#resid_bicgstab_chol = bicgstab.get_residual("bicgstab_chol", matrix, rhs,tol=1e-8,maxiter=1000,M=Cholesky_preconditioner)
#print("lgmres chol resid")
#resid_lgmres_chol = lgmres.get_residual("lgmres_chol", matrix, rhs,tol=1e-8,maxiter=1000, M=Cholesky_preconditioner)
#print("gmres chol resid")
#resid_gmres_chol = gmres.get_residual("gmres_chol", matrix, rhs,tol=1e-8,maxiter=1000, M=Cholesky_preconditioner)

#labels = ['bicgstab_spilu', 'lgmres_spilu', 'gmres_spilu']
#plt.semilogy(resid_bicgstab_chol)
#plt.semilogy(resid_lgmres_chol)
#plt.semilogy(resid_gmres_chol)
#plt.legend(labels)
#plt.show()



