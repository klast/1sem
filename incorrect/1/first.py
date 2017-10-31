# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""
import numpy as np
import pandas as pd
import scipy.optimize as sc
import matplotlib.pyplot as plt


def rho(A, Za,u,ld):
    return np.linalg.norm(np.dot(A,Za) - u)**2 - (ld + sigma + h * np.linalg.norm(Za))**2

def my_solve(A_approx, U_approx, ld, h, sigma):
    B = np.dot(A_approx.T, u_approx)
    A_k = np.dot(A_approx.T, A_approx)
    alpha_0 = 1e-10
    alpha_n = 1000
    tmp = 1
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
    return np.linalg.norm(np.dot(A_approx, z) - u_approx) + h*np.linalg.norm(z) + sigma

def myplot(values, indexes, title, name):
    for i in indexes:
        plt.plot(x, values[i])
    #plt.plot(x, z_all[0],x , z_all[1], x ,z_all[2], x , z_all[3]) #Построение графика
    plt.title(title)
    plt.grid(True) #Сетка
    legends=['1e-4','1e-3','1e-2','1e-1']
    plt.legend(legends)
    #plt.show() #Показать график
    plt.savefig(name)
    plt.clf()
    

n = 10
m = 10
bounds = []
for i in range(10):
    bounds.append((0, 10))

if __name__ == "__main__":
#1
    A = np.random.random((n, m))
    if np.linalg.det(A) < 1000:
        A = A + 3 * np.eye(m)
#2
    z = np.random.random(n)
#3
    u = np.dot(A, z)
#4
    z_vec = [z for i in range(10)]
    theta_A = np.random.uniform(-1, 1, (n, m))
    theta_u = np.random.uniform(-1, 1, m)
    min_value = []
    real_value = []
    h_vec = []
    sigma_vec = []
    norm_1 = []
    norm_2 = []
    epsilon = []
    alpha = []
    Zalpha = []
    for i in reversed(range(4)):
        for j in reversed(range(4)):
            if i==j or i==3 or j==3:
                epsilon.append([10**(-i-1), 10**(-j-1)])
#1 FIRST TASK
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
        lambda_delta = sc.differential_evolution(LambdaDelta, bounds, tol = 0.0001)
        min_value.append(LambdaDelta(lambda_delta.x))
        #print(LambdaDelta(lambda_delta.x))
        real_value.append(LambdaDelta(z))
        #print(LambdaDelta(z))
#LD_val = sc.minimize(LD, x0, method = 'nelder-mead', options={'xtol': 1e-8, 'disp': True}, args=(h, sigma, A_approx, u_approx))
        #print(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        norm_1.append(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), h, sigma)
        alpha.append(ans[0])
        Zalpha.append(ans[1])
        res = np.linalg.norm(z-ans[1])
        #print(res)
        norm_2.append(np.linalg.norm(z-ans[1]))
    df1 = pd.DataFrame({
            'epsilon_A' : [epsilon[i][0] for i in range(10)],
            'epsilon_u' : [epsilon[i][1] for i in range(10)],
            'LD_min' : min_value,
            'LD_real' : real_value,
            'h' : h_vec,
            'sigma' : sigma_vec,
            'alpha' : alpha,
            'z' : z_vec,
            'z_alpha' : Zalpha,
            'norm_1' : norm_1,
            'norm_2' : norm_2
            })
    
    x = np.arange(0,10,1) #Массив значений аргумента
    z_all = [Zalpha[i] for i in range(10)]
    myplot(z_all,[0,1,2,3],'EPS_A=10-4','1_eps_A.png')
    myplot(z_all,[0,4,6,8],'EPS_u=10-4','1_eps_u.png')
    myplot(z_all,[0,5,7,9],'EPS_A=EPS_u','1_eps_Au.png')
    
#1 SECOND TASK
    m = 12
    A0 = A.copy()
    u0 = u.copy()
    A = np.zeros((12,10))
    u = np.zeros(12)
    for i in range(0,10):
        A[i] = A0[i]
        u[i] = u0[i]
    
    A[10] = A[0] + A[9]
    A[11] = A[1] + A[8]
    u[10] = u[0] + u[9]
    u[11] = u[1] + u[8]
    theta_A = np.random.uniform(-1, 1, (12, 10))
    theta_u = np.random.uniform(-1, 1, 12)
    min_value = []
    real_value = []
    h_vec = []
    sigma_vec = []
    norm_1 = []
    norm_2 = []
    alpha = []
    Zalpha = []
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
        lambda_delta = sc.differential_evolution(LambdaDelta, bounds, tol = 0.0001)
        min_value.append(LambdaDelta(lambda_delta.x))
        #print(LambdaDelta(lambda_delta.x))
        real_value.append(LambdaDelta(z))
        #print(LambdaDelta(z))
#LD_val = sc.minimize(LD, x0, method = 'nelder-mead', options={'xtol': 1e-8, 'disp': True}, args=(h, sigma, A_approx, u_approx))
        #print(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        norm_1.append(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), h, sigma)
        alpha.append(ans[0])
        Zalpha.append(ans[1])
        res = np.linalg.norm(z-ans[1])
        #print(res)
        norm_2.append(np.linalg.norm(z-ans[1]))
    df2 = pd.DataFrame({
            'epsilon_A' : [epsilon[i][0] for i in range(10)],
            'epsilon_u' : [epsilon[i][1] for i in range(10)],
            'LD_min' : min_value,
            'LD_real' : real_value,
            'h' : h_vec,
            'sigma' : sigma_vec,
            'alpha' : alpha,
            'z' : z_vec,
            'z_alpha' : Zalpha,
            'norm_1' : norm_1,
            'norm_2' : norm_2
            })
    myplot(z_all,[0,1,2,3],'EPS_A=10-4','2_eps_A.png')
    myplot(z_all,[0,4,6,8],'EPS_u=10-4','2_eps_u.png')
    myplot(z_all,[0,5,7,9],'EPS_A=EPS_u','2_eps_Au.png')

    
    #1 THIRD TASK
    m = 8
    A = np.zeros((8,10))
    u = np.zeros(8)
    for i in range(0,8):
        A[i] = A0[i]
        u[i] = u0[i]
    theta_A = np.random.uniform(-1, 1, (8, 10))
    theta_u = np.random.uniform(-1, 1, 8)
    min_value = []
    real_value = []
    h_vec = []
    sigma_vec = []
    norm_1 = []
    norm_2 = []
    alpha = []
    Zalpha = []
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
        lambda_delta = sc.differential_evolution(LambdaDelta, bounds, tol = 0.0001)
        min_value.append(LambdaDelta(lambda_delta.x))
        #print(LambdaDelta(lambda_delta.x))
        real_value.append(LambdaDelta(z))
        #print(LambdaDelta(z))
#LD_val = sc.minimize(LD, x0, method = 'nelder-mead', options={'xtol': 1e-8, 'disp': True}, args=(h, sigma, A_approx, u_approx))
        #print(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        norm_1.append(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), h, sigma)
        alpha.append(ans[0])
        Zalpha.append(ans[1])
        res = np.linalg.norm(z-ans[1])
        #print(res)
        norm_2.append(np.linalg.norm(z-ans[1]))
    df3 = pd.DataFrame({
            'epsilon_A' : [epsilon[i][0] for i in range(10)],
            'epsilon_u' : [epsilon[i][1] for i in range(10)],
            'LD_min' : min_value,
            'LD_real' : real_value,
            'h' : h_vec,
            'sigma' : sigma_vec,
            'alpha' : alpha,
            'z' : z_vec,
            'z_alpha' : Zalpha,
            'norm_1' : norm_1,
            'norm_2' : norm_2
            })
    
    myplot(z_all,[0,1,2,3],'EPS_A=10-4','3_eps_A.png')
    myplot(z_all,[0,4,6,8],'EPS_u=10-4','3_eps_u.png')
    myplot(z_all,[0,5,7,9],'EPS_A=EPS_u','3_eps_Au.png')
    
#1 FOURTH TASK
    m = 10
    A = np.zeros((10,10))
    A = A0.copy()
    u = u0.copy()
    A[9] = A[8]
    theta_A = np.random.uniform(-1, 1, (10, 10))
    theta_u = np.random.uniform(-1, 1, 10)
    min_value = []
    real_value = []
    h_vec = []
    sigma_vec = []
    norm_1 = []
    norm_2 = []
    alpha = []
    Zalpha = []
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
        lambda_delta = sc.differential_evolution(LambdaDelta, bounds, tol = 0.0001)
        min_value.append(LambdaDelta(lambda_delta.x))
        #print(LambdaDelta(lambda_delta.x))
        real_value.append(LambdaDelta(z))
        #print(LambdaDelta(z))
#LD_val = sc.minimize(LD, x0, method = 'nelder-mead', options={'xtol': 1e-8, 'disp': True}, args=(h, sigma, A_approx, u_approx))
        #print(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        norm_1.append(np.linalg.norm(np.dot(A_approx, lambda_delta.x) - u_approx))
        ans = my_solve(A_approx, u_approx, LambdaDelta(lambda_delta.x), h, sigma)
        alpha.append(ans[0])
        Zalpha.append(ans[1])
        res = np.linalg.norm(z-ans[1])
        #print(res)
        norm_2.append(np.linalg.norm(z-ans[1]))
    df4 = pd.DataFrame({
            'epsilon_A' : [epsilon[i][0] for i in range(10)],
            'epsilon_u' : [epsilon[i][1] for i in range(10)],
            'LD_min' : min_value,
            'LD_real' : real_value,
            'h' : h_vec,
            'sigma' : sigma_vec,
            'alpha' : alpha,
            'z' : z_vec,
            'z_alpha' : Zalpha,
            'norm_1' : norm_1,
            'norm_2' : norm_2
            })
    
    myplot(z_all,[0,1,2,3],'EPS_A=10-4','4_eps_A.png')
    myplot(z_all,[0,4,6,8],'EPS_u=10-4','4_eps_u.png')
    myplot(z_all,[0,5,7,9],'EPS_A=EPS_u','4_eps_Au.png')
    
    writer = pd.ExcelWriter('output.xlsx')
    df1.to_excel(writer,'sheet1')
    df2.to_excel(writer,'sheet2')
    df3.to_excel(writer,'sheet3')
    df4.to_excel(writer,'sheet4')
    dfz = pd.DataFrame({
            'z' : z_vec,
            })
    dfu = pd.DataFrame({
            'u' : u0
            })
    dfz.to_excel(writer, 'z')
    dfu.to_excel(writer, 'u')
    writer.save()
    writer.close()    


    
