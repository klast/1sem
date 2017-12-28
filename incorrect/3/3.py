# -*- coding: utf-8 -*-
import numpy as np
from math import sin, cos, exp, sqrt
import numpy.linalg as ln
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import pandas as pd


target = np.array([1, 2]) # целевая точка
t_start = 0 # начальное время
t_end = 1 # конечное время
total_t = t_end - t_start # время в пути
nt = 100 # количество временных шагов
dt = (t_end - t_start) / (nt - 1) # длина временного шага

eps_d = 0.1 # погрешность по координате
eps_b = 0.1 # погрешность по управлению

delta = 0
h = 0

lambda_delta_arr = []
delta_arr = []
steps_arr = []
X0_arr = []
X1_arr = []


def arg_j_delta(t, x, u):
    a = np.dot(get_d_delta_mat(t), x)
    b = np.dot(get_b_delta_mat(t), u)
    return np.add(a, b)

def arg_psi(t, x, u):
    return np.dot(-np.transpose(get_d_delta_mat(t)), x)

def arg_base_sys(t, x, u):
    a = np.dot(get_d_delta_mat(t), x)
    b = np.dot(get_b_delta_mat(t), u[int(t/dt)])
    return np.add(a, b)

class ODESys:
    def __init__(self, u, arg):
        self._u = np.copy(u)
        self._arg = arg

    def get_value(self, t, x):
        return self._arg(t, x, self._u)

# целевая точка
def get_g_vec():
    return np.array([[1], [2]])
# матрица координат
def get_d_mat(t):
    return np.array([[cos(t), t], [1/(1+t), sin(t)]])
# матрица управления
def get_b_mat(t):
    return np.array([[1-exp(-t), 0], [0, 1+sin(2*t)]])
# возмущенная матрица координат
def get_d_delta_mat(t):
    return np.dot(get_d_mat(t), np.array([[1 + eps_d, 0], [0, 1 - eps_d]]))
# возмущенная матрица управления
def get_b_delta_mat(t):
    return np.dot(get_b_mat(t), np.array([[1 + eps_b * cos(t), -eps_b * sin(t)],
                                          [eps_b * sin(t), 1 + eps_b * cos(t)]]))

def get_d_b_mat_max():
    b_max = 0
    d_max = 0
    for i in range(0, nt):
        d_max = max(ln.norm(get_d_mat(i * dt)), d_max)
        b_max = max(ln.norm(get_b_mat(i * dt)), b_max)
    return d_max * 1.05, b_max * 1.05  # чтобы наверняка

# максимальное отклонение нормы в зависимости от погрешности
def get_delta():
    res = 0
    for i in range(0, nt):
        a = ln.norm(get_d_delta_mat(i * dt) - get_d_mat(i * dt))
        b = ln.norm(get_b_delta_mat(i * dt) - get_b_mat(i * dt))
        res = max(max(a, b), res)
    return res * 1.05  # чтобы наверняка

def get_j_delta(x):
    return ln.norm(x - target) ** 2

def get_psi(x, u):
    return (h*ln.norm(u) + delta)**2 + 2*(h*ln.norm(u) + delta)*(ln.norm(x)+ln.norm(target))

def get_lambda_delta_arg(u):
    s = ODESys(u, arg_j_delta)
    sol = integrate.solve_ivp(s.get_value, [t_start, t_end], np.array([0, 0]), t_eval=[t_end])
    x = sol.y.flatten()
    return get_j_delta(x) + get_psi(x, u)

def calc():
    # считаем lambda_delta один раз до цикла
    lambda_delta = opt.minimize(get_lambda_delta_arg, np.array([1, 1]), method='Powell', tol=1e-8).fun
    #print("LAMBDA_DELTA", lambda_delta)
    lambda_delta_arr.append(lambda_delta)

    u_vec = np.ones((nt, 2))

    step = 0
    while step < 1000:        
        # находим psi в цикле один раз за итерацию, u будем меняться (но не будет оптимизации по нему)
        s = ODESys(u_vec, arg_base_sys)
        sol = integrate.solve_ivp(s.get_value, [t_start, t_end], np.array([0, 0]), t_eval=[t_end])
        x = sol.y.flatten()
        psi_t_start = 2*(x - target)

        s = ODESys(None, arg_psi)

        sol = integrate.solve_ivp(s.get_value, [t_end, t_start], psi_t_start,
                                  t_eval=np.arange(t_end, t_start-1e-10, -dt))
        psi = np.flip(np.transpose(sol.y), 0)

        grad = np.zeros(u_vec.shape) # градиент
        b_arg_vec = np.zeros(nt) # вектор значений функций под интегралом в числителе выражения для step_length
        for i in range(0, nt):
            grad[i] = np.dot(np.transpose(get_b_delta_mat(i * dt)), psi[i])
            b_arg_vec[i] = ln.norm(grad[i])**2

        step_length = 0.5 * integrate.trapz(b_arg_vec, np.arange(t_start, t_end+1e-10, dt))

        # восстанавливаем координату по хитрому управлению
        s = ODESys(u_vec - grad, arg_base_sys)
        sol = integrate.solve_ivp(s.get_value, [t_start, t_end], np.array([0, 0]),
                                  t_eval=np.arange(t_start, t_end + 1e-10, dt))
        delta_x_vec = np.transpose(sol.y)

        # восстанавливаем координату по известному управлению решая ОДУ
        s = ODESys(u_vec, arg_base_sys)
        sol = integrate.solve_ivp(s.get_value, [t_start, t_end], np.array([0, 0]),
                                  t_eval=np.arange(t_start, t_end + 1e-10, dt))
        x_vec = np.transpose(sol.y)

        step_length /= ln.norm(delta_x_vec[-1] - x_vec[-1])**2

        new_u_vec = u_vec - step_length * grad
        tmp = ln.norm(new_u_vec - u_vec)
        if tmp < 1.e-5:
            break
        u_vec = new_u_vec[:]

        rho = get_j_delta(x_vec[-1]) - get_psi(x_vec[-1], u_vec[-1]) - lambda_delta
        if abs(rho) < 1.e-6:
            break        
        
        #print(rho)
        #if ln.norm(x_vec[-1] - target) < 1.e-6:
           #break
        step += 1

    #print("STEPS_NUM",step)
    steps_arr.append(step)
    return u_vec, x_vec

if __name__ == "__main__":
    eps_arr = [0, 0.001, 0.01, 0.1, 1]
    u = list()
    x = list()
    for i in range(0, len(eps_arr)):
        eps_d = eps_arr[i]
        eps_b = eps_arr[i]

        d_mat_max, b_mat_max = get_d_b_mat_max()
        delta = get_delta()
        c_0 = exp(d_mat_max * total_t) * b_mat_max * sqrt(total_t)
        h = delta * exp((d_mat_max * delta) * total_t) * (total_t * c_0 + sqrt(total_t))
        #print(delta)
        delta_arr.append(delta)

        _u, _x = calc()
        u.append(_u)
        x.append(_x)
    
 

    fig = plt.gcf()
    # fig.canvas.set_window_title(exp_name)
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.95, wspace=0.1, hspace=0.3)

    labels = ['0', '0.001', '0.01', '0.1', '1']

    for i in range(0, len(x)):
        tmp = x[i][-1]
        X0_arr.append(tmp[0])
        X1_arr.append(tmp[1])
    df1 = pd.DataFrame({'EPS'         : eps_arr,
                       'DELTA'        : delta_arr,
                       'LAMBDA_DELTA' : lambda_delta_arr,
                       'STEP'         : steps_arr,
                       'X_0'          : X0_arr,
                       'X_1'          : X1_arr})    
    print(df1)
    plt.subplot(2, 2, 1)
    plt.title(r'U_1')
    for res in u:
        plt.plot(res[:,0])
    plt.legend(labels)

    plt.subplot(2, 2, 3)
    plt.title(r'U_2')
    for res in u:
        plt.plot(res[:,1])
    plt.legend(labels)

    plt.subplot(2, 2, 2)
    plt.title(r'X_1')
    for res in x:
        plt.plot(res[:,0])
    plt.legend(labels)

    plt.subplot(2, 2, 4)
    plt.title(r'X_2')
    for res in x:
        plt.plot(res[:,1])
    plt.legend(labels)
    plt.show()