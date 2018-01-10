
# coding: utf-8



import numpy as np
import math
from scipy.optimize import  minimize
from numpy import linalg
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import random




N = 100
a = 0.
b = math.pi
c = 5.
d = 6.
alpha=1e-8
sigma = 0
h=1e-1
#h=0.00012
#h= 2e-4
#ck=1
#ck = 0.02986315562
ck = 0.04761904761
#cn=0.25
cn = 0.02986315562






s = np.arange(a, b+(b-a)/N, (b-a)/N)
solution = []
for i in range(0,N+1):
    solution.append(2*math.cos(s[i]))
t = np.arange(c, d+(d-c)/N, (d-c)/N)
theta = np.random.uniform(-1, 1, t.shape[0])
epsilon = h

def K(s,t,z):
    return (z * math.cos(s) - 1) / (t+z)

def dKdz(s,t,z):
    #result = -math.cos(2 * s) / (t+z)**2
    result = t * math.cos(s) + 1 / (t+z)**2
    return result
best_norm = 1e5
saved_z = np.zeros(N+1)
def u(t):
    return (math.pi / (t*t-4)**0.5) * 0.25*((t*t-4)**0.5 - t)**2 

def dz(z):
    #global s
    dz=np.zeros(N+1)
    for i in range(1,N):
        _s = s[i] - (b-a)/(2*N)
        s_ = s[i] + (b-a)/(2*N)
        _z = z[i-1] + (z[i]-z[i-1])*(_s-s[i-1])/(s[i]-s[i-1])
        z_ = z[i+1] + (z[i]-z[i+1])*(s_-s[i+1])/(s[i]-s[i+1])
        dz[i] = (z_ - _z)
    _s = s[0] - (b-a)/(2*N)
    s_ = s[0] + (b-a)/(2*N)
    _z = z[1] + (z[0]-z[1])*(_s-s[1])/(s[0]-s[1])
    z_ = z[1] + (z[0]-z[1])*(s_-s[1])/(s[0]-s[1])
    dz[0] = (z_ - _z)
    _s = s[N] - (b-a)/(2*N)
    s_ = s[N] + (b-a)/(2*N)
    _z = z[N-1] + (z[N]-z[N-1])*(_s-s[N-1])/(s[N]-s[N-1])
    z_ = z[N-1] + (z[N]-z[N-1])*(s_-s[N-1])/(s[N]-s[N-1])
    dz[N] = (z_ - _z)
    return dz*N/(b - a)

def dz_2(z):
    global s
    s = np.array(s)
    z = np.array(z)
    dz = np.zeros(N+1)
    dz[0:-1] = np.diff(z) * np.diff(s)
    dz[-1] = (z[-1]-z[-2]) * (s[-1] - s[-2])
    return dz

def calc_z(z, beta, alpha):
    z_new =  z-beta*gradM(z, alpha)
    return z_new


def integral(f,a_,b_):
    return trapz(f, dx=(b_-a_)/N)

def Ur(z,t):
    a1 = integral(u(t),a,b) - (1/N)*random.uniform(0,1)*1e-4
    return a1

def omega(z, a, b):
    f = z**2 + dz(z)**2
    #f = z**2
    return integral(f,a,b)

#def helper(z):
    

def psi(z):
    #B = h*(h+2*ck)*(d-c)+(sigma*(h+2*ck)+h*(sigma+2*cn))*(d-c)/(b-a) + sigma*(sigma+2*cn)*(d-c)/((b-a)**2)
    B = (d - c)*(2*h*sigma + 2*h*cn + 2*sigma*ck+4*ck*cn*h*h + 2*h*ck + sigma*sigma + 2*sigma*cn)
    B_ = B*(b-a)*(b-a+1)
    return B_*(1 + omega(z,a,b))

def J(z):
    kk = np.zeros(N+1)
    intk = np.zeros(N+1)
    for i in range(0,N+1):
        for j in range(0,N+1):
            kk[j] = K(s[j], t[i], z[j])
        intk[i] = (integral(kk, a, b) - U[i])**2
    return integral(intk, c, d)

def M(z, alpha):
    return J(z) + alpha*omega(z,a,b)

def L_delta(z):
    return J(z) + psi(z)

def rho(z):
    global lambda_delta
    return J(z) - lambda_delta - psi(z)

def gradM(z, alpha):
    f = z - dz(dz(z))
    #f=z
    #f = z - dz(dz(z))
    kk = np.zeros((N+1,N+1))
    ii = np.zeros(N+1)
    dk = np.zeros(N+1)
    intk = np.zeros(N+1)
    for i in range(0,N+1):
        for j in range(0,N+1):
            kk[i][j] = K(s[j], t[i], z[j])
        ii[i] = integral(kk[i], a, b) - U[i]
    for k in range(0,N+1):
        for i in range(0,N+1):
            intk[i] = dKdz(s[k], t[i], z[k])*ii[i]
        dk[k] = integral(intk, c, d)
    return 2*dk + alpha*f


U = u(t) * (1 + epsilon * theta)
N_arr = [10, 50, 100]
#z_arr = []

for N in N_arr:
    s = np.arange(a, b+(b-a)/N, (b-a)/N)
    #solution = [2*math.cos(i) for i in s]
    if N == 10:
        solution = z_arr[0].copy()
    elif N==50:
        solution = z_arr[1].copy()
    else:
        solution = z_arr[2].copy()
    #for i in range(0,N+1):
        #solution.append(z_arr[0][i])
    t = np.arange(c, d+(d-c)/N, (d-c)/N)
    theta = np.random.uniform(-1, 1, t.shape[0])
    epsilon = h
    U = u(t) * (1 + epsilon * theta)
    i=0
    ii=0
    alpha_begin = 0
    alpha_end = 0.1
    z_init = np.zeros(N+1)
#z_init = np.array(solution.copy())
#lambda_delta = minimize(L_delta, np.ones(N+1), method='Powell', tol=1e-5).fun
    lambda_delta = 0
    print(lambda_delta)
    #z_init = np.array([-4 / math.pi * i + 2 for i in s])
#z_init = np.zeros(N+1)
#z_init = np.array([math.cos(i) for i in s])
    while True:
        alpha = (alpha_begin + alpha_end)* 0.5
        z = z_init.copy()
        i=0
        beta = 10
        while (linalg.norm(gradM(z,alpha))>1e-5):
            print(i, alpha, linalg.norm(gradM(z,alpha)), rho(z))
            #plt.plot(s,z)
            beta = 10
            i = i + 1
            z_ = calc_z(z, beta, alpha)
            while(M(z_, alpha)>M(z, alpha)):
                beta=beta/2
                z_ = calc_z(z, beta, alpha)
                if (beta<1e-10):
                    print("bad beta")
                    break
            if (np.abs(M(z_, alpha)-M(z, alpha))<1e-8):
                break
            z=calc_z(z, beta, alpha)
        if np.abs(alpha_begin-alpha_end) < 1e-5:
            z_arr.append(z)
            break
        print(i," iterations")
        print("rho =", rho(z))       
        temp = rho(z)
        if (np.abs(temp) < 1e-3):
            z_arr.append(z)
            break
        if temp > 0:
            alpha_end = alpha
        else:
            alpha_begin = alpha
    

    
#print (linalg.norm(gradM(z,alpha)))
#print("z = ", z)
#print("s = ", s)


fig = plt.figure()
l1, = plt.plot(s, z_init,  label='initial')
l2, = plt.plot(np.arange(a, b+(b-a)/10, (b-a)/10), z_arr[0],  label='10')
l3, = plt.plot(np.arange(a, b+(b-a)/50, (b-a)/50), z_arr[1],  label='50')
l4, = plt.plot(np.arange(a, b+(b-a)/100, (b-a)/100), z_arr[2],  label='100')
l5, = plt.plot(s, solution,  label='exact')
plt.legend(handles=[l1, l2, l3,l4,l5])
for jj in range(0, 10):
    t=0.2*jj
    #print("Nev = ", integral(z**t, a, b)-Ur(z,t), "t =", t)
plt.show()







