
# coding: utf-8



import numpy as np
import math
from scipy.optimize import  minimize
from numpy import linalg
from scipy.integrate import trapz
from scipy.special import gamma
import matplotlib.pyplot as plt
import random




N = 25
a = 0.
b = math.pi
c = 5.
d = 6.
alpha=1e-8
sigma = 1e-3
#h=1e-2
h=1e-4
#ck=1
#ck = 0.02986315562
ck = 0.04761904761
#cn=0.25
#cn = -0.035714285714
cn = 0.02986315562





s = np.arange(a, b+(b-a)/N, (b-a)/N)
solution = []
for i in range(0,N+1):
    solution.append(2*math.cos(s[i]))
t = np.arange(c, d+(d-c)/N, (d-c)/N)
theta = np.random.uniform(-1, 1, t.shape[0])
epsilon = 1e-1



def K(s,t,z):
    #return (1. - s) / ((t + z) ** 3)
    #return (1.)/((1+t*z)**2)
    #return z**t
    return math.cos(2*s) / (t+z)

def dKdz(s,t,z):
    #return 3 * (s - 1.) / ((t + z) ** 4)
    #return (-2.*t)/((1+t*z)**3)
    #return 3*(s-1.)/((t+z)**4)
    #if (z != 0 ):
        #result = t*(z**(t-1))
    result = -math.cos(2 * s) / (t+z)**2
    #if (z == 0):
        #result = 0
   # print ("t = ", t, "Z=", z, "res = ", result)
    return result

def u(t):
    #return (math.pi)/((1-t**2)**(1.5))
    #return 1./(2*t*t*(1+t))
    #return ((math.pi)/2)*(gamma((t+1)/2)/gamma(t/2 + 1))
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
    #s = np.array(s)
    #z = np.array(z)
    #dz = np.zeros(N+1)
    #dz[0:-1] = np.diff(z) / np.diff(s)
    #dz[-1] = (z[-1]-z[-2]) / (s[-1] - s[-2])
    #return dz*N/(b-a)
def dz_2(z):
    global s
    s = np.array(s)
    z = np.array(z)
    dz = np.zeros(N+1)
    dz[0:-1] = np.diff(z) * np.diff(s)
    dz[-1] = (z[-1]-z[-2]) * (s[-1] - s[-2])
    return dz

def integral(f,a_,b_):
    #s = 0
    #for i in range(0,N):
    #    s = s + f[i+1] + f[i]
    #return s*(b_ - a_)/(2*N)
    return trapz(f, dx=(b_-a_)/N)

def Ur(z,t):
    a1 = integral(u(t),a,b) - (1/N)*random.uniform(0,1)*1e-4
    return a1

def omega(z, a, b):
    #f = z**2 + dz(z)**2
    f = z**2
    return integral(f,a,b)

def psi(z):
    B = h*(h+2*ck)*(d-c)+(sigma*(h+2*ck)+h*(sigma+2*cn))*(d-c)/(b-a) + sigma*(sigma+2*cn)*(d-c)/((b-a)**2)
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

def rho(z):
    return J(z) - psi(z)

def gradM(z, alpha):
    #f = z - dz(dz(z))
    #f=z
    f = z - dz(dz(z))
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

def my_solve(A_approx, U_approx, ld, h, sigma):
    B = np.dot(A_approx.T, u_approx)
    A_k = np.dot(A_approx.T, A_approx)
    alpha_0 = 1e-10
    alpha_n = 1000
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


U = u(t) * (1 + epsilon * theta)
z = np.ones(N+1)
beta=1.
i=0
ii=0
alpha_begin = 0
alpha_end = 10000
alpha = 1000
#z_init = np.ones(N+1)
z_init = np.array([-4 / math.pi * i + 2 for i in s])
while True:
    #alpha = (alpha_begin + alpha_end)* 0.5
    #print("alpha = ", alpha)
    z = z_init.copy()
    i=0
    beta = 0.1
    while (linalg.norm(gradM(z,alpha))>1e-5):
        print(i, alpha, linalg.norm(gradM(z,alpha)), rho(z))
        beta = 0.00001
        i = i + 1
        z_ = z-beta*gradM(z, alpha)
        z_[0] = 2
        z_[N] = -2
        #print(z_)
        #print("z_ = ", z_)
        while(M(z_, alpha)>M(z, alpha)):
            beta=beta/2
            z_ = z-beta*gradM(z, alpha)
            z_[0] = 2
            z_[N] = -2
            if (beta<1e-10):
                break
        if (np.abs(M(z_, alpha)-M(z, alpha))<1e-5):
            break
        z=z-beta*gradM(z,alpha)
        z[0]=2
        z[N]=-2
    print(z)
    print(i," iterations")
    print("rho =", rho(z))
    break
    #temp = rho(z)
    #if (np.abs(temp) < 1e-3):
    #    break
    #if temp > 0:
    #    alpha_end = alpha
    #else:
        #alpha_begin = alpha
    

    
#print (linalg.norm(gradM(z,alpha)))
#print("z = ", z)
#print("s = ", s)

abs = np.arange(a, b+(b-a)/N, (b-a)/N)
fig = plt.figure()
l1, = plt.plot(abs, [-4 / math.pi * i + 2 for i in s],  label='initial')
l2, = plt.plot(abs ,z,  label='z')
l3, = plt.plot(abs ,solution,  label='exact')
plt.legend(handles=[l1, l2, l3])
for jj in range(0, 10):
    t=0.2*jj
    #print("Nev = ", integral(z**t, a, b)-Ur(z,t), "t =", t)
plt.show()







