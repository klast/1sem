# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 15:14:12 2018

@author: spelevova
"""
import numpy as np
from numba import jit
import math
import scipy.optimize as sc
from numpy import linalg
from scipy.integrate import trapz
from scipy.special import gamma
import matplotlib.pyplot as plt
import random

epsilon = 1e-4
N = 10
a = 0.
b = math.pi
c = 5.
d = 6.
numStepT = N
numStepS = N
deltaT = (d-c)/numStepT
deltaS = (b-a)/numStepS
t = np.arange(c, d+deltaT, deltaT)
s = np.arange(a, b+deltaS, deltaS)
theta = np.random.uniform(-1, 1, t.shape[0])
alpha=1e-8
beta = 1.
h = 0.0001
sigma = 0
ck = 0.04761904761
cu = 0.02986315562
#z_k_init = np.array([-4 / math.pi * i + 2 for i in s])
z_k_init = np.ones(N+1)
#z_k_init = np.array([1 for i in s])

def intS1D(Z1D):
    return np.sum(Z1D)*deltaS

def intT1D(K1D):
    return np.sum(K1D)*deltaT

def intT2D(K):
    return np.sum(K)*deltaT

def intS(K):
    return (np.sum(K,1)*deltaS).T

def u(t):
    return (math.pi / (t*t-4)**0.5) * (((t*t-4)**0.5 - t) / 2.)**2

def z(s):
    return 2 * math.cos(s)

@jit
def K(t, s, z_s):
    return math.cos(2*s)/(t + z_s)

def dKdz(t, s, z_s):
    result = (z_s*t-2) / (t+z_s)**2
    return result

def dz(z_s):
    dz=np.zeros(N+1)
    for i in range(1,N):
        _s = s[i] - (b-a)/(2*N)
        s_ = s[i] + (b-a)/(2*N)
        _z = z_s[i-1] + (z_s[i]-z_s[i-1])*(_s-s[i-1])/(s[i]-s[i-1])
        z_ = z_s[i+1] + (z_s[i]-z_s[i+1])*(s_-s[i+1])/(s[i]-s[i+1])
        dz[i] = (z_ - _z)
    _s = s[0] - (b-a)/(2*N)
    s_ = s[0] + (b-a)/(2*N)
    _z = z_s[1] + (z_s[0]-z_s[1])*(_s-s[1])/(s[0]-s[1])
    z_ = z_s[1] + (z_s[0]-z_s[1])*(s_-s[1])/(s[0]-s[1])
    dz[0] = (z_ - _z)
    _s = s[N] - (b-a)/(2*N)
    s_ = s[N] + (b-a)/(2*N)
    _z = z_s[N-1] + (z_s[N]-z_s[N-1])*(_s-s[N-1])/(s[N]-s[N-1])
    z_ = z_s[N-1] + (z_s[N]-z_s[N-1])*(s_-s[N-1])/(s[N]-s[N-1])
    dz[N] = (z_ - _z)
    return dz*N/(b - a)

def myPsi(z_s):
    result = B*(b-a)*(b-a+1)*(1+intS1D(np.power(z_s,2)))
    return result

def myB():
    result = h * h*(h+2*ck)*(d-c)+((sigma*(h+2*ck)+h*(sigma+2*cu))*(d-c))/(b-a)+(sigma*(sigma+2*cu)*(d-c))/(b-a)**2;
    return result

def u_approx(t, theta):
    return u(t) * (1 + epsilon*theta)

def I_approx(z_s):
    uA_t = np.array([u_approx(t[i], theta[i]) for i in range(numStepT+1)])
    KA_ts = np.array([[K(t[i],s[j],z_s[j]) for j in range(numStepS+1)] for i in range(numStepT+1)])
    two = intS(KA_ts)-uA_t
    three = np.power(two, 2)
    result = intT1D(three)
    return result    

def L_delta(z_s):
    return I_approx(z_s) + myPsi(z_s)

def find_lambda_delta(z_s):
    lambda_delta = sc.minimize(L_delta, np.ones(N+1))
    return L_delta(lambda_delta.x)

def rho(z_alpha):
    uA_t = np.array([u_approx(t[i], theta[i]) for i in range(numStepT+1)])
    KA_ts = np.array([[K(t[i],s[j],z_alpha[j]) for j in range(numStepS+1)] for i in range(numStepT+1)])
    two = intS(KA_ts) - uA_t
    three = np.power(two, 2)
    result = intT1D(three) - lambda_delta - myPsi(z_alpha)
    return result
    
def M_alpha_delta(z_s):
    uA_t = np.array([u_approx(t[i], theta[i]) for i in range(numStepT+1)])
    KA_ts = np.array([[K(t[i],s[j],z_s[j]) for j in range(numStepS+1)] for i in range(numStepT+1)])
    two = intS(KA_ts) - uA_t
    three = np.power(two, 2)
    res1=intT1D(three)
    result = res1 + alpha*intS1D(np.power(z_s,2))
    return result

def DM_alpha_delta(z_s):
    f = z_s - dz(dz(z_s))
    uA_t = np.array([u_approx(t[i], theta[i]) for i in range(numStepT+1)])
    KA_ts = np.array([[K(t[i],s[j],z_s[j]) for j in range(numStepS+1)] for i in range(numStepT+1)])
    DK_txi = np.array([[dKdz(t[i],s[j],z_s[j]) for j in range(numStepS+1)] for i in range(numStepT+1)])
    two = np.array(intS(KA_ts) - uA_t)
    One = np.array([[DK_txi[i][j]*two[i] for j in range(numStepS+1)] for i in range(numStepT+1)])
    result = 2*intT2D(One) + alpha*f
    return result

def IterationProcess():
    global beta
    beta = 1
    z_k = z_k_init.copy()    
    while np.linalg.norm(DM_alpha_delta(z_k)) > 1e-5:
        #print("i ", i, "grad norm = ", np.linalg.norm(DM_alpha_delta(z_k)))
        beta = 1
        z_k_1 = z_k - beta * DM_alpha_delta(z_k)
        z_k_1[0] = 2
        z_k_1[N] = -2
        M_new = M_alpha_delta(z_k_1)
        M_old = M_alpha_delta(z_k)
        while M_new>M_old:
            beta = beta / 2
            z_k_1 = z_k - beta * DM_alpha_delta(z_k)
            z_k_1[0] = 2
            z_k_1[N] = -2
            M_new = M_alpha_delta(z_k_1)
            if beta < 1e-10:
                break
        z_k = z_k - beta * DM_alpha_delta(z_k)
        z_k[0] = 2
        z_k[N] = -2
        if (np.abs(M_old-M_new) < 1e-7):
            break 
    result = z_k
    return result

#начнем искать alpha
B = myB()
alpha = 0.001
alphaBegin = -1e-5
alphaEnd = 1
temp = 1
z_alpha = [z(i) for i in s]
lambda_delta = find_lambda_delta(z_alpha)
#lambda_delta = 0.000828092546412
print("lambda_delta", lambda_delta)

while False:
    alpha = (alphaBegin + alphaEnd) / 2
    z_alpha = IterationProcess()
    temp = rho(z_alpha)
    print("alpha " + str(alpha) + " rho " + str(temp))
    
    if abs(temp) < 1e-6:
        break
    
    if temp>0:
        alphaEnd = alpha
    else:
        alphaBegin = alpha
        
#alpha_arr = np.arange(0.17,0.18, 1e-4)
#temp_arr = []
#for al in alpha_arr:
#    alpha = al
#    z_alpha = IterationProcess()
#    temp = rho(z_alpha)
#    temp_arr.append(temp)
    
#plt.plot(alpha_arr, temp_arr)
#plt.show()
    
alpha = 0
zbrs = [z(i) for i in s]
beta = 1
z_k = z_k_init.copy()
M_new = 0
M_old = M_alpha_delta(z_k)
i=1
while np.linalg.norm(DM_alpha_delta(z_k)) > 1e-5 or (np.abs(M_old - M_new) > 1e-7):
    print("i ", i, "grad norm = ", np.linalg.norm(DM_alpha_delta(z_k)))
    beta = 10
    i = i + 1
    z_k_1 = z_k - beta * DM_alpha_delta(z_k)
    z_k_1[0] = 2
    z_k_1[N] = -2
    M_new = M_alpha_delta(z_k_1)
    M_old = M_alpha_delta(z_k)
    while M_new>M_old:
        beta = beta / 2
        z_k_1 = z_k - beta * DM_alpha_delta(z_k)
        z_k_1[0] = 2
        z_k_1[N] = -2
        M_new = M_alpha_delta(z_k_1)
        if beta < 1e-10:
            break
    z_k = z_k - beta * DM_alpha_delta(z_k)
    z_k[0] = 2
    z_k[N] = -2

    
plt.plot(s, zbrs)
plt.plot(s, z_k)
plt.plot(s, z_k_init)
plt.legend(['exact', 'current', 'init'])
plt.show()
