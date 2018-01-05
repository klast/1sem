# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 15:14:12 2018

@author: spelevova
"""
import numpy as np
import math
import scipy.optimize as sc
from numpy import linalg
from scipy.integrate import trapz
from scipy.special import gamma
import matplotlib.pyplot as plt
import random

epsilon = 1e-3
N = 100
a = 0.
b = math.pi
c = 5
d = 6
numStepT = 100
numStepS = 100
deltaT = (d-c)/numStepT
deltaS = (b-a)/numStepS
t = np.arange(c, d, deltaT)
s = np.arange(a, b, deltaS)
theta = np.random.uniform(-1, 1, t.shape[0])
alpha=1e-8
sigma = 1e-4
h=0
sigma = epsilon*2
ck = 0.02986315562
cu = -0.035714285714

def intS1D(Z1D):
    return np.sum(Z1D)*deltaS

def intT1D(K1D):
    return np.sum(K1D)*deltaT

def intT2D(K):
    return np.sum(K)*deltaT

def intS(K):
    return (np.sum(K,1)*deltaS).transpose()

def u(t):
    return (math.pi / (t*t-4)**0.5) * (((t*t-4)**0.5 - t) / 2.)**2

def z(s):
    return 2 * math.cos(s)

def K(t, s, z_s):
    return math.cos(2*s)/(t + z_s)

def dKdz(t, s, z_s):
    #return 3 * (s - 1.) / ((t + z) ** 4)
    #return (-2.*t)/((1+t*z)**3)
    #return 3*(s-1.)/((t+z)**4)
    if (z_s != 0 ):
        #result = t*(z**(t-1))
        result = (z_s*t-2) / (t+z_s)**2
    if (z_s == 0):
        result = 0
   # print ("t = ", t, "Z=", z, "res = ", result)
    return result

def myPsi(z_s):
    result = B*(b-a)*(b-a+1)*(1+intS1D(np.power(z_s,2)))
    return result

def myB():
    result = h * h*(h+2*ck)*(d-c)+((sigma*(h+2*ck)+h*(sigma+2*cu))*(d-c))/(b-a)+(sigma*(sigma+2*cu)*(d-c))/(b-a)**2;
    return result

def u_approx(t, theta):
    return u(t) * (1 + epsilon*theta)

def I_approx(z_s):
    uA_t = [u_approx(t[i], theta[i]) for i in range(numStepT)]
    KA_ts = [[K(t[i],s[j],z_s[j]) for j in range(numStepS)] for i in range(numStepT)]
    two = intS(KA_ts)-uA_t
    three = np.power(two, 2)
    result = intT1D(three)
    return result    

def L_delta(z_s):
    return I_approx(z_s) + myPsi(z_s)

def find_lambda_delta(z_s):
    bounds = [(-10, 10) for i in range(numStepS)]
    lambda_delta = sc.differential_evolution(L_delta, bounds, tol=1e-5)
    return L_delta(lambda_delta.x)

def rho(z_alpha):
    uA_t = [u_approx(t[i], theta[i]) for i in range(numStepT)]
    KA_ts = [[K(t[i],s[j],z_alpha[j]) for j in range(numStepS)] for i in range(numStepT)]
    two = intS(KA_ts) - uA_t
    three = np.power(two, 2)
    result = intT1D(three) - lambda_delta - myPsi(z_alpha)
    return result
    
def M_alpha_delta(z_s):
    uA_t = [u_approx(t[i], theta[i]) for i in range(numStepT)]
    KA_ts = [[K(t[i],s[j],z_s[j]) for j in range(numStepS)] for i in range(numStepT)]
    two = intS(KA_ts) - uA_t
    three = np.zeros(two, 2)
    res1=intT1D(three)
    result = res1 + alpha*intS1D(np.power(z_s,2))
    return result

def DM_alpha_delta(z_s):
    uA_t = [u_approx(t[i], theta[i]) for i in range(numStepT)]
    KA_ts = [[K(t[i],s[j],z_s[j]) for j in range(numStepS)] for i in range(numStepT)]
    DK_txi = [[dKdz(t[i],s[j],z_s[j]) for j in range(numStepS)] for i in range(numStepT)]
    two = intS(KA_ts) - uA_t
    One = [[DK_txi[i][j]*two[i] for j in range(numStepS)] for i in range(numStepT)]
    result = 2*intT2D(One) + 2*alpha*z_s
    return result

def IterationProcess():
    z_k = 2*s/math.pi
    zbrs = z(s)
    z_k[0] = zbrs[0]
    z_k_init = z_k.copy()
    M_new = 0
    M_old = M_alpha_delta(z_k)
    
    
    for i in range(N):
        z_k_1 = z_k - beta * DM_alpha_delta(z_k)
        M_new = M_alpha_delta(z_k_1)
        z_k = z_k_1
        M_old=M_new
    result = z_k
    return result

#начнем искать alpha
B = myB()
beta = 0
alpha = 0
alphaBegin = 0
alphaEnd = 0.1
temp = 1
z_alpha = [z(i) for i in s]
lambda_delta = find_lambda_delta(z_alpha)
while True:
    alpha = (alphaBegin + alphaEnd) / 2
    z_alpha = IterationProcess()
    temp = rho(z_alpha)
    print("alpha " + alpha + " rho " + temp)
    
    if abs(temp) < 1e-6:
        break
    
    if temp>0:
        alphaEnd = alpha
    else:
        alphaBegin = alpha
        
z_k = 2*s/math.pi
zbrs = z(s)
z_k[0] = zbrs[0]
z_k_init = z_k
M_new = 0
M_old = M_alpha_delta(z_k)
for i in range(N):
    print("i ", i)
    z_k_1 = z_k - beta * DM_alpha_delta(z_k)
    M_new = M_alpha_delta(z_k_1)
    z_k = z_k_1
    M_old = M_new
