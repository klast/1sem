
# coding: utf-8

# In[236]:

import numpy as np
import math
from scipy.optimize import  minimize
from numpy import linalg
from scipy.integrate import trapz
from scipy.special import gamma
import matplotlib.pyplot as plt
import random


# In[237]:

N = 25
a = 0.
b = math.pi/2
c = 0
d = 2
alpha=1e-8
sigma = 0#.001
h=0.1
#h=0.00018
#ck=1
ck = 1
#cn=0.25
cn = 1.571


# In[238]:

s = np.arange(a, b+(b-a)/N, (b-a)/N)
solution = []
for i in range(0,N+1):
    solution.append(math.cos(s[i]))
t = np.arange(c, d+(d-c)/N, (d-c)/N)


# In[239]:

def K(s,t,z):
    #return (1. - s) / ((t + z) ** 3)
    #return (1.)/((1+t*z)**2)
    return z**t

def dKdz(s,t,z):
    #return 3 * (s - 1.) / ((t + z) ** 4)
    #return (-2.*t)/((1+t*z)**3)
    #return 3*(s-1.)/((t+z)**4)
    if (z != 0 ):
        result = t*(z**(t-1))
    if (z == 0):
        result = 0
   # print ("t = ", t, "Z=", z, "res = ", result)
    return result

def u(t):
    #return (math.pi)/((1-t**2)**(1.5))
    #return 1./(2*t*t*(1+t))
    return ((math.pi)/2)*(gamma((t+1)/2)/gamma(t/2 + 1))



def dz(z):
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

def integral(f,a_,b_):
    s = 0
    for i in range(0,N):
        s = s + f[i+1] + f[i]
    return s*(b_ - a_)/(2*N)

def Ur(z,t):
    a1 = integral(z**t,a,b) - (1/N)*random.uniform(0,1)*1e-4
    return a1

def omega(z, a, b):
    f = z**2 + dz(z)**2
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

def gradM(z, alpha):
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
    result = 2*dk + alpha * f
    return result








# In[240]:

U = u(t)
z = np.ones(N+1)
beta=1.
i=0
ii=0
for ii in range(1, 2):
    ii=ii
    print("Numer ",ii)
    alpha = 76.65
    print("alpha = ", alpha)
    z = np.ones(N + 1)
    i=0
    beta = 0.1
    while (linalg.norm(gradM(z,alpha))>1e-5):
        print(linalg.norm(gradM(z,alpha)))
        beta = 0.00001
        i = i + 1
        z_ = z-beta*gradM(z, alpha)
        z_[0] = 1
        z_[N] = 0
        #print(z_)
        #print("z_ = ", z_)
        while(M(z_, alpha)>M(z, alpha)):
            beta=beta/2
            z_ = z-beta*gradM(z, alpha)
            z_[0] = 1
            z_[N] = 0
            if (beta<1e-10):
                break
        if (abs(M(z_, alpha)-M(z, alpha))<1e-5):
            break
        z=z-beta*gradM(z,alpha)
        z[0]=1
        z[N]=0
    #print(z)
    print(i," iterations")
    print("rho =", J(z) - psi(z))


for i in range(0,1):
    h =  0.0175
   # print("rho =", J(z) - psi(z), "  h = ", h)
# In[243]:

#print (linalg.norm(gradM(z,alpha)))
#print("z = ", z)
#print("s = ", s)
# In[244]:
abs = np.arange(a, b+(b-a)/N, (b-a)/N)
fig = plt.figure()
l1, = plt.plot(abs, np.ones(N+1)-0.5,  label='initial')
l2, = plt.plot(abs ,z,  label='z')
l3, = plt.plot(abs ,solution,  label='exact')
plt.legend(handles=[l1, l2, l3])
for jj in range(0, 10):
    t=0.2*jj
    #print("Nev = ", integral(z**t, a, b)-Ur(z,t), "t =", t)
plt.show()







