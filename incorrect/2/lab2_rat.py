
# coding: utf-8

# In[236]:

import numpy as np
from scipy.optimize import fmin, fsolve, minimize
from scipy import integrate
from numpy import linalg, random
import matplotlib.pyplot as plt
from matplotlib import cm

# In[237]:

N = 100
a = 0.
b = 1.
c = 1.
d = 2.
alpha=1e-5
sigma = 1e-6#.001
h=1e-6#.001
ck=1
cn=0.25


# In[238]:

s = np.arange(a, b+(b-a)/N, (b-a)/N)
t = np.arange(c, d+(d-c)/N, (d-c)/N)


# In[239]:

def K(s,t,z):
    return (1.-s)/((t+z)**3)

def dKdz(s,t,z):
    return 3*(s-1.)/((t+z)**4)

def u(t):
    return 1./(2*t*t*(1+t))

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

def omega(z, a, b):
    f = z**2 + dz(z)**2
    return integral(f,a,b)

def psi(z):
    B = h*(h+2*ck)*(d-c)+(sigma*(h+2*ck)+h*(sigma+2*ck))*(d-c)/(b-a) + sigma*(sigma+2*cn)*(d-c)/((b-a)**2)
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

def L(z):
    return J(z) + psi(z)

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
    return 2*dk + alpha*f

def rho(alpha):
    z0=minimize(M, np.ones(N+1), args=(alpha,))
    return J(z0.x) - psi(z0.x) - ldelta


U = u(t)
z = np.ones(N+1)
beta = 1.



ldelta = minimize(L, np.ones(N+1)).fun


i=0
while (linalg.norm(gradM(z,alpha))>1e-5):
    beta=5.
    i=i+1
    z_ = z-beta*gradM(z, alpha)
    while(M(z_, alpha)>M(z, alpha)):
        beta=beta/2
        z_ = z-beta*gradM(z, alpha)
        if (beta<1e-6):
            break
    if (beta<1e-6):
        break
    z=z-beta*gradM(z,alpha)
    z[0]=s[0]
    z[N]=s[N]

print(i)


print(linalg.norm(gradM(z,alpha)))


fig = plt.figure()
l1, = plt.plot(np.ones(N+1),  label='z0')
l2, = plt.plot(z,  label='z')
l3, = plt.plot(s,  label='z*')
plt.legend(handles=[l1, l2, l3])

kk = np.zeros((N + 1, N + 1))
ii = np.zeros(N + 1)
intk = np.zeros(N + 1)
for i in range(0, N + 1):
    for j in range(0, N + 1):
        kk[i][j] = K(s[j], t[i], z[j])
    ii[i] = integral(kk[i], a, b)

fig = plt.figure()
l1, = plt.plot(U-ii,  label='U(t) - int(K)')
plt.legend(handles=[l1])
plt.show()


# In[221]:

print(np.sqrt(np.dot(z-s,z-s)))


# In[222]:

print(np.sqrt(np.dot(z-s,z-s))/np.sqrt(np.dot(s,s)))

print(np.sqrt(np.dot(np.ones(N+1)-s,np.ones(N+1)-s)))
print(np.sqrt(np.dot(np.ones(N+1)-s,np.ones(N+1)-s))/np.sqrt(np.dot(s,s)))



