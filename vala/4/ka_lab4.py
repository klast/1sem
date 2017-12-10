
# coding: utf-8

# In[150]:


import numpy as np
import numpy.linalg as la
np.set_printoptions(edgeitems=7,precision=3) # 3 знако после запятой, 7 элементов в матрице по краям до многоточий
np.core.arrayprint._line_width = 200 # широкие строки в браузере для широкого экрана

# n=100
# dx=float(1/(n-1))
# m=100
# dy float(1/(m-1))
n=30
P=np.zeros((n,n),dtype=np.int32)-1
#P[2:9,3:8]=1

def addRectangle(P,x0,y0,x1,y1):
    i0=round(x0*(n-1))
    j0=round(y0*(n-1))
    i1=round(x1*(n-1))
    j1=round(y1*(n-1))
    P[i0:i1,j0:j1] = 1

def addEllipse(P,x0,y0,r1,r2):
    x=0
    y=0
    for i in range(n):
        for j in range(n):
                x=i/n
                y=j/n
                if ((x-x0)*(x-x0)/r1/r1+(y-y0)*(y-y0)/r2/r2 <=1):
                    P[i,j]=1
    
    
addRectangle(P, 0.1,0.1,0.4,0.4)
addEllipse(P, 0.5,0.5,0.2,0.4)
P




# In[151]:


k=0 # текущий номер переменной
# Nums=np.zeros(P.shape,dtype=np.int32)-1
for i in range(n):
    for j in range(n):
        if P[i,j]==1:
            P[i,j]=k;
            k += 1
            

print("Номера переменных, соответствующих узлам сетки:")
print(P)
m=k


# In[152]:


m


# In[153]:


A=np.zeros((m,m))
for i in range(n):
    for j in range(n):
        if P[i,j]==-1:
            continue
        k = P[i,j]
        A[k,k] = -4
        if P[i+1,j] != -1:
            A[k,P[i+1,j]] = 1
        if P[i-1,j] != -1:
            A[k,P[i-1,j]] = 1
        if P[i,j-1] != -1:
            A[k,P[i,j-1]] = 1
        if P[i,j+1] != -1:
            A[k,P[i,j+1]] = 1
A       
        


# In[154]:


A[0,8]


# In[155]:


A[0,9]


# In[156]:


get_ipython().magic('time k,V=la.eig(A)')
get_ipython().magic('time k,V=la.eigh(A)')

def vectorToGrid(x, Nums):
    n=Nums.shape[0]
    U=np.zeros(Nums.shape)
    for i in range(n):
        for j in range(n):
            if Nums[i,j]== -1:
                U[i,j]=0
            else:
                
                U[i,j]= x[Nums[i,j]]
    #print(U)
    return U

# print(vectorToGrid(np.random.rand(???),Nums))

idx=np.argsort(k) ### idx - номера элементов массива k в порядке возрастания значений
idx=idx[::-1] ### idx в обратном порядке
xx = k[idx] ### отсортированный массив

k[idx[0]] ### ближайшее к нулю собственное значение
k[idx[1]] ### второе от нуля собственное значение


x = V[:,idx[0]] ###-й столбец, собственный вектор, соответствующий lambda_2

U=  vectorToGrid(x,P) 
get_ipython().magic('matplotlib notebook')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

X, Y = np.meshgrid(np.linspace(0,1,n), np.linspace(0,1,n))
fig = plt.figure(figsize=(10,8))
ha = fig.add_subplot(111, projection='3d')
ha.plot_surface(X,Y,U, cmap=cm.hsv)


# In[158]:


x1 = V[:,idx[1]] ###-й столбец, собственный вектор, соответствующий lambda_2

U1=  vectorToGrid(x1,P)
fig1 = plt.figure(figsize=(10,8))
ha1 = fig1.add_subplot(111, projection='3d')
ha1.plot_surface(X,Y,U1, cmap=cm.hsv)

