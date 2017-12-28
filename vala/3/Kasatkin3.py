import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=4)
B=np.random.rand(3,3)-0.5
BT=np.transpose(B)
E=np.identity(3)
k=0.7
A = BT @ B + k*E
print("A=", A)
x0 = np.random.rand(3)
print("exact solution", x0)
b = A @ x0
print("right side = ", b)
# Steepest Decent Algorithm
print("Steepest")
x1=[0,0,0]
r = b - (A @ x1)
p = A @ r
#for i in range(0, 5):
i=0
func = 100
Res = np.array([np.linalg.norm(b-(A @ x1))/np.linalg.norm(b)])
while (np.linalg.norm(b-(A @ x1))/np.linalg.norm(b)>0.000001):
    i=i+1
    alpha = np.dot(r,r)/np.dot(p,r)
    x1=x1 + alpha*r
    r=r - alpha*p
    p = A @ r
    funcprev = func
    func = np.dot(A @x1, x1) - 2*np.dot(b,x1)
  #  print(x1)
    Res = np.append(Res , [np.linalg.norm(b-(A @ x1))/np.linalg.norm(b)])
    print("N = %d , x-x0/x = %f , b-Ax = %f , b-Ax/b = %f" % (i ,np.linalg.norm(x1-x0)/np.linalg.norm(x1), np.linalg.norm(b - (A @ x1)) ,np.linalg.norm(b-(A @ x1))/np.linalg.norm(b)))
    if func>funcprev:
        raise ValueError("Mistake", 100)
   # print("(Ax,x) - 2(b,x) = ", func)

print(x1)
# Residue graphic
N=np.arange(i+1)
print(Res)
#Res = np.log(Res)
print(N)
print(Res)
#plt.semilogy(N,Res)
plt.plot(N,Res)
plt.show()
#plt.semilogy(N,Res)
#plt.show()
# P-G Check
print("PG-Check")
alpha = np.dot(r,r)/np.dot(p,r)
x1=x1 + alpha*r
pg=np.dot(b- (A @ x1), r)
print(pg)



# Minimal Residual Iteration
print("Minimal Residual")
x2=[0,0,0]
r= b - (A @ x2)
p = A @ r
i=0
Residue = np.linalg.norm(b-(A @ x2))
while (np.linalg.norm(b-(A @ x2))/np.linalg.norm(b)>0.0001):
    i+=1
    alpha = np.dot(p,r)/np.dot(p,p)
    x2=x2 + alpha*r
    r=r - alpha*p
    p = A @ r
    ResiduePrev = Residue
    Residue = np.linalg.norm(b-(A @ x2))
    # Проверяем что невязка убывает
    if Residue>ResiduePrev:
        raise ValueError("Mistake",100)
    print("N = %d ,x-x0/x = %f , b-Ax = %f , b-Ax/b = %f" % (i, np.linalg.norm(x2-x0)/np.linalg.norm(x2), np.linalg.norm(b - (A @ x2)) ,np.linalg.norm( b-(A @ x2))/np.linalg.norm(b)))
print(x2)
# P-G Check
print("PG-Check")
alpha = np.dot(p,r)/np.dot(p,p)
x2=x2 + alpha*r
pg=np.dot(b- (A @ x2),A @ r)
print(pg)

#Conjugate Gradient Method
print("Conjugate")
x3=[0,0,0]
r= b - (A @ x3)
p = r
i=0
Residue = np.linalg.norm(b-(A @ x3))
while (np.linalg.norm(b-(A @ x3))/np.linalg.norm(b)>0.0001):
    i+=1
    alpha = np.dot(r,r)/np.dot(A @ p,p)
    x3=x3 + alpha*p
    rprev = r
    r=r - alpha*(A @ p)
    beta = np.dot(r , r)/np.dot(rprev, rprev)
    p = r + beta*p
    ResiduePrev = Residue
    Residue = np.linalg.norm(b-(A @ x3))
    # Проверяем что невязка убывает
    if Residue>ResiduePrev:
        raise ValueError("Mistake",100)
    print("N = %d ,x-x0/x = %f , b-Ax = %f , b-Ax/b = %f" % (i,np.linalg.norm(x3-x0)/np.linalg.norm(x3), np.linalg.norm(b - (A @ x3)) ,np.linalg.norm(b-(A @ x3))/np.linalg.norm(b)))
print(x3)