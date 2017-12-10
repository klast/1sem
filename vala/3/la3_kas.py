import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt
import math
import random
# Следующая строчка для того, чтобы генерировалась одна и та же матрица
np.random.seed(2)
# параметры графика
plt.subplots_adjust(hspace=0.4)
# параметры  вывода матрицы
np.set_printoptions(precision=4, linewidth=200)
# Эпсилоны для программ
eps_gal = 1e-8
eps = 1e-8
eps_3 = 1e-8
# 1 задание
n = 500


k = 10
# Единичная матрица
E = np.identity(n)
# Генерируем случайную матрицу B, где b_ij \in (-0.5:0.5)
B = (2*np.random.random_sample((n, n))-1) / 2
A = B.transpose()*B + k*E

# 2 задание

b = np.random.rand(n)
x_nach = np.zeros(n)
x0 = x_nach.copy()
x = lg.solve(A.copy(), b.copy())
r = b - A.dot(x_nach)
# создаем массив для 3 задания
graph = []

# для 3 задания номер итерации для построения графика
k = 0
# для 3 задания закидываем в graph текущий номер итерации и относительную невязку
graph.append([ lg.norm(A.dot(x_nach) - b)/lg.norm(b)])
# для 5 задания генерируем массив F
F = []
F.append((A.dot(x_nach)).dot(x) - 2*b.dot(x_nach))
while True:
    k += 1
    print(k)
    alpha = r.dot(r)/((A.dot(r)).dot(r))
    x_nach += alpha * r
    # 4 задание
    Petrov_Gal = (b - A.dot(x_nach)).dot(r)/lg.norm(r)
    
    if abs(Petrov_Gal) > eps_gal:
        print(k)
        print("Petrov-Galerkin is not complited")
    # 4
    r = r -alpha*A.dot(r)
    # 3 задание
    graph.append([ lg.norm(A.dot(x_nach) - b) / lg.norm(b)])
    # 3
    #5 задание
    F.append((A.dot(x_nach)).dot(x_nach) - 2 * b.dot(x_nach))
    if F[k] > F[k-1]:
        print("Error!!!!!!")
        exit(0)
    # 5
    flag = True
    for i in range(n):
        if r[i] < eps:
            continue
        else:
            flag = False
            break
    if flag :
        break
x_grad = x_nach
print("Grad_method ", x_grad)
print("Otklonenie ot tochnogo resheniya", lg.norm(x-x_nach))
print("norm nevyazki", lg.norm(b-A.dot(x_nach)))
print("Otklonenie reshenia", lg.norm(x_nach-x0)/lg.norm(x_nach))
print("Otnositelnaya nevyazka", lg.norm(b-A.dot(x_nach))/lg.norm(b))
#3 задание
graph1 = np.array(graph)

print("graph", graph1[:,:])

#5 задание
F1 = np.array(F)
print("F", F1)



#7 задание
k = 0
graph = []
x_nach = np.zeros(n)
r1 = b - A.dot(x_nach)
norm_r = []
# 10 задание
q_prog = []
mu = np.min(lg.eigvals(A.transpose()+A),axis=0)/2
sigma = lg.norm(A)
q = math.sqrt(1-mu*mu/sigma/sigma)
p = A.dot(r1)
graph.append([ lg.norm(A.dot(x_nach) - b) / lg.norm(b)])
while True:
    # 9 задание
    norm_r.append(lg.norm(r1))
    # 9
    k += 1
    alpha = p.dot(r1)/p.dot(p)
    x_nach = x_nach + alpha * r1
    graph.append([ lg.norm(A.dot(x_nach) - b) / lg.norm(b)])
    # 8 задание
    Petrov_Gal = (b - A.dot(x_nach)).dot(r)
    if abs(Petrov_Gal) > eps_gal:
        print("Petrov-Galerkin is not complited")
    # 8
    r_prom = lg.norm(r1)
    r1 = r1 - alpha*p
    q_prog.append( lg.norm(r1)/r_prom)
    p = A.dot(r1)

    flag = True
    for i in range(n):
        if r1[i] < eps:
            continue
        else:
            flag = False
            break
    if flag :
        break
x_min = x_nach
graph2 = np.array(graph)
print("Min nevyazka", x_nach)
print("Norma raznosti grad and min", lg.norm(x_grad-x_min))
# 7
# 9 задание
print("Norma nevyazki", np.array(norm_r))
# 10 задание
print("q_prog",np.array(q_prog))
print("q", q)

# 11, 12 задания
k = 0
x_nach = np.zeros(n)
r1 = b - A.dot(x_nach)
graph = []
graph.append([ lg.norm(A.dot(x_nach) - b)/lg.norm(b)])
mas_r = []
p = r1
mas_p = []
for i in range(n):
    k += 1
    mas_p.append(p)
    mas_r.append(r1)
    alpha = r1.dot(r1)/(A.dot(p)).dot(p)
    x_nach = x_nach + alpha*p
    graph.append([ lg.norm(A.dot(x_nach) - b)/lg.norm(b)])
    r_prom = r1
    r1 = r1 - alpha*A.dot(p)
    betta = r1.dot(r1)/r_prom.dot(r_prom)
    p = r1 + betta*p
    flag = True
    for i in range(n):
        if r1[i] < eps_3:
            continue
        else:
            flag = False
            break
    if flag :
        break
graph3 = np.array(graph)
x_sopr = x_nach
mas_p = np.array(mas_p)
mas_r = np.array(mas_r)
print("x_sopr",x_sopr)
# 13 задание
print("norm x_sopr and x_grad", lg.norm(x_sopr-x_grad))
print("KKKK", k)
number1 = random.randint(0, math.ceil(k/2)-1)
number2 = random.randint(math.ceil(k/2), k-1)
# 14 задание
print("(A*p_i,p_j)",(A.dot(mas_p[number1, :])).dot(mas_p[number2, :]))
print("(A*p_i,p_i+1)", (A.dot(mas_p[number1, :])).dot(mas_p[number1+1, :]))
print("(r_i,r_j)", mas_r[number1, :].dot(mas_r[number2, :]))

# 15 задание
x_nach = np.zeros(n)
r1 = b - A.dot(x_nach)
p = r1
mas_p = []
mas_r = []
for i in range(n):
    mas_p.append( p)
    mas_r.append( r1)
    alpha = r1.dot(p)/p.dot(A.dot(p))
    x_nach = x_nach + alpha*p
    r_prom = r1
    r1 = r1 - alpha*A.dot(p)
    for i in range(n):
        flag = True
        if r1[i] < eps_3:
            continue
        else:
            flag = False
            break
    if flag :
        break
    betta = -(A.dot(r1)).dot(p)/p.dot(A.dot(p))
    p = r1 + betta*p
    

x_sopr_mod = x_nach
print("x_sopr_mod",x_sopr_mod)
# 13 задание
print("norm x_sopr and x_grad", lg.norm(x_sopr-x_sopr_mod))
print("graph2", graph2)
print("graph3", graph3)
graph = np.array([graph1,graph2,graph3])
print("graph", graph)
plt.semilogy(graph1)
plt.semilogy(graph2)
plt.semilogy(graph3)
plt.show()