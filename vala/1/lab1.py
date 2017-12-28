import numpy as np
import numpy.linalg as lg


n = 1000
np.set_printoptions(precision=4, linewidth=400)

b = np.random.rand(n)

A = np.random.rand(n, n)
x = np.zeros(n)
x = lg.solve(A.copy(), b.copy())
print()
for i in range(n):
    if i != n-1:
        index_max = A[i:, i].argmax()+i
        if A[index_max, i] == 0:
            print("Error : det(A)=0!")
            exit(0)
        if i != index_max:
            h = A[index_max, :].copy()
            A[index_max, :] = A[i, :]
            A[i, :] = h

            c = b[index_max].copy()
            b[index_max] = b[i]
            b[i] = c
        b[i] = b[i] / A[i, i]
        A[i, i:] = A[i, i:]/A[i, i]
        for j in range(i+1, n, 1):
            b[j] = b[j] - b[i] * A[j, i]
            A[j, i:] = A[j, i:] - A[j, i]*A[i, i:]
    else:
        b[i] = b[i] / A[i, i]
        A[i, i] = 1

for i in (range(n-2, -1, -1)):
    b[i] = b[i] - A[i, i+1:].dot(b[i+1:])
    A[i, i + 1:] = 0

#print(b)
#print(x)
print(lg.norm(x-b))