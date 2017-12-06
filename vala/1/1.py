import numpy as np

N = 500

def gauss(A, b):

    if abs(np.linalg.det(A)) < 1.0E-6:
        return -1.0

    for column in range(0, N):
        print(column)
        main_elem = A[column, column]
        A[column, :] = A[column, :] / main_elem
        b[column] = b[column] / main_elem

        for row in range(column + 1, N):
            b[row] = b[row] - A[row, column] * b[column]
            A[row, :] = A[row, :] - A[row, column] * A[column, :]

    x_gen = np.zeros(N)
    x_gen[N - 1] = b[N - 1]
    for row in range(N - 2, -1, -1):
        x_gen[row] = b[row] - A[row, (row + 1):] @ x[(row + 1):]

    return x_gen


A = np.random.rand(N, N)
x = np.random.rand(N)
b = A @ x
x_0 = gauss(A, b)
print(np.linalg.norm(x - x_0))
