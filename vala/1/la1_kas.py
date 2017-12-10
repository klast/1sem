import numpy as np
import numpy.linalg as lg


n = 10
np.set_printoptions(precision=2, linewidth=400)

b = np.random.rand(n)

A = np.random.rand(n, n)
x = np.zeros(n)
x = lg.solve(A.copy(), b.copy())
print()
# Спавка: Метод гаусса в выбором максимального элемента в столбце
# Справка: A[i:,:] вернет "cрезку", т.е. все строчки, начиная с i-ой, и все столбы матрицы A
# Если A[i:,i:], то вернет  все строчки, начиная с i-ой, и все столбы, начиная с i-ого, матрицы A
# Бежим по всем строчкам матрицы А
for i in range(n):
    # Если не полседння строчка, то выполняем ( возможно экномия времени будет)
    if i != n-1:
        # Находим максимальный элемент в столбце
        index_max = A[i:, i].argmax()+i
        # Если он 0, то определитель матрицы ), значит он а несовметсна
        if A[index_max, i] == 0:
            print("Error : det(A)=0!")
            exit(0)
        # Если найденный индекс максимального элемента не совпадает с текущим индексом i, то меняем строчки  матрицы и
        # вектора b местами
        if i != index_max:
            h = A[index_max, :].copy()
            A[index_max, :] = A[i, :]
            A[i, :] = h

            c = b[index_max].copy()
            b[index_max] = b[i]
            b[i] = c
        # Выполняем нормировку по i элементу строки i
        b[i] = b[i] / A[i, i]
        A[i, i:] = A[i, i:]/A[i, i]
        # Зануляем все элементы матрицы А под A[i,i]
        for j in range(i+1, n, 1):
            b[j] = b[j] - b[i] * A[j, i]
            A[j, i:] = A[j, i:] - A[j, i]*A[i, i:]
    # Если полседння строчка, то выполняем
    else:
        b[i] = b[i] / A[i, i]
        A[i, i] = 1
# Обратный ход метода Гаусса
for i in (range(n-2, -1, -1)):
    b[i] = b[i] - A[i, i+1:].dot(b[i+1:])
    A[i, i + 1:] = 0

print(b)
print(x)