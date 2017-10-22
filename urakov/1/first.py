# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""
import math

def swap(st, i, j):
    st[i], st[j] = st[j], st[i]
    return st

def next(current):
    i = len(current) - 2
    res = list(current)
    while i >= 0 and current[i] > current[i + 1]:
        i -= 1
    if i < 0:
        return []
    j = len(current) - 1
    while current[j] < current[i]:
        j -= 1
    res = swap(res, i, j)    
    res[i+1:] = res[:i:-1]
    return "".join(res)

def prev(current):
    i = len(current) - 2
    res = list(current)
    while i >= 0 and current[i] < current[i + 1]:
        i -= 1
    if i < 0:
        return []
    j = len(current) - 1
    while current[j] > current[i]:
        j -= 1
    res = swap(res, i, j)
    res[i+1:] = res[:i:-1]
    return "".join(res)

def num_on_str(current):
    res = 0
    n = len(current)
    for i in range(n):
        k = 0
        for elem in current[i+1:]:
            if elem < current[i]:
                k += 1
        res += math.factorial(n - 1 - i) * k
    return res + 1
    
def str_on_num(begin, num):
    num +=1
    res = []
    unused = list(begin)
    n = len(begin)
    while n > 0:
        n -= 1
        k = num // math.factorial(n)
        res.append(unused[k])
        unused.pop(k)
        num %= math.factorial(n)
    return "".join(res)
    
res = next('12354')
print(res)
res = prev(res)
print(res)
res = num_on_str('12354')
print(res)
res = str_on_num('12345', 5)
print(res)
        

