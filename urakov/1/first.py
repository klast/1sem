# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""
import numpy

def swap(st, i, j):
    st = list(st)
    st[i], st[j] = st[j], st[i]
    return ''.join(st)

def reverse(string, i, j):
    res = string[:i]
    res += string[:i-1:-1]
    return res

    

def next(prev):
    n = len(prev)
    for i in range(n-2, 0, -1):
        if prev[i] < prev[i+1]:
            min = i + 1
            for j in range(i + 1, n, 1):
                if(prev[j] < prev[min]) and (prev[j] > prev[i]):
                    min = j
            tmp = swap(prev, i, min)
            res = reverse(tmp, i + 1, n - 1)
            return res
    return None
res = next('13254')
print(res)
        

