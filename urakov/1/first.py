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
    i = len(prev) - 2
    res = prev
    while i>=0 and prev[i]>prev[i + 1]:
        i-=1
    if i < 0:
        return []
    j = len(prev) - 1
    while prev[j] < prev[i]:
        j-=1
    res = swap(res, i, j)
    res = list(res)
    res[i+1:] = res[:i:-1]
    return "".join(res)

res = next('abced')
print(res)
        

