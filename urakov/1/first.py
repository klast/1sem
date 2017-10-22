# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""
import numpy

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
    
res = next('abcd')
print(res)
res = prev(res)
print(res)
        

