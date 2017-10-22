# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 19:11:05 2017

@author: СпелеВова
"""
import numpy as np

def sort_bubble(lst):
    res = lst.copy()
    n = len(res)
    for i in range(n):
        for j in range(n - 1, i , -1):
            if res[j] < res[j - 1]:
                res[j], res[j - 1] = res[j - 1], res[j]
    return res

def sort_insert(lst):
    res = lst.copy()
    n = len(res)
    for i in range(1, n):
        while i > 0 and res[i] < res[i - 1]:
            res[i], res[i - 1] = res[i - 1], res[i]
            i -= 1
    return res

def sort_select(lst):
    res = lst.copy()
    for i in range(len(res)):
        i_min = i
        for j in range(i + 1, len(res)):
            if res[j] < res[i_min]:
                i_min = j
        res[i_min], res[i] = res[i], res[i_min]
    return res

def sort_merge(lst):
    res = lst.copy()
    
    if len(res) > 1:
        mid = len(res) // 2
        left = res[:mid].copy()
        right = res[mid:].copy()
        left = sort_merge(left)
        right = sort_merge(right)
        i, j, k = 0, 0, 0
        while i < len(left) and j < len(right):
            if left[i] < right[j]:
                res[k] = left[i]
                i += 1
            else:
                res[k] = right[j]
                j += 1
            k += 1
        while i < len(left):
            res[k] = left[i]
            i += 1
            k += 1
        
        while j < len(right):
            res[k] = right[j]
            j += 1
            k += 1
    return res

begin = np.random.randint(1,100,10)
print(sort_bubble(begin))
print(sort_insert(begin))
print(sort_select(begin))
print(sort_merge(begin))
