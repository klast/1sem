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

def merge(left, right):
    res = []
    i, j = 0, 0
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:
            res.append(left[i])
            i += 1
        else:
            res.append(right[j])
            j += 1
        res += left[i:]
        res += right[j:]
    return res

def sort_merge(lst):
    res = lst.copy()
    if len(res) <= 1:
        return res
    else:
        mid = len(res) // 2
        left = res[:mid]
        right = res[mid:]
    return merge(sort_merge(left), sort_merge(right))


begin = np.random.randint(1,100,10)
print(sort_bubble(begin))
print(sort_insert(begin))
print(sort_merge(begin))
