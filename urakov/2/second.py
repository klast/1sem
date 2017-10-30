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

def sort_select(lst):
    res = lst.copy()
    n = len(res)
    for i in range(n):
        i_min = i
        for j in range(i + 1, n):
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
        if len(left) > 1:
            print("Beginning to sort left = ",left)
            left = sort_merge(left)
        if len(right) > 1:
            print("Beginning to sort right = ", right)
            right = sort_merge(right)
        i, j, k = 0, 0, 0
        print("Beginning to merge ", left ,'and ', right)
        while i < len(left) and j < len(right):
            if left[i] < right[j]:
                res[k] = left[i]
                i += 1
            else:
                res[k] = right[j]
                j += 1
            k += 1
            [ print(res[i],end=" ") for i in range(k)]
            print("",end = '\n')
            
        while i < len(left):
            res[k] = left[i]
            i += 1
            k += 1
            [ print(res[i],end=" ") for i in range(k)]
            print("",end = '\n')
        while j < len(right):
            res[k] = right[j]
            j += 1
            k += 1
            [ print(res[i],end=" ") for i in range(k)]
            print("",end = '\n')
        if len(left) > 1 or len(right) > 1:
            print('merged = ', res)
    return res

begin = np.random.randint(1,100,8)
print("unsorted = ", begin)
#print(sort_bubble(begin))
#print(sort_select(begin))
print('sorted_list = ', sort_merge(begin))
