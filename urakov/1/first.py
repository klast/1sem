# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""
import math

def next(current):
    i = len(current) - 2
    res = list(current)
    while i >= 0 and current[i] > current[i + 1]:
        i -= 1
    if i < 0:
        return None
    j = len(current) - 1
    while current[j] < current[i]:
        j -= 1
    res[i], res[j] = res[j], res[i] 
    res[i+1:] = res[:i:-1]
    return "".join(res)

def prev(current):
    i = len(current) - 2
    res = list(current)
    while i >= 0 and current[i] < current[i + 1]:
        i -= 1
    if i < 0:
        return None
    j = len(current) - 1
    while current[j] > current[i]:
        j -= 1
    res[i], res[j] = res[j], res[i] 
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
    return res
    
def str_on_num(begin, num):
    res = []
    unused = list(begin)
    n = len(begin)
    if num > math.factorial(n) - 1 or num < 0:
        return None
    while n > 0:
        n -= 1
        k = num // math.factorial(n)
        res.append(unused[k])
        unused.pop(k)
        num %= math.factorial(n)
    return "".join(res)

def print_next(begin):
    while begin != None:
        begin = next(begin)
        print(begin)

def print_prev(begin):
    while begin != None:
        begin = prev(begin)
        print(begin)
        
def print_nums(begin, way='FORWARD'):
    while begin != None:
        print(num_on_str(begin), begin)
        if way == 'FORWARD':
            begin = next(begin)
        elif way == 'BACKWARD':
            begin = prev(begin)
            
def print_str_on_num(begin, num, way = 'FORWARD'):
    current = str_on_num(begin, num)
    print(current, num)
    while True:
        current = str_on_num(begin, num)
        if current == None:
            return        
        print(current, num)
        if way == 'FORWARD':
            num += 1
        elif way == 'BACKWARD':
            num -= 1   
res = next('abcd')
print(res)
res = prev(res)
print(res)
res = num_on_str('bac')
print(res)
res = str_on_num('12345', 5)
print(res)
        

