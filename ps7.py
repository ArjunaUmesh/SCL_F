# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 07:59:10 2022

@author: Arjun
"""

import sympy as sp
x, y = sp.symbols('x, y')


def taylorSeriesMethod(fx, x0, y0, x1):
    arr = [y0]
    for i in range(5):
        arr.append(sp.diff(fx, x, i))
    print(arr)
    
    
taylorSeriesMethod(sp.exp(x ** 2 - 1) / (1 + x * y + x ** 3), 0, -3, 1.5)
    
    
def eulersMethod(fx, x0, y0, x1):
    n = 10
    h = (x1 - x0) / n
    for i in range(n):
        sub = fx.subs([(x, x0), (y, y0)])
        x0 = x0 + h
        y0 = y0 + h * sub
    print(f'y({x1}) = {y0}')
    
eulersMethod((x ** 2 - y ** 2) / (x ** 2 + y ** 2), 1, -2, 4)


def modifiedEulersMethod(fx, x0, y0, x1):
    h = 0.1

    while x1 - x0 > 0.00001:
        y0 += h * fx.subs([(x, x0 + 0.5 * h), (y, y0 + 0.5 * h * fx.subs([(x, x0), (y, y0)]))])
        x0 += h
    print(f'y({x1}) = {y0}')
    
    
modifiedEulersMethod((x + y ** 2 - 2 * x * y) / (x ** 2 - y ** 2 + 2 * x), -1, 0, 2)



def rungeKutta3Method(fx, x0, y0, x1):
    h = 0.1
    while x1 - x0 > 0.00001:
        k1 = h * fx.subs([(x, x0), (y, y0)])
        k2 = h * fx.subs([(x, x0 + 0.5 * h), (y, y0 + 0.5 * k1)])
        kd = h * fx.subs([(x, x0 + h), (y, y0 + k1)])
        k3 = h * fx.subs([(x, x0 + h), (y, y0 + kd)])
        k = (k1 + 4 * k2 + k3) / 6
        y0 += k
        x0 += h
    print(f'y({x1}) = {sp.N(y0)}')
        

rungeKutta3Method(sp.tan(x) / sp.exp(x + y), -0.5, -1.5, 0.75)



def rungeKutta4Method(fx, x0, y0, x1):
    h = 0.1
    while x1 - x0 > 0.00001:
        k1 = h * fx.subs([(x, x0), (y, y0)])
        k2 = h * fx.subs([(x, x0 + 0.5 * h), (y, y0 + 0.5 * k1)])
        k3 = h * fx.subs([(x, x0 + 0.5 * h), (y, y0 + 0.5 * k2)])
        k4 = h * fx.subs([(x, x0 + h), (y, y0 + k3)])
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        y0 += k
        x0 += h
    print(f'y({x1}) = {sp.N(y0)}')   
    
    

rungeKutta4Method(sp.sin(x) + sp.exp(x ** 2 - 2 * y), 2, 10, 3)




