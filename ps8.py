# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:02:29 2022

@author: Arjun
"""

import sympy as sp
x, y = sp.symbols('x, y')


def factorial(n):
    f = 1
    for i in range(2, n + 1):
        f *= i
    return f

def taylorSeriesMethod(fx, x0, y0, x1):
    arr = [y0]
    for i in range(5):
        arr.append(sp.diff(fx, x, i))
    print(arr)
    
    
def eulersMethod(fx, x0, y0, x1):
    n = 10
    h = (x1 - x0) / n
    for i in range(n):
        sub = fx.subs([(x, x0), (y, y0)])
        x0 = x0 + h
        y0 = y0 + h * sub
    return y0


def modifiedEulersMethod(fx, x0, y0, x1):
    h = 0.1

    while x1 - x0 > 0.00001:
        y0 += h * fx.subs([(x, x0 + 0.5 * h), (y, y0 + 0.5 * h * fx.subs([(x, x0), (y, y0)]))])
        x0 += h
    return y0


def rungeKuttaMethod(fx, x0, y0, x1):
    h = 0.1
    while x1 - x0 > 0.00001:
        k1 = h * fx.subs([(x, x0), (y, y0)])
        k2 = h * fx.subs([(x, x0 + 0.5 * h), (y, y0 + 0.5 * k1)])
        k3 = h * fx.subs([(x, x0 + 0.5 * h), (y, y0 + 0.5 * k2)])
        k4 = h * fx.subs([(x, x0 + h), (y, y0 + k3)])
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        y0 += k
        x0 += h
    return sp.N(y0)


def milnesPredictorCorrectorMethod(fx, X, Y, x0):
    h = X[1] - X[0]
    n = len(X) - 1
    
    pValue = Y[n - 3] + 4 * h * (2 * fx.subs([(x, X[n - 2]), (y, Y[n - 2])]) - fx.subs([(x, X[n - 1]), (y, Y[n - 1])]) + 2 * fx.subs([(x, X[n]), (y, Y[n])])) / 3
    cValue = Y[n - 1] + h * (fx.subs([(x, X[n - 1]), (y, Y[n - 1])]) + fx.subs([(x, X[n] + h), (y, pValue)]) + 4 * fx.subs([(x, X[n]), (y, Y[n])])) / 3
    print(f'y({x0}) = {cValue}')
    
    
    
def adamBashforthMethod(fx, X, Y, x0):
    h = X[1] - X[0]
    n = len(X) - 1

    pValue = Y[n] + h * (55 * fx.subs([(x, X[n]), (y, Y[n])]) - 59 * fx.subs([(x, X[n - 1]), (y, Y[n - 1])]) + 37 * fx.subs([(x, X[n - 2]), (y, Y[n - 2])]) - 9 * fx.subs([(x, X[n - 3]), (y, Y[n - 3])])) / 24
    cValue = Y[n] + h * (9 * fx.subs([(x, X[n] + h), (y, pValue)]) + 19 * fx.subs([(x, X[n]), (y, Y[n])]) - 5 * fx.subs([(x, X[n - 1]), (y, Y[n - 1])]) + fx.subs([(x, X[n - 2]), (y, Y[n - 2])])) / 24
    print(f'y({x0}) = {cValue}')
    
    
X = [0.1, 0.2, 0.3]
Y = [eulersMethod(2 * x * y / (1 + x ** 2), 0, 0, i) for i in X]
Y = [taylorSeriesMethod(x * y + x ** 2, 0, 1, i) for i in X]
Y = [rungeKuttaMethod(2 * x * y / (1 + x ** 2), 0, 0, i) for i in X]   
Y = [modifiedEulersMethod((2 - y ** 2) / (5 * x), 4, 1, i) for i in X]

milnesPredictorCorrectorMethod((1 + x ** 2) * y ** 2,X,Y, 0.4)
adamBashforthMethod((1 + x ** 2) * y ** 2, [0, 0.1, 0.2, 0.3], [1, 1.06, 1.12, 1.21], 0.4)