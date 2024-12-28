import sympy
import numpy as np
from sympy import *
from Parameters import *
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.optimize import fsolve
import math


def Calculate_New_Capacitance(Number):
    R, L, C, r, m, S, n, x, y = symbols('R, L, C, r, m, S, n, x, y', real = True)

    a = (R**2*C)/(m*r)
    b = a + (np.pi**2*r**3*n)/(1-y)
    c = (2/3)*(np.pi**2*r**3*n) + 2*S
    d = (2*L)/(m*r)
    e = sqrt(b**2 - 2*a*(d + c))
    g = m*r*C

    omega_of_mu_1 = (sqrt(2))/(sqrt(-g*b + g*d + g*c + g*e)) #для первого промежутка на графике (до резонанса + резонанс)
    omega_of_mu_2 = (sqrt(2))/(sqrt(-g*b + g*d + g*c - g*e)) #для второго промежутка на графике (после резонанса)

    omega_subs = omega_of_mu_2.subs({R: Params['R'], L: Params['L'], C: Params['C'], r: Params['Radius'], m: Params['mu_0'], S: Params['Sigma'], n: 1/(Params['Dz'])**3})
    omega_func = lambdify(y, omega_subs, "math")

    mu_value_1 = 0
    value_omega = omega_func(mu_value_1)/2/np.pi/10**6

    print("Значение частоты при магнитной проницаемости", mu_value_1, "равно", value_omega, "МГц")

    h = (2*L)/(m*r) + (2/3)*np.pi**2*r**3*n + 2*S - (np.pi**2*r**3*n)/(1-y)
    s = x**2*m*r

    c_1 = 2/(s*(h + sqrt(((np.pi**2*r**3*n)/(1-y))**2 - ((2*R)/(x*m*r))**2))) #значения для промежутка частот до резонанса
    c_2 = 2/(s*(h - sqrt(((np.pi**2*r**3*n)/(1-y))**2 - ((2*R)/(x*m*r))**2))) #после резонанса

    C_subs = c_1.subs({R: Params['R'], L: Params['L'], r: Params['Radius'], m: Params['mu_0'], S: Params['Sigma'], n: 1/(Params['Dz'])**3, x: value_omega*2*np.pi*10**6})

    C_func = lambdify(y, C_subs, 'math')

    mu_value_2 = 0
    C_value = C_func(mu_value_2)


    Number = Params['N']['z']['ny']
    print('n=', Number)

    step_for_grad = abs(Params['C'] - C_value)/15

    print(step_for_grad)
    print(C_value)

    return step_for_grad, C_value
