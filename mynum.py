from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Patch


def set_func():
    pass


# f в y’ = f(x,y), y(x0)=y0
def f(x, y):
    return (x - 3.2) * y + 8 * x * exp((x - 3.2) ** 2 / 2) * cos(4 * x ** 2)


# analytic solution to the IVP y’ = f(t,y), y(x0)=y0
def sol(x, x0, y0):
    C = y0 * exp(-((x0 - 3.2) ** 2) / 2) - sin(4 * x0 ** 2)
    return exp((x - 3.2) ** 2 / 2) * (sin(4 * x ** 2) + C)


# Метод Рунге-Кутты
def runge_kutta(x0, xn, n, y0):
    h = abs(xn - x0) / n
    x = linspace(x0, xn, n + 1)
    y = zeros(n + 1)
    y[0] = y0
    for i in range(0, n):
        K1 = f(x[i], y[i])
        K2 = f(x[i] + h / 2, y[i] + K1 * h / 2)
        K3 = f(x[i] + h / 2, y[i] + K2 * h / 2)
        K4 = f(x[i] + h, y[i] + K3 * h)
        y[i + 1] = y[i] + h * (K1 + 2 * K2 + 2 * K3 + K4) / 6
    return y


# Метод Адамса с предиктором и корректором
def adams_precor(x0, xn, n, y0):
    h = abs(xn - x0) / n
    x = linspace(x0, xn, n + 1)
    y = zeros(n + 1)
    # Вычисляем начальные шаги с помощью метода РК
    y[0:3] = runge_kutta(x0, x0 + 2 * h, 2, y0)
    K1 = f(x[1], y[1])
    K2 = f(x[0], y[0])
    for i in range(2, n):
        K3 = K2
        K2 = K1
        K1 = f(x[i], y[i])
        # Предиктор
        y[i + 1] = y[i] + h * (23 * K1 - 16 * K2 + 5 * K3) / 12
        K0 = f(x[i + 1], y[i + 1])
        # Корректор
        y[i + 1] = y[i] + h * (9 * K0 + 19 * K1 - 5 * K2 + K3) / 24
    return y

