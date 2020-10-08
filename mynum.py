from numpy import *
from matplotlib.pyplot import *
from matplotlib.patches import Patch

# f in the IVP y’ = f(t,y), y(x0)=y0
def f(x, y):
    return (x - 3.2) * y + 8 * x * exp((x - 3.2) ** 2 / 2) * cos(4 * x ** 2)


def dfy(x, y):
    return x - 3.2


# analytic solution to the IVP y’ = f(t,y), y(x0)=y0
def sol(x, x0, y0):
    C = y0 * exp(-((x0 - 3.2) ** 2) / 2) - sin(4 * x0 ** 2)
    return exp((x - 3.2) ** 2 / 2) * (sin(4 * x ** 2) + C)


# Runge-Kutta "Classic" Order 4 method
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


# Adams-Bashforth 3 Step Method
def adams_3(x0, xn, n, y0):
    step = abs(xn - x0) / n
    x = linspace(x0, xn, n + 1)
    y = zeros(n + 1)
    y[0:3] = runge_kutta(x0, x0 + 2 * step, 2, y0)
    K1 = f(x[1], y[1])
    K2 = f(x[0], y[0])
    for i in range(2, n):
        K3 = K2
        K2 = K1
        K1 = f(x[i], y[i])
        y[i + 1] = y[i] + step * (23 * K1 - 16 * K2 + 5 * K3) / 12
    return y


# Adams-Bashforth 3/Moulton 4 Step Predictor/Corrector
def adams_precor(x0, xn, n, y0):
    step = abs(xn - x0) / n
    x = linspace(x0, xn, n + 1)
    y = zeros(n + 1)
    # Calculate initial steps with Runge-Kutta 4
    y[0:3] = runge_kutta(x0, x0 + 2 * step, 2, y0)
    K1 = f(x[1], y[1])
    K2 = f(x[0], y[0])
    for i in range(2, n):
        K3 = K2
        K2 = K1
        K1 = f(x[i], y[i])
        # Adams-Bashforth Predictor
        y[i + 1] = y[i] + step * (23 * K1 - 16 * K2 + 5 * K3) / 12
        K0 = f(x[i + 1], y[i + 1])
        # Adams-Moulton Corrector
        y[i + 1] = y[i] + step * (9 * K0 + 19 * K1 - 5 * K2 + K3) / 24
    return y

