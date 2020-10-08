from data import function_strings
from util import print_selection_of, read_int
from mynum import *

print_selection_of(function_strings)
function = read_int("Choose a function: ")

fg = 1
n = 300
x0 = 0
xn = 6
y0 = 0.75
x = linspace(x0, xn, n + 1)
ypc = adams_precor(x0, xn, n, y0)
figure(fg)
plot(x, ypc, "cyan", label="Predictor/Corrector 3/4")
x = linspace(x0, xn, 401)
ysol = sol(x, x0, y0)
plot(x, ysol, color="green", label="Exact")
title("n = %d" % n)
axis([0, xn, -60, 40])
legend(loc="lower left")
show()
