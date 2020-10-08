from numpy import sin

function_strings = ["y' = y + (1 + x) * y ^ 2", "y' = y * sin(x) + x", "y' = x^2 + y^2"]

functions = [
    lambda x, y: y + (1 + x) * y * y,
    lambda x, y: y * sin(x) + x,
    lambda x, y: x * x + y * y,
]
