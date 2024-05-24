import math
import random

def horner(coeffs, x):
    result = coeffs[0]
    for coeff in coeffs[1:]:
        result = result * x + coeff
    return result

def calculate_abc(coefficients, x):
    x_0, x_1, x_2 = x
    h_0 = x_1 - x_0
    h_1 = x_2 - x_1
    delta_0 = calculate_delta(coefficients, h_0, x_1, x_0)
    delta_1 = calculate_delta(coefficients, h_1, x_2, x_1)
    a = calculate_div((delta_1 - delta_0), (h_1 + h_0))
    b = a * h_1 + delta_1
    c = horner(coefficients, x_2)
    return {"x_0": x_0, "x_1": x_1, "x_2": x_2, "h_0": h_0, "h_1": h_1, "delta_0": delta_0, "delta_1": delta_1, "a": a,
            "b": b, "c": c}

def calculate_polynomial(coefficients, x):
    result = coefficients[0]
    for coeff in coefficients[1:]:
        result = result * x + coeff
    return result

def calculate_delta(coefficients, h, x, y):
    return calculate_div((calculate_polynomial(coefficients, x) - calculate_polynomial(coefficients, y)), h)

def calculate_div(a, b):
    if b == 0:
        return 0
    return a / b

def signum(x):
    if x >= 0:
        return 1
    else:
        return -1

def muller_method(f, x0, x1, x2, tol=1e-6, max_iter=100):
    iteration_count = 0
    while iteration_count < max_iter:
        iteration_count += 1
        result = calculate_abc(polynomial_coefficients, [x0, x1, x2])
        delta = result['b'] ** 2 - 4 * result['a'] * result['c']
        if delta < 0:
            return None
        delta_x = calculate_div(2 * result['c'], result['b'] + signum(result['b']) * math.sqrt(delta))
        x_3 = x2 - delta_x
        if abs(delta_x) < tol:
            return x_3
        x0, x1, x2 = x1, x2, x_3
    return None

def find_roots_muller(coeffs, start_points, tol=1e-6, max_iter=100):
    roots = set()
    for i in range(len(start_points)):
        x0, x1, x2 = start_points[i], start_points[(i+1)%len(start_points)], start_points[(i+2)%len(start_points)]
        root = muller_method(lambda x: horner(coeffs, x), x0, x1, x2, tol, max_iter)
        if root is not None:
            roots.add(round(root, 6))
    return roots

def calculate_roots_interval(coeffs):
    A = max(abs(coeff) for coeff in coeffs[1:])
    R = (abs(coeffs[0]) + A) / abs(coeffs[0])
    return -R, R

polynomial_coefficients = [1.0, -6.0, 11.0, -6.0]
roots_interval = calculate_roots_interval(polynomial_coefficients)
print(roots_interval)

roots = set()
for i in range(100):
    start_points = [random.uniform(roots_interval[0], roots_interval[1]) for _ in range(3)]
    roots.update(find_roots_muller(polynomial_coefficients, start_points))

with open("roots.txt", "w") as file:
    for root in roots:
        file.write(str(root) + "\n")
