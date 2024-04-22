import numpy as np
import matplotlib.pyplot as plt

def aitken(x, y, target_x):
    n = len(x)
    F = np.zeros((n, n))
    for i in range(n):
        F[i, 0] = y[i]
    for j in range(1, n):
        for i in range(n - j):
            F[i, j] = F[i + 1, j - 1] + (F[i + 1, j - 1] - F[i, j - 1]) / (x[i + j] - x[i])
    for k in range(n):
        if x[k] == target_x:
            return F[0, k]
    return None

def function(x):
    return x ** 4 - 12 * x ** 3 + 30 * x ** 2 + 12

def least_squares(x, y, target_x):
    n = len(x)
    for m in range(1, min(6, n)):
        A = np.vander(x, m + 1, increasing=True)
        coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        Pm_x = np.polyval(coeffs, target_x)
        return Pm_x

if __name__ == '__main__':
    n = int(input("Enter number of points: "))
    x0 = float(input("x0: "))
    xn = float(input("xn: "))
    h = (xn - x0) / n

    x = [x0 + i * h for i in range(n + 1)]
    y = [function(var) for var in x]
    print(x)
    print(y)
    new_x = float(input("Enter the new value of x: "))

    Ln_x = aitken(x, y, new_x)
    Pm_x = least_squares(x, y, new_x)

    print("Ln(x): ", Ln_x)
    print("Pm(x): ", Pm_x)
    if Ln_x is not None:
        print("Ln(x): ", Ln_x)
        print("|Ln(x) - f(x)|: ", abs(Ln_x - function(new_x)))

    plt.plot(x, y, label='f(x)', color='black', linestyle='dashed')
    plt.plot(new_x, Ln_x, 'bo', label='Ln(x)')
    plt.plot(new_x, Pm_x, 'ro', label='Pm(x)')

    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.show()
