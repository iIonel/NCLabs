import numpy as np
from numpy.linalg import det

try:
    n = int(input("Introduceti n: "))
    t = int(input("Introduceti t: "))
    if t < 5:
        raise ValueError("t trebuie sa fie mai mare decat 4")
    if n < 1:
        raise ValueError("n trebuie sa fie mai mare decat 0")
except ValueError as exception:
    print(exception)
    exit()

e = 10 ** (-t)

def crout(A):
    n = len(A)
    L = np.zeros((n, n)) # Matrice de n*n cu elementele egale cu 0
    U = np.eye(n) # Matrice de n*n cu 1 pe diagonala principala
    try:
        for p1 in range(n):
            for i in range(p1, n):
                L[i, p1] = A[i, p1] - np.sum(L[i, :p1] * U[:p1, p1])
            for i in range(p1 + 1, n):
                if np.abs(L[p1, p1]) > e:
                    U[p1, i] = (A[p1, i] - np.sum(L[p1, :i] * U[:i, i])) / L[p1, p1]
                else:
                    raise Exception(f"LU nu se poate descompune, det(A_{p1 + 1}) = 0")
    except Exception as exception:
        print(exception)
        exit()
    return L, U


def crout_restricted(A):
    n = len(A)
    try:
        for p in range(n):
            for i in range(p, n):
                A[i, p] = A_init[i, p] - np.sum(A[i, :p] * A[:p, p])
            for i in range(p + 1, n):
                if np.abs(A[p, p]) > e:
                    A[p, i] = (A_init[p, i] - np.sum(A[p, :p] * A[:p, i])) / A[p, p]
                else:
                    raise Exception(f"LU nu se poate descompune, det(A_{p + 1}) = 0")
    except Exception as exception:
        print(exception)
        exit()
    return A


print(f"Parametrii programului sunt: n = {n}, t = {t}, precession e={e}\n")

#A = np.random.rand(n, n)
#b = np.random.rand(n)
A = np.array([[2, 0, 2], [1, 2, 5], [1, 1, 7]])
b = np.array([4, 10, 10])

A_init = A.copy()
b_init = b.copy()

L, U = crout(A)
print("LU folosind crout: ")
print("A = \n", A)
print("L = \n", L)
print("U = \n", U)

A = crout_restricted(A)
print("\nDescompunerea lui A = \n", A)

print("\ndet(A) = det(A_init)")
det_A = 1
for i in range(len(A)):
    det_A = det_A * A[i, i]
print(det_A, "=", det(A_init))

# Substitutie directa, Ly = b
y = np.zeros_like(b, dtype=float) # Copie a lui b doar ca toate elem. sunt 0
for i in range(n):
    y[i] = b[i] - np.sum(A[i, :i] * y[:i])
    if np.abs(A[i, i]) > e:
        y[i] /= A[i, i]
    else:
        print(f"Impartire la 0 la i={i}, L[i, i]={L[i, i]}")
        exit()

#Substitutie inversa, Ux = y
x_LU = np.zeros_like(y, dtype=float)
for i in range(n - 1, -1, -1):
    x_LU[i] = y[i] - np.sum(A[i, i + 1:] * x_LU[i + 1:])
    if np.abs(A[i, i]) > e:
        x_LU[i] /= 1
    else:
        print(f"Impartire la 0 la i={i}, U[i, i]={U[i, i]}")
        exit()

print("\nSolutia sistemului: ")
print(x_LU)

print("\nVerificarea solutiei ")
norm = np.linalg.norm(np.dot(A_init, x_LU) - b_init) # Norma Euclidiana
print(f"||A*x_LU - b||2 = {norm}")
print(f"Norma < 1e-8:", norm < 10 ** (-8))

print("\nNorme:")

x_lib = np.linalg.solve(A_init, b)
A_inv = np.linalg.inv(A_init) # Inversa matricei

norm1 = np.linalg.norm(x_LU - x_lib)
norm2 = np.linalg.norm(x_LU - np.dot(A_inv, b_init))

print(f"||x_LU - x_lib||2 = {norm1}")
print(f"||x_LU - A^(-1)*b_init||2 = {norm2}")

#Bonus

def indexare(i, j):
    if i >= j:
        return i * (i + 1) // 2 + j
    else:
        return j * (j + 1) // 2 + i

def LU_optimized(A):
    n = len(A)
    L = np.ones(n * (n + 1) // 2)
    U = np.ones(n * (n + 1) // 2)
    try:
        for j in range(n):
            for i in range(j, n):
                L[indexare(i, j)] = A[i, j] - sum(L[indexare(i, k)] * U[indexare(k, j)] for k in range(j))
            for i in range(j + 1, n):
                if np.abs(L[indexare(j, j)]) > e:
                    U[indexare(j, i)] = (A[j, i] - sum((L[indexare(j, k)] * U[indexare(k, i)]) for k in range(j))) / L[
                        indexare(j, j)]
                else:
                    raise Exception(f"LU nu se poate descompune, det(A_{j + 1}) = 0")
    except Exception as exception:
        print(exception)
        exit()
    return L, U


L, U = LU_optimized(A_init)

#print("\nVectorul L:", L)
#print("Vectorul U:", U)

#Metoda substitutiei directe, Ly = b
y = np.zeros_like(b, dtype=float)
for i in range(n):
    y[i] = b[i] - sum(L[indexare(i, k)] * y[k] for k in range(i))
    if np.abs(L[indexare(i, i)]) > e:
        y[i] /= L[indexare(i, i)]
    else:
        print(f"Impartire la 0 la i={i}, L[i, i]={L[i, i]}")
        exit()

#Metoda substitutiei inverse, Ux = y
x_LU = np.zeros_like(y, dtype=float)
for i in range(n - 1, -1, -1):
    x_LU[i] = y[i] - sum(U[indexare(i, k)] * x_LU[k] for k in range(i + 1, n))
    if np.abs(U[indexare(i, i)]) > e:
        x_LU[i] /= U[indexare(i, i)]
    else:
        print(f"Impartire la 0 la i={i}, U[i, i]={U[i, i]}")
        exit()

print("\nSolutia folosind 2 vectori L and U: \n", x_LU)