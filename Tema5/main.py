import math as m
from cmath import sqrt

import numpy as np
from numpy import matrix
from scipy.linalg import lu  #descompunere Cholesky
from scipy.linalg import svd

# p = n
def Jacobi(A, epsilon):
    n = len(A)
    k = 0
    k_max = 1000
    U = np.identity(n)
    p, q = calcul_p_q(A)
    alfa = calcul_alfa(A, p, q)
    s, c, t = calcul_TSC(alfa)
    print("T este:")
    print(t)
    while abs(A[p][q]) > epsilon and k <= k_max:
        A = rotatie_A(A, p, q, c, s, t)
        if k == 0:
            print(A)
        U = calcul_U(U, p, q, c, s)
        p, q = calcul_p_q(A)
        alfa = calcul_alfa(A, p, q)
        s, c, t = calcul_TSC(alfa)
        k = k + 1
    print("Lambda: ")
    print(np.round(A, 6))
    return (A.diagonal(), U)


def calcul_p_q(A):
    maximum_value = 0
    p = -1
    q = -1
    for i in range(0, len(A)):
        for j in range(0, i):
            if abs(A[i][j]) > maximum_value:
                maximum_value = abs(A[i][j])
                p = i
                q = j

    return (p, q)


def calcul_alfa(A, p, q):
    return ((A[p][p] - A[q][q]) / (2 * A[p][q]))


def calcul_TSC(alfa):
    if alfa >= 0:
        t = -alfa + m.sqrt(alfa * alfa + 1)
    else:
        t = -alfa - m.sqrt(alfa * alfa + 1)
    c = 1 / m.sqrt((1 + t * t))
    s = t / m.sqrt((1 + t * t))

    return (s, c, t)

def rotatie_A(A, p, q, c, s, t):
    for j in range(0, len(A)):
        if j != p and j != q:
            A[p][j] = (c * A[p][j]) + (s * A[q][j])
            A[j][q] = (-s * A[j][p]) + (c * A[q][j])
            A[q][j] = A[j][q]
            A[j][p] = A[p][j]
    A[p][p] = A[p][p] + (t * A[p][q])
    A[q][q] = A[q][q] - (t * A[p][q])
    A[p][q] = 0
    A[q][p] = 0
    return A


def calcul_U(U, p, q, c, s):
    for i in range(0, len(U)):
        U_ip_vechi = U[i][p]
        U[i][p] = c * U[i][p] + s * U[i][q]
        U[i][q] = -s * U_ip_vechi + c * U[i][q]
    return U

def actualizare_norma(A, lambda_diagonal, U):
    norma = abs(A * U - U * lambda_diagonal)
    return norma

def definire_SI(S, p, n):
    SI = [[0 for i in range(p)] for j in range(n)]
    for i in range(p):
        for j in range(n):
            if i == j:
                SI[i][j] = 1 / S[i]
    return SI
def Jacobi_test(A, epsilon):
    U, lambda_diagonal = Jacobi(A, epsilon)
    norma = actualizare_norma(A, lambda_diagonal, U)
    print("Norma matriciala este:")
    print(np.linalg.norm(norma))



def diferenta_matriciala(A, epsilon):
    k = 0
    k_max = 100000
    P, L, U = lu(A)
    L = np.matrix(L)
    index = 1
    diferenta = 1
    while np.linalg.norm(diferenta) > epsilon and k < k_max:
        for i in range(index):
            L = L * L
        L_transpus = L.transpose()
        A_prim = L * L_transpus
        A_secund = L_transpus * L
        A_prima_matrice = np.matrix(A_prim)
        A_matricea_secunda = np.matrix(A_secund)
        diferenta = np.linalg.norm(abs(A_prima_matrice - A_matricea_secunda))
        index += index
    print("A(k): ")
    print(A_prima_matrice)
    print("Diferenta dintre A(k) si A(k+1):")
    print(diferenta)


def SVD(A, p, n):
    U, s, VT = svd(A)
    print("Valorile singulare ale matricei A:")
    print(s)

    print("Rangul matricei A este:")
    A_prim = matrix(A)
    print(np.linalg.matrix_rank(A_prim))

    print("Numarul de conditionare al matricei A:")
    min_value = min(s)
    max_value = max(s)
    print(max_value / min_value)

    print("Pseudoinversa Moore-Penrose a matricei A:")
    SI = definire_SI(s, p, n)
    SI = np.array(SI, dtype= float)
    UT = U.transpose()
    V = VT.transpose()
    UT = np.array(UT, dtype = float)
    V = np.array(V, dtype = float)
    AI = V.dot(SI).dot(UT)
    print(AI)

    print("Matricea pseudo-inversa in sensul celor mai mici patrate:")
    AT = A.transpose()
    AJ = np.linalg.inv((AT.dot(A))).dot(AT)
    print(AJ)
    norma = np.linalg.norm(abs(AI - AJ))
    print("Norma este:", float(norma))

A1 = np.array([
        [1, 1, 2],
        [1, 1, 2],
        [2, 2, 2]
    ], dtype = float)
A2 = np.array([
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [1 ,0, 1, 0],
        [0, 1, 0, 1]
    ], dtype = float)
A3 = np.array([
        [1, 2, 3, 4],
        [2, 3, 4, 5],
        [3, 4, 5, 6],
        [4, 5, 6, 7]
    ], dtype = float)

A4 = np.array([
        [1, 2, 3, 4],
        [2, 30, 4, 5],
        [3, 4, 5, 6],
        [4, 5, 6, 7]
    ], dtype = float)
epsilon = 0.0000001
print("Exemplu 1")
Jacobi_test(A1, epsilon)
print("Exemplu 2")
Jacobi_test(A2, epsilon)
print("Exemplu 3")
Jacobi_test(A3, epsilon)
print("Exemplu 4")
Jacobi_test(A4, epsilon)
diferenta_matriciala(A1, epsilon)
p = 3
n = 3
#SVD(A1, p, n)
p = 4
n = 4
SVD(A4, p ,n)
