import math
import random
import numpy as np

def generate_matrix(n):
    random_matrix = []
    for i in range(0,n):
        random_line = []
        for j in range(0,n):
            random_line.append(random.randint(0,100))
        random_matrix.append(random_line)
    return random_matrix

def generate_vector(n):
    random_vector = []
    for i in range(0,n):
        random_vector.append(random.randint(0,100))
    return random_vector

def calculate_b(A,s):
    b = [0 for i in range(0,len(A))]

    for i in range(0,len(A)):
        b[i] = sum([s[j] * A[i][j] for j in range(0, len(A))])
    return b

def calculare_QR(A, b):
    dimensiune = len(A)
    Q_Tilda = [[0 for j in range(dimensiune)] for i in range(dimensiune)]
    for i in range(dimensiune):
        Q_Tilda[i][i] = 1
    for r in range(dimensiune - 1):
        suma_pătrate = sum([A[i][r] ** 2 for i in range(r, dimensiune)])
        if suma_pătrate <= epsilon:
            return None
        k = math.sqrt(suma_pătrate)
        if A[r][r] > 0:
            k = -k
        beta = suma_pătrate - k * A[r][r]
        u = [0] * dimensiune
        u[r] = A[r][r] - k
        for i in range(r + 1, dimensiune):
            u[i] = A[i][r]
        for j in range(r + 1, dimensiune):
            gamma = sum([u[i] * A[i][j] for i in range(r, dimensiune)]) / beta
            for i in range(r, dimensiune):
                A[i][j] -= gamma * u[i]
        A[r][r] = k
        for i in range(r + 1, dimensiune):
            A[i][r] = 0
        gamma = sum([u[i] * b[i] for i in range(r, dimensiune)]) / beta
        for i in range(r, dimensiune):
            b[i] -= gamma * u[i]
        for j in range(dimensiune):
            gamma = sum([u[i] * Q_Tilda[i][j] for i in range(r, dimensiune)]) / beta
            for i in range(r, dimensiune):
                Q_Tilda[i][j] -= gamma * u[i]
    R = A
    Q_T = Q_Tilda
    Q_T_b = b
    return R, Q_T, Q_T_b

def calculate_x(A, b):
    x = [0 for i in range(0, len(A))]
    for i in range(len(A) - 1, -1, -1):
        s = 0
        for j in range(i, len(A)):
            s += A[i][j] * x[j]
        x[i] = (b[i] - s) / A[i][i]
    return x


n = int(input("Size of matrix: "))
A = np.array([[0,0,4],[1,2,3],[0,1,2]])
s = np.array([3,2,1])
epsilon = 0.00000000000000001

#ex1
b = calculate_b(A,s)

#ex2
R,Qt,Qtb = calculare_QR(A,b)

#ex3
A_matrix = np.array(A)
b_vector = np.array(b)
QQR, RQR = np.linalg.qr(A_matrix)
scalar = np.dot(QQR.T, b_vector)
xQR = np.linalg.solve(RQR, scalar)
xQR = xQR.tolist()
xHouseholder = calculate_x(R, Qtb)
norma = np.linalg.norm(np.array(xQR) - np.array(xHouseholder))

#ex4
norma1 = np.linalg.norm(np.dot(A, xHouseholder) - b)
norma2 = np.linalg.norm(np.dot(A,xQR) - b)
norma3 = np.linalg.norm(np.array(xHouseholder) - np.array(s)) / np.linalg.norm(s)
norma4 = np.linalg.norm(np.array(xQR) - np.array(s)) / np.linalg.norm(s)

print(A)
print(s)
print(b)
print(R)
print("Q: ",Qt)
print(Qtb)
print(norma)
print(norma1)
print(norma2)
print(norma3)
print(norma4)