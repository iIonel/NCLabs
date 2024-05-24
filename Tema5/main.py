import math as m
from cmath import sqrt

import numpy as np
from numpy import matrix
from scipy.linalg import lu
from scipy.linalg import svd

class JacobiMethod:
    def __init__(self, matrix, epsilon):
        self.matrix = matrix
        self.epsilon = epsilon
        self.size = len(matrix)
        self.max_iterations = 1000

    def run(self):
        iterations = 0
        eigenvectors = np.identity(self.size)
        p, q = self._find_max_off_diagonal()
        alfa = self._calculate_alpha(p, q)
        s, c, t = self._calculate_TSC(alfa)
        print(t)
        while abs(self.matrix[p][q]) > self.epsilon and iterations <= self.max_iterations:
            self._rotate(p, q, c, s, t)
            if iterations == 0:
                print(self.matrix)
            eigenvectors = self._update_eigenvectors(eigenvectors, p, q, c, s)
            p, q = self._find_max_off_diagonal()
            alfa = self._calculate_alpha(p, q)
            s, c, t = self._calculate_TSC(alfa)
            iterations += 1
        print(np.round(self.matrix, 6))
        return (self.matrix.diagonal(), eigenvectors)

    def _find_max_off_diagonal(self):
        max_value = 0
        p = -1
        q = -1
        for i in range(self.size):
            for j in range(i):
                if abs(self.matrix[i][j]) > max_value:
                    max_value = abs(self.matrix[i][j])
                    p = i
                    q = j
        return (p, q)

    def _calculate_alpha(self, p, q):
        return (self.matrix[p][p] - self.matrix[q][q]) / (2 * self.matrix[p][q])

    def _calculate_TSC(self, alfa):
        if alfa >= 0:
            t = -alfa + m.sqrt(alfa * alfa + 1)
        else:
            t = -alfa - m.sqrt(alfa * alfa + 1)
        c = 1 / m.sqrt(1 + t * t)
        s = t / m.sqrt(1 + t * t)
        return (s, c, t)

    def _rotate(self, p, q, c, s, t):
        for j in range(self.size):
            if j != p and j != q:
                self.matrix[p][j] = (c * self.matrix[p][j]) + (s * self.matrix[q][j])
                self.matrix[j][q] = (-s * self.matrix[j][p]) + (c * self.matrix[q][j])
                self.matrix[q][j] = self.matrix[j][q]
                self.matrix[j][p] = self.matrix[p][j]
        self.matrix[p][p] += t * self.matrix[p][q]
        self.matrix[q][q] -= t * self.matrix[p][q]
        self.matrix[p][q] = 0
        self.matrix[q][p] = 0

    def _update_eigenvectors(self, eigenvectors, p, q, c, s):
        for i in range(self.size):
            u_old = eigenvectors[i][p]
            eigenvectors[i][p] = c * eigenvectors[i][p] + s * eigenvectors[i][q]
            eigenvectors[i][q] = -s * u_old + c * eigenvectors[i][q]
        return eigenvectors


class MatrixDifference:
    def __init__(self, matrix, epsilon):
        self.matrix = matrix
        self.epsilon = epsilon
        self.max_iterations = 100000

    def compute_difference(self):
        P, L, U = lu(self.matrix)
        L = np.matrix(L)
        index = 1
        difference = 1
        iterations = 0
        while np.linalg.norm(difference) > self.epsilon and iterations < self.max_iterations:
            for _ in range(index):
                L = L * L
            L_transpose = L.transpose()
            A_prime = L * L_transpose
            A_secondary = L_transpose * L
            difference = np.linalg.norm(abs(A_prime - A_secondary))
            index += index
            iterations += 1
        print("A(k): ")
        print(A_prime)
        print(difference)


class SVDDecomposition:
    def __init__(self, matrix):
        self.matrix = matrix

    def run(self):
        U, s, VT = svd(self.matrix)
        print(s)

        A_prime = matrix(self.matrix)
        print(np.linalg.matrix_rank(A_prime))

        min_value = min(s)
        max_value = max(s)
        print(max_value / min_value)

        SI = self._define_SI(s, len(s), len(self.matrix[0]))
        UT = U.transpose()
        V = VT.transpose()
        AI = V.dot(SI).dot(UT)
        print(AI)

        AT = self.matrix.transpose()
        AJ = np.linalg.inv(AT.dot(self.matrix)).dot(AT)
        print(AJ)
        norm = np.linalg.norm(abs(AI - AJ))
        print("Norm is:", float(norm))

    def _define_SI(self, S, p, n):
        SI = np.zeros((p, n))
        for i in range(p):
            if i < len(S):
                SI[i][i] = 1 / S[i]
        return SI


def main():
    A1 = np.array([
        [1, 1, 2],
        [1, 1, 2],
        [2, 2, 2]
    ], dtype=float)
    A2 = np.array([
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [1, 0, 1, 0],
        [0, 1, 0, 1]
    ], dtype=float)
    A3 = np.array([
        [1, 2, 3, 4],
        [2, 3, 4, 5],
        [3, 4, 5, 6],
        [4, 5, 6, 7]
    ], dtype=float)
    A4 = np.array([
        [1, 2, 3, 4],
        [2, 30, 4, 5],
        [3, 4, 5, 6],
        [4, 5, 6, 7]
    ], dtype=float)
    epsilon = 1e-7

    jacobi1 = JacobiMethod(A1, epsilon)
    jacobi1.run()
    jacobi2 = JacobiMethod(A2, epsilon)
    jacobi2.run()
    jacobi3 = JacobiMethod(A3, epsilon)
    jacobi3.run()
    jacobi4 = JacobiMethod(A4, epsilon)
    jacobi4.run()

    matrix_diff = MatrixDifference(A1, epsilon)
    matrix_diff.compute_difference()

    svd_decomp = SVDDecomposition(A4)
    svd_decomp.run()

if __name__ == "__main__":
    main()
