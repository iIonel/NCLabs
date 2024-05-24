import math
import random
import numpy as np

class MatrixOperations:
    def __init__(self, size, epsilon):
        self.size = size
        self.epsilon = epsilon

    def generate_random_matrix(self):
        matrix = []
        for i in range(self.size):
            row = [random.randint(0, 100) for _ in range(self.size)]
            matrix.append(row)
        return matrix

    def generate_random_vector(self):
        return [random.randint(0, 100) for _ in range(self.size)]

    def calculate_b_vector(self, matrix, vector):
        b_vector = [sum(vector[j] * matrix[i][j] for j in range(self.size)) for i in range(self.size)]
        return b_vector

    def qr_decomposition(self, matrix, b):
        dim = len(matrix)
        Q_tilda = [[1 if i == j else 0 for j in range(dim)] for i in range(dim)]
        for r in range(dim - 1):
            sum_squares = sum(matrix[i][r] ** 2 for i in range(r, dim))
            if sum_squares <= self.epsilon:
                return None
            k = -math.sqrt(sum_squares) if matrix[r][r] > 0 else math.sqrt(sum_squares)
            beta = sum_squares - k * matrix[r][r]
            u = [0] * dim
            u[r] = matrix[r][r] - k
            for i in range(r + 1, dim):
                u[i] = matrix[i][r]
            for j in range(r + 1, dim):
                gamma = sum(u[i] * matrix[i][j] for i in range(r, dim)) / beta
                for i in range(r, dim):
                    matrix[i][j] -= gamma * u[i]
            matrix[r][r] = k
            for i in range(r + 1, dim):
                matrix[i][r] = 0
            gamma = sum(u[i] * b[i] for i in range(r, dim)) / beta
            for i in range(r, dim):
                b[i] -= gamma * u[i]
            for j in range(dim):
                gamma = sum(u[i] * Q_tilda[i][j] for i in range(r, dim)) / beta
                for i in range(r, dim):
                    Q_tilda[i][j] -= gamma * u[i]
        return matrix, Q_tilda, b

    def back_substitution(self, R, Qt_b):
        x = [0] * len(R)
        for i in range(len(R) - 1, -1, -1):
            s = sum(R[i][j] * x[j] for j in range(i, len(R)))
            x[i] = (Qt_b[i] - s) / R[i][i]
        return x

class QRComparison:
    def __init__(self, matrix, b_vector):
        self.matrix = matrix
        self.b_vector = b_vector

    def compute_with_numpy(self):
        A_matrix = np.array(self.matrix)
        b_vec = np.array(self.b_vector)
        Q, R = np.linalg.qr(A_matrix)
        Q_transpose_b = np.dot(Q.T, b_vec)
        x_numpy = np.linalg.solve(R, Q_transpose_b)
        return x_numpy.tolist()

    def compare_solutions(self, x_qr_custom, x_qr_numpy, s):
        norm_diff_solutions = np.linalg.norm(np.array(x_qr_numpy) - np.array(x_qr_custom))
        norm_residual_custom = np.linalg.norm(np.dot(self.matrix, x_qr_custom) - self.b_vector)
        norm_residual_numpy = np.linalg.norm(np.dot(self.matrix, x_qr_numpy) - self.b_vector)
        relative_error_custom = np.linalg.norm(np.array(x_qr_custom) - np.array(s)) / np.linalg.norm(s)
        relative_error_numpy = np.linalg.norm(np.array(x_qr_numpy) - np.array(s)) / np.linalg.norm(s)

        return norm_diff_solutions, norm_residual_custom, norm_residual_numpy, relative_error_custom, relative_error_numpy

if __name__ == "__main__":
    n = int(input("Size of matrix: "))
    epsilon = 1e-16

    matrix_ops = MatrixOperations(n, epsilon)
    matrix = np.array([[0,0,4],[1,2,3],[0,1,2]])
    s = np.array([3, 2, 1])

    b = matrix_ops.calculate_b_vector(matrix, s)

    R, Q_t, Q_t_b = matrix_ops.qr_decomposition(matrix.copy(), b.copy())

    qr_comp = QRComparison(matrix, b)
    x_numpy = qr_comp.compute_with_numpy()
    x_custom = matrix_ops.back_substitution(R, Q_t_b)
    norm_diff, norm_res_custom, norm_res_numpy, rel_err_custom, rel_err_numpy = qr_comp.compare_solutions(x_custom, x_numpy, s)


    print(matrix)
    print(s)
    print(b)
    print( R)
    print(Q_t)
    print(Q_t_b)
    print(norm_diff)
    print(norm_res_custom)
    print(norm_res_numpy)
    print(rel_err_custom)
    print(rel_err_numpy)
