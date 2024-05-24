import numpy as np
from numpy.linalg import det

class LUDecomposition:
    def __init__(self, dim, precision):
        if precision < 5:
            raise ValueError("eroare")
        if dim < 1:
            raise ValueError("eroare")

        self.dim = dim
        self.precision = precision
        self.threshold = 10 ** (-precision)
        self.matrix = np.array([[2, 0, 2], [1, 2, 5], [1, 1, 7]])
        self.vector = np.array([4, 10, 10])
        self.matrix_initial = self.matrix.copy()
        self.vector_initial = self.vector.copy()

    def crout_decomposition(self, matrix):
        n = len(matrix)
        L = np.zeros((n, n))
        U = np.eye(n)
        try:
            for k in range(n):
                for i in range(k, n):
                    L[i, k] = matrix[i, k] - np.sum(L[i, :k] * U[:k, k])
                for j in range(k + 1, n):
                    if np.abs(L[k, k]) > self.threshold:
                        U[k, j] = (matrix[k, j] - np.sum(L[k, :j] * U[:j, j])) / L[k, k]
                    else:
                        raise Exception(f"LU nu se poate descompune, det(A_{k + 1}) = 0")
        except Exception as exception:
            print(exception)
            exit()
        return L, U

    def crout_decomposition_restricted(self, matrix):
        n = len(matrix)
        try:
            for k in range(n):
                for i in range(k, n):
                    matrix[i, k] = self.matrix_initial[i, k] - np.sum(matrix[i, :k] * matrix[:k, k])
                for j in range(k + 1, n):
                    if np.abs(matrix[k, k]) > self.threshold:
                        matrix[k, j] = (self.matrix_initial[k, j] - np.sum(matrix[k, :k] * matrix[:k, j])) / matrix[k, k]
                    else:
                        raise Exception(f"LU nu se poate descompune, det(A_{k + 1}) = 0")
        except Exception as exception:
            print(exception)
            exit()
        return matrix

    def solve_system(self, L, U, b):
        y = np.zeros_like(b, dtype=float)
        for i in range(self.dim):
            y[i] = b[i] - np.sum(L[i, :i] * y[:i])
            if np.abs(L[i, i]) > self.threshold:
                y[i] /= L[i, i]
            else:
                print(f"Impartire la 0 la i={i}, L[i, i]={L[i, i]}")
                exit()

        x = np.zeros_like(y, dtype=float)
        for i in range(self.dim - 1, -1, -1):
            x[i] = y[i] - np.sum(U[i, i + 1:] * x[i + 1:])
            if np.abs(U[i, i]) > self.threshold:
                x[i] /= U[i, i]
            else:
                print(f"Impartire la 0 la i={i}, U[i, i]={U[i, i]}")
                exit()
        return x

    def index(self, i, j):
        if i >= j:
            return i * (i + 1) // 2 + j
        else:
            return j * (j + 1) // 2 + i

    def optimized_lu(self, matrix):
        n = len(matrix)
        L = np.ones(n * (n + 1) // 2)
        U = np.ones(n * (n + 1) // 2)
        try:
            for k in range(n):
                for i in range(k, n):
                    L[self.index(i, k)] = matrix[i, k] - sum(L[self.index(i, m)] * U[self.index(m, k)] for m in range(k))
                for j in range(k + 1, n):
                    if np.abs(L[self.index(k, k)]) > self.threshold:
                        U[self.index(k, j)] = (matrix[k, j] - sum((L[self.index(k, m)] * U[self.index(m, j)]) for m in range(k))) / L[self.index(k, k)]
                    else:
                        raise Exception(f"LU nu se poate descompune, det(A_{k + 1}) = 0")
        except Exception as exception:
            print(exception)
            exit()
        return L, U

    def solve_optimized(self, L, U, b):
        y = np.zeros_like(b, dtype=float)
        for i in range(self.dim):
            y[i] = b[i] - sum(L[self.index(i, k)] * y[k] for k in range(i))
            if np.abs(L[self.index(i, i)]) > self.threshold:
                y[i] /= L[self.index(i, i)]
            else:
                print(f"Impartire la 0 la i={i}, L[i, i]={L[i, i]}")
                exit()

        x = np.zeros_like(y, dtype=float)
        for i in range(self.dim - 1, -1, -1):
            x[i] = y[i] - sum(U[self.index(i, k)] * x[k] for k in range(i + 1, self.dim))
            if np.abs(U[self.index(i, i)]) > self.threshold:
                x[i] /= U[self.index(i, i)]
            else:
                print(f"Impartire la 0 la i={i}, U[i, i]={U[i, i]}")
                exit()
        return x

    def run(self):
        print(f"Parametrii programului sunt: dim = {self.dim}, precision = {self.precision}, threshold = {self.threshold}\n")

        L, U = self.crout_decomposition(self.matrix)
        print("LU folosind crout: ")
        print("Matrix = \n", self.matrix)
        print("L = \n", L)
        print("U = \n", U)

        A_restricted = self.crout_decomposition_restricted(self.matrix.copy())
        print("\nDescompunerea lui Matrix = \n", A_restricted)

        print("\ndet(Matrix) = det(matrix_initial)")
        det_A = 1
        for i in range(len(A_restricted)):
            det_A *= A_restricted[i, i]
        print(det_A, "=", det(self.matrix_initial))

        x_LU = self.solve_system(L, U, self.vector)
        print("\nSolutia sistemului: ")
        print(x_LU)

        print("\nVerificarea solutiei ")
        norm = np.linalg.norm(np.dot(self.matrix_initial, x_LU) - self.vector_initial)
        print(f"||Matrix*x_LU - vector||2 = {norm}")
        print(f"Norma < 1e-8:", norm < 10 ** (-8))

        print("\nNorme:")
        x_lib = np.linalg.solve(self.matrix_initial, self.vector)
        matrix_inv = np.linalg.inv(self.matrix_initial)
        norm1 = np.linalg.norm(x_LU - x_lib)
        norm2 = np.linalg.norm(x_LU - np.dot(matrix_inv, self.vector_initial))

        print(f"||x_LU - x_lib||2 = {norm1}")
        print(f"||x_LU - matrix^(-1)*vector_initial||2 = {norm2}")

        L_optimized, U_optimized = self.optimized_lu(self.matrix_initial)
        x_optimized = self.solve_optimized(L_optimized, U_optimized, self.vector)
        print("\nSolutia folosind 2 vectori L and U: \n", x_optimized)


if __name__ == "__main__":
    try:
        dim = int(input())
        precision = int(input())
        decomposition = LUDecomposition(dim, precision)
        decomposition.run()
    except ValueError as exception:
        print(exception)
