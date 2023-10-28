from matrix import Matrix


# Lab 3
test1 = (3,2,[2, 1,
              1, 0,
              1, 1,])

test2 = (2, 3, [3, 2, 6,
                0, 2, 2])

test3 = (3, 3, [1, 2, 3, 
                4, 5, 6,
                5, 6, 6])


m1 = Matrix(*test1)

B, C = m1.get_skeleton_decomposition()

print(B, end='\n\n')
print(C, end='\n\n')
print(B * C, end='\n\n')


m1 = Matrix(*test2)

B, C = m1.get_skeleton_decomposition()

print(B, end='\n\n')
print(C, end='\n\n')
print(B * C, end='\n\n')

m1 = Matrix(*test3)

B, C = m1.get_skeleton_decomposition()

print(B, end='\n\n')
print(C, end='\n\n')
print(B * C, end='\n\n')

def dot_product(v1, v2):
    # Вычисление скалярного произведения двух векторов
    return sum(x * y for x, y in zip(v1, v2))

def vector_subtraction(v1, v2):
    # Вычитание векторов
    return [x - y for x, y in zip(v1, v2)]

def qr_algorithm(matrix, max_iterations):
    size = len(matrix)
    eigenvalues = [0] * size  # Массив для хранения собственных значений
    eigenvectors = [[0] * size for _ in range(size)]  # Массив для хранения собственных векторов

    for _ in range(max_iterations):
        # QR-разложение матрицы
        Q, R = qr_decomposition(matrix)
        matrix = matrix_multiplication(R, Q)  # Умножение матрицы R на Q

    # Вычисление собственных значений и векторов
    for i in range(size):
        eigenvalues[i] = matrix[i][i]
        eigenvector = [0] * size
        eigenvector[i] = 1

        # Умножение вектора на Q
        for _ in range(max_iterations):
            eigenvector = matrix_vector_multiplication(Q, eigenvector)
        eigenvectors[i] = eigenvector

    return eigenvalues, eigenvectors

def qr_decomposition(matrix):
    size = len(matrix)
    Q = [[0] * size for _ in range(size)]
    R = [[0] * size for _ in range(size)]

    # Вычисление QR-разложения
    for j in range(size):
        v = matrix[j]
        for i in range(j):
            R[i][j] = dot_product(Q[i], matrix[j])
            v = vector_subtraction(v, scalar_vector_multiplication(R[i][j], Q[i]))
        R[j][j] = magnitude(v)
        Q[j] = normalize(v)

    return Q, R

def matrix_multiplication(matrix1, matrix2):
    result = [[0] * len(matrix2[0]) for _ in range(len(matrix1))]

    # Умножение матриц
    for i in range(len(matrix1)):
        for j in range(len(matrix2[0])):
            for k in range(len(matrix2)):
                result[i][j] += matrix1[i][k] * matrix2[k][j]

    return result

def matrix_vector_multiplication(matrix, vector):
    result = [0] * len(vector)

    # Умножение матрицы на вектор
    for i in range(len(matrix)):
        for j in range(len(vector)):
            result[i] += matrix[i][j] * vector[j]

    return result

def scalar_vector_multiplication(scalar, vector):
    # Умножение скаляра на вектор
    return [scalar * x for x in vector]

def magnitude(vector):
    # Вычисление длины вектора
    return sum(x**2 for x in vector)**0.5

def normalize(vector):
    # Нормализация вектора
    length = magnitude(vector)
    return [x / length for x in vector]



# Создание симметричной матрицы
matrix = [[5, -2],
          [-2, 8],]

# Запуск QR-алгоритма
eigenvalues, eigenvectors = qr_algorithm(matrix, max_iterations=100)

# Вывод собственных значений и векторов
print("Собственные значения:", eigenvalues)
print("Собственные векторы:")
for i, eigenvector in enumerate(eigenvectors):
    print("Собственное значение:", eigenvalues[i])
    print("Собственный вектор:", eigenvector)
    print()