from fractions import Fraction
import numpy as np
import time
import random


def benchmark(func):
    def wrapper(*args, **kwargs):
        start = round(time.time() * 1000000)
        return_value = func(*args, **kwargs)
        end = round(time.time() * 1000000)
        print('[*] Время выполнения: {} микросекунд:'.format(end-start))
        return return_value
    return wrapper

class ElementaryTransformationsMixin:
    def swap_rows(self, first:int, second:int) -> None:
        self._elements[first], self._elements[second] = self._elements[second], self._elements[first]
    
    def swap_columns(self, first:int, second:int) -> None:
        for row in self._elements:
            row[first], row[second] = row[second], row[first]
    
    def multiply_by_row(self, row:int, target) -> None:
        for i in range(self.column):
            self._elements[row][i] *= target
    
    def multiply_by_column(self, column:int, target) -> None:
        for i in range(self.row):
            self._elements[i][column] *= target
            
    def divide_by_num(self, target):
        #target = Fraction(target)
        if target == 0:
            raise Exception('divison by zero')
        for i in range(self.row):
            for j in range(self.column):
                self[i, j] = self[i, j] / target
        
    def mul_by_num(self, target):
        #target = Fraction(target)
        for i in range(self.row):
            for j in range(self.column):
                self._elements[i][j] *= target
    
    def sum_rows(self, source_row:int, target_row:int) -> None:
        for i in range(self.column):
            self._elements[target_row][i] += self._elements[source_row][i]
        
    def sum_columns(self, source_column:int, target_column:int) -> None:
        for i in range(self.row):
            self._elements[i][target_column] += self._elements[i][source_column]
    
    def sum_num(self, num):
        for i in range(self.row):
            for j in range(self.column):
                self._elements[i][j] += num
    
    def sub_num(self, num):
        for i in range(self.row):
            for j in range(self.column):
                self._elements[i][j] -= num


# Метод Гаусса
class GaussMethodMixin:
    def to_fractions(self) -> None:
        for i in range(self.row):
            for j in range(self.column):
                self._elements[i][j] = Fraction(self._elements[i][j])
                        
    def _set_max_el_in_row(self, i:int, count_swap:int) -> int:
        max_el = abs(self._elements[i][i])
        max_row = i
        
        for j in range(i+1, self.row):
            if max_el < abs(self._elements[j][i]):
                max_el = abs(self._elements[j][i])
                max_row = j
        
        if i != max_row:
            self.swap_rows(i, max_row)
            count_swap += 1
        
        return count_swap
    
    def _set_max_el_in_column(self, i:int, count_swap:int) -> int:
        max_el = abs(self._elements[i][i])
        max_column = i
        
        for j in range(i+1, self.column):
            if max_el < abs(self._elements[i][j]):
                max_el = abs(self._elements[i][j])
                max_column = j
        
        if i != max_column:
            self.swap_columns(i, max_column)
            count_swap += 1
        
        return count_swap
    
    def get_triangle(self):
        matr = self.copy()
        count_swap = 0
        
        for i in range(matr.row):
            if i == matr.column:
                continue
            if matr._elements[i][i] == 0:
                count_swap += matr._set_max_el_in_row(i, count_swap)
            
            if matr._elements[i][i] == 0:
                count_swap += matr._set_max_el_in_column(i, count_swap)

            if matr._elements[i][i] == 0:
                return matr, None
            
            for j in range(i+1, matr.row):
                c = -(matr._elements[j][i] / matr._elements[i][i])
                for k in range(i, matr.column):
                    matr._elements[j][k] += c * matr._elements[i][k]
            
        return matr, count_swap


class GrevilleMethod:
    def is_zeros_matr(self):
        for i in range(self.row):
            for j in range(self.column):
                if self._elements[i][j] != 0:
                    return False
        return True
    
    def get_pseudoinverse_matrix(self):
        list_A = []
        k = 0
        while k < self.column:
            a_k = self.get_column(k, 1)
            k += 1
            if k == 1:
                a_kT = a_k.T
                list_A.append(((a_kT * a_k) ** -1) * a_kT)
            else:
                d_k = list_A[-1] * a_k
                c_k = a_k - (self.get_column(range(k-1), k-1) * d_k)
                if c_k.is_zeros_matr():
                    temp_d = d_k.T * d_k
                    temp_d.sum_num(1)
                    b_k = temp_d ** -1 * d_k.T * list_A[-1]
                else:
                    b_k = (c_k.T * c_k) ** -1 * c_k.T
                B = list_A[-1] - d_k * b_k
                B.add_row(b_k)
                list_A.append(B)
        return list_A[-1]


class EigenvalueEigenvectorMixin:
    def power_method(self, epsilon=1e-6, max_iter=100):
        """
        return: Largest Eigenvalue and Eigenvector
        """
        tol = epsilon
        max_iter = max_iter
        x = Matrix(self.row,  1, [i + 1 for i in range(self.row)])
        lam_prev = 0
        
        for i in range(max_iter): 
            norm = (self * x).norm()
            temp_x = self * x
            temp_x.divide_by_num(norm)
            x = temp_x
        
            lam = (x.T * self * x)._elements[0][0] / (x.T * x)._elements[0][0]
        
            if np.abs(lam - lam_prev) < tol: 
                break
        
            lam_prev = lam
        
        return float(lam), x 
    
    def norm(self):
        norm = 0
        for i in range(self.row):
            for j in range(self.column):
                norm += self._elements[i][j] ** 2
        return norm ** 0.5

    
    def dot_sum(self, right):
        if self.row == right.row and self.column == right.column:
            temp = []
            for i in range(self.row):
                for j in range(self.column):
                    temp.append(self._elements[i][j] * right._elements[i][j])
            return sum(temp)
        return None
    
    def set_column(self, vec, column):
        for i in range(self.row):
            self._elements[i][column] = vec._elements[i][0]

    def qr_decomposition(self):
        n = self.column
        
        Q = Matrix(self.row, self.column, [0 for _ in range(self.row*self.column)])
        R = Matrix(self.row, self.column, [0 for _ in range(self.row*self.column)])

        for i in range(n):
            # Begin the Gram-Schmidt process
            v = self.get_column(i, 1)
            
            for j in range(i):
                q_col = Q.get_column(j, 1)
                temp = q_col.dot_sum(self.get_column(i, 1))
                R._elements[j][i] = temp
                q_col.mul_by_num(temp)
                v = v - q_col
            
            R._elements[i][i] = v.norm()
            v.divide_by_num(R._elements[i][i])
            Q.set_column(v, i)
        
        return Q, R

    def qr_algorithm(self, max_iter):
        q, r = self.qr_decomposition()
        for _ in range(max_iter):
            A_i = r * q
            q, r = A_i.qr_decomposition()
            
        eigenvalues = [A_i._elements[i][i] for i in range(A_i.row)]
        
        
        vectors = []
        
        for i in range(q.column):
            e = Matrix(q.column, 1, [0 if i!=j else 1 for j in range(q.column)])
            print(e)
            vectors.append(q * e)
        
        return eigenvalues, vectors
        

class SVDMixin:
    @staticmethod
    def F(X, a, b):
        err = 0
        for i in range(X.row):
            for j in range(X.column):
                err += (X[i, j] - b[i, 0] * a[j, 0]) ** 2
        return err / 2
    
    def svd(self):
        n = self.row
        m = self.column

        X = self.copy()
        u = Matrix(m, 1, [random.randint(0,1) for i in range(m)])
        v = Matrix(n, 1, [0]*n)

        U = Matrix(m, 0)
        V = Matrix(n, 0)
        s = Matrix(n, m, [0] * (n*m))

        u_list = []
        v_list = []
        s_list = []
        last = 1
        curr = 0

        for _ in range(max(n,m)):
            
            u = Matrix(m, 1, [random.randint(1,1) for i in range(m)])
            last = 1
            curr = 0
            
            while (last-curr) > 0.0000001:
                last = F(X, u, v)
                
                for i in range(n):
                    sum_x_u = 0
                    sum_u2 = 0
                    for j in range(m):
                        sum_x_u += X[i, j] * u[j, 0] 
                        sum_u2 += u[j, 0] ** 2
                    v[i, 0] = sum_x_u / sum_u2

                for i in range(m):
                    sum_x_v = 0
                    sum_v2 = 0
                    for j in range(n):
                        sum_x_v += X[j, i] * v[j, 0] 
                        sum_v2 +=  v[j, 0] ** 2
                    u[i, 0] = sum_x_v / sum_v2
                
                curr = F(X, u, v)
            
            X = X - (v * u.T)
            
            s_list.append(u.norm() * v.norm())
            u.divide_by_num(u.norm())
            v.divide_by_num(v.norm())
            u_list.append(u.copy())
            v_list.append(v.copy())
             
        for i in range(m):
            U.add_column(u_list[i])

        for i in range(n):
            V.add_column(v_list[i])

        for i in range(min(n,m)):
            s[i, i] = s_list[i]
    
        return U, s, V


class Matrix(SVDMixin, EigenvalueEigenvectorMixin, GrevilleMethod, GaussMethodMixin, ElementaryTransformationsMixin):
    def __init__(self, row:int, column:int, elements:list|None = None) -> None:
        self._row = row
        self._column = column
        self._elements = []
        if elements:
            index = 0
            for i in range(row):
                self._elements.append([])
                for _ in range(column):
                    self._elements[i].append(elements[index])
                    index += 1
                    
            #self.to_fractions()
            
    @property
    def row(matrix) -> int:
        return matrix._row
    
    @property
    def column(matrix) -> int:
        return matrix._column
    
    def _is_zeros(self, row:list) -> bool:
        for el in row:
            if el != 0:
                return False
        return True
    
    @benchmark     
    def rank(self) -> int:
        triangle_matr, _ = self.get_triangle()
        print(triangle_matr)
        rank = 0
        for i in range(self.row):
            if triangle_matr._elements[i][i] != 0 or not triangle_matr._is_zeros(triangle_matr._elements[i]):
                rank += 1
        return rank        
                
    def copy(self):
        matr = Matrix(self.row, self.column,)
        matr._elements = [row.copy() for row in self._elements]
        return matr

    def get_column(self, columns, num_columns):
        temp = []
        if isinstance(columns, int):
            for row in range(self.row):
                temp.append(self._elements[row][columns])
        else:
            for row in range(self.row):
                for c in columns:
                    temp.append(self._elements[row][c])
        return Matrix(self.row, num_columns, temp)
    
    def get_row(self, row):
        return Matrix(1, self.column, self._elements[row])
    
    def add_row(self, row):
        self._row += 1
        self._elements.append(row._elements[0])
    
    def add_column(self, column):
        self._column += 1
        if not self._elements:
            for i in range(self.row):
                self._elements.append([])
                self._elements[i].append(column[i, 0])
            return
        for i in range(self.row):
            self._elements[i].append(column[i, 0])
    
    @benchmark
    def get_skeleton_decomposition(self):
        triangle_matr, _ = self.get_triangle()
        temp = []
        rank_A = 0
        for row in range(triangle_matr.row):
            if not self._is_zeros(triangle_matr._elements[row]):
                rank_A += 1
                temp.extend(self._elements[row])
        C = Matrix(rank_A, self.column, temp)
        C_inv = C.get_pseudoinverse_matrix()
        B = self * C_inv
        return B, C
    
    @property
    def T(self):
        temp = []
        for i in range(self.column):
            for j in range(self.row):
                temp.append(self._elements[j][i])
        return Matrix(self.column, self.row, temp)
        
    def __mul__(self, matr):
        if self.column != matr.row:
            raise Exception('column != row')
        temp = []
        for i in range(self.row):
            for j in range(matr.column):
                el = 0
                for k in range(self.column):
                    el += self._elements[i][k] * matr._elements[k][j]
                temp.append(el)
        return Matrix(self.row, matr.column, temp)
    
    def __sub__(self, a):
        if self.row != a.row or self.column != a.column:
            return None
        temp = []
        for i in range(self.row):
            for j in range(self.column):
                temp.append(self._elements[i][j] - a._elements[i][j])
        return Matrix(self.row, self.column, temp)
    
    def __pow__(self, a):
        temp = []
        for i in range(self.row):
            for j in range(self.column):
                temp.append(self._elements[i][j] ** a)
        return Matrix(self.row, self.column, temp)
    
    def __str__(self) -> str:
        for i in range(self.row):
            for j in range(self.column):
                self[i, j] = round(self[i, j], 3) 
        return '\n'.join(' | '.join(map(str, row)) for row in self._elements) 
    
    def __getitem__(self, pos):
        i, j = pos
        return self._elements[i][j]
    
    def __setitem__(self, pos, v):
        i, j = pos
        self._elements[i][j] = v
                    

n = 2
m = 3

X = Matrix(n, m, [2, -1, 0, 
                4, 3, -2,])

U, s, V = X.svd()

print(U)
print()
print(V.T)
print()
print(s)

print((U * s.get_pseudoinverse_matrix() * V.T))

matr = np.array([[2, -1, 0], 
                [4, 3, -2],])

q,w,e = np.linalg.svd(matr)

print(q,w,e, sep='\n')

print(np.linalg.pinv(matr))