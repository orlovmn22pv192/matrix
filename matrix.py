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


class Matrix(GrevilleMethod, GaussMethodMixin, ElementaryTransformationsMixin):
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
                    
            self.to_fractions()
            
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
    
    def add_row(self, row):
        self._row += 1
        self._elements.append(row._elements[0])
    
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
        return '\n'.join(' | '.join(map(str, row)) for row in self._elements)  
                    

# rank 2
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

