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

class ElementaryЕransformationsMixin:
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
        matr.to_fractions()
        count_swap = 0
        
        for i in range(matr.row):
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


class Matrix(GaussMethodMixin, ElementaryЕransformationsMixin):
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
        rank = 0
        for i in range(self.row):
            if triangle_matr._elements[i][i] != 0 or not triangle_matr._is_zeros(triangle_matr._elements[i]):
                rank += 1
        return rank        
                
    def copy(self):
        matr = Matrix(self.row, self.column,)
        matr._elements = [row.copy() for row in self._elements]
        return matr

    def __str__(self) -> str:
        return '\n'.join(' | '.join(map(str, row)) for row in self._elements)  
                    

# rank 2
test1 = (3,3,[1,2,3,4,5,6,7,8,9]) 
# rank 3
test2 = (3,3,[2,1,-1,-3,-1,2,-2,1,2])
# rank 4
test3 = (4,4,[1,2,3,14,2,3,4,15,3,4,51,6,4,5,16,7])
# rank 5
test4 = (5,5,[1,20,3,4,51,34,3,56,-2,3,43,0,1,67,78,1,34,56,87,12,11,22,33,44,55])

tests = zip([test1, test2, test3, test4], [2,3,4,5])

for test, answer in tests:
    print('res=', Matrix(*test).rank(), ' answer=', answer)

@benchmark
def test_np(matr) -> int:
    m = np.array(matr)    
    rank = np.linalg.matrix_rank(m)
    return rank

n = 100
m = 100
    
huge_matrix = [[random.randint(1,100) for _ in range(m)] for _ in range(n)]
print('Numpy:')
print('rank = ', test_np(huge_matrix))
m = Matrix(n,m)
m._elements = huge_matrix
print('Gauss method:')
print('rank = ',m.rank())


