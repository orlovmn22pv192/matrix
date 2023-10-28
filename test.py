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
