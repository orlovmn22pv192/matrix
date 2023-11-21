from matrix import Matrix
import math
import matplotlib.pyplot as plt
import numpy as np

#Lab 7

class AffineTransformations:
    def get_rotate_matrix_x(a):
        return Matrix(4, 4, [1, 0, 0, 0,
                             0, math.cos(a), -math.sin(a), 0,
                             0, math.sin(a), math.cos(a), 0,
                             0, 0, 0, 1])
    
    def get_rotate_matrix_y(a):
        return Matrix(4, 4, [math.cos(a), 0, -math.sin(a), 0,
                             0, 1, 0, 0,
                             math.sin(a), 0, math.cos(a), 0,
                             0, 0, 0, 1])
        
    def get_rotate_matrix_z(a):
        return Matrix(4, 4, [math.cos(a), math.sin(a), 0, 0,
                             -math.sin(a), math.cos(a), 0, 0,
                             0, 0, 1, 0,
                             0, 0, 0, 1])
    
    def get_scale_matrix(k_x, k_y, k_z):
        return Matrix(4, 4, [k_x, 0, 0, 0,
                             0, k_y, 0, 0,
                             0, 0, k_z, 0,
                             0, 0,   0, 1])
    
    def get_transfer_matrix(x, y, z):
        return Matrix(4, 4, [1, 0, 0, x,
                             0, 1, 0, y,
                             0, 0, 1, z,
                             0, 0, 0, 1]) 
        

'''fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

x = [0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0,   ]
y = [0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1,  ]
z = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,  ]
ax.plot(x,y,z, color='r')
plt.show()'''

import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
 
vertices = np.zeros([3,8],dtype=int)
vertices[0,:] = [1, 7, 5, 8, 2, 4, 6, 3]
vertices[1,:] = [1, 7, 4, 6, 8, 2, 5, 3]
vertices[2,:] = [6, 1, 5, 2, 8, 3, 7, 4]
vertices = vertices - 1 #(adjust the indices by one since python starts with zero indexing)

# Define an array with dimensions 8 by 3
# 8 for each vertex
# -> indices will be vertex1=0, v2=1, v3=2 ...
# 3 for each coordinate
# -> indices will be x=0,y=1,z=1
cube = np.zeros([8,3])

# Define x values
cube[:,0] = [0, 0, 0, 0, 1, 1, 1, 1]
# Define y values
cube[:,1] = [0, 1, 0, 1, 0, 1, 0, 1]
# Define z values
cube[:,2] = [0, 0, 1, 1, 0, 0, 1, 1]

# First initialize the fig variable to a figure
fig = plt.figure()
# Add a 3d axis to the figure
ax = fig.add_subplot(111, projection='3d')

# plotting cube
# Initialize a list of vertex coordinates for each face
# faces = [np.zeros([5,3])]*3
faces = []
faces.append(np.zeros([5,3]))
faces.append(np.zeros([5,3]))
faces.append(np.zeros([5,3]))
faces.append(np.zeros([5,3]))
faces.append(np.zeros([5,3]))
faces.append(np.zeros([5,3]))
# Bottom face
faces[0][:,0] = [0,0,1,1,0]
faces[0][:,1] = [0,1,1,0,0]
faces[0][:,2] = [0,0,0,0,0]
# Top face
faces[1][:,0] = [0,0,1,1,0]
faces[1][:,1] = [0,1,1,0,0]
faces[1][:,2] = [1,1,1,1,1]
# Left Face
faces[2][:,0] = [0,0,0,0,0]
faces[2][:,1] = [0,1,1,0,0]
faces[2][:,2] = [0,0,1,1,0]
# Left Face
faces[3][:,0] = [1,1,1,1,1]
faces[3][:,1] = [0,1,1,0,0]
faces[3][:,2] = [0,0,1,1,0]
# front face
faces[4][:,0] = [0,1,1,0,0]
faces[4][:,1] = [0,0,0,0,0]
faces[4][:,2] = [0,0,1,1,0]
# front face
faces[5][:,0] = [0,1,1,0,0]
faces[5][:,1] = [1,1,1,1,1]
faces[5][:,2] = [0,0,1,1,0]
ax.add_collection3d(Poly3DCollection(faces, facecolors='cyan', linewidths=1, edgecolors='k', alpha=.25))

# plotting lines


plt.show()