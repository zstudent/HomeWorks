import numpy as np
import os, time
os.system('cls')

class Matrix:
    # This class with its methods is going to help us while multiplying matrices

    def __init__(self,matrix):
        self.matrix = matrix
    
    # This method returns the specific row of a matrix
    def get_row(self,index):
        try:
            return self.matrix[index][:]
        except:
            return "Index is out of bounds"

    # This method returns the specific column of a matrix
    def get_column(self,index):
        try:
            column = []
            for row in self.matrix:
                column.append(row[index])
            return column
        except:
            return "Index is out of bounds"
    
    # This method is responsible for performing dot product
    # Example of dot product: a = [1,2], b = [3,4] a*b = 1*3 + 2*4 = 3 + 8 = 11
    def dot_product(a='',b=''):
        try:
            return sum([a_i*b_i for a_i, b_i in zip(a,b)])
        except:
            return "Both of vectors must be provided"
    
    # This method returns a boolean value that indicates whether matrix multiplication can be performed or not
    def is_multiplicable(a,b):
        first_matrix, second_matrix = Matrix(a), Matrix(b)
        for i in range(len(a)):
            row, column = first_matrix.get_row(i), second_matrix.get_column(i)
            if len(row) != len(column): 
                return False
        return True

# This function is responsible for multiplying two matrices
def matrix_multiplication(A,B):
    if not Matrix.is_multiplicable(A,B):
        return 'Matrix shapes are inconsistent'
    row_count, column_count = len(A), len(B[:][0])   
    first_matrix, second_matrix = Matrix(A), Matrix(B)
    result_matrix = []
    for i in range(row_count):
        row_results = []
        for j in range(column_count):
            row, column = first_matrix.get_row(i), second_matrix.get_column(j)
            row_results.append(Matrix.dot_product(row,column))
        result_matrix.append(row_results)
    return result_matrix

A = [[100,22,35,433,54], # 3x5 matrix
    [33,63,34,-130,100],
    [5,7,10,-50,1000]]

B = [[1,200,3,4,100,200], # 5x6 matrix
    [5,6,7,8,-1000,-2000],
    [10,-205,30,40,5,8],
    [-5,-6,-7,-8,1,5],
    [0.1,-0.1,2,4,6,-240]]

# the result should be 3 x 6 matrix

my_function_time = time.time()
print('My result: {0}'.format(matrix_multiplication(A,B)))
print('My function performance: ',"--- %s seconds ---" % (time.time() - my_function_time))

numpy_method_time = time.time()
print('NumPy result: {0}'.format(np.dot(A,B)))
print("NumPy method performance","--- %s seconds ---" % (time.time() - numpy_method_time))