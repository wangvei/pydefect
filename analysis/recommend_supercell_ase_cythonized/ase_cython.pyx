from distutils.core import setup
from Cython.Build import cythonize
from cpython cimport array
import array
import numpy as np
cimport cython
from ase.build import get_deviation_from_optimal_cell_shape
import math

cdef c_matrix_from_python_list(list python_list, double c_array[3][3]):
    for i in range(3):
        for j in range(3):
            c_array[i][j] = python_list[i][j]

cdef addition_c_matrix(double array1[3][3], double array2[3][3], double answer[3][3]):
    for i in range(3):
        for j in range(3):
            answer[i][j] = array1[i][j] + array2[i][j]

cdef python_list_from_c_matrix(double c_array[3][3]):
    python_list = [[0,0,0],[0,0,0],[0,0,0]]
    for i in range(3):
        for j in range(3):
            python_list[i][j] = c_array[i][j]
    return python_list

cdef print_c_matrix(double array[3][3]):
    for i in range(3):
        for j in range(3):
            print("i, j, array[i][j] = ", i, j, array[i][j])

cdef void _find_optimal_cell_shape_cdef(int lower_limit, int upper_limit, 
                                        double h[3][3], # norm_cell
                                        double *best_score,
                                        double optimal_P[3][3], 
                                        double target_metric[3][3],
                                        int target_size,
                                        double starting_P[3][3]):
    """
    Mimic ase code by Cython.
    """

    # Set up target metric
    cdef (double[3][3]) P, d
    cdef double det_P = 0
    cdef double m
    cdef double current_score
    cdef int dxx, dxy, dxz, dyx, dyy, dyz, dzx, dzy, dzz
    best_score[0] = 1e+10

    for dxx in range(lower_limit, upper_limit + 1):
        for dxy in range(lower_limit, upper_limit + 1):
            for dxz in range(lower_limit, upper_limit + 1):
                for dyx in range(lower_limit, upper_limit + 1):
                    for dyy in range(lower_limit, upper_limit + 1):
                        for dyz in range(lower_limit, upper_limit + 1):
                            for dzx in range(lower_limit, upper_limit + 1):
                                for dzy in range(lower_limit, upper_limit + 1):
                                    for dzz in range(lower_limit, upper_limit + 1):
                                        d[0][0] = dxx
                                        d[0][1] = dxy
                                        d[0][2] = dxz
                                        d[1][0] = dyx
                                        d[1][1] = dyy
                                        d[1][2] = dyz
                                        d[2][0] = dzx
                                        d[2][1] = dzy
                                        d[2][2] = dzz
                                        addition_c_matrix(starting_P, d, P)
                                        det_P = 0
                                        det_P += P[0][0] * P[1][1] * P[2][2]                                                                                                    
                                        det_P += P[0][1] * P[1][2] * P[2][0]                                                                                                    
                                        det_P += P[0][2] * P[1][0] * P[2][1]                                                                                                    
                                        det_P -= P[0][0] * P[1][2] * P[2][1]                                                                                                    
                                        det_P -= P[0][1] * P[1][0] * P[2][2]                                                                                                    
                                        det_P -= P[0][2] * P[1][1] * P[2][0]

                                        if det_P == target_size:
                                            current_score = 0.0
                                            current_score += pow(P[0][0] * h[0][0] + P[0][1] * h[1][0] + P[0][2] * h[2][0] - target_metric[0][0], 2)
                                            current_score += pow(P[0][0] * h[0][1] + P[0][1] * h[1][1] + P[0][2] * h[2][1] - target_metric[0][1], 2)
                                            current_score += pow(P[0][0] * h[0][2] + P[0][1] * h[1][2] + P[0][2] * h[2][2] - target_metric[0][2], 2)
                                            current_score += pow(P[1][0] * h[0][0] + P[1][1] * h[1][0] + P[1][2] * h[2][0] - target_metric[1][0], 2)
                                            current_score += pow(P[1][0] * h[0][1] + P[1][1] * h[1][1] + P[1][2] * h[2][1] - target_metric[1][1], 2)
                                            current_score += pow(P[1][0] * h[0][2] + P[1][1] * h[1][2] + P[1][2] * h[2][2] - target_metric[1][2], 2)
                                            current_score += pow(P[2][0] * h[0][0] + P[2][1] * h[1][0] + P[2][2] * h[2][0] - target_metric[2][0], 2)
                                            current_score += pow(P[2][0] * h[0][1] + P[2][1] * h[1][1] + P[2][2] * h[2][1] - target_metric[2][1], 2)
                                            current_score += pow(P[2][0] * h[0][2] + P[2][1] * h[1][2] + P[2][2] * h[2][2] - target_metric[2][2], 2)


                                            if current_score < best_score[0]:
                                                best_score[0] = current_score
                                                for i in range(3):
                                                    for j in range(3):
                                                        optimal_P[i][j] = P[i][j]

def find_optimal_cell_shape(cell, target_size,
                            lower_limit=-2, upper_limit=2):
    target_metric = np.eye(3)
    norm = (target_size * np.linalg.det(cell) /
                        np.linalg.det(target_metric))**(-1.0 / 3)
    norm_cell = norm * cell
    ideal_P = np.dot(target_metric, np.linalg.inv(norm_cell))
    starting_P = np.array(np.around(ideal_P, 0), dtype=float)
    optimal_P = np.zeros(9).reshape(3,3)
    cdef (double[3][3]) c_norm_cell, c_optimal_P, c_target_metric, c_starting_P
    cdef double[1] best_score

    c_matrix_from_python_list(norm_cell.tolist(), c_norm_cell)
    c_matrix_from_python_list(target_metric.tolist(), c_target_metric)
    c_matrix_from_python_list(optimal_P.tolist(), c_optimal_P)
    c_matrix_from_python_list(starting_P.tolist(), c_starting_P)

    lower_limit = int(lower_limit)
    upper_limit = int(upper_limit)
    target_size = int(target_size)

    _find_optimal_cell_shape_cdef(lower_limit, upper_limit, 
                                  c_norm_cell, best_score,
                                  c_optimal_P, c_target_metric,
                                  target_size, c_starting_P)
    optimal_P = python_list_from_c_matrix(c_optimal_P)
    return optimal_P

