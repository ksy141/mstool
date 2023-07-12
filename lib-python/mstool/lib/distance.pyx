import   numpy as np
cimport  numpy
cimport  cython
from libc.math cimport round as cround

ctypedef numpy.float64_t DTYPE_t
# https://cython-docs2.readthedocs.io/en/latest/src/tutorial/numpy.html

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def distance_matrix(numpy.ndarray[DTYPE_t, ndim=2] pos1,
                    numpy.ndarray[DTYPE_t, ndim=2] pos2,
                    numpy.ndarray[DTYPE_t, ndim=1] dimensions,
                    bint pbc):
    
    # Py_ssize_t is a signed integer type provided by Python 
    # which covers the same range of values as is supported as NumPy array indices. 
    # It is the preferred type to use for loops over arrays.
    
    cdef Py_ssize_t x_max = pos1.shape[0]
    cdef Py_ssize_t y_max = pos2.shape[0]
    cdef Py_ssize_t x = 0
    cdef Py_ssize_t y = 0
    cdef Py_ssize_t i = 0
    cdef double    dx = 0
    cdef double    dy = 0
    cdef double    dz = 0
    cdef double   ddx = 0
    cdef double   ddy = 0
    cdef double   ddz = 0
    
    cdef numpy.ndarray[DTYPE_t, ndim=2] result = np.zeros([x_max, y_max], dtype=np.float64)
    
    cdef DTYPE_t value
    for x in range(x_max):
        for y in range(y_max):

            dx = pos1[x,0] - pos2[y,0]
            dy = pos1[x,1] - pos2[y,1]
            dz = pos1[x,2] - pos2[y,2]
            
            if pbc:
                ddx = dimensions[0] * cround(dx/dimensions[0])
                ddy = dimensions[1] * cround(dy/dimensions[1])
                ddz = dimensions[2] * cround(dz/dimensions[2])
                
            dx -= ddx
            dy -= ddy
            dz -= ddz

            value = (dx**2 + dy**2 + dz**2) ** 0.5
            result[x, y] = value
            
    return result

