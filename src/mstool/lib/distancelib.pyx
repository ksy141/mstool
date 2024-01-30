# cython: language_level=3
import   numpy as np
cimport  numpy
cimport  cython
from libc.math cimport round as cround

ctypedef numpy.float64_t DTYPE_t
ctypedef numpy.int64_t   INT_t
# https://cython-docs2.readthedocs.io/en/latest/src/tutorial/numpy.html

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def distance_matrix_lib(numpy.ndarray[DTYPE_t, ndim=2] pos1,
                        numpy.ndarray[DTYPE_t, ndim=2] pos2,
                        numpy.ndarray[DTYPE_t, ndim=1] dimensions,
                        bint pbc):
    '''Calculate all the possible distances between pos1 and pos2.
    If pbc=True, periodic boundary condition (PBC) will be taken into account.
    Please note that the input arrays must be in ``np.float64``
    because it uses cython.

    Parameters
    ----------
    pos1 : numpy.ndarray(dtype=np.float64)
        array of shape ``(n, 3)``
    pos2 : numpy.ndarray(dtype=np.float64)
        array of shape ``(m, 3)``
    dimensions : numpy.aray(dtype=np.float64)
        array of shape ``(6,)`` or (3,)``. 
        Only the first three dimensions will be used.
    pbc : bool
        If ``pbc=True``, PBC will be taken into account.

    Returns
    -------
    result : numpy.ndarray(dtype=np.float64)
        distance matrix of shape ``(n, m)``

    Examples
    --------
    >>> from   MDAnalysis.analysis.distances import distance_array
    >>> from   mstool.lib.distance           import distance_matrix
    >>> pos1 = np.random.rand(5,  3)
    >>> pos2 = np.random.rand(10, 3)
    >>> dim  = np.array([0.8, 0.8, 0.8, 90, 90, 90])
    >>> dm1 = distance_array(pos1, pos2,  box=dim)
    >>> dm2 = distance_matrix(pos1, pos2, dimensions=dim, pbc=True)
    >>> np.all(np.isclose(dm1, dm2, atol=1e-3))
    '''

    # Py_ssize_t is a signed integer type provided by Python 
    # which covers the same range of values as is supported as NumPy array indices. 
    # It is the preferred type to use for loops over arrays.
    
    cdef Py_ssize_t x_max = pos1.shape[0]
    cdef Py_ssize_t y_max = pos2.shape[0]
    cdef Py_ssize_t x = 0
    cdef Py_ssize_t y = 0
    cdef Py_ssize_t i = 0
    cdef double    dx = 0.0
    cdef double    dy = 0.0
    cdef double    dz = 0.0
    cdef double   ddx = 0.0
    cdef double   ddy = 0.0
    cdef double   ddz = 0.0
    
    cdef numpy.ndarray[DTYPE_t, ndim=2] result = np.zeros([x_max, y_max], dtype=np.float64)
    
    cdef DTYPE_t value = 0.0
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


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def distance_overlap_lib(numpy.ndarray[DTYPE_t, ndim=2] pos1,
                         numpy.ndarray[DTYPE_t, ndim=2] pos2,
                         double rcut,
                         numpy.ndarray[DTYPE_t, ndim=1] dimensions,
                         bint pbc):
    '''Calculate overlap between pos1 and pos2 within a provided cutoff distance.
    pos1 is the coordinates of atoms that you do not want to keep if overlapped with pos2 (e.g., water atoms)
    pos2 is the coordinates of atoms that you do want to keep even if overlapped with pos1 (e.g., protein atoms)
    If pbc=True, periodic boundary condition (PBC) will be taken into account.
    Please note that the input arrays must be in ``np.float64``
    because it uses cython.

    Parameters
    ----------
    pos1 : numpy.ndarray(dtype=np.float64)
        array of shape ``(n, 3)`` e.g., positions of water atoms
    pos2 : numpy.ndarray(dtype=np.float64)
        array of shape ``(m, 3)`` e.g., positions of protein atoms
    rcut : float
        cutoff distance to define overlap 
    dimensions : numpy.aray(dtype=np.float64)
        array of shape ``(6,)`` or (3,)``. 
        Only the first three dimensions will be used.
    pbc : bool
        If ``pbc=True``, PBC will be taken into account.

    Returns
    -------
    result : numpy.ndarray(dtype=np.bool)
        bool array of ``(n,)``.
        ``pos1[result]`` are the water atoms that you want to save.

    Examples
    --------
    >>> from   mstool.lib.distance_overlap import distance_overlap
    >>> pos1 = np.random.rand(50,  3) # solvent
    >>> pos2 = np.random.rand(10, 3)  # solute
    >>> dim  = np.array([0.2, 0.2, 0.2, 90, 90, 90])
    >>> bA   = distance_overlap(pos1, pos2, rcut=0.05, dim=dim, pbc=True)
    >>> pos1[bA] # solvent atoms w/o overlap with solute
    '''

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
    
    cdef numpy.ndarray[INT_t, ndim=1] result = np.zeros([x_max], dtype=np.int64) + 1
    
    cdef DTYPE_t value2 = 0.0
    cdef DTYPE_t rcut2  = rcut ** 2
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

            value2 = dx**2 + dy**2 + dz**2

            if value2 < rcut2:
                result[x] = 0
                break
    
    return result != 0 

