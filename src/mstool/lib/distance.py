import numpy as np
from   .distancelib import distance_matrix_lib, distance_overlap_lib

def distance_matrix(pos1, pos2, dimensions=None):
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
    >>> dm2 = distance_matrix(pos1, pos2, dimensions=dim)
    >>> np.all(np.isclose(dm1, dm2, atol=1e-3))
    '''

    pos1 = np.asarray(pos1, dtype=np.float64)
    pos2 = np.asarray(pos2, dtype=np.float64)
    
    if dimensions is not None:
        dimensions = np.asarray(dimensions, dtype=np.float64)
        dm = distance_matrix_lib(pos1, pos2, dimensions, True)
    else:
        dimensions = np.zeros(6, dtype=np.float64)
        dm = distance_matrix_lib(pos1, pos2, dimensions, False)

    return dm



def distance_overlap(pos1, pos2, rcut, dimensions=None):
    '''Calculate overlap between pos1 and pos2 within a provided cutoff distance.
    pos1 is the coordinates of atoms that you do not want to keep if overlapped with pos2 (e.g., water atoms)
    pos2 is the coordinates of atoms that you do want to keep even if overlapped with pos1 (e.g., protein atoms)
    If dimensions is given, PBC will be taken int account.

    Parameters
    ----------
    pos1 : numpy.ndarray
        array of shape ``(n, 3)`` e.g., positions of water atoms
    pos2 : numpy.ndarray
        array of shape ``(m, 3)`` e.g., positions of protein atoms
    rcut : float
        cutoff distance to define overlap 
    dimensions : numpy.array
        array of shape ``(6,)`` or (3,)``. 
        Only the first three values will be used.
        The default is None (PBC will not be taken into account).

    Returns
    -------
    result : numpy.ndarray(dtype=np.bool)
        bool array of ``(n,)``.
        ``pos1[result]`` are the water atoms that you want to save.

    Examples
    --------
    >>> from   mstool.lib.distance import distance_overlap
    >>> pos1 = np.random.rand(50,  3) # solvent
    >>> pos2 = np.random.rand(10, 3)  # solute
    >>> dim  = np.array([0.2, 0.2, 0.2, 90, 90, 90])
    >>> bA   = distance_overlap(pos1, pos2, rcut=0.05, dimensions=dim)
    >>> pos1[bA] # solvent atoms w/o overlap with solute
    '''

    pos1 = np.asarray(pos1, dtype=np.float64)
    pos2 = np.asarray(pos2, dtype=np.float64)
    rcut = float(rcut)

    if dimensions is not None:
        dimensions = np.asarray(dimensions, dtype=np.float64)
        bA = distance_overlap_lib(pos1, pos2, rcut, dimensions, True)
    else:
        dimensions = np.zeros(6, dtype=np.float64)
        bA = distance_overlap_lib(pos1, pos2, rcut, dimensions, False)

    return bA
