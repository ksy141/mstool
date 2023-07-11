import numpy as np

def distance_matrix(pos1, pos2, dimensions=None):
    """Calculate all possible distances between pos1 and pos2.
    If dimensions is not given, the returned value should be the same as 
    ``scipy.spatial.distance_matrix(pos1, pos2)``.
    This is memory-inefficient. I should update it to cython.
    
    Parameters
    ----------
    pos1 : numpy.ndarray
        array of shape ``(n, 3)``
    pos2 : numpy.ndarray
        array of shape ``(m, 3)``
    dimensions : numpy.ndarray or list, optional
        array or list in the format of [pbcx, pbcy, pbcz, 90, 90, 90].
        Periodic bounday condition will be taken into account if dimensions is provided.
    Returns
    -------
    dis : numpy.ndarray
        a distance matrix of shape ``(n, m)``

    Examples
    --------
    >>> import MDAnalysis as mda
    >>> from   MDAnalysis.analysis.distances import distance_array
    >>> from   mstool.utils.distance import distance_matrix
    >>> u = mda.Universe('cg.pdb')
    >>> pos = u.atoms.positions
    >>> dim = u.dimensions
    >>> dm1 = distance_array(pos, pos, box=dim)
    >>> dm2 = distance_matrix(pos, pos, dimensions=dim)
    >>> np.all(np.isclose(dm1, dm2), atol=1e-3)
    """

    # if any of the values of dimensions is 0, pbc will not be taken into account
    if np.all(dimensions):
        print("PBC takein into account")
        dimensions = np.array(dimensions, dtype=np.float64)
        arr = np.repeat(pos1, len(pos2), axis=0) - np.tile(pos2, (len(pos1), 1))
        arr = arr - dimensions[0:3] * np.round(arr / dimensions[0:3])
        dis = np.sqrt(np.sum(arr**2, axis=1)).reshape(len(pos1), len(pos2))
        print(dis.shape)
    
    else:
        print("PBC not takein into account")
        print(pos1.shape, pos2.shape)
        import scipy
        dis  = scipy.spatial.distance_matrix(pos1, pos2)
        print(dis.shape)

    return dis

