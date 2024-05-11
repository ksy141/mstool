import numpy as np

def isfloat(string):
    try:
        float(string.strip())
        return True
    except:
        return False

### the below is copied from MDAnalysis 2.4.1
def ltruncate_int(value, ndigits):
    """Truncate an integer, retaining least significant digits

    Parameters
    ----------
    value : int
      value to truncate
    ndigits : int
      number of digits to keep

    Returns
    -------
    truncated : int
      only the `ndigits` least significant digits from `value`

    Examples
    --------
    >>> ltruncate_int(123, 2)
    23
    >>> ltruncate_int(1234, 5)
    1234
    """
    return int(str(value)[-ndigits:])

def triclinic_box(x, y, z):
    """Convert the three triclinic box vectors to
    ``[lx, ly, lz, alpha, beta, gamma]``.

    If the resulting box is invalid, i.e., any box length is zero or negative,
    or any angle is outside the open interval ``(0, 180)``, a zero vector will
    be returned.

    All angles are in degrees and defined as follows:

    * ``alpha = angle(y,z)``
    * ``beta  = angle(x,z)``
    * ``gamma = angle(x,y)``

    Parameters
    ----------
    x : array_like
        Array of shape ``(3,)`` representing the first box vector
    y : array_like
        Array of shape ``(3,)`` representing the second box vector
    z : array_like
        Array of shape ``(3,)`` representing the third box vector

    Returns
    -------
    numpy.ndarray
        A numpy array of shape ``(6,)`` and dtype ``np.float32`` providing the
        unitcell dimensions in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:\n
        ``[lx, ly, lz, alpha, beta, gamma]``.\n
        Invalid boxes are returned as a zero vector.

    Note
    ----
    Definition of angles: http://en.wikipedia.org/wiki/Lattice_constant

    See Also
    --------
    :func:`~MDAnalysis.lib.mdamath.triclinic_vectors`


    .. versionchanged:: 0.20.0
       Calculations are performed in double precision and invalid box vectors
       result in an all-zero box.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    lx = np.linalg.norm(x)
    ly = np.linalg.norm(y)
    lz = np.linalg.norm(z)
    if lx == 0 or ly ==0 or lz == 0:
        return np.zeros(6, dtype=np.float32)
    alpha = np.rad2deg(np.arccos(np.dot(y, z) / (ly * lz)))
    beta = np.rad2deg(np.arccos(np.dot(x, z) / (lx * lz)))
    gamma = np.rad2deg(np.arccos(np.dot(x, y) / (lx * ly)))
    box = np.array([lx, ly, lz, alpha, beta, gamma], dtype=np.float32)
    # Only positive edge lengths and angles in (0, 180) are allowed:
    if np.all(box > 0.0) and alpha < 180.0 and beta < 180.0 and gamma < 180.0:
        return box
    # invalid box, return zero vector:
    return np.zeros(6, dtype=np.float32)

def triclinic_vectors(dimensions, dtype = np.float32):
    """Convert ``[lx, ly, lz, alpha, beta, gamma]`` to a triclinic matrix
    representation.

    Original `code by Tsjerk Wassenaar`_ posted on the Gromacs mailinglist.

    If `dimensions` indicates a non-periodic system (i.e., all lengths are
    zero), zero vectors are returned. The same applies for invalid `dimensions`,
    i.e., any box length is zero or negative, or any angle is outside the open
    interval ``(0, 180)``.

    .. _code by Tsjerk Wassenaar:
       http://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html

    Parameters
    ----------
    dimensions : array_like
        Unitcell dimensions provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:\n
        ``[lx, ly, lz, alpha, beta, gamma]``.
    dtype: numpy.dtype
        The data type of the returned box matrix.

    Returns
    -------
    box_matrix : numpy.ndarray
        A numpy array of shape ``(3, 3)`` and dtype `dtype`,
        with ``box_matrix[0]`` containing the first, ``box_matrix[1]`` the
        second, and ``box_matrix[2]`` the third box vector.

    Notes
    -----
    * The first vector is guaranteed to point along the x-axis, i.e., it has the
      form ``(lx, 0, 0)``.
    * The second vector is guaranteed to lie in the x/y-plane, i.e., its
      z-component is guaranteed to be zero.
    * If any box length is negative or zero, or if any box angle is zero, the
      box is treated as invalid and an all-zero-matrix is returned.


    .. versionchanged:: 0.7.6
       Null-vectors are returned for non-periodic (or missing) unit cell.
    .. versionchanged:: 0.20.0
       * Calculations are performed in double precision and zero vectors are
         also returned for invalid boxes.
       * Added optional output dtype parameter.
    """
    dim = np.asarray(dimensions, dtype=np.float64)
    lx, ly, lz, alpha, beta, gamma = dim
    # Only positive edge lengths and angles in (0, 180) are allowed:
    if not (np.all(dim > 0.0) and
            alpha < 180.0 and beta < 180.0 and gamma < 180.0):
        # invalid box, return zero vectors:
        box_matrix = np.zeros((3, 3), dtype=dtype)
    # detect orthogonal boxes:
    elif alpha == beta == gamma == 90.0:
        # box is orthogonal, return a diagonal matrix:
        box_matrix = np.diag(dim[:3].astype(dtype, copy=False))
    # we have a triclinic box:
    else:
        box_matrix = np.zeros((3, 3), dtype=np.float64)
        box_matrix[0, 0] = lx
        # Use exact trigonometric values for right angles:
        if alpha == 90.0:
            cos_alpha = 0.0
        else:
            cos_alpha = np.cos(np.deg2rad(alpha))
        if beta == 90.0:
            cos_beta = 0.0
        else:
            cos_beta = np.cos(np.deg2rad(beta))
        if gamma == 90.0:
            cos_gamma = 0.0
            sin_gamma = 1.0
        else:
            gamma = np.deg2rad(gamma)
            cos_gamma = np.cos(gamma)
            sin_gamma = np.sin(gamma)
        box_matrix[1, 0] = ly * cos_gamma
        box_matrix[1, 1] = ly * sin_gamma
        box_matrix[2, 0] = lz * cos_beta
        box_matrix[2, 1] = lz * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        box_matrix[2, 2] = np.sqrt(lz * lz - box_matrix[2, 0] ** 2 -
                                   box_matrix[2, 1] ** 2)
        # The discriminant of the above square root is only negative or zero for
        # triplets of box angles that lead to an invalid box (i.e., the sum of
        # any two angles is less than or equal to the third).
        # We don't need to explicitly test for np.nan here since checking for a
        # positive value already covers that.
        if box_matrix[2, 2] > 0.0:
            # all good, convert to correct dtype:
            box_matrix = box_matrix.astype(dtype, copy=False)
        else:
            # invalid box, return zero vectors:
            box_matrix = np.zeros((3, 3), dtype=dtype)
    return box_matrix


