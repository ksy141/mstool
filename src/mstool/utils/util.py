import numpy as np

def fibo(r=None, N=None, APL=None, alpha=0, beta=0, gamma=0, 
        verbose=True, plot='plot.pdf'):

    ''' Make a Fibonacci sphere.
    Parameters
    ----------
    r : float
        radius in A
    N : int
        the number of points. Provide either N or APL.
    APL : float
        area per lipid in A^2. Provide either N or APL.
    alpha : float
        angle in degree for yaw. The default is 0.0.
    beta : float
        angle in degree for pitch. The default is 0.0.
    gamma : float
        angle in degree for roll. The default is 0.0.
    verbose : bool
        make it verbose. The default is True.
    plot : str
        plot the points. Provide the filename. The default is plot.pdf

    Returns
    -------
    points : array
        a numpy array of ``(N, 3)``. ``N`` is the provided N or 
        the calculated N (if you provide APL). 
    '''
    
    alpha *= np.pi / 180
    beta  *= np.pi / 180
    gamma *= np.pi / 180


    if r is not None and N is not None:
        APL = 4 * np.pi * r**2 / N
    elif r is not None and APL is not None:
        N = 4 * np.pi * r**2 / APL
    elif N is not None and APL is not None:
        r = np.sqrt(APL * N / 4 / np.pi)
    else:
        raise ValueError('sepcify two out of three args: r, N, APL')

    N = int(N)
    if verbose: print('r: %8.3f A\nN: %8d\nAPL: %6.3f A^2\n\n' %(r, N, APL))

    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(N):

        if N != 1:
            y = r * (1 - (i / (N - 1)) * 2)  # y goes from r to -r
            radius = np.sqrt(r ** 2 - y ** 2)  # radius at y
            theta = phi * i  # golden angle increment
            x = np.cos(theta) * radius
            z = np.sin(theta) * radius

        else:
            x = 0.0
            y = 0.0
            z = r

        if not (alpha == 0 and beta == 0 and gamma == 0):
            t00 = np.cos(beta)  * np.cos(gamma)
            t01 = np.sin(alpha) * np.sin(beta) * np.cos(gamma) - np.cos(alpha) * np.sin(gamma)
            t02 = np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma)
            t10 = np.cos(beta)  * np.sin(gamma)
            t11 = np.sin(alpha) * np.sin(beta) * np.sin(gamma) + np.cos(alpha) * np.cos(gamma)
            t12 = np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma)
            t20 = -np.sin(beta)
            t21 = np.sin(alpha) * np.cos(beta)
            t22 = np.cos(alpha) * np.cos(beta)
            RM3D = np.array([[t00, t01, t02], [t10, t11, t12], [t20, t21, t22]])
            points.append(RM3D @ np.array([x, y, z]))

        else:
            points.append([x, y, z])

    points = np.array(points)

    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        # aspect ratio is 1:1:1 in data space
        ax.set_box_aspect((np.ptp(points[:,0]), np.ptp(points[:,1]), np.ptp(points[:,2])))
        ax.scatter(points[:,0], points[:,1], points[:,2])
        fig.savefig(plot)

    return points


def align_a_b(a, b):
    ### a, b are unit vector
    ### return a rotation matrix that aligns a into b.
    # u.atoms.positions = np.matmul(R, u.atoms.positions.T).T + shift

    a = np.array(a) / np.linalg.norm(a)
    b = np.array(b) / np.linalg.norm(b)

    if np.all(a == b):
        R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        return R

    elif np.all(a == -b):
        R = -np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        return R

    else:
        v = np.cross(a, b)
        s = np.sqrt(np.sum(v**2))
        c = np.sum(a * b)

        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R  = np.identity(3) + vx + np.matmul(vx, vx) * (1 - c) / s**2
        return R


def completeTetra(v1,v2):
    # from CG2AA
    # https://www.ic.fcen.uba.ar/cg2aa/cg2aa.py
    # finds v3,v4 the other vertices of the tetra centered at zero
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    nv = np.cross(v1,v2)
    nv = 0.8164965809277259 * nv / np.linalg.norm(nv)
    v = (-v1-v2)/2
    w1 = v + nv
    w2 = v - nv
    w1 = w1 / np.linalg.norm(w1)
    w2 = w2 / np.linalg.norm(w2)
    return	w1, w2

