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



