from .universe import Universe
from ..lib.align import rotation_matrix
import numpy as np

def calInertialTensor(pos):
    Ixx = np.sum(pos[:,1]**2 + pos[:,2]**2)
    Iyy = np.sum(pos[:,0]**2 + pos[:,2]**2)
    Izz = np.sum(pos[:,0]**2 + pos[:,1]**2)

    Ixy = -np.sum(pos[:,0] * pos[:,1])
    Ixz = -np.sum(pos[:,0] * pos[:,2])
    Iyz = -np.sum(pos[:,1] * pos[:,2])

    T = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    w, v = np.linalg.eig(T)

    return w, v

def align_a_to_b(a, b):
    ''' 
    a, b are unit vector
    return a rotation matrix that aligns a into b
    
    R = align_a_b(a, b)
    xyz = np.matmul(R, pos.T).T
    '''
    a /= np.linalg.norm(a)
    b /= np.linalg.norm(b)

    v = np.cross(a, b)
    s = np.sqrt(np.sum(v**2))
    c = np.sum(a * b)
    
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R  = np.identity(3) + vx + np.matmul(vx, vx) * (1 - c) / s**2
    return R

def Orient(structure, out=None, select='@CA,BB', translate=[0,0,0], axis='z', rotation=True):
    if isinstance(structure, Universe):
        u  = structure
    else:
        u  = Universe(structure)
    
    allpos = u.atoms[['x','y','z']].to_numpy()
    subpos = u.select(select)[['x','y','z']].to_numpy()

    if axis == 'x':
        axis = [1, 0, 0]
    elif axis == 'y':
        axis = [0, 1, 0]
    elif axis == 'z':
        axis = [0, 0, 1]
    else:
        assert len(axis) == 3, 'provide a vector'
        axis = axis
    
    if rotation:
        w, v = calInertialTensor(subpos - subpos.mean(axis=0))
        R = align_a_to_b(v[:,0], axis)
        newpos = np.matmul(R, np.transpose(allpos - subpos.mean(axis=0))).T
        u.atoms[['x','y','z']] = newpos
    
    u.atoms[['x','y','z']] -= u.select(select)[['x','y','z']].mean(axis=0)
    u.atoms[['x','y','z']] += translate

    if out: u.write(out)
    return u

