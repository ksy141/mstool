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


def align_principal_component(points, subset_points, target_vector):
    """
    Rotate a 3D object so that its principal component aligns with a specified vector,
    calculated from a subset of the object's coordinates.

    Parameters:
        points (numpy.ndarray): Nx3 array of 3D coordinates representing the object.
        subset_points (numpy.ndarray): Mx3 array of 3D coordinates of selected atoms.
        target_vector (numpy.ndarray): 3D vector to which the principal component will be aligned.

    Returns:
        numpy.ndarray: Nx3 array of rotated 3D coordinates.
    """
    # Step 1: Calculate the principal component from the subset of points
    centroid = np.mean(subset_points, axis=0)
    centered_points = subset_points - centroid
    _, _, vh = np.linalg.svd(centered_points)
    principal_component = vh[0]

    # Step 2: Determine the rotation needed to align the principal component with the target vector
    rotation_axis = np.cross(principal_component, target_vector)
    rotation_axis /= np.linalg.norm(rotation_axis)
    dot_product = np.dot(principal_component, target_vector)
    rotation_angle = np.arccos(dot_product)

    # Step 3: Apply the rotation to the object's coordinates
    rotation_matrix = rotation_matrix_from_axis_angle(rotation_axis, rotation_angle)
    rotated_points = np.dot(points, rotation_matrix.T) + centroid

    return rotated_points

def rotation_matrix_from_axis_angle(axis, angle):
    """
    Generate a rotation matrix from an axis and an angle.

    Parameters:
        axis (numpy.ndarray): 3D vector representing the rotation axis.
        angle (float): Angle of rotation in radians.

    Returns:
        numpy.ndarray: 3x3 rotation matrix.
    """
    axis = axis / np.linalg.norm(axis)
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    rotation_matrix = np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                                [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                                [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])
    return rotation_matrix


def Orient(structure, out=None, select='@CA,BB', translate=[0,0,0], axis='z', rotation=True):
    if isinstance(structure, Universe):
        u  = structure
    else:
        u  = Universe(structure)
    
    print(u.select(select))
    allpos = u.atoms[['x','y','z']].to_numpy()
    subpos = u.select(select)[['x','y','z']].to_numpy()

    if axis == 'x':
        axis = [1, 0, 0]
    elif axis == 'y':
        axis = [0, 1, 0]
    elif axis == 'z':
        axis = [0, 0, 1]
    else:
        axis = axis
    
    if rotation:
        newpos = align_principal_component(allpos, subpos, axis/np.linalg.norm(axis))
        u.atoms[['x','y','z']] = newpos
    u.atoms[['x','y','z']] -= u.select(select)[['x','y','z']].mean(axis=0)
    u.atoms[['x','y','z']] += translate

    if out: u.write(out)
    return u

