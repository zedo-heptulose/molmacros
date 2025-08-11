import numpy as np

reflection_xy = np.array([[1,0,0],
                          [0,1,0],
                          [0,0,-1]])

reflection_yz = np.array([[-1,0,0],
                          [0,1,0],
                          [0,0,1]])

reflection_xz = np.array([[1,0,0],
                          [0,-1,0],
                          [0,0,1]])

def close(val1, val2):
    '''
    USE THIS INSTEAD OF == WITH float VALUES
    '''
    tolerance = 1e-4
    return abs(val1 - val2) < tolerance

def rotation_x(angle):
    """
    chatgpt code
    Returns the rotation matrix for rotation about the x-axis by the given angle (in radians).
    """
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    return np.array([[1, 0, 0],
                     [0, cos_theta, -sin_theta],
                     [0, sin_theta, cos_theta]])

def rotation_y(angle):
    """
    chatgpt code
    Returns the rotation matrix for rotation about the y-axis by the given angle (in radians).
    """
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    return np.array([[cos_theta, 0, sin_theta],
                     [0, 1, 0],
                     [-sin_theta, 0, cos_theta]])

def rotation_z(angle):
    """
    chatgpt code
    Returns the rotation matrix for rotation about the z-axis by the given angle (in radians).
    """
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    return np.array([[cos_theta, -sin_theta, 0],
                     [sin_theta, cos_theta, 0],
                     [0, 0, 1]])

#need to unit test all of these functions

# def angle_between_vectors(vector_1, vector_2):
#     '''
#     this doesn't work on principle because it doesn't know directions
#     '''
#     v1 = vector_1.copy()
#     v2 = vector_2.copy()
    
#     if np.allclose(v1, [0,0,0]) or np.allclose(v2, [0,0,0]):
#         raise ValueError("zero vector in angle_between_vectors")
    
#     if np.linalg.norm(v1) == 0:
#         raise ValueError("floating point failure")
#     if np.linalg.norm(v1) == 0:
#         raise ValueError("floating point failure")
    
#     v1 = v1 / np.linalg.norm(v1)
#     v2 = v2 / np.linalg.norm(v2)
    
#     if np.allclose(v1, [0,0,0]) or np.allclose(v2, [0,0,0]):
#         raise ValueError("zero vector in angle_between_vectors")
    
#     #THESE ARE NECESSARY BECAUSE OF FLOATING POINT IMPRECISION
#     if np.allclose(v1, v2):
#         return 0
#     if np.allclose(v1, -v2):
#         return np.pi
    
    
#     angle = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
#     #print (angle)
    
#     return angle


import numpy as np

def angle_between_vectors(_v1, _v2):
    '''
    copilot wrote this
    check floating point errors around 
    theta = 0, pi
    Calculate the angle between two vectors in 3D space, taking into account the direction.
    '''
    # Copy vectors
    v1 = _v1.copy()
    v2 = _v2.copy()

    # Check for zero vectors
    if np.allclose(v1, [0,0,0]) or np.allclose(v2, [0,0,0]):
        raise ValueError("zero vector in angle_between_vectors")

    # Normalize vectors
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    # Calculate the cross and dot products
    cross_product = np.cross(v1, v2)
    dot_product = np.dot(v1, v2)

    # Calculate the angle
    angle = np.arctan2(np.linalg.norm(cross_product), dot_product)

    return angle

#new set of dihedral functions: just use vectors, not points.

def calculate_dihedral(p0, p1, p2, p3):
    """
    Calculate the dihedral angle between the planes formed by the points p0, p1, p2, and p3.

    Args:
    p0, p1, p2, p3: numpy arrays representing the coordinates of the points.

    Returns:
    The dihedral angle in radians.
    """
    # Calculate the vectors between the points
    v0 = p1 - p0
    v1 = p2 - p1
    v2 = p3 - p2
    # Calculate the normals to the planes formed by the points
    n1 = np.cross(v0, v1)
    n2 = np.cross(v1, v2)
    # Calculate the angle between the normals 
    return angle_between_vectors(n1, n2)

def rotation_about_axis(_axis, angle):
    """
    unit tests.... these linear transformations would be a great thing to learn and apply unit testing for, no?
    
    Calculate the rotation matrix from axis and angle using Rodrigues' rotation formula.

    Args:
    axis: The rotation axis as a numpy array.
    angle: The rotation angle in radians.

    Returns:
    The rotation matrix.
    """
    # Normalize the axis
    ###print(_axis)
    
    axis = _axis.copy()
    if np.allclose(axis, np.array([0,0,0])):
        raise ValueError("zero vector forbidden")
    
    axis /= np.linalg.norm(axis)
    
    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)  # Normalize the axis vector
    
    # Rodrigues' rotation formula
    cos_theta = np.cos(angle)
    sin_theta = np.sin(angle)
    cross_matrix = np.array([[0, -axis[2], axis[1]],
                             [axis[2], 0, -axis[0]],
                             [-axis[1], axis[0], 0]])
    rotation_matrix = np.eye(3) * cos_theta + sin_theta * cross_matrix + (1 - cos_theta) * np.outer(axis, axis)
    
    return rotation_matrix

def align_matrix(v1, v2):
    raise NotImplementedError()

def antialign_matrix(v1, v2, **kwargs):
    raise NotImplementedError()


def dihedral_rotation_matrix(p1, p2, dihedral_angle):
    """
    Rotate points p2 and p3 to have an arbitrary dihedral angle with points p0 and p1.

    Args:
    p0, p1, p2, p3: numpy arrays representing the coordinates of the points.
    desired_dihedral: The desired dihedral angle in radians.

    Returns:
    The rotated coordinates of points p2 and p3.
    """
    # Calculate the axis of rotation
    axis = p2 - p1

    # Normalize the axis
    axis /= np.linalg.norm(axis)

    # Calculate the rotation matrix
    rotation_matrix = rotation_about_axis(axis, dihedral_angle)

    return rotation_matrix
    

def transform(molecule, matrix): #pass atom_coords type as atom_coord_list
    '''
    apply matrix multiplication to each atomic coordinate in molecule
    '''
    
    new_molecule = {atom : matrix @ molecule[atom] for
                    atom in molecule}
    return new_molecule


def translate(molecule, vector):
    '''
    do vector addition for each coord
    '''
    new_molecule = {atom : molecule[atom] + vector for
                    atom in molecule}
    return new_molecule
    
