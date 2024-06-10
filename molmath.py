import numpy as np

def rotation_about_axis(axis, angle):
    """
    Calculate the rotation matrix from axis and angle.

    Args:
    axis: The rotation axis as a numpy array.
    angle: The rotation angle in radians.

    Returns:
    The rotation matrix.
    """
    # Normalize the axis
    axis = axis / np.linalg.norm(axis)

    # If the angle is zero (or very close), the vectors are already parallel
    if np.isclose(angle, 0):
        return np.eye(3)

    # If the angle is pi (or very close), the vectors are antiparallel
    if np.isclose(angle, np.pi):
        return -np.eye(3)

    # Calculate the quaternion components
    w = np.cos(angle / 2)
    x, y, z = axis * np.sin(angle / 2)

    # Calculate the rotation matrix
    rotation_matrix = np.array([
        [1 - 2*y*y - 2*z*z, 2*x*y - 2*z*w, 2*x*z + 2*y*w],
        [2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z, 2*y*z - 2*x*w],
        [2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - 2*x*x - 2*y*y]
    ])

    return rotation_matrix


def align_matrix(v1, v2, axis=None):
    """
    Returns a matrix that will align v1 to v2.
    """
    # Normalize the input vectors
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    # If axis is not provided, calculate the cross product of the vectors
    if axis is None:
        axis = np.cross(v1, v2)

    # If the cross product is a zero vector, the vectors are parallel or antiparallel
    if np.allclose(axis, np.array([0,0,0])):
        # If the vectors are parallel, return the identity matrix
        if np.allclose(v1, v2):
            return np.eye(3)
        # If the vectors are antiparallel, return a 180 degree rotation matrix around an arbitrary axis
        else:
            return rotation_about_axis(np.array([1,0,0]), np.pi)

    # Calculate the dot product of the vectors
    dot = np.dot(v1, v2)

    # Calculate the elements of the skew-symmetric cross-product matrix of v1
    skew = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])

    # Calculate the rotation matrix using the formula R = I + skew + skew^2/(1-dot)
    rotation_matrix = np.eye(3) + skew + np.dot(skew, skew)/(1-dot)

    return rotation_matrix

def antialign_matrix(v1, v2, axis=None):
    """
    Returns a matrix that will align v1 antiparallel to v2.
    """
    matrix = align_matrix(v1, v2,axis)
    matrix *= -1
    return matrix