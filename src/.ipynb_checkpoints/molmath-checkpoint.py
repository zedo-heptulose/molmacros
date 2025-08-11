import numpy as np

import structure_math as smt

try:
    import plotly.graph_objects as go
except:
    print('could not import plotly')

def debug_print_vectors(*args):
    data = []
    colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black']  # Add more colors if needed
    max_range = [0, 0,0]  # Initialize max range

    for i, arg in enumerate(args):
        try:
            label, vector = arg
        except ValueError:
            print('invalid arg:', arg)
            continue
        
        if all(np.issubdtype(type(v), np.number) for v in vector):
            color = colors[i % len(colors)]
            if len(vector) == 2:
                data.append(go.Scatter3d(x=[0, vector[0]], y=[0, vector[1]], z=[0, 0], mode='lines', line=dict(color=color), name=label))
                data.append(go.Scatter3d(x=[vector[0]], y=[vector[1]], z=[0], mode='markers+text', text=[label], textposition="bottom center", marker=dict(color=color)))
                max_range = [max(max_range[0], abs(vector[0])), max(max_range[1], abs(vector[1])),max_range[2]]
            elif len(vector) == 3:
                data.append(go.Scatter3d(x=[0, vector[0]], y=[0, vector[1]], z=[0, vector[2]], mode='lines', line=dict(color=color), name=label))
                data.append(go.Scatter3d(x=[vector[0]], y=[vector[1]], z=[vector[2]], mode='markers+text', text=[label], textposition="bottom center", marker=dict(color=color)))
                max_range = [max(max_range[0], abs(vector[0])), max(max_range[1], abs(vector[1])), max(max_range[2], abs(vector[2])), max(max_range[1], abs(vector[1])), max(max_range[2], abs(vector[2]))]
            else:
                print('vector dimension not supported:', vector)
        else:
            print('vector not a number:', vector)


    max_max = max(max_range)
    if data:
        fig = go.Figure(data=data)
        fig.update_layout(
            scene=dict(
                xaxis=dict(range=[-max_max, max_max]),
                yaxis=dict(range=[-max_max, max_max]),
                zaxis=dict(range=[-max_max, max_max]),
            ),
            showlegend=False
        )
        fig.show()
    else:
        print('no valid vectors found')


def align_matrix(_v1, _v2, **kwargs):
    '''
    ALIGNS v2 to v1 !!!!
    accepts:
    v1, v2: numpy arrays
    returns:
    align_matrix: numpy array
    '''
    debug = kwargs.get('debug', False)
    
    v1 = _v1 / np.linalg.norm(_v1)
    v2 = _v2 / np.linalg.norm(_v2)
    
    if np.allclose(v1,v2,atol=1e-2):
        if debug:
            print('vectors are the same')
            
        return np.eye(3)
    if np.allclose(v1,-v2,atol=1e-2):
        if debug:
            print('vectors are antiparallel')
        #copilot code check later
        #construct rotation about arbitrary axis
        z_axis = np.array([0,0,1])
        y_axis = np.array([0,1,0])
        rot_ax = np.cross(v1, z_axis)
        if np.allclose(rot_ax, np.array([0,0,0])):
            rot_ax = np.cross(v1, y_axis)
            if np.allclose(rot_ax, np.array([0,0,0])):
                raise ValueError("bad vector?")
        return rotation_about_axis(rot_ax, np.pi)
    
    if debug:   
        print('arbitrary angle')
    
    axis = kwargs.get('axis', None)
    if axis is None:
        axis = np.cross(v1, v2)
        
    if debug:
        tuple1 = ('v1', v1)
        tuple2 = ('v2', v2)
        tuple3 = ('axis', axis)
        
    
    angle_between = angle_about_axis(v1, v2, axis)
    #HERE behaves unpredictably. returns align where it should antialign, etc.
    par_angle = - angle_between
    ##print (axis)
    par_matrix = rotation_about_axis(axis, par_angle)
    
    #print ('gets to end')
    return par_matrix



def angle_about_axis(_v1, _v2, _axis, **kwargs):
    '''
    accepts:
    v1, v2, axis: numpy arrays
    returns:
    angle: float
    '''
    debug = kwargs.get('debug', False)
    
    v1 = _v1 / np.linalg.norm(_v1)
    v2 = _v2 / np.linalg.norm(_v2)
    axis = _axis / np.linalg.norm(_axis)
    
    if np.allclose(v1, v2,atol=1e-2):
        return 0
    if np.allclose(v1, -v2,atol=1e-2):
        return np.pi
    
    dot = np.dot(v1, v2)
    cross = np.cross(v1, v2)
    c_n = cross / np.linalg.norm(cross)
    
    if debug:
        debug_print_vectors(('v1', v1), ('v2', v2), ('axis', axis), ('cross', cross), ('c_n', c_n))
    
    if np.allclose(c_n, axis,atol = 1e-1):
        return np.arccos(dot)
    elif np.allclose(c_n, -axis,atol = 1e-1):
        return 2*np.pi - np.arccos(dot)
    else:
        raise ValueError('Axis not parallel to cross product of v1 and v2')


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

