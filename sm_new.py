import numpy as np
import molmath as mm
import structure_math as smt
import structuremaker as sm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import re
from scipy.spatial.distance import pdist, squareform


color_map = {
    'C': 'black',
    'H': 'white',
    'O': 'red',
    'N': 'blue',
    'F': 'lime',
    'Cl': 'green',
    'Br': 'brown',
    'I': 'purple',
    'S': 'yellow',
    'P': 'orange',
    'B': 'pink',
    'R': 'grey'
}


def show(mol_coords, **kwargs):
    '''
    shows the molecule
    '''
    show_hydrogens = kwargs.get('show_hydrogens', True)
    if not show_hydrogens:
        mol_coords = {key: mol_coords[key] for key in mol_coords if not key.startswith('H')}
    x = [mol_coords[atom][0] for atom in mol_coords]
    y = [mol_coords[atom][1] for atom in mol_coords]
    z = [mol_coords[atom][2] for atom in mol_coords]
    atom_types = [re.match(r'[A-Za-z]+', key).group() for key in mol_coords]

    colors = [color_map[atom_type] for atom_type in atom_types]

    atom_labels = [key for key in mol_coords]
    
    data = [
        go.Scatter3d(
            x=x, 
            y=y, 
            z=z, 
            mode='markers',
            text=atom_labels,
            hoverinfo='text',
            marker=dict(
                size=3,
                color=colors,                # set color to an array/list of desired values
                colorscale='Viridis',   # choose a colorscale
                opacity=0.8
            )
        )
    ]
    
    coordinates = list(zip(x,y,z))
    distances = squareform(pdist(coordinates))
    
    for i in range(len(distances)):
        for j in range(i+1, len(distances)):
            if distances[i][j] < 1.6 and not (atom_types[i] == 'H' and atom_types[j] == 'H'):
                data.append(
                    go.Scatter3d(
                        x=[coordinates[i][0], coordinates[j][0]],
                        y=[coordinates[i][1], coordinates[j][1]],
                        z=[coordinates[i][2], coordinates[j][2]],
                        mode='lines',
                        line=dict(
                            color='black', 
                            width=3
                        )
                    )
                )
    fig = go.Figure(data=data)
    
    min_coord = min(min(x),min(y),min(z))
    max_coord = max(max(x),max(y),max(z))
    
    
    show_geom = kwargs.get('show_geom', False)
    fig.update_layout(
        showlegend=False,
        scene = dict(
            xaxis=dict(showgrid=show_geom, showline=show_geom, showbackground =show_geom, range=[min_coord, max_coord]),
            yaxis=dict(showgrid=show_geom, showline=show_geom, showbackground =show_geom, range=[min_coord, max_coord]),
            zaxis=dict(showgrid=show_geom, showline=show_geom, showbackground =show_geom, range=[min_coord, max_coord]),

        ),
        width = 500,
        height= 500,
        
        autosize=False,
        margin=dict(
            l=50,  # left margin
            r=50,  # right margin
            b=100,  # bottom margin
            t=100,  # top margin
            pad=10  # padding
                ),
        paper_bgcolor='rgba(140,140,140,0.5)',
        plot_bgcolor='rgba(0,100,140,0.5)',
    ),
    fig.show()
    
def replace_r_atoms(molecule):
    '''
    '''
    new_molecule = molecule.copy()
    assert type(new_molecule) is type(molecule)
    new_atoms = []
    
    for atom in molecule:
        symbol = re.match(r'[A-Za-z]+', atom).group()
        if symbol == 'R':
            del new_molecule[atom]
            new_atoms.append(('H', molecule[atom]))
            
    for atom in new_atoms:
        new_molecule.add_atom(atom[0], atom[1])
        
    del molecule
    return new_molecule
            
    
# def debug_print_vectors(*args):
#     for arg in args:
#         label, vector = arg
#         if np.issubdtype(type(vector[0]), np.number) and np.issubdtype(type(vector[1]), np.number):
#             plt.quiver(0, 0, vector[0], vector[1], angles='xy', scale_units='xy', scale=1, label=label)
#         else:
#             print('vector not a number:', vector)
    
#     plt.xlim(-10, 10)
#     plt.ylim(-10, 10)
#     plt.grid()
#     plt.legend()
#     plt.title('Debug Vectors')
#     plt.show()
    

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


def add_group(_m1, keys1, _m2, keys2, **kwargs):
    '''
    '''
    debug = kwargs.get('debug', False)
    bond_length = kwargs.get('bond_length', 1.2)
    dihedral_angle = kwargs.get('dihedral_angle', 0)
    if debug:
        print('\n\n\n\n\n////////////////////in add_group////////////////////////\n\n\n\n\n')
        print('first molecule:')
        _m1.show()
        print('second molecule:')
        _m2.show()
        print('making copies of both molecules')
    
    m1 = _m1.copy()
    m2 = _m2.copy()
    
    m1_r = keys1[0]
    m1_c = keys1[1]
    m1_d = keys1[2] if len(keys1) > 2 else None
    
    m2_r = keys2[0]
    m2_c = keys2[1]
    m2_d = keys2[2] if len(keys2) > 2 else None
    
    m1_cr = m1[m1_r] - m1[m1_c]
    m2_cr = m2[m2_r] - m2[m2_c]
    
    if debug: 
        tuple1 = ('m1_cr:', m1_cr)
        tuple2 = ('m2_cr:', m2_cr)
        debug_print_vectors(tuple1, tuple2)
    
    #align bonds...
    
    m1 = m1.translate(-m1[m1_c])
    m2 = m2.translate(-m2[m2_c])
    
    assert np.allclose(m1[m1_c], np.array([0,0,0]))
    assert np.allclose(m2[m2_c], np.array([0,0,0]))
    
    anti_mat = mm.align_matrix(-m1_cr,m2_cr)
    
    m2 = m2.transform(anti_mat)
    
    if debug:
        print('molecule before dihedral transform:',end='\n\n')
        m3 = m1.union_with(m2)
        m3.show()
        del m3
    #only now that bonds are aligned...
    
    if m1_d and m2_d:
        m1_cr = m1[m1_r] - m1[m1_c]
        m2_cr = m2[m2_r] - m2[m2_c]
        
        m1_cd = m1[m1_d] - m1[m1_c]
        m2_cd = m2[m2_d] - m2[m2_c]

        
        if debug: 
            tuple1 = ('m1_cr:', m1_cr)
            tuple2 = ('m2_cr:', m2_cr)
            debug_print_vectors(tuple1, tuple2)
    

        m1_plane_norm = np.cross(m1_cr, m1_cd)
        m2_plane_norm = np.cross(m2_cr, m2_cd)
        
        
        if debug:
            tuple1 = ('m1_cd:', m1_cd)
            tuple2 = ('m2_cd:', m2_cd)
            tuple3 = ('m1_cr:', m1_cr)
            tuple4 = ('m2_cr:', m2_cr)
            tuple5 = ('m1_plane_norm:', m1_plane_norm)
            tuple6 = ('m2_plane_norm:', m2_plane_norm)
            debug_print_vectors(tuple1, tuple2,tuple3, tuple4, tuple5, tuple6)
        
        angle_between = mm.angle_about_axis(m1_plane_norm, m2_plane_norm, m1_cr,debug=debug)
        
        trans_angle = dihedral_angle - angle_between
        
        if debug:
            print('angle_between:', angle_between)
            print('diheral_angle:', dihedral_angle)
            print('trans_angle:', trans_angle  )
        #might just scrap and rewrite smt
        rotation_mat = smt.rotation_about_axis(m1_cr, trans_angle)
        
        #need to see if this goes the way I want it to..
        m2 = m2.transform(rotation_mat)
        if debug:
            print('molecule after dihedral transform:',end='\n\n')
            m3 = m1.union_with(m2)
            m3.show()
            del m3
            
    
    bond_vector = (m1[m1_r] - m1[m1_c])/np.linalg.norm(m1[m1_r] - m1[m1_c]) * bond_length
    
    
    if debug: 
        print('bond_vector:', bond_vector)
        print('length: ' + str(np.linalg.norm(bond_vector)))
        tuple1 = ('bond_vector:', bond_vector)
        debug_print_vectors(tuple1)
    
    m2 = m2.translate(bond_vector)
    
    del m1[m1_r]
    del m2[m2_r]
    
    m1 = m1.union_with(m2)
    
    if debug:
        print('molecule after bond translation and union:',end='\n\n')
        m1.show()
    
    return m1
    

def add_across_bond(_m1, keys1, _m2, keys2, **kwargs):
    '''
    for appending rings together along a shared bond,
    and the like.
    the keys for this one should be two C's (mandatory) and two R's (optional)
    '''
    # line up molecules along shared bond
    # find the shared bond
    debug = kwargs.get('debug', False)    
    if debug:
        print('\n\n\n\n\n////////////////////in add_across_bond////////////////////////\n\n\n\n\n')
        print('first molecule:')
        show(_m1)
        print('second molecule:')
        show(_m2)
        print('making copies of both molecules')

    m1 = _m1.copy()
    m2 = _m2.copy()
    
    if debug:
        debug_args = []
        for key in keys1:
            if m1.get(key,None) is None:
                print('key not found:', key)
            else:
                debug_args.append(tuple([f'mol1{key}', m1[key]]))
        for key in keys2:
            if m2.get(key,None) is None:
                print('key not found:', key)
            else:
                debug_args.append(tuple([f'mol2{key}', m2[key]]))
        debug_print_vectors(*debug_args)
    
    m1_c1 = keys1[0]
    m1_c2 = keys1[1]
    m1_r1 = keys1[2] if len(keys1) > 2  else None 
    m1_r2 = keys1[3] if len(keys1) >3 else None
        
    m2_c1 = keys2[0]
    m2_c2 = keys2[1]
    m2_r1 = keys2[2] if len(keys2) > 2 else None
    m2_r2 = keys2[3] if len(keys2) > 3 else None
    
    # find the shared bond
    m1_c1c2 = m1[m1_c2] - m1[m1_c1]
    m2_c1c2 = m2[m2_c2] - m2[m2_c1]
    if debug:
        tuple_1 = ('m1_c1c2:', m1_c1c2)
        tuple_2 = ('m2_c1c2:', m2_c1c2)
        print('showing the two bonds to align:')
        debug_print_vectors(tuple_1, tuple_2)
    
    #align mol2 to mol1
    align_mat = smt.align_matrix(m1_c1c2, m2_c1c2)
    m2 = smt.translate(m2, -m2[m2_c1])
    m2 = smt.transform(m2, align_mat)
    if debug:
        print('align_mat:', align_mat)
        print('m2 after align:')
        show(m2)
    
    #find cross product...
    cond1 = m1.get(m1_r1,None) is not None or m2.get(m1_r2,None) is not None
    cond2 = m2.get(m2_r1,None) is not None or m2.get(m2_r2,None) is not None
    
    ignore_dihedral = kwargs.get('ignore_dihedral', False)
    
    if (cond1 and cond2) and not ignore_dihedral:
        m1_r = m1_r1 if m1.get(m1_r1,None) is not None else m1_r2
        m2_r = m2_r1 if m2.get(m2_r1,None) is not None else m2_r2
        
        m1_c1 = keys1[0]
        m1_c1r = m1[m1_r] - m1[m1_c1]
        
        if debug:
            tuple1 = ('m1_c1r1:', m1_c1r)
            tuple2 = ('m1_c1:', m1[m1_c1])
            tuple3 = ('m1_r:', m1[m1_r])
            debug_print_vectors(tuple1, tuple2, tuple3)
        
        # m2_c2 = keys2[1]
        m2_c1r = m2[m2_r] - m2[m2_c1]
        m2_c1c2 = m2[m2_c2] - m2[m2_c1]
        
        
        if debug:
            tuple1 = ('m2_c1r:', m2_c1r)
            tuple2 = ('m2_c1:', m2[m2_c1])
            tuple3 = ('m2_r:', m2[m2_r])
            debug_print_vectors(tuple1, tuple2, tuple3)
            print(f'm2_c1r: {m2_c1r}')
            print(f'm1_c1r: {m1_c1r}')
            print(f'm2_c1: {m2[m2_c1]}')
            print(f'm1_c1: {m1[m1_c1]}')
            print(f'm2_r: {m2[m2_r]}')
            print(f'm1_r: {m1[m1_r]}')
            
        m1_plane_norm = np.cross(m1_c1c2, m1_c1r)
        m2_plane_norm = np.cross(m2_c1c2, m2_c1r)
        
        if debug:
            tuple1 = ('m1_plane_norm:', m1_plane_norm)
            tuple2 = ('m2_plane_norm:', m2_plane_norm)
            debug_print_vectors(tuple1, tuple2)

        anti_align_mat = smt.antialign_matrix(m1_plane_norm, m2_plane_norm, axis = m1_c1c2)
        #need to be able to set axis for anti_align_mat
        m2 = smt.transform(m2, anti_align_mat)
        
        if debug:
            print('m2 after anti_align:')
            show(m2)
            print('anti_align_mat:', anti_align_mat)
    
    m2 = smt.translate(m2, m1[m1_c1] - m2[m2_c1])
    
    for key in keys2[:2]:
        H = find_close_H(m2, key)
        if H is not None:
            del m2[H]
        if key is not None and m2.get(key,None) is not None:
            if debug: print('deleting key:', key)
            del m2[key]
    
    for key in keys1[:2]:
        H = find_close_H(m1, key)
        if H is not None:
            del m1[H]
    
    if debug:
        print('m2 after deleting keys:')
        show(m2)
        print('making union:')
    
    
    m1 = sm.make_molecule_union(m1, m2)
    if debug:
        print('m1 after union:')
        show(m1)
    
    if debug: print('returning m1')
    return m1
    # if third atoms are present, antialign dihedrals
    
def find_close_H(molecule, atom, threshold = 1.3, **kwargs):
    '''
    find the closest H to an atom in the molecule
    '''
    debug = kwargs.get('debug', False)
    if debug:
        print('\n\n\n\n\n////////////////////in find_close_H////////////////////////\n\n\n\n\n')
        print('molecule:')
        show(molecule)
        print('atom:', atom)
    
    min_dist = np.inf
    closest_H = None
    for key in molecule:
        if key.startswith('H'):
            dist = np.linalg.norm(molecule[key] - molecule[atom])
            if dist < min_dist:
                min_dist = dist
                closest_H = key
    if debug:
        print('closest_H:', closest_H)
    if min_dist > threshold:
        return None
    return closest_H
    
    
def distort(molecule, function, **kwargs):
    '''
    apply a distortion function to the molecule
    '''
    debug = kwargs.get('debug', False)
    exempt = kwargs.get('exempt', [])
    if debug:
        print('\n\n\n\n\n////////////////////in distort////////////////////////\n\n\n\n\n')
        print('molecule:')
        show(molecule)
    
    for key in molecule:
        if key not in exempt:
            molecule[key] = function(molecule[key])
    
    if debug:
        print('molecule after distortion:')
        show(molecule)
    return molecule


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a grid of points in the x-y plane
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
x, y = np.meshgrid(x, y)

# Compute the distance from the origin for each point
r = np.sqrt(x**2 + y**2)

# Define the radial function z = r^2
z = r**2

# Create a 3D plot
fig = plt.figure();
ax = fig.add_subplot(111, projection='3d');
ax.plot_surface(x, y, z, cmap='viridis')

plt.show()