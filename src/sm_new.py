import numpy as np
import molmath as mm
import structure_math as smt
import structuremaker as sm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import re
from scipy.spatial.distance import pdist, squareform

import io
import py3Dmol

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

def show_py3Dmol(mol_coords,**kwargs):
    xyz = to_xyz_string(mol_coords)
    view = py3Dmol.view(width=400,height=400)
    view.addModel(xyz,'xyz')
    view.setStyle({'sphere':{'radius' : 0.3}, 'stick' : {'radius': 0.2}})
    for element, color in color_map.items():
        view.setStyle({'element': element}, {'sphere': {'color': color, 'radius': 0.3}})
    excluded_symbols = ['Bq', 'DA']
    for symbol in excluded_symbols:
        view.setStyle({'elem': symbol}, 
                      {'sphere': {'colorscheme': 'element', 'radius': 0.3}, 'stick': {'hidden': True}})
    view.zoomTo()
    view.show()

def show(mol_coords, unitcell=None, **kwargs):
    '''
    Shows a molecule or multiple molecules along with an optional unit cell.
    
    Parameters:
    - mol_coords: Either a single molecule's coordinate dictionary or a list of such dictionaries.
    - lattice_vectors: A 3x3 numpy array representing the unit cell lattice vectors (optional).
    - show_hydrogens: Whether to display hydrogen atoms (default True).
    - show_geom: Whether to show geometry lines (default False).
    '''
    show_hydrogens = kwargs.get('show_hydrogens', True)
    show_geom = kwargs.get('show_geom', False)
    data = []
    
    # Check if we're dealing with a single molecule or multiple molecules
    if isinstance(mol_coords, dict):
        # Convert single molecule into a list
        mol_coords_list = [mol_coords]
    elif isinstance(mol_coords, list):
        mol_coords_list = mol_coords
    else:
        raise ValueError("mol_coords must be a dictionary (for a single molecule) or a list of dictionaries (for multiple molecules).")
    
    # Loop through the list of molecule coordinates
    for mol_index, mol_coords in enumerate(mol_coords_list):
        if not show_hydrogens:
            mol_coords = {key: mol_coords[key] for key in mol_coords if not key.startswith('H')}

        x = [mol_coords[atom][0] for atom in mol_coords]
        y = [mol_coords[atom][1] for atom in mol_coords]
        z = [mol_coords[atom][2] for atom in mol_coords]
        atom_types = [re.match(r'[A-Za-z]+', key).group() for key in mol_coords]

        # Assign different colors to each molecule
        colors = [color_map[atom_type] for atom_type in atom_types]

        atom_labels = [key for key in mol_coords]
        
        # Plot atoms as scatter points
        data.append(
            go.Scatter3d(
                x=x, 
                y=y, 
                z=z, 
                mode='markers' + ('+text' if kwargs.get('show_labels',False) else ''),
                text=[f'{label}' for label in atom_labels], 
                hoverinfo='text',
                marker=dict(
                    size=3,
                    color=colors,   # Atom colors based on their types
                    colorscale='Viridis',
                    opacity=0.8
                )
            )
        )
        
        # Calculate distances between atoms and plot bonds
        coordinates = list(zip(x, y, z))
        distances = squareform(pdist(coordinates))
        
        for i in range(len(distances)):
            for j in range(i+1, len(distances)):
                if distances[i][j] < 1.8 and not (atom_types[i] == 'H' and atom_types[j] == 'H'):
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

    if unitcell:
        print('USING UNIT CELL:')
        #display vertices connected by lines
        origin = np.array([0,0,0])
        corners = [
            origin,
            unitcell.tv1,
            unitcell.tv2,
            unitcell.tv3,
            unitcell.tv1 + unitcell.tv2,
            unitcell.tv1 + unitcell.tv3,
            unitcell.tv2 + unitcell.tv3,
            unitcell.tv1 + unitcell.tv2 + unitcell.tv3,
        ]

        edges = [
            (0, 1), (0, 2), (0, 3),
            (1, 4), (1, 5), (2, 4), (2, 6),
            (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)
        ]

        for start, end in edges:
            data.append(
                go.Scatter3d(
                    x=[corners[start][0], corners[end][0]],
                    y=[corners[start][1], corners[end][1]],
                    z=[corners[start][2], corners[end][2]],
                    mode='lines',
                    line=dict(color='blue', width=2),
                    name='Unit Cell Edge'
                )
            )

    # Set up the figure layout
    fig = go.Figure(data=data)
    
    ax_style = dict(showbackground=False,
                    showgrid=True,
                    zeroline=True,
                    gridcolor='red',
                    zerolinecolor='green',
                    )

    fig.update_layout(
        scene=dict(xaxis=ax_style,
                   yaxis=ax_style,
                   zaxis=ax_style,
                   ),

        showlegend=False,
        width=500,
        height=500,
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
        scene_aspectmode='data',
    )
    
    show = kwargs.get('show', True)
    if show:
        fig.show()
    else:
        return fig

    # fig.show()






def replace_r_atoms(molecule, **kwargs):
    '''
    '''
    print ('AAAA')
    debug = kwargs.get('debug', False)
    
    new_molecule = molecule.copy()
    assert type(new_molecule) is type(molecule)
    coords = []
    
    if debug:
        print('molecule before replacing R atoms:')
        new_molecule.show()
    
    for atom in molecule: #don't really get it.
        if atom.startswith('R'):
            del new_molecule[atom] # this also works.
            if debug:
                print('deleting atom:', atom)
                new_molecule.show() 
            assert type(molecule[atom]) is np.ndarray
            coords.append(molecule[atom]) # so should this?
            
    for coord in coords:
        new_molecule = new_molecule.add_atom('H', coord) # I know this works.
        if debug: 
            print('adding atom:', 'H', coord)
            new_molecule.show()
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

    try:
        m1_cr = m1[m1_r] - m1[m1_c]
        m2_cr = m2[m2_r] - m2[m2_c]
    #TODO: fix tight coupling of sm_new and molecule
    except:
        print("Invalid keys for one of the atoms.")
        print(f"m1 keys:{m1_r} {m1_c} {m1_d}")
        print("m1 atom coords:")
        print(m1.atom_coords)
        print(f"m2 keys:{m2_r} {m2_c} {m2_d}")
        print("m2 atom coords:")
        print(m2.atom_coords)
    
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
    
    try:
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
    
    except:        
        print("FINDING DIHEDRAL ANGLE FAILED")
    
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
    m1_r1 = keys1[2] if len(keys1) > 2 else None 
    m1_r2 = keys1[3] if len(keys1) > 3 else None
        
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
    align_mat = mm.align_matrix(m1_c1c2, m2_c1c2)
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
        
        m1_c1r = m1[m1_r] - m1[m1_c1]

        new_m2_c1c2 = m2[m2_c2] - m2[m2_c1]
        m2_c1r =m2[m2_r] - m2[m2_c1]
        
        m1_plnorm = np.cross(m1_c1c2,m1_c1r)
        m1_plnorm /= np.linalg.norm(m1_plnorm)
        m2_plnorm = np.cross(new_m2_c1c2,m2_c1r)
        m2_plnorm /= np.linalg.norm(m1_plnorm)

        intersect_angle = kwargs.get('angle',np.pi)

        current_angle = mm.angle_about_axis(m1_plnorm,m2_plnorm,new_m2_c1c2)

        difference_angle = intersect_angle - current_angle

        rotation_matrix = smt.rotation_about_axis(new_m2_c1c2,difference_angle)

        #it is still properly translated to do this
        m2 = smt.transform(m2,rotation_matrix)
        
    
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
    
    m1 = make_molecule_union(m1, m2)
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
    atom_mask = kwargs.get('atom_mask',None)
    anti_mask = kwargs.get('anti_mask', [])
    
    if debug:
        print('\n\n\n\n\n////////////////////in distort////////////////////////\n\n\n\n\n')
        print('molecule:')
        show(molecule)
    
    if atom_mask:
        for key in atom_mask:
            molecule[key] = function(molecule[key])
    else:
        for key in molecule:
            if key not in anti_mask:
                molecule[key] = function(molecule[key])
    
    if debug:
        print('molecule after distortion:')
        show(molecule)
    return molecule


def make_molecule_union(molecule1, molecule2, **kwargs):
    mol1 = molecule1.copy()
    mol2 = molecule2.copy()

    molecule_union = {}
    atom_count = 1
    #count keys in mol1 to see what highest number is
    #count from there when adding keys from mol2
    #index starts at 1
    num_pattern = re.compile(r'\d+')    
    for key in molecule1:
        #look for highest index
        num = re.search(num_pattern,key)
        if num:
            num = int(num.group(0))
        else:
            raise ValueError('atom without number in molecule')
        if num > atom_count:
            atom_count = num
        #add atom
        coords = molecule1[key]
        molecule_union[key] = coords
    
    for key in molecule2:
        atom_count += 1
        symbol = re.match(r'([A-Za-z]{1,2})',key).group(1)
        new_key = symbol + str(atom_count)
        coords = molecule2[key]
        molecule_union[new_key] = coords

    return molecule_union
    
def sort_keys(_molecule):
    '''
    collapses numbers in keys so that there 
    are no gaps counting from 1
    '''
    molecule = {}
    keys = [key for key in _molecule]
    keys = sorted(keys,key=lambda x: int(re.search(r'\d+',x).group(0)))
    atom_count = 1
    for key in keys:
        symbol = re.match(r'([A-Za-z]{1,2})',key).group(1)
        new_key = symbol + str(atom_count)
        coords = _molecule[key]
        molecule[new_key] = coords
        atom_count += 1
    return molecule

def rotate_keys(_molecule,key1,position):
    '''
    Swaps the numbers associated with two keys.
    Expects a molecule with collapsed atom numbers.
    Unpredictable results if this is not so.
    '''
    key1_num = re.search(r'\d+',key1).group(0)
    molecule = _molecule.copy()
    coords1 = _molecule[key1]
    key2 = None
    for key in _molecule:
        num = int(re.search(r'\d+',key).group(0))
        if num == position:
            key2 = key
            coords2 = _molecule[key2]
            
            key2_symbol = re.match(r'([A-Za-z]{1,2})',key1).group(1)
            key2_num = re.search(r'\d+',key).group(0)
            key1_symbol = re.match(r'([A-Za-z]{1,2})',key1).group(1)
            new_key1 = key1_symbol + key2_num
            new_key2 = key2_symbol + key1_num
            break
    if not key2:
        raise ValueError('bad position number')
    
    del molecule[key1]
    del molecule[key2]
    molecule[new_key1] = coords1
    molecule[new_key2] = coords2

    return molecule

#writes coordinate output
def write_xyz(filename, atom_dict, comment = None):
    with open(filename, 'w') as coordinate_output:
        coordinate_output.write(to_xyz_string(atom_dict,comment))

def to_xyz_string(atom_dict, comment=None):
    coordinate_output = io.StringIO()
    coordinate_output.write(f'{str(len(atom_dict))}\n')
    if comment is not None:
        coordinate_output.write(f'{comment}\n')
    else:
        coordinate_output.write('\n')
    atom_list = [(key,atom_dict[key]) for key in atom_dict]
    atom_list = sorted(atom_list,key=lambda x: int(re.search(r'\d+',x[0]).group(0)))
    
    for atom, coords in atom_list:
        symbol = re.match(r'[A-Za-z]{1,2}', atom)[0]
        coordinate_output.write(f'{symbol} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}\n')
    
    contents = coordinate_output.getvalue()
    coordinate_output.close()

    return contents

