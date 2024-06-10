#want to have: lotta stuff. Honestly, big problem here. I don't know what I'm doing.
import re
import structure_math as smt
import numpy as np
import re
import matplotlib.pyplot as plt
import itertools
from itertools import combinations

#need to modify unitcell to support this
#might need connectivity table for larger molecules. maybe can use inferred for now
#....
#first, need to be able to read and write coordinates in an appropriate format
#Gaussian uses the same format as .xyz, only with charge/multiplicity before coordinates
#so xyz is an appropriate format to read to and write from
#the very first functions I need are coordinate readers and writers, because I can't do or test anything without i/o


#molecules are dicts of numbered atomic symbols : coordinates (np.array)
#'atom' and 'key' are used interchangably in these functions

def FormatError():
    print('bad .xyz format')


#reads coordinate output
def read_xyz(filename):
    with open(filename, 'r') as coordinate_file:
        lines = coordinate_file.readlines()
    
    match = re.match(r'(?:\s*)?(\d+)', lines[0])
    if match == False:
        raise ValueError()
    num_atoms = float(match[1])
       
    atom_dict = {}   
        
    atom_counter = 0
    r_counter = 0    
    for line in lines[2:]:
        if atom_counter >= num_atoms:
            break
        
        atom_symbol = re.search(r'[A-Za-z]{1,2}', line)[0]
        coord_match = re.findall(r'[-0-9.]+', line)
        
        coords = np.array([float(coord_match[0]), float(coord_match[1]), float(coord_match[2])])
        
        if (atom_symbol == False or coord_match == False):
            raise ValueError()       
        
        r_group = re.search(r'[R]{1,2}', line)
        if r_group:
            r_counter += 1
            atom_name = 'R'+ str(r_counter)
        else:
            atom_counter += 1
            atom_name = atom_symbol + str(atom_counter)
        
        atom_dict.update({atom_name : coords})
        
    return atom_dict, atom_counter, r_counter

#manipulations happen here
def read_gauss_out(out_filename):
    '''
    
    '''
    #read all the lines of the output file into list of strings
    with open(out_filename, 'r') as out:
        lines = out.readlines()

    #define REs for finding the start and end
    #section of outfile with coordinates
    start_pattern = r'(\|\d,\d\|)'
    end_pattern = r'(\|@)'
    #list to store lines w/ coordinate info
    coord_lines = []
    #flag to tell if in this desired section
    read = False
    #for each line:
    for line in lines:
        #if we find the start pattern, set read flag to true
        if re.search(start_pattern, line):
            #print('flag is true')
            read = True
        # if read flag is true, record the current line
        if read:
            #print(line)
            coord_lines.append(line)
        # if read flag is true and line with end pattern reached,
        # stop reading and break from for loop
        if re.search(end_pattern, line) and read == True:
            break
    #pattern that all coords follow: atom symbol (non cap ','), 
    coord_pattern = r'([A-Za-z]{1,2})(?:,)(-?\d+\.\d+)(?:,)(-?\d+\.\d+)(?:,)(-?\d+\.\d+)'
    #for each line from section with coordinate info, strip \n character
    # and append to very long string
    long_string = ' '
    for line in coord_lines: 
        #print(line)
        line = line.strip()
        long_string += line
    #split very long string into the strings delimited by '|'
    #print(long_string)
    info_strings = long_string.split('|')
    #make list to hold internal data read from coord string
    atom_dict = {}
    # for each 
    index = 1
    for info_string in info_strings:
        
        match = re.match(coord_pattern, info_string)
        if match:
            atom_symbol = match.group(1) + str(index)
            x = float(match.group(2))
            y = float(match.group(3))
            z = float(match.group(4))
            coords = np.array([x,y,z])
            atom_dict.update({atom_symbol : coords})
            index+= 1
        
    return atom_dict
            
            

#organize thoughts
#display 3d structure with bonds and labels
#^first function to code, makes everything after much easier
def show(molecule, **kwargs):
    '''
    displays the molecule in 3d.
    
    accepts atom dictionary string : np.array as an argument.
    also accepts flags, which are display hydrogens, label atoms
    respectively.
    
    plots all the points with their
    '''
    _label = kwargs.get('label', None)
    _flags = kwargs.get('flags', [True, True, True])
    _polygons_colors = kwargs.get('polygons_colors', [None, None])
    display_hydrogens = _flags[0]
    label_atoms = _flags[1]
    show_axes = _flags[2]
    
    bond_threshold = 1.7

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    keys = [key for key in molecule if not re.match('H\d+', key)]
    
    for key in keys:
        if not isinstance(molecule[key], (list, np.ndarray)):
            raise TypeError(f'Expected list or array for key {key}, got {type(molecule[key])}')

    x_points = [molecule[key][0] for key in keys]
    y_points = [molecule[key][1] for key in keys]
    z_points = [molecule[key][2] for key in keys]
    
    #if display hydrogens flag is true, appends h keys to keys
    #and H points to points
    if display_hydrogens:
        H_keys = [key for key in molecule if re.match('H\d+', key)]
        
        H_x_points = [molecule[key][0] for key in H_keys]
        H_y_points = [molecule[key][1] for key in H_keys]
        H_z_points = [molecule[key][2] for key in H_keys]
        
        keys = keys + H_keys
        x_points = x_points + H_x_points
        y_points = y_points + H_y_points
        z_points = z_points + H_z_points
        
    ax.scatter(x_points, y_points, z_points)
        
    #plots atoms at their coordinates, if label atoms flag is true
    if(label_atoms):
        for i, txt in enumerate(keys):
            ax.text(x_points[i], y_points[i], z_points[i], txt)
    
    
    #draws bonds
    for atom1, atom2 in combinations(molecule.keys(), 2):
        distance = np.linalg.norm(molecule[atom1] - molecule[atom2])
        if distance < bond_threshold and not (atom1[0] == 'H' and atom2[0] == 'H'):
            ax.plot([molecule[atom1][0], molecule[atom2][0]],
                    [molecule[atom1][1], molecule[atom2][1]],
                    [molecule[atom1][2], molecule[atom2][2]], color ='gray')
    
    ax.view_init(elev=20., azim=30)
    
    if not show_axes:
        plt.axis('off')
        
    if _label:
        plt.title(_label)
        
    plt.show()
    
def rotate_view(atom_dict, angle = 0, **kwargs):
    '''
    returns view of other angle without changing original coordinates
    '''
    axis = kwargs.get('axis', 'z')
    get = kwargs.get('get', False)
    label = kwargs.get('label', None)
    flags = kwargs.get('flags', [True, True, True])
    polygons_colors = kwargs.get('polygons_colors', [None, None])
    
    if axis == 'z':
        copy = smt.transform(atom_dict.copy(), smt.rotation_z(angle))
    elif axis == 'y':
        copy = smt.transform(atom_dict.copy(), smt.rotation_y(angle))
    elif axis == 'x':
        copy = smt.transform(atom_dict.copy(), smt.rotation_x(angle))
    
    show(copy, label, flags, polygons_colors)
    
    if get:
        return copy
        


#add r-groups
    #make an atom with symbol "R", "R1, R2," etc... use this as placeholder for each one
    #specify dihedral, make planar besides for now (MM opt will take care of the rest)
    
    #steps: align C1-R1, C2-R2 bonds,
    #       translate to position 
    #       measure dihedral, 
    #       rotation to desired dihedral
    
def make_molecule_union(_molecule1, _molecule2):
    '''
    accepts two molcules as an argument
    and returns a union, !with differently named internal coordinates!
    if adding r groups in a loop, add from the highest numbered atoms down.
    '''
    molecule1 = _molecule1.copy()
    molecule2 = _molecule2.copy()
    
    molecule_union = dict()
    
    atom_counter = 1
    r_counter = 1
    
    for key in molecule1:
        r_group = re.match(r'[R]{1,2}', key)
        if r_group:
            new_key = 'R' + str(r_counter)
            molecule_union[new_key] = molecule1[key]
            r_counter += 1
            continue
        symbol = re.match(r'([A-Za-z]{1,2})', key)[0]
        new_key = symbol + str(atom_counter)
        molecule_union[new_key] = molecule1[key]
        atom_counter += 1
    
    for key in molecule2:
        r_group = re.match(r'[R]{1,2}', key)
        if r_group:
            new_key = 'R' + str(r_counter)
            molecule_union[new_key] = molecule2[key]
            r_counter += 1
            continue
        symbol = re.match(r'([A-Za-z]{1,2})', key)[0]
        new_key = symbol + str(atom_counter)
        molecule_union[new_key] = molecule2[key]
        atom_counter += 1
        
    return molecule_union
    

def add_r_group(_molecule, keys_r_c_d, _r_group, r_keys_r_c_d, spatial_arguments = [0, 1.2]):
    '''
    adds one molecule to another, as specified by arguments.
    accepts molecule, keys for replaced atom, C it's bonded to, atom to measure dihedral by
    and another molecule, called _r_group, and a similar set of keys
    accepts also a list of the following spatial arguments: dihedral angle, bond_length
    
    dihedral angle uses molecule c to r group bond as its axis (directional) so angles are uniquely specified by this
    '''
    
    #ensure all arguments are copied and not references to avoid weird bugs
    #but floats do not have .copy() functionality...
    dihedral_angle = spatial_arguments[0]
    bond_length = spatial_arguments[1]
    
    molecule = _molecule.copy()
    r_group = _r_group.copy()
    
    #display_molecule(make_molecule_union(molecule, r_group))
    
    #display_molecule(molecule)
    #display_molecule(r_group)
    
    m_r = keys_r_c_d[0]
    m_c = keys_r_c_d[1]
    m_d = keys_r_c_d[2] if len(keys_r_c_d) == 3 else None
    
    r_r = r_keys_r_c_d[0]
    r_c = r_keys_r_c_d[1]
    r_d = r_keys_r_c_d[2] if len(r_keys_r_c_d) == 3 else None
    
    #end copy arguments
    
    #translate c atom of r group to origin for linear transformations
    r_group = smt.translate(r_group, -r_group[r_c])
    assert np.allclose(r_group[r_c], np.array([0,0,0]))
    
    #this works, passes the assertion
    
    # make vectors from the bond C to r atoms for each species
    #similarly, make vectors for the c to dihedral atom
    v_mcr = molecule[m_r] - molecule[m_c]
    assert not np.allclose(v_mcr, np.array([0,0,0]))
    if m_d:
        v_mcd = molecule[m_d] - molecule[m_c]
        if np.allclose(v_mcd, np.array([0,0,0])):
            show(molecule)
            print(f'dihedral key: {m_d}  carbon key: {m_c}')
            raise ValueError('molecule dihedral atom is same as molecule carbon atom')
    v_rcr = r_group[r_r] - r_group[r_c]
    assert not np.allclose(v_rcr, np.array([0,0,0]))
    if r_d:
        v_rcd = r_group[r_d] - r_group[r_c]
        assert not np.allclose(v_rcd, np.array([0,0,0]))

    #all good so far
    
    # find angle between vectors, make matrix to make them antiparallel (180deg angle)
    antipar_matrix = smt.antialign_matrix(v_rcr, v_mcr)

    #apply antiparallel matrix to r group
    r_group = smt.transform(r_group, antipar_matrix)
    #display_molecule(make_molecule_union(molecule, r_group))
    #now the molecules are algined for a bond; need to correct the dihedral.
    #first get vectors normal to the plane formed by the atoms near the bond
    #since r_group was transformed, need to re-make v_rcr and v_rcd
    if len(keys_r_c_d) == 3  and len(r_keys_r_c_d) == 3:
        v_rcr = r_group[r_r] - r_group[r_c]
        v_rcd = r_group[r_d] - r_group[r_c]
        
        r_norm = np.cross(v_rcr, v_rcd)
        assert not np.allclose(r_norm, np.array([0,0,0]))

        m_norm = np.cross(v_mcd, v_mcr)
        assert not np.allclose(m_norm, np.array([0,0,0]))
        
        #then set the dihedral angle to zero (if dihedral atoms specified)
        
        align_dihedral = smt.align_matrix(r_norm, m_norm)
        smt.transform(r_group, align_dihedral)
        #display_molecule(make_molecule_union(molecule, r_group))
        #then rotate the r_group about the v_mcr axis by the dihedral_angle
        dihedral_matrix = smt.rotation_about_axis(v_mcr, dihedral_angle)
        r_group = smt.transform(r_group, dihedral_matrix)

    #put r_group c on top of moleucle c (with error checking that this is in fact what is happening)
    assert np.allclose(r_group[r_c], np.array([0,0,0]))
    translation = molecule[m_c].copy() #np.arrays are passed by reference...
    #translate r group into place, in the molecular bond
    new_bond_vector = bond_length * (v_mcr / np.linalg.norm(v_mcr))
    translation += new_bond_vector
    #apply this translation

    r_group = smt.translate(r_group, translation)
    #display_molecule(make_molecule_union(molecule, r_group))
    #delete the atoms given as references for the bond

    del molecule[m_r]

    del r_group[r_r]

    #make a union of the molecule with group bonded to it (molecule must come first so its keys don't change)
    
    return make_molecule_union(molecule, r_group)


def mirror_radially(fragment):
    '''
    mirrors a fragment radially
    '''
    frag1 = fragment.copy()
    frag2 = smt.transform(frag1, smt.rotation_z(np.pi*2/3))
    frag3 = smt.transform(frag1, smt.rotation_z(np.pi*4/3))
    big_frag = make_molecule_union(frag1, frag2)
    big_frag = make_molecule_union(big_frag, frag3)
    return big_frag




def make_radial_group(l_of_l_of_fs):
    list_of_molecules = []
    for fragments in l_of_l_of_fs:
        frag1 = fragments[0].copy()
        frag2 = fragments[1].copy()
        frag3 = fragments[2].copy()
        frag2 = smt.transform(frag2.copy(), smt.rotation_z(np.pi*2/3))
        frag3 = smt.transform(frag3.copy(), smt.rotation_z(np.pi*4/3))
        frag12 = make_molecule_union(frag1.copy(), frag2.copy())
        full_mol = make_molecule_union(frag12.copy(), frag3.copy())
        list_of_molecules.append(full_mol)
        
    return list_of_molecules

def iterate_radial_groups(list_of_radial_fragments):
    list_of_molecules = []
    for combination in itertools.combinations_with_replacement(list_of_radial_fragments, 3):
        frag1 = combination[0].copy()
        frag2 = combination[1].copy()
        frag3 = combination[2].copy()
        frag2 = smt.transform(frag2.copy(), smt.rotation_z(np.pi*2/3))
        frag3 = smt.transform(frag3.copy(), smt.rotation_z(np.pi*4/3))
        frag12 = make_molecule_union(frag1.copy(), frag2.copy())
        full_mol = make_molecule_union(frag12.copy(), frag3.copy())
        list_of_molecules.append(full_mol)
        
    return list_of_molecules

#iterate through r-groups
def iterate_r_groups(list_of_r_groups, rules):
    '''
    argument: list of r-groups with list of rules
    
    returns: list of generated molecules with new r-groups
    '''

#combine radial components in every unique way

#find generic solutions to all of these

#convert hexagonal coordinates to xyz file

def crude_distance(c1,c2):
    #don't use a distance function.
    #do something much cruder.
    '''
    c1 is coords 1 c2 is coords 2
    '''
    return np.abs((c1[0] - c2[0]) + (c1[1] - c2[1]) + (c1[2]-c2[2]))

def find_far_apart_coords(molecule, pivot_key):
    #much better, went from 30s to 2s using this way.
    #check the farthest things in the x and y slice planes.
    pivot = np.array(molecule[pivot_key])
    vectors = {key: np.array(coord) - pivot for key, coord in molecule.items() if key != pivot_key}
    sorted_keys = sorted(vectors, key=lambda key: crude_distance(vectors[key],[0,0,0]))
    return sorted_keys[0], sorted_keys[-1]
    
    
def condense_dict(_atom_dict):
    '''
    fixes missing indices.
    assumes no R-atoms
    '''
    atom_dict = {}
    index = 1
    r_index = 1
    for atom in _atom_dict:
        new_key = re.match('[A-Za-z]{1,2}', atom).group() + str(index)
        if new_key.startswith('R'):
            new_key = 'R' + str(r_index)
            r_index += 1
        else:
            index += 1
        
        atom_dict[new_key] = _atom_dict[atom]

    return atom_dict

def prune_close_atoms(_atom_dict, threshold,**kwargs):
    debug = kwargs.get('debug',False)
    atom_dict = _atom_dict.copy()
    deleted_atoms = []
    for key1 in _atom_dict:
        for key2 in _atom_dict:
            deleted = False
            same_key = key1 == key2
            dist_cond = np.linalg.norm(_atom_dict[key1] - _atom_dict[key2]) < threshold 
            redundant = key2 in deleted_atoms or key1 in deleted_atoms
            same_atom = re.match(r'[A-Za-z]{1,2}', key1)[0] == re.match(r'[A-Za-z]{1,2}', key2)[0]
            if dist_cond and same_atom and not redundant and not same_key:
                del atom_dict[key2]
                deleted_atoms.append(key2)
                deleted = True
            second_is_h_atom = re.match(r'H', key2) and not re.match(r'H', key1)
            if dist_cond and second_is_h_atom and not redundant and not same_key:
                del atom_dict[key2]
                deleted_atoms.append(key2)
                deleted = True
            if deleted and debug:
                show(atom_dict)
                
    return atom_dict


def count_atoms(atom_dict):
    atom_count = 0
    r_count = 0
    for key in atom_dict:
        if key.startswith('R'):
            r_count += 1
        else:
            atom_count += 1
    return atom_count, r_count


def similar(mol1, mol2, threshold):
    for key in mol1:
        if mol2.get(key) is None:
            return False
    first_key = list(mol1.keys())[0]
    key_set = find_far_apart_coords(mol1, first_key)
    
    m1v1 = mol1[key_set[0]]-mol1[first_key]
    m1v2 = mol1[key_set[1]]-mol1[first_key]
    
    m2v1 = mol2[key_set[0]]-mol2[first_key]
    m2v2 = mol2[key_set[1]]-mol2[first_key]
    
    m1_pn = np.cross(m1v1,m1v2)
    m2_pn = np.cross(m2v1,m2v2)
    
    align_planes = smt.align_matrix(m1_pn, m2_pn)
    
    mol2c = mol2.copy()
    
    mol2c = smt.transform(mol2c,align_planes)
    
    m2cv1 = mol2c[key_set[0]]-mol2c[first_key]
    m2cv2 = mol2c[key_set[1]]-mol2c[first_key]
    
    align_directions = smt.align_matrix(m1v1,m2cv1)
    
    mol2c = smt.transform(mol2c,align_directions)
    
    for key in mol1:
        #can make this faster by checking boxes
        #instead of spheres if it's slow
        if np.linalg.norm((mol1[key]-mol2c[key])) > threshold:
            return False
        
    #if got through every key and all within threshold:
    return True


#writes coordinate output
def write_xyz(filename, atom_dict, comment = None):
    with open(filename, 'w') as coordinate_output:
        coordinate_output.write(str(len(atom_dict)))
        if comment is not None:
            coordinate_output.write('\n' + comment + '\n')
        else:
            coordinate_output.write('\n\n')
        for key in atom_dict:
            symbol = re.match(r'[A-Za-z]{1,2}', key)[0]
            
            coordinates = atom_dict[key].copy()
            
            coordinate_string = ' '
            
            for value in coordinates:
                coordinate_string += '   ' + str(value)
                
            line = symbol + coordinate_string
            
            coordinate_output.write(line + '\n')



def write_cif_coords(molecule, filename):
    #incomplete implementation, need to change to fractional basis
    symbol_pattern = r'([A-Za-z]{1,2})'
    with open(filename, 'w') as file:
        index = 0
        for atom in molecule:
            match = re.match(symbol_pattern, atom)
            write_string = match.group(0) + ' '
            write_string += match.group(0) + str(index) + ' '
            write_string += str(molecule[atom][0]) + ' '
            write_string += str(molecule[atom][1]) + ' '
            write_string += str(molecule[atom][2]) + ' '
            write_string += '1' + '\n'
            file.write(write_string)
            index += 1
            
            
            
def convert_out_to_xyz(out_filename, xyz_filename):
    atom_dict = read_out(out_filename)
    write_xyz(xyz_filename, atom_dict, 'converted from ' + out_filename)