###
'''
goal-
decouple components of molecule
make a molecule class
clean up weird stuff
'''
###

###again.... need to fix the basic dihedral function.
#it has floating point error issues at the very least
##### DIHEDRAL IS BROKEN, CHECK IT OUT #####

import structuremaker as sm
import sm_new
import structure_math as smt
import numpy as np
import os
import copy


class Atom:
    def __init__(self):
        self.type = None
        self.isotope_mass = None
        self.charge = 0
        self.spin = 0

INPLACE_DEFAULT = False        

class Molecule:
    def __init__(self, **kwargs):
        atom_dict = kwargs.get('atom_dict', None)
        filename = kwargs.get('filename', None)
        if atom_dict is not None:
            self.atom_coords = atom_dict
            self.num_atoms = len(atom_dict)
            self.num_r_atoms = 0
        elif filename is not None:
            extension = os.path.splitext(filename)[1]
            if extension == '.xyz':
                self.atom_coords, self.num_atoms, self.num_r_atoms = sm.read_xyz(filename)
            else:
                raise ValueError('Invalid file type. Must be .xyz')
        else:
            self.atom_coords = {}  # maps keys to numpy arrays
            self.num_atoms = 0  # number of non-R atoms
            self.num_r_atoms = 0  # number of R atoms
            
        self.atom_info = {}  # maps keys to Atom objects
        self.sites = {}  # maps keys to lists of keys
        self.name = None
        self.active_site = None
        self.inplace = INPLACE_DEFAULT
    
    def __getitem__(self, key):
        return self.atom_coords[key]
    
    def __setitem__(self, key, value):
        self.atom_coords[key] = value
    
    def __delitem__(self, key):
        del self.atom_coords[key]
        
    def __iter__(self):
        return iter(self.atom_coords)
    
    def print_sites(self):
        for key in self.sites:
            print(key + ':' + str(self.sites[key]) )
    
    def site (self, key):
        return self.sites[key]
    
    def show(self, **kwargs):
        sm_new.show(self.atom_coords, **kwargs)
        show_sites = kwargs.get('show_sites', False)
        if show_sites:
            self.print_sites()
        
    def get_atom_info(self, atom):
        return self.atom_info[atom]
    
    def get_coords(self, atom):
        return self.atom_coords[atom]
    
    def get_vector(self, atom1, atom2):
        return self.atom_coords[atom2] - self.atom_coords[atom1]
    
    def get_plane_normal(self, atom1, atom2, atom3):
        v1 = self.get_vector(atom1, atom2)
        v2 = self.get_vector(atom2, atom3)
        if np.allclose(np.cross(v1, v2), np.array([0,0,0])):
            return np.array([0,0,0])
        return np.cross(v1, v2) / np.linalg.norm(np.cross(v1, v2))
    
    def transform(self, matrix, **kwargs):
        inplace = kwargs.get('inplace', self.inplace)
        
        _molecule = self.instance(inplace)
        
        _molecule.atom_coords = smt.transform(self.atom_coords, matrix)
        
        return _molecule
    
    def translate(self, vector, **kwargs):
        inplace = kwargs.get('inplace', self.inplace)

        _molecule = self.instance(inplace)
        _molecule.atom_coords = smt.translate(self.atom_coords, vector)
        return _molecule
  
    
    def add_atom(self, atom_symbol, coords,**kwargs):
        '''
        not working at the moment
        '''
        #atom_data = kwargs.get('atom_data', None)
        inplace = kwargs.get('inplace', self.inplace)
        _molecule = self.instance(inplace)
        
        psuedo_atom = {atom_symbol: coords}
        _molecule.atom_coords = sm.make_molecule_union(_molecule.atom_coords, psuedo_atom)
        
        if atom_symbol.startswith('R'):
            _molecule.num_r_atoms += 1
        else:
            _molecule.num_atoms += 1
        return _molecule
    
    def instance(self, is_inplace):
        if is_inplace:
            return self
        else:
            return self.copy()
    
    def remove_atom(self, atom, **kwargs):
        inplace = kwargs.get('inplace', self.inplace)
        _molecule = self.instance(inplace)
            
        del _molecule[atom]
        _molecule.atom_coords = sm.condense_dict(_molecule.atom_coords)
        if atom.startswith('R'):
            _molecule.num_r_atoms -= 1
        else:  
            _molecule.num_atoms -= 1
        return _molecule
    
    def union_with(self, other_molecule, **kwargs):
        inplace = kwargs.get('inplace', self.inplace)
        _molecule = self.instance(inplace)
        
        _molecule.atom_coords = sm.make_molecule_union(_molecule.atom_coords, other_molecule.atom_coords)
        _molecule.num_atoms += other_molecule.num_atoms
        _molecule.num_r_atoms += other_molecule.num_r_atoms
        return _molecule
        
    def add_group(self, group, **kwargs):
        '''
        either accepts name of set of atoms for substitution,
        or the actual atoms as a list or tuple
        '''
        key = kwargs.get('key', None)
        if key is None:
            key = self.active_site
        if key is None:
            raise ValueError('add_group called with no key')
        
        g_key = kwargs.get('g_key', None)
        if g_key is None:
            g_key = group.active_site
        if g_key is None:
            raise ValueError('add_group called with no key')
        
        inplace = kwargs.get('inplace', self.inplace)
        
        _molecule = self.instance(inplace)
        
        if not (type(key) is list or type(key) is tuple):
            key = _molecule.site(key)
        if not (type(g_key) is list or type(g_key) is tuple):
            g_key = group.site(g_key)
        
        _molecule = sm_new.add_group(_molecule, key, group, g_key, **kwargs)
    
        _molecule.num_atoms, _molecule.num_r_atoms = sm.count_atoms(_molecule.atom_coords)
        
        return _molecule
        


    def add_across_bond(self, group,**kwargs):
        # '''
        # accepts either a tuple of four atom keys
        # or a single variable as a key from the atom's dict
        # adds rings across bonds, deletes nearby hydrogens
        # '''
        # simplest way to do this with fewest rabbit holes and least confusion
        # just accept dicts of up to four keys for this one
        # where the last two get pruned, if they are there,
        # and the first two are the bond to add across
        
        inplace = kwargs.get('inplace', self.inplace)
        _molecule = self.instance(inplace)
        
        key = self.active_site
        g_key = group.active_site
        
        if not (type(key) is list or type(key) is tuple):
            key = _molecule.sites[key]

        if not (type(g_key) is list or type(g_key) is tuple):
            g_key = group.sites[g_key]
            
        _molecule.atom_coords = sm_new.add_across_bond(_molecule.atom_coords, key, group.atom_coords, g_key, **kwargs)
        
        db = kwargs.get('debug', False)
        clean = kwargs.get('clean', True)
        if clean:
            _molecule = _molecule.prune_close_atoms(**kwargs)        
            _molecule.atom_coords = sm.condense_dict(_molecule.atom_coords)
        return _molecule
       
       
       
    def replace_atom(self, key,new_key,**kwargs):
        '''
        #NEXT: MAKE THIS STABLE.
        THEN CAN FREELY TEST DOPING
        '''   
        inplace = kwargs.get('inplace',self.inplace)
        _molecule = self.instance(inplace)
        
        _molecule.atom_coords[new_key] = _molecule[key]
        del _molecule[key]
        
        debug = kwargs.get('debug',False)
        if debug: print(f'Atom Coords:\n{_molecule.atom_coords}')
        _molecule.atom_coords = sm.condense_dict(_molecule.atom_coords)
    
        return _molecule

    def replace_r_atoms(self, **kwargs):
        '''
        accepts no arguments, replaces r atoms with H atoms.
        '''
        inplace = kwargs.get('inplace', self.inplace)
        _molecule = self.instance(inplace)
        _molecule = sm_new.replace_r_atoms(_molecule)
        _molecule.num_r_atoms = 0
        _molecule.num_atoms = len(_molecule.atom_coords)
        return _molecule
        
    
    def prune_close_atoms(self, **kwargs):
        threshold = kwargs.get('threshold', 0.5)
        db = kwargs.get('debug', False)
        inplace = kwargs.get('inplace', self.inplace)
        _molecule = self.instance(inplace)
        _molecule.atom_coords = sm.prune_close_atoms(_molecule.atom_coords, threshold, debug=db)
        _molecule.num_atoms, _molecule.num_r_atoms = sm.count_atoms(_molecule.atom_coords)
        return _molecule

    def distort(self, function, **kwargs):
        inplace = kwargs.get('inplace', self.inplace)
        _molecule = self.instance(inplace)
        _molecule.atom_coords = sm_new.distort(_molecule.atom_coords, function, **kwargs)
        return _molecule
        
        
    def is_similar(self, other_molecule):
        return sm.similar(self.atom_coords, other_molecule.atom_dict)
    
    def write_xyz(self, filename, comment = None):
        sm.write_xyz(filename, self.atom_coords, comment=None)
        
    def copy(self):
        new_molecule = Molecule()
        if type(self.atom_coords) is dict:
            new_molecule.atom_coords = copy.deepcopy(self.atom_coords)
        if type(self.sites) is dict:
            new_molecule.sites = copy.deepcopy(self.sites)
        if type(self.atom_info) is dict:
            new_molecule.atom_info = copy.deepcopy(self.atom_info)
            
        new_molecule.num_atoms = self.num_atoms
        new_molecule.num_r_atoms = self.num_r_atoms
        new_molecule.sites = self.sites
        new_molecule.name = self.name
        new_molecule.active_site = self.active_site
        new_molecule.inplace = self.inplace
        return new_molecule