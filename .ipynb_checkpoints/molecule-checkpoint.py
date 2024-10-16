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
import molmath as mm
import structure_math as smt
import numpy as np
import os
import shutil
import copy

#persistence and storing keys
import pickle
import json

#xtb harness integration
import sys


path = os.path.abspath('../batch-manager')
if path not in sys.path:
    sys.path.append(path)

import input_files
import job_harness

class Atom:
    def __init__(self):
        self.type = None
        self.isotope_mass = None
        self.charge = 0
        self.spin = 0

INPLACE_DEFAULT = False        

class Molecule:
    def __init__(self, **kwargs):
        self.atom_coords = {}  # maps keys to numpy arrays
        self.num_atoms = 0  # number of non-R atoms
        self.num_r_atoms = 0  # number of R atoms
        self.atom_info = {}  # maps keys to Atom objects
        self.sites = {}  # maps keys to lists of keys
        self.name = None
        self.active_site = None
        self.inplace = INPLACE_DEFAULT
    
        atom_dict = kwargs.get('atom_dict', None)
        filename = kwargs.get('filename', None)
        if atom_dict is not None:
            self.atom_coords = atom_dict
            self.num_atoms = len(atom_dict)
            self.num_r_atoms = 0
        elif filename is not None:
            self.name = os.path.basename(filename).split('.')[0]
            extension = os.path.splitext(filename)[1]
            if extension == '.xyz':
                self.atom_coords, self.num_atoms, self.num_r_atoms = sm.read_xyz(filename)
            elif extension == '.json':
                with open(f'{filename}', 'r') as file:
                    self = self.deserialize_json(json.load(file))
                
            elif extension == '.pkl':
                raise NotImplementedError('.pkl file extension does not work with this class')
                try:
                    with open(filename, 'rb') as file:
                        self = pickle.load(file)
                except Exception as e:
                    print(f"Error during loading: {e}")           
            else:
                raise ValueError('Invalid file type. Must be .xyz or .json')   

    def deserialize_json(self,json_dict):
        jd = json_dict
        self.atom_coords = {key: np.array(jd['atom_coords'][key]) for key in jd['atom_coords']}.copy()
        self.num_atoms = jd['num_atoms']
        self.atom_info = jd['atom_info'].copy()
        self.sites = jd['sites'].copy()
        self.name = jd['name']
        self.active_site = jd['active_site']
        self.inplace = jd['inplace']
        return self
        
    def serialize_dict(self):
        #TODO: IMPLEMENT ATOM INFO
            return {
                'atom_coords': {key : list(self.atom_coords[key]) for key in self.atom_coords},  # maps keys to numpy arrays
                'num_atoms': self.num_atoms,  # number of non-R atoms
                'atom_info': {}, #NOT IMPLEMENTED # maps keys to Atom objects
                'sites': self.sites,  # maps keys to lists of keys
                'name': self.name,
                'active_site': self.active_site,
                'inplace': self.inplace
            }
        
    def write_xyz(self,**kwargs):
        directory = kwargs.get('directory','../data')
        basename = kwargs.get('basename',self.name)
        if basename.endswith('.xyz'):
            basename = basename[:-4]
        full_path = os.path.join(directory,basename)+'.xyz'
        sm_new.write_xyz(full_path,self.atom_coords, comment=None)

    def write_json(self,**kwargs):
        directory=kwargs.get('directory','../data/lib/gen')
        if directory.endswith('/'):
            directory = directory[:-1]
        name = kwargs.get('name',None)
        if not name:
            name = self.name
        if not name:
            raise ValueError('Must have name to write object')
        if '.' in name:
            name = os.path.splitext(name)[0]
        with open(f'{directory}/{name}.json', 'w') as file:
            json.dump(self.serialize_dict(), file)
    

    def __getitem__(self, key):
        return self.atom_coords[key]
    
    def __setitem__(self, key, value):
        self.atom_coords[key] = value
    
    def __delitem__(self, key):
        del self.atom_coords[key]
        
    def __iter__(self):
        return iter(self.atom_coords)
    
    
    def _instance(self, is_inplace):
        if is_inplace:
            return self
        else:
            return self.copy()
            
    def _inplace_or_copy(method):
        def wrapper(self, *args, **kwargs):
            inplace = kwargs.get('inplace', self.inplace)
            temp = self._instance(inplace)
            result = method(temp, *args, **kwargs)
            return result
        return wrapper

    def add_bonding_site(self,name,atoms,**kwargs):
        '''
        accepts a name for a site and a list of atom keys
        expects atoms in this order:
        atom to be replaced
        atom to bond to
        atom to use for dihedral angle
        (second atom for dihedral angle, for add_across_bond sites)

        Functions concerning bonding sites are inplace by default
        '''
        if self.sites.get(name,None) and not kwargs.get('update',True):
            raise ValueError('Tried to overwrite site with update=False')
        self.sites[name] = atoms
        return self
        
    @_inplace_or_copy
    def gfn2opt(self,**kwargs):
        scratch_dir = '.scratch'
        os.mkdir(scratch_dir)
        job_basename = self.name
        self.write_xyz(directory=scratch_dir)

        xtb_builder = input_files.xTBInputBuilder()
        xtb_job = xtb_builder.change_params({
            'xyz_directory': scratch_dir,
            'xyz_file': job_basename+'.xyz',
            'job_basename': job_basename,
            'write_directory': scratch_dir,
            }).build()
        xtb_job.debug = True
        xtb_job.create_directory()
        
        xtb_harness = job_harness.xTBHarness()
        #xtb_harness.mode = 'direct'
        xtb_harness.job_name = job_basename
        xtb_harness.directory = scratch_dir
        ret_val = xtb_harness.MainLoop()
        if ret_val == 0:
            print('Successful optimization')

            self.atom_coords, self.num_atoms, self.num_r_atoms = sm.read_xyz(os.path.join(scratch_dir,'xtbopt.xyz'))
            shutil.rmtree(scratch_dir)
            return self

        else:
            raise RuntimeError('Geometry optimization failed, check input')



    def remove_bonding_site(self,name,**kwargs):
        del self.sites[name]
        return self
        
    def reset_sites(self,**kwargs):
        self.sites = {}
        return self
    
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

    @_inplace_or_copy
    def transform(self, matrix, **kwargs):
        self.atom_coords = smt.transform(self.atom_coords, matrix)
        return self

    @_inplace_or_copy
    def translate(self, vector, **kwargs):
        self.atom_coords = smt.translate(self.atom_coords, vector)
        return self

    @_inplace_or_copy
    def add_atom(self, atom_symbol, coords,**kwargs):
        '''
        not working at the moment
        '''
        psuedo_atom = {atom_symbol: coords}
        self.atom_coords = sm.make_molecule_union(self.atom_coords, psuedo_atom)
        
        if atom_symbol.startswith('R'):
            self.num_r_atoms += 1
        else:
            self.num_atoms += 1
        return self

    @_inplace_or_copy
    def remove_atoms(self, atoms, **kwargs):
        if type(atoms) is list:
            for atom in atoms:
                del self[atom]
                if atom.startswith('R'):
                    self.num_r_atoms -= 1
                else:  
                    self.num_atoms -= 1
        elif type(atoms) is str:
            del self[atoms]

        else:
            raise ValueError('invalid type. remove_atoms accepts str or list(str)')
        if kwargs.get('condense',False):
            self.atom_coords = sm.condense_dict(self.atom_coords)
 
        return self
    
    @_inplace_or_copy
    def union_with(self, other_molecule, **kwargs):
        self.atom_coords = sm_new.make_molecule_union(self.atom_coords, other_molecule.atom_coords)
        self.num_atoms += other_molecule.num_atoms
        self.num_r_atoms += other_molecule.num_r_atoms
        return self

    
    @_inplace_or_copy
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
        if not (type(key) is list or type(key) is tuple):
            key = self.site(key)
        if not (type(g_key) is list or type(g_key) is tuple):
            g_key = group.site(g_key)
        
        self = sm_new.add_group(self, key, group, g_key, **kwargs)
        self.num_atoms, self.num_r_atoms = sm.count_atoms(self.atom_coords)
        
        return self
        
    @_inplace_or_copy
    def add_across_bond(self, group,**kwargs):
        # '''
        # accepts either a tuple of four atom keys
        # or a single variable as a key from the atom's dict
        # adds rings across bonds, deletes nearby hydrogens
        # '''
        key = self.active_site
        g_key = group.active_site
        
        if not (type(key) is list or type(key) is tuple):
            key = self.sites[key]

        if not (type(g_key) is list or type(g_key) is tuple):
            g_key = group.sites[g_key]

        if kwargs.get('reverse',False):
            temp = g_key[0]
            g_key[0] = g_key[1]
            g_key[1] = temp
            
        self.atom_coords = sm_new.add_across_bond(self.atom_coords, key, group.atom_coords, g_key, **kwargs)
        
        db = kwargs.get('debug', False)
        clean = kwargs.get('clean', True)
        if clean:
            self = self.prune_close_atoms(**kwargs)        
            #self.atom_coords = sm.condense_dict(self.atom_coords)
        return self
       
       
    #TODO: this shouldn't change the order of atoms.
    @_inplace_or_copy
    def replace_atom(self, key,new_key,**kwargs):
        '''
        #NEXT: MAKE THIS STABLE.
        THEN CAN FREELY TEST DOPING
        '''   
        self.atom_coords[new_key] = self[key]
        del self[key]
        
        debug = kwargs.get('debug',False)
        if debug: print(f'Atom Coords:\n{self.atom_coords}')
        self.atom_coords = sm.condense_dict(self.atom_coords)
    
        return self

    @_inplace_or_copy
    def replace_r_atoms(self, **kwargs):
        '''
        accepts no arguments, replaces r atoms with H atoms.
        '''
        self = sm_new.replace_r_atoms(self,**kwargs)
        self.num_r_atoms = 0
        self.num_atoms = len(self.atom_coords)
        return self
        
    @_inplace_or_copy
    def prune_close_atoms(self, **kwargs):
        '''
        accepts the following kwargs:
        threshold (default 0.5)
        debug (default False)
        inplace (default self.inplace, usually False)
        '''
        threshold = kwargs.get('threshold', 0.5)
        db = kwargs.get('debug', False)
        self.atom_coords = sm.prune_close_atoms(self.atom_coords, threshold, debug=db)
        self.num_atoms, self.num_r_atoms = sm.count_atoms(self.atom_coords)
        return self

    @_inplace_or_copy
    def distort(self, function, **kwargs):
        self.atom_coords = sm_new.distort(self.atom_coords, function, **kwargs)
        return self

    @_inplace_or_copy
    def collapse_atom_numbers(self,**kwargs):
        '''
        gets rid of gaps in atom numbering, while
        preserving order
        '''
        self.atom_coords = sm_new.sort_keys(self.atom_coords)
        return self

    @_inplace_or_copy
    def align_to_plane(self,atoms,**kwargs):
        '''
        args here should be a list of atoms to align to a plane.
        If one atom, do nothing
        if two atoms, align those to the x axis
        if three atoms, align those to the xy axis
        kwargs:
        axis1,
        axis2
        ='x','y','z'
        plane='xy','xz','yz'
        '''
        debug = kwargs.get('debug',False)
        
        if debug: 
            print('molecule before edit:')
            self.show()
        #deine axes based on kwargs. xy is default
        #MAN this syntax highlighting is ugly.
        ax1 = np.array([1,0,0])
        ax2 = np.array([0,1,0])

        if len(atoms) == 0 or len(atoms) == 1:
            print('Zero or one atoms passed to align_to_plane')
            return self
        if len(atoms) == 2:
            #align to ax1
            #TODO: IMPLEMENT THIS
            pass

        if len(atoms) == 3:
            #align plane normals of molecu;e and plane
            a1coord = self[atoms[0]]
            a2coord = self[atoms[1]]
            a3coord = self[atoms[2]]
            mol_v1 = a2coord - a1coord
            mol_v2 = a3coord - a1coord
            mol_plnorm = np.cross(mol_v1,mol_v2)
            ext_plnorm = np.cross(ax1,ax2)

            alignment_mat = mm.align_matrix(ext_plnorm,mol_plnorm)

            self = self.transform(alignment_mat)
            
            if debug:
                print(f'a1coord : {a1coord}')
                print(f'a2coord : {a2coord}')
                print(f'a3coord : {a3coord}')
                print(f'mol_v1 : {mol_v1}')
                print(f'mol_v2 : {mol_v2}')
                print(f'mol_plnorm : {mol_plnorm}')
                print(f'ext_plnorm : {ext_plnorm}')
                print(f'alignment_mat : {alignment_mat}')
                print(f'molecule after transformation:')
                self.show()      
        return self

    @_inplace_or_copy
    def rotate_atoms(self,atom1,index,**kwargs):
        '''
        rotates atom1 to specified index
        '''
        self.atom_coords = sm_new.rotate_keys(self.atom_coords,atom1,index)
        return self
        
    
    def is_similar(self, other_molecule):
        return sm.similar(self.atom_coords, other_molecule.atom_dict)

    @_inplace_or_copy
    def bs(self,bonding_site,**kwargs):
        self.active_site = bonding_site
        return self
        
    def copy(self):
        new_molecule = Molecule()
        
        assert type(self.atom_coords) is dict
        new_molecule.atom_coords = copy.deepcopy(self.atom_coords)
        
        assert type(self.sites) is dict or self.sites is None
        new_molecule.sites = copy.deepcopy(self.sites)
        
        assert type(self.atom_info) is dict or self.atom_info is None
        new_molecule.atom_info = copy.deepcopy(self.atom_info)
        
        assert type(self.num_atoms) is int
        new_molecule.num_atoms = self.num_atoms
        
        assert type(self.num_r_atoms) is int or self.num_r_atoms is None
        new_molecule.num_r_atoms = self.num_r_atoms
        
        assert type(self.name) is str or self.name is None
        new_molecule.name = self.name

        assert type(self.active_site) is str or self.active_site is None
        new_molecule.active_site = self.active_site

        assert type(self.inplace) is bool or self.inplace is None
        new_molecule.inplace = self.inplace
        
        return new_molecule
