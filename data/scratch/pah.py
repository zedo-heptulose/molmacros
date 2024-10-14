#GOALS FOR THIS PROGRAM-
#MOLECULE MAKER, but optimized for PAHs

import numpy as np
import itertools
import structure_math as smt
import structuremaker as sm
import molecule as mol
import os
import re

class Substituent:
    def __init__(self):
        self.atom_dict = {}
        self.bond_length = 1.4
        self.dihedral = 0
        self.site_atom = None
        self.dihedral_atom = None
        self.r_atom = None
        self.bond_vector = np.array([0,0,0])

rotate_dict = lambda d, n: {k : d[(k + n) % len(d)] for k in d}

def rotate_hex(coords, n): 
    for i in range(n%6):
        coord_copy = np.array([coords[1], coords[2], -coords[0]])
        return coord_copy

hex_const = -np.pi*np.sqrt(3)/2

hex_to_cartesian = np.array([[-hex_const,0,hex_const],
                             [0.5,       1,      0.5],
                             [0,         0,        0],])


class RingNode:
    #can't implement this way. use nodes and an adjacency dict
    #or use the coordinate system.
    #or use a graph library
    #hex coords are good but only allow 6 membered rings interally...
    #yeah, let's do a hex grid.
    #can store nodes at each coordinate.

    def __init__(self, ring_size):
        self.ring_size = ring_size
        self.subst_dict = {i : None for i in range(self.ring_size)}
        self.atoms_dict = {i : 'C' for i in range(self.ring_size)}

    def set_subst(self, subst_dict):
        for key in subst_dict:
            self.subst_dict[key] = subst_dict[key]
            
    def set_atoms(self, atoms_dict):
        for key in atoms_dict:
            self.atoms_dict[key] = atoms_dict[key]
            
    def remove_subst(self, key):
        self.subst_dict[key] = None
    
    def remove_atom(self, key):
        self.heteroatoms_dict[key] = None
    
    def rotate(self, angle):
        self.subst_dict = rotate_dict(self.subst_dict, angle)
        self.heteroatoms_dict = rotate_dict(self.heteroatoms_dict, angle)
    
    
    
class PAH:
    def __init__(self):
        self.rings = {} #key is coordinates value is node
    
    def add_ring(self, coords, node):
        self.rings[coords] = node
        add_coords = np.array([0,0,-1])
        for i in range(6):
            if self.rings[coords + rotate_hex(add_coords,i)]:
                self.rings[coords + rotate_hex(add_coords,i)].remove_subst(i%6)
                self.rings[coords + rotate_hex(add_coords,i)].remove_subst((i-1)%6)
        
    def remove_ring(self, coords):
        del self.rings[coords]
        
    def rotate(self, angle):
        self.rings = {coords : rotate_hex(coords, angle) for coords in self.rings}
        for coords in self.rings:
            self.rings[coords].set_subst(rotate_dict(self.rings[coords].subst_dict, angle))
            self.rings[coords].set_atoms(rotate_dict(self.rings[coords].heteroatoms_dict, angle))
            
    def set_subst(self, coords, subst_dict):
        self.rings[coords].set_subst(subst_dict)
    
    def set_atoms(self, coords, heteroatoms_dict):
        self.rings[coords].set_atoms(heteroatoms_dict)
        add_coords = np.array([0,0,-1])

        left_dict = {i + 2 : atom for i, atom in heteroatoms_dict}
        right_dict = {i + 3 : atom for i, atom in heteroatoms_dict}
        
        for i in heteroatoms_dict:
            if self.rings[coords + rotate_hex(add_coords, i%6)]:
                self.rings[coords + rotate_hex(add_coords, i%6)].set_atoms(left_dict)
                
            if self.rings[coords + rotate_hex(add_coords, (i + 1)%6)]:
                self.rings[coords + rotate_hex(add_coords, (i + 1)%6)].set_atoms(right_dict)
                
                
    def fix_hydrogens(self):
        for coordinate in self.rings:
            add_coord = np.array([0,0,-1])
            other_coord = np.array([0,-1,-1])
            for i in range(6):
                c1 = self.rings[coordinate].subst_dict[i] is None
                c2 = self.rings[coordinate].atoms_dict[i] is None
                c3 = self.rings[coordinate + rotate_hex(add_coord, i%6)] is None
                c4 = self.rings[coordinate + rotate_hex(add_coord, (i + 1)%6)] is None
                c5 = self.rings[coordinate + rotate_hex(other_coord, i%6)] is None
                c6 = self.rings[coordinate + rotate_hex(other_coord, (i + 1)%6)] is None
                
                if c1 and c2 and c3 and c4 and c5 and c6:
                    self.rings[coordinate].subst_dict[i] = 'H'
                    
                    
    def make_atom_dict(self):
        index = 0
        molecule = {}
        visited_coords = {}
        for coords in self.rings:
            index += 1
            add_coord = np.array([0,-0.5,-0.5])
            current_keys = {}
            skip = {}
            for i in range(6):
                skip[i] = False
                if self.rings[coords].atoms_dict[i] is None:
                    skip[i] = True
                    continue
                
                atom = self.rings[coords].atoms_dict[i]
                coord = hex_to_cartesian @ (coords + rotate_hex(add_coord, i))
                self.rings[coords].subst_dict[i]
                key = atom + str(index)
                current_keys[i] = key
                molecule[key] = coord
                
            for i in range(6):
                if skip[i] is True or self.rings[coords].subst_dict[i] is None:
                    continue
                
                sub = self.rings[coords].subst_dict[i]                
                molecule = sm.add_r_group(molecule, sub.atom_dict)
                
        prune_close_atoms(molecule, 0.3)
        return molecule
    
    def make_molecule(self):
        molecule = mol.Molecule()
        for coords in self.rings:
            skip = {}
            ring = self.rings[coords]
            for i in range(6):
                skip[i] = False
                if ring.atoms_dict[i] is None:
                    skip[i] = True
                    continue
                atom_symbol = ring.atoms_dict[i]
                coord = hex_to_cartesian @ (coords + rotate_hex(np.array([0,-0.5,-0.5]), i))
                molecule.add_atom(atom_symbol, coord)
                
            for i in range(6):
                if skip[i] is True or ring.subst_dict[i] is None:
                    continue
                mol_keys = [ring.atoms_dict[i], ring.atoms_dict[(i+1)%6]]
                molecule.add_group()
                
    
    def make_molecule(self):
        return mol.Molecule(self.make_atom_dict())
    
    def show(self):
        self.make_molecule().show()
                
                