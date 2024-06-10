
# 

# ///////////////OLD CODE FROM PAHITER. Not likely needed anymore////////////
#  def make_atom_dict(self):
#         index = 0
#         molecule = {}
#         visited_coords = {}
#         for coords in self.rings:
#             index += 1
#             add_coord = np.array([0,-0.5,-0.5])
#             current_keys = {}
#             skip = {}
#             for i in range(6):
#                 skip[i] = False
#                 if self.rings[coords].atoms_dict[i] is None:
#                     skip[i] = True
#                     continue
                
#                 atom = self.rings[coords].atoms_dict[i]
#                 coord = hex_to_cartesian @ (coords + rotate_hex(add_coord, i))
#                 self.rings[coords].subst_dict[i]
#                 key = atom + str(index)
#                 current_keys[i] = key
#                 molecule[key] = coord
                
#             for i in range(6):
#                 if skip[i] is True or self.rings[coords].subst_dict[i] is None:
#                     continue
                
#                 sub = self.rings[coords].subst_dict[i]                
#                 molecule = sm.add_r_group(molecule, sub.atom_dict)
                
#         prune_close_atoms(molecule, 0.3)
#         return molecule
    
    
    




















# /////////////////////////////////////////


import numpy as np
benzene.active_site = '1'
alkyne = t['alkyne'].copy()
alkyne.active_site = '1'
phenylacetylene = benzene.add_group(alkyne)
phenylacetylene.active_site = '2'
biphenyl = phenylacetylene.add_group(benzene,spatial_args=[0.2,1.4])

#need to skew the outer ring up
#and the alkyne down


def distort_y(array):
    x = array[0]
    y = array[1]
    z = array[2]
    
    new_array = [x,y + 0.1*(x**2)*np.exp(-5*z**2),z]
    
    return np.array(new_array)



biphenyl.active_site = '3b'
biphenyl = biphenyl.translate(-biphenyl['C3'])
biphenyl = biphenyl.distort(distort_y)
biphenyl.show()

core = benzene.copy()
core.active_site = '1b'
core=core.add_across_bond(biphenyl)
core.active_site = '3b'
core=core.add_across_bond(biphenyl)
core.active_site = '5b'
core=core.add_across_bond(biphenyl)

core.write_xyz('test/core.xyz')
core.show()

#



#for dimers.
#could just about push these.

#wonder if can optimize with PBCS.

#worth checking.
import numpy as np

import molecule as mol
import utilities as util
import mol_tree as mt

t = util.get_templates('lib/frag')

monomer = util.get_templates('lib/graphenes')
monomer_r = monomer['monomer_R']
monomer_r.active_site = '1'
dimer1 = monomer_r.add_group(monomer_r,spatial_args = [0-np.pi/6,1.2])
dimer2 = monomer_r.add_group(monomer_r,spatial_args =[np.pi-np.pi/6,1.2])
dimer1.show()
dimer1.write_xyz('test/dimer1.xyz')
dimer2.show()
dimer2.write_xyz('test/dimer2.xyz')