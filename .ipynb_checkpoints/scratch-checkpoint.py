# for jupyter scratch work

import numpy as np

import molecule as mol
import utilities as util
import mol_tree as mt

t = util.get_templates('../data/lib/frag')

monomer = util.get_templates('../data/lib/graphenes')
monomer_r = monomer['monomer_R']
monomer_r.active_site = '1'
dimer1 = monomer_r.add_group(monomer_r,dihedral_angle=0)
dimer2 = monomer_r.add_group(monomer_r,dihedral_angle=np.pi)
dimer1.show()
#dimer1.write_xyz('../data/test/dimer1.xyz')
dimer2.show()
#dimer2.write_xyz('../data/test/dimer2.xyz')

import numpy as np
t = util.get_templates('../data/lib/frag')
benzene = t['benzene'].copy()
benzene.active_site = '1'
alkyne = t['alkyne'].copy()
alkyne.active_site = '1'
phenylacetylene = benzene.add_group(alkyne)
phenylacetylene.active_site = '2'
biphenyl = phenylacetylene.add_group(benzene,debug=True)

#need to skew the outer ring up
#and the alkyne down
#lyfe is pane



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

core.write_xyz('../data/test/core.xyz')
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

t = util.get_templates('../data/lib/frag')

monomer = util.get_templates('../data/lib/graphenes')
monomer_r = monomer['monomer_R']
monomer_r.active_site = '1'
dimer1 = monomer_r.add_group(monomer_r,spatial_args = [0-np.pi/6,1.2])
dimer2 = monomer_r.add_group(monomer_r,spatial_args =[np.pi-np.pi/6,1.2])
dimer1.show()
dimer1.write_xyz('../data/test/dimer1.xyz')
dimer2.show()
dimer2.write_xyz('../data/test/dimer2.xyz')

















import molecule as mol
import utilities as util

g = util.get_templates('../data/lib/graphenes')
cmr = g['closed_monomer_R'].copy()
cmr.active_site = '1'

import structure_math as smt
cmr2 = cmr.transform(smt.reflection_xz)

cmr.show()
cmr2.show()

import numpy as np 
closed_dimer1 = cmr.add_group(cmr,dihedral_angle=np.pi,bond_length=1.6,debug=True)
closed_dimer1.show()
closed_dimer1.write_xyz('../data/test/closed_dimer1.xyz')

closed_dimer2 = cmr.add_group(cmr2,dihedral_angle = 0, bond_length=1.6,debug=True)
closed_dimer2.show()
closed_dimer2.write_xyz('../data/test/closed_dimer2.xyz')













import molecule as mol
import utilities as util
import mol_tree as mt #moltree is a pain. I'll probably just make a function to build these manually.
#moltree would, however, be useful for iterating groups. I really just need to change some of the more annoying things about it and it would be nice. eh.

def n_chain(molecule, site1,site2, length,debug=True):
    '''
    this makes the side chains for the big ol molecules.
    '''
    if debug: print('////////////////In n_chain: /////////////////\n\n')
    if length == 1:
        return molecule

    old_site = molecule.active_site
    
    mol2 = molecule.copy()
    if debug:
        print('molecule\n')
        molecule.show()
        molecule.print_sites()
        print('\n\nmol2:\n')
        mol2.show()
        mol2.print_sites()
        
    for i in range(length-1):       
        molecule.active_site = site2
        mol2.active_site = site1
        mol2 = molecule.add_group(mol2)
        if debug:
            print(f'added group {i+1}')
            mol2.show()
    
    molecule.active_site = old_site
    return mol2
    
def chain_squared(molecule, sites1, sites2, length, debug=True):
    '''
    this makes the side chains for the big ol molecules.
    '''
    if length == 1:
        return molecule

    old_site = molecule.active_site
    
    site1 = sites1[0]
    site2 = sites1[1]
    
    side_chain = n_chain(molecule,site1,site2,length,debug=debug)
    
    n_site1 = sites2[0]
    n_site2 = sites2[1]
    
    mol2 = molecule.copy()
    for i in range(length-1):       
        #trial and error lol let's goooo
        mol2 = n_chain(side_chain,n_site1,n_site2,length,debug=debug)
        if debug:
            print(f'added group {i+1}')
            mol2.show()
    
    molecule.active_site = old_site
    return mol2

g = util.get_templates('../data/lib/graphenes',using_keys=True)
cd2 = g['closed_dimer2'].copy()
#cd2.show()


four_chain = n_chain(cd2,'2','4',4,debug=False)
#four_chain.show()
four_chain.write_xyz('../data/test/four_chain.xyz')

sites1 = ('1','3')
sites2 = ('2','4')

for i in range(1,4):
    sixteen_chain = chain_squared(cd2,sites2,sites1,i,debug=False)
    sixteen_chain.replace_r_atoms()
    sixteen_chain.prune_close_atoms()
    sixteen_chain.write_xyz(f'../data/lib/graphenes/closed_dimer_{i}x{i}.xyz')





#/////////////////////////////////////////////////
# THIS IS THE WORK I DID JULY 11TH.
# IT WAS USED TO CREATE THE OPEN POLYMERS.





import molecule as mol
import utilities as util
import mol_tree as mt #moltree is a pain. I'll probably just make a function to build these manually.
#moltree would, however, be useful for iterating groups. I really just need to change some of the more annoying things about it and it would be nice. eh.

def n_chain(molecule, site1,site2, length,dihedral_angle,debug=True):
    '''
    this makes the side chains for the big ol molecules.
    '''
    if debug: print('////////////////In n_chain: /////////////////\n\n')
    if length == 1:
        return molecule

    old_site = molecule.active_site
    
    mol2 = molecule.copy()
    if debug:
        print('molecule\n')
        molecule.show()
        molecule.print_sites()
        print('\n\nmol2:\n')
        mol2.show()
        mol2.print_sites()
        
    for i in range(length-1):       
        molecule.active_site = site2
        mol2.active_site = site1
        mol2 = molecule.add_group(mol2,dihedral_angle=dihedral_angle)
        if debug:
            print(f'added group {i+1}')
            mol2.show()
    
    molecule.active_site = old_site
    return mol2
    
def chain_squared(molecule, sites1, sites2, length, dihedral_angle,debug=True):
    '''
    this makes the side chains for the big ol molecules.
    '''
    if length == 1:
        return molecule

    old_site = molecule.active_site
    
    site1 = sites1[0]
    site2 = sites1[1]
    
    side_chain = n_chain(molecule,site1,site2,length,dihedral_angle,debug=debug)
    
    n_site1 = sites2[0]
    n_site2 = sites2[1]
    
    mol2 = molecule.copy()
    for i in range(length-1):       
        #trial and error lol let's goooo
        mol2 = n_chain(side_chain,n_site1,n_site2,length,np.pi/100,debug=debug)
        if debug:
            print(f'added group {i+1}')
            mol2.show()
    
    molecule.active_site = old_site
    return mol2

g = util.get_templates('../data/lib/graphenes',using_keys=True)
od = g['open_dimer'].copy()
#cd2.show()


# four_chain = n_chain(od,'2','4',4,debug=False)
# #four_chain.show()
# four_chain.write_xyz('../data/test/four_chain.xyz')

sites1 = ('1','3')
sites2 = ('2','4')

for i in range(1,4):
    sixteen_chain = chain_squared(od,sites2,sites1,i,np.pi,debug=False)
    sixteen_chain = sixteen_chain.replace_r_atoms()
    sixteen_chain = sixteen_chain.prune_close_atoms(threshold = 1.0)
    sixteen_chain.write_xyz(f'../data/graphenes/open_dimer_{i}x{i}.xyz')



################ CODE USED TO GENERATE THE TETRARADICAL ##########################def distortion(coords): 
    x = coords[0]
    z = coords[2]
    y = coords[1] + 0.2 * (np.cos(x) + np.cos(z))
    return np.array([x,y,z])

newcore = core.distort(distortion)

newcore = newcore.distort(lambda c : np.array([c[0],c[1]+1,c[2]]), anti_mask=['H7','H69','H120','H30'])

newcore = newcore.distort(lambda c : np.array([c[0],c[1]-1,c[2]]),anti_mask=['H19','H98','H108','H61'])
g = util.get_templates('../data/lib/graphenes',using_keys=True)
core = g['core'].copy()

sp2 = temp['sp2'].copy()
sp2.active_site = '1'

core.active_site = '1'
core = core.add_group(sp2,dihedral_angle=np.pi/4)

core.active_site = '2'
core = core.add_group(sp2,dihedral_angle=2*np.pi/3)

core.active_site = '3'
core = core.add_group(sp2,dihedral_angle=np.pi/4)

core.active_site = '4'
core = core.add_group(sp2,dihedral_angle=np.pi/4)


core.show()


newcore.show()
newcore.write_xyz(f'../data/lib/graphenes/core.xyz')