import molecule as mol
import os
import numpy as np

def get_templates(directory, **kwargs):
    '''
    '''
    newtemp = {}
    for file in [file for file in os.listdir(directory) if not file.startswith('.')]:
        # print(file)
        #newtemp[file.split('.')[0]] = mol.Molecule(filename=f'{path}/{file}').copy()
        mol_obj = mol.Molecule(filename=f'{directory}/{file}')
        copied_mol = mol_obj.copy()
        newtemp[file.split('.')[0]] = copied_mol
    return newtemp

def get_all_xyz(directory,**kwargs):
    '''
    '''
    mols = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.xyz'):
                basename = os.path.basename(file)
                path = os.path.join(root,file)
                mols.append((basename, mol.Molecule(filename=path)))
    return mols





##### structure building utilities
def chain_groups(monomer, key1, key2, repeats, **kwargs):
    '''
    chain group as a polymer of itself, repeats times
    key1 is the safe key
    key2 is the bad key
    '''
    spargs = kwargs.get('spargs', None)
    cap = kwargs.get('cap', None)
    capkey = kwargs.get('capkey', None)
    cap1 = kwargs.get('cap1', None)
    cap1key = kwargs.get('cap1key', None)
    cap2 = kwargs.get('cap2', None)
    cap2key = kwargs.get('cap2key', None)
    
    chain = monomer.copy()
    chain.active_site = key1
    
    the_cap = None
    if cap1 is not None and cap1key is not None:
        cap1.active_site = cap1key
        the_cap = cap1
    elif cap is not None and capkey is not None:
        cap.active_site = capkey
        the_cap = cap
    if the_cap is not None:
        chain.add_group(the_cap,spatial_args=spargs)

    for i in range(repeats-1):
        new = monomer.copy()
        new.active_site = key2
        new.add_group(chain,spatial_args=spargs)
        chain = new.copy()
    
    chain.active_site = key2
    if cap2 is not None and cap2key is not None:
        the_cap = cap2
        cap2.active_site = cap2key
    elif cap is not None and cap is not None:
        the_cap = cap
        cap.active_site = capkey
    
    chain.add_group(the_cap, spatial_args=spargs)

    return chain



def iterate_chains(monomer, key1, key2, repeats, **kwargs):
    '''
    like chain_groups, but iterates from 1 up to repeats and returns a list
    iterate chains of a monomer, key1 is the safe key, key2 is the bad key
    '''
    chains = []
    for i in range(repeats):
        temp = chain_groups(monomer, key1, key2, i+1, **kwargs)
        chains.append(temp)
        
    return chains

def iterate_groups(cores, groups, **kwargs):
    '''
    accepts a dict of cores,
    a dict of groups to iterate on them,
    and the following optional keyword arguments:
    spatial_args =(float,float) or [float,float]
    key1 
    key2
    verbose (True/False)
    file_out (True/False)
    '''
    verbose = kwargs.get('verbose',False)
    file_out = kwargs.get('file_out',False)
    directory = kwargs.get('directory','.')
    iterations = {}
    for core in cores:
        for group in groups:
            temp = cores[core].copy()
            key = core + '_' + str(temp.active_site) + '_' + group + '_' + str(groups[group].active_site)
            
            if verbose:
                print(f"Generating Molecule With: {key}")
            if file_out:
                temp.write_xyz(os.path.join(directory,key+'.xyz'))
                
            temp.add_group(groups[group],**kwargs)
            iterations[key] = temp
            
    return iterations


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
    
def chain_squared(molecule, sites1, sites2, length1, length2, dihedral_angle,debug=True):
    '''
    this makes the side chains for the big ol molecules.
    '''
    if length1 == 1 and length2 == 1:
        return molecule

    old_site = molecule.active_site
    
    site1 = sites1[0]
    site2 = sites1[1]
    
    side_chain = n_chain(molecule,site1,site2,length1,dihedral_angle,debug=debug)
    
    n_site1 = sites2[0]
    n_site2 = sites2[1]
    
    mol2 = molecule.copy()
    for i in range(length1-1):       
        #trial and error lol let's goooo
        mol2 = n_chain(side_chain,n_site1,n_site2,length2,np.pi/100,debug=debug)
        if debug:
            print(f'added group {i+1}')
            mol2.show()
    
    molecule.active_site = old_site
    return mol2

