import molecule as mol
import os

def get_templates(directory, **kwargs):
    '''
    get all template molecules in the directory
    as a dict. 
    if present, reads keys.txt to get sites of the molecules; 
    these are appended as a dict.
    '''
    using_keys = kwargs.get('using_keys', True)
    if using_keys:
        keys_file = kwargs.get('keys_file', 'keys.txt')
        keys = read_keys_file(os.path.join(directory, keys_file))
    molecules = {}
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.xyz'):
                basename = os.path.splitext(file)[0]
                molecule = mol.Molecule(filename=os.path.join(root, file))
                molecules[basename] = (molecule)
                if using_keys:
                    molecule.sites = keys[basename]
    return molecules


def read_keys_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = {}
    current_key = None

    for line in lines:
        line = line.strip()
        if line.startswith('%'):
            current_key = line[1:].strip()
            data[current_key] = {}
        elif line == 'end':
            current_key = None
        elif current_key is not None and ':' in line:
            name, keys = line.split(":", 1)
            data[current_key][name.strip()] = tuple(keys.strip().split())

    return data



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

