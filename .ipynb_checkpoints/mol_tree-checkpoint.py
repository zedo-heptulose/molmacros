import molecule as mol
import re

#how to make this work?
#there can be a list of child keys to add to
#and the group_nodes can be stored in a dict,
#where each child is stored with the key that 
# used to add it to the upper

#or the group_nodes can be a tuple with the site they're meant to bond to

#one way tree
#no need to know depth etc

class MolTree:
    class MolNode:
        def __init__(self, molecule = None, lower_keys = None, upper_key = None):
            if molecule:
                self.molecule = molecule.copy()
            self.upper_key = upper_key #key used to add this node to upper's group_nodes
            self.lower_keys = lower_keys #keys of bonding sites to add group_nodes to 
            self.group_nodes = {}

        def copy(self):
            return MolTree.MolNode(self.molecule, self.lower_keys, self.upper_key)

    def __init__(self, **kwargs):
        root = kwargs.get('root', None)
        if root is not None:
            self.root = root
            return root
        
        molecule = kwargs.get('molecule', None)
        if molecule is not None:
            lower_keys = kwargs.get('lower_keys', None)
            upper_key = kwargs.get('upper_key', None)
            self.root = self.MolNode(molecule,lower_keys,upper_key)
            return root
            
        self.root = self.MolNode()
        return root

    def add_molecule(self, molecule, upper=None, lower_keys = None, upper_key=None):
        node = self.MolNode(molecule, lower_keys, upper_key)
        if upper:
            upper.group_nodes.append(node)
            node.upper = upper
        else:
            self.root = node
        return node
    
    def connect(self, upper, lower_key, child, upper_key):
        upper.group_nodes[lower_key] = child
        child.upper_key = upper_key


    def build(self, node=None, **kwargs):
        '''
        This function builds 
        '''
        debug = kwargs.get('debug', False)

        working_mol = node.molecule.copy()

        if debug:
            print(f"Initial working molecule: {working_mol}")

        if node is None:
            return working_mol

        sites = []
        if node.group_nodes:
            for lower_key in node.group_nodes:
                result = self.build(node.group_nodes[lower_key])
                if len(result) != 2:
                    print(f"Error: Expected build() to return 2 values, but got {len(result)}")
                    continue
                lower_mol, low_upper_key = result
                sites.append((lower_key, lower_mol, low_upper_key))
        
        if debug:
            print(f"Sites before ordering: {sites}")

        sites = sorted(sites, key=lambda x : int(re.match(r'(\d+)',x[0]).group()))

        if debug:
            print(f"Sites after ordering: {sites}")

        for site in sites:
            low_key = site[0]
            low_mol = site[1]
            low_upper_key = site[2]
            
            working_mol.active_site = low_key
            low_mol.active_site = low_upper_key
            
            if re.search(r'[A-Za-z]',low_key):
                working_mol = working_mol.add_across_bond(low_mol)
            else:
                working_mol = working_mol.add_group(low_mol)

            if debug:
                print(f"Working molecule after adding site {low_key}: {working_mol}")

        if debug:
            print("Finished building molecule")

        return working_mol, node.upper_key

            #psuedocode:
            #this function returns a molecule.
            #that molecule is made of all the group_nodes joined to the upper.
            
            #first, we call this function on the group_nodes, because it is postorder
            #if there are no group_nodes, simply return the molecule
            #we compile the result of each call into a list of child molecules
            #we add each child molecule to the upper in order of increasing key
            #we return the upper molecule
            
            
# import mol_tree as mt
# import molecule as mol
# import utilities as util
# t = util.get_templates('lib',using_keys = True)
# del t['biphenyl']
# mytree = mt.MolTree()
# mytree.add_molecule(t['benzene'])
# mytree.root.molecule.show()
# mytree.root.lower_keys = ['1','2','3','4','5']
# sub = mytree.root.copy()

# mytree.connect(mytree.root,'1',sub,'6')
# mytree.build(mytree.root,debug=True).show()